import argparse
from collections import OrderedDict
from multiprocessing import Process, Queue
import time
from scipy.sparse import csr_matrix
import traceback
from icecream import ic
import pyranges as pr
import networkx as nx

import gzip
import pandas as pd
import numpy as np

import hicmatrix.HiCMatrix as hm
import scipy.sparse as sp

import logging
log = logging.getLogger(__name__)

def parse_arguments():
    parser = argparse.ArgumentParser(description='epigenetic_graph_analyzer')
    parser.add_argument('--interactionData', '-id', type=str, nargs='+', help='Path to the washU or cool file containing the interaction data. List the replicates in the order they should be combined.', required=True)
    parser.add_argument('--interactionFormat', '-if',choices=['washU', 'cool'], help='Format of the interaction data (washU or cool)', required=True)

    parser.add_argument('--chromosomeSizes','-cs', type=str, help='Path to the file containing chromosome sizes', required=True)
    parser.add_argument('--digestionSites', '-ds', type=str, help='Path to the file containing digestion sites, alternative the resolution of the matrix as integer.', required=True)
    parser.add_argument('--promoterFile', '-p', type=str, help='Path to the file containing promoter data', required=True)
    parser.add_argument('--outputFile', '-o', type=str, help='Path to the output file where the graph will be saved', required=True)
    parser.add_argument('--interChromosomalInteractions', '-ic', action='store_true', help='Allow inter-chromosomal interactions (default: False)')
    parser.add_argument('--threshold', '-th', type=float, default=5.0, help='Threshold for interaction to be considered (default: 5.0)')

    parser.add_argument('-t', '--threads', type=int, default=8, help='Number of threads to use (default: 1)')
    parser.add_argument('-v', '--version', action='version', version='epigenetic_graph_analyzer 1.0')
    return parser.parse_args()


def createCutIntervals(chromosome_size_file, resolution):
    chrom_sizes = OrderedDict()
    size_genome = 0
    with open(chromosome_size_file, 'r') as file:
        file_ = True
        while file_:
            file_ = file.readline().strip()
            if file_ != '':
                line_split = file_.split('\t')
                chrom_sizes[line_split[0]] = int(line_split[1])
                size_genome += int(line_split[1])

    cut_intervals = []
    if isinstance(resolution, int):
        for chromosome in chrom_sizes:
            for interval in range(0, chromosome[1], resolution):
                cut_intervals.append(tuple([chromosome[0], interval,
                                            min(chromosome[1], interval + resolution), 1.0]))
    else:
        restriction_cut_sites = pd.read_table(resolution, sep='\t')
        for row in restriction_cut_sites.iterrows():
            
            if row[1]['Chromosome'] in chrom_sizes:
                if int(row[1]['Start']) != int(row[1]['End']):
                    cut_intervals.append(tuple([str(row[1]['Chromosome']), int(row[1]['Start']), int(row[1]['End']), 1.0]))
    
    return cut_intervals

def importInteractionsFromWashU(list_of_washU_files, cut_intervals, allow_inter_chr_interactions=False, compression='gzip', pNumberOfProcesses=8):
    
    hic_matrix_list = []
    matrix_merged_csr = None
                
    for i, path_interactions in enumerate(list_of_washU_files):
        dtype_dict = {0:str, 1:int, 2:int, 3:str, 4:int, 5:int, 6:float}
        pd_interactions = pd.read_table(path_interactions, header=None, sep='\t|,', compression=compression, dtype=dtype_dict)
        
        hic_matrix = createMatrix(pd_interactions, cut_intervals, allow_inter_chr_interactions, pNumberOfProcesses)
        hic_matrix_list.append(hic_matrix)

        # merge matrices
        if matrix_merged_csr is None:
            matrix_merged_csr = hic_matrix_list[0].matrix.tocsr()
        else:
            matrix_merged_csr += hic_matrix_list[i].matrix.tocsr()
    
    # compute an intersection matrix to make sure that the values are present in all matrices i.e. all replicates
    intersection_result_matrix = hic_matrix_list[0].matrix
            
    if len(hic_matrix_list) > 1:
        for hic_matrix in hic_matrix_list[1:]:
            intersection_matrix = intersection_result_matrix.multiply(hic_matrix.matrix).tocsr()
            intersection_matrix.eliminate_zeros()
            non_zero_ids = intersection_matrix.nonzero()
            data_merged_matrix = np.asarray(intersection_result_matrix[non_zero_ids]).flatten()
            data_hic_matrix = np.asarray(hic_matrix.matrix[non_zero_ids]).flatten()
            data_intersection = data_merged_matrix + data_hic_matrix
            intersection_result_matrix = csr_matrix((data_intersection, non_zero_ids),shape=intersection_result_matrix.shape)

    hic_matrix_list.append(matrix_merged_csr)
    hic_matrix_list.append(intersection_result_matrix / len(list_of_washU_files))

    return hic_matrix_list

def process_interactions(interactions, hic_matrix, allow_inter_chr_interactions, pQueue):
    hic_lil_matrix = sp.lil_matrix(hic_matrix.matrix.shape)
    for chromosome_1, start_1, end_1, chromosome_2, start_2, end_2, value in zip(interactions[0], interactions[1], interactions[2], interactions[3], interactions[4], interactions[5], interactions[6]):
        if value == 0:
            continue
        if not allow_inter_chr_interactions:
            if chromosome_1 != chromosome_2:
                continue
        start_1 = start_1 + 1
        end_1 = end_1 - 1
        start_2 = start_2 + 1
        end_2 = end_2 - 1
        try:
            bin_id_1 = hic_matrix.getRegionBinRange(chromosome_1, start_1, end_1)
            bin_id_2 = hic_matrix.getRegionBinRange(chromosome_2, start_2, end_2)
            hic_lil_matrix[bin_id_1, bin_id_2] = value
        except Exception as exp:
            print(traceback.print_exc())
            ic(chromosome_1, start_1, end_1, value)
            ic(chromosome_2, start_2, end_2, value)
            pQueue.put('Fail: ' + str(exp) + traceback.format_exc())
            return
    
    pQueue.put(hic_lil_matrix)
    return

def run_parallel(interactions, hic_matrix, allow_inter_chr_interactions, pNumberOfProcesses=8):
    chunk_size = len(interactions) // pNumberOfProcesses
    queue = [None] * pNumberOfProcesses
    process = [None] * pNumberOfProcesses
    all_data_collected = False
    thread_done = [False] * pNumberOfProcesses
    interactions_thread = [None] * pNumberOfProcesses

    for i in range(pNumberOfProcesses):

        if i < pNumberOfProcesses - 1:
            interactions_thread[i] =  interactions.iloc[(i * chunk_size) : ((i + 1) * chunk_size)]
           
        else:
            interactions_thread[i] = interactions.iloc[(i * chunk_size):]
           
        queue[i] = Queue()
        process[i] = Process(target=process_interactions, kwargs=dict(
                interactions=interactions_thread[i],
                hic_matrix=hic_matrix,
                allow_inter_chr_interactions=allow_inter_chr_interactions,
                pQueue=queue[i]
            )
        )
        process[i].start()

    fail_flag = False
    fail_message = ''
    while not all_data_collected:
        for i in range(pNumberOfProcesses):
            if queue[i] is not None and not queue[i].empty():
                queue_data = queue[i].get()
                if (type(queue_data) == str) and ('Fail' in queue_data):
                    fail_flag = True
                    fail_message = queue_data
                else:
                    hic_matrix.matrix += queue_data

                queue[i] = None
                process[i].join()
                process[i].terminate()
                process[i] = None
                thread_done[i] = True
                interactions_thread[i] = None
        all_data_collected = True
        for thread in thread_done:
            if not thread:
                all_data_collected = False
        time.sleep(1)

    if fail_flag:
        return fail_message, None
    
    hic_matrix.matrix = hic_matrix.matrix.tocsr()
    return hic_matrix

def createMatrix(pd_interactions, cut_intervals, allow_inter_chr_interactions, pNumberOfProcesses=8):
    hic_matrix_csr = sp.lil_matrix((len(cut_intervals), len(cut_intervals)))
    hic_matrix = hm.hiCMatrix()
    hic_matrix.setMatrix(hic_matrix_csr, cut_intervals)
    hic_matrix = run_parallel(pd_interactions, hic_matrix, allow_inter_chr_interactions, pNumberOfProcesses)

    return hic_matrix

def createSampleInteractionMatrix(path_promoter, chic_matrix):
    pyrange_promoter = pr.read_bed(path_promoter)

    matrix_list = []
    # matrix_list_replicates = [[] for _ in range(len(chic_matrix[0:-2]))]
    bait_list = []
    for index, row in pyrange_promoter.as_df().iterrows():
        bait_location = row[['Chromosome', 'Start', 'End']]
        
        try:
            view_point_start = chic_matrix[0].getRegionBinRange(bait_location[0], bait_location[1], bait_location[1])
            if view_point_start is None:
                
                view_point_start = chic_matrix[0].getRegionBinRange(bait_location[0], bait_location[1]+1, bait_location[1]+1)
                
                if view_point_start is None:
                    ic(traceback.print_exc())
                    ic(bait_location, row)
                    break
            # -1 to get the last matrix which is the mean of all replicates
            matrix_list.append(chic_matrix[-1][view_point_start[0], :])
            # for i, matrix in enumerate(chic_matrix[0:-2]):
                # matrix_list_replicates[i].append(matrix.matrix[view_point_start[0], :])
            bait_list.append([row[['Name']].str.split(',')[0][0], str(row['Chromosome']) + ":" + str(row['Start']) + "-"+ str(row['End'])])
        except IndexError as exp:
            pass
    # Assuming you have a list of CSR matrices: matrix_list
    sample_matrix = sp.vstack(matrix_list)
    # for i, matrix in enumerate(chic_matrix[0:-2]):
        # matrix_list_replicates[i] = sp.vstack(matrix_list_replicates[i])
    return sample_matrix, bait_list

def createGraphAndFindCommunities(sample_matrix, bait_list, cut_intervals, threshold):
    G = nx.Graph()
    number_of_samples = sample_matrix.shape[0]
    bait_ids, OE_ids = sample_matrix.nonzero()
    for i, (bait, OE) in enumerate(zip(bait_ids, OE_ids)):
        if sample_matrix.data[i] < threshold:
            continue
        if (bait != OE):
            G.add_weighted_edges_from([(bait, OE+number_of_samples, sample_matrix.data[i])])
            G.nodes[OE + number_of_samples]['colorType'] = "orange"
            G.nodes[OE + number_of_samples]['position'] = cut_intervals[OE][0] + ":" + str(cut_intervals[OE][1]) + "-" +  str(cut_intervals[OE][2])
       
    for bait in bait_ids:
        if bait in G.nodes:
            G.nodes[bait]['name'] = bait_list[bait][0]
            G.nodes[bait]['colorType'] = "green"
            G.nodes[bait]['position'] = bait_list[bait][1]

    partition_list = nx.community.louvain_communities(G, seed=42)

    # enumerate partion creating a mapping from gene name to partition
    gene_name_to_partition_mapping = {}
    for partition_id, partition in enumerate(partition_list):
        for node in partition:
            G.nodes[node]['partition'] = partition_id
            if "name" in G.nodes[node]:
                gene_name_to_partition_mapping[G.nodes[node]['name']] = partition_id

    return G, partition_list, gene_name_to_partition_mapping

def compressedFile(pFilename):
    """
    Determines if a file is compressed or not
    """
    with open(pFilename, 'rb') as f:
        if f.read(2) == b'\x1f\x8b':
            return 'gzip'
        return None

def main():
    args = parse_arguments()
    log.info("Reading interaction data...")
    cut_intervals = createCutIntervals(args.chromosomeSizes, resolution=args.digestionSites)

    if args.interactionFormat == 'washU':

        interaction_matrices = importInteractionsFromWashU(args.interactionData, cut_intervals, args.interChromosomalInteractions, compression=compressedFile(args.interactionData[0]), pNumberOfProcesses=args.threads)

    else:
        pass
    
    log.info("Create interaction matrix...")

    sample_matrix, bait_list = createSampleInteractionMatrix(args.promoterFile, interaction_matrices)

    log.info("Create graph...")

    graph, partition, gene_name_to_partition_mapping = createGraphAndFindCommunities(sample_matrix, bait_list, cut_intervals=cut_intervals, threshold=args.threshold)

    nx.write_graphml(graph, args.outputFile + "_graph.graphml")
    df = pd.DataFrame(list(gene_name_to_partition_mapping.items()), columns=['Gene', 'Community'])
    df.to_csv(args.outputFile + "_communities.csv", index=False)

    partition_df = pd.DataFrame({
        'Index': range(len(partition)),
        'Content': [','.join(map(str, partition)) for partition in partition]
    })

    # Set the 'Index' column as the index of the DataFrame
    partition_df.set_index('Index', inplace=True)

    # Write the DataFrame to a CSV file
    partition_df.to_csv(args.outputFile + "_partitions.csv")
    
