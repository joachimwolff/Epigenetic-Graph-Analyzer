# main.py


# import necessary libraries
import numpy
import pandas
import hicexplorer
import hicmatrix
import networkx
import argparse
import gzip
import pyBigWig
import pandas as pd
import pyranges as pr
from scipy.sparse import dok_matrix
# This is the main script for the epigenetic_graph_analyzer project.



def parse_arguments():
    parser = argparse.ArgumentParser(description='epigenetic_graph_analyzer')
    parser.add_argument('interactions_sample_1', type=str, nargs='+', help='Path to the text file with interactions')
    parser.add_argument('interactions_sample_2', type=str, nargs='+', help='Path to the text file with interactions')

    parser.add_argument('histone_data_sample_1', type=str, nargs='+', help='Path(s) to the histone data file(s) (bigwig)')
    parser.add_argument('histone_data_sample_2', type=str, nargs='+', help='Path(s) to the histone data file(s) (bigwig)')

    parser.add_argument('rna_expression_data_sample_1', type=str, nargs='+', help='Path(s) to the RNA expression data file(s) (bigwig)')
    parser.add_argument('rna_expression_data_sample_2', type=str, nargs='+', help='Path(s) to the RNA expression data file(s) (bigwig)')

    parser.add_argument('atac_expression_data_sample_1', type=str, nargs='+', help='Path(s) to the ATAC data file(s) (bigwig)')
    parser.add_argument('atac_expression_data_sample_2', type=str, nargs='+', help='Path(s) to the ATAC data file(s) (bigwig)')

    parser.add_argument('promoter_region_file', type=str, help='Path to the promoter region file')
    parser.add_argument('enhancer_region_file', type=str, help='Path to the enhancer region file')

    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use (default: 1)')
    parser.add_argument('-v', '--version', action='version', version='epigenetic_graph_analyzer 1.0')
    return parser.parse_args()

# provide a function to read bigwig files
def read_bigwig_file(file_path):
    # read the bigwig file
    bw = pyBigWig.open(file_path)
    values = bw.values()
    bw.close()
    return values




def create_interaction_matrix(interactions_path, promoter_regions_path, enhancer_regions_path):
    # create a matrix of the interactions
    # the matrix must be of the form (promoter, enhancer, interaction)
    # each row is a promoter which was enriched in the capture Hi-C experiment
    # each column is a promoter or enhancer from the enhancer region file and promoter region file
    # the interaction is the number of interactions between the promoter and enhancer
    pyrange_promoter = pr.read_bed(promoter_regions_path)
    pyrange_enhancer = pr.read_bed(enhancer_regions_path)
    pd_interactions = pd.read_table(interactions_path, header=None, sep='\t|,', compression='gzip')
    pd_interactions.columns = ["Chromosome", "Start", "End", "Chromosome_OE", "Start_OE", "End_OE", "interactions"]

    # match promoter locations with bait locations
    pyrange_interactions_promoters = pyrange_interactions.intersect(pyrange_promoter)

    # match enhancer locations with OE locations
    pyrange_interactions_promoters.columns = ["Chromosome_Bait", "Start_Bait", "End_Bait", "Chromosome", "Start", "End", "interactions"]
    pyrange_interactions_promoters_enhancer = pyrange_interactions_promoters.intersect(pyrange_enhancer)

    # match promoter locations with OE locations
    pyrange_interactions_promoters_bait_promoters_OE = pyrange_interactions_promoters.intersect(pyrange_promoter)

    # stack both pyranges of enhancer OE and promoter OE
    stacked_pyranges = pr.concat([pyrange_interactions_promoters_enhancer, pyrange_interactions_promoters_bait_promoters_OE])
    stacked_pyranges.columns = ["Chromosome", "Start", "End", "Chromosome_OE", "Start_OE", "End_OE", "interactions"]

    # get unique locations of baits
    unique_bait_locations = stacked_pyranges.as_df()[['Chromosome', 'Start', 'End']].drop_duplicates().sort_values(by=['Chromosome', 'Start', 'End']).reset_index(drop=True)

    # get unique locations of OE
    unique_OE_locations = stacked_pyranges.as_df()[['Chromosome_OE', 'Start_OE', 'End_OE']].drop_duplicates().sort_values(by=['Chromosome_OE', 'Start_OE', 'End_OE']).reset_index(drop=True)

    length_baits = len(unique_bait_locations)
    length_OE = len(unique_OE_locations)
    matrix = dok_matrix((length_baits, length_OE), dtype=float)

    # Create a dictionary to map bait locations to matrix indices
    bait_index_map = {tuple(row): index for index, row in unique_bait_locations.iterrows()}

    # Create a dictionary to map OE locations to matrix indices
    OE_index_map = {tuple(row): index for index, row in unique_OE_locations.iterrows()}

    for index, row in stacked_pyranges.as_df().iterrows():
        bait_location = tuple(row[['Chromosome', 'Start', 'End']])
        OE_location = tuple(row[['Chromosome_OE', 'Start_OE', 'End_OE']])
        
        bait_index = bait_index_map[bait_location]
        OE_index = OE_index_map[OE_location]
        
        matrix[bait_index, OE_index] = row["interactions"]

    return matrix 

def main():
    args = parse_arguments()

    # TODO: Add your code here
    # read promoter and enhancer regions
    # the promoter and enhancer together form the nodes of the graph aka the entries of a matrix
    # there must be a mapping of position in the matrix to the promoter and enhancer regions
    # read the interaction file and 

    # read the promoter region file
    promoter_regions = read_bed_file(args.promoter_region_file)
    enhancer_regions = read_bed_file(args.enhancer_region_file) 

    # read the interactions file
    interactions_sample_1 = read_interaction_file(args.interactions_sample_1)
    interactions_sample_2 = read_interaction_file(args.interactions_sample_2)

    # read the histone data files
    histone_data_sample_1 = read_bigwig_file(args.histone_data_sample_1)
    histone_data_sample_2 = read_bigwig_file(args.histone_data_sample_2)

    # read the rna expression data files
    rna_expression_data_sample_1 = read_bigwig_file(args.rna_expression_data_sample_1)
    rna_expression_data_sample_2 = read_bigwig_file(args.rna_expression_data_sample_2)

    # read the atac expression data files
    atac_expression_data_sample_1 = read_bigwig_file(args.atac_expression_data_sample_1)
    atac_expression_data_sample_2 = read_bigwig_file(args.atac_expression_data_sample_2)

    # create a matrix of the interactions
    # the matrix must be of the form (promoter, enhancer, interaction)
    # each row is a promoter which was enriched in the capture Hi-C experiment
    # each column is a promoter or enhancer from the enhancer region file and promoter region file
    # the interaction is the number of interactions between the promoter and enhancer






if __name__ == '__main__':
    main()
