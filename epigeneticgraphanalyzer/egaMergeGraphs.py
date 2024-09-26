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

import os
import plotly.graph_objects as go
# import plotly.express as px

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from adjustText import adjust_text 
import matplotlib.patheffects as path_effects
import copy

import logging
log = logging.getLogger(__name__)

def parse_arguments():
    parser = argparse.ArgumentParser(description='epigenetic_graph_analyzer')
    parser.add_argument('-g1', '--graph1', type=str, required=True, help='Path to the first graphml file')
    parser.add_argument('-g2', '--graph2', type=str, required=True, help='Path to the second graphml file')

    parser.add_argument('-c1', '--communities1', type=str, required=True, help='Path to the first communities file')
    parser.add_argument('-c2', '--communities2', type=str, required=True, help='Path to the second communities file')

    parser.add_argument('-p1', '--partitions1', type=str, required=True, help='Path to the first partitions file')
    parser.add_argument('-p2', '--partitions2', type=str, required=True, help='Path to the second partitions file')

    parser.add_argument('-pr', '--promoters', type=str, required=True, help='Path to the promoters file')
    parser.add_argument('-r', '--rnaseq', type=str, required=True, help='Path to the RNAseq DESeq file')
    parser.add_argument('-m', '--manorm', type=str, nargs='+', required=True, help='List of paths to the MANorm ChIP-seq files')
    parser.add_argument('-mn', '--manorm_names', type=str, nargs='+', required=True, help='List of names for the MANorm ChIP-seq files')

    parser.add_argument('-p', '--pvalue', type=float, default=0.05, help='P-value threshold for RNA-seq analysis')
    parser.add_argument('-l', '--level', type=int, default=100, help='Number of up and downregulated genes to be considered')
    parser.add_argument('-o', '--outputFolder', type=str, required=True, help='Path to the output folder')

    parser.add_argument('-t', '--threads', type=int, default=8, help='Number of threads to use (default: 1)')
    parser.add_argument('-v', '--version', action='version', version='epigenetic_graph_analyzer 1.0')
    return parser.parse_args()

def getCommunituesFromDESeq2Ranks(df_deseq_list, gene_name_to_partition_mapping, name_extension='_1'):
    cluster = {}
    for df_deseq in df_deseq_list:
        for gene_name in df_deseq['Gene']:
            partition_id = None
            if gene_name in gene_name_to_partition_mapping:
                partition_id = gene_name_to_partition_mapping[gene_name]
            elif gene_name + name_extension in gene_name_to_partition_mapping:
                partition_id = gene_name_to_partition_mapping[gene_name + name_extension]
            if partition_id:
                if partition_id not in cluster:
                    cluster[partition_id] = []
                data_rna = df_deseq[df_deseq['Gene'] == gene_name].values.tolist()
                if len(data_rna) != 0:
                    data_rna = data_rna[0]
                cluster[partition_id].append((gene_name, data_rna))
    return cluster

def jaccardSimiliarityCommunities(community_1, community_2, number_of_baits):
    intersection_communities_nodes = []
    intersection_communities_jaccard = np.zeros((len(community_1), len(community_2)))
    intersection_communities_bait_jaccard = np.zeros((len(community_1), len(community_2)))
    for i in range(len(community_1)):
        intersection_communities_nodes.append({})
        for j in range(len(community_2)):
            intersection = community_1[i].intersection(community_2[j])

            if len(intersection) > 0:
                intersection_communities_nodes[i][j] = intersection
                intersection_communities_jaccard[i, j] = len(intersection) / len(community_1[i].union(community_2[j]))
                bait_intersection = {x for x in intersection if x <= number_of_baits}
                bait_community_1 = {x for x in community_1[i] if x <= number_of_baits}
                bait_community_2 = {x for x in community_2[j] if x <= number_of_baits}
                intersection_communities_bait_jaccard[i, j] = len(bait_intersection) / len(bait_community_1.union(bait_community_2))

    return intersection_communities_nodes, intersection_communities_jaccard, intersection_communities_bait_jaccard

# iterate throught a graph of networkx and copy the nodes. Check for the nodes if they have been existing in the other graph too. if yes, compute the log2Foldchange between them, if not, simply color them green if new in subgraph2, and purple if only in subgraph1
def mergeSubgraphs(pSubgraph_1, pSubgraph_2, RNA_genes, chip_data_list_pyrange):
    G = nx.Graph()
    G.add_nodes_from(copy.deepcopy(list(pSubgraph_1.nodes(data=True))))
    for node in pSubgraph_2.nodes:
        if node in G.nodes:
            degree_node_1 = pSubgraph_1.degree[node]
            degree_node_2 = pSubgraph_2.degree[node]

            G.nodes[node]["degree_log2FC"] = np.nan_to_num(np.log2(degree_node_2/degree_node_1))
            G.nodes[node]["degree_log2FC_binary"] = 1 if np.nan_to_num(np.log2(degree_node_2/degree_node_1)) >= 0 else -1

        else:
            # node is in condition 2 but not in condition 1
            G.add_nodes_from([(node, copy.deepcopy(pSubgraph_2.nodes[node]))])
            # G.nodes[node] = pSubgraph_2.nodes[node]
            G.nodes[node]["degree_log2FC"] = np.nan_to_num(np.log2(pSubgraph_2.degree[node]/1)) 
            G.nodes[node]["degree_log2FC_binary"] = 1 if np.nan_to_num(np.log2(pSubgraph_2.degree[node]/1)) >= 0 else -1

        if 'position' in pSubgraph_2.nodes[node]:
            position = pSubgraph_2.nodes[node]['position']
            for chip_data, chip_name in chip_data_list_pyrange:
                intersections = intersect(chip_data, position)
                if len(intersections) > 0:
                    G.nodes[node][chip_name] = intersections.df.values[0][4]

        if 'name' in G.nodes[node]:
            for genes in RNA_genes:
                # ic(genes)
                # ic(G.nodes[node]['name'], genes)
                if G.nodes[node]['name'].split('_')[0]  == genes[0]:
                    G.nodes[node]['RNA log2FC'] = genes[1][2]
    for node in pSubgraph_1.nodes:
        if node not in pSubgraph_2.nodes:
            G.nodes[node]["degree_log2FC"] = np.nan_to_num(np.log2(1/pSubgraph_1.degree[node]))
            G.nodes[node]["degree_log2FC_binary"] = 1 if np.nan_to_num(np.log2(1/pSubgraph_1.degree[node])) >= 0 else -1

            if 'name' in G.nodes[node]:
                for genes in RNA_genes:
                    if G.nodes[node]['name'].split('_')[0] == genes[0]:
                        G.nodes[node]['RNA log2FC'] = genes[1][2]

            if 'position' in pSubgraph_1.nodes[node]:
                position = pSubgraph_1.nodes[node]['position']
                for chip_data, chip_name in chip_data_list_pyrange:
                    intersections = intersect(chip_data, position)
                    if len(intersections) > 0:
                        G.nodes[node][chip_name] = intersections.df.values[0][4]

    for edge in pSubgraph_1.edges:
        if edge in pSubgraph_2.edges:
            weight_1 = copy.deepcopy(pSubgraph_1.get_edge_data(edge[0], edge[1])['weight'])
            weight_2 = copy.deepcopy(pSubgraph_2.get_edge_data(edge[0], edge[1])['weight'])
            edge = (edge[0], edge[1], weight_2 / weight_1)
            G.add_weighted_edges_from([edge])
            G[edge[0]][edge[1]]['log2FC'] = np.log2(weight_2 / weight_1)
            G[edge[0]][edge[1]]['log2FC_binary'] = 1 if np.log2(weight_2 / weight_1) >= 0 else -1

        else:
            edge = (edge[0], edge[1], 1 / pSubgraph_1.get_edge_data(edge[0], edge[1])['weight'])
            G.add_weighted_edges_from([edge])
            G[edge[0]][edge[1]]['log2FC'] = np.log2(1 / pSubgraph_1.get_edge_data(edge[0], edge[1])['weight'])
            G[edge[0]][edge[1]]['log2FC_binary'] = 1 if np.log2(1 / pSubgraph_1.get_edge_data(edge[0], edge[1])['weight']) >= 0 else -1
            
    for edge in pSubgraph_2.edges:
        if edge not in pSubgraph_1.edges:
            edge = (edge[0], edge[1], pSubgraph_2.get_edge_data(edge[0], edge[1])['weight'] / 1)
            G.add_weighted_edges_from([edge])
            G[edge[0]][edge[1]]['log2FC'] = np.log2(pSubgraph_2.get_edge_data(edge[0], edge[1])['weight'] / 1)
            G[edge[0]][edge[1]]['log2FC_binary'] = 1 if np.log2(pSubgraph_2.get_edge_data(edge[0], edge[1])['weight'] / 1) >= 0 else -1

    return G



def writePlotlyGraphs(cluster_of_regulated_genes, output_path, graph_sample_1, partition_of_graph_sample_1, graph_sample_2, partition_of_graph_sample_2, jaccardSimilarityCommunities, sample_1_name_prefix, sample_2_name_prefix, chip_data_list_pyrange, pBaitSimilarity=0.5):
    analyse_result_df_bait_similarity_high = None
    analyse_result_df_bait_similarity_low = None
    
    for partition_id, genes in cluster_of_regulated_genes.items():
        # Check if the path exists
        if not os.path.exists(os.path.join(output_path, str(partition_id))):
            # Create the path
            os.makedirs(os.path.join(output_path, str(partition_id)))
        # ic(partition_id)
        # ic(genes)
        # ic(graph_sample_1.nodes)
        # ic(partition_of_graph_sample_1[partition_id])
        # ic(partition_of_graph_sample_2[partition_id])

        writePlotlyGraph(partition_id, os.path.join(output_path, str(partition_id)), graph_sample_1.subgraph(partition_of_graph_sample_1[partition_id]), genes, sample_1_name_prefix)
        intersection = jaccardSimilarityCommunities[0][partition_id]
        
        for id in intersection:
            # writePlotlyGraph(id, os.path.join(output_path, str(partition_id)), graph_sample_2.subgraph(partition_of_graph_sample_2[id]), genes, sample_2_name_prefix)
           
            merged_graph = mergeSubgraphs(graph_sample_1.subgraph(partition_of_graph_sample_1[partition_id]), graph_sample_2.subgraph(partition_of_graph_sample_2[id]), genes, chip_data_list_pyrange)
            ic(merged_graph)
            ic(merged_graph.nodes)
            analyseMergedGraphResult =  createSiginficantTestingTable(merged_graph, output_path, 3, chip_data_list_pyrange, partition_id, id)
            createGraphsBasedOnRNA(merged_graph, 3, os.path.join(output_path, str(partition_id)), str(partition_id) + '_' + str(id) , [x[1] for x in chip_data_list_pyrange], jaccardSimilarityCommunities[1][partition_id, id], jaccardSimilarityCommunities[2][partition_id, id])
            if jaccardSimilarityCommunities[2][partition_id, id] >= pBaitSimilarity:
                if analyse_result_df_bait_similarity_high is None:
                    analyse_result_df_bait_similarity_high = analyseMergedGraphResult
                else:
                    analyse_result_df_bait_similarity_high = pd.concat([analyse_result_df_bait_similarity_high, analyseMergedGraphResult], axis=0)
            else:
                if analyse_result_df_bait_similarity_low is None:
                    analyse_result_df_bait_similarity_low = analyseMergedGraphResult
                else:
                    analyse_result_df_bait_similarity_low = pd.concat([analyse_result_df_bait_similarity_low, analyseMergedGraphResult], axis=0)
            writePlotlyGraphMerged(partition_id, id, os.path.join(output_path, str(partition_id)), merged_graph, genes, sample_1_name_prefix + '_' + sample_1_name_prefix)
        
    # return analyse_result_df_bait_similarity_high, analyse_result_df_bait_similarity_low
    if analyse_result_df_bait_similarity_high is not None:
        analyse_result_df_bait_similarity_high.to_csv(os.path.join(output_path, "results_merged_high.csv"))
    if analyse_result_df_bait_similarity_low is not None:
        analyse_result_df_bait_similarity_low.to_csv(os.path.join(output_path, "results_merged_low.csv"))
    
    return analyse_result_df_bait_similarity_high, analyse_result_df_bait_similarity_low

def createPyRangeObject(path, name, log2FCThreshold=1, p_value_threshold=0.05):
    df_MAnorm = pd.read_csv(path, sep='\t')

    # Display the dataframe
    df_MAnorm_pvalue = df_MAnorm[df_MAnorm['P_value'] <= p_value_threshold]
    df_MAnorm_pvalue_mvalue = df_MAnorm_pvalue[(df_MAnorm_pvalue['M_value'] >= log2FCThreshold) | (df_MAnorm_pvalue['M_value'] <= -log2FCThreshold)]
    df_MAnorm_pvalue_mvalue.rename(columns={"chr":'Chromosome', "start" : "Start", 'end':"End"}, inplace=True)
    py_MAnorm = pr.PyRanges(df_MAnorm_pvalue_mvalue)
    return (py_MAnorm, name)

def position_to_pyranges(chrom, start, end):
    return pr.PyRanges(chromosomes=[chrom], starts=[start], ends=[end])

def createSiginficantTestingTable(pMergeGraph, outputPath, numberOfSteps, chip_data_list_pyrange, id_sample_1, id_sample_2):
    # search for nodes with logRNA field
    # get all connected 1..n steps nodes until numberOfSteps is reached
    
    connected_nodes, connected_edges = getConnectedNodes(pMergeGraph, numberOfSteps)

    results = []

    visited_nodes = set()
    result_for_dfferential_dict = {}
    # Collect data for each bait node
    for bait_node in connected_edges:
        # ic(bait_node)
        # result_for_dfferential_dict = {}


        for i, level_edges in enumerate(connected_edges[bait_node]):
            visited_nodes = set()
            interaction_logFC = []
            m_values = {}
            m_values_bait = {}


            for edge in level_edges:
                # ic(edge)
                node1, node2 = edge
                connected_node = node2 #if node1 == bait_node else node1
                

                if connected_node not in visited_nodes:
                    visited_nodes.add(connected_node)
                

                    for chip_name in chip_data_list_pyrange:
                        if chip_name[1] in pMergeGraph.nodes[bait_node]:
                            if chip_name[1] not in m_values_bait:
                                m_values_bait[chip_name[1]] = []
                            m_values_bait[chip_name[1]].append(pMergeGraph.nodes[bait_node][chip_name[1]])

                            if 'log2FC' in pMergeGraph[node1][node2]:

                                if 'positive interaction_' + chip_name[1] not in m_values_bait:
                                    m_values_bait['positive interaction_' + chip_name[1]] = []
                                if 'negative interaction_' + chip_name[1] not in m_values_bait:
                                    m_values_bait['negative interaction_' + chip_name[1]] = []
                                
                                if pMergeGraph[node1][node2]['log2FC'] > 0:
                                    m_values_bait['positive interaction_' + chip_name[1]].append(pMergeGraph.nodes[bait_node][chip_name[1]])
                                else:
                                    m_values_bait['negative interaction_' + chip_name[1]].append(pMergeGraph.nodes[bait_node][chip_name[1]])
                    # Collect chip data
                    for chip_name in chip_data_list_pyrange:
                        if chip_name[1] in pMergeGraph.nodes[connected_node]:
                                # If there is an intersection, add the m-value to m_values
                            if chip_name[1] not in m_values:
                                m_values[chip_name[1]] = []
                            m_values[chip_name[1]].append(pMergeGraph.nodes[connected_node][chip_name[1]])
                            if 'log2FC' in pMergeGraph[node1][node2]:
                                if 'positive interaction_' + chip_name[1] not in m_values:
                                    m_values['positive interaction_' + chip_name[1]] = []
                                if 'negative interaction_' + chip_name[1] not in m_values:
                                    m_values['negative interaction_' + chip_name[1]] = []
                                
                                if pMergeGraph[node1][node2]['log2FC'] > 0:
                                    m_values['positive interaction_' + chip_name[1]].append(pMergeGraph.nodes[connected_node][chip_name[1]])
                                else:
                                    m_values['negative interaction_' + chip_name[1]].append(pMergeGraph.nodes[connected_node][chip_name[1]])
                    # Collect interaction logFC
                    if 'log2FC' in pMergeGraph[node1][node2]:
                        interaction_logFC.append(pMergeGraph[node1][node2]['log2FC'])

            # Compute statistics
            if interaction_logFC:
                positive_values = [x for x in interaction_logFC if x > 0]
                negative_values = [x for x in interaction_logFC if x < 0]

                mean_positive = sum(positive_values) / len(positive_values) if positive_values else 0
                mean_negative = sum(negative_values) / len(negative_values) if negative_values else 0
                count_positive = len(positive_values)
                count_negative = len(negative_values)

                count_negative_compute = count_negative if count_negative > 0 else count_negative + 0.00001
                count_positive_compute = count_positive if count_positive > 0 else count_positive + 0.00001


                logFC_negative_positive = np.log2(count_positive_compute / count_negative_compute)
            

                # Collect the result
                result = {
                    'bait_node': bait_node,
                    'name' : pMergeGraph.nodes[bait_node]['name'],
                    'step': i,
                    'mean_positive_logFC': mean_positive,
                    'mean_negative_logFC': mean_negative,
                    'count_positive_logFC': count_positive,
                    'count_negative_logFC': count_negative,
                    'logFC_count' : logFC_negative_positive,
                    'cluster_id' : str(id_sample_1) + '_' + str(id_sample_2),
                    'RNA_log2FC' : pMergeGraph.nodes[bait_node]['RNA log2FC'],
                    'position' : pMergeGraph.nodes[bait_node]['position']
                }

                # print(m_values)
                if m_values:
                    for chip_name, m_value in m_values.items():
                        # ic(chip_name, m_value)
                        positive_values = [x for x in m_value if x > 0]
                        negative_values = [x for x in m_value if x < 0]
                        
                        mean_positive = sum(positive_values) / len(positive_values) if positive_values else 0
                        mean_negative = sum(negative_values) / len(negative_values) if negative_values else 0
                        count_positive = len(positive_values)
                        count_negative = len(negative_values)

                        count_negative_compute = count_negative if count_negative > 0 else count_negative + 0.00001
                        count_positive_compute = count_positive if count_positive > 0 else count_positive + 0.00001


                        logFC_negative_positive = np.log2(count_positive_compute / count_negative_compute)
                                
                        # Add the sum of m_value to the result dictionary with chip_name as the key
                        result[chip_name + '_mean_positive_OE'] = mean_positive
                        result[chip_name + '_count_positive_OE'] = count_positive
                        result[chip_name + '_mean_negative_OE'] = mean_negative
                        result[chip_name + '_count_negative_OE'] = count_negative
                        result[chip_name + '_logFC_count_OE'] = logFC_negative_positive

                if m_values_bait:
                    for chip_name, m_value in m_values_bait.items():
                        # ic(chip_name, m_value)
                        positive_values = [x for x in m_value if x > 0]
                        negative_values = [x for x in m_value if x < 0]
                        
                        mean_positive = sum(positive_values) / len(positive_values) if positive_values else 0
                        mean_negative = sum(negative_values) / len(negative_values) if negative_values else 0
                        count_positive = len(positive_values)
                        count_negative = len(negative_values)

                        count_negative_compute = count_negative if count_negative > 0 else count_negative + 0.00001
                        count_positive_compute = count_positive if count_positive > 0 else count_positive + 0.00001


                        logFC_negative_positive = np.log2(count_positive_compute / count_negative_compute)
                                
                        # Add the sum of m_value to the result dictionary with chip_name as the key
                        result[chip_name + '_mean_positive_bait'] = mean_positive
                        result[chip_name + '_count_positive_bait'] = count_positive
                        result[chip_name + '_mean_negative_bait'] = mean_negative
                        result[chip_name + '_count_negative_bait'] = count_negative
                        result[chip_name + '_logFC_count_bait'] = logFC_negative_positive
                results.append(result)
      

    # Convert results to a pandas DataFrame
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(outputPath ,str(id_sample_1) ,'dataframe_' + str(id_sample_1) + '_' + str(id_sample_2) +  '.csv'), index=False)

    return df

# def normalize_and_color(nodes, attribute):
#     values = [data[attribute] for _, data in nodes]
#     norm_values = [np.nan_to_num((value - min(values)) / (max(values) - min(values))) for value in values]
#     colors = px.colors.sequential.Plasma
#     node_color = [colors[int(norm * (len(colors) - 1))] for norm in norm_values]
#     return node_color

def plotGraph(pGraph, pNodeColor, pEdgeColor, pBaitName, pOutputPath, pPartitionString, pJaccardSimilarityCommunities=None, pJaccardSimilarityBaitCommunities=None, pFigureSize=(8, 6), pFileFormat='pdf', pSeed=42):
    G=pGraph
    # levels = nx.get_node_attributes(G, pNodeColor)

    # Set a default value for nodes without the 'level' attribute
    labels = {node: G.nodes[node].get('name', "") for node in G.nodes}

    # Step 3: Use a fixed seed for the layout
    fixed_seed = pSeed
    pos = nx.spring_layout(G, seed=fixed_seed)

    # Step 4: Normalize the node attribute 'level' using the seismic colormap and TwoSlopeNorm
    levels = nx.get_node_attributes(G, pNodeColor)
    default_node_value = 0.0
    all_levels = {node: levels.get(node, default_node_value) for node in G.nodes}

    vmin_nodes = min(all_levels.values())
    vmax_nodes = max(all_levels.values())

    # Ensure vmin, vcenter, and vmax are in ascending order
    if vmin_nodes >= 0:
        vmin_nodes = 0-0.001
        vmax_nodes = vmax_nodes if vmax_nodes > 0 else 0 + 0.001  # Set vmax to a positive value if it's 0 or less
    elif vmax_nodes <= 0:
        vmax_nodes = 0+0.001
        vmin_nodes = vmin_nodes if vmin_nodes < 0 else 0 - 0.001 # Set vmin to a negative value if it's 0 or more

    if not vmin_nodes < 0 < vmax_nodes:
        ic(vmin_nodes, vmax_nodes)
    norm_nodes = mcolors.TwoSlopeNorm(vmin=vmin_nodes, vcenter=0, vmax=vmax_nodes)
    node_colors = [cm.seismic(norm_nodes(all_levels[node])) for node in G.nodes]

    # Step 5: Normalize the edge attribute 'weight' using the seismic colormap and TwoSlopeNorm
    weights = nx.get_edge_attributes(G, pEdgeColor)
    default_edge_value = 0.0
    all_weights = {edge: weights.get(edge, default_edge_value) for edge in G.edges}

    vmin_edges = min(all_weights.values())
    vmax_edges = max(all_weights.values())

    # Ensure vmin, vcenter, and vmax are in ascending order
    if vmin_edges >= 0:
        vmin_edges = 0-0.001
        vmax_edges = vmax_edges if vmax_edges > 0 else 0 + 0.001  # Set vmax to a positive value if it's 0 or less
    elif vmax_edges <= 0:
        vmax_edges =  0+0.001
        vmin_edges = vmin_edges if vmin_edges < 0 else 0 - 0.001  # Set vmin to a negative value if it's 0 or more


    norm_edges = mcolors.TwoSlopeNorm(vmin=vmin_edges, vcenter=0, vmax=vmax_edges)
    edge_colors = [cm.seismic(norm_edges(all_weights[edge])) for edge in G.edges]

    # Step 6: Draw the graph with the fixed seed-based positions and the seismic colormap-based colors
    edge_widths = [(abs(weights[edge]) + 0.001) * 2 for edge in G.edges]

    # Step 7: Draw the graph with black edge lines on nodes and custom edge widths
    plt.figure(figsize=pFigureSize)
    nx.draw(
        G, pos=pos, node_color=node_colors, edge_color=edge_colors, node_size=200, font_size=15,
        edgecolors='black', linewidths=2, width=edge_widths  # Set edge widths based on absolute weights
    )
    # Step 8: Draw labels with `adjustText` and apply a white outline to black text
    texts = []

    for node, (x, y) in pos.items():
        # Apply the offset to the text position
        text = plt.text(x , y, labels[node], fontsize=15, ha='center', va='center', color='black')
        text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='white'), path_effects.Normal()])
        texts.append(text)
    # Automatically adjust text positions to avoid overlap
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

    # Add smaller colorbars with titles
    sm_nodes = plt.cm.ScalarMappable(cmap=cm.seismic, norm=norm_nodes)
    sm_nodes.set_array([])
    cbar_nodes = plt.colorbar(sm_nodes, ax=plt.gca(), orientation='vertical', shrink=0.4, pad=0.02, aspect=20)
    cbar_nodes.set_label('Nodes: log2(FC) {}'.format(pNodeColor), labelpad=10)

    sm_edges = plt.cm.ScalarMappable(cmap=cm.seismic, norm=norm_edges)
    sm_edges.set_array([])
    cbar_edges = plt.colorbar(sm_edges, ax=plt.gca(), orientation='vertical', shrink=0.4, pad=0.02, aspect=20)
    cbar_edges.set_label('Edges: log2(FC) cHi-C', labelpad=10)

    plt.text(0.9, 0.8, 'Jaccard Index between clusters {} {}:\n    Full community: {:.5f}\n    Only baits: {:.5f}'.format(pPartitionString.split("_")[0], pPartitionString.split("_")[1], pJaccardSimilarityCommunities, pJaccardSimilarityBaitCommunities), transform=plt.gca().transAxes,
         fontsize=14, verticalalignment='bottom', horizontalalignment='left',
         bbox=dict(facecolor='white', alpha=0.6, edgecolor='black'))

    plt.title("Interaction vs {} change. Bait {}. Partition: {} {} ".format(pNodeColor, pBaitName, pPartitionString.split("_")[0], pPartitionString.split("_")[1]))
    plt.savefig(pOutputPath + '.' + pFileFormat, bbox_inches='tight')
    plt.close()

def createGraphsBasedOnRNA(pMergedGraph, pNumberOfSteps, pOutputPath, pPartitionsIdString, pChipNamesToPlotList, pJaccardSimilarityCommunities, pJaccardSimilarityBaitCommunities, pPlot=True):

    connected_nodes, connected_edges = getConnectedNodes(pMergedGraph, pNumberOfSteps)

    visited_nodes = set()
    result_for_dfferential_dict = {}
    # Collect data for each bait node
    for bait_node in connected_edges:
        ic(bait_node)
        G = nx.Graph()

        node_attributes = pMergedGraph.nodes[bait_node]

        # Add the node to the new graph with the same attributes
        G.add_node(bait_node, **node_attributes)
        
        for i, level_edges in enumerate(connected_edges[bait_node]):
            visited_nodes = set()
            for edge in level_edges:
                # ic(edge)
                node1, node2 = edge
                if node2 not in visited_nodes:
                    
                    node_attributes = pMergedGraph.nodes[node2]

                    # Add the node to the new graph with the same attributes
                    G.add_node(node2, **node_attributes)
                    edge_attributes = pMergedGraph.get_edge_data(*edge)

                    # Add the edge to the new graph with the same attributes
                    G.add_edge(*edge, **edge_attributes)
                    visited_nodes.add(node2)

            for chip_name in pChipNamesToPlotList:
                if pPlot:
                    plotGraph(G, chip_name, "log2FC", G.nodes[bait_node].get('name', ''), 
                          os.path.join(pOutputPath, pPartitionsIdString + '_' + str(bait_node) + "_" + chip_name + "_" + str(i)), 
                          pPartitionsIdString, pFigureSize=(15,10), pSeed=39,
                         pJaccardSimilarityCommunities=pJaccardSimilarityCommunities, pJaccardSimilarityBaitCommunities=pJaccardSimilarityBaitCommunities)
                nx.write_graphml(G, os.path.join(os.path.join(pOutputPath, pPartitionsIdString + '_' + str(bait_node) + "_" + chip_name + "_" + str(i)) + ".graphml"))

    # plot merged graph
    for chip_name in pChipNamesToPlotList:
        if pPlot:
            plotGraph(pMergedGraph, chip_name, "log2FC", 'Full partition.', os.path.join(pOutputPath, pPartitionsIdString  + "_full_partition_" + chip_name),
                  pPartitionsIdString, pFigureSize=(15,10), pSeed=39, 
                  pJaccardSimilarityCommunities=pJaccardSimilarityCommunities, pJaccardSimilarityBaitCommunities=pJaccardSimilarityBaitCommunities)
        nx.write_graphml(pMergedGraph, os.path.join(os.path.join(pOutputPath, pPartitionsIdString + "_full_partition_" + chip_name) + ".graphml"))


def getConnectedNodes(pMergedGraph, numberOfSteps):
    RNA_logFC_nodes = []
    for node in pMergedGraph.nodes:
        if 'RNA log2FC' in pMergedGraph.nodes[node]:
            RNA_logFC_nodes.append(node)
    
    connected_nodes = {}
    connected_edges = {}

    for node in RNA_logFC_nodes:
        connected_nodes[node] = [set() for _ in range(numberOfSteps)]
        connected_edges[node] = [set() for _ in range(numberOfSteps)]

        # Initialize with direct neighbors and edges
        connected_nodes[node][0] = set(pMergedGraph.neighbors(node))
        connected_edges[node][0] = set(pMergedGraph.edges(node))

        visited_nodes = set()
        visited_nodes.add(node)

        for i in range(1, numberOfSteps):
            for connected_node in connected_nodes[node][i - 1]:
                if connected_node not in visited_nodes:
                    visited_nodes.add(connected_node)
                    neighbor_nodes = pMergedGraph.neighbors(connected_node)
                    neighbor_edges = pMergedGraph.edges(connected_node)

                    neighbor_nodes_update = []
                    neighbor_edges_update = []

                    for neighbor_node in neighbor_nodes:
                        if neighbor_node not in visited_nodes:
                            neighbor_nodes_update.append(neighbor_node)
                    for neighbor_edge in neighbor_edges:
                        if neighbor_edge[1] not in visited_nodes:
                            neighbor_edges_update.append(neighbor_edge)

                    connected_nodes[node][i].update(neighbor_nodes_update)
                    connected_edges[node][i].update(neighbor_edges_update)
    return connected_nodes, connected_edges

def writePlotlyGraph(partition_id, output_path, graph, genes_to_highlight, file_name_prefix):

    
    pos = nx.spring_layout(graph)
    edge_x = []
    edge_y = []
    for edge in graph.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')
    node_x = []
    node_y = []
    node_text = []
    for node in graph.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_text.append(str(node))
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))
    node_adjacencies = []
    node_text = []
    for node, adjacencies in enumerate(graph.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append('# of connections: '+str(len(adjacencies[1])))
    node_trace.marker.color = node_adjacencies
    node_trace.text = node_text
    fig = go.Figure(data=[edge_trace, node_trace],
                layout=go.Layout(
                title='Network graph made with Plotly',
                titlefont_size=16,
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                annotations=[ dict(
                    text="Python code: <a href='https://plotly.com/python/'> https://plotly.com/python/</a>",
                    showarrow=False,
                    xref="paper", yref="paper",
                    x=0.005, y=-0.002 ) ],
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                )
    fig.write_html(os.path.join(output_path, file_name_prefix + str(partition_id) + ".html"))

def intersect(pyRangeObject, position):
    chromosome, position_start_end = position.split(":")
    start, end = position_start_end.split("-")
    query_data = {'Chromosome': [chromosome], 'Start': [int(start)], 'End': [int(end)]}
    query_df = pd.DataFrame(query_data)
    query = pr.PyRanges(query_df)

    return pyRangeObject.intersect(query)

def writePlotlyGraphMerged(id_sample_1, id_sample_2, output_path, graph, genes_to_highlight, sample_name):
    nx.write_graphml(graph, os.path.join(output_path, 'Merged_' + str(id_sample_1) + "_" + str(id_sample_2) + ".graphml"))


def main():
    args = parse_arguments()
    # Check if the output folder exists, if not, create it
    if not os.path.exists(args.outputFolder):
        os.makedirs(args.outputFolder)
    # Load the first graph
    graph1 = nx.read_graphml(args.graph1, node_type=int, edge_key_type=int)
    # ic(graph1)
    # ic(graph1.nodes(data=True))
    # Load the second graph
    graph2 = nx.read_graphml(args.graph2)
    # ic(graph2)
    # ic(graph2.nodes(data=True))
    # Load the first communities file
    communities1 = pd.read_csv(args.communities1)
    # Load the second communities file
    communities2 = pd.read_csv(args.communities2)

    # Transform communities into a dictionary where Gene field is the key and Community the value
    gene_name_to_partition_mapping_sample_1 = {row['Gene']: row['Community'] for _, row in communities1.iterrows()}
    gene_name_to_partition_mapping_sample_2 = {row['Gene']: row['Community'] for _, row in communities2.iterrows()}
    

    # Load the first partitions file
    partitions1 = pd.read_csv(args.partitions1)
    # Load the second partitions file
    partitions2 = pd.read_csv(args.partitions2)

    # Transform partitions into a list of lists
    partition_sample_1 = [set(map(int, entry.split(','))) for entry in partitions1['Content']]
    partition_sample_2 = [set(map(int, entry.split(','))) for entry in partitions2['Content']]

    # Load the RNAseq DESeq file
    rnaseq = pd.read_csv(args.rnaseq,  sep='\t', names=['Gene', 'BaseMean', 'Log2FoldChange', 'LFCSE', 'Stat', 'PValue', 'PAdjust'])
    df_deseq_significant = rnaseq[rnaseq['PAdjust'] <= args.pvalue]
    df_deseq_down = df_deseq_significant.sort_values(by='Log2FoldChange', ascending=True).head(args.level)
    df_deseq_up = df_deseq_significant.sort_values(by='Log2FoldChange', ascending=False).head(args.level)
    log.info(f"Loaded rnaseq from {args.rnaseq}")

    promoters = pd.read_csv(args.promoters, sep='\t')

    ma_norm_pyrange = []

    for i, file in enumerate(args.manorm):
        ma_norm_pyrange.append(createPyRangeObject(file, args.manorm_names[i]))


    cluster_of_regulated_genes = getCommunituesFromDESeq2Ranks([df_deseq_up, df_deseq_down], gene_name_to_partition_mapping_sample_1, name_extension='_1')

    # ic(len(cluster_of_regulated_genes))
    # ic(cluster_of_regulated_genes)

    # ic(partition_sample_1)
    # ic(partition_sample_2)
    jaccardSimilarity = jaccardSimiliarityCommunities(partition_sample_1, partition_sample_2, len(promoters))

    # ic(jaccardSimilarity)
    analyseMergedGraphResult = writePlotlyGraphs(cluster_of_regulated_genes, args.outputFolder, 
                                                 graph1, partition_sample_1, graph2, partition_sample_2, jaccardSimilarity, 
                                                 "Md_", "PA2_", ma_norm_pyrange)

    analyseMergedGraphResultHigh, analyseMergedGraphResultLow = analyseMergedGraphResult

    if analyseMergedGraphResultHigh is not None:
        analyseMergedGraphResultHigh.to_csv(os.path.join(args.outputFolder, "results_high.csv"), index=False)
    if analyseMergedGraphResultLow is not None:
        analyseMergedGraphResultLow.to_csv(os.path.join(args.outputFolder, "results_low.csv"), index=False)