import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import os
import pandas as pd

import icecream as ic
import pyBigWig

from epigeneticgraphanalyzer.egaMergeGraphs import getConnectedNodes

def parse_arguments():
    parser = argparse.ArgumentParser(description="Epigenetic Graph Analyzer")
    parser.add_argument(
        '-g', '--graph', 
        type=str, 
        required=True, 
        help='Path to the input graph file'
    )

    parser.add_argument(
        '-n', '--conditions', 
        type=str, 
        nargs='+', 
        required=True, 
        help='List of names corresponding to the input bigwig files'
    )
    parser.add_argument(
        '-b', '--bigwigFiles', 
        type=str, 
        nargs='+', 
        required=True, 
        help='List of paths to the input bigwig files'
    )

    parser.add_argument(
        '-d', '--deseqCounts', 
        type=str, 
        required=True, 
        help='Path to the input file with normalized counts from DESeq2'
    )
    parser.add_argument(
        '-o', '--outputFolder', 
        type=str, 
        required=True, 
        help='Path to the output folder'
    )
    return parser.parse_args()

def annotate(pGraph, pConditions, pRNA):
    '''
    Annotate the graph with the bigwig files
    '''


    for name in pConditions:
        for i, replicate in enumerate(pConditions[name]):
            bw = pyBigWig.open(replicate)
            for node in pGraph.nodes:
                position = node['position']
                chromosome, value = position.split(':')
                start, end = value.split('-')
                values = bw.values(chromosome, int(start), int(end))
                node[name + '_R' + str(i)] = values

                if pRNA:
                    
        bw.close()

def createRLMTable(pGraph, pSteps, pConditionNames):
    '''
    Create the table with the RLM values
    '''
    
    for node in pGraph.nodes:
        connected_nodes, connected_edges = getConnectedNodes(pMergeGraph, numberOfSteps)

        results = []

        visited_nodes = set()
        result_for_dfferential_dict = {}
        # Collect data for each bait node
        for bait_node in connected_edges:

            for i, level_edges in enumerate(connected_edges[bait_node]):
                visited_nodes = set()
                interaction_logFC = []
                m_values = {}
                m_values_bait = {}


                for edge in level_edges:
                    node1, node2 = edge
                    connected_node = node2 #if node1 == bait_node else node1

                    if connected_node not in visited_nodes:
                        visited_nodes.add(connected_node)

                        for name in pConditionNames:
                            pass

                            # if /

def main():
    args = parse_arguments()
    if not os.path.exists(args.outputFolder):
        os.makedirs(args.outputFolder)
    # Load the first graph
    graph = nx.read_graphml(args.graph, node_type=int, edge_key_type=int)

    if (len(args.conditions) % len(args.bigwigFiles)) == 0:
        raise ValueError('The number of conditions must be multiple of the number of bigwig files')
    
    conditions = {}
    for i, condition in enumerate(args.conditions):
        conditions[condition] = args.bigwigFiles[i].split(',')
    
    annotated_graph = annotate(graph, conditions)
    
