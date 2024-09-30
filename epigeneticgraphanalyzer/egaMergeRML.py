import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import os
import pandas as pd

import icecream as ic
import pyBigWig

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
        '-o', '--outputFolder', 
        type=str, 
        required=True, 
        help='Path to the output folder'
    )
    return parser.parse_args()

def annotate(pGraph, pConditions):
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
        bw.close()

def main():
    args = parse_arguments()
    if not os.path.exists(args.outputFolder):
        os.makedirs(args.outputFolder)
    # Load the first graph
    graph = nx.read_graphml(args.graph, node_type=int, edge_key_type=int)