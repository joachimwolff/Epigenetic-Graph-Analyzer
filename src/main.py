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
# This is the main script for the Epigenetic Graph Analyzer project.
# The specific functionalities will be added as the project progresses.

# main.py



# provide a function to read bigwig files
def read_bigwig_file(file_path):
    # read the bigwig file
    bw = pyBigWig.open(file_path)
    values = bw.values()
    bw.close()
    return values


def read_interaction_file(file_path):
    interactions = []
    with gzip.open(file_path, 'rt') as file:
        for line in file:
            line = line.strip()
            if line:
                parts = line.split('\t')
                interaction_x = parts.split(',')
                interactions.append(interaction)
                interaction_y = parts.split(',')
                interactions.append(interaction)
    return interactions

def parse_arguments():
    parser = argparse.ArgumentParser(description='Epigenetic Graph Analyzer')
    parser.add_argument('interactions', type=str, nargs='+', help='Path to the text file with interactions')
    parser.add_argument('histone_data', type=str, nargs='+', help='Path(s) to the histone data file(s) (bigwig)')
    parser.add_argument('rna_expression_data', type=str, nargs='+', help='Path(s) to the RNA expression data file(s) (bigwig)')
    parser.add_argument('promoter_region_file', type=str, help='Path to the promoter region file')
    parser.add_argument('enhancer_region_file', type=str, help='Path to the enhancer region file')

    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use (default: 1)')
    parser.add_argument('-v', '--version', action='version', version='Epigenetic Graph Analyzer 1.0')
    return parser.parse_args()



def main():
    args = parse_arguments()

    # TODO: Add your code here
    # read promoter and enhancer regions
    # the promoter and enhancer together form the nodes of the graph aka the entries of a matrix
    # there must be a mapping of position in the matrix to the promoter and enhancer regions
    # read the interaction file and 

if __name__ == '__main__':
    main()
