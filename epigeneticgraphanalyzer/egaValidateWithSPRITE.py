import numpy as np
import scipy
import argparse
import gzip
import pandas as pd

from hicexplorer.hicBuildMatrix import intervalListToIntervalTree, get_bins
from scipy.sparse import lil_matrix

def parse_arguments():
    parser = argparse.ArgumentParser(description='Validate with SPRITE')
    parser.add_argument('--sprite', required=True, type=str, help='Path to the SPRITE input file (must be a .gz file)')
    parser.add_argument('--graph', required=True, help='Path to the graph input file')
    parser.add_argument('--genomeSize', required=True, type=str, help='Path to the genome size file')
    parser.add_argument('--binSize', required=False, type=int, default=10000, help='Size of the bins for genomic intervals')
    return parser.parse_args()

def createSprite3DMap(pSpriteFile, pBinIntervalTree):
    '''
    Open the SPRITE file
    '''

    number_of_bins = 0
    for chromosome in pBinIntervalTree:
        number_of_bins += len(pBinIntervalTree[chromosome])
    
    sprite_3d_map = [lil_matrix((number_of_bins, number_of_bins)) for _ in range(number_of_bins)]
    with gzip.open(pSpriteFile, 'rt') as f:
        for line in f:
            line_split = line.strip().strip().split('\t')
            cluster = line_split[0]
            data = sorted(line_split[1:])
            if len(data) == 1:
                continue
            if len(data) == 2:
                chrom1, start1 = data[0].split(':')
                start1 = int(start1)
                chrom2, start2 = data[1].split(':')
                start2 = int(start2)
                bin1 = sorted(pBinIntervalTree[chrom1][start1]).data
                bin2 = sorted(pBinIntervalTree[chrom2][start2]).data
             
                sprite_3d_map[bin1][bin2] = 1
                sprite_3d_map[bin2][bin1] = 1

            for i in range(len(data)):
                chrom1, start1 = data[0].split(':')
                start1 = int(start1)
                if chrom1 not in pBinIntervalTree:
                    continue
                bin1 = sorted(pBinIntervalTree[chrom1][start1]).data
                for j in range(i+1, len(data)):
                    chrom2, start2 = data[j].split(':')
                    start2 = int(start2)
                    if chrom2 not in pBinIntervalTree:
                        continue
                    bin2 = sorted(pBinIntervalTree[chrom2][start2]).data
                    for k in range(j+1, len(data)):
                        chrom3, start3 = data[k].split(':')
                        start3 = int(start3)
                        if chrom3 not in pBinIntervalTree:
                            continue                                                          
                        bin3 = sorted(pBinIntervalTree[chrom3][start3]).data
                        sprite_3d_map[bin1][bin2, bin3] += 1
                        sprite_3d_map[bin1][bin3, bin2] += 1
                        sprite_3d_map[bin2][bin1, bin3] += 1
                        sprite_3d_map[bin2][bin3, bin1] += 1
                        sprite_3d_map[bin3][bin1, bin2] += 1
                        sprite_3d_map[bin3][bin2, bin1] += 1
    return sprite_3d_map


def main():
    args = parse_arguments()
    graph = nx.read_graphml(args.graph, node_type=int, edge_key_type=int)

    # Read the genome size file with pandas
    genome_size_df = pd.read_csv(args.genomeSize, sep='\t', header=None, names=['chromosome', 'size'])
    bin_intervals = get_bins(args.binSize, genome_size_df.set_index('chromosome')['size'].to_dict())
    bin_interval_tree = intervalListToIntervalTree(bin_intervals)
    three_dimensional_map = createSprite3DMap(args.sprite, bin_interval_tree)


    args.graph
