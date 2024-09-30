import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import os
import pandas as pd

import icecream as ic

def parse_arguments():
    parser = argparse.ArgumentParser(description="Epigenetic Graph Analyzer")
    parser.add_argument('--regulationFile', type=str, required=True, help='Path to the regulation file')
    parser.add_argument('--step', type=int, required=True, help='Step value for the analysis')
    parser.add_argument('--outputFolder', type=str, required=True, help='Path to the output file')
    parser.add_argument('--format', type=str, choices=['pdf', 'svg', 'eps', 'png'], default='pdf', help='Output format for the plots')
    parser.add_argument('--colormap', type=str, default='seismic', help='Colormap to use for the plots')
    parser.add_argument('--figsize', default="10,30", help='Figure size for the plots, format: width,height')
    parser.add_argument('--skipLabels', action='store_false', help='Skip y-axis labels in the plots')
    return parser.parse_args()

def plotClustermap(pData, pColumns, pStepLevel, pOutputPath, pFigsize=(10, 30), pYticklabels=True, colormap='seismic', format='pdf'):

    cmap = plt.get_cmap(colormap)
    pDataStep = pData
    pDataStep.fillna(0, inplace=True)

    pDataStep['Combined_Index'] = pDataStep.cluster_id.astype(str) + '-' + pDataStep['name']

    maximum = np.max(pDataStep[pColumns].values)
    minimum = np.min(pDataStep[pColumns].values)

    norm = mcolors.TwoSlopeNorm(vmin=minimum, vcenter=0, vmax=maximum)
    # Set this new combined index as the index of the DataFrame
    pDataStep.set_index('Combined_Index', inplace=True)
    clusterGrid = sns.clustermap(pDataStep[pColumns], cmap=cmap, norm=norm, figsize=pFigsize, yticklabels=pYticklabels)
    plt.title('Step size {}'.format(pStepLevel + 1))
    plt.savefig(os.path.join(pOutputPath, 'clustermap_step_{}.' + format).format(pStepLevel + 1), bbox_inches='tight')

    return clusterGrid

def identifyRegulation(pDataFrame, pSuffixBait, pSuffixOE):
    '''
    Identify the regulation of the genes in the cluster. Is the gene regulated by bait region histone modification change, OE change, both or none.
    '''

    if pDataFrame.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    
    bait_list = [col for col in pDataFrame.columns if col.endswith(pSuffixBait)]
    OE_list = [col for col in pDataFrame.columns if col.endswith(pSuffixOE)]

    dataFrameBait = pDataFrame[bait_list]
    dataFrameOE = pDataFrame[OE_list]

    dataFrameBait['Sum'] = (dataFrameBait != 0).sum(axis=1)
    dataFrameOE['Sum'] = (dataFrameOE != 0).sum(axis=1)

    baitRegulatedFilter = dataFrameBait['Sum'] != 0
    OERegulatedFilter = dataFrameOE['Sum'] != 0

    baitRegulated = pDataFrame[baitRegulatedFilter]
    OERegulated = pDataFrame[OERegulatedFilter]

    bothRegulated = pDataFrame[baitRegulatedFilter | OERegulatedFilter]

    noneRegulated = pDataFrame[~(baitRegulatedFilter | OERegulatedFilter)]

    print(len(baitRegulated))
    print(len(OERegulated))

    print(len(bothRegulated))
    print(len(noneRegulated))
    print(len(pDataFrame))


    return baitRegulated, OERegulated, bothRegulated, noneRegulated

def main():
    args = parse_arguments()
    # Add your main logic here
    print(args.outputFolder)
    if not os.path.exists(args.outputFolder):
        if args.outputFolder:
            os.makedirs(args.outputFolder, exist_ok=True)
    
    width, height = map(float, args.figsize.split(','))

    dataFrame = pd.read_csv(args.regulationFile)

    dataFrame.fillna(0, inplace=True)
    dataFrame['step'] = dataFrame['step'].astype(int)
    suffixBait = '_bait'
    suffixOE = '_OE'
    negative_list = [col for col in dataFrame.columns if 'mean' in col and col.startswith('negative interaction')]
    positive_list = [col for col in dataFrame.columns if 'mean' in col and col.startswith('positive interaction')]

    rna_list = [col for col in dataFrame.columns if col.endswith('RNA_log2FC')]

    interaction_list = [col for col in dataFrame.columns if ('mean' in col) and col.endswith('logFC')]

    columnsListIdenticication = rna_list + interaction_list + negative_list + positive_list + ['name', 'cluster_id']
    columnsListPlot = rna_list + interaction_list + negative_list + positive_list

    # dataFrameBait = dataFrame[bait_list]
    # dataFrameOE = dataFrame[OE_list]
    noneRegulated = None
    for i in range(0, args.step):
        print( i+1)

        if noneRegulated is None:
            dataFrameLocal = dataFrame[dataFrame['step'] == i]
        else:
            names = noneRegulated['name']
            dataFrameLocal = dataFrame[dataFrame['step'] == i]
            dataFrameLocal = dataFrameLocal[dataFrameLocal['name'].isin(names)]
        if dataFrameLocal is None:
            break
        print( i+1)
        plotClustermap(dataFrameLocal, columnsListPlot, i, args.outputFolder, pFigsize=(width, height), colormap=args.colormap, pYticklabels=args.skipLabels, format=args.format)
        
        # if noneRegulated is None:
        #     break
        baitRegulated, OERegulated, bothRegulated, noneRegulated = identifyRegulation(noneRegulated[columnsListIdenticication], suffixBait, suffixOE)       

        baitRegulated.to_csv(os.path.join(args.outputFolder, f'baitRegulated_step_{i}.csv'), index=False)
        OERegulated.to_csv(os.path.join(args.outputFolder, f'OERegulated_step_{i}.csv'), index=False)
        bothRegulated.to_csv(os.path.join(args.outputFolder, f'bothRegulated_step_{i}.csv'), index=False)
        noneRegulated.to_csv(os.path.join(args.outputFolder, f'noneRegulated_step_{i}.csv'), index=False)
