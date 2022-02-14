#!/usr/bin/env python
###############################################################################################################
#
# Clusters_plot.py
# Author: Gao Yueyao
# Python 3.6.10
# Requires the following Python packages:
# numpy(=1.18.1), pandas(1.0.3), matplotlib(3.2.1), sklearn(=0.41)
#
###############################################################################################################
#
# Import dependencies
#
###############################################################################################################
import os
import math
import numpy as np
import pandas as pd
import collections
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
from sklearn.neighbors.nearest_centroid import NearestCentroid
import argparse
###############################################################################################################
#
# Description of script
#
###############################################################################################################

parser = argparse.ArgumentParser(description="""
Clusters_plot.py takes a gene expression matrix (GEM) and a clustering results file.
This script plot gene expression over a time course with a panel for each cluster. 
Each panel contains black lines for the expression of each individual gene(or replicate) 
within the cluster and a red line represent the cluster centroid""")

###############################################################################################################
#
# Required arguments
#
###############################################################################################################
parser.add_argument('-exm','--gene_expression_matrix',dest='gene_expression_matrix',action="store",required=True,help="""
The format of the gene_expression_matrix.csv is:
gene    1     2    3    ...    time_t
gene_1  10    20   5    ...    8
gene_2  3     2    50   ...    8
gene_3  18    100  10   ...    22
...
gene_n  45    22   15   ...    60
""")
parser.add_argument('-c','--cluster_label',dest='cluster_label',action='store',required=True,help="""
The format of Clustering_results.csv is:
gene,cluster
gene1,cluster1
gene2,cluster3
gene3,cluster5
...
gene_n,clusterN
The order of index have to be the same as gene_expression_matrix.txt.
""")
parser.add_argument('-o','--output',dest="output_path_prefix",action='store')
args = parser.parse_args()
###############################################################################################################
#
# Import GEM and cluster label file
#
###############################################################################################################

# Reading Input GEM and apply standard scaler
GEMinput_df = pd.read_csv(args.gene_expression_matrix,index_col=[0])
GEM_arr = GEMinput_df.values
print('Input Gene Expression Matrix has {} entries with {} time points'.format(GEM_arr.shape[0],GEM_arr.shape[1]))
# Normalization with standard scaler
gene_expression_matrix = GEM_arr
gene_expression_matrix -= np.vstack(np.nanmean(gene_expression_matrix, axis=1))
gene_expression_matrix /= np.vstack(np.nanstd(gene_expression_matrix, axis=1))

cluster_df = pd.read_csv(args.cluster_label,index_col='gene')
print('Reading clustering results')
ClstNum = len(set(cluster_df['Cluster'].tolist()))
print('There are {} clusters in total'.format(ClstNum))

###############################################################################################################
#
# Plot
#
###############################################################################################################
# Time Points for Plotting
t_label = list(map(int,GEMinput_df.columns.tolist()))
t = t_label
t /= np.mean(np.diff(t))

# Compile the clustering results and compute the Nearest Centroid 
clusters = cluster_df['Cluster'].tolist()
cluster_list = set(cluster_df['Cluster'].tolist())
print('Computing the nearest centroid for each cluster')
X = gene_expression_matrix
y = np.array(clusters)
clf = NearestCentroid()
clf.fit(X, y)

# The gene expression trajectory plots will be generated in several subplots
IDs = list(cluster_list)
# one panel per cluster:
total_subplots = len(IDs)
# max of 6 panels per figure or page
subplots_per_fig = 6
total_no_of_figs = int(np.ceil(total_subplots/float(subplots_per_fig)))
total_cols = 2 # generate this many columns of subplots in the figure.
total_rows = np.ceil(subplots_per_fig/total_cols) # each figure generate will have this many rows.
print('total_no_of_figs',total_no_of_figs)
print('total_rows',total_rows)

IDs_split = [IDs[i:i+subplots_per_fig] for i in range(0, len(IDs), subplots_per_fig)]

clf_Dict = {}
for num,cluster in enumerate(sorted(cluster_list)):
    clf_Dict[cluster] = clf.centroids_[num]
clf_Dict

for c, IDs in enumerate(IDs_split):
    fig = plt.figure(num=None, figsize=(8,12), dpi=300, facecolor='w', edgecolor='k') #figsize=(12,8),
    for i, ID in enumerate(IDs):
        ax = fig.add_subplot(total_rows, total_cols, i+1)
        genes = cluster_df[cluster_df['Cluster']==ID].index
        for gene in genes:
            gene_exp = np.array(GEMinput_df.loc[gene])
            ax.plot(t,gene_exp.ravel(),'k-',alpha=0.15)
        ax.plot(t,clf_Dict[ID], "r-")
        ax.set_xticks(t)
        ax.set_xticklabels(t_label)
        ax.set_ylabel('Gene expression')
        ax.set_title('{} ({})'.format(ID,len(cluster_df[cluster_df['Cluster']==ID])))
    plt.tight_layout()
    plt.savefig('{}_gene_exp_fig_{}.png'.format(args.output_path_prefix,c+1),dpi=600)



