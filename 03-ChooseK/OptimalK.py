#!/usr/bin/env python
########################################################################################
#
# OptimalK.py
# Author: Gao Yueyao
# Python 3.6.10
# Requires the following Python packages:
# numpy(=1.18.1), pandas(1.0.3), matplotlib(3.2.1), scikit-learn(0.23.2)
#
########################################################################################
#
# Import depenencies
#
########################################################################################
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt 
from sklearn.metrics import silhouette_samples, silhouette_score, calinski_harabasz_score,davies_bouldin_score
########################################################################################
#
# Description of script
#
########################################################################################
parser = argparse.ArgumentParser(description="""Choose the optimal K using three statistical measurements""")
#########################################################################################
#
# Required arguments
#
##########################################################################################
parser.add_argument('-i','--inputDIR',dest="inputdir",action='store',required=True,help="Input Directory")
parser.add_argument('-d','--outdir',dest="dpgp_outdir",action='store',required=True,help="DPGP output Directory")
parser.add_argument('-emx','--EMX',dest='gene_expression_matrix',action="store",required=True,help="the path of OFF-removed GEM")
args = parser.parse_args()
#########################################################################################

##########################################################################################
#Load the expression file and the location of DPGP clustering results
os.chdir(args.inputdir)
DTW_DPGP_output = '03-ChooseK/ClusterSumPBS/DTW-K{}-DPGP_ResultsSum.csv'
GEMinput_df = pd.read_csv(args.gene_expression_matrix,index_col=[0])
GEM_arr = GEMinput_df.values
print('Input Gene Expression Matrix has {} entries with {} time points'.format(GEM_arr.shape[0],GEM_arr.shape[1]))

# Normalization with standard scaler
gene_expression_matrix = GEM_arr
gene_expression_matrix -= np.vstack(np.nanmean(gene_expression_matrix, axis=1))
gene_expression_matrix /= np.vstack(np.nanstd(gene_expression_matrix, axis=1))

# Calculate the calinski harabasz index, silhouette score, davies bouldin index of the various K values and Visualize the results
Kval_list = []
for dir in os.listdir(args.dpgp_outdir):
    if "output" in dir:
        Kval = int(dir.split('K')[1].split('D')[0])
        Kval_list.append(Kval)
ch_scorelist = []
sil_scorelist = []
db_scorelist = []
for i in sorted(Kval_list):
    df = pd.read_csv(DTW_DPGP_output.format(i),index_col='gene')
    cluster_labels = df['Cluster'].tolist()
    ch_score = calinski_harabasz_score(gene_expression_matrix, cluster_labels)
    ch_scorelist.append(ch_score)
    print("For DTW_n_clusters =", i,
          "The calinski_harabasz_score is :", ch_score)
    silhouette_avg = silhouette_score(gene_expression_matrix, cluster_labels)
    sil_scorelist.append(silhouette_avg)
    print("For DTW_n_clusters =", i,
          "The average silhouette_score is :", silhouette_avg)
    db_index = davies_bouldin_score(gene_expression_matrix, cluster_labels)
    db_scorelist.append(db_index)
    print("For DTW_n_clusters =", i,
          "The Davies-Bouldin index is :", db_index)
    
    
fig, (ax0,ax1,ax2) = plt.subplots(3,figsize=(7, 21))
ax0.plot(sorted(Kval_list), db_scorelist, 'bo-')
ax0.set_xlabel('DTW_n_clusters')
ax0.set_ylabel('db_score')
ax0.set_title('Davies-Bouldin Index')
ax0.grid(False)

ax1.plot(sorted(Kval_list), ch_scorelist, 'g.-')
ax1.set_xlabel('DTW_n_clusters')
ax1.set_ylabel('ch_score')
ax1.set_title('Calinski and Harabaz Index')
ax1.grid(False)

ax2.plot(sorted(Kval_list), sil_scorelist, 'r*-')
ax2.set_xlabel('DTW_n_clusters')
ax2.set_ylabel('Averge Silhouette Score')
ax2.set_title('Averge Silhouette Coefficient')
ax2.grid(False)

plt.tight_layout()
plt.savefig('03-ChooseK/DTW-DPGP_diffKs_PERFs.png',dpi=200)



        

