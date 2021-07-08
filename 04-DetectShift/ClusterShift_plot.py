#!~/.conda/envs/my_Python_3.7.3/bin/python3.7

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import argparse 

parser = argparse.ArgumentParser(description="""
ClusterShift_plot.py takes a gene expression matrix (GEM) and two clustering results files.
This script plot clusters of gene expression over a time course """)

parser.add_argument('-exm','--gene_expression_matrix',dest='gene_expression_matrix',action="store",required=True)
parser.add_argument('-c','--cluster_label',dest='cluster_label',action='store',required=True,help="[prefix]DTW_DPGP_ClusterShift.csv")
parser.add_argument('-s','--cluster_sum',dest='cluster_sum',action='store',required=True,help="[prefix]ClusterShift_Summary.csv")
parser.add_argument('-o','--output',dest="output_path",action='store')
args = parser.parse_args()

# Reading Input GEM 
#GEMinput_df = pd.read_csv(args.gene_expression_matrix,index_col='Unnamed: 0')
GEMinput_df = pd.read_csv(args.gene_expression_matrix,index_col=0)
GEM_arr = GEMinput_df.values
print('Input Gene Expression Matrix has {} entries with {} time points'.format(GEM_arr.shape[0],GEM_arr.shape[1]))

# Import Cluster summary and cluster label files
GeneCluster_shift_df = pd.read_csv(args.cluster_label,index_col='gene')
Cluster_summary_df = pd.read_csv(args.cluster_sum,index_col='Cluster')
#allgenes = GEMinput_df.index.tolist()

#Plotting without normalization
ALLGENES = GEMinput_df.index.tolist()
#cluster = 'A_0_B_Mtr-Shoot-DTW-K50-028-DPGP001'
def Cluster_shift_FPKMplot(cluster):
    cluster_genes = GeneCluster_shift_df[GeneCluster_shift_df['Cluster']==cluster].index.tolist()
    A_clustergenes = [i for i in ALLGENES if i.split('_')[1] in cluster_genes and i.split('_')[0]=='A']
    B_clustergenes = [i for i in ALLGENES if i.split('_')[1] in cluster_genes and i.split('_')[0]=='B']

    t_label = list(map(int,GEMinput_df.columns.tolist()))
    t = t_label
    t /= np.mean(np.diff(t))

    A_cluster = cluster.split('A_')[1].split('_B')[0]
    B_cluster = cluster.split('B_')[1]
    if A_cluster == B_cluster:
        title = 'NOT SHIFTED'
    elif A_cluster != B_cluster:
        title ='SHIFTED'
    #fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 6),constrained_layout=True, sharey=True)
    for gene in A_clustergenes:
            gene_exp = np.array(GEMinput_df.loc[gene])
            ax1.plot(t,gene_exp.ravel(),'b-',alpha=0.2,label='A')
            ax1.set_xticks(t)
            ax1.set_xticklabels(t_label)
            ax1.set_ylabel('log2(FPKM+1)')
            ax1.set_xlabel('Time point')
            ax1.set_title('A: {}'.format(A_cluster))

        
    for gene in B_clustergenes:
        gene_exp = np.array(GEMinput_df.loc[gene])
        ax2.plot(t,gene_exp.ravel(),'g-',alpha=0.2,label='B')
        ax2.set_xticks(t)
        ax2.set_xticklabels(t_label)
        ax2.set_xlabel('Time point')
        ax2.set_title('B: {}'.format(B_cluster))

    #plt.tight_layout()
    fig.suptitle('{} ({})'.format(title,len(cluster_genes)))
    plt.savefig('FPKM_cluster_shift_{}.png'.format(cluster),dpi=600)
#Cluster_shift_FPKMplot(cluster)


# Create Figures of each cluster shift

#Warning Control
plt.rcParams.update({'figure.max_open_warning': 0})

os.chdir(args.output_path)
os.mkdir('GroupA')
os.chdir('GroupA')
GroupA_clusters = Cluster_summary_df[Cluster_summary_df['Category']=='GroupA'].index.tolist()
for cluster in GroupA_clusters:
    Cluster_shift_FPKMplot(cluster)
os.chdir(args.output_path)
os.mkdir('GroupB')
os.chdir('GroupB')
GroupB_clusters = Cluster_summary_df[Cluster_summary_df['Category']=='GroupB'].index.tolist()
for cluster in GroupB_clusters:
    Cluster_shift_FPKMplot(cluster)
os.chdir(args.output_path)
os.mkdir('GroupC')
os.chdir('GroupC')
GroupC_clusters = Cluster_summary_df[Cluster_summary_df['Category']=='GroupC'].index.tolist()
for cluster in GroupC_clusters:
    Cluster_shift_FPKMplot(cluster)

# #Plotting with Normalization
# #cluster = 'A_Mtr-Shoot-DTW-K35-007-DPGP001_B_Mtr-Shoot-DTW-K35-007-DPGP001'
# def Cluster_shift_plot(cluster):
#     cluster_genes = GeneCluster_shift_df[GeneCluster_shift_df['Cluster']==cluster].index.tolist()
#     A_clustergenes = [i for i in allgenes if i.split('_')[1] in cluster_genes and i.split('_')[0]=='A']
#     B_clustergenes = [i for i in allgenes if i.split('_')[1] in cluster_genes and i.split('_')[0]=='B']

#     t_label = list(map(int,GEMinput_df.columns.tolist()))
#     t = t_label
#     t /= np.mean(np.diff(t))

#     A_cluster = cluster.split('C_')[1].split('_B')[0]
#     B_cluster = cluster.split('B_')[1]
#     if A_cluster == B_cluster:
#         title = 'NOT SHIFTED'
#     elif A_cluster != B_cluster:
#         title ='SHIFTED'

#     #fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))
#     fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 6),constrained_layout=True, sharey=True)

#     for gene in A_clustergenes:
#         gene_exp = np.array(GEMinput_df.loc[gene])
#         ax1.plot(t,gene_exp.ravel(),'b-',alpha=0.15)
#         ax1.set_xticks(t)
#         ax1.set_xticklabels(t_label)
#         ax1.set_ylabel('Gene expression')
#         ax1.set_title('A: {}'.format(A_cluster))

    
#     for gene in B_clustergenes:
#         gene_exp = np.array(GEMinput_df.loc[gene])
#         ax2.plot(t,gene_exp.ravel(),'r-',alpha=0.15)
#         ax2.set_xticks(t)
#         ax2.set_xticklabels(t_label)
#         ax2.set_ylabel('Gene expression')
#         ax2.set_title('B: {}'.format(B_cluster))

#     fig.suptitle('{} ({})'.format(title,len(cluster_genes)))  
#     #plt.title('{} ({})'.format(cluster,len(cluster_genes)),loc='center')
#     plt.tight_layout()
#     plt.savefig('Normed_cluster_shift_{}.png'.format(cluster),dpi=100)


# os.chdir('/scratch2/yueyaog/TimeCourse/DP_GP_cluster/DTW-8LOOPs/DTW_K35-DPGP/GeneShift_Plots/Normed')
# GroupC_clusters = Cluster_summary_df[Cluster_summary_df['Category']=='GroupC'].index.tolist()
# for cluster in GroupC_clusters:
#     Cluster_shift_plot(cluster)

