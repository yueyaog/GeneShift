#!~/.conda/envs/my_Python_3.7.3/bin/python3.7

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import seaborn as sns 
from itertools import repeat
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

GEMinput_df = pd.read_csv(args.gene_expression_matrix,index_col=0)
GEM_arr = GEMinput_df.values
print('Input Gene Expression Matrix has {} entries with {} time points'.format(GEM_arr.shape[0],GEM_arr.shape[1]))

# Import Cluster summary and cluster label files
GeneCluster_shift_df = pd.read_csv(args.cluster_label,index_col='gene')
Cond_A = GeneCluster_shift_df['Cluster'][0].split('_DTW')[0]
Cond_B = GeneCluster_shift_df['Cluster'][0].split('_DTW')[1].split("_")[1]
Cluster_summary_df = pd.read_csv(args.cluster_sum,index_col='Cluster')
#allgenes = GEMinput_df.index.tolist()
GeneCluster_shift_df['Cluster']
#Plotting without normalization
ALLGENES = GEMinput_df.index.tolist()

# LINE PLOT
#cluster = 'A_0_B_Mtr-Shoot-DTW-K50-028-DPGP001'
def LnPt_Plot(cluster):
    cluster_genes = GeneCluster_shift_df[GeneCluster_shift_df['Cluster']==cluster].index.tolist()
    A_clustergenes = [i for i in ALLGENES if i.split('_')[1] in cluster_genes and i.split('_')[0]==Cond_A]
    B_clustergenes = [i for i in ALLGENES if i.split('_')[1] in cluster_genes and i.split('_')[0]==Cond_B]

    t_label = list(map(int,GEMinput_df.columns.tolist()))
    t = t_label
    t /= np.mean(np.diff(t))

    A_cluster = cluster.split('{}_'.format(Cond_A))[1].split('_{}'.format(Cond_B))[0]
    B_cluster = cluster.split('{}_'.format(Cond_B))[1]
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
            ax1.set_title('{}: {}'.format(Cond_A, A_cluster))

        
    for gene in B_clustergenes:
        gene_exp = np.array(GEMinput_df.loc[gene])
        ax2.plot(t,gene_exp.ravel(),'g-',alpha=0.2,label='B')
        ax2.set_xticks(t)
        ax2.set_xticklabels(t_label)
        ax2.set_xlabel('Time point')
        ax2.set_title('{}: {}'.format(Cond_B, B_cluster))

    #plt.tight_layout()
    fig.suptitle('{} ({})'.format(title,len(cluster_genes)))
    plt.savefig('LnPt_FPKM_TrajecotrySet_{}.png'.format(cluster),dpi=600)


#BOXPOINT PLOT
# Parse the input GEM into BoxPlot_df
All_Exp_list = []
Time_list = []
for i in range(len(GEMinput_df.columns)):
    Exp_list = GEMinput_df.iloc[:,i].tolist()
    Time_list.extend(repeat(int(GEMinput_df.columns[i]),len(Exp_list)))
    All_Exp_list.extend(Exp_list)
    
Gene_index = [i.split("_")[1] for i in GEMinput_df.index.tolist()]*len(GEMinput_df.columns)
Condition_list = [i.split("_")[0] for i in GEMinput_df.index.tolist()]*len(GEMinput_df.columns)
Rep_list = [i.split("_")[2] for i in GEMinput_df.index.tolist()]*len(GEMinput_df.columns)

BoxPlot_df = pd.DataFrame()
BoxPlot_df['Time'] =  Time_list
BoxPlot_df['Exp'] =  All_Exp_list
BoxPlot_df['Treatment'] =  Condition_list
BoxPlot_df['Rep'] =  Rep_list
BoxPlot_df.index =  Gene_index

# Plot the trajectory set into box point plot
def BxPt_Plot(cluster):
    genelist = GeneCluster_shift_df[GeneCluster_shift_df['Cluster']==cluster].index.tolist()

    fig, (ax0,ax1) = plt.subplots(nrows=1,ncols=2,figsize=(18,6))
    t_label = list(map(int,GEMinput_df.columns.tolist()))
    sns.boxplot(x="Time", y="Exp", data=BoxPlot_df[BoxPlot_df['Treatment']==Cond_A].loc[genelist],order= t_label,ax=ax0)
    sns.boxplot(x="Time", y="Exp", data=BoxPlot_df[BoxPlot_df['Treatment']==Cond_B].loc[genelist],order= t_label,ax=ax1)

    medians_A = BoxPlot_df[BoxPlot_df['Treatment']==Cond_A].loc[genelist].groupby('Time')['Exp'].median().reset_index()
    medians_B = BoxPlot_df[BoxPlot_df['Treatment']==Cond_B].loc[genelist].groupby('Time')['Exp'].median().reset_index()
    sns.pointplot(x="Time", y="Exp", data=medians_A, linestyles='--', color='k', ax=ax0)
    sns.pointplot(x="Time", y="Exp", data=medians_B, linestyles='--', color='k', ax=ax1)

    fig.suptitle('TrajectorySet_{}_Expression'.format(cluster))
    ax0.set_ylabel('log2(FPKM+1)')
    ax1.set_ylabel('log2(FPKM+1)')
    ax0.set_title(Cond_A)
    ax1.set_title(Cond_B)
    plt.savefig('BxPt_FPKM_TrajecotrySet_{}.png'.format(cluster),dpi=600)


# Create Figures of each cluster shift

#Warning Control
plt.rcParams.update({'figure.max_open_warning': 0})

os.chdir(args.output_path)
os.mkdir('GroupA')
os.chdir('GroupA')
GroupA_clusters = Cluster_summary_df[Cluster_summary_df['Category']=='GroupA'].index.tolist()
for cluster in GroupA_clusters:
    LnPt_Plot(cluster)
    BxPt_Plot(cluster)
os.chdir(args.output_path)
os.mkdir('GroupB')
os.chdir('GroupB')
GroupB_clusters = Cluster_summary_df[Cluster_summary_df['Category']=='GroupB'].index.tolist()
for cluster in GroupB_clusters:
    LnPt_Plot(cluster)
    BxPt_Plot(cluster)
os.chdir(args.output_path)
os.mkdir('GroupC')
os.chdir('GroupC')
GroupC_clusters = Cluster_summary_df[Cluster_summary_df['Category']=='GroupC'].index.tolist()
for cluster in GroupC_clusters:
    LnPt_Plot(cluster)
    BxPt_Plot(cluster)


