#!~/.conda/envs/my_Python_3.7.3/bin/python3.7

import pandas as pd 
import numpy as np 
import os
import argparse

parser = argparse.ArgumentParser(description="""
Cluster_change.py is the final step of GeneShift workflow""")

parser.add_argument('-off','--OFF',dest='off_emx',action='store',help='required, The OFF gene expression matrix')
parser.add_argument('-c','--cluster_label',dest='cluster_label',action='store',required=True,help="""
The format of Clustering_results.csv is:
gene,cluster
gene1,cluster1
gene2,cluster3
gene3,cluster5
...
gene_n,clusterN
""")
parser.add_argument('-r','--rep',dest="rep",type=int,action='store',help="The replicate threshold")
parser.add_argument('-o','--output',dest="output_path_prefix",action='store',help='required, e.g. /path/to/PERF_results')
args = parser.parse_args()

#Obtain OFF gene in A and B in root
OFF_genes_df = pd.read_csv(args.off_emx,index_col=[0])
AllOFF_genes = OFF_genes_df.index.tolist()
A_OFF_genes = set([i.split('_r')[0] for i in AllOFF_genes if 'A' in i])
B_OFF_genes = set([i.split('_r')[0] for i in AllOFF_genes if 'B' in i])
print('##################### Cluster Change Output by GeneShift#######################')
print('')
print('#A OFF genes: ',len(A_OFF_genes))
print('#B OFF genes: ',len(B_OFF_genes))


""" Import the optimal DTW-DPGP clustering results on input Data 
    Remove the samples are not consistent cross replicates
    Set the shared cluster as final cluster"""

Cluster_df = pd.read_csv(args.cluster_label,index_col='gene')
sample_genelist = list(set([i.split('_r')[0] for i in Cluster_df.index.tolist()]))
 
rep1_Dict={}
rep2_Dict={}
rep3_Dict={}
for index,row in Cluster_df.iterrows():
    for i in sample_genelist:
        if i in index and 'rep1' in index:
            rep1_Dict[i]=row['Cluster']
        elif i in index and 'rep2' in index:
            rep2_Dict[i]=row['Cluster']
        elif i in index and 'rep3' in index:
            rep3_Dict[i]=row['Cluster']


 
smpl_GeneCluster_df = pd.DataFrame()
smpl_GeneCluster_df['gene'] = sample_genelist
smpl_GeneCluster_df['rep1'] = smpl_GeneCluster_df['gene'].map(rep1_Dict)
smpl_GeneCluster_df['rep2'] = smpl_GeneCluster_df['gene'].map(rep2_Dict)
smpl_GeneCluster_df['rep3'] = smpl_GeneCluster_df['gene'].map(rep3_Dict)
 
Consist_Dict = {}
Final_Cluster_Dict = {}
for index,row in smpl_GeneCluster_df.iterrows():
    if row['rep1'] == row['rep2'] ==row['rep3']:
        Consist_Dict[row['gene']] = 3
        Final_Cluster_Dict[row['gene']] = row['rep1']
    elif row['rep1'] == row['rep2'] != row['rep3']:
        Consist_Dict[row['gene']] = 2
        Final_Cluster_Dict[row['gene']] = row['rep1']
    elif row['rep1'] == row['rep3'] != row['rep2']:
        Consist_Dict[row['gene']] = 2
        Final_Cluster_Dict[row['gene']] = row['rep1']
    elif row['rep2'] == row['rep3'] != row['rep1']:
        Consist_Dict[row['gene']] = 2
        Final_Cluster_Dict[row['gene']] = row['rep2']
    else:
        Consist_Dict[row['gene']] = 1


smpl_GeneCluster_df['consist_num'] = smpl_GeneCluster_df['gene'].map(Consist_Dict)
smpl_GeneCluster_df['Final_Cluster'] = smpl_GeneCluster_df['gene'].map(Final_Cluster_Dict)

print('The threshold for consistent rep is ',args.rep)
if args.rep==2:
    Consist_GeneCluster_df = smpl_GeneCluster_df[smpl_GeneCluster_df['consist_num']>1]
elif args.rep==3:
    Consist_GeneCluster_df = smpl_GeneCluster_df[smpl_GeneCluster_df['consist_num']==3]

cond_genes = Consist_GeneCluster_df['gene'].tolist()
allgenes = list(set([i.split('_')[1] for i in Consist_GeneCluster_df['gene'].tolist()]))

first_list = []
second_list = []
third_list = []
for gene in allgenes:
    if 'B_'+gene in cond_genes and 'A_'+gene in cond_genes:
        first_list.append(gene)
    elif 'B_'+gene in cond_genes and 'A_'+gene in A_OFF_genes:
        second_list.append(gene)
    elif 'A_'+gene in cond_genes and 'B_'+gene in B_OFF_genes:
        third_list.append(gene)
qualified_genes = first_list+second_list+third_list
print('There are {} valid(consistent) genes in total\n'.format(len(allgenes)))
print('There are {} genes are qualified for downstream analysis\n'.format(len(qualified_genes)))
print('Pattern to Pattern: ',len(first_list))
print('A OFF to B Pattern: ',len(second_list))
print('A Pattern to B OFF: ',len(third_list))

def overlap_check(list):
   if len(list)>len(set(list)):
       print('{} elements are overlapping'.format(len(list)-len(set(list))))
   elif len(list) == len(set(list)):
       print('No overlapping elements being found')
print('Pattern to Pattern: {}'.format(overlap_check(first_list)))
print('A OFF to B Pattern: {}'.format(overlap_check(second_list)))
print('A Pattern to B OFF: {}'.format(overlap_check(third_list)))


"""For the selected genes, I seperate them into clusters based on A->B
    There are three groups based on the cluster shift"""
A_Dict = {}
B_Dict = {}
for index,row in Consist_GeneCluster_df.iterrows():
    for i in qualified_genes:
        if i in row['gene'] and 'A' in row['gene']:
            A_Dict[i] = row['Final_Cluster']
        elif i in row['gene'] and 'B' in row['gene']:
            B_Dict[i] = row['Final_Cluster']

Cluster_Switch_df = pd.DataFrame()
Cluster_Switch_df['gene'] = qualified_genes
Cluster_Switch_df['A_Cluster'] = Cluster_Switch_df['gene'].map(A_Dict)
Cluster_Switch_df['B_Cluster'] = Cluster_Switch_df['gene'].map(B_Dict)


Cluster_Switch_df = Cluster_Switch_df.fillna(0)
"""Analyze the selected genes based on the gene shift between A and B"""
Cluster_Dict = {}
for index,row in Cluster_Switch_df.iterrows():
    Cluster_Dict[row['gene']] = 'A_'+str(row['A_Cluster'])+'_B_'+str(row['B_Cluster'])
Cluster_Switch_df['Cluster'] = Cluster_Switch_df['gene'].map(Cluster_Dict)

 
Cluster_Summary_df = pd.DataFrame(Cluster_Switch_df.groupby(['Cluster']).count()['gene'])
 

"""Genes that shift between conditions and not shift between conditions"""
Non_shifts_list = []
for index,row in Cluster_Summary_df.iterrows():
    A_cluster = index.split('_')[1]
    B_cluster = index.split('_')[-1]
    if A_cluster == B_cluster: 
        Non_shifts_list.append(index)
print(Cluster_Summary_df[Cluster_Summary_df.index.isin(Non_shifts_list)])
Nonshifts_genes = Cluster_Summary_df[Cluster_Summary_df.index.isin(Non_shifts_list)]['gene'].sum()
Shifts_genes = Cluster_Summary_df['gene'].sum()-Nonshifts_genes

Shift_Dict = {}
for index,row in Cluster_Switch_df.iterrows():
    if row['Cluster'] in Non_shifts_list:
        Shift_Dict[row['gene']] = 0
    else:
        Shift_Dict[row['gene']] = 1

Cluster_Switch_df['Shift'] = Cluster_Switch_df['gene'].map(Shift_Dict)

print('For {} genes under our analysis.\n{} genes shift between A and B.\n{} genes behave the same pattern under A and B.'.format(len(qualified_genes),Shifts_genes,Nonshifts_genes))
print('For {} genes under condition A and B, they have {} ways of responding'.format(len(qualified_genes),len(set(Cluster_Switch_df['Cluster'].tolist()))))


 
Category_Dict = {}
for index,row in Cluster_Summary_df.iterrows():
    if 'A_0' in index:
        Category_Dict[index] = 'GroupA'
    elif 'B_0' in index:
        Category_Dict[index] = 'GroupB'
    else:
        Category_Dict[index] = 'GroupC'

Cluster_Summary_df['Category'] = Cluster_Summary_df.index.map(Category_Dict.get)
Cluster_Switch_df['Category'] = Cluster_Switch_df['Cluster'].map(Category_Dict)
 
#os.chdir('/scratch2/yueyaog/TimeCourse/DP_GP_cluster/Phase2_LOOP/DTW_DPGP_Results/Optimal_K50')

smpl_GeneCluster_df.to_csv('{}GeneCluster_Consist_Filter.csv'.format(args.output_path_prefix))
Cluster_Switch_df.to_csv('{}_DTW_DPGP_ClusterShift.csv'.format(args.output_path_prefix),index=False)
Cluster_Summary_df = Cluster_Summary_df.sort_values(by=['gene'],ascending=False)
Cluster_Summary_df.to_csv('{}ClusterShift_Summary.csv'.format(args.output_path_prefix))
print('The results have saved to csv files')
