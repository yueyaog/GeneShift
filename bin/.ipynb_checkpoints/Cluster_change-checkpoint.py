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

#Obtain OFF gene in conditionA and conditionB
OFF_genes_df = pd.read_csv(args.off_emx,index_col=[0])
AllOFF_genes = OFF_genes_df.index.tolist()
Cond_A = list(set([i.split('_')[0] for i in AllOFF_genes]))[0]
Cond_B = list(set([i.split('_')[0] for i in AllOFF_genes]))[1]
OFF_genes_df = pd.read_csv(args.off_emx,index_col=[0])
AllOFF_genes = OFF_genes_df.index.tolist()

A_OFF_genes = set([i.split('_r')[0] for i in AllOFF_genes if Cond_A in i])
B_OFF_genes = set([i.split('_r')[0] for i in AllOFF_genes if Cond_B in i])
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
    if '{}_'.format(Cond_B)+gene in cond_genes and '{}_'.format(Cond_A)+gene in cond_genes:
        first_list.append(gene)
    elif '{}_'.format(Cond_B)+gene in cond_genes and '{}_'.format(Cond_A)+gene in A_OFF_genes:
        second_list.append(gene)
    elif '{}_'.format(Cond_A)+gene in cond_genes and '{}_'.format(Cond_B)+gene in B_OFF_genes:
        third_list.append(gene)
qualified_genes = first_list+second_list+third_list
print('There are {} valid(consistent) genes in total\n'.format(len(allgenes)))
print('There are {} genes are qualified for downstream analysis\n'.format(len(qualified_genes)))
print('Pattern to Pattern: ',len(first_list))
print('{} Not Expressed to {} Dynamic Pattern: '.format(Cond_A,Cond_B),len(second_list))
print('{} Dynamic Pattern to {} Not Expressed: '.format(Cond_A,Cond_B),len(third_list))

def overlap_check(list):
   if len(list)>len(set(list)):
       print('{} elements are overlapping'.format(len(list)-len(set(list))))
   elif len(list) == len(set(list)):
       print('No overlapping elements being found')
print('Pattern to Pattern: {}'.format(overlap_check(first_list)))
print('{} Not Expressed to {} Dynamic Pattern: {}'.format(Cond_A,Cond_B,overlap_check(second_list)))
print('{} Dynamic Pattern to {} Not Expressed: {}'.format(Cond_A,Cond_B,overlap_check(third_list)))


"""For the selected genes, I seperate them into clusters based on A->B
    There are three groups based on the cluster shift"""
A_Dict = {}
B_Dict = {}
for index,row in Consist_GeneCluster_df.iterrows():
    for i in qualified_genes:
        if i in row['gene'] and Cond_A in row['gene']:
            A_Dict[i] = row['Final_Cluster']
        elif i in row['gene'] and Cond_B in row['gene']:
            B_Dict[i] = row['Final_Cluster']

Cluster_Switch_df = pd.DataFrame()
Cluster_Switch_df['gene'] = qualified_genes
Cluster_Switch_df['{}_Cluster'.format(Cond_A)] = Cluster_Switch_df['gene'].map(A_Dict)
Cluster_Switch_df['{}_Cluster'.format(Cond_B)] = Cluster_Switch_df['gene'].map(B_Dict)


Cluster_Switch_df = Cluster_Switch_df.fillna(0)
"""Analyze the selected genes based on the gene shift between two conditions"""
Cluster_Dict = {}
for index,row in Cluster_Switch_df.iterrows():
    Cluster_Dict[row['gene']] = '{}_'.format(Cond_A)+str(row['{}_Cluster'.format(Cond_A)])+'_{}_'.format(Cond_B)+str(row['{}_Cluster'.format(Cond_B)])
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

print('For {} genes under our analysis.\n{} genes shift between {} and {}.\n{} genes behave the same pattern under {} and {}.'.format(len(qualified_genes),Shifts_genes,Cond_A,Cond_B,Nonshifts_genes,Cond_A,Cond_B))
print('For {} genes under condition {} and {}, they have {} ways of responding'.format(len(qualified_genes),Cond_A,Cond_B,len(set(Cluster_Switch_df['Cluster'].tolist()))))


####### Feb6, EDITSSSSSS############# 
Category_Dict = {}
for index,row in Cluster_Summary_df.iterrows():
    if '{}_0'.format(Cond_A) in index:
        Category_Dict[index] = 'GroupA'
    elif '{}_0'.format(Cond_B) in index:
        Category_Dict[index] = 'GroupB'
    else:
        Category_Dict[index] = 'GroupC'

Cluster_Summary_df['Category'] = Cluster_Summary_df.index.map(Category_Dict.get)
Cluster_Switch_df['Category'] = Cluster_Switch_df['Cluster'].map(Category_Dict)
 
smpl_GeneCluster_df.to_csv('{}GeneCluster_Consist_Filter.csv'.format(args.output_path_prefix))
Cluster_Switch_df.to_csv('{}_DTW_DPGP_ClusterShift.csv'.format(args.output_path_prefix),index=False)
Cluster_Summary_df = Cluster_Summary_df.sort_values(by=['gene'],ascending=False)
Cluster_Summary_df.to_csv('{}ClusterShift_Summary.csv'.format(args.output_path_prefix))
print('The results have saved to csv files')
