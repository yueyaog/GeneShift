#!/usr/bin/env python
########################################################################################
#
# Cluster_Compile.py
# Author: Gao Yueyao
# Python 3.6.10
# Requires the following Python packages:
# numpy(=1.18.1), pandas(1.0.3)
#
########################################################################################
#
# Import depenencies
#
########################################################################################
import pandas as pd
import numpy as np
import os
import argparse
########################################################################################
#
# Description of script
#
########################################################################################
parser = argparse.ArgumentParser(description="""Summarize the Cluster results from DTW-DPGP. The script use a argv parser \
instead of a for loop approach to save computational time. The pbs jobs will be run in parallel""")
#########################################################################################
#
# Required arguments
#
##########################################################################################
parser.add_argument('-i','--inputDIR',dest="inputdir",action='store',required=True,help="Input Directory")
parser.add_argument('-emx','--EMX',dest='gene_expression_matrix',action="store",required=True,help="the path of OFF-removed GEM")
parser.add_argument('-k','--Kval',dest="Kval",type=int,action='store',help="The DTW K value")
args = parser.parse_args()
#########################################################################################

##########################################################################################
os.chdir(args.inputdir)
DPGP_output = '02-DP_GP/DTW-K{}DPGP_output/DTW-Cluster{}_optimal_clustering.txt'
Exp_df = pd.read_csv(args.gene_expression_matrix,index_col=[0])

DTWDPGP_Results_df = pd.DataFrame(Exp_df.index.tolist())
DTWDPGP_Results_df.set_index(0,inplace=True)
DTWDPGP_Results_df.rename_axis('gene',inplace=True)
DTWDPGP_Results_df.head()

ClusterDict={}
for i in range(1,args.Kval+1):
    df = pd.read_csv(DPGP_output.format(args.Kval,i),sep='\t',index_col='cluster')
    print('DTWK{}_Cluster{} DPGP fine clustering output exists.\nReading DPGP Fine Clustering Results.'.format(args.Kval,i))
    for index,row in df.iterrows():
        key = row['gene']
        val = 'DTW-K{}-{}-DPGP{}'.format(str(args.Kval),str(i).rjust(3,'0'),str(index).rjust(3,'0'))
        ClusterDict[key] = val
        DTWDPGP_Results_df['Cluster'] = DTWDPGP_Results_df.index.map(ClusterDict.get)
        DTWDPGP_Results_df.to_csv('03-ChooseK/ClusterSummary/DTW-K{}-DPGP_ResultsSum.csv'.format(str(args.Kval)))