#!/usr/bin/env python

# Import dependencies
import numpy as np 
import pandas as pd 
import argparse

parser = argparse.ArgumentParser(description="""
Prepare_inputs.py takes a gene expression matrix file and seperate it
into two files: OFF and OFFremoved. If any replicate in A_geneX is complete zero cross all time points, A_geneX
will be categorized into OFF """)

# Required arguments
parser.add_argument('-i','--input',dest='gene_expression_matrix',action="store",required=True,help="""
The input gene_expression_matrix need to be in csv format.""")
parser.add_argument('-o','--output_prefix',dest="output_file_prefix",action='store')
args = parser.parse_args()

# Import data 
input_df = pd.read_csv(args.gene_expression_matrix,index_col=[0])
Cond_index = input_df.index.tolist()
sum_df = pd.DataFrame(input_df.sum(axis=1))
# The gene expression are 0 at every time point
AllZero_index = sum_df[sum_df[0]==0].index.tolist()

# Include the replicates of the off genes
off_gene_index = []
for i in AllZero_index:
    off_gene_index.append(i.split('_r')[0])

OFFgenes = []
for i in set(off_gene_index):
    OFFgenes.append(i+'_rep1')
    OFFgenes.append(i+'_rep2')
    OFFgenes.append(i+'_rep3')

off_gene_df = input_df.reindex(OFFgenes)
off_gene_df.to_csv(args.output_file_prefix+'_OFF_exp.csv')

off_rm_df = input_df.reindex(set(Cond_index)-set(OFFgenes))
off_rm_df.to_csv(args.output_file_prefix+'_OFFremoved_exp.csv')



