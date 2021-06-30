#%% 
import numpy as np 
import pandas as pd 
import random
# %%
df = pd.read_csv('/scratch2/yueyaog/TimeCourse/TSLEARN_Kmean/Mtr_Shoot_Prep/logquan_DATA/MtrShoot_FPKM_rep.csv',index_col=[0])
df.head()
# %%
Mtr_genelist = set([i.split('_')[1] for i in df.index.tolist()])
# Select 500 random genes from the Mtr_genelist
sample_list = random.sample(Mtr_genelist,500)

# %%
sample_index_list = []
for index,row in df.iterrows():
    for i in sample_list:
        if i in index:
            sample_index_list.append(index)
#%%
GeneRef_Dict = {}
for count, value in enumerate(sample_list, start=1):
    GeneRef_Dict[value] = 'Gene'+str(count)
GeneRef_Dict

#%%
new_index=[]
for i in sample_index_list:
    new_index.append(i.replace(i.split("_")[1],GeneRef_Dict[i.split("_")[1]]))

Cond_index = []
for i in new_index:
    if 'A17C' in i:
        Cond_index.append(i.replace('A17C','A'))
        
    elif 'A17R' in i:
        Cond_index.append(i.replace('A17R','B'))
    
Cond_index
    
#%%
exp_df = df.reindex(sample_index_list)
exp_df.index = Cond_index

# Add of a column containing a numbered version of the index
exp_df['indexNumber'] = [int(i.split('_')[1].split('Gene')[1]) for i in Cond_index]
exp_df['Condition'] = [i.split("_")[0] for i in Cond_index]
exp_df['Rep'] = [int(i.split("_")[2].split("rep")[1]) for i in Cond_index]
# Perform sort of the rows
exp_df.sort_values(by=['indexNumber','Condition','Rep'],inplace=True)
#%%


# Deletion of the added column
exp_df.drop(['indexNumber','Condition','Rep'], axis=1, inplace = True)

#%%
exp_df.to_csv('/zfs/lasernode/feltuslab/gaoyy/GeneShift/Test/Input/Test_exp.csv')
# %%
#exp_df = df.reindex(sample_index_list)
sum_df = pd.DataFrame(exp_df.sum(axis=1))
# The gene expression are 0 at every time point
AllZero_index = sum_df[sum_df[0]==0].index.tolist()
#AllZero_df = exp_df.reindex(AllZero_index)

# Include the replicates of the off genes
off_gene_index = []
for i in AllZero_index:
    off_gene_index.append(i.split('_r')[0])

OFFgenes = []
for i in set(off_gene_index):
    OFFgenes.append(i+'_rep1')
    OFFgenes.append(i+'_rep2')
    OFFgenes.append(i+'_rep3')

#sum_df.reindex(All_OffGenes).sort_values(by=[0],ascending=False).head()
OFFgenes_df = sum_df.reindex(OFFgenes)
OFFgenes_df.sort_values(by=[0],ascending=False,inplace=True)
OFFgenes_df.describe()
#MtrRoot_OFFgenes =exp_df.reindex(OFFgenes)
# %%
off_gene_df = exp_df.reindex(OFFgenes)
off_gene_df.to_csv('/zfs/lasernode/feltuslab/gaoyy/GeneShift/Test/Input/Test_OFF_exp.csv')
# %%
off_rm_df = exp_df.reindex(set(Cond_index)-set(OFFgenes))
# %%
off_rm_df.to_csv('/zfs/lasernode/feltuslab/gaoyy/GeneShift/Test/Input/Test_OFFremoved_exp.csv')

# %%
