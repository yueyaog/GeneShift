# GeneShift

This workflow is designed to detect pattern change of time-series gene expression data. 
![Image of GeneShift](https://github.com/yueyaog/GeneShift/blob/master/Auxiliary/GeneShift_repo.png)

## Workflow Summary
 GeneShift workflow performs the following tasks:
  1. Generate input dataset
  2. Compute initial clustering using [soft-DTW-KMeans](https://arxiv.org/abs/1703.01541)
  3. Compute fine clustering with [DP_GP_cluster](https://github.com/PrincetonUniversity/DP_GP_cluster/tree/master/DP_GP)
  4. Determine the optial number of clusters
  5. Using deep learning model(RNN LSTM) to test the robusticity of clustering
  6. Post-clustering analysis (replicate sorting) 
  
## Installation
All of GeneShift dependencies can be installed through Anaconda3. The following installation commands are specified to create Anaconda environments using Clemson's Palmetto cluster. But the dependencies can be easily install to any other HPC system.

Two anaconda environments need to created for GeneShift
```
module load anaconda3/5.1.0-gcc/8.3.1

conda create -n GeneShift_env python=3.6 matplotlib numpy pandas tslearn scikit-learn 

conda create -n DPGP_env python=2.7 GPy pandas numpy scipy matplotlib
```
Once the two anaconda environments have been created, simply clone GeneShift repository to use GeneShift
```
git clone https://github.com/yueyaog/GeneShift.git
```

## Input Data
GeneShift takes ```*_exp.csv``` of the format:
|               | time_1 | time_2 | time_3 | ... | time_n |
|---------------|--------|--------|--------|-----|--------|
| A_gene_1_rep1 |  6.128 |  3.564 |  1.232 | ... |  4.217 |
| A_gene_1_rep2 |  5.940 |  2.817 |  0.715 | ... |  3.829 |
| A_gene_1_rep3 |  6.591 |  3.902 |  1.594 | ... |  4.336 |
| B_gene_1_rep1 |  5.412 |  2.781 |  0.790 | ... |  8.772 |
| B_gene_1_rep2 |  6.195 |  2.066 |  0.815 | ... |  8.891 |
| B_gene_1_rep3 |  5.836 |  3.097 |  0.836 | ... |  9.096 |
| A_gene_2_rep1 |  0.734 |  1.236 |  4.849 | ... |  6.110 |
|      ...      |   ...  |   ...  |   ...  | ... |   ...  |
| B_gene_n_rep3 |  7.889 |  13.206|  11.192| ... |  9.761 |

where the first row is a header containing the time points and the first column is an index containing condition, geneID, and rep. Entries are seperated by comma. 
## Execute the Workflow
The workflow contains a test gene expression matrix for testing. First, user need to execute ```initiate.sh``` to make sure the output from each step will be stored properly.
```
$ ./initiate.sh
```

### Prepare Input
To avoid noises in gene expression data clustering, input data will be seperated into ```OFF_exp.csv``` and ```OFFremoved_exp.csv```. Clustering will only be performed on ```OFFremoved_exp.csv``` . ```OFF_exp.csv``` will be analyzed in post-clustering analysis.
```
$ cd 00-DataPrep/
$ ./00-DataPrep.sh
```

### Initial Clustering (DTW-KMeans)
[Soft-DTW-KMeans](https://arxiv.org/abs/1703.01541) with a range of K values will be appied to the ```OFFremoved_exp.csv```. 
```
$ cd 01-DTWKMeans/
$ ./01-DTWKMeans.sh
```

### Fine Clustering (DP_GP_Cluster)
The initial clustering results will be fine clustered by [Dirichlet process Gaussian process model](https://github.com/PrincetonUniversity/DP_GP_cluster/tree/master/DP_GP).
```
$ cd 02-DP_GP/
$ qsub dpgp_prep.pbs
$ ./02-DP
```

### Choose Optimal K (ch index, db index, silhouette coefficient)
Before determine the optimal K from a range of k values we tested, the clustering results (from 01-DTWKMeans and 02-DP-GP) need to be compiled.
```
$ cd 03-ChooseK
$ ./03-1-ClusterSum.sh
```
Once the clustering results are being summarized into several csv files, three analysis will be used to choose an optimal K value. [Calinski harabasz index](https://doi.org/10.1080/03610927408827101), [silhouette score](https://doi.org/10.1016/0377-0427(87)90125-7), [davies bouldin index](https://doi.org/10.1109/TPAMI.1979.4766909) will be calculated of various K values. The performance of different K values will be visualized in the output plot. 
- Silhouette score is bounded between -1 for incorrect clustering and +1 for highly dense clustering. Scores around zero indicate overlapping clusters.
- Calinski harabasz index is higher when clusters are dense and well separated, which relates to a standard concept of a cluster.
- Davies bouldin index closer to zero indicate a better partition.

### Post-clustering analysis (replicate sorting)


## Classification

__Overview:__

The classifcation script sorts clusters using time series data and supervised learning techniques. The script trains 3 models: a 1-D CNN, an MLP, and an LSTM. The LSTM is the main model while the other two are used for comparison. Each model is trained using a 70-30 train test split of the data provided and outputs a confusion matrix and an F1 score to show the accuracy of each model. Each cluster must contain at least two samples or else the confusion matrix will not be alligned due to the train test split. This can be controlled by the threshold parameter in the pbs script if the user wants to exclude clusters that only have a few samples. 

__Dependencies:__

- pandas
- numpy
- scikit-learn
- matplotlib
- seaborn
- tensorflow

__Data:__

The script takes a single csv file of RNA expression changes over time. The first row is the header containing the gene column, the different time steps, and the cluster column. The other rows contain the gene, RNA expressions at each time step, and the cluster it belongs to. The csv file must be in the 'data' directory. Multiple csv files can be in this directory and the script will output a results directory for each file. An example is shown below, this example has one feature, the RNA expression, and five time steps.
```       
gene                      0             12            24            48            72            clusters    
A17C_Medtr2g026050_rep2	  5.221980164   5.157623152   5.294473715   5.29495152    5.37787568    Mtr-Shoot-DTW-K50-045-DPGP001
A17R_Medtr7g059515_rep2	  0             0.403085897   0.091756243   0.140124224   0.125210248   Mtr-Shoot-DTW-K50-044-DPGP002
```

__PBS Script:__

The PBS script is used to submit a job on the palmetto cluster. The script loops through all the csv files in the 'data' and runs the classification on each file. The command line arguments can be changed in the pbs script and are as follows:
- file_name - The file name of the csv file. This is automatically grabbed in the script so it does not need to be changed.
- threshold - The minimum number of samples a cluster must have or else it is removed from the data. The default value is 2 and it can not be lower than 2.
- test_split - The percentage of data to split into a train-test set. The default value is 0.3 which splits the data into 70% train and 30% test. 
- epochs - The number of passes through the dataset to train the model. Default is 100.

__Output:__

The main output of the script is a confusion matrix and an F1 score for each model, these are created in a 'results' directory. The confusion matrix shows how well the trained model classified each sample from the test split. The expected values are on the y-axis and the predicted values are on the x-axis. This matrix is normalized to show a percentage of samples classfied. So if there are five genes in cluster A but the model predicted two of these genes are in cluster B, the cluster A row would could contain a 60 in column A and a 40 in column B. This normalization can be disabled in the confusion matrix function. The F1 scored is written to a text file

__Running:__
1. Verify all the dependencies are installed in an anaconda virual environment and the correct environment is set in the PBS script
2. Adjust command line arguments in the PBS script if desired
   - Additional parameters can be adjusted within the script itsel
4. Correctly format data move it to the 'data' directory
5. On a palmetto login node enter the following command to submit the job
   - ```qsub cluster-sort.pbs``` 
 
