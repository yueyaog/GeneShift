# GeneShift

This workflow is designed to detect pattern change of time-series gene expression data.   GeneShift workflow will perform following tasks:
  * OFF gene check
  * Cluster expression of individual replicate over time using [DTW-KMeans](https://arxiv.org/abs/1703.01541) and then fine clustering with [DP_GP_cluster](https://github.com/PrincetonUniversity/DP_GP_cluster/tree/master/DP_GP)
  * Using deep learning model(RNN LSTM) to test the robusticity of clustering 
  * Cluster finalization (remove poor quality gene based on replicate consistency) 
  * Conduct pattern and degree change capture and visualize the performance


## Classification

__Overview:__

The classifcation script sorts clusters using time series data and supervised learning techniques. The script trains 3 models: a 1-D CNN, an MLP, and an LSTM. The LSTM is the main model while the other two are used for comparison. Each model is trained using a 70-30 train test split of the data provided and outputs a confusion matrix and an F1 score to show the accuracy of each model. Each cluster must contain at least two samples or else the confusion matrix will not be alligned due to the train test split. This can be controlled by the threshold parameter in the pbs script if the user wants to exclude clusters that only have a few samples. 

__Data:__

The script takes a single csv file of RNA expression changes over time. The first row is the header containing the gene column, the different time steps, and the cluster column. The other rows contain the gene, RNA expressions at each time step, and the cluster it belongs to. The csv file must be in the 'data' directory. Multiple csv files can be in this directory and the script will output a results directory for each file. An example is shown below, this example has one feature, the RNA expression, and five time steps.
```       
gene                      0             12            24            48            72            clusters    
A17C_Medtr2g026050_rep2	  5.221980164   5.157623152   5.294473715   5.29495152    5.37787568    Mtr-Shoot-DTW-K50-045-DPGP001
A17R_Medtr7g059515_rep2	  0             0.403085897   0.091756243   0.140124224   0.125210248   Mtr-Shoot-DTW-K50-044-DPGP002
```


__Output:__

The main output of the script is a confusion matrix and an F1 score for each model, these are created in a 'results' directory. The confusion matrix shows how well the trained model classified each sample from the test split. The expected values are on the y-axis and the predicted values are on the x-axis. This matrix is normalized to show a percentage of samples classfied. So if there are five genes in cluster A but the model predicted two of these genes are in cluster B, the cluster A row would could contain a 60 in column A and a 40 in column B. This normalization can be disabled in the confusion matrix function. The F1 scored is written to a text file

TODO: Write in insctructions to run script
  
