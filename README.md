# GeneShift
This workflow will capture **pattern and degree change over time** based on user input data. GeneShift workflow will perform following tasks:
  * OFF gene check
  * Cluster expression of individual replicate over time using [DTW-KMeans](https://arxiv.org/abs/1703.01541) and then fine clustering with [DP_GP_cluster](https://github.com/PrincetonUniversity/DP_GP_cluster/tree/master/DP_GP)
  * Using deep learning model(RNN LSTM) to test the robusticity of clustering 
  * Cluster finalization (remove poor quality gene based on replicate consistency) 
  * Conduct pattern and degree change capture and visualize the performance
  