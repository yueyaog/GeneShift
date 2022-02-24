nextflow.enable.dsl=2

println """\
====================================
G E N E S H I F T   W O R K F L O W
====================================

Workflow Information:
---------------------
  Launch Directory:   ${workflow.launchDir}
  Work Directory:     ${workflow.workDir}
  Config Files:       ${workflow.configFiles}
  Profiles:           ${workflow.profile}
  
Execution Parameters:
---------------------
input:
    gem_file: ${params.gem_file}
    Kmin: ${params.Kmin}
    Kmax: ${params.Kmax}
    StepSize: ${params.StepSize}
    iteration: ${params.iteration}
    

    

"""





workflow {  // Create a K range list 
            def emptyList = []
            params.Kmin.step(params.Kmax, params.StepSize){emptyList.add("$it")}
            //Load Input
            gem_files = Channel.fromList([ params.gem_file ])
            
           //Prep Input Files 
           data_prep(gem_files)
           off_gem_file = data_prep.out.off_gem_file
           off_removed_gem_file = data_prep.out.off_removed_gem_file
           
           //Perform Initial Clustering Using DTWKMeans
           ch = Channel.fromList(emptyList)   
           DTWKMeans(off_removed_gem_file, ch)
           DTW_Clustered_files = DTWKMeans.out.DTW_Clustered_files
           
           //Perform Fine Clustering with DP_GP
           DTW_files = Channel.fromPath("DTWKMeans/DTW-K*/DTW-K*DPGP_input/DTW-K*-Cluster*.txt").flatten()
           DP_GP(DTW_Clustered_files.flatten())
           DPGP_Clustered_files = DP_GP.out.DPGP_Clustered_files.collect()
           
           //Compiling Clustering Results for Each K Value
           k = Channel.fromList(emptyList) 
           Compile(k,off_removed_gem_file,DPGP_Clustered_files)
           Compiled_Clustered_files = Compile.out.Compiled_Clustered_files.collect()
           
           //Choosing the Optimal Number of Clusters 
           OptimalK(off_removed_gem_file,Compiled_Clustered_files)
           Optimal_Kval_file = OptimalK.out.Optimal_Kval_file
           
           //Detecting Cluster Changes between Conditions
           DetectShift(Optimal_Kval_file,off_removed_gem_file,off_gem_file)
           Cluster_Shift_File = DetectShift.out.Cluster_Shift_File
           Shift_Summary_File  = DetectShift.out.Shift_Summary_File 
           
           //Visualization of GeneShift Outputs
            gem_file = Channel.fromList([ params.gem_file ])
            Visualization(gem_file,Cluster_Shift_File,Shift_Summary_File)
            GeneShift_LinePlots = Visualization.out.Line_Plots
            GeneShift_BoxPlots = Visualization.out.Box_Plots

          
           
}

/**
 * Prepare input for GeneShift
 * Seperate the user input GEM into OFF_exp.csv and OFFremoved_exp.csv
 */
 process data_prep {
    publishDir "${params.output_dir}"

    input:
        val(gem_file)
        
    output:
        path("*OFFremoved_exp.csv"),emit: off_removed_gem_file
        path("*OFF_exp.csv"),emit: off_gem_file
        //path("*OFF*_exp.csv"), emit: off_gem_file

    script:
    """
    module load anaconda3/5.1.0-gcc 
    source activate ${params.conda_env}
    python ${GeneShift_PATH}/InputPrep.py\
            -i ${gem_file} \
            -o ${params.data_prefix} \

    """
}

/**
 * Perform initial clustering on OFF_removed gem file. 
 */
           
        
process DTWKMeans {
     tag "k${ch}"
     publishDir "${params.output_dir}"
     
     input: 

         path(off_removed_gem_file)
         each ch
         
     output:
         path("DTWKMeans/DTW-K*/DTW-K*DPGP_input/DTW-K*-Cluster*.txt"),emit: DTW_Clustered_files
         path("DTWKMeans/DTW-K*/DTW-K*_ClusteringResults.csv"),emit: DTW_Clustering_results
         path("DTWKMeans/DTW-K*/DTW-K*_ClusteringSummary.csv"),emit: DTW_Clustering_sum
         
     script:
     """
     
        module add anaconda3/5.1.0-gcc/8.3.1

        source activate deep-learning
        
        cd \$PBS_O_WORKDIR
        mkdir \$(pwd)/DTWKMeans/
        
        python ${GeneShift_PATH}/DTWKMeans.py -i ${off_removed_gem_file} -k ${ch} -o \$(pwd)/DTWKMeans/DTW-K${ch} -p DTW-K${ch}

    """
        
 
 }
 


/**
  * Perform fine clustering on DTW_Clustered_files
  */
process DP_GP{
    tag "${DTW_file.baseName}"
    publishDir "${params.output_dir}"
     
    input:
        file DTW_file
        
     
    output:
        path("DPGP/*optimal_clustering.txt"), emit: DPGP_Clustered_files
        
    script:
    """
    module load anaconda3/5.1.0-gcc 
    source activate ${params.dp_gp_env}
    
    cd \$PBS_O_WORKDIR
    mkdir \$(pwd)/DPGP/
 
    python ${GeneShift_PATH}/DP_GP_cluster.py\
            -i ${DTW_file} \
            -o \$(pwd)/DPGP/${DTW_file.baseName} \
            -n ${params.iteration} \
            --fast \
            --true_times
    """
}

/**
  * Compiling Clustering Results for Each K Value
  */

process Compile{
     tag "Compile${k}"
     publishDir "${params.output_dir}"
     
     input: 
         each k 
         file off_removed_gem_file
         file DPGP_Clustered_file

       
     output:
         path("ClusterSummary/DTW-K*-DPGP_ResultsSum.csv"),emit: Compiled_Clustered_files
         
     
     script:
     """
        module add anaconda3/5.1.0-gcc/8.3.1

        source activate deep-learning
        
        cd \$PBS_O_WORKDIR
        mkdir \$(pwd)/ClusterSummary/
        #echo ${DPGP_Clustered_file}
        
        python ${GeneShift_PATH}/Cluster_Compile.py\
        -i ${params.output_path}/DPGP\
        -emx ${off_removed_gem_file}\
        -k ${k}\
        -o \$(pwd)/ClusterSummary
 
        
        

    """
        
}

/**
  * Choosing the Optimal Number of Clusters
  */
process OptimalK{
     publishDir "${params.output_dir}"
     
     input: 
         file off_removed_gem_file
         file Compiled_Clustered_file

         

     output:
         path("OptimalK/DTW-DPGP_diffKs_PERFs.png"),emit: diffK_plot
         path("OptimalK/Optimal_K_value.txt"), emit: Optimal_Kval_file
         
     
     script:
     """
        module add anaconda3/5.1.0-gcc/8.3.1

        source activate deep-learning
        
        cd \$PBS_O_WORKDIR
        mkdir \$(pwd)/OptimalK
        echo ${Compiled_Clustered_file}
        
        python ${GeneShift_PATH}/OptimalK.py\
        -i ${params.output_path}/ClusterSummary\
        -o \$(pwd)/OptimalK\
        -emx ${off_removed_gem_file}\
        -kmin ${params.Kmin}\
        -kmax ${params.Kmax}\
        -step ${params.StepSize}

    """
        
}



/**
  * Post-clustering analysis (replicate sorting)
  */
  
process DetectShift{
     publishDir "${params.output_dir}"
     
     input: 
         file Optimal_Kval_file
         file off_removed_gem_file
         file off_gem_file

     output:
         path("OptimalK/OptimalK_plots/*_gene_exp_fig_*.png"),emit: OptimalK_ClusterPlots
         path("OptimalK/REP3_Outputs/*GeneCluster_Consist_Filter.csv"), emit: Replicate_Consistence_File
         path("OptimalK/REP3_Outputs/*_DTW_DPGP_ClusterShift.csv"), emit: Cluster_Shift_File
         path("OptimalK/REP3_Outputs/*ClusterShift_Summary.csv"), emit: Shift_Summary_File         

     script:
     """
        module add anaconda3/5.1.0-gcc/8.3.1

        source activate deep-learning
        
        cd \$PBS_O_WORKDIR
        mkdir -p \$(pwd)/OptimalK/OptimalK_plots
        mkdir -p \$(pwd)/OptimalK/REP3_Outputs
        
        
        k=`cat ${Optimal_Kval_file}`


        #Plot the optimal cluster results
        python ${GeneShift_PATH}/Clusters_plot.py\
        -exm ${off_removed_gem_file}\
        -c ${params.output_path}/ClusterSummary/DTW-K\${k}-DPGP_ResultsSum.csv\
        -o \$(pwd)/OptimalK/OptimalK_plots/K\${k}

        #Detect the shift under two conditions
        #Replicate threshold 3

        python ${GeneShift_PATH}/Cluster_change.py\
        -r 3\
        -off ${off_gem_file}\
        -c ${params.output_path}/ClusterSummary/DTW-K\${k}-DPGP_ResultsSum.csv\
        -o \$(pwd)/OptimalK/REP3_Outputs/ClusterShiftREP3
        

    """
        
}

/**
  * Visualization of GeneShift Outputs
  */
process Visualization{
     publishDir "${params.output_dir}"
     
     input: 
         val gem_file
         file Cluster_Shift_File
         file Shift_Summary_File

     output:
         path("GeneShift_Plots/Group*/LnPt_FPKM_TrajecotrySet_*.png"),emit: Line_Plots
         path("GeneShift_Plots/Group*/BxPt_FPKM_TrajecotrySet_*.png"),emit: Box_Plots
         
     
     script:
     """
        module add anaconda3/5.1.0-gcc/8.3.1

        source activate deep-learning
        
        cd \$PBS_O_WORKDIR
        mkdir \$(pwd)/GeneShift_Plots
        
        python ${GeneShift_PATH}/ClusterShift_plot.py\
        -exm ${gem_file} \
        -c ${Cluster_Shift_File} \
        -s ${Shift_Summary_File} \
        -o \$(pwd)/GeneShift_Plots

    """
        
}
