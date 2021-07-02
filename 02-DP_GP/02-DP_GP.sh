#!/bin/bash
 
BASEDIR=$(pwd)
mkdir ${BASEDIR}/Logs


# Create Cluster List, individual pbs directory, output directory for each K
for i in `cat ${BASEDIR}/Loops.txt`
do 
    mkdir -p ${BASEDIR}/PBS/${i}-PBS
    cp ${BASEDIR}/DP_GP_cluster.py ${BASEDIR}/PBS/${i}-PBS
    cp ${BASEDIR}/basedir.txt  ${BASEDIR}/PBS/${i}-PBS/
    mkdir ${BASEDIR}/${i}DPGP_output
    cat ${BASEDIR}/Inputs/DTW-Loops-Results/${i}_ClusteringSummary.csv | awk -F, '{print $2}' | tail -n +2 > ${BASEDIR}/PBS/${i}-PBS/${i}-Clusters-List.txt
done

# Create PBS script for each DPGP job
for i in `cat ${BASEDIR}/Loops.txt`
do
    for j in `cat ${BASEDIR}/PBS/${i}-PBS/${i}-Clusters-List.txt`
    do cat ${BASEDIR}/DP_GP_cluster.template | sed s/KVal/${i}/g | sed s/clusterX/DTW-Cluster${j}/g > ${BASEDIR}/PBS/${i}-PBS/${i}_Cluster${j}.DPGP.pbs
    done
done

# Submit PBS jobs and write things in the log files
for i in `cat ${BASEDIR}/Loops.txt`
do 
    for j in ${BASEDIR}/PBS/${i}-PBS/*DPGP.pbs
    do
        jobName=$(echo $j | awk -F/ '{print $NF}')
        qsub -o ${BASEDIR}/Logs/${jobName}.out -e ${BASEDIR}/Logs/${jobName}.err ${j}
    done
done



