#!/bin/bash
 
BASEDIR=$(pwd)
INPTDIR=/scratch2/yueyaog/02-DP_GP/Inputs

# Create a list for K looping
#for i in `ls ${BASEDIR}/DTW-Loops-Results/DTW*` ; do echo $i | awk -F DP '{print$1}';done > ${BASEDIR}/Loops.txt
# Create Cluster List, individual pbs directory, output directory for each K
for i in `cat ${BASEDIR}/Loops.txt`
do 
    mkdir -p ${BASEDIR}/PBS/${i}-PBS
    cp ${BASEDIR}/DP_GP_cluster.py ${BASEDIR}/PBS/${i}-PBS
    mkdir ${BASEDIR}/${i}DPGP_output
    cat ${INPTDIR}/DTW-Loops-Results/${i}_ClusteringSummary.csv | awk -F, '{print $2}' | tail -n +2 > ${BASEDIR}/PBS/${i}-PBS/${i}-Clusters-List.txt
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
        #qsub -o ${BASEDIR}/Logs/${j}.out -e ${BASEDIR}/Logs/${j}.err $j
        qsub $j
    done
done

#mv ${BASEDIR}/*DTWdpgp.o* ${BASEDIR}/Logs/
#mv ${BASEDIR}/*DTWdpgp.e* ${BASEDIR}/Logs/

