#!/bin/bash
 
BASEDIR=$(pwd)

mkdir -p ${BASEDIR}/ClusterSummary/Logs
mkdir -p ${BASEDIR}/ClusterSummary/PBS

cp ${BASEDIR}/Cluster_Compile.py ${BASEDIR}/ClusterSummary/PBS

cp ${BASEDIR}/basedir.txt ${BASEDIR}/ClusterSummary/PBS

for i in {10..50..5} ; do cat ${BASEDIR}/Cluster_Summary.template | sed s/KX/${i}/g > ${BASEDIR}/ClusterSummary/PBS/K${i}.Cluster_Sum.pbs ; done

cd ${BASEDIR}/ClusterSummary/PBS

for i in *Cluster_Sum.pbs ; do qsub -o ${BASEDIR}/ClusterSummary/Logs/${i}.out -e ${BASEDIR}/ClusterSummary/Logs/$i.err ./$i ; done