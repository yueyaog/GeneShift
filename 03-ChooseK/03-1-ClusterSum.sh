#!/bin/bash
 
BASEDIR=$(pwd)

mkdir ${BASEDIR}/ClusterSumPBS

cp ${BASEDIR}/Cluster_Compile.py ${BASEDIR}/ClusterSumPBS/

for i in {35..90..5} ; do cat ${BASEDIR}/Cluster_Summary.template | sed s/KX/${i}/g > ${BASEDIR}/ClusterSumPBS/K${i}.Cluster_Sum.pbs ; done

cd ${BASEDIR}/ClusterSumPBS/

for i in *Cluster_Sum.pbs ; do qsub ./$i ; done