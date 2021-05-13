#!/bin/bash
 
BASEDIR=$(pwd)

mkdir ${BASEDIR}/MtrRoot_DTWKMeans/PBS
mkdir ${BASEDIR}/MtrRoot_DTWKMeans/Logs

cp ${BASEDIR}/DTWKMeans.py ${BASEDIR}/MtrRoot_DTWKMeans/PBS

for i in {35..90..5} ; do cat ${BASEDIR}/DTWKMeans.template | sed s/KX/${i}/g > ${BASEDIR}/MtrRoot_DTWKMeans/PBS/K${i}.DTWKMeans.pbs ; done

cd ${BASEDIR}/MtrRoot_DTWKMeans/PBS

for i in *DTWKMeans.pbs ; do qsub -o ${BASEDIR}/MtrRoot_DTWKMeans/Logs/$i.out \
-e ${BASEDIR}/MtrRoot_DTWKMeans/Logs/$i.err ./$i ; done