#!/bin/bash
 
BASEDIR=$(pwd)

mkdir -p ${BASEDIR}/DTWKMeans/PBS
mkdir -p ${BASEDIR}/DTWKMeans/Logs

cp ${BASEDIR}/DTWKMeans.py ${BASEDIR}/DTWKMeans/PBS

for i in {10..50..5} ; do cat ${BASEDIR}/DTWKMeans.template | sed s/KX/${i}/g > ${BASEDIR}/DTWKMeans/PBS/K${i}.DTWKMeans.pbs ; done

cd ${BASEDIR}/DTWKMeans/PBS

for i in *DTWKMeans.pbs ; do qsub -o ${BASEDIR}/DTWKMeans/Logs/$i.out \
-e ${BASEDIR}/DTWKMeans/Logs/$i.err ./$i ; done