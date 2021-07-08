#!/bin/bash
 
BASEDIR=$(pwd)

Kmin=$1
Kmax=$2
StepSize=$3

mkdir -p ${BASEDIR}/DTWKMeans/PBS
mkdir -p ${BASEDIR}/DTWKMeans/Logs
cp ${BASEDIR}/basedir.txt ${BASEDIR}/DTWKMeans/PBS/

cp ${BASEDIR}/DTWKMeans.py ${BASEDIR}/DTWKMeans/PBS

for i in {$Kmin..$Kmax..$StepSize} ; do cat ${BASEDIR}/DTWKMeans.template | sed s/KX/${i}/g > ${BASEDIR}/DTWKMeans/PBS/K${i}.DTWKMeans.pbs ; done

cd ${BASEDIR}/DTWKMeans/PBS

for i in *DTWKMeans.pbs ; do qsub -o ${BASEDIR}/DTWKMeans/Logs/$i.out \
-e ${BASEDIR}/DTWKMeans/Logs/$i.err ./$i ; done