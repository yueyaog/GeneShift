#!/bin/bash
 
BASEDIR=$(pwd)

Kmin=$1
Kmax=$2
StepSize=$3


mkdir -p ${BASEDIR}/DTWKMeans/PBS
mkdir -p ${BASEDIR}/DTWKMeans/Logs
cp ${BASEDIR}/basedir.txt ${BASEDIR}/DTWKMeans/PBS/

cp ${BASEDIR}/DTWKMeans.py ${BASEDIR}/DTWKMeans/PBS


# UPDATE Dec27: The following command is not working in palmetto shell env anymore. The input arguments are not treated as integers in shell env. I currently update it into seq function.
#for i in {$Kmin..$Kmax..$StepSize} ; do cat ${BASEDIR}/DTWKMeans.pbs.template | sed s/KX/${i}/g > ${BASEDIR}/DTWKMeans/PBS/K${i}.DTWKMeans.pbs ; done

for i in `seq $Kmin $StepSize $Kmax` ; do cat ${BASEDIR}/DTWKMeans.pbs.template | sed s/KX/${i}/g > ${BASEDIR}/DTWKMeans/PBS/K${i}.DTWKMeans.pbs ; done

cd ${BASEDIR}/DTWKMeans/PBS

for i in *DTWKMeans.pbs ; do qsub -o ${BASEDIR}/DTWKMeans/Logs/$i.out \
-e ${BASEDIR}/DTWKMeans/Logs/$i.err ./$i ; done