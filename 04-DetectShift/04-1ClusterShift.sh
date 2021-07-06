#!/bin/bash

BASEDIR=$(pwd)


mkdir -p ${BASEDIR}/OptimalK$1/ClusterPlots
mkdir -p ${BASEDIR}/OptimalK$1/Logs
mkdir -p ${BASEDIR}/OptimalK$1/PBS

cat ${BASEDIR}/ClusterShift.template | sed s/KX/K$1/g > ${BASEDIR}/OptimalK$1/PBS/04-1ClusterShift.pbs

qsub -o ${BASEDIR}/OptimalK$1/Logs/04-1ClusterShift.out -e ${BASEDIR}/OptimalK$1/Logs/04-1ClusterShift.err ${BASEDIR}/OptimalK$1/PBS/04-1ClusterShift.pbs