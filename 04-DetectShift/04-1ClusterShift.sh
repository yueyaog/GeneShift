#!/bin/bash

BASEDIR=$(pwd)

OPTIMAL_K=$1

mkdir -p ${BASEDIR}/OptimalK$OPTIMAL_K/ClusterPlots
mkdir -p ${BASEDIR}/OptimalK$OPTIMAL_K/Logs
mkdir -p ${BASEDIR}/OptimalK$OPTIMAL_K/PBS
mkdir -p ${BASEDIR}/OptimalK$OPTIMAL_K/REP2
mkdir -p ${BASEDIR}/OptimalK$OPTIMAL_K/REP3

cat ${BASEDIR}/ClusterShift.pbs.template | sed s/KX/K$OPTIMAL_K/g > ${BASEDIR}/OptimalK$OPTIMAL_K/PBS/04-1ClusterShift.pbs

qsub -o ${BASEDIR}/OptimalK$OPTIMAL_K/Logs/04-1ClusterShift.out -e ${BASEDIR}/OptimalK$OPTIMAL_K/Logs/04-1ClusterShift.err ${BASEDIR}/OptimalK$OPTIMAL_K/PBS/04-1ClusterShift.pbs