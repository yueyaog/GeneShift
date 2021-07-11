#!/bin/bash

BASEDIR=$(pwd)

OPTIMAL_K=$1

mkdir -p  $BASEDIR/OptimalK$OPTIMAL_K/GeneShift_Plots/REP2
mkdir -p  $BASEDIR/OptimalK$OPTIMAL_K/GeneShift_Plots/REP3

cat ${BASEDIR}/ClusterShift_Visu.pbs.template | sed s/KX/K$OPTIMAL_K/g > ${BASEDIR}/OptimalK$OPTIMAL_K/PBS/04-2ClusterShift_Visu.pbs

qsub -o ${BASEDIR}/OptimalK$OPTIMAL_K/Logs/04-2ClusterShift_Visu.out -e ${BASEDIR}/OptimalK$OPTIMAL_K/Logs/04-2ClusterShift_Visu.err ${BASEDIR}/OptimalK$OPTIMAL_K/PBS/04-2ClusterShift_Visu.pbs