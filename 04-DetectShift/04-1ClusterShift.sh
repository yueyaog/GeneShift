#!/bin/bash

BASEDIR=$(pwd)
mkdir $BASEDIR/Logs

mkdir -p OptimalK/ClusterPlots

qsub -o ${BASEDIR}/Logs/04-1ClusterShift.out -e ${BASEDIR}/Logs/04-1ClusterShift.err ./ClusterShift.template