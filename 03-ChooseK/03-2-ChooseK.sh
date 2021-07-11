#!/bin/bash

BASEDIR=$(pwd)
mkdir $BASEDIR/Logs

qsub -o ${BASEDIR}/Logs/ChooseK.out -e ${BASEDIR}/Logs/ChooseK.err ./ChooseK.pbs.template