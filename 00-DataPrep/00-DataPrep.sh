#!/bin/bash

BASEDIR=$(pwd)
mkdir $BASEDIR/Logs

qsub -o ${BASEDIR}/Logs/00InputPrep.out -e ${BASEDIR}/Logs/00InputPrep.err ./Input_prep.template