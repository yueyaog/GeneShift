#!/bin/bash

pwd > basedir.txt
BASEDIR=$(pwd)

cp basedir.txt $BASEDIR/00-DataPrep
cp basedir.txt $BASEDIR/01-DTWKMeans
cp basedir.txt $BASEDIR/02-DP_GP

rm $BASEDIR/basedir.txt