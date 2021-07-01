#!/bin/bash

pwd > basedir.txt
BASEDIR=$(pwd)

cp basedir.txt $BASEDIR/00-DataPrep

rm $BASEDIR/basedir.txt