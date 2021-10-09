#!/usr/bin/bash



# SPEC=$1
TUNNELS=/home/rxz/dev/ribxz/ciftools/TUNNELS/$SPEC/*

parallel 'python3 driver.py -pdbid {1}' \
	::: $(for folder in $TUNNELS; do echo $(basename $folder); done) \


