#!/bin/bash

#PBS -N RNA_Assembly
#PBS -q WitsLong
#PBS -l walltime=100:00:00,mem=24G
#PBS -l nodes=1:ppn=7
#PBS -o 
#PBS -e 

#Folders needed
WORK_DIR=

cd $WORK_DIR

Trinity --seqType fq --max_memory 24G --left $left_reads --right $right_reads --SS_lib_type RF --CPU 6