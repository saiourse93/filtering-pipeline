#!/bin/bash

#PBS -N RNA_Assembly
#PBS -q WitsLong
#PBS -l walltime=100:00:00,mem=100gb
#PBS -l nodes=1:ppn=7

#Folders needed
WORK_DIR=/home/phelelani/scleroderma_analysis/jobs/assembly-pipeline/tools # Change to your working directory!

cd $WORK_DIR

Trinity --seqType fq --max_memory 24G --left CASE_01_ARM_001_filtered_R1.fastq --right CASE_01_ARM_001_filtered_R2.fastq --SS_lib_type RF --CPU 6
