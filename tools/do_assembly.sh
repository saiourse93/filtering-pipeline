#!/bin/bash                                                                                                                                     

#PBS -N RNA_Assembly                                                                                                                            
#PBS -q WitsLong                                                                                                                                
#PBS -l walltime=100:00:00,mem=100gb                                                                                                            
#PBS -l nodes=1:ppn=7                                                                                                                           

#Folders needed                                                                                                                                 
WORK_DIR=/home/mkhari/git/assembly-pipeline/tools 
DATA_DIR=/global/blast/test_data/assembly

cd $WORK_DIR
samples=sample_dirs

while IFS= read -r file
do
        if [ ! -d $WORK_DIR/TRINITY_RESULTS ]
        then
            mkdir $WORK_DIR/TRINITY_RESULTS
        fi

        if [ ! -d $WORK_DIR/TRINITY_RESULTS/$file ]
        then
            mkdir $WORK_DIR/TRINITY_RESULTS/$file
        fi

        if [ -s $WORK_DIR/TRINITY_RESULTS/$file/trinity_out_dir/Trinity.fasta ]
        then
            rm $WORK_DIR/TRINITY_RESULTS/$file/trinity_out_dir/Trinity.fasta
        fi

	cp $DATA_DIR/$file/*filtered_R1.fastq $WORK_DIR/TRINITY_RESULTS/$file
	cp $DATA_DIR/$file/*filtered_R2.fastq $WORK_DIR/TRINITY_RESULTS/$file

	cd $WORK_DIR/TRINITY_RESULTS/$file

        /opt/exp_soft/bioinf/trinity/Trinity --seqType fq --max_memory 24G --left *filtered_R1.fastq \
        --right *filtered_R2.fastq --SS_lib_type RF --CPU 6

done <"$samples"