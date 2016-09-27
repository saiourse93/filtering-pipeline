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
        rm $WORK_DIR/trinity_out_dir/Trinity.fasta

        if [ -s $WORK_DIR/trinity_out_dir/RESULTS/$file/Trinity.fasta ]
        then
            rm $WORK_DIR/trinity_out_dir/RESULTS/$file/Trinity.fasta
        fi

        /opt/exp_soft/bioinf/trinity/Trinity --seqType fq --max_memory 24G --left $DATA_DIR/$file/*filtered_R1.fastq \
        --right $DATA_DIR/$file/*filtered_R2.fastq --SS_lib_type RF --CPU 6

        if [ ! -d $WORK_DIR/trinity_out_dir/RESULTS ]
        then
            mkdir $WORK_DIR/trinity_out_dir/RESULTS
        fi

        if [ ! -d $WORK_DIR/trinity_out_dir/RESULTS/$file ]
        then
            mkdir $WORK_DIR/trinity_out_dir/RESULTS/$file
        fi

        cat $WORK_DIR/trinity_out_dir/Trinity.fasta > $WORK_DIR/trinity_out_dir/RESULTS/$file/Trinity.fasta

done <"$samples"