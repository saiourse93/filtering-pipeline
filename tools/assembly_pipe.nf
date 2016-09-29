#!/usr/bin/env nextflow

sample_dir = Channel.from(
           'CASE_01_ARM_001','CASE_02_BAK_001','CASE_04_ARM_001','CASE_06_ARM_001','CASE_07_BAK_001','CASE_09_BAK_001','CTRL_03_BAK_001',
           'CTRL_08_BAK_001','CASE_01_BAK_001','CASE_03_ARM_001','CASE_04_BAK_001','CASE_06_ARM_002','CASE_08_ARM_001','CTRL_01_BAK_001',
           'CTRL_05_BAK_001','CASE_01_BAK_002','CASE_03_ARM_002','CASE_05_ARM_001','CASE_06_BAK_001','CASE_08_BAK_001','CTRL_01_BAK_002',
           'CTRL_06_BAK_001','CASE_02_ARM_001','CASE_03_BAK_001','CASE_05_BAK_001','CASE_07_ARM_001','CASE_09_ARM_001','CTRL_02_BAK_001',
           'CTRL_07_BAK_001')

WORK_DIR='/home/mkhari/git/assembly-pipeline/tools'
DATA_DIR='/global/blast/test_data/assembly'
Trinity='/opt/exp_soft/bioinf/trinity/Trinity'

process trinity {

	executor = 'pbs'
	queue = 'WitsLong'
	cpus = 5
	memory = '2GB'
	time = '1h'
	penv 'smp'
               
	input:
	val sample_dir

	output:
	file 'Trinity.fasta' into results
 
	"""
	cd $WORK_DIR/work
	mkdir $sample_dir
	cd $sample_dir

	cp $DATA_DIR/$sample_dir/*filtered_R1.fastq $WORK_DIR/work/$sample_dir
	cp $DATA_DIR/$sample_dir/*filtered_R2.fastq $WORK_DIR/work/$sample_dir

	$Trinity --seqType fq --max_memory 24G --left *filtered_R1.fastq \
	--right *filtered_R2.fastq --SS_lib_type RF --CPU 6
	"""

}
