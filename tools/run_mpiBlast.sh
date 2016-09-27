#!/bin/bash                                                                                                                                     
#PBS -N run_mpiBlast_bench
#PBS -q WitsLong
#PBS -l walltime=100:00:00,mem=250gb
#PBS -l nodes=7:ppn=2 
#PBS -o /home/mkhari/git/assembly-pipeline/tools/mpi_run.out
#PBS -e /home/mkhari/git/assembly-pipeline/tools/mpi_run.err

WORK_DIR=/home/mkhari/git/assembly-pipeline/tools

cd $WORK_DIR

DATA_DIR=/home/mkhari/git/mpi-blast-benchmarking/tools/trinity_out_dir/RESULTS

samples=sample_dirs
params=parameters_pipe

# This function will run mpiBlast on the test sample with the parameters specified
# Needs input :
# -input fastq file
#
#

do_mpi_blast () {
    module load mpich-x86_64
    /usr/bin/time mpi-start -t mpich2 -np $1 mpiblast -d $2 \
        -i $3 \
        -p blastn \
        -e 1e-15 -m 8 -b 200 \
        -o $4 \
        --time-profile=time.txt \
        --partition-size=$5 \
        --replica-group-size=$6 \
        --use-virtual-frags \
        --query-segment-size=$7 \
        --use-parallel-write \
        --predistribute-db
}

while IFS= read -r file
do
    cd $DATA_DIR/$file
    # Get the file with the parameters for each mpi-blast
    while IFS='|' read -r num_proc db in_file out_file p_size rg_size qs_size
    do
        # Run mpiBlast if the final output file does not exist!
        if [ ! -s $out_file ] 
        then
	    # run mpiBlast
	    do_mpi_blast $num_proc $db $in_file $out_file $p_size $rg_size $qs_size
	
	    # Join the temporary output files into the final output file
	    cat *.part* > $out_file
        fi
    done < "$params"
done < "samples"