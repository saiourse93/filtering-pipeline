#!/usr/bin/env nextflow

process trinity {

    """
    qsub /home/mkhari/git/assembly-pipeline/tools/do_assembly.sh 
    """

}

process mpi {

    """
    qsub /home/mkhari/git/assembly-pipeline/tools/run_mpiBlast.sh
    """

}

