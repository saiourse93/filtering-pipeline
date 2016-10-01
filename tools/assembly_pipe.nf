#!/usr/bin/env nextflow

the_read_pair = Channel.fromFilePairs('/global/blast/test_data/assembly/{CASE,CRTL}_*/*_filtered_R{1,2}.fastq')

Trinity='/opt/exp_soft/bioinf/trinity/Trinity'

OUT_DIR="$HOME/RESULTS"

process trinity {
	executor = 'pbs'
	queue = 'WitsLong'
	cpus = 9
	memory = '100GB'
	time = '3h'

	publishDir "$OUT_DIR" , mode:'symlink', overwrite: true               

	input:
		set sample, file(reads) from the_read_pair

		output:
			set sample, "trinity_${sample}/Trinity.fasta" into assemblies
 
 """
 $Trinity --seqType fq --max_memory 100G --left ${reads[0]} --right ${reads[1]} --SS_lib_type RF --CPU 8

 mkdir trinity_${sample}
 cp trinity_out_dir/Trinity.fasta trinity_${sample}
 """
}

//assemblies.subscribe{ println "$it" }

process mpi_blast {
	executor = 'pbs'
	queue = 'WitsLong'
	cpus = 9
	memory = '100GB'
	time = '3h'

	publishDir "$OUT_DIR" , mode:'symlink', overwrite: true

	input:
		set sample, file(fasta) from assemblies
		output:
			set sample, "blast_${sample}/results.bln" into blast_hits

    """
    blastn -query ${fasta} -db pathogens_32 -out results.bln -evalue 1e-15 -outfmt 6 -num_alignments 200 -num_threads 8

    mkdir blast_${sample}
    cp results.bln blast_${sample}
    """

}

blast_hits.subscribe{ println "$it" }