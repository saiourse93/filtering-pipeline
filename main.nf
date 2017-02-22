#!/usr/bin/env nextflow

params.imsa = "$PWD/imsa" // The IMSA pipeline will come packaged! User will edit the config file.
params.data = "/spaces/phelelani/ssc_data/data_trimmed/inflated" // User to provide folder with own data. 
params.actions = "$PWD/actions.txt" // Action file will come packaged. User will edit.
params.out = "/spaces/phelelani/ssc_data/filtering/filtering_results"

imsa_path = params.imsa
out_path = file(params.out)
data_path = params.data
actions = params.actions
out_path.mkdir()

read_pair = Channel.fromFilePairs("${data_path}/*_R{1,2}.fastq", type: 'file')
//read_pair = Channel.fromFilePairs("${data_path}/caskiSubset_01_R{1,2}.fq", type: 'file')
//read_pair = Channel.fromFilePairs("${data_path}/CASE_01_ARM_001_R{1,2}.fastq", type: 'file')

//Initial pre-processing of the data before running the pipeline
process imsaMaster_process {
    cache true
    scratch '$HOME/tmp'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false, pattern: '*.py'

    input:
    set val(sample), file(reads) from read_pair
    
    output:
    set val(sample), file("${sample}_pipelineScript_bowtie2.py"), file("${sample}_pipelineScript_blat.py"), file("${sample}_pipelineScript_blast.py"), file(reads) into bowtie_in
    
    """
    python ${imsa_path}/master.py -i ${reads.get(0)} -j ${reads.get(1)} -p ${actions} -a 10 -o 33

    mv pipelineScript_bowtie2.py ${sample}_pipelineScript_bowtie2.py
    mv pipelineScript_blat.py ${sample}_pipelineScript_blat.py
    mv pipelineScript_blast.py ${sample}_pipelineScript_blast.py

    for item in \$(ls *.py)
    do
        sed -i 's@'"\$PWD"'@.@g' \$item
    done
    """
}

// Process 1: Mapping of reads to the reference genome
process runBowtie_process {
    cache true
    scratch '$HOME/tmp'
    executor 'pbs'
    queue 'WitsLong'
    cpus 11
    memory '50 GB'
    time '50h'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false, pattern: '*{.*.fq,.sam}'
    
    input:
    set val(sample), file(bowtie_script), file(blat_script), file(blast_script), file(reads) from bowtie_in
    
    output:
    set val(sample), file(blat_script), file(blast_script), file(reads), file("${sample}*{.fq,.sam}") into blat_in
    
    """
    python ${bowtie_script.getName()}
    """
}

// Process 2: 
process runBlat_process {
    cache true
    scratch '$HOME/tmp'
    executor 'pbs'
    queue 'WitsLong'
    cpus 1
    memory '20 GB'
    time '50h'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false, pattern: '*{.psl,.fa}'
    
    input:
    set val(sample), file(blat_script), file(blast_script), file(reads), file(bowtie_results) from blat_in
    
    output:
    set sample, file(blast_script), file(reads), file("${sample}*{.psl,.fa}") into blast_in
    
    """
    python ${blat_script.getName()}
    """
}

// Process 3:
process runBlastRef_process {
    cache true
    scratch '$HOME/tmp'
    executor 'pbs'
    queue 'WitsLong'
    cpus 11
    memory '200 GB'
    time '50h'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false, pattern: '*{.bln,.fa}'
    
    input:
    set sample, file(blast_script), file(reads), file(blat_results) from blast_in
    
    output:
    set sample,
    file(reads), 
    file("${sample}*{.bln,.fa}") into getFastq_in
    
    """
    python ${blast_script.getName()}
    """
}

// Process 4:
process getFastq_process {
    cache true
    scratch '$HOME/tmp'
    executor 'pbs'
    queue 'WitsLong'
    cpus 1
    memory '5 GB'
    time '50h'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false
    
    input:
    set sample, file(reads), file(blast_results) from getFastq_in
    
    output:
    set sample,
    file("${sample}_filtered_R{1,2}.fq") into trinity_in
    
    """
    grep '>' ${(blast_results as List).findAll { it =~ '.fa' }.max()} | sort -u > ${sample}_unique_names
    sed -i 's|>|@|g' ${sample}_unique_names
    sed -i 's|\$| |g' ${sample}_unique_names
    LC_ALL=C fgrep -F -f ${sample}_unique_names -A 3 --no-group-separator ${reads.get(0)} > ${sample}_filtered_R1.fq
    LC_ALL=C fgrep -F -f ${sample}_unique_names -A 3 --no-group-separator ${reads.get(1)} > ${sample}_filtered_R2.fq
    """
}

// Process 5:
process runTrinity_process {
    cache true
    scratch '$HOME/tmp'
    executor 'pbs'
    queue 'WitsLong'
    cpus 11
    memory '200 GB'
    time '50h'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false
    
    input:
    set sample, 
    file(reads) from trinity_in
    
    output:
    set sample, "trinity_${sample}/Trinity.fasta" into assemblies
    
    """
    Trinity --seqType fq --max_memory 150G --left ${reads.get(0)} --right ${reads.get(1)} --SS_lib_type RF --CPU 12 --output trinity_${sample}
    """
}

// Process 6:
process runBlastNt_process {
    cache true
    scratch '$HOME/tmp'
    executor 'pbs'
    queue 'WitsLong'
    cpus 11
    memory '200 GB'
    time '50h'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false
    
    input:
    set sample, trinity_assembly from assemblies
    
    output:
    set sample, "${sample}_blast_hits.bln" into blast_hits
    
    """
    blastn -query $trinity_assembly -db pathogens_32 -evalue 1e-15 -outfmt 6 -num_alignments 200 > ${sample}_blast_hits.bln
    """
}

// Process 7:
process runPostProcess_process {
    cache true
    scratch '$HOME/tmp'
    executor 'pbs'
    queue 'WitsLong'
    cpus 1
    memory '20 GB'
    time '50h'
    tag { sample }
    stageInMode 'symlink'
    stageOutMode 'rsync'
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false
    
    input:
    set sample, file(hits) from blast_hits
    
    output:
    set sample, file("${sample}*") into post_process
    
    """
    python ${imsa_path}/postprocess.py -b $hits
    """
}

// On completion
workflow.onComplete {
    println """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
}

