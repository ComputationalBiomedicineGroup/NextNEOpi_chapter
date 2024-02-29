#!/usr/bin/env nextflow

params.input = "/path/to/SRR_Acc_List.txt"

input_ch = Channel.fromPath(params.input).splitCsv()

process download_SRA_files {
    errorStrategy 'ignore'
    publishDir "/path/to/Sharma_fastq_out/", mode: "copy"

    input:
    tuple val(input) from input_ch

    output:
    path("**.fastq") into fastq_channel

    script:
    """
    /path/to/sratoolkit.3.0.10-centos_linux64/bin/fasterq-dump \
    -e 16 --split-files $input -t SRA_tmp
    """
}