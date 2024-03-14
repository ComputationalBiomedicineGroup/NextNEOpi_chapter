# Repository for the nextNEOpi Chapter scripts

Within this chapter, we utilize the ```nextNEOpi``` pipeline to conduct an analysis of genomic and transcriptomic data derived from several biopsy sites of three patients with lung tumors (```Sharma et al. 2019```). This investigation illustrates how nextNEOpi aids in detecting potential canonical neoantigens, neoantigens arising from fusions, and other cancer-related attributes crucial for immunotherapy.

This repository houses a selection of scripts essential for replicating the analysis outlined in the chapter:

The list of SRA sample accessions (```SRR_Acc_List.txt```) and the associated metadata (```SraRunTable.txt```).

A Nextflow script (```SRA_download.nf```), designed to download the raw sequencing data using the SRA toolkit, with input from the ```SRR_Acc_List.txt``` file, and its accompanying Nextflow configuration file (```nextflow.config```). You may need to modify the ```nextflow.config``` to accommodate your computational environment.

An R script (```Sharma_nextneopi_samplesheet.R```), which processes the ```SraRunTable.txt``` file to generate the CSV input sample sheet for nextNEOpi, plus an example sample sheet (```Sharma_nextneopi_samplesheet.csv```). Installation of the ```dplyer```, ```readr```, and ```stringr``` R packages is required for this script.

An additional R script for visualizing the analysis outcomes (```Sharma_nextneopi_analysis.R```). This script necessitates the installation of several R packages: ```stringr```, ```readr```, ```ggplot2```, ```dplyr```, ```patchwork```, ```vegan```, and ```DiversitySeq``` (available at https://sysbiobig.dei.unipd.it/software/diversityseq/).

For detailed instructions on replicating the analysis, please consult the chapter itself.
