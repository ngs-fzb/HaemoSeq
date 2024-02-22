#!/usr/bin/env bash

###Author###
echo "hello World, this script was written by Margo Diricks!"

###Function###
#This script takes annotated assemblies in GFF3 format (e.g. PROKKA output), calculates the pan genome and a core gene alignment and constructs a phylogenetic tree. 

###Required packages###
#roary: https://github.com/sanger-pathogens/Roary and http://sanger-pathogens.github.io/Roary/
#Fasttree: https://anaconda.org/bioconda/fasttree
	
###Usage###
# bash file_name.sh

###Inputfiles###
#Type: annotated assemblies in GFF3 format (e.g. PROKKA output)
#Path to input files:
PATH_input="$HOME/YourFolder"

###Parameters###
#CorePercentage (percentage of strains that need to have the gene present before this gene is accepted as core gene)
cp="90"
#minimum percentage identity for blastp (default 95)
blastp=60

###Output files###
#1. Directory where the final files will be stored
PATH_output="$HOME/YourFolder"
#2. Name of tree file
Tree="BlastP_"$blastp"core_"$cp"_FastTree_GTR_gamma.newick"

#Amount Of Threads to be used
Threads=1

################################################################################################################################################################

#Pan-genome and core-genome analysis
roary -e -s --mafft -p $Threads -i $blastp -f $PATH_output -cd $cp $PATH_input/*.gff

#CreateTree
FastTree -nt -gtr -gamma < $PATH_output/core_gene_alignment.aln > $PATH_output/$Tree

################################################################################################################################################################
# Usage:   roary [options] *.gff

# Options: -p INT    number of threads [1]
         # -o STR    clusters output filename [clustered_proteins]
         # -f STR    output directory [.]
         # -e        create a multiFASTA alignment of core genes using PRANK
         # -n        fast core gene alignment with MAFFT, use with -e
         # -i        minimum percentage identity for blastp [95]
         # -cd FLOAT percentage of isolates a gene must be in to be core [99]
         # -qc       generate QC report with Kraken
         # -k STR    path to Kraken database for QC, use with -qc
         # -a        check dependancies and print versions
         # -b STR    blastp executable [blastp]
         # -c STR    mcl executable [mcl]
         # -d STR    mcxdeblast executable [mcxdeblast]
         # -g INT    maximum number of clusters [50000]
         # -m STR    makeblastdb executable [makeblastdb]
         # -r        create R plots, requires R and ggplot2
         # -s        dont split paralogs
         # -t INT    translation table [11]
         # -ap       allow paralogs in core alignment
         # -z        dont delete intermediate files
         # -v        verbose output to STDOUT
         # -w        print version and exit
         # -y        add gene inference information to spreadsheet, doesnt work with -e
         # -iv STR   Change the MCL inflation value [1.5]
         # -h        this help message

# Example: Quickly generate a core gene alignment using 8 threads
         # roary -e --mafft -p 8 *.gff 
# Options: -g STR    groups filename [clustered_proteins]
         # -a STR    action (union/intersection/complement/gene_multifasta/difference) [union]
         # -c FLOAT  percentage of isolates a gene must be in to be core [99]
         # -o STR    output filename [pan_genome_results]
         # -n STR    comma separated list of gene names for use with gene_multifasta action
         # -i STR    comma separated list of filenames, comparison set one
         # -t STR    comma separated list of filenames, comparison set two
         # -v        verbose output to STDOUT
         # -h        this help message
		 
