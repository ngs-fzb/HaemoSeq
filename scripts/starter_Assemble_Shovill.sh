#!/usr/bin/env bash

###Author###
echo "hello World, this script was written by Margo Diricks!"

###Function###
#This script runs shovill on a set of fastq files; Parameters (trimming, coverage,...) can be adjusted.

###Required packages###
	#https://githuhttps://github.com/tseemann/shovill

###Usage###
# bash $HOME/Scripts/Name_Script.sh"

###Inputfiles###
#Directory where you have stored the (raw) paired-end illumina fastQ files
PATH_input="$HOME/YourFolder"

#Read notation
FastqType=".fastq.gz" # Extension of your fastq file
Fw="_R1" # Forward read notation
Rv="_R2" # Reverse read notation

###Parameters###
#Expected genome size
egs="2.0M"
#Downsampling to
cov=100
#Assembler (choose skesa velvet megahit or spades)
ass="spades"

###Output files###
#1. Directory where the final files will be stored
PATH_output="$HOME/YourFolder"

#Amount of cpus to be used (0=All)
cpu=0

################################################################################################################################################################
mkdir $PATH_output
mkdir $PATH_output/FinalAssemblies
rm $PATH_output/FinalAssemblies/Failed.txt

#Loop through all the raw fastQ files
for fastq in $PATH_input/*$Fw$FastqType 
do
	SampleName=$(basename $fastq| cut -d '_' -f 1)
	if [ -d $PATH_output/$SampleName ]
	then
		echo "Sample was already analysed"
	else
		shovill --R1 $fastq --R2 $(echo $fastq | sed "s/$Fw$FastqType$/$Rv$FastqType/1") --outdir $PATH_output/$SampleName --gsize $egs --depth $cov --trim --noreadcorr --assembler $ass 

#Move the final assembly to a seperate folder and rename file
		if [ -f "$PATH_output/$SampleName/contigs.fa" ]
		then
			mv $PATH_output/$SampleName/contigs.fa $PATH_output/FinalAssemblies/$SampleName"_Assembly_Shovill_"$ass.fasta
		else
			echo $SampleName >> $PATH_output/FinalAssemblies/Failed.txt
		fi
	fi
done

# Create summary file (as quality check)
#Make header of the summary file 
echo -e "Sample\tEstimated_depth(x)\tRead_max_len\tRead_avg_len\tRead_min_len\tSurviving_Read_Pairs_trimmomatic_perc\tFw_surviving_trimmomatic_perc\tRv_surviving_trimmomatic_perc\tDropped_trimmomatic_perc\tFlash_CombinedPairs_perc\tWalltime\tContigAmount\tMinContigLength\tAssemblyLength\tAssemblyRel%("$egs")" >> $PATH_output/"shovillOutput_summary.txt"

#Loop through all the generated log files and extract necessary info
for log in $PATH_output/*/shovill.log
do
	#ReadStats
	max_len=$(grep 'max_len' $log| cut -d '=' -f 2)
	avg_len=$(grep 'avg_len' $log| cut -d '=' -f 2)
	min_len=$(grep 'min_len' $log| cut -d '=' -f 2)
	EstSeqDepth=$(grep 'Estimated.sequencing.depth' $log| cut -d ':' -f 2 | sed 's/\s//' | sed 's/x//')
	Flash_CombinedPairs_perc=$(grep -o 'Percent.combined.*%' $log| cut -d ':' -f 2 | sed 's/\s//' | sed 's/%//')
	Surviving_Read_Pairs_trimmomatic_perc=$(grep -Po 'Both.Surviving.*?%' $log| grep -Po '.....%' | sed 's/%//') #-P means that you want to use non-greedy expression .*? which means any character multiple times until the first character you put after this expression; -o means output only matching part, not complete line
	Fw_surviving_trimmomatic_perc=$(grep -P -o 'Forward.Only.Surviving.*?%' $log| grep -Po '....%' | sed 's/%//')
	Rv_surviving_trimmomatic_perc=$(grep -P -o 'Reverse.Only.Surviving.*?%' $log| grep -Po '....%' | sed 's/%//')
	Dropped_trimmomatic_perc=$(grep -Po 'Dropped.*?%' $log| grep -Po '....%' | sed 's/%//')
	Walltime=$(grep 'Walltime' $log| cut -d ':' -f 2 | sed 's/\s//')
	Contig_amount=$(grep 'It.contains' $log| cut -d ' ' -f 4)
	min_contig_length=$(grep -P -o '\(min=.*?\)' $log| cut -d '=' -f 2 | sed -e 's/)//')
	Assembly_length=$(grep -P -o 'Assembly.is.*?,' $log| cut -d ' ' -f 3 | sed 's/,//')
	Assembly_rel=$(grep 'Assembly.is' $log| grep -Po '\(.*?\)'| sed -e 's/)//' | sed -e 's/(//' | sed -e 's/%//')
	
	paste --delimiter='\t'  <(basename $(dirname $log))	<(echo $EstSeqDepth) <(echo $max_len) <(echo $avg_len) <(echo $min_len) <(echo $Surviving_Read_Pairs_trimmomatic_perc) <(echo $Fw_surviving_trimmomatic_perc) <(echo $Rv_surviving_trimmomatic_perc) <(echo $Dropped_trimmomatic_perc) <(echo $Flash_CombinedPairs_perc) <(echo $Walltime) <(echo $Contig_amount) <(echo $min_contig_length) <(echo $Assembly_length) <(echo $Assembly_rel) >> $PATH_output/"shovillOutput_summary.txt"
done

################################################################################################################################################################
#QUICK START: shovill --outdir out --R1 test/R1.fq.gz --R2 test/R2.fq.gz

#MAIN STEPS
# Estimate genome size and read length from reads (unless --gsize provided)
# Reduce FASTQ files to a sensible depth (default --depth 100)
# Trim adapters from reads (with --trim only)
# Conservatively correct sequencing errors in reads
# Pre-overlap ("stitch") paired-end reads
# Assemble with SPAdes/SKESA/Megahit with modified kmer range and PE + long SE reads
# Correct minor assembly errors by mapping reads back to contigs
# Remove contigs that are too short, too low coverage, or pure homopolymers
# Produce final FASTA with nicer names and parseable annotations

# SYNOPSIS
  # De novo assembly pipeline for Illumina paired reads
# USAGE
  # shovill [options] --outdir DIR --R1 R1.fq.gz --R2 R2.fq.gz
# GENERAL
  # --help          This help
  # --version       Print version and exit
  # --check         Check dependencies are installed
# INPUT
  # --R1 XXX        Read 1 FASTQ (default: '')
  # --R2 XXX        Read 2 FASTQ (default: '')
  # --depth N       Sub-sample --R1/--R2 to this depth. Disable with --depth 0 (default: 150)
  # --gsize XXX     Estimated genome size eg. 3.2M <blank=AUTODETECT> (default: '')
# OUTPUT
  # --outdir XXX    Output folder (default: '')
  # --force         Force overwite of existing output folder (default: OFF)
  # --minlen N      Minimum contig length <0=AUTO> (default: 0)
  # --mincov n.nn   Minimum contig coverage <0=AUTO> (default: 2)
  # --namefmt XXX   Format of contig FASTA IDs in 'printf' style (default: 'contig%05d')
  # --keepfiles     Keep intermediate files (default: OFF)
# RESOURCES
  # --tmpdir XXX    Fast temporary directory (default: '/tmp/tseemann')
  # --cpus N        Number of CPUs to use (0=ALL) (default: 8)
  # --ram n.nn      Try to keep RAM usage below this many GB (default: 16)
# ASSEMBLER
  # --assembler XXX Assembler: skesa velvet megahit spades (default: 'spades')
  # --opts XXX      Extra assembler options in quotes eg. spades: "--untrusted-contigs locus.fna" ... (default: '')
  # --kmers XXX     K-mers to use <blank=AUTO> (default: '')
# MODULES
  # --trim          Enable adaptor trimming (default: OFF)
  # --noreadcorr    Disable read error correction (default: OFF)
  # --nostitch      Disable read stitching (default: OFF)
  # --nocorr        Disable post-assembly correction (default: OFF)
  
#--trim option:
#[trimmomatic]  -threads 8 -phred33 R1.sub.fq.gz R2.sub.fq.gz R1.fq.gz /dev/null R2.fq.gz /dev/null ILLUMINACLIP:/home/mdiricks/bin/miniconda2/db/trimmomatic.fa:1:30:11 LEADING:3 TRAILING:3 MINLEN:30 TOPHRED33

# This will perform the following in this order
# Remove Illumina adapters provided in the trimmomatic.fa file (provided). Initially
# Trimmomatic will look for seed matches (16 bases) allowing maximally 1
# mismatch. These seeds will be extended and clipped if in the case of paired end
# reads a score of 30 is reached (about 50 bases), or in the case of single ended reads a
# score of 11, (about 17 bases).
# Remove leading low quality or N bases (below quality 3)
# Remove trailing low quality or N bases (below quality 3)
# Drop reads which are less than 30 bases long after these steps
