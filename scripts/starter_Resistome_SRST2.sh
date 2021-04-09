#!/usr/bin/env bash

###Author###
echo "hello World, this script was written by Margo Diricks!"

###Function###
echo "Running this script will reveal the presence of resistance genes included in the resistance database of your choice"
Analysis="Resistance"


###Required packages###
# SRST2 and GNU parallel --> these can be installed e.g. using conda.

###Usage###
# bash $HOME/Scripts/Name_Script.sh"

###Input files###
#Directory where you have stored the (raw) paired-end illumina fastQ files
PATH_input="$HOME/YourFolder"

#Read notation and extension
FastqType=".fastq.gz" # Extension of your fastq file
Fw="_R1" # Forward read notation
Rv="_R2" # Reverse read notation

#Path to nucleotide reference database
geneDB="$HOME/YourFolder/ResistanceDatabase.fasta"

###Parameters###
# Additional parameters can be added if required (see below)
#1. Maximum number of mismatches allowed during mapping - M
GeneMaxMis=50
#2. Maximum %divergence cutoff for gene reporting (in genes report) - D
MaxDiv=15
#3. Minimum %coverage cutoff for gene reporting (in full genes and genes report) - C
MinCov=90
#Note: default values are 50M-15D-90C. Severely truncated genes or homologues with not enough similarity will not be detected using these values!

###Output files###
#1. Directory where the final files will be stored
PATH_output="$HOME/YourFolder/"$Analysis"/Mismatch"$GeneMaxMis"_maxdiv"$MaxDiv"_mincov"$MinCov
#2. Directory where the temporary files will be stored
PATH_temp="$HOME/tmp"

#Amount of samples to be analyzed in parallel 
Jobs=10 #Jobs=0 means you will run as many jobs in parallel as possible. Please note that this can sometimes lead to problems. Preferably use the amount of CPUs/threads you have available.
#Important note: if it is the first time you use a reference database (fasta file defined as geneDB), you need to set this value once at 1 to generate the fai and bt2 files of the reference database. Afterwards, you can set the amount of jobs depending on your computational resources. Please consult the GNU parallel manual for more info. 

################################################################################################################################################################

#Resistome analysis in parallel
source `which env_parallel.bash`
env_parallel -j $Jobs --compress --tmpdir $PATH_temp 'srst2 --input_pe {} $(echo {} | sed "s/$Fw$FastqType$/$Rv$FastqType/1") --forward $Fw --reverse $Rv --log --gene_db $geneDB --gene_max_mismatch $GeneMaxMis --max_divergence $MaxDiv --min_coverage $MinCov --output $PATH_output/$(echo {/}| cut -d '_' -f 1) --use_existing_scores --use_existing_pileup' ::: $PATH_input/*$Fw$FastqType

#Remove pileup and bam files, only keep result files
rm $PATH_output/*.pileup
rm $PATH_output/*.bam
rm $PATH_output/*.sam


#Create a combined file
srst2 --prev_output $PATH_output/*__genes__*__results.txt --output $PATH_output/$Analysis

################################################################################################################################################################
# SRST2 - Short Read Sequence Typer (v2)

# optional arguments:
  # -h, --help            show this help message and exit
  # --version             show program's version number and exit
  # --input_se INPUT_SE [INPUT_SE ...]
                        # Single end read file(s) for analysing (may be gzipped)
  # --input_pe INPUT_PE [INPUT_PE ...]
                        # Paired end read files for analysing (may be gzipped)
  # --merge_paired        Switch on if all the input read sets belong to a
                        # single sample, and you want to merge their data to get
                        # a single result
  # --forward FORWARD     Designator for forward reads (only used if NOT in
                        # MiSeq format sample_S1_L001_R1_001.fastq.gz; otherwise
                        # default is _1, i.e. expect forward reads as
                        # sample_1.fastq.gz)
  # --reverse REVERSE     Designator for reverse reads (only used if NOT in
                        # MiSeq format sample_S1_L001_R2_001.fastq.gz; otherwise
                        # default is _2, i.e. expect forward reads as
                        # sample_2.fastq.gz
  # --read_type {q,qseq,f}
                        # Read file type (for bowtie2; default is q=fastq; other
                        # options: qseq=solexa, f=fasta).
  # --mlst_db MLST_DB     Fasta file of MLST alleles (optional)
  # --mlst_delimiter MLST_DELIMITER
                        # Character(s) separating gene name from allele number
                        # in MLST database (default "-", as in arcc-1)
  # --mlst_definitions MLST_DEFINITIONS
                        # ST definitions for MLST scheme (required if mlst_db
                        # supplied and you want to calculate STs)
  # --mlst_max_mismatch MLST_MAX_MISMATCH
                        # Maximum number of mismatches per read for MLST allele
                        # calling (default 10)
  # --gene_db GENE_DB [GENE_DB ...]
                        # Fasta file/s for gene databases (optional)
  # --no_gene_details     Switch OFF verbose reporting of gene typing
  # --gene_max_mismatch GENE_MAX_MISMATCH
                        # Maximum number of mismatches per read for gene
                        # detection and allele calling (default 10)
  # --min_coverage MIN_COVERAGE
                        # Minimum %coverage cutoff for gene reporting (default
                        # 90)
  # --max_divergence MAX_DIVERGENCE
                        # Maximum %divergence cutoff for gene reporting (default
                        # 10)
  # --min_depth MIN_DEPTH
                        # Minimum mean depth to flag as dubious allele call
                        # (default 5)
  # --min_edge_depth MIN_EDGE_DEPTH
                        # Minimum edge depth to flag as dubious allele call
                        # (default 2)
  # --prob_err PROB_ERR   Probability of sequencing error (default 0.01)
  # --truncation_score_tolerance TRUNCATION_SCORE_TOLERANCE
                        # % increase in score allowed to choose non-truncated
                        # allele
  # --stop_after STOP_AFTER
                        # Stop mapping after this number of reads have been
                        # mapped (otherwise map all). parameter to pass to bowtie2 parameter -u N to stop mapping after the first N reads. Default behaviour remains to map all reads. However, for large read sets (e.g. >100x), extra reads do not help and merely increase the time taken for mapping and scoring, and you may want to limit to the first million reads or read pairs (100x of a 2 Mbp genome (with 100bp PE reads?)) using --stop_after 1000000.
  # --other OTHER         Other arguments to pass to bowtie2 (must be escaped,
                        # e.g. "\--no-mixed".
  # --max_unaligned_overlap MAX_UNALIGNED_OVERLAP
                        # Read discarded from alignment if either of its ends
                        # has unaligned overlap with the reference that is
                        # longer than this value (default 10)
  # --mapq MAPQ           Samtools -q parameter (default 1)
  # --baseq BASEQ         Samtools -Q parameter (default 20)
  # --samtools_args SAMTOOLS_ARGS
                        # Other arguments to pass to samtools mpileup (must be
                        # escaped, e.g. "\-A").
  # --output OUTPUT       Prefix for srst2 output files
  # --log                 Switch ON logging to file (otherwise log to stdout)
  # --save_scores         Switch ON verbose reporting of all scores
  # --report_new_consensus
                        # If a matching alleles is not found, report the
                        # consensus allele. Note, only SNP differences are
                        # considered, not indels.
  # --report_all_consensus
                        # Report the consensus allele for the most likely
                        # allele. Note, only SNP differences are considered, not
                        # indels.
  # --use_existing_bowtie2_sam
                        # Use existing SAM file generated by Bowtie2 if
                        # available, otherwise they will be generated
  # --use_existing_pileup
                        # Use existing pileups if available, otherwise they will
                        # be generated
  # --use_existing_scores
                        # Use existing scores files if available, otherwise they
                        # will be generated
  # --keep_interim_alignment
                        # Keep interim files (sam & unsorted bam), otherwise
                        # they will be deleted after sorted bam is created
  # --threads THREADS     Use multiple threads in Bowtie and Samtools
  # --prev_output PREV_OUTPUT [PREV_OUTPUT ...]
                        # SRST2 results files to compile (any new results from
                        # this run will also be incorporated)
