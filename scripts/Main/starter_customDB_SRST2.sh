#!/usr/bin/env bash

###Usage###
# bash $PathToScript/Scriptname.sh -h

###Author###
#Margo Diricks (mdiricks@fz-borstel.de)!"

###Version###
version="1.0.0"

###Function###
#echo "Running this script will reveal the presence of known plasmids from the PLSDB database"
Analysis="customDB_SRST2"

###Required packages###
# SRST2 and GNU parallel, installed in a conda environment called SRST2

###Required parameters - via command line###
#-i PATH_fastQ=""
#-o PATH_output=""
#-d db_PLSDB_SRST2=""

###Optional parameters that can be changed###
Fw="_R1" #Change if your fastQ files are not SampleName_R1.fastq.gz
Rv="_R2" #Change if your fastQ files are not SampleName_R2.fastq.gz
filetype=".fastq.gz"
cpu=$(lscpu | grep -E '^CPU\(' | awk '{print $(NF)}') #Uses all available threads
set="SampleSet"
conda_env="SRST2"
PATH_tmp="$HOME/tmp"

GeneMaxMis=10 #M - Maximum number of mismatches allowed during mapping (default 10)
MaxDiv=10 #D - Maximum %divergence cutoff for gene reporting (in genes report) (default 10)
MinCov=90 #C - Minimum %coverage cutoff for gene reporting (in full genes and genes report) (default 90)


############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: NameScript.sh [parameters]"
   echo "Required parameters:"
   echo "-i     Full path to folder where fastq files are stored"
   echo "-o     Full path to folder where result files need to be stored"
   echo "-d     Full path to plasmid database ($PATH/db_name.fasta)"
   echo ""
   echo "Optional parameters":
   echo "-c     amount of cpus that need to be used (depends on your hardware configuration); Default: All"
   echo "-e     Name of conda environment; Default: SRST2"
   echo "-f     Forward read notation; Default: _R1"
   echo "-r     Reverse read notation; Default: _R2"
   echo "-s     Name of sample set - used for file naming; Default: Sampleset"
   echo "-t     Full path to directory where temporary files will be stored; Default: $HOME/tmp"
   echo ""
   echo "-v     Display version"
   echo "-h     Display help"
}
############################################################
# Get parameters                                                #
############################################################

while getopts ":hi:o:c:d:f:r:s:t:e::v" option; do #:h does not need argument, f: does need argument
   case $option in
      h) # display Help
         Help
         exit;;
      i) #
         PATH_fastQ=$OPTARG;;
      o) # 
         PATH_output_tmp=$OPTARG
         PATH_output=$PATH_output_tmp/customDB_SRST2 ;;
      d) # 
         db_custom_SRST2=$OPTARG;;
      c) # 
         cpu=$OPTARG;;
      f) # 
         Fw=$OPTARG;;
      r) # 
         Rv=$OPTARG;;
      s) # 
         set=$OPTARG;;
      t) # 
         PATH_tmp=$OPTARG;;
      e) # 
         conda_env=$OPTARG;;
      v) # display Version
         echo $version
         exit;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

###Check if required parameters are provided###
if [[ -z "$PATH_fastQ" ]] || [[ -z "$PATH_output" ]] || [[ -z "$db_custom_SRST2" ]]
then
	echo "Please provide all required arguments (-i PATH_fastQ, -o PATH_output and -d db_custom_SRST2)! use starter_custom_SRST2.sh -h for help on syntax"
	exit
fi

###############################################################################CODE#################################################################################
#Remove previous files
rm $PATH_output/failed.txt

###Create folders###
mkdir -p $PATH_output
mkdir -p $PATH_tmp

###Activate conda environment###

eval "$(conda shell.bash hook)"
conda activate $conda_env

###Create info file###

date > $PATH_output/info.txt
echo "Version script: "$version >> $PATH_output/info.txt
echo "Plasmid database used: "$db_custom_SRST2 >> $PATH_output/info.txt
#echo "SRST2 version: "$(srst2 --version) >> $PATH_output/info.txt #Does not work, outputs in command line
echo "M - Maximum number of mismatches allowed during mapping: "$GeneMaxMis >> $PATH_output/info.txt
echo "D - Maximum %divergence cutoff for gene reporting (in genes report): "$MaxDiv  >> $PATH_output/info.txt
echo "C - Minimum %coverage cutoff for gene reporting (in full genes and genes report): "$MinCov >> $PATH_output/info.txt
echo "Input files: "$PATH_fastQ >> $PATH_output/info.txt
echo "Output files:"$PATH_output >> $PATH_output/info.txt
echo "Temporary directory: "$PATH_tmp >> $PATH_output/info.txt
echo "Amount of threads used: "$cpu >> $PATH_output/info.txt
echo "Sample set: "$set >> $PATH_output/info.txt
echo "Conda environment: "$conda_env >> $PATH_output/info.txt
echo "Amount of samples in input folder:" $(ls $PATH_fastQ/*$Fw$filetype | wc -l) >> $PATH_output/info.txt

###Build index files for SRST2 database (if they don´t exist yet)###
[[ ! -s $db_custom_SRST2".rev.1.bt2" ]] && bowtie2-build $db_custom_SRST2 $db_custom_SRST2; samtools faidx $db_custom_SRST2

#Check if sample was already analysed
ls $PATH_output/*.log > $PATH_output/LogAvailable.txt
ls $PATH_fastQ/*$Fw$filetype > $PATH_output/ToAnalyse.txt

group_toDownload=()

while read line
do
	SampleName=$(basename $line| cut -d '_' -f 1)
	if grep -q "$SampleName" $PATH_output/LogAvailable.txt
	then 
		echo "Sample is already analysed"
	else
		group_toDownload+=( $line )
fi
done <<< "$(cat $PATH_output/ToAnalyse.txt)"

###Plasmid identification###
if [ ${#group_toDownload[@]} -eq 0 ]; then
	echo "Everything is analysed!"
else
source `which env_parallel.bash`
env_parallel -j $cpu --compress --tmpdir $PATH_tmp 'srst2 --input_pe {} $(echo {} | sed "s/$Fw/$Rv/1") --forward $Fw --reverse $Rv --log --gene_db $db_custom_SRST2 --gene_max_mismatch $GeneMaxMis --max_divergence $MaxDiv --min_coverage $MinCov --output $PATH_output/$(echo {/}| cut -d '_' -f 1) --use_existing_scores  --use_existing_bowtie2_sam --use_existing_pileup' ::: "${group_toDownload[@]}"
fi
###Create a combined file###
srst2 --prev_output $PATH_output/*__genes__*__results.txt --output $PATH_output/$Analysis"_"$set


#Check if for each sample if the algorithm finished
COUNTER_f=0
for log in $PATH_output/*.log
do
	SampleName=$(basename $log| cut -d '.' -f 1)
	if grep -q "SRST2 has finished" $log
	then 
		echo "Algorithm finished succesfully for "$SampleName
	else
		echo $SampleName >> $PATH_output/Failed.txt
		COUNTER_f=$((COUNTER_f + 1))
	fi
done

echo "Amount of samples with full genes result file: "$(ls $PATH_output/*__fullgenes__* | wc -l) >> $PATH_output/info.txt #contains only hits with coverage >min coverage parameter set. However, it reports also hits with max_divergence > max_divergence parameter set
echo "Amount of samples with genes result file: "$(ls $PATH_output/*__genes__* | wc -l) >> $PATH_output/info.txt #only contain hits complying to the thresholds 
echo "Amount of samples with (intermediate) pile-up result file: "$(ls $PATH_output/*.pileup | wc -l) >> $PATH_output/info.txt
echo "Amount of samples for which algorithm failed: "$COUNTER_f >> $PATH_output/info.txt

###Create a list of plasmids for which valid hits were found###
head -n 1 $PATH_output/$Analysis"_"$set"__compiledResults.txt" | cut -f 2- | tr '\t' '\n' > $PATH_output/Plasmid_hits.txt

###Naming output files###
Name_itolFile_Summary_plasmids=$set"_custom-SRST2_plasmid_presence_itol.txt" #This file can be dragged and dropped over an itol tree to visualize the rpesence of plasmids

###Make an Itol file for plasmid presence###
echo "DATASET_COLORSTRIP" > $PATH_output/$Name_itolFile_Summary_plasmids
echo "SEPARATOR TAB" >> $PATH_output/$Name_itolFile_Summary_plasmids
echo "DATASET_SCALE	0" >> $PATH_output/$Name_itolFile_Summary_plasmids
#Label is used in the legend table (can be changed later)
echo "DATASET_LABEL	custom-SRST2_plasmid_presence" >> $PATH_output/$Name_itolFile_Summary_plasmids
#Dataset color (can be changed later)
echo "COLOR	#ff0000" >> $PATH_output/$Name_itolFile_Summary_plasmids
#Legend
echo "LEGEND_TITLE	custom-SRST2_plasmid_presence">> $PATH_output/$Name_itolFile_Summary_plasmids
echo "LEGEND_SHAPES	1" >> $PATH_output/$Name_itolFile_Summary_plasmids
echo "LEGEND_COLORS	#808080" >> $PATH_output/$Name_itolFile_Summary_plasmids
echo "LEGEND_LABELS	at least 1 plasmid present" >> $PATH_output/$Name_itolFile_Summary_plasmids
#Data
echo "DATA" >> $PATH_output/$Name_itolFile_Summary_plasmids

COUNTER=0
while read line
do
	COUNTER=$((COUNTER + 1))
	SampleName=$(echo $line | cut -d ' ' -f 1 | cut -d '_' -f 1)
	if [ "$COUNTER" -gt 1 ]
	then
		paste --delimiter='\t'  <(echo $SampleName) <(echo "#808080") >> $PATH_output/$Name_itolFile_Summary_plasmids
	fi
done <<< "$(cat $PATH_output/$Analysis"_"$set"__compiledResults.txt")"


###Closing###
conda deactivate
echo "Script Finished!" >> $PATH_output/info.txt
date >> $PATH_output/info.txt


####################################################################CODE THAT MIGHT BE USED IN ADDITION#################################################################################


###Clean up###
rm $PATH_output/*.pileup
rm $PATH_output/*.bam
rm $PATH_output/*.sam

###Download plasmids for which a hit was found###
#bash starter_download_plasmids.sh 

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
						
exit 