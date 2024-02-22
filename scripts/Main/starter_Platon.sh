#!/usr/bin/env bash

###Usage###
# bash $PathToScript/Scriptname.sh -h

###Author###
#echo "hello World, this script was written by Margo Diricks (mdiricks@fz-borstel.de)!"

###Version###
version="1.0.0"

###Function###
#echo "Running this script will predict the presence of plasmids from draft assemblies using platon"

###Required packages###
# platon, installed in a conda environment called platon

###Required parameters - via command line###
#-i PATH_fastA=""
#-o PATH_output=""
#-d db_platon=""

###Optional parameters that can be changed###
cpu=$(lscpu | grep -E '^CPU\(' | awk '{print $(NF)}') #Uses all available threads
set="SampleSet"
conda_env="platon"


############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: NameScript.sh [parameters]"
   echo "Required parameters:"
   echo "-i     Full path to folder where fastA files (draft assemblies) are stored"
   echo "-o     Full path to folder where result files need to be stored"
   echo "-d     Full path to platon database folder ($PATH/db)"
   echo ""
   echo "Optional parameters":
   echo "-c     amount of cpus that need to be used (depends on your hardware configuration); Default: All"
   echo "-e     Name of conda environment; Default: platon"
   echo "-s     Name of sample set - used for file naming; Default: Sampleset"
   echo ""
   echo "-v     Display version"
   echo "-h     Display help"
}
############################################################
# Get parameters                                                #
############################################################

while getopts ":hi:o:c:d:s:e::v" option; do 
   case $option in
      h) # display Help
         Help
         exit;;
      i) #
         PATH_fastA=$OPTARG;;
      o) # 
         PATH_output_tmp=$OPTARG
         PATH_output=$PATH_output_tmp/Platon ;;
      d) # 
         db_platon=$OPTARG;;
      c) # 
         cpu=$OPTARG;;
      s) # 
         set=$OPTARG;;
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
if [[ -z "$PATH_fastA" ]] || [[ -z "$PATH_output" ]] || [[ -z "$db_platon" ]]
then
	echo "Please provide all required arguments (-i PATH_fastA, -o PATH_output and -d db_platon)! use starter_platon.sh -h for help on syntax"
	exit
fi

###############################################################################CODE#################################################################################

###Create folders###
mkdir -p $PATH_output/FinalAssemblies
rm $PATH_output/AlreadyAnalyzed.txt
rm $PATH_output/Failed.txt

###Activate conda environment###

eval "$(conda shell.bash hook)"
conda activate $conda_env

###Create info file###

date > $PATH_output/info.txt
echo "Version script: "$version >> $PATH_output/info.txt
echo "Plasmid database used: "$db_platon >> $PATH_output/info.txt
echo "platon version: "$(platon --version) >> $PATH_output/info.txt
echo "Input files: "$PATH_fastA >> $PATH_output/info.txt
echo "Amount of threads used: "$cpu >> $PATH_output/info.txt
echo "Sample set: "$set >> $PATH_output/info.txt
echo "Conda environment: "$conda_env >> $PATH_output/info.txt
echo "Amount of samples in input folder:" $(ls $PATH_fastA/*.fasta | wc -l) >> $PATH_output/info.txt


###Start platon###
for file in $PATH_fastA/*.fasta
do
	SampleName=$(basename $file| cut -d '_' -f 1 | cut -d '.' -f 1)
	if [ -d $PATH_output/$SampleName ]
	then
		echo $SampleName >> $PATH_output/AlreadyAnalyzed.txt
	else
		platon --output $PATH_output/$SampleName --prefix $SampleName --db $db_platon "$file" --threads $cpu
		#Clean up
		rm $PATH_output/$SampleName/*.chromosome.fasta
		rm $PATH_output/$SampleName/$SampleName.json
	fi
done



###Check for which samples algorithm failed###
Counter_f=0
for file in $PATH_output/*/*.log
do
	SampleName=$(basename $file| cut -d '.' -f 1)
	tsv=$(echo $file | cut -d '.' -f 1).tsv
	if [[ ! -s "$tsv" ]]
	then
		echo $SampleName >> $PATH_output/Failed.txt
		COUNTER_f=$((COUNTER_f + 1))
	fi
done

echo "Amount of samples for which the algorithm failed (see Failed.txt):" $COUNTER_f >> $PATH_output/info.txt

###Move all plasmid files to FinalAssemblies
mv $PATH_output/*/*.plasmid.fasta $PATH_output/FinalAssemblies

### Make summary and itol files###
#Naming output files#
Name_itolFile_Summary_plasmids=$set"_platon_plasmid_presence_itol.txt" #This file can be dragged and dropped over an itol tree to visualize the rpesence of plasmids

#Make header of the summary files 
echo -e "Sample\tID\tLength\tCoverage\t# ORFs\tRDS\tCircular\tInc Type(s)\t# Replication\t# Mobilization\t# OriT\t# Conjugation\t# AMRs\t# rRNAs\t# Plasmid Hits" > $PATH_output/"Platon_Output_summary_"$set.tsv
echo -e "Sample\t#putative mobile elements (platon)" > $PATH_output/"Platon_Output_summary2_"$set.tsv

###Make an Itol file for plasmid presence###
echo "DATASET_COLORSTRIP" > $PATH_output/$Name_itolFile_Summary_plasmids
echo "SEPARATOR TAB" >> $PATH_output/$Name_itolFile_Summary_plasmids
echo "DATASET_SCALE	0" >> $PATH_output/$Name_itolFile_Summary_plasmids
#Label is used in the legend table (can be changed later)
echo "DATASET_LABEL	platon_plasmid_presence" >> $PATH_output/$Name_itolFile_Summary_plasmids
#Dataset color (can be changed later)
echo "COLOR	#ff0000" >> $PATH_output/$Name_itolFile_Summary_plasmids
#Legend
echo "LEGEND_TITLE	platon_plasmid_presence">> $PATH_output/$Name_itolFile_Summary_plasmids
echo "LEGEND_SHAPES	1" >> $PATH_output/$Name_itolFile_Summary_plasmids
echo "LEGEND_COLORS	#000000" >> $PATH_output/$Name_itolFile_Summary_plasmids
echo "LEGEND_LABELS	at least 1 plasmid present" >> $PATH_output/$Name_itolFile_Summary_plasmids
#Data
echo "DATA" >> $PATH_output/$Name_itolFile_Summary_plasmids

COUNTER2=0
for file in $PATH_output/*/*.tsv
do
	COUNTER=0
	SampleName=$(basename $file| cut -d '.' -f 1)
	while read line
	do
		COUNTER=$((COUNTER + 1))
		if [ "$COUNTER" -gt 1 ]
		then
			paste --delimiter='\t'  <(echo $SampleName) <(echo $line) >> $PATH_output/"Platon_Output_summary_"$set.tsv
		fi
	done <<< "$(cat $file)"
	
	paste --delimiter='\t'  <(echo $SampleName) <(echo $(( $COUNTER - 1 ))) >> $PATH_output/"Platon_Output_summary2_"$set.tsv
	
	if [ "$COUNTER" -gt 1 ] 
	then
		paste --delimiter='\t'  <(echo $SampleName) <(echo "#000000") >> $PATH_output/$Name_itolFile_Summary_plasmids
		COUNTER2=$((COUNTER2 + 1))
	fi
done

echo "Amount of samples with plasmid borne contigs: "$COUNTER2 >> $PATH_output/info.txt

###Clean up###
#Delete empty plasmid files
find $PATH_output/FinalAssemblies -maxdepth 1 -type f -empty -print -delete

###Closing###
conda deactivate
echo "Script Finished!" >> $PATH_output/info.txt
date >> $PATH_output/info.txt
exit 


####################################################################CODE THAT MIGHT BE USED IN ADDITION#################################################################################

##do it again with all filters deactivated and forcing characterization, if nothing was found but you really want it!
# for file in $PATH_input/*.fasta
# do
# SampleName=$(basename $file| cut -d '_' -f 1)
# platon --characterize --output $PATH_output/Forced --prefix $SampleName --db $db "$file"
# done
# exit

################################################################################################################################################################
#https://github.com/oschwengers/platon
# usage: platon [--db DB] [--prefix PREFIX] [--output OUTPUT]
              # [--mode {sensitivity,accuracy,specificity}] [--characterize]
              # [--help] [--verbose] [--threads THREADS] [--version]
              # <genome>

# Identification and characterization of bacterial plasmid contigs from short-read draft assemblies.

# Input / Output:
  # <genome>              draft genome in fasta format
  # --db DB, -d DB        database path (default = <platon_path>/db)
  # --prefix PREFIX, -p PREFIX
                        # Prefix for output files
  # --output OUTPUT, -o OUTPUT
                        # Output directory (default = current working directory)

# Workflow:
  # --mode {sensitivity,accuracy,specificity}, -m {sensitivity,accuracy,specificity}
                        # applied filter mode: sensitivity: RDS only (>= 95%
                        # sensitivity); specificity: RDS only (>=99.9%
                        # specificity); accuracy: RDS & characterization
                        # heuristics (highest accuracy) (default = accuracy)
  # --characterize, -c    deactivate filters; characterize all contigs

# General:
  # --help, -h            Show this help message and exit
  # --verbose, -v         Print verbose information
  # --threads THREADS, -t THREADS
                        # Number of threads to use (default = number of
                        # available CPUs)
  # --version             show program's version number and exit

# Citation:
# Schwengers O., Barth P., Falgenhauer L., Hain T., Chakraborty T., & Goesmann A. (2020).
# Platon: identification and characterization of bacterial plasmid contigs in short-read draft assemblies exploiting protein sequence-based replicon distribution scores.
# Microbial Genomics, 95, 295. https://doi.org/10.1099/mgen.0.000398

# GitHub:
# https://github.com/oschwengers/platon
