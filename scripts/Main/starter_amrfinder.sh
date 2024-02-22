#!/usr/bin/env bash

###Usage###
# bash $PathToScript/Scriptname.sh -h

###Author###
#Margo Diricks (mdiricks@fz-borstel.de)!"

###Version###
version="1.0.0"

###Function###
#This script runs NCBI`s AMRFinderPlus to detect resistance, virulence factors and stress-response genes (https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus)

###Required packages###
#AMRFinder, installed in a conda environment called AMRFinder
	
###Required parameters - via command line###
#-i PATH_input=""
#-o PATH_output=""

###Optional parameters that can be changed###

set="SampleSet"
cpu=$(lscpu | grep -E '^CPU\(' | awk '{print $(NF)}') #Uses all available threads
conda_env="AMRFinder"
ident_min=0.5 #<0-1; Default: 0.9> Minimum identity for a blast-based hit (Methods BLAST or PARTIAL). -1 means use a curated threshold if it exists and 0.9 otherwise. Setting this value to something other than -1 will override any curated similarity cutoffs.
coverage_min=0.5 #<0-1; Default: 0.5> Minimum proportion of reference gene covered for a BLAST-based hit (Methods BLAST or PARTIAL).



############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: NameScript.sh [parameters]"
   echo "Required parameters:"
   echo "-i     Full path to folder where input files (fastA) are stored"
   echo "-o     Full path to folder where result files need to be stored"
   echo ""
   echo "Optional parameters":
   echo "-c     amount of cpus that need to be used (depends on your hardware configuration); Default: 1"
   echo "-e     Name of conda environment; Default: AMRFinder"
   echo "-d     Database with AMRfinder genes"
   echo ""
   echo "-v     Display version"
   echo "-h     Display help"
}
############################################################
# Get parameters                                                #
############################################################

while getopts ":hi:o:c:d:e::v" option; do #:h does not need argument, f: does need argument
   case $option in
      h) # display Help
         Help
         exit;;
      i) #
         PATH_input=$OPTARG;;
      o) # 
         PATH_output_tmp=$OPTARG
		 PATH_output=$PATH_output_tmp/AMRFinder ;;
      c) # 
         cpu=$OPTARG;;
	  d) #
         db_AMRfinder=$OPTARG;;
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
if [[ -z "$PATH_input" ]] || [[ -z "$PATH_output" ]] || [[ -z "$db_AMRfinder" ]]
then
	echo "Please provide all required arguments (-i PATH_input, -o PATH_output, -d db_AMRfinder)! use starter_amrfinder.sh -h for help on syntax"
	exit
fi

###############################################################################CODE#################################################################################

###Create folders###
mkdir -p $PATH_output

###Activate conda environment###

eval "$(conda shell.bash hook)"
conda activate $conda_env

###Create info file###

date > $PATH_output/info.txt
echo "Version script: "$version >> $PATH_output/info.txt
echo "AMRFinder version: "$(amrfinder --database_version) >> $PATH_output/info.txt 
echo "AMRfinder database version: "$db_AMRfinder >> $PATH_output/info.txt 
echo "Minimum identity to detect a gene: " $ident_min >> $PATH_output/info.txt 
echo "Minimum coverage to detect a gene: " $coverage_min >> $PATH_output/info.txt 
echo "Input files: "$PATH_input >> $PATH_output/info.txt
echo "Amount of threads used: "$cpu >> $PATH_output/info.txt
echo "Sample set: "$set >> $PATH_output/info.txt
echo "Conda environment: "$conda_env >> $PATH_output/info.txt
echo "Amount of samples in input folder:" $(ls $PATH_input/*.fasta | wc -l) >> $PATH_output/info.txt

###Actual command###
#Download latest AMRFinder plus database (about 0.2 Gb per database). By default, the database is stored in the conda environment (e.g. /molmyc/miniconda3/envs/AMRFinder/share/amrfinderplus/data/2023-07-13.2). if updated, the older database is not removed!! Remove manual if you want to save space. See info file for the location of database)
#amrfinder -u #To update database in default directory
#amrfinder_update -d <database_directory> #To update database in user-specified directory

for fasta in $PATH_input/*.fasta
do
	SampleName=$(basename $fasta .fasta)
	File=$PATH_output/$SampleName".txt"
	if [ -f "$File" ]
	then
		 echo $SampleName "was already analyzed" #>> $PATH_output/AlreadyAnalyzed.txt 
	else
		amrfinder -n $fasta --threads $cpu --name $SampleName -o $PATH_output/$SampleName".txt" --plus -i $ident_min -c $coverage_min -d $db_AMRfinder
	fi
done

###Making summary file###
rm $PATH_output/AMRfinder_summary.txt
echo -e "Name\tProtein identifier\tContig id\tStart\tStop\tStrand\tGene symbol\tSequence name\tScope\tElement type\tElement subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference sequence\t% Identity to reference sequence\tAlignment length\tAccession of closest sequence\tName of closest sequence\tHMM id\tHMM description" >> $PATH_output/AMRfinder_summary.txt
for file in $PATH_output/*.txt #don´t include info file in summary
do 
if [[ "$file" != "$PATH_output/AMRfinder_summary.txt" ]] && [[ "$file" != "$PATH_output/info.txt" ]] #Very important that you don´t include the summary file, otherwise you get TB files!!
then
	awk 'NR>1' $file >> $PATH_output/AMRfinder_summary.txt
fi
done

###Closing###
date >> $PATH_output/info.txt
conda deactivate
echo "Script Finished!"
exit 
################################################################################################################################################################
