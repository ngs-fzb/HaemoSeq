#!/usr/bin/env bash

###Usage###
# bash $PathToScript/Scriptname.sh -h

###Author###
#Margo Diricks (mdiricks@fz-borstel.de)!"

###Version###
version="1.0.1"

###Function###
#This script performs a quality control of the raw reads using fastQC and summarizes the results with multiqc

###Required packages###
#fastqc and multiqc, installed in a conda environment called multiqc
	
###Required parameters - via command line###
#-i PATH_input=""
#-o PATH_output=""

###Optional parameters that can be changed###
cpu=$(lscpu | grep -E '^CPU\(' | awk '{print $(NF)}')
filetype="fastq.gz"

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: NameScript.sh [parameters]"
   echo "Required parameters:"
   echo "-i     Full path to folder where input files are stored"
   echo "-o     Full path to folder where result files need to be stored"
   echo ""
   echo "Optional parameters":
   echo "-c     amount of cpus that need to be used (depends on your hardware configuration); Default: 1"
   echo "-e     Name of conda environment; Default: mashtree"
   echo ""
   echo "-v     Display version"
   echo "-h     Display help"
}
############################################################
# Get parameters                                                #
############################################################

while getopts ":hi:o:c:d:r:e::v" option; do #:h does not need argument, f: does need argument
   case $option in
      h) # display Help
         Help
         exit;;
      i) #
         PATH_input=$OPTARG;;
      o) # 
         PATH_output_tmp=$OPTARG
         PATH_output=$PATH_output_tmp/MultiQC ;;
      c) # 
         cpu=$OPTARG;;
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
if [[ -z "$PATH_input" ]] || [[ -z "$PATH_output" ]]
then
	echo "Please provide all required arguments (-i PATH_input and -o PATH_output)! use starter_multiqc.sh -h for help on syntax"
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
echo "fastqc version: "$(fastqc -v) >> $PATH_output/info.txt
echo "multiqc version: "$(multiqc --version) >> $PATH_output/info.txt
echo "Input files: "$PATH_input >> $PATH_output/info.txt
echo "Output files:"$PATH_output >> $PATH_output/info.txt
echo "Amount of threads used: "$cpu >> $PATH_output/info.txt
echo "Conda environment: "$conda_env >> $PATH_output/info.txt

###Perform analysis###

for fastq in $PATH_input/*.$filetype
do
	SampleName=$(basename $fastq .$filetype)
	Folder=$PATH_output/$SampleName"_fastqc.zip"
	echo $Folder
	if [ -f $Folder ]
	then
		 echo $SampleName "was already analyzed" #>> $PATH_output/AlreadyAnalyzed.txt 
	else
		fastqc -o $PATH_output --extract -t $cpu $fastq
	fi
done

multiqc --export -o $PATH_output $PATH_output
mv multiqc_data/ $PATH_output

###Clean data###
rm -r $PATH_output/*_fastqc
rm $PATH_output/*_fastqc.html

date >> $PATH_output/info.txt


################################################################################################################################################################
