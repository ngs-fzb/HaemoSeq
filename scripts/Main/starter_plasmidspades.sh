#! /bin/bash

###Usage###
# bash $PathToScript/Scriptname.sh -h

###Author###
#Margo Diricks (mdiricks@fz-borstel.de)!"

###Version###
version="1.0.0"

###Function###
#echo "Running this script will assemble plasmids (and other mobile elements with coverage differing from the genome) using plasmidspades"

###Required packages###
# spades, installed in a conda environment called spades

###Required parameters - via command line###
#-i PATH_fastQ=""
#-o PATH_output=""

###Optional parameters that can be changed###
Fw="_R1" #Change if your fastQ files are not SampleName_R1.fastq.gz
Rv="_R2" #Change if your fastQ files are not SampleName_R2.fastq.gz
filetype=".fastq.gz"
cpu=$(lscpu | grep -E '^CPU\(' | awk '{print $(NF)}') #Uses all available threads
set="SampleSet"
conda_env="spades"
PATH_tmp="$HOME/tmp"

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
   echo ""
   echo "Optional parameters":
   echo "-c     amount of cpus that need to be used (depends on your hardware configuration); Default: All"
   echo "-e     Name of conda environment; Default: spades"
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

while getopts ":hi:o:c:f:r:s:t:e::v" option; do #:h does not need argument, f: does need argument
   case $option in
      h) # display Help
         Help
         exit;;
      i) #
         PATH_fastQ=$OPTARG;;
      o) # 
         PATH_output_tmp=$OPTARG
         PATH_output=$PATH_output_tmp/Plasmidspades ;;
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
if [[ -z "$PATH_fastQ" ]] || [[ -z "$PATH_output" ]]
then
	echo "Please provide all required arguments (-i PATH_fastQ and -o PATH_output)! use starter_plasmidspades.sh -h for help on syntax"
	exit
fi

###############################################################################CODE#################################################################################

###Create folders###
mkdir -p $PATH_output/FinalAssemblies
mkdir -p $PATH_tmp
mkdir $PATH_output/Logfiles

###Activate conda environment###

eval "$(conda shell.bash hook)"
conda activate $conda_env

###Create info file###

date > $PATH_output/info.txt
echo "Version script: "$version >> $PATH_output/info.txt
echo "spades version: "$(spades.py --version) >> $PATH_output/info.txt #Does not work, outputs in command line
echo "Input files: "$PATH_fastQ >> $PATH_output/info.txt
echo "Temporary directory: "$PATH_tmp >> $PATH_output/info.txt
echo "Amount of threads used: "$cpu >> $PATH_output/info.txt
echo "Sample set: "$set >> $PATH_output/info.txt
echo "Conda environment: "$conda_env >> $PATH_output/info.txt
echo "Amount of samples in input folder:" $(ls $PATH_fastQ/*$Fw$filetype | wc -l) >> $PATH_output/info.txt

###Actual command###
for fastq in $PATH_fastQ/*$Fw$filetype
do
	SampleName=$(basename $fastq| cut -d '_' -f 1)
	if [ -f "$PATH_output/Logfiles/$SampleName_spades.log" ]
	then
		 echo $SampleName >> $PATH_output/FinalAssemblies/AlreadyAnalyzed.txt
	else
		spades.py --plasmid --tmp-dir $PATH_tmp -o $PATH_output/$SampleName -1 $fastq -2 $(echo $fastq | sed "s/$Fw$filetype$/$Rv$filetype/1") -t $cpu
		mv $PATH_output/$SampleName/scaffolds.fasta $PATH_output/FinalAssemblies/$SampleName"_pl.fasta"
		mv $PATH_output/$SampleName/spades.log $PATH_output/Logfiles/$SampleName"_spades.log"
		#Clean up, otherwise folders can get big!
		rm -r $PATH_output/$SampleName
	fi
done

###Naming output files###
Name_itolFile_Summary_plasmids=$set"_plasmidspades_plasmid_presence_itol.txt" #This file can be dragged and dropped over an itol tree to visualize the rpesence of plasmids

###Make summary and itol files###
echo "Amount of samples with assembled mobile elements: "$(ls $PATH_output/FinalAssemblies/*_pl.fasta | wc -l) >> $PATH_output/info.txt #contains only hits with coverage >min coverage parameter set. However, it reports also hits with max_divergence > max_divergence parameter set
echo -e "Sample\t#Putative Mobile Elements\t Mobile Elements (node_length-bp_cov)" > $PATH_output/$set"_plasmidspades_summary.txt"

###Make an Itol file for plasmid presence###
echo "DATASET_COLORSTRIP" > $PATH_output/$Name_itolFile_Summary_plasmids
echo "SEPARATOR TAB" >> $PATH_output/$Name_itolFile_Summary_plasmids
echo "DATASET_SCALE	0" >> $PATH_output/$Name_itolFile_Summary_plasmids
#Label is used in the legend table (can be changed later)
echo "DATASET_LABEL	plasmidspades_plasmid_presence" >> $PATH_output/$Name_itolFile_Summary_plasmids
#Dataset color (can be changed later)
echo "COLOR	#ff0000" >> $PATH_output/$Name_itolFile_Summary_plasmids
#Legend
echo "LEGEND_TITLE	plasmidspades_plasmid_presence">> $PATH_output/$Name_itolFile_Summary_plasmids
echo "LEGEND_SHAPES	1" >> $PATH_output/$Name_itolFile_Summary_plasmids
echo "LEGEND_COLORS	#dbd7d2" >> $PATH_output/$Name_itolFile_Summary_plasmids
echo "LEGEND_LABELS	at least 1 plasmid present" >> $PATH_output/$Name_itolFile_Summary_plasmids
#Data
echo "DATA" >> $PATH_output/$Name_itolFile_Summary_plasmids

for fasta in $PATH_output/FinalAssemblies/*_pl.fasta
do
	SampleName=$(basename $fasta| cut -d '_' -f 1)
	MobelCount=$(grep -c ">" $fasta)
	MobilElements=""
	paste --delimiter='\t'  <(echo $SampleName) <(echo "#dbd7d2") >> $PATH_output/$Name_itolFile_Summary_plasmids
	while read line
	do
		if [[ "$line" == *">"* ]]
		then
			Node=$(echo $line | cut -d '_' -f 2)
			Length=$(echo $line | cut -d '_' -f 4)
			Cov=$(echo $line | cut -d '_' -f 6)
			MobilElements+=$Node"_"$Length"_"$Cov";"
		fi
	done <<< "$(cat $fasta)"
	paste --delimiter='\t'  <(echo $SampleName)	<(echo $MobelCount)	<(echo $MobilElements) >> $PATH_output/$set"_plasmidspades_summary.txt"
done

###Closing###
conda deactivate
echo "Script Finished!"
exit 

####################################################################CODE THAT MIGHT BE USED IN ADDITION#################################################################################




#####################################################################HELP##########################################################################################

