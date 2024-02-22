#!/usr/bin/env bash

###Author###
echo "hello World, this script was written by Margo Diricks!"
echo "Don't forget to change/insert the path to the input files in this script!"
echo "Don't forget to activate the right conda environment!"

###Function###
#This script runs hicap, a software tool that looks for the capsule genes of Haemophilus influenzae strains (in silico serotyping). If not all required genes are present, the strains are regarded as being non-typeable (NTHi)

###Required packages###
#Install hicap and it's dependencies https://github.com/scwatts/hicap
#Parallel if you want to run samples in parallel https://www.gnu.org/software/parallel/

###Usage###
# bash file_name.sh

###Inputfiles###
#Type: Fasta file
#Path to input files:
PATHNAME_inputFolder=""
PATHNAME_outputFolder=""

###EstimatedTime###
#20 sec per sample (2 Mbp genome)

################################################################################################################################################################
mkdir $PATHNAME_outputFolder
#Sequential execution
#mkdir -p $PATHNAME_outputFolder
# for assembly_fp in $PATHNAME_inputFolder/*.fasta
# do
	# hicap --query_fp "${assembly_fp}" --output_dir $PATHNAME_outputFolder;
# done

#Parallel execution; IMPORTANT: this will use all your CPU if you dont specify the amount of jobs with -j!
parallel -j 30 hicap --query_fp {} --output_dir $PATHNAME_outputFolder ::: $PATHNAME_inputFolder/*.fasta

#Summarize results in one file
echo "#isolate	predicted_serotype	attributes	genes_identified	locus_location	region_I_genes	region_II_genes	region_III_genes	IS1016_hits" > $PATHNAME_outputFolder/Summary_Hicap.txt
for result in $PATHNAME_outputFolder/*.tsv
do
tail -n 1 $result >> $PATHNAME_outputFolder/Summary_Hicap.txt
done
