#!/usr/bin/env bash

###Usage###
# bash $PathToScript/Scriptname.sh -h

###Author###
#Margo Diricks (mdiricks@fz-borstel.de)!"

###Version###
version="1.0.0"

###Function###
#echo "This script runs kraken2, which identifieds the taxonomic ID of your raw reads or contigs"

###Required packages###
# kraken2 and krona installed in a conda environment called kraken2 (https://anaconda.org/bioconda/kraken2 and https://anaconda.org/bioconda/krona).
# Download kraken2 database and krona database (see below, section HELP!)

###Required parameters - via command line###
#-i PATH_input=""
#-o PATH_output=""
#-d db=""

###Optional parameters that can be changed###
Fw="_R1" #Change if your fastQ files are not SampleName_R1.fastq.gz
Rv="_R2" #Change if your fastQ files are not SampleName_R2.fastq.gz
filetype=".fastq.gz"
cpu=$(lscpu | grep -E '^CPU\(' | awk '{print $(NF)}') #Uses all available threads
db_krona=""
set="SampleSet"
conda_env="kraken2"
#PATH_tmp="$HOME/tmp"
Analysis="kraken2"
#MinBaseQuality=15 #minimum base quality

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: NameScript.sh [parameters]"
   echo "Required parameters:"
   echo "-i     Full path to folder where reads or contigs are stored"
   echo "-o     Full path to folder where result files need to be stored"
   echo "-d     Full path to kraken2 database"
   echo ""
   echo "Optional parameters":
   echo "-c     amount of cpus that need to be used (depends on your hardware configuration); Default: All"
   echo "-e     Name of conda environment; Default: kraken2"
   echo "-f     Forward read notation; Default: _R1"
   echo "-r     Reverse read notation; Default: _R2"
   echo "-k     Full path to krona database (Default: $PATH_to_conda_installation/envs/conda_environment_name/opt/krona/taxonomy)" 
   echo "-s     Name of sample set - used for file naming; Default: Sampleset"
   echo "-t     Full path to directory where temporary files will be stored; Default: $HOME/tmp"
   echo ""
   echo "-v     Display version"
   echo "-h     Display help"
}
############################################################
# Get parameters                                                #
############################################################

while getopts ":hi:o:d:c:f:r:k:p:s:t:e::v" option; do #:h does not need argument, f: does need argument
   case $option in
      h) # display Help
         Help
         exit;;
      i) #
         PATH_input=$OPTARG;;
      o) # 
         PATH_output_tmp=$OPTARG
         PATH_output=$PATH_output_tmp/$Analysis ;;
      d) # 
         db_kraken2=$OPTARG;;
      c) # 
         cpu=$OPTARG;;
      f) # 
         Fw=$OPTARG;;
      r) # 
         Rv=$OPTARG;;
      k) # 
         db_krona=$OPTARG;;
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
if [[ -z "$PATH_input" ]] || [[ -z "$PATH_output" ]] || [[ -z "$db_kraken2" ]]
then
	echo "Please provide all required arguments (-i PATH_input, -o PATH_output, -d db_kraken2 and -k db_krona)! use starter_kraken2.sh -h for help on syntax"
	exit
fi

###############################################################################CODE#################################################################################
#Remove previous files
#rm $PATH_output/failed.txt

###Create folders###
mkdir -p $PATH_output
#mkdir -p $PATH_tmp

###Activate conda environment###

eval "$(conda shell.bash hook)"
conda activate $conda_env

###Initiate krona database###
if [[ ! -z "$db_krona" ]] # if value for krona was given (i.e. other then default db location)
then
	#echo "Symbolic link to krona database will be made"
	PATH_conda_env=$(conda info --envs | grep "*" | cut -d "*" -f 2 | cut -d " " -f 3)
	#echo $PATH_conda_env
	rm -rf $PATH_conda_env/opt/krona/taxonomy
	ln -s $db_krona $PATH_conda_env/opt/krona/taxonomy
else
	echo "Default location will be used" 
# Download/update krona database if necessary
#ktUpdateTaxonomy.sh
#if error see
#https://github.com/bioconda/bioconda-recipes/issues/10959
fi

###Create info file###

date > $PATH_output/info.txt
echo "Version script: "$version >> $PATH_output/info.txt
echo "kraken2 database used: "$db_kraken2 >> $PATH_output/info.txt
echo "kraken2 version: "$(kraken2 --version) >> $PATH_output/info.txt 
echo "Input files: "$PATH_input >> $PATH_output/info.txt
echo "Output files:"$PATH_output >> $PATH_output/info.txt
echo "Temporary directory: "$PATH_tmp >> $PATH_output/info.txt
echo "Amount of threads used: "$cpu >> $PATH_output/info.txt
echo "Sample set: "$set >> $PATH_output/info.txt
echo "Conda environment: "$conda_env >> $PATH_output/info.txt
if [[ "$filetype" == *"fastq"* ]]
then
	echo "Amount of samples in input folder:" $(ls $PATH_input/*$Fw$filetype | wc -l) >> $PATH_output/info.txt
elif [[ "$filetype" == *"fasta"* ]]
then
	echo "Amount of samples in input folder:" $(ls $PATH_input/*$filetype | wc -l) >> $PATH_output/info.txt
fi

echo "Start species identification with kraken2"
###Loop through all the raw fastQ or contig files###
if [[ "$filetype" == *"fastq"* ]]
then
	for input in $PATH_input/*$Fw$filetype 
	do
		SampleName=$(basename $input | cut -d '_' -f 1)
		if [ -s $PATH_output/$SampleName.report ]
		then
			echo "Sample was already analysed"
		else
			kraken2 --threads $cpu --db $db_kraken2 --output $PATH_output/$SampleName.kraken --report $PATH_output/$SampleName.report --paired $input $(echo $input | sed "s/$Fw$filetype$/$Rv$filetype/1") #--minimum-base-quality $MinBaseQuality
			awk -F '\t' '{print $2"	"$3}' $PATH_output/$SampleName.kraken > $PATH_output/$SampleName.krona
			ktImportTaxonomy -o $PATH_output/$SampleName.krona.html $PATH_output/$SampleName.krona
		fi
	done
elif [[ "$filetype" == *"fasta"* ]]
then
	for input in $PATH_input/*$filetype 
	do
		SampleName=$(basename $input | cut -d '_' -f 1 | cut -d '.' -f 1)
		if [ -s $PATH_output/$SampleName.report ]
		then
			echo "Sample was already analysed"
		else	
			kraken2 --threads $threads --db $db_kraken2 --output $PATH_output/$SampleName.kraken --report $PATH_output/$SampleName.report $fasta
			awk -F '\t' '{print $2"	"$3}' $PATH_output/$SampleName.kraken > $PATH_output/$SampleName.krona
			ktImportTaxonomy -o $PATH_output/$SampleName.krona.html $PATH_output/$SampleName.krona
		fi
	done
fi

echo "Amount of samples with krona.html file: "$(ls $PATH_output/*.krona.html | wc -l) >> $PATH_output/info.txt 

###Clean up###
rm -r $PATH_output/*.krona.html.files
rm $PATH_output/*.kraken # if you want to see which ID is assigned to every single contig or read, don´t do this
rm $PATH_output/*.krona

###Closing###
conda deactivate
echo "Script Finished!" >> $PATH_output/info.txt
date >> $PATH_output/info.txt


####################################################################CODE THAT MIGHT BE USED IN ADDITION#################################################################################

###Download kraken2 database###
#https://github.com/DerrickWood/kraken2/wiki/Manual#installation
#After kraken2 is installed with e.g. conda, you need to download the databsae using the kraken build script: $kraken2-build --standard --db $DBNAME.
#Replace "$DBNAME" above with your preferred database name/location. Please note that the database will use approximately 100 GB of disk space during creation, with the majority of that being reference sequences or taxonomy mapping information that can be removed after the build.)

###Download krona database###
#https://github.com/marbl/Krona/wiki/KronaTools
# After Krona is installed, You still need to manually download/update the taxonomy databases before Krona can generate taxonomic reports.
#Use $ktUpdateTaxonomy.sh to download/install.  The default location for storing taxonomic databases is your conda environment e.g. $PATH_to_conda_installation/envs/conda_environment_name/opt/krona/taxonomy.
#If you get an error: check https://github.com/bioconda/bioconda-recipes/issues/10959
#Use $conda info --envs to find path to your conda environment
# If you would like the taxonomic data stored elsewhere, simply replace this directory with a symlink.  For example:
# $rm -rf $PATH_to_conda_installation/envs/conda_environment_name/opt/krona/taxonomy
# $mkdir /path/on/big/disk/taxonomy
# $ln -s /path/on/big/disk/taxonomy $PATH_to_conda_installation/envs/conda_environment_name/opt/krona/taxonomy
# $ktUpdateTaxonomy.sh
#If you already downloaded the krona database: add the /path/on/big/disk/taxonomy folder as a parameter -k to the script




##################################################################################HELP##################################################################################################
###INSTALLATION###






###Program parameters###

						
exit 