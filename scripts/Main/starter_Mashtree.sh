#!/usr/bin/env bash

###Usage###
# bash $PathToScript/Scriptname.sh -h

###Author###
#Margo Diricks (mdiricks@fz-borstel.de)!"

###Version###
version="1.0.0"

###Function###
#This script calculates mash distances between set of samples and creates a tree out of it

###Required packages###
#mashtree, installed in a conda environment called mashtree
	
###Required parameters - via command line###
#-i PATH_input=""
#-o PATH_output=""
#-n filetype=""

###Optional parameters that can be changed###
cpu=1 #Note that >1 cpu can lead to problems on some devices; All CPUS: $(lscpu | grep -E '^CPU\(' | awk '{print $(NF)}')
set="SampleSet"
conda_env="mashtree"
Bootstrap="No" #Yes or No; Default: No. Note that performing bootstrapping increases calculation time considerably. Bootstrap functionality is not yet implemented in this script
PATH_tmp="" #add --tempdir $PATH_tmp to command if you want to use this!
#Fw="_R1" #Change if your fastQ files are not SampleName_R1.fastq.gz 
MinDepth=5 #Default 5, D; Use 0 if you want to do bootstrap or jacknife
SketchSize=10000 #Default, 10000 S 
kmerLength=21 #Default, 21 k
reps=100 #Default, 100; For confidence values/bootstrapping

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
   echo "-n     Input file type: .fastq, .fastq.gz, .fasta, .gbk or .msh"
   echo ""
   echo "Optional parameters":
   echo "-b     perform bootstrapping; Default: No, because this functionality is not yet implemented!"
   echo "-c     amount of cpus that need to be used (depends on your hardware configuration); Default: 1"
   echo "-e     Name of conda environment; Default: mashtree"
   echo "-s     Name of sample set - used for file naming; Default: Sampleset"
   echo ""
   echo "-v     Display version"
   echo "-h     Display help"
}
############################################################
# Get parameters                                                #
############################################################

while getopts ":hi:o:b:c:d:n:r:s:e::v" option; do #:h does not need argument, f: does need argument
   case $option in
      h) # display Help
         Help
         exit;;
      i) #
         PATH_input=$OPTARG;;
      n) # 
         filetype=$OPTARG;;
      o) # 
         PATH_output_tmp=$OPTARG ;;
      b) # 
         Bootstrap=$OPTARG;;
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
if [[ -z "$PATH_input" ]] || [[ -z "$PATH_output_tmp" ]] || [[ -z "$filetype" ]]
then
	echo "Please provide all required arguments (-i PATH_input, -o PATH_output and -n filetype)! use starter_Mashtree.sh -h for help on syntax"
	exit
fi

###############################################################################CODE#################################################################################

###Create folders###
filetype_c=$(echo $filetype | cut -d '.' -f 2)
PATH_output=$PATH_output_tmp/Mashtree/$filetype_c
mkdir -p $PATH_output
#mkdir -p $PATH_tmp

###Activate conda environment###

eval "$(conda shell.bash hook)"
conda activate $conda_env

###Create info file###

date > $PATH_output/info.txt
echo "Version script: "$version >> $PATH_output/info.txt
echo "mashtree version: "$(mashtree --version) >> $PATH_output/info.txt
echo "Sketch size: "$SketchSize >> $PATH_output/info.txt
echo "kmer length: "$kmerLength >> $PATH_output/info.txt
echo "Bootstrap performed: "$Bootstrap >> $PATH_output/info.txt
echo "replications (for bootstrap): "$reps >> $PATH_output/info.txt
echo "Input files: "$PATH_input >> $PATH_output/info.txt
echo "Output files:"$PATH_output >> $PATH_output/info.txt
#echo "Temporary directory: "$PATH_tmp >> $PATH_output/info.txt
echo "Amount of threads used: "$cpu >> $PATH_output/info.txt
echo "Sample set: "$set >> $PATH_output/info.txt
echo "Conda environment: "$conda_env >> $PATH_output/info.txt

if [[ "$filetype" == ".fastq" ]] || [[ "$filetype" == ".fastq.gz" ]]
then
	echo "Amount of samples in input folder:" $(ls $PATH_input/*1$filetype | wc -l) >> $PATH_output/info.txt
else
	echo "Amount of samples in input folder:" $(ls $PATH_input/*$filetype | wc -l) >> $PATH_output/info.txt
fi

if [[ "$Bootstrap" == "No" ]]
then
	mashtree --mindepth $MinDepth --kmerlength $kmerLength --sketch-size $SketchSize --outmatrix $PATH_output/$set"_"distance --numcpus $cpu $PATH_input/*$filetype > $PATH_output/$set"_mash.dnd"
elif [[ "$Bootstrap" == "Yes" ]]
then
	#find $PATH_input -name "*1$filetype" > $PATH_input/fastas.fofn
	#mashtree_bootstrap.pl --reps $reps --outmatrix $PATH_output/$set"_"distance --numcpus $cpu --mindepth 0 --kmerlength $kmerLength --sketch-size $SketchSize --tempdir $PATH_tmp -- --file-of-files $PATH_input/fastas.fofn > $PATH_output/$set"_mash_bootstrap.dnd"
	echo "This functionality does not work yet"
fi
date >> $PATH_output/info.txt
#mv $PATH_tmp/* $PATH_output/$set

################################################################################################################################################################
# Usage: mashtree [options] *.fastq *.fasta *.gbk *.msh > tree.dnd
# NOTE: fastq files are read as raw reads;
      # fasta, gbk, and embl files are read as assemblies;
      # Input files can be gzipped.
# --tempdir            ''   If specified, this directory will not be
                          # removed at the end of the script and can
                          # be used to cache results for future
                          # analyses.
                          # If not specified, a dir will be made for you
                          # and then deleted at the end of this script.
# --numcpus            1    This script uses Perl threads.
# --outmatrix          ''   If specified, will write a distance matrix
                          # in tab-delimited format
# --file-of-files           If specified, mashtree will try to read
                          # filenames from each input file. The file of
                          # files format is one filename per line. This
                          # file of files cannot be compressed.
# --outtree                 If specified, the tree will be written to
                          # this file and not to stdout. Log messages
                          # will still go to stderr.
# --version                 Display the version and exit

# TREE OPTIONS
# --truncLength        250  How many characters to keep in a filename
# --sort-order         ABC  For neighbor-joining, the sort order can
                          # make a difference. Options include:
                          # ABC (alphabetical), random, input-order

# MASH SKETCH OPTIONS
# --genomesize         5000000
# --mindepth           5    If mindepth is zero, then it will be
                          # chosen in a smart but slower method,
                          # to discard lower-abundance kmers.
# --kmerlength         21
# --sketch-size        10000