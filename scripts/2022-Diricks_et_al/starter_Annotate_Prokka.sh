#!/usr/bin/env bash

###Author###
echo "hello World, this script was written by Margo Diricks!"

###Function###
#This script annotates fasta files with Prokka

###Required packages###
#prokka: https://github.com/tseemann/prokka
	#
	
###Usage###
# bash $HOME/Scripts/Name_Script.sh"

###Inputfiles###
#Type: fasta files
#Path to input files:
PATH_input="$HOME/YourFolder"

###Parameters###
#Minimum contig length
col="200"
#Centre
cr="YourCenter"
#Genus
ge="Haemophilus"
#Species
sp=""
#Strain
str=""

###Output files###
#1. Directory where the final files will be stored
PATH_output="$HOME/YourFolder"

#Amount of cpus to be used(0=all)
cpu=0

################################################################################################################################################################
mkdir $PATH_output
for fa in $PATH_input/*.fasta
do
	SampleName=$(basename $fa| cut -d '_' -f 1)
	if [ -d $PATH_output/$SampleName ]
	then
		echo "Sample was already analysed"
	else
		prokka --outdir $PATH_output/$SampleName --prefix $SampleName --addgenes --mincontiglen $col --centre $cr --genus $ge --usegenus --species $sp --strain $str --cpus $cpu $fa
	fi
done
mkdir $PATH_output/FinalAnnotatedFiles
mv $PATH_output/*/*.gff $PATH_output/FinalAnnotatedFiles

################################################################################################################################################################
# General:
  # --help            This help
  # --version         Print version and exit
  # --citation        Print citation for referencing Prokka
  # --quiet           No screen output (default OFF)
  # --debug           Debug mode: keep all temporary files (default OFF)
# Setup:
  # --listdb          List all configured databases
  # --setupdb         Index all installed databases
  # --cleandb         Remove all database indices
  # --depends         List all software dependencies
# Outputs:
  # --outdir [X]      Output folder [auto] (default '')
  # --force           Force overwriting existing output folder (default OFF)
  # --prefix [X]      Filename output prefix [auto] (default '')
  # --addgenes        Add 'gene' features for each 'CDS' feature (default OFF)
  # --locustag [X]    Locus tag prefix (default 'PROKKA')
  # --increment [N]   Locus tag counter increment (default '1')
  # --gffver [N]      GFF version (default '3')
  # --compliant       Force Genbank/ENA/DDJB compliance: --genes --mincontiglen 200 --centre XXX (default OFF)
  # --centre [X]      Sequencing centre ID. (default '')
# Organism details:
  # --genus [X]       Genus name (default 'Genus')
  # --species [X]     Species name (default 'species')
  # --strain [X]      Strain name (default 'strain')
  # --plasmid [X]     Plasmid name or identifier (default '')
# Annotations:
  # --kingdom [X]     Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
  # --gcode [N]       Genetic code / Translation table (set if --kingdom is set) (default '0')
  # --prodigaltf [X]  Prodigal training file (default '')
  # --gram [X]        Gram: -/neg +/pos (default '')
  # --usegenus        Use genus-specific BLAST databases (needs --genus) (default OFF)
  # --proteins [X]    Fasta file of trusted proteins to first annotate from (default '')
  # --hmms [X]        Trusted HMM to first annotate from (default '')
  # --metagenome      Improve gene predictions for highly fragmented genomes (default OFF)
  # --rawproduct      Do not clean up /product annotation (default OFF)
# Computation:
  # --fast            Fast mode - skip CDS /product searching (default OFF)
  # --cpus [N]        Number of CPUs to use [0=all] (default '8')
  # --mincontiglen [N] Minimum contig size [NCBI needs 200] (default '1')
  # --evalue [n.n]    Similarity e-value cut-off (default '1e-06')
  # --rfam            Enable searching for ncRNAs with Infernal+Rfam (SLOW!) (default '0')
  # --norrna          Don't run rRNA search (default OFF)
  # --notrna          Don't run tRNA search (default OFF)
  # --rnammer         Prefer RNAmmer over Barrnap for rRNA prediction (default OFF)