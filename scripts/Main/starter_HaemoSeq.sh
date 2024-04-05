#! /bin/bash

###Usage###
# bash $PathToScript/starter_HaemoSeq.sh

###Author###
echo "hello World, this script was written by Margo Diricks (mdiricks@fz-borstel.de)!"

###Function###
echo "This is a starter script for HaemoSeq, a pipeline for analyis of WGS data from Haemophilus bacteria"

#######################
###REQUIRED PACKAGES###
#######################

#Create conda environments with the names in between "". Alternatively, create one environment. However, this can lead to dependency issues. 

conda_env_SRST2="SRST2" #Required packages: SRST2 (tested with v0.2.0) and GNU parallel (tested with v. 20210222)
conda_env_spades="spades" #Required packages: spades (tested with v3.15.5)
conda_env_shovill="Shovill" #Required packages: shovill (tested with v1.1.0)
conda_env_platon="platon" #Required packages: platon (tested with v1.6)
conda_env_mashtree="mashtree" #Required packages: mashtree (tested with v.1.2.0)
conda_env_kraken2="kraken2" #Required packages: Kraken (tested with v.2.1.2) and ktImportText. Note that you also need to install an appropriate kraken2 database and krona database!
conda_env_seqkit="seqkit" #Required packages: seqkit (tested with v2.3.1) and pigz tested with (v2.6)
conda_env_multiqc="multiqc" #Required packages: fastqc (tested with v0.11.4) and multiqc (tested with v1.13.dev0)
conda_env_fastp="fastp" #Required packages: fastp (tested with v.0.23.2)
conda_env_AMRfinder="AMRFinder" #Required packages: AMRfinder (tested with v.3.11.2 and database version 2022-12-19.1)

##############################################
###!!!PARAMATERS THAT NEED TO BE CHANGED!!!###
##############################################

##Paths to Input/Output##
#PATH_scripts: Path to folder where scripts are stored
PATH_scripts="$PATH/to/scripts"
#PATH_fastQ:Path to folder where fastQ files are stored
PATH_fastQ="$PATH/to/FastQ"
#PATH_output: Path to folder where result files will be stored
PATH_output="$PATH/to/HaemoSeq"

#Paths to databases
#db_classification_serotyping: Path to folder where classification/serotyping database is stored (Diricks et al., 2022). Download from https://github.com/ngs-fzb/HaemoSeq/tree/main/Databases.
db_classification_serotyping="$Path_To/Class_Haemophilus_and_Serotyping-v2.0.fasta"
#db_platon: Path to folder where platon database files are stored. Download from https://github.com/oschwengers/platon#database or https://zenodo.org/search?page=1&size=20&q=conceptrecid:3349651&all_versions&sort=-version
db_platon="$PATH/to/platon_db_folder"
#db_PLSDB_SRST2: Path to a PLSDB plasmid database. Download from https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/
db_PLSDB_SRST2="$PATH/to/plsdb.fna"
#db_custom_SRST2: Path to a custom database of genes/sequences that you want to look for. 
db_custom_SRST2="$PATH/to/custom_db.fna"
#db_MLST_SRST2: Path to folder where pubMLST databases (genes and allele nomenclature) will be stored. These databases will be automatically downloaded by HaemoSeq.
db_MLST_SRST2="$PATH/to/pubMLST_database_folder"
#db_kraken2: Path to folder where kraken2 database files are stored. Download from https://benlangmead.github.io/aws-indexes/k2
db_kraken2="$PATH/to/Kraken_db_folder"
#db_krona: Path to folder where krona files are stored. If left empty, the default folder will be used (i.e. conda environment). Download/update using ktUpdateTaxonomy.sh --> if error see https://github.com/bioconda/bioconda-recipes/issues/10959
db_krona="$PATH/to/krona_db_folder"
#db_AMRfinder: Path to folder where AMRfinder database files are stored. Download with amrfinder_update -d $Path/to/AMRfinder_db_folder
db_AMRfinder="$PATH/to/AMRfinder_db_folder" 


species="Haemophilus influenzae" #Required for MLST, which is currently only available for H. influenzae
genome_size="1.8M" #Expected genome size

###########################################
###!!!CHOOSE ANALYSES YOU WANT TO DO!!! ###
###########################################

#[options: Yes or No]
#Estimated time for 30 HI samples (simulated from public assemblies using dwgsim; 70-90 Mb read files) using Intel® Xeon® processor (E5-2650 v4 @ 2.2 GHz) with 48 Gb RAM and 16 CPUs
##Read quality control##
Do_multiqc="No"
##Assembly##
Do_assembly="No" #Using shovill and spades
#Estimate time: 4 hours 
##Taxonomy and phylogeny##
Do_classification_serotyping="Yes" #Based on Haemophilus species specific markers and capsule genes (Diricks et al. 2022)
#Estimated time: 3 min
Do_kraken2="No"
#Estimated time: 324 min ( 5.7 hours)
Do_mash_fastQ="No" #Can be used as a quick quality check of the fastQ files!
#Estimated time: 7 min
Do_mash_fastA="No"
#Estimated time: 10 sec
Do_MLST_fastQ="No" # If results are unclear: compare with mlst_tseemann or pubMLST web-page?
#Estimated time: 25 min
##Plasmid prediction##
Do_SRST2_PLSDB="No" #Against PLSDB database
Do_SRST2_customDB="No" # Against custom database
#Estimated time: Depends on size of the database
Do_plasmidspades="No"
#Estimated time: ~ min
Do_platon="No" # Requires assemblies!
#Estimated time: ~ min
Do_amrfinder="No" #Requires assemblies!
#Estimated time: 

#############################################
###OPTIONAL PARAMATERS THAT YOU CAN CHANGE###
#############################################

set="Test" #Name of sample set - used for file naming"
#cpu: Resources to be used (max. amount depends on your hardware configuration)
cpu=$(lscpu | grep -E '^CPU\(' | awk '{print $(NF)}') # $(lscpu | grep -E '^CPU\(' | awk '{print $(NF)}') Uses all available threads

Fw="_R1" #Change if your fastQ files are not SampleName_R1.fastq.gz
Rv="_R2" #Change if your fastQ files are not SampleName_R2.fastq.gz
PATH_tmp="$HOME/tmp" #Directory were temporary files are stored; Default: $HOME/tmp
#Note: It is important that PATH_tmp has enough space, otherwise some tools might fail!! 
PATH_fastA="$PATH_output/Assemblies/FinalAssemblies" # Path to folder where fastA files are stored; $PATH_output/Assemblies/FinalAssemblies = Default output if Do_assembly = Yes
cov=100 #reads will be downsampled to reach this theoretic coverage; Default: 100"
ass="spades" #assembler that needs to be used (choose skesa velvet megahit or spades); Default: spades"

Do_extract_plasmids="No" #Extract fasta sequences for which a hit was found (in SRST2_PLSDB analysis); Default: No
#Estimated time: ~0.5 min per plasmid
conda_env_seqkit="seqkit" #Required packages: seqkit


##########################################
###NO CHANGES REQUIRED BELOW THIS POINT###
##########################################

###############################################################################CODE#################################################################################

###QUALITY CONTROL READS###
if [[ "$Do_multiqc" == "Yes" ]]
then
	if [[ ! -z "$PATH_fastQ" ]] || [[ ! -z "$PATH_output" ]]
	then
		echo "Starting to perform quality control on reads"
		bash $PATH_scripts/starter_multiqc.sh -i $PATH_fastQ -o $PATH_output -c $cpu -e $conda_env_multiqc
	else
		echo "Please provide all required arguments (-i PATH_fastQ and -o PATH_output )! use starter_multiqc.sh -h for help on syntax"
	fi
elif [[ "$Do_multiqc" == "No" ]]
then
	echo "quality control on reads skipped."
else 
	echo "Please decide on whether you want to do quality control on reads or not!"
fi

###ASSEMBLY###
if [[ "$Do_assembly" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]] && [[ ! -z "$genome_size" ]]
	then
		echo "Starting to assemble"
		bash $PATH_scripts/starter_Assemble_usingShovill.sh -i $PATH_fastQ -o $PATH_output -g $genome_size -d $cov -c $cpu -e $conda_env_shovill -f $Fw -r $Rv -s $set -a $ass
	else
		echo "Please provide all required arguments: -i PATH_fastQ, -o PATH_output and -d db_classification_serotyping!"
	fi
elif [[ "$Do_assembly" == "No" ]]
then
	echo "Assembly skipped."
else 
	echo "Please decide on whether you want to make an assembly or not!"
fi

###TAXONOMY AND PHYLOGENY###
if [[ "$Do_classification_serotyping" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]] && [[ ! -z "$db_classification_serotyping" ]]
	then 
		echo "Starting species classification and serotyping"
		bash $PATH_scripts/starter_classification-serotyping_Haemophilus_SRST2.sh -i $PATH_fastQ -o $PATH_output -d $db_classification_serotyping -c $cpu -e $conda_env_SRST2 -f $Fw -r $Rv -s $set -t $PATH_tmp
	else
		echo "Please provide all required arguments: -i PATH_fastQ, -o PATH_output and -d db_classification_serotyping!"
	fi
elif [[ "$Do_classification_serotyping" == "No" ]]
then
	echo "species classification and serotyping skipped."
else 
	echo "Please decide on whether you want to do species classification and serotyping or not (Choose Yes or No)!"
fi

if [[ "$Do_mash_fastQ" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]]
	then 
		echo "Starting mash on FastQ files"
		bash $PATH_scripts/starter_Mashtree.sh -i $PATH_fastQ -o $PATH_output -n .fastq.gz -e $conda_env_mashtree
	else
		echo "Please provide all required arguments: -i PATH_fastQ and -o PATH_output!"
	fi
elif [[ "$Do_mash_fastQ" == "No" ]]
then
	echo "mash on FastQ files skipped."
else 
	echo "Please decide on whether you want to do mash on FastQ files or not (Choose Yes or No)!"
fi

if [[ "$Do_mash_fastA" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastA" ]]
	then 
		echo "Starting mash on FastA files"
		bash $PATH_scripts/starter_Mashtree.sh -i $PATH_fastA -o $PATH_output -n .fasta -e $conda_env_mashtree
	else
		echo "Please provide all required arguments: -i PATH_fastA and -o PATH_output!"
	fi
elif [[ "$Do_mash_fastA" == "No" ]]
then
	echo "mash on FastA files skipped."
else 
	echo "Please decide on whether you want to do mash on FastA files or not (Choose Yes or No)!"
fi

if [[ "$Do_kraken2" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]] && [[ ! -z "$db_kraken2" ]] && [[ ! -z "$db_krona" ]]
	then 
		echo "Starting kraken2 on FastQ files"
		bash $PATH_scripts/starter_kraken2.sh -i $PATH_fastQ -o $PATH_output -d $db_kraken2 -k $db_krona -e $conda_env_kraken2
		perl $PATH_scripts/kraken_parse_results.v2.0.0.pl -s "$species" $PATH_output/kraken2/*.report
	else
		echo "Please provide all required arguments: -i PATH_fastQ, -o PATH_output, db_krona and db_kraken2!"
	fi
elif [[ "$Do_kraken2" == "No" ]]
then
	echo "kraken2 on FastQ files skipped."
else 
	echo "Please decide on whether you want to do kraken2 on FastQ files or not (Choose Yes or No)!"
fi

if [[ "$Do_MLST_fastQ" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]] && [[ ! -z "$db_MLST_SRST2" ]] && [[ ! -z "$species" ]]
	then 
		echo "Starting MLST on FastQ files"
		bash $PATH_scripts/starter_MLST_SRST2.sh -i $PATH_fastQ -o $PATH_output -d $db_MLST_SRST2 -e $conda_env_SRST2 -p "$species"
	else
		echo "Please provide all required arguments: -i PATH_fastQ, -o PATH_output, db_MLST_SRST2 and species!"
	fi
elif [[ "$Do_MLST_fastQ" == "No" ]]
then
	echo "MLST on FastQ files skipped."
else 
	echo "Please decide on whether you want to do MLST on FastQ files or not (Choose Yes or No)!"
fi

###TAXONOMY AND RESISTANCE PREDICTION###

if [[ "$Do_NTMprofiler" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]]
	then 
		echo "Starting running NTMprofiler"
		bash $PATH_scripts/starter_NTMprofiler.sh -i $PATH_fastQ -o $PATH_output -n .fastq.gz -c $cpu -e $conda_env_NTMprofiler -f $Fw -r $Rv -s $set
	else
		echo "Please provide all required arguments: -i PATH_fastQ and -o PATH_output!"
	fi
elif [[ "$Do_NTMprofiler" == "No" ]]
then
	echo "running NTMprofiler skipped."
else 
	echo "Please decide on whether you want to run NTMprofiler or not (Choose Yes or No)!"
fi

####PLASMID DETECTION###
###PART1
if [[ "$Do_plasmidspades" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]]
	then
		echo "Starting plasmid assembly using plasmidspades"
		bash $PATH_scripts/starter_plasmidspades.sh -i $PATH_fastQ -o $PATH_output -c $cpu -e $conda_env_spades -f $Fw -r $Rv -s $set -t $PATH_tmp
		bash $PATH_scripts/starter_Platon.sh -i $PATH_output/Plasmidspades/FinalAssemblies -o $PATH_output/Plasmidspades/FinalAssemblies -d $db_platon -c $cpu -e $conda_env_platon -s $set
	else
		echo "Please provide all required arguments: -i PATH_fastQ -o PATH_output!"
	fi
elif [[ "$Do_plasmidspades" == "No" ]]
then
	echo "Plasmid assembly using plasmidspades skipped."
else 
	echo "Please decide on whether you want to do assembly using plasmidspades or not (Choose Yes or No)!"
fi

###PART2b
if [[ "$Do_SRST2_customDB" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]] && [[ ! -z "$db_custom_SRST2" ]]
	then
		customDB=$(basename $db_custom_SRST2| cut -d '.' -f 1)
		customDB_SRST2=$customDB"_SRST2.fasta"
		if [ ! -f "$(dirname $db_custom_SRST2)/$customDB_SRST2" ]
		then
			echo "Starting converting custom database into SRST2 compatible database"
			bash $PATH_scripts/customDB_to_SRST2_db.sh -d $db_custom_SRST2
		else
			echo "SRST2 compatible custom database found"
		fi
		echo "Starting plasmid detection with SRST2 using custom database"
		bash $PATH_scripts/starter_customDB_SRST2.sh -i $PATH_fastQ -o $PATH_output -d $(dirname $db_custom_SRST2)/$customDB_SRST2 -c $cpu -e $conda_env_SRST2 -f $Fw -r $Rv -s $set -t $PATH_tmp
	else
		echo "Please provide all required arguments: -i PATH_fastQ -o PATH_output and -d db_custom_SRST2!"
	fi
elif [[ "$Do_SRST2_customDB" == "No" ]]
then
	echo "Plasmid detection with SRST2 using custom database skipped."
else 
	echo "Please decide on whether you want to do plasmid detection with SRST2 using custom database or not (Choose Yes or No)!"
fi

###PART2
if [[ "$Do_SRST2_PLSDB" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastQ" ]] && [[ ! -z "$db_PLSDB_SRST2" ]]
	then
		if [ ! -f "$(dirname $db_PLSDB_SRST2)/plsdb_SRST2.fasta" ]
		then
			echo "Starting converting plsdb database into SRST2 compatible database"
			bash $PATH_scripts/PLSDB_to_SRST2_db.sh -d $db_PLSDB_SRST2
		else
			echo "SRST2 compatible plsdb database found"
		fi
		echo "Starting plasmid detection with SRST2 using PLSDB database"
		bash $PATH_scripts/starter_PLSDB_SRST2.sh -i $PATH_fastQ -o $PATH_output -d $(dirname $db_PLSDB_SRST2)/plsdb_SRST2.fasta -c $cpu -e $conda_env_SRST2 -f $Fw -r $Rv -s $set -t $PATH_tmp
	else
		echo "Please provide all required arguments: -i PATH_fastQ -o PATH_output and -d db_PLSDB_SRST2!"
	fi
elif [[ "$Do_SRST2_PLSDB" == "No" ]]
then
	echo "Plasmid detection with SRST2 using PLSDB database skipped."
else 
	echo "Please decide on whether you want to do plasmid detection with SRST2 using PLSDB database or not (Choose Yes or No)!"
fi

if [[ "$Do_extract_plasmids" == "Yes" ]]
then
	bash $PATH_scripts/Extract_fasta_from_multifasta.sh -i $PATH_output/PLSDB_SRST2/Plasmid_hits.txt -o $PATH_output/PLSDB_SRST2/Plasmid_hits -e $conda_env_seqkit -m $(dirname $db_PLSDB_SRST2)/plsdb_SRST2.fasta
	bash $PATH_scripts/starter_Mashtree.sh -i $PATH_output/PLSDB_SRST2/Plasmid_hits -o $PATH_output/PLSDB_SRST2/Plasmid_hits -e $conda_env_mashtree -n .fasta
fi

###PART3
if [[ "$Do_platon" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastA" ]] && [[ ! -z "$db_platon" ]]
	then
		echo "Plasmid detection with platon"
		bash $PATH_scripts/starter_platon.sh -i $PATH_fastA -o $PATH_output -d $db_platon -c $cpu -e $conda_env_platon -s $set
	else
		echo "Please provide all required arguments: -i PATH_fastA, -o PATH_output and -d db_platon!"
	fi
elif [[ "$Do_platon" == "No" ]]
then
	echo "Plasmid detection with platon skipped."
else 
	echo "Please decide on whether you want to do plasmid detection with platon or not (Choose Yes or No)!"
fi

echo "Script Finished!"


####RESISTANCE AND VIRULENCE DETECTION###
if [[ "$Do_amrfinder" == "Yes" ]]
then
	if [[ ! -z "$PATH_output" ]] && [[ ! -z "$PATH_fastA" ]] && [[ ! -z "$db_AMRfinder" ]]
	then
		echo "Resistance gene and virulence detection with AMRFinderplus"
		bash $PATH_scripts/starter_amrfinder.sh -i $PATH_fastA -o $PATH_output -c $cpu -e $conda_env_AMRfinder -d $db_AMRfinder
	else
		echo "Please provide all required arguments: -i PATH_fastA, -d db_AMRfinder and -o PATH_output!"
	fi
elif [[ "$Do_amrfinder" == "No" ]]
then
	echo "Resistance gene and virulence detection with AMRFinderplus skipped."
else 
	echo "Please decide on whether you want to do resistance gene and virulence detection with AMRFinderplus or not (Choose Yes or No)!"
fi

echo "Script Finished!"
exit


####################################################################CODE THAT MIGHT BE USED IN ADDITION#################################################################################
###Download and extract PLSDB database### 
#https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/
#db_PLSDB_SRST2: Download FASTA archive. Use bzip2 -d db_name to extract
#Convert fna file to SRST2 compatible database: automated in NTMseq pipeline

###Download and extract latest platon database### 
#Info: Platon depends on a custom database based on MPS, RDS, RefSeq Plasmid database, PlasmidFinder db as well as manually curated MOB HMM models from MOBscan, custom conjugation and replication HMM models and oriT sequences from MOB-suite.
#Download page: https://github.com/oschwengers/platon#database or https://zenodo.org/search?page=1&size=20&q=conceptrecid:3349651&all_versions&sort=-version

