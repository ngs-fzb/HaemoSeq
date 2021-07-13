#!/usr/bin/env bash

###Author###
echo "hello World, this script was written by Margo!"

###Function###
echo "Running this script will reveal the identity of the Haemophilus spp. that is growing in your culture! In addition, the serotype will be determined if applicable"
echo "This script detects HIvsHH marker genes, haemin synthesis genes and capsule loci genes!"
echo "Don´t forget to install SRST2 and GNU parallel, e.g. with conda, before you run this script!"
echo "Don´t forget to activate the right conda environment before you run this script!"
echo "Make sure your FastQ file are named according to following convention: SampleName_*_R1.fastq.gz or change the script!"
echo "Don´t forget to change the parameters PATH_input, PATH_output, jobs and others if necessary!"
Analysis="Classification_Serotyping"

###Required packages###
# SRST2 and GNU parallel --> these can be installed e.g. using conda.

###Usage###
# bash $HOME/Scripts/Name_Script.sh"

###Input files###
#Directory where you have stored the (raw) paired-end illumina fastQ files
PATH_input=
#PATH_input="/mnt/c/Users/diric/SciData/Genomes_Haemophilus/FastQ" #Example Path

#Read notation and extension
FastqType=".fastq.gz" # Extension of your fastq file
Fw="_R1" # Forward read notation, e.g. R1
Rv="_R2" # Reverse read notation, e.g. R2

#Path to nucleotide reference database
geneDB=
#geneDB="/mnt/c/Users/diric/Desktop/Bashscripts2/Class_Haemophilus_and_Serotyping.fasta" #Example Path

###Parameters###
# Additional parameters can be added if required (see below)
#1. Maximum number of mismatches allowed during mapping - M
GeneMaxMis=50
#2. Maximum %divergence cutoff for gene reporting (in genes report) - D
MaxDiv=15
#3. Minimum %coverage cutoff for gene reporting (in full genes and genes report) - C
MinCov=90
#Note: default values are 50M-15D-90C and decision rules are based on these thresholds. Severely truncated genes or homologues with not enough similarity will not be detected!

#Amount of samples to be analyzed in parallel 
Jobs=1 #Jobs=0 means you will run as many jobs in parallel as possible. Please note that this can sometimes lead to problems. Preferably use the amount of CPUs/threads you have available.
#Important note: if it is the first time you use a reference database (fasta file defined as geneDB), you need to set this value once at 1 to generate the fai and bt2 files of the reference database. Afterwards, you can set the amount of jobs depending on your computational resources. Please consult the GNU parallel manual for more info. 

###Output files###
#1. Directory where the final files will be stored
PATH_output=
#PATH_output="/mnt/c/Users/diric/SciData/Genomes_Haemophilus/"$Analysis"/Revised_Script/Mismatch"$GeneMaxMis"_maxdiv"$MaxDiv"_mincov"$MinCov"_Jobs"$Jobs #Example Path
#2. Directory where the temporary files will be stored
PATH_temp="$HOME/tmp"
#3. Name of the decision and itol output files that will be generated automatically
Name_itolFile_Decision_Species=$Analysis"_SRST2_itol_Decision_Species" #This file can be dragged and dropped over an itol tree to visualize the species detected
Name_itolFile_Binary=$Analysis"_SRST2_itol_Binary" #This file can be dragged and dropped over an itol tree to visualize the presence/absence of the genes
Name_DecisionFile=$Analysis"_finalDecision" #Contains the identity of your strain in a simple text file
Name_itolFile_Decision_Serotype=$Analysis"_SRST2_itol_Decision_Serotype"
Name_itolFile_Decision_Subspecies=$Analysis"_SRST2_itol_Decision_Subspecies"


###EstimatedTime###
#<2 min per sample (100x coverage)
#Analysis of 200 samples with 100x coverage in less then 10 min on a machine with 16 cores (32 threads)

################################################################################################################################################################
#Species identification in parallel
source `which env_parallel.bash`
env_parallel -j $Jobs --compress --tmpdir $PATH_temp 'srst2 --input_pe {} $(echo {} | sed "s/$Fw$FastqType$/$Rv$FastqType/1") --forward $Fw --reverse $Rv --log --gene_db $geneDB --gene_max_mismatch $GeneMaxMis --max_divergence $MaxDiv --min_coverage $MinCov --output $PATH_output/$(echo {/}| cut -d '_' -f 1) --use_existing_scores --use_existing_pileup' ::: $PATH_input/*$Fw$FastqType

#Remove pileup and bam files, only keep result files
rm $PATH_output/*.pileup
rm $PATH_output/*.bam
rm $PATH_output/*.sam


#Create a combined file
srst2 --prev_output $PATH_output/*__genes__*__results.txt --output $PATH_output/$Analysis

#Make a Decision file with (sub)species and serotype determination
echo "###Species determination based on presence of marker genes###" > $PATH_output/$Name_DecisionFile.txt
echo "Only samples with a valid hit for at least one of the genes are listed here." >> $PATH_output/$Name_DecisionFile.txt
echo "Valid hit = hit with coverage > $MinCov and divergence < $MaxDiv towards one of the reference alleles." >> $PATH_output/$Name_DecisionFile.txt
echo "Minor contaminations (with incomplete marker patterns) are not reported here." >> $PATH_output/$Name_DecisionFile.txt
echo "#Legend#" >> $PATH_output/$Name_DecisionFile.txt
echo "HD: Haemophilus ducreyi: strains with a valid hit for nadV" >> $PATH_output/$Name_DecisionFile.txt
echo "HpI: Haemophilus parainfluenzae: strains with valid hit for seven haemin genes hemA/C/D/E/G/H/L/N-HpI" >> $PATH_output/$Name_DecisionFile.txt
echo "HpH: Haemophilus parahaemolyticus: strains with valid hit for at least seven haemin genes hemA/C/D/E/G/H/L/N-HpH" >> $PATH_output/$Name_DecisionFile.txt
echo "HHi: Haemophilus haemolyticus subsp. intermedius: strains with valid hit for nine haemin genes hemA/B/C/D/E/G/H/L/N-HH" >> $PATH_output/$Name_DecisionFile.txt
echo "HP: Haemophilus pittmaniae: strains with valid hit for at least seven haemin genes hemA/C/D/E/G/H/L/N-HP" >> $PATH_output/$Name_DecisionFile.txt
echo "HS: Haemophilus sputorum: strains with valid hit for at least seven haemin genes hemA/C/D/E/G/H/L/N-HS" >> $PATH_output/$Name_DecisionFile.txt
echo "HH: Haemophilus haemolyticus: strains with valid hit for at least four HH marker genes" >> $PATH_output/$Name_DecisionFile.txt
echo "HI: Haemophilus influenzae: strains with valid hit for at least four HI marker genes" >> $PATH_output/$Name_DecisionFile.txt
echo "Serotype a-f: valid hits for all capsule loci were found" >> $PATH_output/$Name_DecisionFile.txt
echo "Multiple serotypes: valid hits were found for multiple region II clusters" >> $PATH_output/$Name_DecisionFile.txt
echo "Capsule-deficient: not for all capsule loci genes valid hits were found" >> $PATH_output/$Name_DecisionFile.txt
echo "NTHi: none of the capsule genes were found" >> $PATH_output/$Name_DecisionFile.txt
echo "Mix: gene patterns from multiple Haemophilus spp. are detected. Please investigate further." >> $PATH_output/$Name_DecisionFile.txt
echo "Unknown: a pattern of genes was detected that has not yet been linked with a specific Haemophilus species. Please investigate further." >> $PATH_output/$Name_DecisionFile.txt
echo "Details about hits can be found in *__fullgenes__* files" >> $PATH_output/$Name_DecisionFile.txt
echo "#Results#" >> $PATH_output/$Name_DecisionFile.txt

#Make a Decision Itol file for species determination
echo "DATASET_COLORSTRIP" > $PATH_output/$Name_itolFile_Decision_Species.txt 
echo "SEPARATOR TAB" >> $PATH_output/$Name_itolFile_Decision_Species.txt
echo "DATASET_SCALE	0" >> $PATH_output/$Name_itolFile_Decision_Species.txt
#Label is used in the legend table (can be changed later)
echo "DATASET_LABEL	"$Analysis"_MaxMismatch"$GeneMaxMis"_MaxDiv"$MaxDiv"_MinCov"$MinCov"_Decision_Species" >> $PATH_output/$Name_itolFile_Decision_Species.txt
#Dataset color (can be changed later)
echo "COLOR	#ff0000" >> $PATH_output/$Name_itolFile_Decision_Species.txt
#Define colors for each individual field column
echo "FIELD_COLORS	#808080	#9bddff	#4cbb17	#eee600	#ae0c00	#ff00ff	#ff9933	#dbd7d2	#fdd5b1" >> $PATH_output/$Name_itolFile_Decision_Species.txt
echo "FIELD_LABELS	HD	HH	HI	HpI	HpH	HP	HS	Unknown	Mixed" >> $PATH_output/$Name_itolFile_Decision_Species.txt
echo "FIELD_SHAPES	1	1	1	1	1	1	1	1	1" >> $PATH_output/$Name_itolFile_Decision_Species.txt
#Legend
echo "LEGEND_TITLE	Species" >> $PATH_output/$Name_itolFile_Decision_Species.txt
echo "LEGEND_SHAPES	1	1	1	1	1	1	1	1	1" >> $PATH_output/$Name_itolFile_Decision_Species.txt
echo "LEGEND_COLORS	#808080	#9bddff	#4cbb17	#eee600	#ae0c00	#ff00ff	#ff9933	#dbd7d2	#fdd5b1" >> $PATH_output/$Name_itolFile_Decision_Species.txt
echo "LEGEND_LABELS	HD	HH	HI	HpI	HpH	HP	HS	Unknown	Mixed" >> $PATH_output/$Name_itolFile_Decision_Species.txt
#Data
echo "DATA" >> $PATH_output/$Name_itolFile_Decision_Species.txt

#Make a Decision Itol file for subspecies determination and serotyping results
echo "DATASET_COLORSTRIP" > $PATH_output/$Name_itolFile_Decision_Serotype.txt 
echo "SEPARATOR TAB" >> $PATH_output/$Name_itolFile_Decision_Serotype.txt
echo "DATASET_SCALE	0" >> $PATH_output/$Name_itolFile_Decision_Serotype.txt
#Label is used in the legend table (can be changed later)
echo "DATASET_LABEL	"$Analysis"_MaxMismatch"$GeneMaxMis"_MaxDiv"$MaxDiv"_MinCov"$MinCov"_Decision_Serotype" >> $PATH_output/$Name_itolFile_Decision_Serotype.txt
#Dataset color (can be changed later)
echo "COLOR	#ff0000" >> $PATH_output/$Name_itolFile_Decision_Serotype.txt
#Legend
echo "LEGEND_TITLE	Serotype" >> $PATH_output/$Name_itolFile_Decision_Serotype.txt
echo "LEGEND_SHAPES	1	1	1	1	1	1	1	1	1" >> $PATH_output/$Name_itolFile_Decision_Serotype.txt
echo "LEGEND_COLORS	#000000	#fd0e35	#800080	#eee600	#d2691e	#808080	#dbd7d2	#fdd5b1	#006600" >> $PATH_output/$Name_itolFile_Decision_Serotype.txt
echo "LEGEND_LABELS	HIa	HIb	HIc	HId	HIe	HIf	Capsule-deficient	Mix of serotypeable strains	NTHi" >> $PATH_output/$Name_itolFile_Decision_Serotype.txt
#Data
echo "DATA" >> $PATH_output/$Name_itolFile_Decision_Serotype.txt

#Make a Decision Itol file for subspecies determination and serotyping results
echo "DATASET_COLORSTRIP" > $PATH_output/$Name_itolFile_Decision_Subspecies.txt 
echo "SEPARATOR TAB" >> $PATH_output/$Name_itolFile_Decision_Subspecies.txt
echo "DATASET_SCALE	0" >> $PATH_output/$Name_itolFile_Decision_Subspecies.txt
#Label is used in the legend table (can be changed later)
echo "DATASET_LABEL	"$Analysis"_MaxMismatch"$GeneMaxMis"_MaxDiv"$MaxDiv"_MinCov"$MinCov"_Decision_Subspecies" >> $PATH_output/$Name_itolFile_Decision_Subspecies.txt
#Dataset color (can be changed later)
echo "COLOR	#ff0000" >> $PATH_output/$Name_itolFile_Decision_Subspecies.txt
#Legend
echo "LEGEND_TITLE	Subspecies">> $PATH_output/$Name_itolFile_Decision_Subspecies.txt
echo "LEGEND_SHAPES	1" >> $PATH_output/$Name_itolFile_Decision_Subspecies.txt
echo "LEGEND_COLORS	#93ccea" >> $PATH_output/$Name_itolFile_Decision_Subspecies.txt
echo "LEGEND_LABELS	HH subsp. intermedius" >> $PATH_output/$Name_itolFile_Decision_Subspecies.txt
#Data
echo "DATA" >> $PATH_output/$Name_itolFile_Decision_Subspecies.txt

#Make a Binary Itol file
echo "DATASET_BINARY" > $PATH_output/$Name_itolFile_Binary.txt 
echo "SEPARATOR TAB" >> $PATH_output/$Name_itolFile_Binary.txt
echo "DATASET_SCALE	0" >> $PATH_output/$Name_itolFile_Binary.txt
#label is used in the legend table (can be changed later)
echo "DATASET_LABEL	"$Analysis"_MaxMismatch"$GeneMaxMis"_MaxDiv"$MaxDiv"_MinCov"$MinCov"_Binary" >> $PATH_output/$Name_itolFile_Binary.txt
#dataset color (can be changed later)
echo "COLOR	#ff0000" >> $PATH_output/$Name_itolFile_Binary.txt
#define colors for each individual field columns
echo "FIELD_COLORS	#808080	#9bddff	#9bddff	#9bddff	#9bddff	#9bddff	#4cbb17	#4cbb17	#4cbb17	#4cbb17	#4cbb17	#556b2f	#556b2f	#556b2f	#556b2f	#000000	#fd0e35	#800080	#eee600	#d2691e	#808080	#556b2f	#556b2f	#ae0c00	#ae0c00	#ae0c00	#ae0c00	#ae0c00	#ae0c00	#ae0c00	#ae0c00	#eee600	#eee600	#eee600	#eee600	#eee600	#eee600	#eee600	#eee600	#9bddff	#9bddff	#9bddff	#9bddff	#9bddff	#9bddff	#9bddff	#9bddff	#9bddff	#ff9933	#ff9933	#ff9933	#ff9933	#ff9933	#ff9933	#ff9933	#ff9933	#ff00ff	#ff00ff	#ff00ff	#ff00ff	#ff00ff	#ff00ff	#ff00ff	#ff00ff" >> $PATH_output/$Name_itolFile_Binary.txt
echo "FIELD_LABELS	nadV	hypD	fklB	ndhI	tpd	resA	phoB	pdxT	dat	oppC	hxuB	bexA	bexB	bexC	bexD	acs1234	bcs1234	ccs1234	dcs12345	ecs12345678	fcs123	hcsA	hcsB	hemA-HpH	hemC-HpH	hemD-HpH	hemE-HpH	hemG-HpH	hemH-HpH	hemL-HpH	hemN-HpH	hemA-HpI	hemC-HpI	hemD-HpI	hemE-HpI	hemG-HpI	hemH-HpI	hemL-HpI	hemN-HpI	hemA-HH	hemB-HH	hemC-HH	hemD-HH	hemE-HH	hemG-HH	hemH-HH	hemL-HH	hemN-HH	hemA-HS	hemC-HS	hemD-HS	hemE-HS	hemG-HS	hemH-HS	hemL-HS	hemN-HS	hemA-HP	hemC-HP	hemD-HP	hemE-HP	hemG-HP	hemH-HP	hemL-HP	hemN-HP" >> $PATH_output/$Name_itolFile_Binary.txt
echo "FIELD_SHAPES	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1" >> $PATH_output/$Name_itolFile_Binary.txt
#Legend
echo "LEGEND_TITLE	Marker genes" >> $PATH_output/$Name_itolFile_Binary.txt
echo "LEGEND_SHAPES	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1" >> $PATH_output/$Name_itolFile_Binary.txt
echo "LEGEND_COLORS	#808080	#9bddff	#9bddff	#9bddff	#9bddff	#9bddff	#4cbb17	#4cbb17	#4cbb17	#4cbb17	#4cbb17	#556b2f	#556b2f	#556b2f	#556b2f	#000000	#fd0e35	#800080	#eee600	#d2691e	#808080	#556b2f	#556b2f	#ae0c00	#ae0c00	#ae0c00	#ae0c00	#ae0c00	#ae0c00	#ae0c00	#ae0c00	#eee600	#eee600	#eee600	#eee600	#eee600	#eee600	#eee600	#eee600	#9bddff	#9bddff	#9bddff	#9bddff	#9bddff	#9bddff	#9bddff	#9bddff	#9bddff	#ff9933	#ff9933	#ff9933	#ff9933	#ff9933	#ff9933	#ff9933	#ff9933	#ff00ff	#ff00ff	#ff00ff	#ff00ff	#ff00ff	#ff00ff	#ff00ff	#ff00ff" >> $PATH_output/$Name_itolFile_Binary.txt
echo "LEGEND_LABELS	nadV	hypD	fklB	ndhI	tpd	resA	phoB	pdxT	dat	oppC	hxuB	bexA	bexB	bexC	bexD	acs1234	bcs1234	ccs1234	dcs12345	ecs12345678	fcs123	hcsA	hcsB	hemA-HpH	hemC-HpH	hemD-HpH	hemE-HpH	hemG-HpH	hemH-HpH	hemL-HpH	hemN-HpH	hemA-HpI	hemC-HpI	hemD-HpI	hemE-HpI	hemG-HpI	hemH-HpI	hemL-HpI	hemN-HpI	hemA-HH	hemB-HH	hemC-HH	hemD-HH	hemE-HH	hemG-HH	hemH-HH	hemL-HH	hemN-HH	hemA-HS	hemC-HS	hemD-HS	hemE-HS	hemG-HS	hemH-HS	hemL-HS	hemN-HS	hemA-HP	hemC-HP	hemD-HP	hemE-HP	hemG-HP	hemH-HP	hemL-HP	hemN-HP" >> $PATH_output/$Name_itolFile_Binary.txt
#DATA
echo "DATA" >> $PATH_output/$Name_itolFile_Binary.txt

for fullgenes in $PATH_output/*fullgenes*.txt
do
	SampleName=$(basename $fullgenes| cut -d '_' -f 1 | cut -d '.' -f 1) #Depends on your fastq file names!
	hemAHpH=-1
	hemCHpH=-1
	hemDHpH=-1
	hemEHpH=-1
	hemGHpH=-1
	hemHHpH=-1
	hemLHpH=-1
	hemNHpH=-1
	hemAHpI=-1
	hemCHpI=-1
	hemDHpI=-1
	hemEHpI=-1
	hemGHpI=-1
	hemHHpI=-1
	hemLHpI=-1
	hemNHpI=-1
	hemAHH=-1
	hemBHH=-1
	hemCHH=-1
	hemDHH=-1
	hemEHH=-1
	hemGHH=-1
	hemHHH=-1
	hemLHH=-1
	hemNHH=-1
	hemAHS=-1
	hemCHS=-1
	hemDHS=-1
	hemEHS=-1
	hemGHS=-1
	hemHHS=-1
	hemLHS=-1
	hemNHS=-1
	hemAHP=-1
	hemCHP=-1
	hemDHP=-1
	hemEHP=-1
	hemGHP=-1
	hemHHP=-1
	hemLHP=-1
	hemNHP=-1
	fklB=-1
	ndhI=-1
	tpd=-1
	resA=-1
	hypD=-1
	phoB=-1
	pdxT=-1
	dat=-1
	oppC=-1
	hxuB=-1
	nadV=-1
	bexA=-1
	bexB=-1
	bexC=-1
	bexD=-1
	acs1234=-1
	bcs1234=-1
	ccs1234=-1
	dcs12345=-1
	ecs12345678=-1
	fcs123=-1
	hcsA=-1
	hcsB=-1
	
	HS=0
	HP=0
	HpI=0
	HI=0
	HpH=0
	HH=0
	HD=0
	
	Subspecies=""
	Color_Subspecies=""
	Serotype=""
	Color_Serotype=""
	Species=""
	Color_Species=""
	while read line
	do
		Divergence=$(echo "[$line]" | cut -f 9 | cut -d '.' -f 1) #round down
#HpH haemin synthesis genes
		if [[ "$line" == *"hemA-HpH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemAHpH=1
		fi
		if [[ "$line" == *"hemC-HpH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemCHpH=1
		fi
		if [[ "$line" == *"hemE-HpH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemEHpH=1
		fi
		if [[ "$line" == *"hemG-HpH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemGHpH=1
		fi
		if [[ "$line" == *"hemH-HpH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemHHpH=1
		fi
		if [[ "$line" == *"hemL-HpH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemLHpH=1
		fi
		if [[ "$line" == *"hemN-HpH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemNHpH=1
		fi
		if [[ "$line" == *"hemD-HpH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemDHpH=1
		fi
#HpI haemin synthesis genes
		if [[ "$line" == *"hemA-HpI"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemAHpI=1
		fi
		if [[ "$line" == *"hemC-HpI"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemCHpI=1
		fi
		if [[ "$line" == *"hemE-HpI"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemEHpI=1
		fi
		if [[ "$line" == *"hemG-HpI"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemGHpI=1
		fi
		if [[ "$line" == *"hemH-HpI"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemHHpI=1
		fi
		if [[ "$line" == *"hemL-HpI"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemLHpI=1
		fi
		if [[ "$line" == *"hemN-HpI"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemNHpI=1
		fi
		if [[ "$line" == *"hemD-HpI"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemDHpI=1
		fi
#HH haemin synthesis
		if [[ "$line" == *"hemA-HH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemAHH=1
		fi
		if [[ "$line" == *"hemB-HH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemBHH=1
		fi
		if [[ "$line" == *"hemC-HH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemCHH=1
		fi
		if [[ "$line" == *"hemE-HH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemEHH=1
		fi
		if [[ "$line" == *"hemG-HH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemGHH=1
		fi
		if [[ "$line" == *"hemH-HH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemHHH=1
		fi
		if [[ "$line" == *"hemL-HH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemLHH=1
		fi
		if [[ "$line" == *"hemN-HH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemNHH=1
		fi
		if [[ "$line" == *"hemD-HH"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemDHH=1
		fi
#HP haemin synthesis genes
		if [[ "$line" == *"hemA-HP"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemAHP=1
		fi
		if [[ "$line" == *"hemC-HP"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemCHP=1
		fi
		if [[ "$line" == *"hemE-HP"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemEHP=1
		fi
		if [[ "$line" == *"hemG-HP"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemGHP=1
		fi
		if [[ "$line" == *"hemH-HP"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemHHP=1
		fi
		if [[ "$line" == *"hemL-HP"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemLHP=1
		fi
		if [[ "$line" == *"hemN-HP"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemNHP=1
		fi
		if [[ "$line" == *"hemD-HP"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemDHP=1
		fi
#HS haemin synthesis genes
		if [[ "$line" == *"hemA-HS"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemAHS=1
		fi
		if [[ "$line" == *"hemC-HS"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemCHS=1
		fi
		if [[ "$line" == *"hemE-HS"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemEHS=1
		fi
		if [[ "$line" == *"hemG-HS"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemGHS=1
		fi
		if [[ "$line" == *"hemH-HS"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemHHS=1
		fi
		if [[ "$line" == *"hemL-HS"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemLHS=1
		fi
		if [[ "$line" == *"hemN-HS"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemNHS=1
		fi
		if [[ "$line" == *"hemD-HS"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hemDHS=1
		fi
#HH marker genes
		if [[ "$line" == *"fklB"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			fklB=1
		fi
		if [[ "$line" == *"ndhI"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			ndhI=1
		fi
		if [[ "$line" == *"tpd"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			tpd=1
		fi
		if [[ "$line" == *"resA"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			resA=1
		fi
		if [[ "$line" == *"hypD"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hypD=1
		fi
##HI marker genes
		if [[ "$line" == *"dat"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			dat=1
		fi
		if [[ "$line" == *"pdxT"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			pdxT=1
		fi
		if [[ "$line" == *"phoB"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			phoB=1
		fi
		if [[ "$line" == *"oppC"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			oppC=1
		fi
		if [[ "$line" == *"hxuB"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hxuB=1
		fi
#HD
		if [[ "$line" == *"nadV"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			nadV=1
		fi
#Capsule loci
		if [[ "$line" == *"bexA"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			bexA=1
		fi
		
		if [[ "$line" == *"bexB"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			bexB=1
		fi
		
		if [[ "$line" == *"bexC"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			bexC=1
		fi
		
		if [[ "$line" == *"bexD"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			bexD=1
		fi
		
		if [[ "$line" == *"acs1234"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			acs1234=1
		fi
		
		if [[ "$line" == *"bcs1234"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			bcs1234=1
		fi
		
		if [[ "$line" == *"ccs1234"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			ccs1234=1
		fi
		
		if [[ "$line" == *"dcs12345"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			dcs12345=1
		fi
		
		if [[ "$line" == *"ecs12345678"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			ecs12345678=1
		fi
		
		if [[ "$line" == *"fcs123"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			fcs123=1
		fi
		
		if [[ "$line" == *"hcsA"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hcsA=1
		fi
		
		if [[ "$line" == *"hcsB"* ]] && [[ "$Divergence" -lt "$MaxDiv" ]]
		then
			hcsB=1
		fi
	done <<< "$(cat $fullgenes)"

##Decision rules
#Species determination
	HpH_calc=$(($hemAHpH+$hemCHpH+$hemDHpH+$hemEHpH+$hemGHpH+$hemHHpH+$hemLHpH+$hemNHpH))
	if [[ "$HpH_calc" -ge 6 ]]
	then
		HpH=1
	fi
	HpI_calc=$(($hemAHpI+$hemCHpI+$hemDHpI+$hemEHpI+$hemGHpI+$hemHHpI+$hemLHpI+$hemNHpI))
	if [[ "$HpI_calc" -ge 6 ]]
	then
		HpI=1
	fi
	HP_calc=$(($hemAHP+$hemCHP+$hemDHP+$hemEHP+$hemGHP+$hemHHP+$hemLHP+$hemNHP))
	if [[ "$HP_calc" -ge 6 ]]
	then
		HP=1
	fi
	HS_calc=$(($hemAHS+$hemCHS+$hemDHS+$hemEHS+$hemGHS+$hemHHS+$hemLHS+$hemNHS))
	if [[ "$HS_calc" -ge 6 ]]
	then
		HS=1
	fi
	HH_calc=$(($hypD+$fklB+$ndhI+$tpd+$resA))
	if [[ "$HH_calc" -ge 3 ]]
	then
		HH=1
	fi
	HI_calc=$(($phoB+$pdxT+$dat+$hxuB+$oppC))
	if [[ "$HI_calc" -ge 3 ]]
	then
		HI=1
	fi
	
	if [[ "$nadV" -eq 1 ]]
	then
		HD=1
	fi
	H_calc=$(($HI+$HH+$HD+$HpH+$HpI+$HS+$HP))
	if [[ "$H_calc" -ge 2 ]]
	then
		Species="Mix"
		Color_Species="#fdd5b1"
	elif [[ "$H_calc" -eq 0 ]]
	then 
		Species="Unknown"
		Color_Species="#dbd7d2"
	elif [[ "$HD" -eq 1 ]]
	then
		Species="HD"
		Color_Species="#808080"
	elif [[ "$HI" -eq 1 ]]
	then
		Species="HI"
		Color_Species="#4cbb17"
	elif [[ "$HH" -eq 1 ]]
	then
		Species="HH"
		Color_Species="#9bddff"
	elif [[ "$HpH" -eq 1 ]]
	then
		Species="HpH"
		Color_Species="#ae0c00"
	elif [[ "$HpI" -eq 1 ]]
	then
		Species="HpI"
		Color_Species="#eee600"
	elif [[ "$HP" -eq 1 ]]
	then
		Species="HP"
		Color_Species="#ff00ff"
	elif [[ "$HS" -eq 1 ]]
	then
		Species="HS"
		Color_Species="#ff9933"
	fi
#Subspecies determination
	HHi_calc=$(($hemAHH+$hemBHH+$hemCHH+$hemDHH+$hemEHH+$hemGHH+$hemHHH+$hemLHH+$hemNHH))
	if [[ "$HHi_calc" -eq 9 ]]
	then
		Subspecies="subsp. intermedius"
		Color_Subspecies="#93ccea"
	fi
#In silico Serotyping
	Serotype_calc=$(($acs1234+$bcs1234+$ccs1234+$dcs12345+$ecs12345678+$fcs123))
	Serotype_calc2=$(($bexA+$bexB+$bexC+$bexD+$hcsA+$hcsB))
	Serotype_calc3=$(($Serotype_calc+$Serotype_calc2))
	if [[ "$Serotype_calc" -ge -2 ]]
	then
		Serotype="Multiple serotypes"
		Color_Serotype="#fdd5b1"
	elif [[ "$Serotype_calc3" -lt 2 ]] && [[ "$Serotype_calc3" -ge -10 ]]
	then
		Serotype="Capsule-deficient"
		Color_Serotype="#dbd7d2"
	elif [[ "$Serotype_calc3" -eq -12 ]] && [[ "$HI" -eq 1 ]]
	then
		Serotype="NTHi"
		Color_Serotype="#006600" 
	elif [[ "$Serotype_calc3" -eq 2 ]]
	then
		if [[ "$acs1234" -eq 1 ]]
		then
			Serotype="Serotype a"
			Color_Serotype="#000000"
		elif [[ "$bcs1234" -eq 1 ]]
		then
			Serotype="Serotype b"
			Color_Serotype="#fd0e35"
		elif [[ "$ccs1234" -eq 1 ]]
		then
			Serotype="Serotype c"
			Color_Serotype="#800080"
		elif [[ "$dcs12345" -eq 1 ]]
		then
			Serotype="Serotype d"
			Color_Serotype="#eee600"
		elif [[ "$ecs12345678" -eq 1 ]]
		then
			Serotype="Serotype e"
			Color_Serotype="#d2691e"
		elif [[ "$fcs123" -eq 1 ]]
		then
			Serotype="Serotype f"
			Color_Serotype="#808080"
		fi
	fi
#Print results to decision and itol files
	echo $SampleName"	"$Species"	"$Subspecies"	"$Serotype >> $PATH_output/$Name_DecisionFile.txt
	echo $SampleName"	"$Color_Species >> $PATH_output/$Name_itolFile_Decision_Species.txt
	echo $SampleName"	"$Color_Serotype >> $PATH_output/$Name_itolFile_Decision_Serotype.txt
	echo $SampleName"	"$Color_Subspecies >> $PATH_output/$Name_itolFile_Decision_Subspecies.txt
	echo $SampleName"	"$nadV"	"$hypD"	"$fklB"	"$ndhI"	"$tpd"	"$resA"	"$phoB"	"$pdxT"	"$dat"	"$oppC"	"$hxuB"	"$bexA"	"$bexB"	"$bexC"	"$bexD"	"$acs1234"	"$bcs1234"	"$ccs1234"	"$dcs12345"	"$ecs12345678"	"$fcs123"	"$hcsA"	"$hcsB"	"$hemAHpH"	"$hemCHpH"	"$hemDHpH"	"$hemEHpH"	"$hemGHpH"	"$hemHHpH"	"$hemLHpH"	"$hemNHpH"	"$hemAHpI"	"$hemCHpI"	"$hemDHpI"	"$hemEHpI"	"$hemGHpI"	"$hemHHpI"	"$hemLHpI"	"$hemNHpI"	"$hemAHH"	"$hemBHH"	"$hemCHH"	"$hemDHH"	"$hemEHH"	"$hemGHH"	"$hemHHH"	"$hemLHH"	"$hemNHH"	"$hemAHS"	"$hemCHS"	"$hemDHS"	"$hemEHS"	"$hemGHS"	"$hemHHS"	"$hemLHS"	"$hemNHS"	"$hemAHP"	"$hemCHP"	"$hemDHP"	"$hemEHP"	"$hemGHP"	"$hemHHP"	"$hemLHP"	"$hemNHP >> $PATH_output/$Name_itolFile_Binary.txt
done

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
