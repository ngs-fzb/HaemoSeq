#!/usr/bin/env bash

###Usage###
# bash $PathToScript/Scriptname.sh -h

###Author###
#echo "hello World, this script was written by Margo Diricks (mdiricks@fz-borstel.de)!"

###Version###
version="1.0.0"

###Function###
#To convert the custom plasmid database file into a SRST2-compatible database file

###Required packages###
#None

###Required parameters - via command line###
#-d db_custom_SRST2=""

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Syntax: NameScript.sh [parameters]"
   echo "Required parameters:"
   echo "d     Full path to plasmid database ($PATH/custom.fasta)"
   echo "v     Display version"
   echo "h     Display help"
}
############################################################
# Get parameters                                                #
############################################################

while getopts ":ho:d::v" option; do #:h does not need argument, f: does need argument
   case $option in
      h) # display Help
         Help
         exit;;
      d) # 
         db_custom_SRST2=$OPTARG;;
      v) # display Version
         echo $version
         exit;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done

###Check if required parameters are provided###
if [[ -z "db_custom_SRST2" ]]
then
	echo "Please provide all required arguments (-d db_custom_SRST2)! use custom_to_SRST2_db.sh -h for help on syntax"
	exit
fi

###############################################################################CODE#################################################################################
customDB=$(basename $db_custom_SRST2| cut -d '.' -f 1)
customDB_SRST2=$customDB"_SRST2.fasta"
COUNTER1=0
cat $db_custom_SRST2 | while read line
do
	if [[ "$line" == *">"* ]]
	then
		COUNTER1=$((COUNTER1 + 1))
		Part1=$(echo $line | cut -d ' ' -f 1 | cut -d '>' -f 2)
		#Part2=$(echo $line | cut -d ' ' -f 2)
		echo ">"$COUNTER1"__"$Part1"__"$Part1"__"$COUNTER1 >> $(dirname $db_custom_SRST2)/$customDB_SRST2
	else
		echo $line >> $(dirname $db_custom_SRST2)/$customDB_SRST2
	fi
done
