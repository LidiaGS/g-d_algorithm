#!/bin/bash

# gamma-delta-workflow.sh
#  
# Pipeline for the implementation of the gamma-delta algorithm 
# This file is part of the https://github.com/LidiaGS/g-d_algorithm.
# Copyright (c) 2019 Universitat Autònoma de Barcelona
# 
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
# 
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## FUNCTIONS
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
PROGNAME=$(basename $0)

error_no_exit()
{
	echo ${blue}"${PROGNAME}: $1"${white} 1>&2
}

helpFunction()
{
   echo "-----------------------------------------------------------------"
   echo $red"Usage: "$white 
   echo $blue"       $0 <single-end.fastq> "$white "(Single-end reads in FASTA/FASTQ format)"
   echo $blue"       $0 <forward-R1.fastq> <reverse-R2.fastq> "$white "(Paired-end reads in FASTA/FASTQ format)"
   exit 1 # Exit script after printing help
}

red=$( tput setaf 1 )
blue=$( tput setaf 4 )
green=$( tput setaf 2 )
white=$( tput setaf 7)

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## EXPORT PATHS AND INPUT DATA
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Notice that paths must be adapted to your personal setup. 
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Tools
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TRIMMOMATIC_PATH=/path/to/Trimmomatic ## WRITE YOUR OWN PATH HERE!
BWA_PATH=/path/to/BWA ## WRITE YOUR OWN PATH HERE!
ST_PATH=/path/to/SAMtools ## WRITE YOUR OWN PATH HERE!
gd_PATH=/path/to/gamma-delta_algorithm_script ## WRITE YOUR OWN PATH HERE!

# Input data
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
REF_PATH=/path/to/references  ## WRITE YOUR OWN PATH HERE!
INDEX_PATH=${REF_PATH}/indexes 


if [ -z ${1} ] && [ -z ${2} ]; then 
	helpFunction;
	exit 0;
else
if [-z ${2} ]; then 
	R1=$1
	echo "Processing single-end reads file: "; 
	echo "> $R1"; 
else
	R2=$2
	echo "Processing paired-end reads files: "; 
	echo "> R1: ${R1} ";
	echo "> R2: ${R2} "; 
fi
fi

# Output data
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
CWD=$(pwd) ## Output folder will be created in the same folder than the one were the user is.
output_folder=${R1##*/}
OUTPUT_PATH=${CWD}/${output_folder%%.*}

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## CREATE FOLDERS
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mkdir ${INDEX_PATH} || error_no_exit "Directory already exist! May be overwritting..."
mkdir ${OUTPUT_PATH} || error_no_exit "Directory already exist! May be overwritting..."
mkdir ${OUTPUT_PATH}/filtered_reads || error_no_exit "Directory already exist! May be overwritting..."
mkdir ${OUTPUT_PATH}/mapped_reads || error_no_exit "Directory already exist! May be overwritting..."
mkdir ${OUTPUT_PATH}/filtered_map || error_no_exit "Directory already exist! May be overwritting..."

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## CREATE REFERENCES' INDEXES
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Notice that this step should only be performed once!
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
echo "Indexing starts at: `date +%Y/%m/%d-%H:%M:%S`";
for ref in ${REF_PATH}/*.fna; do
 	echo ">" $ref
 	${BWA_PATH}/bwa index ${ref};
 	mv ${ref}.* ${INDEX_PATH}/.
done

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## FILTER SAMPLE DATA 
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Trim at 150bp length and remove reads shorter than 140bp while keeping paired-end reads
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
echo "Filtering samples starts at: `date +%Y/%m/%d-%H:%M:%S`";

if [[ -z ${R2} ]]; then  
	echo "HERE /1"
	## Single-end case:
	java -jar ${TRIMMOMATIC_PATH}/trimmomatic-0.36.jar SE -threads 24 -phred33 ${R1} ${OUTPUT_PATH}/filtered_reads/filtered_${R1##*/} MINLEN:140 CROP:150
else 
	echo "HERE /2"
	## Paired-end case:
	java -jar ${TRIMMOMATIC_PATH}/trimmomatic-0.36.jar PE -threads 24 -phred33 ${R1} ${R2} ${OUTPUT_PATH}/filtered_reads/paired_${R1##*/} ${OUTPUT_PATH}/filtered_reads/unpaired_${R1##*/} ${OUTPUT_PATH}/filtered_reads/paired_${R2##*/} ${OUTPUT_PATH}/filtered_reads/unpaired_${R2##*/} MINLEN:140 CROP:150
fi

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## MAPPING SAMPLE TO REFERENCES
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
echo "Mapping starts at: `date +%Y/%m/%d-%H:%M:%S`";
if [[ -z ${R2} ]]; then  
	## Single-end case:
	for ref in ${REF_PATH}/*.fna; do
		out_ref=${ref##*/}
		${BWA_PATH}/bwa mem -t 24 ${INDEX_PATH}/${out_ref} ${OUTPUT_PATH}/filtered_reads/filtered_${R1##*/} > ${OUTPUT_PATH}/mapped_reads/${out_ref%%.*}.sam;
	done
else 
	## Paired-end case:
	for ref in ${REF_PATH}/*.fna; do
		out_ref=${ref##*/}
		${BWA_PATH}/bwa mem -t 24 ${INDEX_PATH}/${out_ref} ${OUTPUT_PATH}/filtered_reads/paired__${R1##*/} ${OUTPUT_PATH}/filtered_reads/paired__${R2##*/} > ${OUTPUT_PATH}/mapped_reads/${out_ref%%.*}_PE.sam;
	done
fi

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## FILTER SAM FILE
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Remove supplementary alignments (0x800), not primary alignments (0x100) and not mapped reads (0x4).
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
echo "Filtering the mapped reads starts at: `date +%Y/%m/%d-%H:%M:%S`";
if [[ -z ${R2} ]]; then   
	## Single-end case:
	for sam in ${OUTPUT_PATH}/mapped_reads/*; do
		out_ref=${sam##*/}
		${ST_PATH}/samtools view -@ 24 -S -F2308 ${OUTPUT_PATH}/mapped_reads/${out_ref%%.*}.sam -o ${OUTPUT_PATH}/filtered_map/${out_ref%%.*}.sam;
	done
else 
	## Paired-end case:
	for sam in ${OUTPUT_PATH}/mapped_reads/*; do
		out_ref=${sam##*/}
		${ST_PATH}/samtools view -@ 24 -S -F2308 ${OUTPUT_PATH}/mapped_reads/${out_ref%%.*}_PE.sam -o ${OUTPUT_PATH}/filtered_map/${out_ref%%.*}.sam;
	done
fi

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## ASSING READS TO A SINGLE SPECIES
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Using the gamma-delta algorithm with values of gamma=0.99 and delta=0.98
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
echo "Reads assignment starts at: `date +%Y/%m/%d-%H:%M:%S`";
if [[ -z ${R2} ]]; then  
	## Single-end case:
	python ${gd_PATH}/gamma-delta.py -g 0.99 -d 0.98 -m ${OUTPUT_PATH}/filtered_map -r ${R1} -o ${OUTPUT_PATH}/g-d_assignment.csv;
else 
	## Paired-end case:
	python ${gd_PATH}/gamma-delta.py -g 0.99 -d 0.98 -m ${OUTPUT_PATH}/filtered_map -r1 ${R1} -r2 ${R2} -o ${OUTPUT_PATH}/g-d_assignment.csv;
fi

echo "Pipeline ends at: `date +%Y/%m/%d-%H:%M:%S`";
