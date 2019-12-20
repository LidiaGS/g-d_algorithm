#!/bin/bash

# submit_pipeline.sh
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

blue=$( tput setaf 4 )
white=$( tput setaf 7)

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## EXPORT PATHS AND INPUT DATA
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Notice that paths must be adapted to your personal setup. 
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Tools
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TRIMMOMATIC_PATH=$3 # /path/to/Trimmomatic ## WRITE YOUR OWN PATH HERE!
BWA_PATH=$4 # /path/to/BWA ## WRITE YOUR OWN PATH HERE!
ST_PATH=$5 # /path/to/SAMtools ## WRITE YOUR OWN PATH HERE!
gd_PATH=$(dirname $0)

# Input data
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
REF_PATH=$1 ##/path/to/references  ## WRITE YOUR OWN PATH HERE!
INDEX_PATH=${REF_PATH}/indexes 
R1=$2 ##/path/to/single-end-reads-file.fastq ## WRITE YOUR OWN SAMPLE HERE!

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
java -jar ${TRIMMOMATIC_PATH}/trimmomatic-0.36.jar SE -threads 12 -phred33 ${R1} ${OUTPUT_PATH}/filtered_reads/filtered_${R1##*/} MINLEN:140 CROP:150

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## MAPPING SAMPLE TO REFERENCES
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
echo "Mapping starts at: `date +%Y/%m/%d-%H:%M:%S`";
for ref in ${REF_PATH}/*.fna; do
	out_ref=${ref##*/}
	${BWA_PATH}/bwa mem -t 12 ${INDEX_PATH}/${out_ref} ${OUTPUT_PATH}/filtered_reads/filtered_${R1##*/} > ${OUTPUT_PATH}/mapped_reads/${out_ref%%.*}.sam;
done

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## FILTER SAM FILE
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Remove supplementary alignments (0x800), not primary alignments (0x100) and not mapped reads (0x4).
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
echo "Filtering the mapped reads starts at: `date +%Y/%m/%d-%H:%M:%S`";
for sam in ${OUTPUT_PATH}/mapped_reads/*; do
	out_ref=${sam##*/}
	${ST_PATH}/samtools view -@ 12 -S -F2308 ${OUTPUT_PATH}/mapped_reads/${out_ref%%.*}.sam -o ${OUTPUT_PATH}/filtered_reads/${out_ref%%.*}_F2308.sam;
done

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## ASSING READS TO A SINGLE SPECIES
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  Using the gamma-delta algorithm
##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
echo "Reads assignment starts at: `date +%Y/%m/%d-%H:%M:%S`";
python ${gd_PATH}/g-d_algorithm_script.py -g 0.99 -d 0.98 -m ${OUTPUT_PATH}/filtered_map -s ${R1} -o ${OUTPUT_PATH}/g-d_assignment.csv;

echo "Pipeline ends at: `date +%Y/%m/%d-%H:%M:%S`";
