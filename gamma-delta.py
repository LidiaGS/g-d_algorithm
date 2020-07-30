# gamma-delta.py 
# coding=utf-8
#
# Implementation of the gamma-delta algorithm
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

import csv
import argparse
import os
from os.path import isfile, join
import operator
from decimal import Decimal
import datetime

## ARGUMENTS DEFINITION
parser = argparse.ArgumentParser(prog='ARGUMENTS', usage='%(prog)s [options]')
parser.add_argument("-r",   "--R",              type=str,   help="Single-end read sample (FASTQ/FASTA)")
parser.add_argument("-r1",  "--R1",             type=str,   help="Separeted paired-end reads sample, R1 (forward) file (FASTA/FASTQ format)")
parser.add_argument("-r2",  "--R2",             type=str,   help="Separeted paired-end reads sample, R2 (reverse) file (FASTA/FASTQ format)")
parser.add_argument("-m",   "--SAMfolder",      type=str,   help="Folder containing the mapped reads files in SAM format (e.g: mapped_reads.sam).")
parser.add_argument("-g",   "--gamma",          type=str,   help="Gamma threshold (per one unit) [1-0]")
parser.add_argument("-d",   "--delta",          type=str,   help="Delta threshold (per one unit) [1-0]")
parser.add_argument("-o",   "--output",         type=str,   help="Name of the output file (CSV format)")
parser.add_argument("-O",   "--outInfo",        type=str,   help="Output the assignment of each read (CSV format)")
args = parser.parse_args()

## ABREVIATIONS
# A1:  Highest mapping ratio
# A2:  Second highest mapping ratio

## FUNCTIONS 

## Obtaining the mapping ratio (A) read against a reference by using CIGAR and NM from SAM file
def mapping_ratio(read_CIGAR, read_NM):
    num_cigar = 0
    match_mismatch = 0
    insertion = 0
    deletion = 0
    soft_clipped = 0
    hard_clipped = 0
    padding = 0
    skipped = 0
    mismatch = 0
    match = 0 
    NM_mismatch = 0
    total_cigar = 0
    map_ratio = 0

    ## Trim CIGAR into parts and join same values: 
    for cigar in read_CIGAR: # Looping through elements of CIGAR
        if cigar.isdigit(): # Save numeric values from CIGAR
            num_cigar = num_cigar * 10
            num_cigar = int(cigar) + num_cigar                
        else: # Save numeric values from CIGAR to its cathegory (i.e.: Matching/mismatchig, indels, ...)
            if cigar == 'M': ## M: matches & mismatches
                match_mismatch = num_cigar + match_mismatch
                num_cigar = 0
            elif cigar == 'I': ## I: insertion 
                insertion = num_cigar + insertion
                num_cigar = 0                
            elif cigar == 'D': ## D: deletion
                deletion = num_cigar + deletion
                num_cigar = 0
            elif cigar == 'S': ## S: soft-clipped
                soft_clipped = num_cigar + soft_clipped
                num_cigar = 0
            elif cigar == 'H': ## H: hard-clipped
                hard_clipped = num_cigar + hard_clipped
                num_cigar = 0
            elif cigar == 'P': ## P: Padding
                padding = num_cigar + padding
                num_cigar = 0
            elif cigar == 'N': ## N: Skipped
                skipped = num_cigar + skipped
                num_cigar = 0
            elif cigar == 'X': ## X: Mismatch
                mismatch = num_cigar + mismatch
                num_cigar = 0
            elif cigar == '=': ## =: Match
                match = num_cigar + match
                num_cigar = 0
            elif cigar == '*':
                num_cigar = 0
    
    ## Length of CIGAR (also reads length)
    total_cigar = match_mismatch + insertion + deletion + soft_clipped + hard_clipped + padding + skipped + mismatch + match

    ## Mapping ratio
    map_ratio = (((total_cigar + deletion + skipped) - float(read_NM) - soft_clipped - hard_clipped - mismatch) / (total_cigar + deletion + skipped)) 
    return map_ratio

## Looking for NM information thought columns of SAM rows
def findNM(samLine):
    NM = str()
    for elem in samLine:
        if str(elem[:2])=='NM':
            return elem[5:] ## Keeping NM tag value without "NM:i:"
    return 0

## When working with paried-end reads, add "/1" and "/2" at the end of the R1 and R2 names, respectively
def rename_paired_end(readName, flag):
    if ((int(flag) & 1) == 1): # # Check paired-end read; if true, this are paired-end reads
        if ((int(flag) & 128) == 128): # Check second in pair; if true, return R2
            return(readName+str("/2"))
        else: # R1
            return(readName+str("/1"))
    else: # Single-end read
        return(readName+str(""))

## Check that A1 is above gamma and A2 below delta
def check_distance(ReadDict):
    for ReadId in ReadDict:
        if (float(ReadDict[ReadId][5]) >= float(args.delta)):
            ReadDict[ReadId][9] = "A2-above-delta"
        elif ((float(ReadDict[ReadId][0]) < float(args.gamma))):
            ReadDict[ReadId][9] = "A1-below-gamma"
        else:
            ReadDict[ReadId][9] = "Assigned"
    return ReadDict

## Saving the mapping information for a given read and reference into the dictionary
def save_map_info(mapRow, mapRef, ReadDict):
    mapRow = mapRow.strip('\n')
    row = mapRow.split("\t")
    ReadId =  str(rename_paired_end(row[0],row[1])) # RNAME
    
    # Create read in the dictionary if the key does not exist yet
    if (ReadId not in ReadDict): 
        ReadDict[ReadId] = list([0, '0', 'NA', 'NA', 'NA', 0, 'NA', 'NA', '0', 'NA']) # Create and empty tuple, so A1 and A2 are zero.
    
    # Obtain the mapping ratio of the read in the current row 
    newRatio = mapping_ratio(row[5], findNM(row)) # CIGAR, NM are picked
    
    # New mapping ratio is higher or equal to A1
    if (newRatio >= ReadDict[ReadId][0]): 
        # Changing A2
        ReadDict[ReadId][5] = ReadDict[ReadId][0] # A2 (A2 to old A1)
        ReadDict[ReadId][6] = ReadDict[ReadId][1] # Reference id of A2 (A2 to old A1)
        ReadDict[ReadId][7] = ReadDict[ReadId][2] # RNAME A2 (A2 to old A1)
        ReadDict[ReadId][8] = ReadDict[ReadId][3] # POS A2 (A2 to old A1)
        # Changing A1
        ReadDict[ReadId][0] = newRatio # A1 (old A1 to new A1)
        ReadDict[ReadId][1] = mapRef # Reference id of A1 (old A1 to new A1)
        ReadDict[ReadId][2] = row[2] # RNAME A1 (old A1 to new A1)
        ReadDict[ReadId][3] = row[3] # POS A1 (old A1 to new A1)
        ReadDict[ReadId][4] = row[1] # FLAG A1 (old A1 to new A1)
    # New mapping ratio is higher than A2    
    elif (newRatio > ReadDict[ReadId][5]):
        ReadDict[ReadId][5] = newRatio # A2 (old A2 to new A1)
        ReadDict[ReadId][6] = mapRef # Reference id of A2 (old A2 to new A1)
        ReadDict[ReadId][7] = row[2] # RNAME A2 (old A2 to new A1)
        ReadDict[ReadId][8] = row[3] # POS A2 (old A2 to new A1)
    return ReadDict

## Read and load SAM files. Notice that SAM header is ignored.
def load_sam_data(samDocs):
    ReadDict = dict()
    print "Loading mapping information..."
    for samfile in samDocs: # Loop through files in folder
        refName = samfile.split(".")[0] # The ref code is at the file name
        sam = open(samfile, "r")
        for row in sam:  # Loop through lines in sam file
            if row[0][0] != '@': # Avoiding headers
                ReadDict = save_map_info(row, refName, ReadDict)
        sam.close()
    return ReadDict

## Create a dictionary for storing the number of reads assigned to each reference
def summary_dict(ReadDict):
    SumDic = dict()
    SumDic["A1-below-gamma"] = 0
    SumDic["A2-above-delta"] = 0

    for ReadId in ReadDict:
        ## Check if the read has been assigned 
        if ReadDict[ReadId][9] == "Assigned": 
            if ReadDict[ReadId][1] in SumDic: # Reference already exists
                SumDic[ReadDict[ReadId][1]] += 1
            else: # Reference doesn't exist yet
                newkey = ReadDict[ReadId][1]
                SumDic[newkey] = 1
        ## If the read has not been assigned, save the reason
        elif ReadDict[ReadId][9] == "A1-below-gamma":
            SumDic["A1-below-gamma"] += 1
        elif ReadDict[ReadId][9] == "A2-above-delta":
            SumDic["A2-above-delta"] += 1
    return SumDic

# Pop value "valueKey" from dictionary
def popValue(valueKey, ReadDict):
    try:
        pop_elem = ReadDict.pop(valueKey)
        return (pop_elem, ReadDict) # Return the value associated to the key and the Dictionary without the element
    except KeyError:
        return (0, ReadDict)

# Obtain header column name
def sample_header():
    if (args.R is not None):
        return os.path.basename(args.R)
    elif ((args.R1 is not None) and (args.R2 is not None)):
        return str(os.path.basename(args.R1))+" + "+str(os.path.basename(args.R2))
    elif (args.SAMfolder is not None):
        return args.SAMfolder
    elif (args.SAMfile is not None):
        return args.SAMfile
    else:
        return str("Sample:")+str(datetime.datetime.now().time())

# Count the number of reads in file (FASTQ/FASTA formats)
def count_reads_in_file(fileName):
    if(fileName[-1:]=="a"): # FASTA file
        j = 2
    else: # FASTQ file
        j = 4   
    return(len(open(fileName).readlines( ))/j)

# Return the number of NOT mapped reads
def count_notMappedReads(mappedReads):
    if (args.R is not None):
        totalReads = count_reads_in_file(args.R)
        print ("> Initial reads: "+str(totalReads))
        return(totalReads -int(mappedReads))
    elif ((args.R1 is not None) and (args.R2 is not None)):
        totalReads = count_reads_in_file(args.R1)*2
        print ("> Initial reads: "+str(totalReads))
        return(totalReads -int(mappedReads))
    return(0)

## Save on a CSV file the summary information of the detected species together with the total number of reads assigned and their relative species abundance 
def save_map_info_csv(OldDict):
    totalMapp = sum(OldDict.values())
    totalNotMapp = count_notMappedReads(totalMapp)
    totalReads = totalMapp + totalNotMapp
    print ("> Mapped reads: "+str(totalMapp))

    # Pop top two elements from column
    rm_by_A1, OldDict = popValue('A1-below-gamma', OldDict)
    rm_by_A2, OldDict = popValue('A2-above-delta', OldDict)

    totalAssig = sum(OldDict.values())
    print ("> "+str(totalAssig)+" reads have been assigned")
    
    ## Sort Dictionary by VALUES        
    sorted_Tupla = sorted(OldDict.items(), key=operator.itemgetter(1), reverse=True)
    print sorted_Tupla

    ## Checking if output file already exists, if so, adding a column to the existing file
    exists = os.path.isfile(args.output)   
    if exists:
        # Store configuration file values
        csvinput = open(args.output,'r')
        reader = csv.reader(csvinput)
        lines = csvinput.readlines()
        csvoutput = open(os.path.splitext(args.output)[0]+"AUXoutput.csv", 'w')
        writer = csv.writer(csvoutput, delimiter=";")

        # Write header
        outRow = list(lines[0].rstrip("\r\n").split(";"))
        outRow.append((str(sample_header())))
        writer.writerow(outRow)
        
        # Write Not-mapped, A1-below-gamma and A2-above-delta
        outRow = list(lines[1].rstrip("\r\n").split(";"))
        outRow.append(("Not-Map"+" ("+str(totalNotMapp)+")"))
        writer.writerow(outRow)

        outRow = list(lines[2].rstrip("\r\n").split(";"))
        outRow.append(('A1-below-gamma'+" ("+str(rm_by_A1)+")"))
        writer.writerow(outRow)

        outRow = list(lines[3].rstrip("\r\n").split(";"))
        outRow.append(('A2-above-delta'+" ("+str(rm_by_A2)+")"))
        writer.writerow(outRow)

        # Write assignments
        for i in range(max(len(OldDict),len(lines)-4)): # Loop through lines in max(file/dictionary)
            # Add data from the file
            try:
                outRow = list(lines[i+4].rstrip("\r\n").split(";")) # +4 for avoiding the first 4 lines
            except IndexError: # If the new data has more lines than the old data
                outRow = list([0]*len(lines[0].rstrip("\r\n").split(";")))
            
            # Add new column into an existing file 
            try:
                outRow.append(sorted_Tupla[i][0]+" ("+str(sorted_Tupla[i][1])+"|"+str(round(Decimal(float(sorted_Tupla[i][1])/float(totalAssig)),5))+")")
            except IndexError: # New data has less lines than old data
                outRow.append(0)
            writer.writerow(outRow)
        os.rename(os.path.splitext(args.output)[0]+"AUXoutput.csv", args.output)
    
    # If output file doesn't exist yet, create it and save the mapping information
    else:
        csvoutput = open(args.output, 'w')
        writer = csv.writer(csvoutput)
        writer.writerow(list([str(sample_header())])) # Write sample name
        writer.writerow([str("Not-mapping-reads"+" ("+str(totalNotMapp)+")")]) # Write Not-Map
        writer.writerow([str('A1-below-gamma'+" ("+str(rm_by_A1)+")")]) # Write 'A1-below-gamma'
        writer.writerow([str('A2-above-delta'+" ("+str(rm_by_A2)+")")]) # Write 'A2-above-delta'

        for item in sorted_Tupla:
            writer.writerow([str(item[0]+" ("+str(item[1])+"|"+str(round(Decimal(float(item[1])/float(totalAssig)),5))+")")])
    csvoutput.close()
    return 

# Save on a CSV file the assignment information of each single read
def save_dic(ReadDict):
    csv_columns = ["RNAME","A1","REF_A1","SCAF_A1","POS_A1","FLAG_A1","A2","ref_A2","SCAF_A2","POS_A2","ASSIGNMENT"]
    with open(args.outInfo, 'w') as csvoutput:
        writer = csv.writer(csvoutput, delimiter=";")
        writer.writerow(csv_columns)
        for key in ReadDict.keys():
            o = ReadDict[key]
            o.insert(0,key)
            writer.writerow(o)
    csvoutput.close()     

## MAIN PROGRAMME
def main():

    # Saving mapping information (SAM files)
    wd = os.getcwd() # Get working directory
    if (args.SAMfolder  is not None): # Getting all the SAM files from folder
        os.chdir(args.SAMfolder)
        samfiles = [f for f in os.listdir(args.SAMfolder) if f.endswith(".sam")]
    else:
        samfol = raw_input("Where are the SAM files?:")
        os.chdir(samfol)
        samfiles = [f for f in os.listdir(samfol) if f.endswith(".sam")]
    ReadDict = load_sam_data(samfiles)

    ## Check that the A1 and A2 fulfill the gamma and delta thresholds requirements (e.g. A1 above gamma and A2 below delta)
    ReadDict = check_distance(ReadDict)
    
    print("Summarizing data...")
    os.chdir(wd) # Turn back to the original working directory
    summary = summary_dict(ReadDict)
    save_map_info_csv(summary)

    ## Save each read assingment information
    if (args.outInfo is not None):
        save_dic(ReadDict)
    return

## MAIN 
main()
