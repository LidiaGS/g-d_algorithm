import csv
import argparse
import os
from os.path import isfile, join
import operator
from decimal import Decimal

## ARGUMENTS DEFINITION
parser = argparse.ArgumentParser(prog='ARGUMENTS', usage='%(prog)s [options]')
parser.add_argument("-s",  "--sample",         type=str,   help="Name of sampled data in csv format (e.g.: samples_info.csv)")
parser.add_argument("-r1",  "--R1",            type=str,   help="R1 Sample")
parser.add_argument("-r2",  "--R2",            type=str,   help="R2 Sample")
parser.add_argument("-m",  "--SAMfolder",      type=str,   help="Name of mapped reads file in sam format (e.g: refPM_readsPM.sam). It needs -n argument.")
parser.add_argument("-g",  "--gamma",          type=str,   help="Minimum percentage")
parser.add_argument("-d",  "--delta",          type=str,   help="Maximum percentage")
parser.add_argument("-o",  "--output",         type=str,   help="Output file")
args = parser.parse_args()


## Abreviations
# HMP: Higher Mapping Percentage
# NI:  Non-informative
# A1:  First highest mapping percentage
# A2:  Second highest mapping percentage

## FUNCTIONS 

## Obtaining the mapping percentage read against a reference by using CIGAR and NM from SAM file
def percentage_alignment(read_CIGAR, read_NM, read_len):
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
    aln_percentage = 0

    ## Trim CIGAR into parts and join same values: 
    try: 
        for cigar in read_CIGAR: ## Looping through elements of CIGAR
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
    except ValueError:
        print("Error: ON CIGAR")
    
    ## Length of CIGAR
    total_cigar = match_mismatch + insertion + deletion + soft_clipped + hard_clipped + padding + skipped + mismatch + match

    ## Mapping percentage
    aln_percentage = (((read_len + deletion + skipped) - float(read_NM) - soft_clipped - hard_clipped - mismatch) / (read_len + deletion + skipped)) 
    return aln_percentage

## Looking for NM information thought columns of SAM rows
def findNM(samLine):
    NM = str()
    for elem in samLine:
        if str(elem[:2])=='NM':
            return elem[5:] ## Keeping NM tag value without "NM:i:"
    return 0

## Create a dictionary using read's names as key and adding an empty list of three elements (best hit ref., first HMP and second best hit ref.) to each dictionary element. 
def load_reads_data(sampleDoc): 
    print 'Loading sample data...'
    reads = open(sampleDoc, "r") 
    
    ## Check if input data is in FASTA or FASTQ format
    if(sampleDoc[-1:]=="a"): # FASTA file
        j = 2
    else: # FASTQ file
        j = 4   

    ReadDic = dict()
    i = 0
    for row in reads:
        if (i%j == 0): ## Reading only the headers without the "@"
            readName = row.split(' ')[0][1:]
            ReadDic[readName.replace("\n","")] = list([0, '0', 0]) # Create and empty tuple, so A1 and A2 are zero.
        i+=1         
    reads.close()
    return ReadDic

def check_read_name(readName, ending):
	readName=readName.replace("\n","").replace("/1","").replace("/2","")
	return readName+str(ending)

def load_pairedend_reads_data(): 
    print 'Loading sample data...'

    ## Check if input data is in FASTA or FASTQ format
    if(args.R1[-1:]=="a"): # FASTA file
        j = 2
    else: # FASTQ file
        j = 4   
    i = 0
    ReadDic = dict()
    reads = open(args.R1, "r") 
    for row in reads:
        if (i%j == 0): 
            readName = row.split(' ')[0][1:]
            ReadDic[check_read_name(readName,"/1")] = list([0, '0', 0]) # Create and empty tuple, so A1 and A2 are zero.         
        i += 1
    reads.close()
    i = 0
    reads = open(args.R2, "r") 
    for row in reads:
        if (i%j == 0): 
            readName = row.split(' ')[0][1:]
            ReadDic[check_read_name(readName,"/2")] = list([0, '0', 0]) # Create and empty tuple, so A1 and A2 are zero.         
        i += 1
    reads.close()
    return ReadDic

def rename_paired_end(readName, flag):
    check_R2 = int(flag) & 128
    if check_R2 == 0: # R1
        if (args.sample is not None):  # Single-end read
            return(readName)
        else:
            return(readName+str("/1")) # Paired-end read
    elif check_R2 == 128: # R2
        return(readName+str("/2"))
    else:
        return(readName)

def check_distance(ReadDic): #
    for rName in ReadDic:
        if (ReadDic[rName][0] == 0):
            ReadDic[rName][1] = "Not-Map"
        elif ((ReadDic[rName][0] != 0) and (float(ReadDic[rName][0]) < float(args.gamma))):
            ReadDic[rName][1] = "HMR-below-gamma"
        elif ((float(ReadDic[rName][0]) >= float(args.gamma)) and (float(ReadDic[rName][2]) >= float(args.delta))):
            ReadDic[rName][1] = "SHMR-above-delta"
        ##print(rName,";",ReadDic[rName])
    return ReadDic

## Saving read's mapping information for a given reference into the read's dictionary
def save_map_info(mapRow, mapRef, ReadDic):
    mapRow = mapRow.strip('\n')
    row = mapRow.split("\t")
    rName = str(rename_paired_end(row[0],row[1])) ## Name

    try: # The read has already information
        fHMP, refHMP, sHMP = ReadDic[rName]
        if (refHMP != "SHMR-above-delta"): # Non-informative reads
            #flag = row[1] ## flag
            cigar = row[5]
            readLen = len(row[9])
            NM = findNM(row)
            newPer = percentage_alignment(cigar, NM, readLen)
            if (newPer >= fHMP): # New percentage is higher or equal to the fHMP
                ReadDic[rName][2] = ReadDic[rName][0]
                ReadDic[rName][0] = newPer
                ReadDic[rName][1] = mapRef
                # ReadDic = check_delta(ReadDic, rName)
            elif (newPer > sHMP): # New percentage is higher than the sHMP
                ReadDic[rName][2] = newPer
                # ReadDic = check_delta(ReadDic, rName)
    except KeyError: # The read doesn't exists
        print mapRow
        print rName
        print("ERROR: READ "+str(rName)+" DOESN'T EXISTS")
    #print ReadDic
    return ReadDic

## Reading SAM file and loading information. Notice that SAM header must be removed. 
def load_sam_data(samDocs, ReadDic):
    for samfile in samDocs:
        refName = samfile.split(".")[0] ## The ref code is at the file name
        sam = open(samfile, "r")
        print "Loading mapping information: ", samfile

        for row in sam:    
            if row[0][0] != '@':
                #print(row, refName, ReadDic)
                ReadDic = save_map_info(row, refName, ReadDic)
        sam.close()
    return ReadDic

## Saving a list with the names of all the references using as name the names of the input files (i.e: for "0-PM" we will keep "PM")
def load_ref_data(samDocs):
    listofRefs = list()
    listofRefs.append(("Not-Map"))
    listofRefs.append(("HMR-below-gamma"))
    listofRefs.append(("SHMR-above-delta"))
    for samfile in samDocs:
        refName = samfile.split("-")[0] ## The ref code is at the file name ## change this line if needed for different kinds of names
        listofRefs.append((refName))
    return listofRefs

## Creating a dictionary for references, so saving the number of assigned reads per references
def summary_dic(ReadDic):
    SumDic = dict()
    for rName in ReadDic:
        if ReadDic[rName][1] in SumDic: # Reference already exists
            SumDic[ReadDic[rName][1]] += 1    
        else: # Reference doesn't exist yet
            newkey = ReadDic[rName][1]
            SumDic[newkey] = 1
    # Returning the dictionary sorted by the decreasing number of reads
    return SumDic

def popValue(valueKey, inDict):
    try:
        pop_ele = inDict.pop(valueKey)
        return (pop_ele, inDict) # Return the value associated to the key and the Dictionary without the element
    except KeyError:
        return (0, inDict)

## Saving reference's dictionary on a CSV file
def save_map_info_csv(OldDict):
    totalReads = sum(OldDict.values())
    print ("> "+str(totalReads)+" initial reads")

    # Pop top three elements from column
    # print OldDict.items()
    pop_NM, OldDict = popValue('Not-Map', OldDict)
    pop_HMR, OldDict = popValue('HMR-below-gamma', OldDict)
    pop_SHMR, OldDict = popValue('SHMR-above-delta', OldDict)

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
        #csvoutput = open(os.path.dirname(args.output)+"AUXoutput.csv", 'w')
        writer = csv.writer(csvoutput, delimiter=";")

        # Write header
        outRow = list(lines[0].rstrip("\r\n").split(";"))
        try:
            outRow.append((os.path.basename(args.sample)))
        except AttributeError:
            outRow.append((os.path.basename(args.R1)))
        writer.writerow(outRow)
        
        # Write Not-mapped, HMR-below-gamma and SHMR-above-delta
        outRow = list(lines[1].rstrip("\r\n").split(";"))
        outRow.append(("Not-Map"+" ("+str(pop_NM)+")"))
        writer.writerow(outRow)

        outRow = list(lines[2].rstrip("\r\n").split(";"))
        outRow.append(('HMR-below-gamma'+" ("+str(pop_HMR)+")"))
        writer.writerow(outRow)

        outRow = list(lines[3].rstrip("\r\n").split(";"))
        outRow.append(('SHMR-above-delta'+" ("+str(pop_SHMR)+")"))
        writer.writerow(outRow)

        # Write assignments
        for i in range(max(len(OldDict),len(lines)-4)): # Loop through lines in max(file/dictionary)
            # ADD OLD DATA FROM FILE
            try:
                outRow = list(lines[i+4].rstrip("\r\n").split(";")) # +4 for avoiding the first 4 lines
            except IndexError: # If the new data has more lines than the old data
                outRow = list([0]*len(lines[0].rstrip("\r\n").split(";")))
            
            # ADD NEW DATA ON FILE TO A NEW COLUMN
            try:
                outRow.append(sorted_Tupla[i][0]+" ("+str(sorted_Tupla[i][1])+"|"+str(round(Decimal(float(sorted_Tupla[i][1])/float(totalAssig)),5))+")")
            except IndexError: # New data has less lines than old data
                outRow.append(0)
            writer.writerow(outRow)
        os.rename(os.path.splitext(args.output)[0]+"AUXoutput.csv", args.output)
        #os.rename(os.path.splitext(args.output+"AUXoutput.csv")[0], args.output)
    
    ## If output file doesn't exist yet, create it and save the mapping information
    else:
        csvoutput = open(args.output, 'w')
        writer = csv.writer(csvoutput)
        try:
            writer.writerow(list([str(os.path.basename(args.sample))])) # Write sample name
        except AttributeError:
            writer.writerow(list([str(os.path.basename(args.R1))])) # Write sample name
        writer.writerow([str("Not-Map"+" ("+str(pop_NM)+")")]) # Write Not-Map
        writer.writerow([str('HMR-below-gamma'+" ("+str(pop_HMR)+")")]) # Write 'HMR-below-gamma'
        writer.writerow([str('SHMR-above-delta'+" ("+str(pop_SHMR)+")")]) # Write 'SHMR-above-delta'

        for item in sorted_Tupla:
            writer.writerow([str(item[0]+" ("+str(item[1])+"|"+str(round(Decimal(float(item[1])/float(totalAssig)),5))+")")])
    csvoutput.close()
    return 

def inputSamples():
    if (args.sample is not None):
        if (args.sample[-6:] != '.fasta') and (args.sample[-6:] != '.fastq'):
            print 'Input sample data is in a wrong format'
            sampleFile = raw_input("Which is your target file (FASTA/FASTQ)?:")
        else:
            sampleFile = args.sample
        return(load_reads_data(sampleFile), "R1")
    elif (args.R1 is not None) and (args.R2 is not None):
        return(load_pairedend_reads_data(), "PE")
    else:
        print "ERROR: specified samples are wrong"
    return()

## MAIN PROGRAMME

def main():

    ## Saving read's name in a dictionary
    ReadDic, end = inputSamples()
    ## Saving mapping information
    # Gettin all the SAM files from folder
    if (args.SAMfolder  is not None):
        os.chdir(args.SAMfolder)
        samfiles = [f for f in os.listdir(args.SAMfolder) if f.endswith(".sam")]
    else:
        samfol = raw_input("Where are the SAM files?:")
        os.chdir(samfol)
        samfiles = [f for f in os.listdir(samfol) if f.endswith(".sam")]
    print(samfiles)

    ReadDic = load_sam_data(samfiles, ReadDic)
    #RefList = load_ref_data(samfiles)
    ReadDic = check_distance(ReadDic)
    
    print("Summarizing data...")
    summary = summary_dic(ReadDic)
    save_map_info_csv(summary)
    return

## MAIN 
main()
