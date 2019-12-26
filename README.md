## 1. Introduction

#### What's the *gamma-delta* workflow?

* is an automated pipeline for analyses of DNA samples that provides a quantitative estimate of the species that are part of such samples. 
* Given a DNA sample and a set of reference genomes, corresponding to the possible species included in the sample, the algorithm generates an output file that includes all the species identified in the sample and their relative abundance.
* Input sample might consist on a set of single- or paired-end reads in FASTA or FASTAQ format. Workflow output is a simple text file in csv format.
* The core of the workflow is the *gamma-delta* algorithm that classifies reads of the sample by using a series of thresholds (gamma and delta parameters) to ensure the accuracy of the quantitative estimate. Details about this process are described at the paper mentioned at the Citation section.

The *gamma-delta* algorithm aims to identify reads that provide taxonomical information at species level. In particular, it is designed to retain only those reads that help for identifying species. Briefly, what gamma-delta algorithm does is for each mapping of a read *r* against a reference, it obtains a mapping ratio *A*, which is calculated by dividing the number of matching nucleotides from the query read *r* to the target sequences among the total number of nucleotides involved in the alignment. Then, a read *r* will be assigned to species *i* when the mapping ratio *A* against species *i* is higher than gamma and the alternative species’ mapping ratios are bellow delta (Garrido-Sanz et al. 2019, MBMG). This algorithm has been written in python2.7 language and runs on command line under Linux.

## 2. Setup

#### 2.1. Tools
The following tools need to be installed in the system to run the pipeline: Trimmomatic, BWA aligner and SAMtools. 
The *gamma-delta* algorithm requires Python 2.7 as well as the csv, argparse, os, operator and decimal libraries.

#### 2.2. Paths
For the correct execution of the pipeline, different paths have to be at the “gamma-delta_workflow.sh” script:

  TRIMMOMATIC_PATH=/path/to/Trimmomatic<br>
  BWA_PATH=/path/to/BWA<br>
  ST_PATH=/path/to/SAMtools<br>
  gd_PATH=/path/to/gamma-delta_algorithm_script<br>
  REF_PATH=/path/to/references<br>
  R1=/path/to/sample<br>
 
 #### 2.3. BWA indexes
The *gamma-delta* workflow uses BWA as a read mapper. This requires the existence of the indices of each of the references against which the sample is compared to. Index generation only needs to be performed once and, therefore, the workflow script can be modified to avoid index recomputation when new samples are analyzed and reference indexes already exit. 

## 3. Command-line and options

cd /path/to/script<br>

#### For single-end reads: <br>
  ./gamma-delta-workflow.sh reads.fastq <br>
#### For paired-end reads:<br>
  ./gamma-delta-workflow.sh forward_reads_R1.fastq reverse_reads_R2.fastq <br>

## 4. Output format

**Column header:** Query name of the sample<br>
**Not-Map:** Number of reads that did not map to any reference<br>
**HMR-below-gamma:** Number of reads that were removed by *gamma* threshold<br>
**SHMR-above-delta:** Number of reads that were removed by *delta* threshold<br>
**List of recovered species:** Name of the reference as the name of the SAM file (Number of reads \| Relative proportion of reads)<br>

**Example:**<br>

| Sample 1  | Sample 2 |
| --- | ---|
| Not-Map (50)  | Not-Map (60)  |
| HMR-below-gamma (30)  | HMR-below-gamma (38)   |
| SHMR-above-delta (20) | SHMR-above-delta (2)  |
| Reference 1 (90 \| 0.9)  | Reference 4 (55 \| 0.55)  |
| Reference 2 (5 \| 0.05)  | Reference 2 (30 \| 0.30)  |
| Reference 3 (4 \| 0.04)  | Reference 1 (10 \| 0.10)  |
| Reference 4 (1 \| 0.01)  | 0  |

## 5. Authors
* Lidia Garrido-Sanz [lidia.garrido@uab.cat] 
* Miquel Àngel Senar
* Josep Piñol

## 6. Reporting bugs
All reports and feedbacks are highly appreciate. Please report any suggestion on github or by email to lidia.garrido@uab.cat. 

## 7. Disclaimer
The authors provided the information and software in good faith. Under no circumstance shall authors and the Universitat Autònoma de Barcelona have any liability for any loss or damage of any kind incurred as a result of the use of the information and software provided. The use of this tool is solely at your own risk.

## 8. Citation

Comming soon...
