## 1. Introduction

#### 1.1. What's the *gamma-delta* algorithm?

The *gamma-delta* algorithm aims to identify reads that provide taxonomical information at species level. For each mapping of a read *r* against a reference we obtain a mapping ratio *A*, which is calculated by dividing the number of matching nucleotides from the query read *r* to the target sequences (*nm*) among the total number of nucleotides involved in the alignment (*nt*). Then, a read *r* will be assigned to species *i* when the mapping ratio *A* against species *i* is higher than *gamma* and the alternative species’ mapping ratios are bellow *delta* (Garrido-Sanz et al. 2019, *MBMG*).

#### 1.2. Setup
For the proper performance of the pipeline the user should modify the below listed paths in the "submit_pipeline.sh" script to its own:

  TRIMMOMATIC_PATH=/path/to/Trimmomatic<br>
  BWA_PATH=/path/to/BWA<br>
  ST_PATH=/path/to/SAMtools<br>
  gd_PATH=/path/to/g-d_algorithm_script<br>
  REF_PATH=/path/to/references<br>
  R1=/path/to/sample<br>

## 2. Command-line and options

cd /path/to/script<br>
./submit_pipeline.sh

#### Singled-end reads assingment
python g-d_algorithm_sort_and_sv.py -g 0.99 -d 0.98 -m /path/to/sams/folder -s reads.se.fastq -o assignment_SE-reads.csv;<br>

#### Paired-end reads assingment
python g-d_algorithm_sort_and_sv.py -g 0.99 -d 0.98 -m /path/to/sams/folder -r1 forward.pe.fastq -r2 reverse.pe.fastq -o assignment_PE-reads.csv;<br>

## 3. Output format

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

## 4. Authors
* Lidia Garrido-Sanz [lidia.garrido@uab.cat] 
* Miquel Àngel Senar
* Josep Piñol

## 5. Reporting bugs
All reports and feedbacks are highly appreciate. Please report any suggestion on github or by email to lidia.garrido@uab.cat. 

## 6. Disclaimer
The authors provided the information in good faith. Under no circumstance shall authors and the Autonomous University of Barcelona have any liability to you for any loss or damage of any kind incurred as a result of the use of the information provided. The use of this information is solely at your own risk. 
