# uditaspreprocess
## Illumina demultiplexing for uditas on miniseq:
1.	Input barcodes for i7 and i5 on a CSV file named SampleSheetUditas.csv
2.	Make sure i5 indexes are reverse complement (e.g. i501 == AGGCTATA). i7 is sequenced on the sense strand, so leave it be. 
3.	“Sample_ID” and “Sample_Name” columns should both be named "i7##_i5##" (this is important later)
4.	Make sure cell A10 (below [Reads]) is equal to R1 + R2. If you sequence 145+145, this number is 290. In your Illumina run folder, this information is available in “RunInfo.xml”
5.	Eliminate spaces (ctrl+g , replace space with nothing)
6.	Use dos2unix to convert line breaks to something readable by bcl2fastq
7.	Rename sample sheet to anything but “SampleSheet.csv”. This is necessary to avoid issues later.
8.	Upload your sample sheet to the cluster and put it into the Illumina run folder↓. You can transfer your runfolder directly from basespace (fast) using the BaseSpace CLI portable software on the cluster (bin/bs download runs -i <RUN ID> -o <OUTPUT DIR>). Alternatively, you can upload your reads t
 
9.	Now, use bcl2fastq to demultiplex. This is usually performed on a cluster due to the >32G memory requirements of bcl2fastq. Load this on the GHPCC using: module load bcl2fastq2/2.20.0
a.	The number of index reads indicated in the command MUST match the number of nucleotides used in your sample sheet.
b.	You must “mask” UMI nucleotides on the i5 read and any extra nucleotides on the i7 index 
c.	For fast demultiplexing, use “-p” flag + the number of cores available to you
d.	Must have >32G RAM on your machine or allocated to your cluster node/job
e.	Bcl2fastq allows for 1 index nucleotide mismatch. If you use indexes that have <2 mismatches, you must set “--barcode-mismatches” to 0.
f.	Successful demultiplexing will populate a message, “Processing completed with 0 errors and 0 warnings.”

```
bcl2fastq --use-bases-mask Y*,I6N2,N9I8,Y* --sample-sheet SampleSheetUditas.csv -r 1 -p 15 -w 1 --runfolder-dir . -o fastq
```

10.	Rename your output files to i7##_i5##_R#.fq.gz and gunzip
```
for i in *.fastq.gz; do mv $i ${i%_S*}.fq.gz;done
gunzip *
```
Getting UMIs:
1.	Use bcl2fastq to get first 9 nucleotides, without demultiplexing:
bcl2fastq -r 1 -p 15 -w 1 --use-bases-mask Y1N*,N*,I9N8,Y1N* --create-fastq-for-index-reads -o index
2.	Delete undetermined read 1 and 2 files and rename I2 undetermined file to index.fq.gz, placing it in the fastq folder with your samples and unzipping
cd index && rm Undet*_R* Undet*_I1*
mv Undetermined_S0_L001_I2_001.fastq.gz ../fastq/index.fq.gz
cd ../fastq && gunzip index.fq.gz 
At this point, you’ll need to upload your “fastq” data folder to a computer with docker. We could not get the uditas software on the cluster, unfortunately. I recommend rclone for the transfer. Pairing UMIs and R1/R2 FASTQs:
1.	Use fastq_pair utility available in umasstr/uditas container to find index.fq reads that overlap with each sample’s R1 read.
2.	Rename index.fq.paired.fq to i7##_i5##.umi.fq and delete extra files
```
for i in *R1.fq; do echo $i && time fastq_pair $i index.fq && mv index.fq.paired.fq ${i%R1*}umi.fq && rm *single* *paired*;done
```
Trimming FASTQ headers and making sample folders:
1.	I create a folder called “uditas01” (attempt 1). 
2.	Move your Bowtie2 index files and whole genome FASTA here
3.	Make a folder for each sample (must me i7##_i5##)
4.	Make a folder inside those called “fastq_files” 
5.	Copy each sample’s R1, R2, and umi file into this fastq folder, renaming with fastq extension and gzipping, and trimming off the superfluous info from FASTQ header that comes after the space.
```
mkdir uditas01
for i in i7*R1*.fq; do echo $i && mkdir uditas01/${i%_R*}/ && mkdir uditas01/${i%_R*}/fastq_files/;done
for i in i7*R*.fq; do echo $i && time cut -f1 -d' ' $i | gzip > uditas01/${i%_R*}/fastq_files/${i%.fq}.fastq.gz;done

for i in i7*umi.fq; do echo $i && time cut -f1 -d' ' $i | gzip > uditas01/${i%_umi*}/fastq_files/${i%.fq}.fastq.gz;done
```
6.	Structure within uditas01 folder should look like this for each sample:
 
7.	In each sample fastq_files folder, check to see that the files have the same # lines
8.	Before running, make sure computer can’t go to sleep due to inactivity.
Final product:
 
Reference Genomes:
Make a folder containing your whole-genome FASTA and associated indexes. All of the below files are necessary. You can also download the hg38 and mm10 folders at their respective hyperlinks which contain all files needed. 
 
If you’re editing human loci in a model organism (human gene in a mouse e.g.), you’ll need to make a custom reference. 
1.	Make a new FASTA 
2.	index (samtools) 
3.	bowtie index (bowtie2)
4.	2bit (kentutils)
Coordinates proof:
Check your sample_info.csv coordinates in the umasstr/uditas docker container. Mount instructions below. 
export BOWTIE2_INDEXES=/DATA/Bowtie2Index
export GENOMES_2BIT=/DATA/Bowtie2Index

python
import matplotlib
import matplotlib.pyplot as plt
import pylab
import pandas as pd
import numpy as np
import os
import sys
import gzip
import itertools
import operator
import subprocess
import twobitreader

genome = twobitreader.TwoBitFile("/DATA/Bowtie2Index/mm10.2bit")

def reverse_complement(seq):
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return "".join([seq_dict[base] for base in reversed(seq)])

chr5	123182085	123182107	-	2-22_FWD GGAATTTCCCAATTCTGATCCT	

#Examples
#>>> reverse_complement(genome['chr5'][ 123182085:123182107 ])
#'GGAATTTCCCAATTCTGATCCT'

#chr5	123182002	123182026	+	sg1	#ACCTCCGAGAATGGGACTCACCCT
#>>> genome['chr5'][ 123182002:123182026 ]
#'ACCTCCGAGAATGGGACTCACCCT'
#chr5	123181414	123181438	-	sg2
#CATCAAGCAAGGGAAAGTGAGTAC
#>>> reverse_complement(genome['chr5'][ 123181414:123181438 ])
#'CATCAAGCAAGGGAAAGTGAGTAC'

#i505
#chr5	123181316	123181337	+	4-10_REV
#AGCCATACATCTTGGAACCAG
#>>> genome['chr5'][ 123181316:123181337 ]
#'AGCCATACATCTTGGAACCAG'

Running the pipeline:
1.	start umasstr/uditas interactive docker container, mounting your uditas01 directory to “/DATA” 
docker run -it -v <your fastq directory>:/DATA -w /DATA umasstr/uditas
2.	Define reference file environment variables 

export BOWTIE2_INDEXES=/DATA/Bowtie2Index

export GENOMES_2BIT=/DATA/Bowtie2Index

3.	start uditas, skipping demultiplexing on “DATA” directory

uditas -ncpu 14 -skip_demultiplexing 1 -folder_genome_2bit <FOLDER_GENOME_2BIT without trailing / > ./

Analyzing The Data
When the UDiTaS analysis pipeline is finished, there will be two key excel files containing read numbers for different editing events: “DATA_big_results” and “results_combined_detailed”. The “big_results” file has the most simplified categories and the “results_combined_detailed” file contains more expanded categories.
Important things to note about the results files: 
1.	Column Headers labeled “Collapsed” have had the read numbers collapsed based on the UMIs present in each read. Therefore, these categories are counting each unique repair event.
2.	The difference between collapsed and uncollapsed categories will give you an idea of your library complexity and whether you sequenced them with enough depth. Divide the collapsed number by the uncollapsed number and it the percentage is low, that’s an indicator that you’ve sequenced most of the unique UMIs in your library (ie: you had enough depth to sequence most of the unique molecules that were present in the library)
3.	In the “results_combined_detailed” file, the sub-categories of a repair event are not cumulative since some events can fit into multiple categories. For example: an event can be counted in both the large_deletion_cut1_total_indels_collapsed and the large_deletion_cut1_total_deletions_collapsed categories. The large_deletion_cut1_total_reads_collapsed category, and similar categories, will count all events including indels and perfect repair etc. 
4.	Categories like the following: “1a_1a_cut1_total_reads_collapsed” measure chromosomal transversions. See Figure 2C in the original UDiTaS publication (Giannoukos G, Ciulla DM, Marco E, Abdulkerim HS, Barrera LA, Bothmer A, Dhanapal V, Gloskowski SW, Jayaram H, Maeder ML, Skor MN, Wang T, Myer VE, Wilson CJ (2018) UDiTaSTM, a genome editing detection method for indels and genome rearrangements. BMC genomics, 19(1):212)
5.	From Eugenio Marco of Editas: ““target plus plasmid reads” means that those are reads that had one end on the target region (where the UDiTaS primer binds) and the other one on the plasmid/AAV. If your primer were to misprime from the plasmid and your reads have no sequence from the target region, they would be counted in the plasmid_only_total_reads” 
6.	In the “Big_results” file the total_aligned_junctions_collapsed (column BU) are the total number of events counted. To find a percentage of a certain type of repair event, take the number from the repair event category (ex: large_deletion_cut1_total_reads_collapsed) and divide by the number in total_aligned_junctions_collapsed.

There are many ways to represent the data like a pie chart detailing the percent of reads for different repair outcomes for a single sample or group (See original paper FIG2A for an example Giannoukos G, Ciulla DM, Marco E, Abdulkerim HS, Barrera LA, Bothmer A, Dhanapal V, Gloskowski SW, Jayaram H, Maeder ML, Skor MN, Wang T, Myer VE, Wilson CJ (2018) UDiTaSTM, a genome editing detection method for indels and genome rearrangements. BMC genomics, 19(1):212). Stacking bar graphs to show the percent of reads for repair outcomes for different groups side by side (See LCA10 study publication FIG1E Maeder ML, Stefanidakis M, Wilson CJ, Baral R, Barrera LA, Bounoutas GS, Bumcrot D, Chao H, Ciulla DM, DaSilva JA, Dass A, Dhanapal V, Fennell TJ, Friedland AE, Giannoukos G, Gloskowski SW, Glucksmann A, Gotta GM, Jayaram H, Haskett SJ, Hopkins B, Horng JE, Joshi S, Marco E, Mepani R, Reyon D, Ta T, Tabbaa DG, Samuelsson SJ, Shen S, Skor MN, Stetkiewicz P, Wang T, Yudkoff C, Myer VE, Albright CF, Jiang H (2019) Development of a gene-editing approach to restore vision loss in Leber congenital amaurosis type 10. Nature medicine, 25(2):229–233.)
I think it is important to include the total read numbers for each category you’re plotting just to show that you have enough depth. 
Calculating the ratio of collapsed:uncollapsed reads is also important to show the quality. Also know the percent aligned reads (percent_aligned_all_amplicons category, Column BR in “big_results”). A low percentage of aligned reads can indicate a quality problem with the experiment (in one case, I had a primer for UDiTaS that performed very poorly and I was able to identify this because libraries prepped with this primer had a very low percentage of aligned reads.) 
In “big_results”, columns BW, BX, BY, BZ, and CA are the percentages of the listed repair event calculated by dividing the repair event reads by the total reads in column BU (these are the most simplified results, you can use “results_combined_detailed” to further break down these categories if desired.

