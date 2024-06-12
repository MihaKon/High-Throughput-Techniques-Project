#!/bin/bash
set -x

echo "PRE-PROCESSING"
echo "-----------------"
echo "Map to Reference; Generate Mapping Stats | Start"

ProcCount = $(nproc)

bwa mem -t $ProcCount -o "NA12878.sam" ./ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta NA12878_chr20_1.fastq.gz NA12878_chr20_2.fastq.gz
samtools stats -@ $ProcCount "NA12878.sam" > "NA12878_stats.txt" 

echo "Map to Reference; Generate Mappings Stats | End"
echo "Mark Duplicates - Start"
echo "Mark Duplicates - End"
echo "VARIANT DISCOVERY"
echo "-----------------"
