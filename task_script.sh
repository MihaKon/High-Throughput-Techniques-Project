#!/bin/bash

ProcCount=$(nproc)

echo "PRE-PROCESSING"
echo "-----------------"

bwa mem -t $ProcCount -o "NA12878.bam" ./ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta NA12878_chr20_1.fastq.gz NA12878_chr20_2.fastq.gz
samtools stats -@ $ProcCount "NA12878.bam" > "NA12878_stats.txt" 
samtools sort -@ $ProcCount -o "NA12878_sorted.bam" "NA12878.bam" 
gatk AddOrReplaceReadGroups -I NA12878_sorted.bam -O NA12878_tagged.bam -LB 'cos' -PL 'ILLUMINA' -PU 'cos' -SM 'NA12878' -ID "NA12878.1"
gatk MarkDuplicates -I "NA12878_tagged.bam" -O "NA12878_marked_dupl.bam" -M "md_metrics.txt"
gatk BaseRecalibrator -I NA12878_marked_dupl.bam --known-sites /db/dbsnp_138.hg38.vcf.gz -O NA12878_bl.rt -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta
gatk ApplyBQSR -I "NA12878_tagged.bam" -bqsr "NA12878_bl.rt"  -O "NA12878_bqsr.bam"
gatk BaseRecalibrator -I "NA12878_bqsr.bam" -O "NA12878_bqsr_br.rt" --known-sites /db/dbsnp_138.hg38.vcf.gz -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta
gatk AnalyzeCoveriants -before "NA12878_bl.rt" -after "NA12878_bqsr_bl.rt" -plots "analyze_base_bqsr.pdf"



echo "VARIANT DISCOVERY"
echo "-----------------"

gatk HaplotypeCaller -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta -I NA12878_bqsr.bam -O NA12878_variants.vcf.gz
gatk CollectVariantCallingMetrics -I NA12878_variants.vcf.gz -O NA12878_variant_metrics --DBSNP /db/dbsnp_138.hg38.vcf.gz

