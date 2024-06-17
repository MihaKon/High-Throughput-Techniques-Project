#!/usr/bin/env bash

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

echo "CALLSET REFINEMENT"
echo "-----------------"

gatk SelectVariants -V NA12878_variants.vcf.gz -select-type SNP -O NA12878.snp.vcf.gz
gatk VariantRecalibrator -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta \
    -V NA12878.snp.vcf.gz \
    -O NA12878.snp.recal.gz \
    --tranches-file NA12878.snp.tranches \
    --mode SNP \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /db/hapmap_3.3.hg38.vcf.gz \
    -resource:omni,known=false,training=true,truth=false,prior=12.0 /db/1000G_omni2.5.hg38.vcf.gz \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 /db/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /db/dbsnp_138.hg38.vcf.gz \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    --rscript-file NA12878.snp.plots.R
gatk ApplyVQSR -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta \
    -V NA12878.snp.vcf.gz \
    --tranches-file NA12878.snp.tranches \
    --recal-file NA12878.snp.recal.gz \
    --truth-sensitivity-filter-level 99.5 \
    -mode SNP \
    -O NA12878.snp.filtered.vcf.gz

gatk SelectVariants -V NA12878_variants.vcf.gz -select-type INDEL -O NA12878.indels.vcf.gz
gatk VariantRecalibrator -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta \
    -V NA12878.indel.vcf.gz \
    -O NA12878.indel.recal.gz \
    --tranches-file NA12878.indel.tranches \
    --mode INDEL \
    -resource:mills,known=false,training=true,truth=true,prior=10.0 /db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /db/dbsnp_138.hg38.vcf.gz \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    --rscript-file NA12878.indel.plots.R
gatk ApplyVQSR -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta -V NA12878.indel.vcf.gz -O NA12878.indel.vqsr.vcf.gz --tranches-file NA12878.indel.tranches --recal-file NA12878.indel.recal.gz --truth-sensitivity-filter-level 99.5 -mode INDEL

gatk MergeVcfs -I NA12878.snp.filtered.vcf.gz -I NA12878.indel.vqsr.vcf.gz -O NA12878.filtered.vcf.gz
gatk CNNScoreVariants -V NA12878.filtered.vcf.gz -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta -O NA12878.CNN1D.vcf
gatk CNNScoreVariants -V NA12878.filtered.vcf.gz -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta -O NA12878.CNN2D.vcf -I NA12878_bqsr.bam -tensor-type read_tensor
gatk FilterVariantTranches --output NA12878.CNN1D.filtered.vcf --snp-tranche 99.95 --indel-tranche 99.4 --resource /db/hapmap_3.3.hg38.vcf.gz --resource /db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --info-key CNN_1D --variant NA12878.CNN1D.vcf
gatk FilterVariantTranches --output NA12878.CNN2D.filtered.vcf --snp-tranche 99.95 --indel-tranche 99.4 --resource /db/hapmap_3.3.hg38.vcf.gz --resource /db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --info-key CNN_2D --variant NA12878.CNN2D.vcf 

gatk SelectVariants -V NA12878.CNN1D.filtered.vcf -O NA12878.CNN1D.filtered_out.vcf --exclude-filtered true
gatk SelectVariants -V NA12878.CNN2D.filtered.vcf -O NA12878.CNN2D.filtered_out.vcf --exclude-filtered true

