#!/usr/bin/env bash

cp /data/NA12878_chr20_1.fastq.gz /task_data/
cp /data/NA12878_chr20_2.fastq.gz /task_data/

ProcCount=$(nproc)

echo "PRE-PROCESSING"
echo "-----------------"

bwa mem -t $ProcCount -o "NA12878.bam" /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta NA12878_chr20_1.fastq.gz NA12878_chr20_2.fastq.gz
samtools stats -@ $ProcCount "NA12878.bam" > "NA12878_stats.txt" 
samtools sort -@ $ProcCount -o "NA12878_sorted.bam" "NA12878.bam" 
gatk AddOrReplaceReadGroups -I NA12878_sorted.bam -O NA12878_tagged.bam -LB 'cos' -PL 'ILLUMINA' -PU 'cos' -SM 'NA12878' -ID "NA12878.1"
gatk MarkDuplicates -I "NA12878_tagged.bam" -O "NA12878_marked_dupl.bam" -M "md_metrics.txt"
gatk BaseRecalibrator -I NA12878_marked_dupl.bam \
    --known-sites /db/dbsnp_138.hg38.vcf.gz \
    -O NA12878_br_hcr.rt \
    -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta
gatk ApplyBQSR -I "NA12878_tagged.bam" -bqsr "NA12878_br_hcr.rt"  -O "NA12878_bqsr_hcr.bam"
gatk BaseRecalibrator -I "NA12878_bqsr_hcr.bam" \
    -O "NA12878_bqsr_br.rt" \
    --known-sites /db/dbsnp_138.hg38.vcf.gz \
    -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta
gatk AnalyzeCovariates -before "NA12878_br_hcr.rt" -after "NA12878_bqsr_br.rt" -plots "analyze_base_bqsr.pdf"

echo "VARIANT DISCOVERY"
echo "-----------------"

gatk HaplotypeCaller -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta.fai -I NA12878_br_hcr.rt -O NA12878_br.vcf.gz
exit 0
gatk HaplotypeCaller -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta.fai -I NA12878_bqsr_hcr.bam -O NA12878_bqsr.vcf.gz

echo "CALLSET REFINEMENT"
echo "-----------------"

function get_analyze_ready() {
    
    local filename=$1
    local variant_type=$2 

    gatk SelectVariants -V "${filename}.vcf.gz" -select-type ${variant_type} -O "${filename}.${variant_type}.vcf.gz"
    gatk VariantRecalibrator -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta.fai \
        -V "${filename}.${variant_type}.vcf.gz" \
        -O "${filename}.${variant_type}.recal.gz" \
        --tranches-file "${filename}.${variant_type}.tranches" \
        --mode ${variant_type} \
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /db/hapmap_3.3.hg38.vcf.gz \
        -resource:omni,known=false,training=true,truth=false,prior=12.0 /db/1000G_omni2.5.hg38.vcf.gz \
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 /db/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /db/dbsnp_138.hg38.vcf.gz \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        --rscript-file "${filename}.${variant_type}.plots.R"
    gatk ApplyVQSR -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta.fai \
        -V "${filename}.${variant_type}.vcf.gz" \
        --tranches-file "${filename}.${variant_type}.tranches" \
        --recal-file "${filename}.${variant_type}.recal.gz" \
        --truth-sensitivity-filter-level 99.5 \
        -mode ${variant_type} \
        -O "${filename}.${variant_type}.filtered.vcf.gz"
}

function add_score() {

    local filename=$1

    gatk MergeVcfs -I "${filename}.snp.filtered.vcf.gz" -I "${filename}.indel.vqsr.vcf.gz" -O "${filename}.filtered.vcf.gz"
    gatk CNNScoreVariants -V "${filename}.filtered.vcf.gz" \
        -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta.fai \
        -O "${filename}.CNN1D.vcf"
    gatk CNNScoreVariants -V "${filename}.filtered.vcf.gz" \
        -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta.fai \
        -O "${filename}.CNN2D.vcf" \
        -I "${filename}_bqsr.bam" \
        -tensor-type read_tensor
    gatk FilterVariantTranches --output "${filename}.CNN1D.filtered.vcf" \
        --snp-tranche 99.95 \
        --indel-tranche 99.4 \
        --resource /db/hapmap_3.3.hg38.vcf.gz \
        --resource /db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        --info-key CNN_1D \
        --variant "${filename}.CNN1D.vcf"
    gatk FilterVariantTranches --output "${filename}.CNN2D.filtered.vcf" \
        --snp-tranche 99.95 \
        --indel-tranche 99.4 \
        --resource /db/hapmap_3.3.hg38.vcf.gz \
        --resource /db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        --info-key CNN_2D \
        --variant "${filename}.CNN2D.vcf"

    gatk SelectVariants -V "${filename}.CNN1D.filtered.vcf" -O "${filename}.CNN1D.filtered_out.vcf" --exclude-filtered true
    gatk SelectVariants -V "${filename}.CNN2D.filtered.vcf" -O "${filename}.CNN2D.filtered_out.vcf" --exclude-filtered true
}

select_variants() {
    local filename=$1

    gatk SelectVariants -V "${filename}.CNN1D.filtered.vcf" -O "${filename}.CNN1D.filtered_out.vcf" --exclude-filtered true
    gatk SelectVariants -V "${filename}.CNN2D.filtered.vcf" -O "${filename}.CNN2D.filtered_out.vcf" --exclude-filtered true
}

echo "BQSR ANALYSIS"
get_analyze_ready "NA12878_bqsr" "SNP"
get_analyze_ready "NA12878_bqsr" "INDEL"

echo "BQSR SCORE"
add_score "NA12878_bqsr"

echo "BQSR SELECT VARIANTS"
select_variants "NA12878_bqsr"
