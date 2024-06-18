#!/usr/bin/env bash

if [ ! -f /task_data/NA12878_chr20_1.fastq.gz ]; then
    cp /data/NA12878_chr20_1.fastq.gz /task_data/
fi

if [ ! -f /task_data/NA12878_chr20_2.fastq.gz ]; then
    cp /data/NA12878_chr20_2.fastq.gz /task_data/
fi

ProcCount=$(nproc)

echo "PRE-PROCESSING"
echo "-----------------"

if [ ! -f "NA12878_bqsr_hcr.bam" ]; then
    bwa mem -t $ProcCount -o "NA12878.bam" /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta NA12878_chr20_1.fastq.gz NA12878_chr20_2.fastq.gz
    samtools stats -@ $ProcCount "NA12878.bam" > "NA12878_stats.txt" 
    samtools sort -@ $ProcCount -o "NA12878_sorted.bam" "NA12878.bam" 
    gatk AddOrReplaceReadGroups -I NA12878_sorted.bam -O NA12878_tagged.bam -LB 'cos' -PL 'ILLUMINA' -PU 'cos' -SM 'NA12878' -ID "NA12878.1"
    gatk MarkDuplicates -I "NA12878_tagged.bam" -O "NA12878_marked_dupl.bam" -M "md_metrics.txt" -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta
    gatk BaseRecalibrator -I NA12878_marked_dupl.bam \
        --known-sites /db/dbsnp_138.hg38.vcf.gz \
        -O NA12878_br_hcr.rt \
        -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta
    gatk ApplyBQSR -I "NA12878_marked_dupl.bam" -bqsr "NA12878_br_hcr.rt"  -O "NA12878_bqsr_hcr.bam"
    gatk BaseRecalibrator -I "NA12878_bqsr_hcr.bam" \
        -O "NA12878_bqsr_br.rt" \
        --known-sites /db/dbsnp_138.hg38.vcf.gz \
        -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta
    gatk AnalyzeCovariates -before "NA12878_br_hcr.rt" -after "NA12878_bqsr_br.rt" -plots "analyze_base_bqsr.pdf"
    gatk AnalyzeCovariates -before "NA12878_br_hcr.rt" -after "NA12878_bqsr_br.rt" -csv "analyze_base_bqsr.csv"
fi

echo "VARIANT DISCOVERY"
echo "-----------------"

if [ ! -f "NA12878_br.vcf" ]; then
    gatk HaplotypeCaller -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta -I NA12878_marked_dupl.bam -O NA12878_br.vcf
fi
if [ ! -f "NA12878_bqsr.vcf" ]; then
    gatk HaplotypeCaller -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta -I NA12878_bqsr_hcr.bam -O NA12878_bqsr.vcf
fi

echo "CALLSET REFINEMENT"
echo "-----------------"

function get_analyze_ready() {
    
    local filename=$1
    local variant_type=$2
    local truth_sensitivity_level=${3:-"99.5"}    


    gatk SelectVariants -V "${filename}.vcf" -select-type ${variant_type} -O "${filename}.${variant_type,,}.vcf"
    if [ "${variant_type,,}" == "indel" ]; then
        gatk VariantRecalibrator -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta \
            -V "${filename}.indel.vcf" \
            -O "${filename}.indel.recal" \
            --tranches-file "${filename}.indel.tranches" \
            --mode INDEL \
            -resource:mills,known=false,training=true,truth=true,prior=10.0 /db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /db/dbsnp_138.hg38.vcf.gz \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            --rscript-file "${filename}.indel.plots.R"
    else
        gatk VariantRecalibrator -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta \
            -V "${filename}.${variant_type,,}.vcf" \
            -O "${filename}.${variant_type,,}.recal" \
            --tranches-file "${filename}.${variant_type,,}.tranches" \
            --mode SNP \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /db/hapmap_3.3.hg38.vcf.gz \
            -resource:omni,known=false,training=true,truth=false,prior=12.0 /db/1000G_omni2.5.hg38.vcf.gz \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 /db/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /db/dbsnp_138.hg38.vcf.gz \
            -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
            --rscript-file "${filename}.${variant_type,,}.plots.R"
    fi
    gatk ApplyVQSR -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta \
        -V "${filename}.${variant_type,,}.vcf" \
        --tranches-file "${filename}.${variant_type,,}.tranches" \
        --recal-file "${filename}.${variant_type,,}.recal" \
        --truth-sensitivity-filter-level $truth_sensitivity_level \
        -mode ${variant_type} \
        -O "${filename}.${variant_type,,}.cnnready.vcf"
}

function add_score() {

    local filename=$1
    local bam_filename=$2
    local snp_tranche=$3
    local indel_tranche=$4

    gatk MergeVcfs -I "${filename}.snp.cnnready.vcf" -I "${filename}.indel.cnnready.vcf" -O "${filename}.filtered.vcf"
    gatk CNNScoreVariants -V "${filename}.filtered.vcf" \
        -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta \
        -O "${filename}.CNN1D.vcf"
    gatk CNNScoreVariants -V "${filename}.filtered.vcf" \
        -R /ref/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fasta \
        -O "${filename}.CNN2D.vcf" \
        -I "${bam_filename}.bam" \
        -tensor-type read_tensor
    gatk FilterVariantTranches --output "${filename}.CNN1D.filtered.vcf" \
        --snp-tranche $snp_tranche \
        --indel-tranche $indel_tranche \
        --resource /db/hapmap_3.3.hg38.vcf.gz \
        --resource /db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        --info-key CNN_1D \
        --variant "${filename}.CNN1D.vcf"
    gatk FilterVariantTranches --output "${filename}.CNN2D.filtered.vcf" \
        --snp-tranche $snp_tranche \
        --indel-tranche $indel_tranche \
        --resource /db/hapmap_3.3.hg38.vcf.gz \
        --resource /db/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
        --info-key CNN_2D \
        --variant "${filename}.CNN2D.vcf"
}

select_variants() {
    local filename=$1

    gatk SelectVariants -V "${filename}.CNN1D.filtered.vcf" -O "${filename}.CNN1D.filtered_out.vcf" --exclude-filtered true
    gatk SelectVariants -V "${filename}.CNN2D.filtered.vcf" -O "${filename}.CNN2D.filtered_out.vcf" --exclude-filtered true
}

if [ ! -f "NA12878_bqsr.CNN2D.filtered_out.vcf" ]; then
    echo "BQSR ANALYSIS"
    get_analyze_ready "NA12878_bqsr" "SNP"
    get_analyze_ready "NA12878_bqsr" "INDEL"

    echo "BQSR SCORE"
    add_score "NA12878_bqsr" "NA12878_bqsr_hcr" 99.95 99.5

    echo "BQSR SELECT VARIANTS"
    select_variants "NA12878_bqsr"
fi

echo "BR ANALYSIS"
get_analyze_ready "NA12878_br" "SNP" 0 
get_analyze_ready "NA12878_br" "INDEL" 0 

echo "BR SCORE"
add_score "NA12878_br" "NA12878_marked_dupl" 98.00 92.0

echo "BR SELECT VARIANTS"
select_variants "NA12878_br"
