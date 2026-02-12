#!/bin/bash

# Configuration - Placeholders
# Replace these with actual paths to your data
TUMOUR_BAM="/nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/ALPACA-WGS/data/bams/test.bqsr.bam"
NORMAL_BAM="/nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/ALPACA-WGS/data/bams/normal.bqsr.bam"
TUMOUR_NAME="LTX1231_SU_T1-R1--8c625f964d"
NORMAL_NAME="LTX1231_BS_GL--59c4c25703"

# Reference files (1000 Genomes)
# Battenberg requires 1000 Genomes loci and allele files split by chromosome
# Provide the prefix for these files
BEAGLE_BASEDIR="/nemo/project/proj-tracerx-wgs/working/picho/software/Battenberg/data/GRCh38/dropbox_2022/"
G1000_LOCI_PREFIX="$BEAGLE_BASEDIR/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr"
G1000_ALLELES_PREFIX="$BEAGLE_BASEDIR/1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_allele_index_chr"

# Path to alleleCounter executable
ALLELECOUNTER="alleleCounter"

# Output directory
OUTPUT_DIR="/nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/ALPACA-WGS/data/output/test_bb_output"

# Script execution
echo "Running extract_snp_baf_logr.R..."
echo "Tumour BAM: ${TUMOUR_BAM}"
echo "Normal BAM: ${NORMAL_BAM}"
echo "Output Dir: ${OUTPUT_DIR}"
echo "AlleleCounter: ${ALLELECOUNTER}"
echo "1000 Genomes Loci Prefix: ${G1000_LOCI_PREFIX}"
echo "1000 Genomes Alleles Prefix: ${G1000_ALLELES_PREFIX}"


# GC Content and Replication Timing files
# Standard Battenberg reference files
GCCORRECT_PREFIX="$BEAGLE_BASEDIR/1000G_loci_hg38/1000_genomes_GC_corr_chr_"
REPLICCORRECT_PREFIX="$BEAGLE_BASEDIR/1000G_loci_hg38/1000_genomes_replication_timing_chr_"

Rscript extract_snp_baf_logr.R \
    --tumour_bam "${TUMOUR_BAM}" \
    --normal_bam "${NORMAL_BAM}" \
    --tumour_name "${TUMOUR_NAME}" \
    --normal_name "${NORMAL_NAME}" \
    --output_dir "${OUTPUT_DIR}" \
    --g1000_loci_prefix "${G1000_LOCI_PREFIX}" \
    --g1000_alleles_prefix "${G1000_ALLELES_PREFIX}" \
    --gccorrect_prefix "${GCCORRECT_PREFIX}" \
    --repliccorrect_prefix "${REPLICCORRECT_PREFIX}" \
    --allelecounter "${ALLELECOUNTER}" \
    --chromosomes "21" \
    --nthreads 4

if [ $? -eq 0 ]; then
    echo "Success! Output is in ${OUTPUT_DIR}"
else
    echo "Error: script failed (expected if paths are placeholders)"
fi
