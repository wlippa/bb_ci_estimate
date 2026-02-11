#!/bin/bash

# Configuration - Placeholders
# Replace these with actual paths to your data
TUMOUR_BAM="/path/to/tumour.bam"
NORMAL_BAM="/path/to/normal.bam"
TUMOUR_NAME="TumourSample"
NORMAL_NAME="NormalSample"

# Reference files (1000 Genomes)
# Battenberg requires 1000 Genomes loci and allele files split by chromosome
# Provide the prefix for these files
G1000_LOCI_PREFIX="/path/to/battenberg_reference/1000genomesloci2012_chr"
G1000_ALLELES_PREFIX="/path/to/battenberg_reference/1000genomesAlleles2012_chr"

# Path to alleleCounter executable
ALLELE_COUNTER="/path/to/alleleCounter"

# Output directory
OUTPUT_DIR="./battenberg_snp_baf_logr_output"

# Script execution
echo "Running extract_snp_baf_logr.R..."
echo "Tumour BAM: ${TUMOUR_BAM}"
echo "Normal BAM: ${NORMAL_BAM}"
echo "Output Dir: ${OUTPUT_DIR}"

Rscript extract_snp_baf_logr.R \
    --tumour_bam "${TUMOUR_BAM}" \
    --normal_bam "${NORMAL_BAM}" \
    --tumour_name "${TUMOUR_NAME}" \
    --normal_name "${NORMAL_NAME}" \
    --output_dir "${OUTPUT_DIR}" \
    --g1000_loci_prefix "${G1000_LOCI_PREFIX}" \
    --g1000_alleles_prefix "${G1000_ALLELES_PREFIX}" \
    --allelecounter "${ALLELE_COUNTER}" \
    --chromosomes "21,22" \
    --nthreads 4

if [ $? -eq 0 ]; then
    echo "Success! Output is in ${OUTPUT_DIR}"
else
    echo "Error: script failed (expected if paths are placeholders)"
fi
