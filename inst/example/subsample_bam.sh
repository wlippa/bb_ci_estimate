#!/bin/bash

# subsample_bam.sh
# Subsample a BAM file by chromosome and optional coordinate range.
# Usage: ./subsample_bam.sh -i <input.bam> -o <output.bam> -c <chromosome> [-s <start>] [-e <end>]
# bash /nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/ALPACA-WGS/repos/bb_ci_estimate/inst/example/subsample_bam.sh -i /nemo/project/proj-tracerx-wgs/working/release/ALIGNMENT/LTX1231_SU_T1-R1--8c625f964d/LTX1231_SU_T1-R1--8c625f964d.bqsr.bam -o /nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/ALPACA-WGS/data/bams/test.bqsr.bam -c chr21 -s 20000000 -e 21000000
# bash /nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/ALPACA-WGS/repos/bb_ci_estimate/inst/example/subsample_bam.sh -i /nemo/project/proj-tracerx-wgs/working/release/ALIGNMENT/LTX1231_BS_GL--59c4c25703/LTX1231_BS_GL--59c4c25703.bqsr.bam -o /nemo/project/proj-tracerx-lung/tctProjects/CN-CCF/ALPACA-WGS/data/bams/normal.bqsr.bam -c chr21 -s 20000000 -e 21000000



usage() {
    echo "Usage: $0 -i <input.bam> -o <output.bam> -c <chromosome> [-s <start>] [-e <end>]"
    echo "  -i  Input BAM file"
    echo "  -o  Output BAM file"
    echo "  -c  Chromosome to extract"
    echo "  -s  Start position (optional)"
    echo "  -e  End position (optional)"
    exit 1
}

INPUT_BAM=""
OUTPUT_BAM=""
CHROMOSOME=""
START_POS=""
END_POS=""

while getopts "i:o:c:s:e:" opt; do
    case $opt in
        i) INPUT_BAM=$OPTARG ;;
        o) OUTPUT_BAM=$OPTARG ;;
        c) CHROMOSOME=$OPTARG ;;
        s) START_POS=$OPTARG ;;
        e) END_POS=$OPTARG ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT_BAM" || -z "$OUTPUT_BAM" || -z "$CHROMOSOME" ]]; then
    echo "Error: Missing required arguments."
    usage
fi

REGION="$CHROMOSOME"
if [[ -n "$START_POS" ]]; then
    REGION="${REGION}:${START_POS}"
    if [[ -n "$END_POS" ]]; then
        REGION="${REGION}-${END_POS}"
    fi
elif [[ -n "$END_POS" ]]; then
     # samtools view region format is chr:start-end. If start is missing but end is present, it's ambiguous/invalid in this simple construction scheme without start.
     # Assuming if end is provided, start should be provided or default to 1.
     # Let's enforce start if end is present for simplicity, or default start to 1.
     REGION="${REGION}:1-${END_POS}"
fi

echo "Subsampling BAM..."
echo "Input: $INPUT_BAM"
echo "Output: $OUTPUT_BAM"
echo "Region: $REGION"

if samtools view -b "$INPUT_BAM" "$REGION" > "$OUTPUT_BAM"; then
    echo "Subsampling complete."
    echo "Indexing output BAM..."
    samtools index "$OUTPUT_BAM"
    echo "Done."
else
    echo "Error: samtools view failed."
    exit 1
fi
