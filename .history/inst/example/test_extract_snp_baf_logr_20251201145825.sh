#!/bin/bash
#
# Test script for extract_snp_baf_logr.R
# 
# This script creates mock data to test the BAF/LogR extraction
# without needing real BAM files or 1000 Genomes reference data.
#
# Usage: bash test_extract_snp_baf_logr.sh
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="${SCRIPT_DIR}/test_snp_extraction"
MOCK_DATA_DIR="${TEST_DIR}/mock_data"
OUTPUT_DIR="${TEST_DIR}/output"

echo "=== Testing BAF/LogR SNP Extraction ==="
echo "Test directory: ${TEST_DIR}"
echo ""

# Clean up previous test runs
rm -rf "${TEST_DIR}"
mkdir -p "${MOCK_DATA_DIR}"
mkdir -p "${OUTPUT_DIR}"

# ============================================================================
# Generate mock 1000 Genomes reference files
# Format: position ref_allele_index alt_allele_index
# Indices: A=1, C=2, G=3, T=4
# ============================================================================
echo "Creating mock 1000 Genomes reference files..."

# Chromosome 1 - 100 SNPs
cat > "${MOCK_DATA_DIR}/1000genomesAlleles_chr1.txt" << 'EOF'
1000	1	2
2000	3	4
3000	1	3
4000	2	4
5000	1	4
6000	2	3
7000	1	2
8000	3	4
9000	1	3
10000	2	4
11000	1	4
12000	2	3
13000	1	2
14000	3	4
15000	1	3
16000	2	4
17000	1	4
18000	2	3
19000	1	2
20000	3	4
EOF

# Chromosome 2 - 100 SNPs
cat > "${MOCK_DATA_DIR}/1000genomesAlleles_chr2.txt" << 'EOF'
1500	1	2
2500	3	4
3500	1	3
4500	2	4
5500	1	4
6500	2	3
7500	1	2
8500	3	4
9500	1	3
10500	2	4
11500	1	4
12500	2	3
13500	1	2
14500	3	4
15500	1	3
16500	2	4
17500	1	4
18500	2	3
19500	1	2
20500	3	4
EOF

# 1000 Genomes loci files (for alleleCounter - just positions)
cat > "${MOCK_DATA_DIR}/1000genomesloci_chr1.txt" << 'EOF'
1	1000
1	2000
1	3000
1	4000
1	5000
1	6000
1	7000
1	8000
1	9000
1	10000
1	11000
1	12000
1	13000
1	14000
1	15000
1	16000
1	17000
1	18000
1	19000
1	20000
EOF

cat > "${MOCK_DATA_DIR}/1000genomesloci_chr2.txt" << 'EOF'
2	1500
2	2500
2	3500
2	4500
2	5500
2	6500
2	7500
2	8500
2	9500
2	10500
2	11500
2	12500
2	13500
2	14500
2	15500
2	16500
2	17500
2	18500
2	19500
2	20500
EOF

# ============================================================================
# Generate mock allele count files (simulating alleleCounter output)
# Format: CHR POS Count_A Count_C Count_G Count_T Depth
# This allows us to skip the BAM processing step
# ============================================================================
echo "Creating mock allele count files..."

# Function to generate mock allele counts
# Simulating a tumour with some copy number changes:
# - Chr1:1-10000: Normal (diploid, BAF ~0.5)
# - Chr1:10001-20000: CN gain (BAF shift towards 0.33 or 0.67)
# - Chr2: LOH region (BAF ~0 or ~1)

# Normal sample - mostly heterozygous (50/50)
cat > "${MOCK_DATA_DIR}/normal_alleleFrequencies_chr1.txt" << 'EOF'
1	1000	25	27	0	0	52
1	2000	0	0	24	28	52
1	3000	26	0	25	0	51
1	4000	0	24	0	26	50
1	5000	27	0	0	25	52
1	6000	0	26	24	0	50
1	7000	24	28	0	0	52
1	8000	0	0	25	26	51
1	9000	27	0	24	0	51
1	10000	0	25	0	27	52
1	11000	26	0	0	25	51
1	12000	0	27	24	0	51
1	13000	25	26	0	0	51
1	14000	0	0	26	25	51
1	15000	24	0	27	0	51
1	16000	0	26	0	25	51
1	17000	25	0	0	26	51
1	18000	0	25	27	0	52
1	19000	26	25	0	0	51
1	20000	0	0	25	26	51
EOF

cat > "${MOCK_DATA_DIR}/normal_alleleFrequencies_chr2.txt" << 'EOF'
2	1500	25	26	0	0	51
2	2500	0	0	24	27	51
2	3500	26	0	25	0	51
2	4500	0	25	0	26	51
2	5500	27	0	0	24	51
2	6500	0	26	25	0	51
2	7500	25	26	0	0	51
2	8500	0	0	26	25	51
2	9500	25	0	26	0	51
2	10500	0	26	0	25	51
2	11500	25	0	0	26	51
2	12500	0	25	26	0	51
2	13500	26	25	0	0	51
2	14500	0	0	25	26	51
2	15500	25	0	26	0	51
2	16500	0	26	0	25	51
2	17500	26	0	0	25	51
2	18500	0	25	26	0	51
2	19500	25	26	0	0	51
2	20500	0	0	26	25	51
EOF

# Tumour sample - shows copy number changes
# Chr1:1-10000: diploid (similar to normal)
# Chr1:10001-20000: gain with BAF shift (e.g., 3+1 copies)
# Chr2: LOH (allelic imbalance)
cat > "${MOCK_DATA_DIR}/tumour_alleleFrequencies_chr1.txt" << 'EOF'
1	1000	48	52	0	0	100
1	2000	0	0	47	53	100
1	3000	51	0	49	0	100
1	4000	0	48	0	52	100
1	5000	53	0	0	47	100
1	6000	0	51	49	0	100
1	7000	47	53	0	0	100
1	8000	0	0	52	48	100
1	9000	49	0	51	0	100
1	10000	0	50	0	50	100
1	11000	75	0	0	25	100
1	12000	0	72	28	0	100
1	13000	74	26	0	0	100
1	14000	0	0	73	27	100
1	15000	71	0	29	0	100
1	16000	0	76	0	24	100
1	17000	73	0	0	27	100
1	18000	0	74	26	0	100
1	19000	72	28	0	0	100
1	20000	0	0	75	25	100
EOF

# Chr2: Strong LOH (one allele dominant)
cat > "${MOCK_DATA_DIR}/tumour_alleleFrequencies_chr2.txt" << 'EOF'
2	1500	92	8	0	0	100
2	2500	0	0	7	93	100
2	3500	91	0	9	0	100
2	4500	0	6	0	94	100
2	5500	93	0	0	7	100
2	6500	0	8	92	0	100
2	7500	94	6	0	0	100
2	8500	0	0	91	9	100
2	9500	90	0	10	0	100
2	10500	0	93	0	7	100
2	11500	95	0	0	5	100
2	12500	0	8	92	0	100
2	13500	91	9	0	0	100
2	14500	0	0	6	94	100
2	15500	93	0	7	0	100
2	16500	0	92	0	8	100
2	17500	90	0	0	10	100
2	18500	0	9	91	0	100
2	19500	94	6	0	0	100
2	20500	0	0	92	8	100
EOF

echo "Mock data created!"
echo ""

# ============================================================================
# Run the extraction script with mock data
# ============================================================================
echo "Running BAF/LogR extraction..."
echo ""

cd "${OUTPUT_DIR}"

Rscript "${SCRIPT_DIR}/extract_snp_baf_logr.R" \
    --skip_allele_counting \
    --tumour_name "tumour" \
    --normal_name "normal" \
    --g1000_loci_prefix "${MOCK_DATA_DIR}/1000genomesloci_chr" \
    --g1000_alleles_prefix "${MOCK_DATA_DIR}/1000genomesAlleles_chr" \
    --chromosomes "1,2" \
    --min_depth 10 \
    --output_dir "${OUTPUT_DIR}" \
    --output_prefix "test_sample" \
    --seed 42

echo ""
echo "=== Output Files ==="
ls -la "${OUTPUT_DIR}"

echo ""
echo "=== Preview of combined output (first 25 lines) ==="
head -25 "${OUTPUT_DIR}/test_sample_per_snp.tsv"

echo ""
echo "=== Preview of BAF file ==="
head -15 "${OUTPUT_DIR}/test_sample_tumour_BAF.tsv"

echo ""
echo "=== Preview of LogR file ==="
head -15 "${OUTPUT_DIR}/test_sample_tumour_LogR.tsv"

echo ""
echo "=== Validation ==="
# Check that output files exist and have expected structure
SNP_COUNT=$(wc -l < "${OUTPUT_DIR}/test_sample_per_snp.tsv")
echo "Total lines in output (including header): ${SNP_COUNT}"

# Check BAF values for chr2 (should show LOH - values near 0 or 1)
echo ""
echo "Chr2 BAF values (should show LOH pattern - near 0 or near 1):"
grep "^2" "${OUTPUT_DIR}/test_sample_per_snp.tsv" | awk '{print $1, $2, $4}' | head -5

# Check chr1 first half (should be ~0.5) vs second half (should show shift)
echo ""
echo "Chr1 early positions (should be ~0.5 BAF):"
grep "^1" "${OUTPUT_DIR}/test_sample_per_snp.tsv" | head -5 | awk '{print $1, $2, $4}'

echo ""
echo "Chr1 later positions (should show BAF shift due to CN gain):"
grep "^1" "${OUTPUT_DIR}/test_sample_per_snp.tsv" | tail -5 | awk '{print $1, $2, $4}'

echo ""
echo "=== Test Complete ==="
echo "All output files are in: ${OUTPUT_DIR}"
