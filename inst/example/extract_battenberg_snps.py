import pandas as pd
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(description="Extract per-SNP LogR values for Battenberg segments.")
    parser.add_argument("--subclones", required=True, help="Path to Battenberg subclones.txt file")
    parser.add_argument("--baf_segmented", required=True, help="Path to Battenberg .BAFsegmented.txt file")
    parser.add_argument("--logr_file", required=True, help="Path to Battenberg _mutantLogR_gcCorrected.tab file (or _mutantLogR.tab)")
    parser.add_argument("--output", required=True, help="Path to output tab-separated file")
    
    args = parser.parse_args()

    print(f"Reading subclones from: {args.subclones}")
    try:
        df_subclones = pd.read_csv(args.subclones, sep='\t')
    except Exception as e:
        print(f"Error reading subclones file: {e}")
        sys.exit(1)

    print(f"Reading BAF segmented from: {args.baf_segmented}")
    try:
        # BAFsegmented often has header: Chromosome Position BAF BAFphased BAFseg
        df_baf = pd.read_csv(args.baf_segmented, sep='\t')
    except Exception as e:
        print(f"Error reading BAF segmented file: {e}")
        sys.exit(1)

    print(f"Reading LogR from: {args.logr_file}")
    try:
        # LogR file usually has header: Chromosome Position LogR
        df_logr = pd.read_csv(args.logr_file, sep='\t')
    except Exception as e:
        print(f"Error reading LogR file: {e}")
        sys.exit(1)

    # Standardize column names if necessary (remove potential quotes/spaces)
    df_subclones.columns = df_subclones.columns.str.replace('"', '').str.strip()
    df_baf.columns = df_baf.columns.str.replace('"', '').str.strip()
    df_logr.columns = df_logr.columns.str.replace('"', '').str.strip()

    # 1. Filter LogR to only include SNPs that are in BAFsegmented (the Heterozygous ones used by Battenberg)
    print("Filtering LogR data to match BAFsegmented SNPs...")
    # Ensure keys are compatible (str vs int for Chromosome)
    df_baf['Chromosome'] = df_baf['Chromosome'].astype(str)
    df_baf['Position'] = df_baf['Position'].astype(int)
    
    df_logr['Chromosome'] = df_logr['Chromosome'].astype(str)
    df_logr['Position'] = df_logr['Position'].astype(int)

    # Inner join to keep only matching SNPs
    # We want columns: Chromosome, Position, LogR from df_logr, and maybe BAF/BAFseg from df_baf
    df_snps = pd.merge(
        df_logr[['Chromosome', 'Position', df_logr.columns[2]]], # Assume 3rd col is the sample LogR
        df_baf,
        on=['Chromosome', 'Position'],
        how='inner'
    )
    
    # Rename the LogR column to 'LogR' if it's not already
    logr_col = df_snps.columns[2] # The one from df_logr
    df_snps.rename(columns={logr_col: 'LogR'}, inplace=True)

    print(f"Retained {len(df_snps)} SNPs after filtering.")

    # 2. Map SNPs to segments from subclones.txt
    print("Mapping SNPs to segments...")
    
    # Create a list to store segment info for each SNP
    # Using a simple approach: For each chromosome, iterate segments and assign
    
    # Add segment columns
    df_snps['Segment_Start'] = -1
    df_snps['Segment_End'] = -1
    df_snps['Segment_LogR'] = float('nan')
    df_snps['Segment_CopyNumber'] = float('nan') # derived or placeholder

    # Optimize: Process by chromosome
    chromosomes = df_snps['Chromosome'].unique()
    
    result_dfs = []

    for chrom in chromosomes:
        # Get segments for this chrom
        # subclones 'chr' column might differ in type/formatting
        sub_chrom = df_subclones[df_subclones['chr'].astype(str) == chrom]
        
        if sub_chrom.empty:
            print(f"Warning: No segments found for chromosome {chrom} in subclones file.")
            continue
            
        snps_chrom = df_snps[df_snps['Chromosome'] == chrom].copy()
        
        # Sort for potential optimization (or just use simple checking)
        snps_chrom = snps_chrom.sort_values('Position')
        
        # IntervalIndex is fast
        # Create IntervalIndex from segments
        try:
            intervals = pd.IntervalIndex.from_arrays(
                sub_chrom['startpos'], 
                sub_chrom['endpos'], 
                closed='both' # Battenberg segments are inclusive? Usually yes.
            )
            
            # Map positions to intervals
            # pd.cut matches bins, but IntervalIndex.get_indexer is for matching
            # logic: snps_chrom['Position'] -> which interval?
            
            # Since intervals might be disjoint or not cover everything, we need safe handling
            # Using a apply approach is slow for millions of SNPs.
            # Using searchsorted or IntervalIndex get_indexer
            
            # Note: overlapping segments shouldn't happen in subclones usually
            if not intervals.is_non_overlapping_monotonic:
                 # Sort intervals
                 sub_chrom = sub_chrom.sort_values('startpos')
                 intervals = pd.IntervalIndex.from_arrays(sub_chrom['startpos'], sub_chrom['endpos'], closed='both')
            
            if pd.__version__ >= '1.0.0':
                 # Use get_indexer to find which interval each SNP pos falls into
                 indexer = intervals.get_indexer(snps_chrom['Position'])
                 
                 # indexer has -1 for no match
                 mask = indexer != -1
                 valid_indices = indexer[mask]
                 
                 # Assign values
                 snps_chrom.loc[mask, 'Segment_Start'] = sub_chrom.iloc[valid_indices]['startpos'].values
                 snps_chrom.loc[mask, 'Segment_End'] = sub_chrom.iloc[valid_indices]['endpos'].values
                 snps_chrom.loc[mask, 'Segment_LogR'] = sub_chrom.iloc[valid_indices]['LogR'].values
                 # Add ntot as basic copy number info
                 snps_chrom.loc[mask, 'Segment_CopyNumber'] = sub_chrom.iloc[valid_indices]['ntot'].values
                 
        except Exception as e:
            print(f"Error processing chromosome {chrom}: {e}")
            # Fallback to slower method if needed, but for now just print error
            pass
            
        result_dfs.append(snps_chrom)

    # Concatenate results
    if result_dfs:
        df_final = pd.concat(result_dfs)
        
        # Filter out SNPs that didn't match a segment (optional, but requested "SNPs... used by battenberg")
        matched_count = len(df_final[df_final['Segment_Start'] != -1])
        print(f"Mapped {matched_count} SNPs to segments.")
        
        df_final = df_final[df_final['Segment_Start'] != -1]
        
        print(f"Writing output to {args.output}")
        df_final.to_csv(args.output, sep='\t', index=False)
    else:
        print("No matching SNPs/Segments found.")

if __name__ == "__main__":
    main()
