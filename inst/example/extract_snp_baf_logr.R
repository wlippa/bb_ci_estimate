#!/usr/bin/env Rscript
#' Extract per-SNP BAF and LogR values from WGS data
#' 
#' This script extracts B-Allele Frequency (BAF) and Log R Ratio (LogR) values
#' at the SNP level, matching Battenberg's internal. Battenberg does not output these
#' values, but we need to them to estimate confidence intervals for copy  number calls.
#' 
#' The script performs the following steps:
#' 1. Counts alleles at 1000 Genomes loci using alleleCounter
#' 2. Calculates BAF and LogR per SNP position
#' 3. Outputs a combined file with per-SNP values
#'
#' Usage:
#'   Rscript extract_snp_baf_logr.R [options]
#'
#' Author: Based on Battenberg (Wedge-lab)
#' License: GPL-3
#' 
#' Required libraries: optparse, readr, parallel, splines

library(optparse)
library(splines)

option_list <- list(
  make_option(c("-t", "--tumour_bam"), type = "character", default = NULL,
              help = "Path to tumour BAM file", metavar = "FILE"),
  make_option(c("-n", "--normal_bam"), type = "character", default = NULL,
              help = "Path to normal BAM file", metavar = "FILE"),
  make_option(c("--tumour_name"), type = "character", default = "tumour",
              help = "Sample name for tumour [default: %default]", metavar = "STRING"),
  make_option(c("--normal_name"), type = "character", default = "normal",
              help = "Sample name for normal [default: %default]", metavar = "STRING"),
  make_option(c("-o", "--output_dir"), type = "character", default = ".",
              help = "Output directory [default: %default]", metavar = "DIR"),
  make_option(c("--output_prefix"), type = "character", default = "snp_baf_logr",
              help = "Prefix for output files [default: %default]", metavar = "STRING"),
  make_option(c("-g", "--g1000_loci_prefix"), type = "character", default = NULL,
              help = "Prefix to 1000 Genomes loci files (e.g., path/to/1000genomesloci2012_chr)", metavar = "PATH"),
  make_option(c("-a", "--g1000_alleles_prefix"), type = "character", default = NULL,
              help = "Prefix to 1000 Genomes allele files (e.g., path/to/1000genomesAlleles2012_chr)", metavar = "PATH"),
  make_option(c("-c", "--chromosomes"), type = "character", default = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22",
              help = "Comma-separated list of chromosomes [default: %default]", metavar = "STRING"),
  make_option(c("--gccorrect_prefix"), type = "character", default = NULL,
              help = "Prefix to GC content files [default: %default]", metavar = "PATH"),
  make_option(c("--repliccorrect_prefix"), type = "character", default = NULL,
              help = "Prefix to replication timing files [default: %default]", metavar = "PATH"),
  make_option(c("--min_base_qual"), type = "integer", default = 20,
              help = "Minimum base quality [default: %default]", metavar = "INT"),
  make_option(c("--min_map_qual"), type = "integer", default = 35,
              help = "Minimum mapping quality [default: %default]", metavar = "INT"),
  make_option(c("--min_depth"), type = "integer", default = 10,
              help = "Minimum read depth in normal [default: %default]", metavar = "INT"),
  make_option(c("--allelecounter"), type = "character", default = "alleleCounter",
              help = "Path to alleleCounter executable [default: %default]", metavar = "PATH"),
  make_option(c("--skip_allele_counting"), action = "store_true", default = FALSE,
              help = "Skip allele counting step (expects files to exist)"),
  make_option(c("--nthreads"), type = "integer", default = 1,
              help = "Number of threads for parallel processing [default: %default]", metavar = "INT"),
  make_option(c("--seed"), type = "integer", default = NULL,
              help = "Random seed for reproducibility [default: current time]", metavar = "INT")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Extract per-SNP BAF and LogR values from WGS tumour/normal pairs")
opt <- parse_args(opt_parser)

# Validate required arguments
if (!opt$skip_allele_counting) {
  if (is.null(opt$tumour_bam) || is.null(opt$normal_bam)) {
    stop("Both --tumour_bam and --normal_bam are required unless --skip_allele_counting is set")
  }
}
if (is.null(opt$g1000_loci_prefix) || is.null(opt$g1000_alleles_prefix)) {
  stop("Both --g1000_loci_prefix and --g1000_alleles_prefix are required")
}

# Parse chromosomes
chr_names <- strsplit(opt$chromosomes, ",")[[1]]

# Set seed
if (is.null(opt$seed)) {
  opt$seed <- as.integer(Sys.time())
}
set.seed(opt$seed)

# Create output directory
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}
setwd(opt$output_dir)

cat("=== Extract SNP-level BAF and LogR ===\n")
cat("Tumour BAM:", opt$tumour_bam, "\n")
cat("Normal BAM:", opt$normal_bam, "\n")
cat("Chromosomes:", opt$chromosomes, "\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Random seed:", opt$seed, "\n")
cat("======================================\n\n")


#' Generic table reader (from Battenberg util.R)
read_table_generic <- function(file, header = TRUE, sep = "\t", chrom_col = 1, skip = 0) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    # Fallback to base R
    d <- read.table(file, header = header, sep = sep, stringsAsFactors = FALSE,
                    skip = skip, comment.char = "")
    d[, chrom_col] <- as.character(d[, chrom_col])
    return(d)
  }
  
  d <- readr::read_delim(file = file, delim = sep, col_names = header, n_max = 1, 
                         skip = skip, col_types = readr::cols())
  col_types <- list()
  for (i in chrom_col) {
    first_colname <- colnames(d)[i]
    col_types[[first_colname]] <- readr::col_character()
  }
  d <- readr::read_delim(file = file, delim = sep, col_names = header, 
                         col_types = col_types, skip = skip)
  colnames(d) <- gsub(" ", ".", colnames(d))
  return(as.data.frame(d))
}


#' Parser for GC content reference data (from Battenberg util.R)
read_gccontent <- function(filename) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    # Fallback to base R if readr is not available (though it should be)
    # GC content format: chr, position, then many columns of GC content
    # We skip 1 line (header is usually not standard)
    d <- read.table(filename, skip = 1, header = FALSE, stringsAsFactors = FALSE)
    return(d)
  }
  return(readr::read_tsv(file=filename, skip = 1, col_names = FALSE, col_types = "-cinnnnnnnnnnnn------"))
}

#' Parser for replication timing reference data (from Battenberg util.R)
read_replication <- function(filename) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    d <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
    return(d)
  }
  return(readr::read_tsv(file=filename, col_types = paste0("ci", paste0(rep("n", 15), collapse = ""))))
}

#' Function to correct LogR for waviness that correlates with GC content
#' @param snp_data Data frame containing Tumour_LogR, Chromosome, and Position
#' @param gc_content_file_prefix Prefix to GC content files
#' @param replic_timing_file_prefix Prefix to replication timing files (optional)
#' @param chr_names Vector of chromosome names
#' @return Corrected snp_data
gc_correct <- function(snp_data, gc_content_file_prefix, replic_timing_file_prefix, chr_names) {
  
  if (is.null(gc_content_file_prefix)) {
    cat("Warning: No GC content prefix provided. Skipping GC correction.\n")
    return(snp_data)
  }
  
  cat("Performing GC correction...\n")
  
  Tumor_LogR <- data.frame(Chromosome = snp_data$Chromosome, Position = snp_data$Position, LogR = snp_data$Tumour_LogR, stringsAsFactors = FALSE)
  
  # Read GC content
  cat("  Reading GC content data...\n")
  gc_files <- paste0(gc_content_file_prefix, chr_names, ".txt.gz")
  
  # Check if files exist
  missing_files <- gc_files[!file.exists(gc_files)]
  if (length(missing_files) > 0) {
    # Try without .gz
    gc_files_nogz <- paste0(gc_content_file_prefix, chr_names, ".txt")
    if (all(file.exists(gc_files_nogz))) {
      gc_files <- gc_files_nogz
    } else {
      stop(paste("GC content files not found. Checked .txt.gz and .txt. Missing:", head(missing_files, 1)))
    }
  }
  
  GC_data <- do.call(rbind, lapply(gc_files, read_gccontent))
  # Assuming standard Battenberg GC file format
  colnames(GC_data) <- c("chr", "Position", paste0(c(25,50,100,200,500), "bp"),
                         paste0(c(1,2,5,10,20,50,100), "kb"))
  
  replic_data <- NULL
  if (!is.null(replic_timing_file_prefix)) {
    cat("  Reading replication timing data...\n")
    replic_files <- paste0(replic_timing_file_prefix, chr_names, ".txt.gz")
    
    # Check if files exist
    missing_files <- replic_files[!file.exists(replic_files)]
    if (length(missing_files) > 0) {
      # Try without .gz
      replic_files_nogz <- paste0(replic_timing_file_prefix, chr_names, ".txt")
      if (all(file.exists(replic_files_nogz))) {
        replic_files <- replic_files_nogz
      } else {
         cat("Warning: Replication timing files not found. Skipping replication correction.\n")
         replic_timing_file_prefix <- NULL
      }
    }
    
    if (!is.null(replic_timing_file_prefix)) {
      replic_data <- do.call(rbind, lapply(replic_files, read_replication))
    }
  }
  
  # Match loci
  # Note: Battenberg's gc.correct.wgs filters Tumor_LogR to match GC_data. 
  # Here we want to keep all SNPs, but we can only correct those that overlap.
  # However, typically GC data covers the whole genome or at least the 1000G loci.
  
  cat("  Matching loci...\n")
  # Use keys for matching
  key_logr <- paste(Tumor_LogR$Chromosome, Tumor_LogR$Position, sep = "_")
  key_gc <- paste(GC_data$chr, GC_data$Position, sep = "_")
  
  common_keys <- intersect(key_logr, key_gc)
  cat("  Overlapping loci:", length(common_keys), "of", nrow(Tumor_LogR), "\n")
  
  if (length(common_keys) < 100) {
    warning("Too few overlapping loci for GC correction. Skipping.")
    return(snp_data)
  }
  
  # Subset to common loci for model fitting
  idx_logr <- match(common_keys, key_logr)
  idx_gc <- match(common_keys, key_gc)
  
  sub_logr <- Tumor_LogR[idx_logr, ]
  sub_gc <- GC_data[idx_gc, ]
  
  if (!is.null(replic_data)) {
    key_rep <- paste(replic_data$chr, replic_data$Position, sep = "_")
    idx_rep <- match(common_keys, key_rep)
    sub_rep <- replic_data[idx_rep, ]
    # Ensure no NAs in replication data if we are using it
    valid_rep <- !is.na(idx_rep)
    sub_logr <- sub_logr[valid_rep, ]
    sub_gc <- sub_gc[valid_rep, ]
    sub_rep <- sub_rep[valid_rep[valid_rep], ]
  }
  
  # Calculate correlations
  corr <- abs(cor(sub_gc[, 3:ncol(sub_gc)], sub_logr[,3], use="complete.obs")[,1])
  
  index_1kb <- which(names(corr)=="1kb")
  maxGCcol_insert <- names(which.max(corr[1:index_1kb]))
  index_100kb <- which(names(corr)=="100kb")
  maxGCcol_amplic <- names(which.max(corr[(index_1kb+2):index_100kb]))
  
  cat("  Max GC correlation (short):", maxGCcol_insert, "-", format(max(corr[1:index_1kb]), digits=3), "\n")
  cat("  Max GC correlation (long):", maxGCcol_amplic, "-", format(max(corr[(index_1kb+2):index_100kb]), digits=3), "\n")
  
  corrdata <- data.frame(logr = sub_logr[,3],
                         GC_insert = sub_gc[,maxGCcol_insert],
                         GC_amplic = sub_gc[,maxGCcol_amplic])
  
  model_formula <- logr ~ ns(x = GC_insert, df = 5, intercept = T) + ns(x = GC_amplic, df = 5, intercept = T)
  
  if (!is.null(replic_data)) {
    corr_rep <- abs(cor(sub_rep[, 3:ncol(sub_rep)], sub_logr[,3], use="complete.obs")[,1])
    maxreplic <- names(which.max(corr_rep))
    cat("  Max Replication correlation:", maxreplic, "-", format(max(corr_rep), digits=3), "\n")
    
    corrdata$replic <- sub_rep[, maxreplic]
    model_formula <- logr ~ ns(x = GC_insert, df = 5, intercept = T) + ns(x = GC_amplic, df = 5, intercept = T) + ns(x = replic, df = 5, intercept = T)
  }
  
  # Fit model
  cat("  Fitting regression model...\n")
  model <- lm(model_formula, data = corrdata, na.action="na.exclude")
  
  # Predict correction for all sites (residuals = observed - predicted)
  # predicting on the original full dataset (where possible)
  # But we need GC values for ALL sites in snp_data to predict.
  
  # We will map the GC values back to the full snp_data
  # If a SNP is missing from GC data, we can't correct it easily (or set specific NA).
  # We'll use constraints of common_keys for updating.
  
  residuals <- residuals(model)
  
  # Update the LogR values in the main object
  # Only for the ones we used in the model
  
  # Corrected LogR = Residuals
  # Map back to original indices
  snp_data$Tumour_LogR[idx_logr] <- residuals
  
  # Mark non-corrected as NA or keep raw? Battenberg typically outputs NA if not corrected or filtered.
  # But here we might have filtered down to common keys. existing code likely wants values.
  # Let's keep raw values for those we couldn't correct (if any), but warn.
  
  if (nrow(snp_data) > length(idx_logr)) {
    cat("  Warning:", nrow(snp_data) - length(idx_logr), "SNPs could not be GC corrected (missing in reference files). Keeping raw LogR for these.\n")
  }
  
  return(snp_data)
}

#' Concatenate allele count files (from Battenberg util.R)
concatenateAlleleCountFiles <- function(inputStart, inputEnd, chr_names) {
  infiles <- c()
  for (chrom in chr_names) {
    filename <- paste0(inputStart, chrom, inputEnd)
    if (file.exists(filename) && file.info(filename)$size > 0) {
      infiles <- c(infiles, filename)
    }
  }
  return(as.data.frame(do.call(rbind, lapply(infiles, FUN = function(x) {
    read_table_generic(x)
  }))))
}


#' Concatenate 1000 Genomes SNP reference files (from Battenberg util.R)
concatenateG1000SnpFiles <- function(inputStart, inputEnd, chr_names) {
  data <- list()
  for (chrom in chr_names) {
    filename <- paste0(inputStart, chrom, inputEnd)
    if (file.exists(filename) && file.info(filename)$size > 0) {
      data[[chrom]] <- cbind(chromosome = chrom, read_table_generic(filename))
    }
  }
  return(as.data.frame(do.call(rbind, data)))
}


#' Get allele counts using alleleCounter (from Battenberg prepare_wgs.R)
getAlleleCounts <- function(bam.file, output.file, g1000.loci, 
                            min.base.qual = 20, min.map.qual = 35, 
                            allelecounter.exe = "alleleCounter") {
  cmd <- paste(allelecounter.exe,
               "-b", bam.file,
               "-l", g1000.loci,
               "-o", output.file,
               "-m", min.base.qual,
               "-q", min.map.qual)
  
  # Check alleleCounter version for dense-snp mode
  counter_version <- tryCatch({
    system(paste(allelecounter.exe, "--version"), intern = TRUE)
  }, error = function(e) "0")
  
  if (length(counter_version) > 0 && as.integer(substr(counter_version, 1, 1)) >= 4) {
    cmd <- paste(cmd, "--dense-snps")
  }
  
  cat("Running:", cmd, "\n")
  exit_code <- system(cmd, wait = TRUE)
  if (exit_code != 0) {
    stop(paste("alleleCounter failed with exit code", exit_code))
  }
}


#' Calculate BAF and LogR from allele counts (adapted from Battenberg getBAFsAndLogRs)
#' This is the core function that extracts per-SNP BAF and LogR values
calculateBAFandLogR <- function(tumourAlleleCountsFile.prefix,
                                 normalAlleleCountsFile.prefix,
                                 g1000file.prefix,
                                 chr_names,
                                 minCounts = NA,
                                 samplename = "sample1",
                                 seed = as.integer(Sys.time())) {
  
  set.seed(seed)
  
  cat("Reading allele count files...\n")
  input_data <- concatenateAlleleCountFiles(tumourAlleleCountsFile.prefix, ".txt", chr_names)
  normal_input_data <- concatenateAlleleCountFiles(normalAlleleCountsFile.prefix, ".txt", chr_names)
  allele_data <- concatenateG1000SnpFiles(g1000file.prefix, ".txt", chr_names)
  
  cat("  Tumour SNPs:", nrow(input_data), "\n")
  cat("  Normal SNPs:", nrow(normal_input_data), "\n")
  cat("  Reference SNPs:", nrow(allele_data), "\n")
  
  # Strip "chr" prefix for matching
  allele_data[, 1] <- gsub("chr", "", allele_data[, 1])
  normal_input_data[, 1] <- gsub("chr", "", normal_input_data[, 1])
  input_data[, 1] <- gsub("chr", "", input_data[, 1])
  
  # Synchronise all data frames
  chrpos_allele <- paste(allele_data[, 1], "_", allele_data[, 2], sep = "")
  chrpos_normal <- paste(normal_input_data[, 1], "_", normal_input_data[, 2], sep = "")
  chrpos_tumour <- paste(input_data[, 1], "_", input_data[, 2], sep = "")
  matched_data <- Reduce(intersect, list(chrpos_allele, chrpos_normal, chrpos_tumour))
  
  cat("  Matched SNPs:", length(matched_data), "\n")
  
  allele_data <- allele_data[chrpos_allele %in% matched_data, ]
  normal_input_data <- normal_input_data[chrpos_normal %in% matched_data, ]
  input_data <- input_data[chrpos_tumour %in% matched_data, ]
  
  # Clean up column names
  names(input_data)[1] <- "CHR"
  names(normal_input_data)[1] <- "CHR"
  
  normal_data <- normal_input_data[, 3:6]
  mutant_data <- input_data[, 3:6]
  
  # Calculate allele depths
  len <- nrow(normal_data)
  normCount1 <- normal_data[cbind(1:len, allele_data[, 3])]
  normCount2 <- normal_data[cbind(1:len, allele_data[, 4])]
  totalNormal <- normCount1 + normCount2
  mutCount1 <- mutant_data[cbind(1:len, allele_data[, 3])]
  mutCount2 <- mutant_data[cbind(1:len, allele_data[, 4])]
  totalMutant <- mutCount1 + mutCount2
  
  # Clean up
  rm(normal_data, mutant_data, allele_data, normal_input_data)
  
  # Filter by minimum depth
  indices <- 1:nrow(input_data)
  if (!is.na(minCounts)) {
    cat("Filtering by minimum depth:", minCounts, "\n")
    indices <- which(totalNormal >= minCounts & totalMutant >= 1)
    
    totalNormal <- totalNormal[indices]
    totalMutant <- totalMutant[indices]
    normCount1 <- normCount1[indices]
    normCount2 <- normCount2[indices]
    mutCount1 <- mutCount1[indices]
    mutCount2 <- mutCount2[indices]
    cat("  SNPs after filtering:", length(indices), "\n")
  }
  n <- length(indices)
  
  # Allocate vectors
  normalBAF <- vector(length = n, mode = "numeric")
  mutantBAF <- vector(length = n, mode = "numeric")
  normalLogR <- vector(length = n, mode = "numeric")
  mutantLogR <- vector(length = n, mode = "numeric")
  
  # Randomise A and B alleles (as in Battenberg)
  selector <- round(runif(n))
  normalBAF[selector == 0] <- normCount1[selector == 0] / totalNormal[selector == 0]
  normalBAF[selector == 1] <- normCount2[selector == 1] / totalNormal[selector == 1]
  mutantBAF[selector == 0] <- mutCount1[selector == 0] / totalMutant[selector == 0]
  mutantBAF[selector == 1] <- mutCount2[selector == 1] / totalMutant[selector == 1]
  
  # LogR calculation (normalised to normal)
  normalLogR <- vector(length = n, mode = "integer")  # Normal LogR = 0
  mutantLogR_raw <- totalMutant / totalNormal
  mutantLogR <- log2(mutantLogR_raw / mean(mutantLogR_raw, na.rm = TRUE))
  
  # Create output data frame
  result <- data.frame(
    Chromosome = input_data$CHR[indices],
    Position = input_data$POS[indices],
    Normal_BAF = normalBAF,
    Tumour_BAF = mutantBAF,
    Normal_LogR = normalLogR,
    Tumour_LogR = mutantLogR,
    Normal_Depth = totalNormal,
    Tumour_Depth = totalMutant,
    Normal_Count_Allele1 = normCount1,
    Normal_Count_Allele2 = normCount2,
    Tumour_Count_Allele1 = mutCount1,
    Tumour_Count_Allele2 = mutCount2
  )
  
  return(result)
}


# ==============================================================================
# Main execution
# ==============================================================================
# Step 1: Run allele counting (if not skipped)
browser()
if (!opt$skip_allele_counting) {
  cat("Step 1: Running allele counting...\n")
  if (opt$nthreads > 1 && requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(opt$nthreads)
    parallel::clusterExport(cl, varlist = c("getAlleleCounts", "opt"), envir = environment())
    parallel::parLapply(cl, chr_names, function(chrom) {
      # Tumour
      getAlleleCounts(
        bam.file = opt$tumour_bam,
        output.file = paste0(opt$tumour_name, "_alleleFrequencies_chr", chrom, ".txt"),
        g1000.loci = paste0(opt$g1000_loci_prefix, chrom, ".txt"),
        min.base.qual = opt$min_base_qual,
        min.map.qual = opt$min_map_qual,
        allelecounter.exe = opt$allelecounter
      )
      # Normal
      getAlleleCounts(
        bam.file = opt$normal_bam,
        output.file = paste0(opt$normal_name, "_alleleFrequencies_chr", chrom, ".txt"),
        g1000.loci = paste0(opt$g1000_loci_prefix, chrom, ".txt"),
        min.base.qual = opt$min_base_qual,
        min.map.qual = opt$min_map_qual,
        allelecounter.exe = opt$allelecounter
      )
    })
    parallel::stopCluster(cl)
  } else {
    for (chrom in chr_names) {
      cat("  Chromosome", chrom, "...\n")
      # Tumour
      getAlleleCounts(
        bam.file = opt$tumour_bam,
        output.file = paste0(opt$tumour_name, "_alleleFrequencies_chr", chrom, ".txt"),
        g1000.loci = paste0(opt$g1000_loci_prefix, chrom, ".txt"),
        min.base.qual = opt$min_base_qual,
        min.map.qual = opt$min_map_qual,
        allelecounter.exe = opt$allelecounter
      )
      # Normal
      getAlleleCounts(
        bam.file = opt$normal_bam,
        output.file = paste0(opt$normal_name, "_alleleFrequencies_chr", chrom, ".txt"),
        g1000.loci = paste0(opt$g1000_loci_prefix, chrom, ".txt"),
        min.base.qual = opt$min_base_qual,
        min.map.qual = opt$min_map_qual,
        allelecounter.exe = opt$allelecounter
      )
    }
  }
} else {
  cat("Step 1: Skipping allele counting (using existing files)\n")
}

# Step 2: Calculate BAF and LogR
cat("\nStep 2: Calculating BAF and LogR per SNP...\n")

snp_data <- calculateBAFandLogR(
  tumourAlleleCountsFile.prefix = paste0(opt$tumour_name, "_alleleFrequencies_chr"),
  normalAlleleCountsFile.prefix = paste0(opt$normal_name, "_alleleFrequencies_chr"),
  g1000file.prefix = opt$g1000_alleles_prefix,
  chr_names = chr_names,
  minCounts = opt$min_depth,
  samplename = opt$tumour_name,
  seed = opt$seed
)

# Step 2b: perform GC correction
if (!is.null(opt$gccorrect_prefix)) {
  snp_data <- gc_correct(
    snp_data = snp_data,
    gc_content_file_prefix = opt$gccorrect_prefix,
    replic_timing_file_prefix = opt$repliccorrect_prefix,
    chr_names = chr_names
  )
}

# Step 3: Write output
cat("\nStep 3: Writing output files...\n")

# Combined output file with all data
output_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_per_snp.tsv"))
write.table(snp_data, file = output_file, row.names = FALSE, quote = FALSE, sep = "\t")
cat("  Combined output:", output_file, "\n")

# Also write separate BAF and LogR files (matching Battenberg format)
baf_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_tumour_BAF.tsv"))
logr_file <- file.path(opt$output_dir, paste0(opt$output_prefix, "_tumour_LogR.tsv"))

write.table(
  data.frame(Chromosome = snp_data$Chromosome, Position = snp_data$Position, 
             BAF = snp_data$Tumour_BAF),
  file = baf_file, row.names = FALSE, quote = FALSE, sep = "\t",
  col.names = c("Chromosome", "Position", opt$tumour_name)
)

write.table(
  data.frame(Chromosome = snp_data$Chromosome, Position = snp_data$Position,
             LogR = snp_data$Tumour_LogR),
  file = logr_file, row.names = FALSE, quote = FALSE, sep = "\t",
  col.names = c("Chromosome", "Position", opt$tumour_name)
)

cat("  BAF output:", baf_file, "\n")
cat("  LogR output:", logr_file, "\n")

# Summary statistics
cat("\n=== Summary ===\n")
cat("Total SNPs:", nrow(snp_data), "\n")
cat("Mean Tumour BAF:", round(mean(snp_data$Tumour_BAF, na.rm = TRUE), 4), "\n")
cat("Mean Tumour LogR:", round(mean(snp_data$Tumour_LogR, na.rm = TRUE), 4), "\n")
cat("Mean Tumour Depth:", round(mean(snp_data$Tumour_Depth, na.rm = TRUE), 1), "\n")
cat("Mean Normal Depth:", round(mean(snp_data$Normal_Depth, na.rm = TRUE), 1), "\n")

cat("\nDone!\n")
