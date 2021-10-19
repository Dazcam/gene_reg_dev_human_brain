#--------------------------------------------------------------------------------------
#
#     Map PGC3 index SNPs and proxies with R2 > 0.8 to snATACseq peaks
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# LDlinkR - requires a personal access token
# Github - https://github.com/CBIIT/LDlinkR
# Vingette: https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html

## Info  ------------------------------------------------------------------------------

# STATUS: Incomplete: BiomaRt choking one large batches even in split batch data 
# corrupted may need to download dsSNP VCF file and cross ref manually

# Extract proxy SNPs of PGC3 index SNPs with r2 > 0.8
# 287 SNPs in PGC3_SZ_index_SNPs.xlsx
# 9 SNPs failed 

    # Not encoded as an rsID in xslx file: 8:4180090_T_A 
    # Not in 1000G reference panel: rs2494638, rs12514660, rs2914025, rs55858161
    # Not a not a biallelic variant: rs62152282, rs35026989, rs11241041, rs650520

# 287 proxy SNP files
    
##  Load Packages  --------------------------------------------------------------------
library(LDlinkR)  
library(readxl)
library(tidyverse)
library(Repitools)

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
OUT_DIR <- "~/Desktop/single_nuclei/snATACseq/"
REGIONS <- c("FC", "GE")
CELL_TYPES <- c("FC-ExN", "FC-InN", "FC-RG", "FC-MG", "FC-N.undef", "GE-InN", "GE-RG")
TOKEN <- ""

##  Read in PGC3 index SNPs one is not an rsID just removed it  ------------------------
snps_0.1 <- read_excel("Desktop/single_nuclei/snRNAseq/tests/PGC3_SZ_index_SNPs_0.1.xlsx") %>%
  dplyr::select(`top-index`) %>%
  filter(!grepl('8:4180090_T_A', `top-index`)) %>%
  base::as.data.frame(snps)
  
# Extract proxy SNPs in LD - can't specify threshold of 0.8 for this so need to go
# With default of 0.01. LD taken from 1000Gs phase 3 - CEU population
snps_in_LD <- LDlinkR::LDproxy_batch(snp = snps, token = TOKEN)

## Remove SNPs that failed proxy batch call from snp list -----------------------------
error_snps <- c('rs2494638', 'rs12514660', 'rs2914025', 'rs55858161',
                'rs62152282', 'rs35026989', 'rs11241041', 'rs650520')
snps_no_errors <- setdiff(snps, error_snps)

##  Get proxy SNPs with R2 > 0.8 for all index SNPs. ----------------------------------
for (SNP in snps$`top-index`) {
  
  skip_to_next <- FALSE
  cat(paste0('\nLoading proxy SNPs for ', SNP, '... \n'))
  
  tryCatch(
    
    
    proxies <<- read_delim(paste0("~/Desktop/single_nuclei/snRNAseq/tests/PGC3_proxy_SNPs/", 
                                  SNP, ".txt"), delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE),
    
    error = function(e) {
      
      cat(paste0("No file for: ", SNP, "\n"))
      skip_to_next <<- TRUE })
  
  if (skip_to_next) { next }     
  
  proxies_r2_0.8 <- proxies %>%
    filter(R2 > 0.8)
  
  # Generate stats
  index_in_proxies <- SNP %in% proxies$Coord
  index_in_proxies_r2_0.8 <- SNP %in% proxies_r2_0.8$Coord
  snps_in_proxies <- nrow(proxies)
  snps_in_proxies_r2_0.8 <- nrow(proxies_r2_0.8)
  cat(paste0('\n', nrow(proxies), ' snps in ', SNP, ' proxies table ... \n'))
  cat(paste0('\nChecking if ', SNP, ' index SNP is in proxies table: ', index_in_proxies, '\n'))
  cat(paste0('\nChecking if ', SNP, ' index SNP is in proxies r2 > 0.8 table: ', index_in_proxies_r2_0.8, '\n'))
  
  if (exists("snp_proxies_df")) {
    
    snp_proxies_df <- rbind(snp_proxies_df, proxies_r2_0.8) 
    snp_stats <- cbind(SNP, index_in_proxies, snps_in_proxies, snps_in_proxies_r2_0.8)
    snp_proxies_summary_df <- rbind(snp_proxies_summary_df, snp_stats)
      
     } else {
    
    snp_proxies_df <- proxies_r2_0.8
    snp_proxies_summary_df <- cbind(SNP, index_in_proxies, snps_in_proxies, snps_in_proxies_r2_0.8)
    
  }
  
  
}

# Retain only unique SNPs in list - reduces from 237K to 1.7K SNPs!
cat('\nRetain only unique SNPs  ... \n')
snp_proxies_unique_df <- snp_proxies_df %>% distinct(Coord, .keep_all = TRUE)
scz_snps <- as.vector(snp_proxies_unique_df$Coord)
scz_snps

# Remove weird IDs: "ALU_umary_ALU_20" etc.
scz_snps <- scz_snps[grepl('rs', scz_snps)]

cat(paste0(length(scz_snps), ' SNPs retained.\n'))

library(httr)    
set_config(use_proxy(url="10.3.100.207",port=8080))
SNPs <- ncbi_snp_query(scz_snps)

## Get hg38 base postions for rsIDs using biomaRt -------------------------------------
cat('\nUsing BiomaRt to get hg38 base postions for SNP rsIDs  ... \n')
ensembl <- useEnsembl("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
# Need to batch the query: https://support.bioconductor.org/p/23684/
SNPs_a <- getBM(attributes=c("refsnp_id",
                           "chr_name",
                           "chrom_start",
                           "chrom_end"),
              filters ="snp_filter", 
              values = scz_snps[1:50000], 
              mart = ensembl, 
              uniqueRows=TRUE)
Sys.sleep(1)
SNPs_b <- getBM(attributes=c("refsnp_id",
                             "chr_name",
                             "chrom_start",
                             "chrom_end"),
                filters ="snp_filter", 
                values = scz_snps[50001:100000], 
                mart = ensembl, 
                uniqueRows=TRUE)
Sys.sleep(1)
SNPs_c <- getBM(attributes=c("refsnp_id",
                             "chr_name",
                             "chrom_start",
                             "chrom_end"),
                filters ="snp_filter", 
                values = scz_snps[100001:150000], 
                mart = ensembl, 
                uniqueRows=TRUE)
Sys.sleep(1)
SNPs_d <- getBM(attributes=c("refsnp_id",
                             "chr_name",
                             "chrom_start",
                             "chrom_end"),
                filters ="snp_filter", 
                values = scz_snps[150001:length(scz_snps)], 
                mart = ensembl, 
                uniqueRows=TRUE)

SNPs <- rbind(SNPs_a, SNPs_b, SNPs_c, SNPs_d)

cat(paste0(nrow(SNPs), ' SNPs retained.\n'))

# Some SNPs include duplicates due to chr patches - I removed these
cat('\nRemoving SNPs on CHR patches ... \n')
snps_no_patches <- SNPs %>%
  filter(!grepl('_', chr_name)) %>%
  dplyr::select(-chrom_end) %>%
  dplyr::rename('snpID' = refsnp_id, 'hg38_base_position' = chrom_start)
cat(paste0(nrow(snps_no_patches), ' SNPs retained. \n'))

## Check for overlap of SNPs in snATACseq peaks of all cell types  --------------------
for (CELL_TYPE in CELL_TYPES) {
  
  REGION <- word(CELL_TYPE, 1, sep = "-")
  CELL_ID <- word(CELL_TYPE, 2, sep = "-")

  cat(paste0('\nLoading peaks for ', CELL_TYPE, ' ... \n'))
  peaks <- readRDS(paste0(OUT_DIR, 'PeakCalls/', REGION, '/', CELL_ID,'-reproduciblePeaks.gr.rds'))
  peaks_df <- Repitools::annoGR2DF(peaks)
  
  cell_overlaps <- data.frame()
  
  cat(paste0('\nChecking for SNP overlaps in ', CELL_TYPE, ' ... \n'))
  for (i in 1:nrow(snps_no_patches)) {
    
    BASE_POSITION <- snps_no_patches$hg38_base_position[i]
    CHR <- snps_no_patches$chr_name[i]
    cat(paste0('\nSNP: ', 
               snps_no_patches$snpID[i], ', position: ',
               BASE_POSITION, '... \n'))
    
    overlaps <- filter(peaks_df, start <= BASE_POSITION, end >= BASE_POSITION, chr == paste0('chr', CHR))
    
    if (nrow(overlaps) > 0) { 
      
      overlaps <- cbind(overlaps, snps_no_patches[i,])
      print(overlaps) 
      
    }
    
    cell_overlaps <- rbind(cell_overlaps, overlaps)
    
  } 
  
  cat(paste0('\nAll ', nrow(snps_no_patches), ' SNPs checked in ', CELL_TYPE, ' ... \n'))
  
  cat(paste0('\nWriting overlapping SNPs to file ... \n'))
  cell_overlaps <- cell_overlaps[, c(1:3, 16:18, 4:15)]
  write_tsv(cell_overlaps, paste0(OUT_DIR, 'PeakCalls/', REGION, '/', REGION, '-', CELL_ID,'_PGC3_r2_0.8_SNP_overlaps.tsv'))
  assign(paste0(REGION, '_', CELL_ID, '_SNP_overlaps'), cell_overlaps)
  
}

cat('Done.\n')
  
  ## Need to cross reference this with scripts/R/snATACseq_map_PGC3_snps_to_peaks.R on cluster
  
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

