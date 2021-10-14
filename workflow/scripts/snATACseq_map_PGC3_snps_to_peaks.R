#--------------------------------------------------------------------------------------
#
#    Map PGC3 SNPs to peaks 
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  1. Load PGC3 fine-mapped SNPs (591) PP 0.1 - These are mapped to hg19
#  2. Get hg38 base positions for each SNP using BiomaRt
#  3. Remove duplicate SNPs mapped to scaffolds
#  4. Check for SNPs overlapping peaks in FC and GE cell types

##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(biomaRt)
library(readxl)
library(Repitools)

##  Set variables  --------------------------------------------------------------------
cat('\nDefining variables ... \n')
OUT_DIR <- 'results/ARCHR/'
CELL_TYPES <- c('FC-ExN.1', 'FC-ExN.2', 'FC-ExN.4', 'FC-ExN.5', 'FC-InN.1', 'FC-InN.2', 'FC-InN.3',
                   'FC-MG', 'FC-N.undef', 'FC-RG.1', 'FC-RG.2', 'GE-InN.1', 'GE-InN.3', 'GE-InN.6', 
                   'GE-InN.7', 'GE-RG.1', 'GE-RG.2', 'GE-RG.3')

## Load and prep data  ----------------------------------------------------------------
cat('\nLoading SNPs for fine-mapping  ... \n')
scz_snps <- read_excel('input_files/PGC3_finemapped_SNPs.xlsx')
scz_snps <- as.vector(scz_snps$rsid)
scz_snps

## Get hg38 base postions for rsIDs using biomaRt -------------------------------------
cat('\nUsing BiomaRt to get hg38 base postions for SNP rsIDs   ... \n')
ensembl <- useEnsembl("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
SNPs <- getBM(attributes=c("refsnp_id",
                           "chr_name",
                           "chrom_start",
                           "chrom_end"),
              filters ="snp_filter", 
              values = scz_snps, 
              mart = ensembl, 
              uniqueRows=TRUE)

# listAttributes(mart = ensembl)

# Some SNPs include duplicates due to chr patches - I removed these
snps_no_patches <- SNPs %>%
  filter(!grepl('_', chr_name)) %>%
  dplyr::select(-chrom_end) %>%
  rename(snpID = refsnp_id, 'hg38_base_position' = chrom_start)

## Check for overlap of SNPs in snATACseq peaks of all cell types  --------------------
for (CELL_TYPE in CELL_TYPES) {
  
  REGION <- word(CELL_TYPE, 1, sep = "-")
  CELL_ID <- word(CELL_TYPE, 2, sep = "-")
  
  cat(paste0('\nLoading peaks for ', CELL_TYPE, ' ... \n'))
  peaks <- readRDS(paste0(OUT_DIR, REGION, '/PeakCalls/', CELL_ID,'-reproduciblePeaks.gr.rds'))
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
  
  cat(paste0('\nAll 591 SNPs checked in ', CELL_TYPE, ' ... \n'))
  
  cat(paste0('\nWriting overlapping SNPs to file ... \n'))
  cell_overlaps <- cell_overlaps[, c(1:3, 16:18, 4:15)]
  write_tsv(cell_overlaps, paste0('results/fine_mapping/', REGION, '_', CELL_ID, '_SNP_overlaps.tsv'))
  assign(paste0(REGION, '_', CELL_ID, '_SNP_overlaps'), cell_overlaps)
  
}

cat('Done.\n')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
