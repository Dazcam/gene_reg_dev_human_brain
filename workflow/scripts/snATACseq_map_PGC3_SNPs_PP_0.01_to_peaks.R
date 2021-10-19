#--------------------------------------------------------------------------------------
#
#   Map PGC3 index SNPs and SNPs with posterior probabilty > 0.01 to snATACseq peaks
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# 1. Merge rsID and index SNP column and remove duplicates from SNPs PP > 0.01 file
# 2. Get hg38 base positions for ~5K SNPs from biomaRt
# 3. Remove SNPs that map to chr patches
# 4. Check if SNPs intersect snATACseq peaks for each cell type

##  Load Packages  --------------------------------------------------------------------
library(readxl)
library(tidyverse)
library(biomaRt)

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
OUT_DIR <- "~/Desktop/single_nuclei/snATACseq/"
REGIONS <- c("FC", "GE")
CELL_TYPES <- c("FC-ExN", "FC-InN", "FC-RG", "FC-MG", "FC-N.undef", "GE-InN", "GE-RG")

# Read in PGC3 index SNPs one is not an rsID just removed it
snps <- read_excel("Desktop/single_nuclei/snRNAseq/tests/PGC3_SZ_index_SNPs_0.01.xlsx") %>%
  filter(!grepl('8:4180090_T_A', index_snp))

# Merge rsID and index SNP column and remove duplicates
cat('\nRetain only unique SNPs  ... \n')
scz_snps <- union(as.vector(snps$index_snp), as.vector(snps$rsid))


## Get hg38 base postions for rsIDs using biomaRt -------------------------------------
cat('\nUsing BiomaRt to get hg38 base postions for SNP rsIDs  ... \n')
ensembl <- useEnsembl("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
SNPs <- getBM(attributes=c("refsnp_id",
                             "chr_name",
                             "chrom_start",
                             "chrom_end"),
              filters ="snp_filter", 
              values = scz_snps, 
              mart = ensembl, 
              uniqueRows=TRUE)

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
  write_tsv(cell_overlaps, paste0(OUT_DIR, 'PeakCalls/', REGION, '/', REGION, '-', 
                                  CELL_ID,'_PGC3_PP_0.01_snATACseq_peak_overlaps.tsv'))
  assign(paste0(REGION, '_', CELL_ID, '_SNP_overlaps'), cell_overlaps)
  
}


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
