# -------------------------------------------------------------------------------------
#
#    snATACseq overlap of OCRs containing fine-mapped SNPs with peaks from
#    peak2gene linkage analysis (FC and GE)
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Checking for peak overlaps overlaps in OCRs containing fine-mapped SNPs and OCRs
#  with a peak to gene linkage in either FC or GE. Issue with duplicates here as 
#  'duplicate' peaks don't have identical boundries so aren't identical. Also had 
#  to deal with a linkage from the same peak linking to multiple genes.
 
##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)
library(readxl)

##  Initialise variables  -------------------------------------------------------------
DATA_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/data_for_figures/"
PLOT_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/"
DROPBOX_DIR <- "Dropbox/BRAY_sc_analysis/files_for_paper/ATAC/archR_hawk_broad_cell_types/"

fine_mapped_peaks <- read_xlsx(paste0(DATA_DIR, "PGE3 SZ finemapped SNPs with PP at least 0.01 in broad cell OCRs.xlsx")) %>%
  rename(start = `hg38 OCR start`,
         stop = `hg38 OCR end`,
         snp_ID = `SNP ID`) %>%
  dplyr::mutate(chr = paste("chr" , chr , sep=""))

unique(fine_mapped_peaks$snp_ID)

FC_peak2gene_peaks <- read_delim(paste0(DROPBOX_DIR, "peak2gene_table_FC.tsv"), 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE) %>%
  rename(start = peak_start,
         stop = peak_end)


GE_peak2gene_peaks <- read_delim(paste0(DROPBOX_DIR, "peak2gene_table_GE.tsv"), 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE) %>%
  rename(start = peak_start,
         stop = peak_end)


# Get peaks into correct format for bedr 
fine_mapped_peaks$chr_start <- str_c(fine_mapped_peaks$chr, ':', fine_mapped_peaks$start)
fine_mapped_peaks$bed <- str_c(fine_mapped_peaks$chr_start, '-', fine_mapped_peaks$stop)
fine_mapped_peaks_vect <- fine_mapped_peaks %>% pull(bed)

FC_peak2gene_peaks$chr_start <- str_c(FC_peak2gene_peaks$chr, ':', FC_peak2gene_peaks$start)
FC_peak2gene_peaks$bed <- str_c(FC_peak2gene_peaks$chr_start, '-', FC_peak2gene_peaks$stop)
FC_peak2gene_peaks_vect <- FC_peak2gene_peaks %>% pull(bed) %>% gsub("1.41e+08", "141000000")
FC_peak2gene_peaks_vect[20639] <- sub("1.41e\\+08", "141000000", FC_peak2gene_peaks_vect[20639])

GE_peak2gene_peaks$chr_start <- str_c(GE_peak2gene_peaks$chr, ':', GE_peak2gene_peaks$start)
GE_peak2gene_peaks$bed <- str_c(GE_peak2gene_peaks$chr_start, '-', GE_peak2gene_peaks$stop)
GE_peak2gene_peaks_vect <- GE_peak2gene_peaks %>% pull(bed) 


is.valid.region(GE_peak2gene_peaks_vect)
fine_mapped_srtd <- bedr.sort.region(fine_mapped_peaks_vect)
FC_p2G_peaks_srtd <- bedr.sort.region(FC_peak2gene_peaks_vect)
GE_p2G_peaks_srtd <- bedr.sort.region(GE_peak2gene_peaks_vect)

for (i in 1:length(fine_mapped_srtd))  {
  
  print(i)
  
  cat("\nChecking overlaps in the following peak: ")
  print(fine_mapped_peaks[i,])
  fc_overlaps <- in.region(FC_p2G_peaks_srtd, fine_mapped_srtd[i])
  ge_overlaps <- in.region(GE_p2G_peaks_srtd, fine_mapped_srtd[i])
  cat(paste0('\nNumber of overlaps fc found: ', sum(fc_overlaps)))
  cat(paste0('\nNumber of overlaps ge found: ', sum(ge_overlaps)))
  
  # Get indices for overlaps
  fine_mapped_idx <- i
  
  if (sum(fc_overlaps) == 0) {
    
    fc_p2G_idx <- 'no_overlap'
  
  
  } else {
    
    fc_p2G_idx <- paste(which(fc_overlaps), collapse=",")
    
  }
  
  if (sum(ge_overlaps) == 0) {
    
    ge_p2G_idx <- 'no_overlap'
    
    
  } else {
    
    ge_p2G_idx <- paste(which(ge_overlaps), collapse=",")
    
  }

  if(exists('overlap_df'))  {
    
    overlap_entry <- base::as.data.frame(cbind(fine_mapped_idx, fc_p2G_idx, ge_p2G_idx))
    overlap_df <- rbind(overlap_df, overlap_entry)
    
  } else {
    
    overlap_df <- base::as.data.frame(cbind(fine_mapped_idx, fc_p2G_idx, ge_p2G_idx))
    
  }
  
}

# Pull out entries with overlapping peaks
overlap_peaks <- overlap_df %>% filter(fc_p2G_idx != "no_overlap" |  ge_p2G_idx != "no_overlap") 

for (i in 1:nrow(overlap_peaks)) {
  
  cat('\n\n|----- Extracting genes for the following entry: -----|\n')
  overlap_peaks_entry <- overlap_peaks[i,]
  print(overlap_peaks_entry)
  
  fine_mapped_idx <- as.integer(overlap_peaks_entry$fine_mapped_idx)
  cat('\nFine mapped index is: ', fine_mapped_idx)
  
  ## Section to extract FC gene IDs
  cat('\n\n|-----  Extracting FC gene IDS  -----|\n')
  if (overlap_peaks_entry$fc_p2G_idx != "no_overlap" ) {
    
    fc_idx <- unlist(strsplit(overlap_peaks_entry$fc_p2G_idx, split=","))
    cat(paste0('\nFC_idx is: ', fc_idx))
    
    if (length(fc_idx) > 1) {
      
      fc_bed_ID <- unique(FC_p2G_peaks_srtd[min(fc_idx):max(fc_idx)]) # Note this can be >1 and still not dealt with extras yet 3 in fc
      cat(paste0('\nFC_bed_ID is: ', fc_bed_ID))
      
      if (length(fc_bed_ID) > 1) {
        
        fc_bed_ID <- fc_bed_ID[1] 
        cat(paste0('\nAs FC_bed_ID > 1 final bed ID for this entry is: ', fc_bed_ID))
        
      }
      
    } else {
      
      fc_bed_ID <- FC_p2G_peaks_srtd[as.integer(fc_idx)] # FC idx can also be > 1 
      cat(paste0('\nAs FC_bed_ID = 1 final bed ID for this entry is: ', fc_bed_ID))
      
    }
    
    fc_overlap_entries <- FC_peak2gene_peaks %>% filter(bed == fc_bed_ID) %>%
      pull(bed)
    fc_gene_overlap_entries <- FC_peak2gene_peaks %>% filter(bed == fc_bed_ID) 
    fc_gene_overlap_entries <- fc_gene_overlap_entries %>% pull(gene_ID)
      
    
    
    cat(paste0('\nFC_overlap_entries is/are: ', fc_overlap_entries))
    fc_overlap_gene_IDs <- paste(fc_gene_overlap_entries, collapse=",")
    cat(paste0('\nFC_overlap_gene IDs is/are: ', fc_overlap_gene_IDs))
    
  } else {
    
    fc_bed_ID <- 'no_overlaps'
    fc_overlap_gene_IDs <- 'no_overlaps'
    
  }
  
  ## Section to extract GE gene IDs
  cat('\n\n|-----  Extracting GE gene IDS  -----|\n\n')
  if (overlap_peaks_entry$ge_p2G_idx != "no_overlap" ) {
    
    ge_idx <- unlist(strsplit(overlap_peaks_entry$ge_p2G_idx, split=","))
    cat(paste0('\nge_idx is: ', ge_idx))
    
    if (length(ge_idx) > 1) {
      
      ge_bed_ID <- unique(GE_p2G_peaks_srtd[min(ge_idx):max(ge_idx)]) # Note this can be >1 and still not dealt with extras yet 3 in ge
      cat(paste0('\nGE_bed_ID is: ', ge_bed_ID))
      
      if (length(ge_bed_ID) > 1) {
        
        ge_bed_ID <- ge_bed_ID[1] 
        cat(paste0('\nAs GE_bed_ID > 1 final bed ID for this entry is: ', ge_bed_ID))
        
      }
      
    } else {
      
      ge_bed_ID <- GE_p2G_peaks_srtd[as.integer(ge_idx)] # GE idx can also be > 1 
      cat(paste0('\nAs GE_bed_ID = 1 final bed ID for this entry is: ', ge_bed_ID))
      
    }
    
    ge_overlap_entries <- GE_peak2gene_peaks %>% filter(bed == ge_bed_ID) %>%
      pull(bed)
    ge_gene_overlap_entries <- GE_peak2gene_peaks %>% filter(bed == ge_bed_ID) 
    ge_gene_overlap_entries <- ge_gene_overlap_entries %>% pull(gene_ID)
    
    
    
    cat(paste0('\nGE_overlap_entries is/are: ', ge_overlap_entries))
    ge_overlap_gene_IDs <- paste(ge_gene_overlap_entries, collapse=",")
    cat(paste0('\nGE_overlap_gene IDs is/are: ', ge_overlap_gene_IDs))
    
  } else {
    
    ge_bed_ID <- 'no_overlaps'
    ge_overlap_gene_IDs <- 'no_overlaps'
    
  }
  
  
  
  
  
  # Extract fine mapped peak entry    
  fine_mapped_bed_ID <- fine_mapped_srtd[fine_mapped_idx]
  cat('\nfine_mapped_bed_ID is: ', fine_mapped_bed_ID)
  fine_mapped_peak <- fine_mapped_peaks %>% filter(bed == fine_mapped_bed_ID)
  cat('\nfine_mapped_peak is: ')
  
  if (exists('final_overlap_df')) {
    
    overlap_final_entry <- base::as.data.frame(cbind(fine_mapped_peak, fc_bed_ID, fc_overlap_gene_IDs, ge_bed_ID, ge_overlap_gene_IDs))
#    cat('\noverlap_final_entry is: ', overlap_final_entry )
    final_overlap_df <- rbind(final_overlap_df, overlap_final_entry)
    
  } else {
    
    final_overlap_df <- base::as.data.frame(cbind(fine_mapped_peak, fc_bed_ID, fc_overlap_gene_IDs, ge_bed_ID, ge_overlap_gene_IDs))
    
  }
 
 
}

# Check for duplicated rows
cat('\nTotal entries in final overlap df: ', nrow(final_overlap_df))
cat('\nDuplicated rows in final overlap df: ', sum(duplicated(final_overlap_df)))

# Remove duplicate rows
final_overlap_df <- final_overlap_df[!duplicated(final_overlap_df), ]
cat('\nTotal entries after duplicate removal in final overlap df: ', nrow(final_overlap_df))

write_tsv(final_overlap_df, paste0(DROPBOX_DIR, 'overlap_fineMapSNP_peaks_with_peak2gene_peaks.tsv'))
  

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
