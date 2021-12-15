#--------------------------------------------------------------------------------------
#
#    snRNAseq create differentially expressed gene tables - laptop
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Supplementary tables 1-5

##  Load Packages  --------------------------------------------------------------------
library(cowplot)
library(tidyverse)
library(Seurat)

## Set variables  ---------------------------------------------------------------------
DATA_DIR <- '~/Desktop/single_cell/scRNAseq/batch2_CR5_200121/r_objects/final/'
FIG_DIR <- '~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/'
REGIONS <- c('pfc', 'wge', 'hip',  'tha', 'cer')

for (REGION in REGIONS) {
  
  OBJ <- readRDS(paste0(DATA_DIR, 'seurat.', REGION, '.final.rds'))
  
  assign(paste0('seurat.', REGION), OBJ)
                  
}

for (REGION in c('cer', 'hip', 'pfc', 'tha', 'wge')) {
  
  cat(paste0('\n\n\nGet the ', REGION, 
             ' top DE genes per cluster vs. all other clusters ...\n\n\n'))
  
  OBJ <- get(paste0('seurat.', REGION))
  
  OBJ.markers <- FindAllMarkers(OBJ, only.pos = TRUE, 
                                min.pct = 0.25, logfc.threshold = 0.25) 
  assign(paste0(REGION, '_DE_markers'), OBJ.markers, .GlobalEnv)
  
}
  
  # Split by cluster into a list of dfs
  OBJ_SPLIT <- OBJ %>% 
    group_split(cluster)
  
  # Write top30 DE genes to xlsx file with a single cluster assigned to each tab
  for (i in seq_along(OBJ_SPLIT)) {
    
    xlsx::write.xlsx(x = as.data.frame(OBJ_SPLIT[[i]]), 
                     file = paste0(DATA_DIR, 'pass1/', REGION, 'top30_DEgenes_pass1.xlsx'), 
                     sheetName = paste0('C', unique(OBJ$cluster)[i]),
                     append = TRUE)
  }
  
  # Assign marker files
  assign(paste0(REGION, '_DE_markers'), OBJ.markers, .GlobalEnv)
  assign(paste0(REGION, '_top30_DE_markers'), TOP30, .GlobalEnv)
  
}

for (REGION in c('cer', 'hip', 'pfc', 'tha', 'wge')) {
  
  cat(paste0('\n\n\nGet the ', REGION, 
             ' top DE genes per cluster vs. all other clusters ...\n\n\n'))
  
  OBJ <- get(paste0('seurat.', REGION))
  
  OBJ.markers <- FindAllMarkers(OBJ, only.pos = TRUE, 
                                min.pct = 0.25, logfc.threshold = 0.25) 
  assign(paste0(REGION, '_DE_markers'), OBJ.markers, .GlobalEnv)
  
}

# Save tables
for (INDEX in 1:length(REGIONS)) {
  
  OBJ <- get(paste0(REGIONS[INDEX], '_DE_markers'))
  print(head(OBJ))
  write.table(OBJ, paste0(FIG_DIR, "Supp_table_", INDEX, ".txt"), , sep = '\t', 
              row.names = FALSE, quote = FALSE, col.names = TRUE)

}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
