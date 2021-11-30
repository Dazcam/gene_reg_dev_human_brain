#--------------------------------------------------------------------------------------
#
#    snRNAseq stacked violin plots - laptop
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Draft for stacked volin plots may be better done in python scanpy

##  Load Packages  --------------------------------------------------------------------
library(cowplot)
library(tidyverse)
library(Seurat)



## Set variables  ---------------------------------------------------------------------
DATA_DIR <- '~/Desktop/single_cell/scRNAseq/batch2_CR5_200121/r_objects/final/'
FIG_DIR <- '~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/'
REGIONS <- c('cer', 'hip', 'pfc', 'tha', 'wge')

## Load Data --------------------------------------------------------------------------
# Seurat objects  ----
for (REGION in REGIONS) { 
  
  seurat.obj <- readRDS(paste0(DATA_DIR, 'seurat.', REGION, '.final.rds'))
  features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                "MKI67", "C3", "ITM2A")
  seurat.plot <- VlnPlot(seurat.obj, features, stack = TRUE, sort = TRUE, flip = TRUE) +
    theme(legend.position = "none") + ggtitle(REGION)
  
  assign(paste0('seurat.', REGION), seurat.obj, .GlobalEnv)
  assign(paste0('seurat.', REGION, '_plot'), seurat.plot, .GlobalEnv)

}







features <- c("EOMES", "MEF2C", "RBFOX3", "PTN", "SLC17A6", 
              "SLC17A7", "ITM2A", "GAD1", "GAD2", "C3", "OLIG1", 
              "MKI67", "SLC1A3", "PAX3", "TNC", "GLI3")

features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
              "MKI67", "C3", "ITM2A")

GAD1
SLC17A7
EOMES
GLI3
OLIG1
MKI67
C3
ITM2A

a <- VlnPlot(seurat.pfc, features, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("FC")
