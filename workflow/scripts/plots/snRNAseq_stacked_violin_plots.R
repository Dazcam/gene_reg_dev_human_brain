#--------------------------------------------------------------------------------------
#
#    snRNAseq stacked violin plots - laptop
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Figures for stacked violin plots (may be better done in python scanpy?)
#  Extended data figures 1-5

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
  assign(paste0('seurat.', REGION), seurat.obj, .GlobalEnv)
  
  
}

# Set features
fc_features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                 "MKI67", "C3", "ITM2A", "LHX6", "SST", "CALB2")
ge_features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                 "MKI67", "C3", "ITM2A", "LHX6", "MEIS2", 
                 "PROX1", "FOXP2")
hip_features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                  "MKI67", "C3", "ITM2A", "SLC17A6", "NEUROD2",
                  "GAD2", "SLC32A1", "MEF2C", "TNC", "SLC1A3",
                  "PROX1")
tha_features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                  "MKI67", "C3", "ITM2A", "SLC17A6", "LHX9", 
                  "TNC", "GAD2", "SLC32A1", "NEUROD2")
cer_features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                  "MKI67", "C3", "ITM2A", "SLC17A6", "CALB1", "NEUROD1")

# Plots
FC_plot <- VlnPlot(seurat.pfc, fc_features, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Frontal Cortex")
GE_plot <- VlnPlot(seurat.wge, ge_features, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Ganglionic Eminence")
Cer_plot <- VlnPlot(seurat.cer, cer_features, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Cerebellum")
Hip_plot <- VlnPlot(seurat.hip, hip_features, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Hippocampus")
Tha_plot <- VlnPlot(seurat.tha, tha_features, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Thalamus")

# Save
tiff(paste0(FIG_DIR, "Ext_Data_Figure_1.tiff"), height = 30, width = 40, units='cm', 
     compression = "lzw", res = 300)
FC_plot 
dev.off()

tiff(paste0(FIG_DIR, "Ext_Data_Figure_2.tiff"), height = 30, width = 40, units='cm', 
     compression = "lzw", res = 300)
GE_plot 
dev.off()

tiff(paste0(FIG_DIR, "Ext_Data_Figure_3.tiff"), height = 30, width = 40, units='cm', 
     compression = "lzw", res = 300)
Hip_plot 
dev.off()

tiff(paste0(FIG_DIR, "Ext_Data_Figure_4.tiff"), height = 30, width = 40, units='cm', 
     compression = "lzw", res = 300)
Tha_plot 
dev.off()

tiff(paste0(FIG_DIR, "Ext_Data_Figure_5.tiff"), height = 30, width = 40, units='cm', 
     compression = "lzw", res = 300)
Cer_plot 
dev.off()

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
