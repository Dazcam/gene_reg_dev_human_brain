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
                  "MKI67", "C3", "ITM2A", "SLC1A6", "LHX9", 
                  "TNC", "GAD2", "PAX6", "SLC17A6")
cer_features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                  "MKI67", "C3", "ITM2A", "SLC17A6", "CALB1", "NEUROD1")


# Set colours
fc_colours <- c('#DCBEFF', '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                '#CEE5FD', '#CEE5FD', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', 
                '#3CBB75FF', '#D078FF', '#F58231', '#CCCCCC', '#FDE725FF', 
                '#FF5959', '#FF5959')
ge_colours <- c('#DCBEFF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF',
                '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#FF5959', '#FF5959',
                '#FF5959')
cer_colours <- c('#DCBEFF', '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                 '#CEE5FD', '#CEE5FD', '#CEE5FD', '#3CBB75FF', '#3CBB75FF', 
                 '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#F58231', 
                 '#CCCCCC', '#CCCCCC', '#FDE725FF', '#FF5959', '#FF5959',
                 '#FF5959')

tha_colours <- c('#DCBEFF', '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                 '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                 '#CEE5FD', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', 
                 '#3CBB75FF', '#F58231', '#FDE725FF', '#FF5959', '#FF5959',
                 '#FF5959', '#FF5959', '#FF5959')


hip_colours <- c('#DCBEFF', '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                 '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#3CBB75FF', 
                 '#3CBB75FF', '#3CBB75FF', '#F58231', '#CCCCCC', '#CCCCCC', 
                 '#FDE725FF', '#FF5959', '#FF5959', '#FF5959')

# Order Idents
Idents(seurat.pfc) <- factor(x = Idents(seurat.pfc), levels = sort(levels(seurat.pfc)))
Idents(seurat.wge) <- factor(x = Idents(seurat.wge), levels = sort(levels(seurat.wge)))
Idents(seurat.cer) <- factor(x = Idents(seurat.cer), levels = sort(levels(seurat.cer)))
Idents(seurat.tha) <- factor(x = Idents(seurat.tha), levels = sort(levels(seurat.tha)))
Idents(seurat.hip) <- factor(x = Idents(seurat.hip), levels = sort(levels(seurat.hip)))

# Plot
FC_plot <- VlnPlot(seurat.pfc, fc_features, stack = TRUE, flip = TRUE, 
                   cols = fc_colours, same.y.lims = TRUE, fill.by = 'ident') +
  theme(legend.position = "none") + ggtitle("Frontal Cortex")
GE_plot <- VlnPlot(seurat.wge, ge_features, stack = TRUE, flip = TRUE, 
                   cols = ge_colours, same.y.lims = TRUE, fill.by = 'ident') +
  theme(legend.position = "none") + ggtitle("Ganglionic Eminence")
Cer_plot <- VlnPlot(seurat.cer, cer_features, stack = TRUE, flip = TRUE, 
                    cols = cer_colours, same.y.lims = TRUE, fill.by = 'ident') +
  theme(legend.position = "none") + ggtitle("Cerebellum")
Hip_plot <- VlnPlot(seurat.hip, hip_features, stack = TRUE, flip = TRUE, 
                    cols = hip_colours, same.y.lims = TRUE, fill.by = 'ident') +
  theme(legend.position = "none") + ggtitle("Hippocampus")
Tha_plot <- VlnPlot(seurat.tha, tha_features, stack = TRUE, flip = TRUE, 
                    cols = tha_colours, same.y.lims = TRUE, fill.by = 'ident') +
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
