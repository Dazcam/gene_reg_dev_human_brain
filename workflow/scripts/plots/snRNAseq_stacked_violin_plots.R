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
                 "MKI67", "C3", "ITM2A", "LHX6", "SST", 
                 "CALB2", "SCGN", "TLE3", "CUX2", "SOX5",
                 "FEZF2", "CRYM", "LHX2", "UNC5D", "GPR6", "MEF2C", "KCNIP2")
ge_features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                 "MKI67", "C3", "ITM2A", "LHX6", "SIX3", 
                 "PROX1", "TSHZ1", "DLX1", "SCGN", "ISL1")
hip_features <- c("NEUROD1", "GRIK4", "EOMES", "GLI3", "OLIG1", 
                  "MKI67", "C3", "ITM2A", "SLC17A6", "NEUROD2",
                  "GAD2", "MEF2C", "TNC", "PROX1", "RELN", 
                  "PID1", "LHX6")
tha_features <- c("GAD1", "EOMES", "GLI3", "OLIG1", "MKI67", "C3", 
                  "ITM2A", "SLC1A6", "LHX9", "TNC", "GAD2", 
                  "PAX6", "SLC17A6", "SCGN")
cer_features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                  "MKI67", "C3", "ITM2A", "SLC17A6", "CALB1", 
                  "NEUROD1", "SCGN", "CA8", "ITPR1", "RBFOX3", 
                  "RELN", "PAX6", "MEIS2", "PTPRK", "SORCS3")

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

#haem <- c("HBB", "HBD", "HBG1", "HBG2", "HBE1", "HBZ", "HBM", "HBA2", 
#          "HBA1", "HBQ1")

#pouliodakis_genes1 <- c('RGS5', 'CCL3', 'AIF1', 'ITM2A', 'CLDN5', 
#                        'ESAM', 'SOX5', 'DLX1', 'DLX2', 'LHX6', 
#                        'DLX5', 'STMN2', 'BCL11B', 'FEZF2', 'CRYM', 'TFAP2D',
#                          'SEMA3C')
#pouliodakis_genes2 <- c('PPP1R17', 'SSTR2', 'EOMES', 'PENK', 'STMN2',
#                        'POU3F2', 'HMGB2', 'SOX2', 'MKI67', 'PCNA',
#                        'VIM', 'PTPRZ1', 'SLC1A3', 'HES1', 'HOPX')

#molyneaux_genes_1 <- c('LHX5', 'RGS8', 'TLE1', 'TLE3', 'RORB', 
#                     'CYP39A1', 'LHX2', 'DTX4', 'PLXND1',
#                     'MARCKSL1', 'CNTN6', 'FOXO1', 'ECPN', 'LIX1', 'SYT9', 
#                     'OMA1', 'LDB2', 'CRIM1', 'PCP4', 'RAC3', 'DIAP3')

#molyneaux_genes_2 <- c('CTIP2', 'OTX1', 'IGFBP4', 'NR4A3', 'LXN', 
#                       'FOXP2', 'ID2', 'SLITRK1', 'TBR1', "LMO3", 'NURR1')

# Plot
FC_plot <- VlnPlot(seurat.pfc, molyneaux_genes_2, stack = TRUE, flip = TRUE, 
                   cols = fc_colours, same.y.lims = TRUE, fill.by = 'ident') +
  theme(legend.position = "none") + ggtitle("Frontal Cortex")
GE_plot <- VlnPlot(seurat.wge, ge_features, stack = TRUE, flip = TRUE, 
                   cols = ge_colours, same.y.lims = TRUE, fill.by = 'ident') +
  theme(legend.position = "none") + ggtitle("Ganglionic Eminence")
Hip_plot <- VlnPlot(seurat.hip, hip_features, stack = TRUE, flip = TRUE, 
                    cols = hip_colours, same.y.lims = TRUE, fill.by = 'ident') +
  theme(legend.position = "none") + ggtitle("Hippocampus")
Tha_plot <- VlnPlot(seurat.tha, tha_features, stack = TRUE, flip = TRUE, 
                    cols = tha_colours, same.y.lims = TRUE, fill.by = 'ident') +
  theme(legend.position = "none") + ggtitle("Thalamus")
Cer_plot <- VlnPlot(seurat.cer, cer_features, stack = TRUE, flip = TRUE, 
                    cols = cer_colours, same.y.lims = TRUE, fill.by = 'ident') +
  theme(legend.position = "none") + ggtitle("Cerebellum")


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
