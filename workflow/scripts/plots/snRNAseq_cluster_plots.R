#--------------------------------------------------------------------------------------
#
#    snRNAseq cluster plots - laptop
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures 1, S1-5 - laptop

#  Nowakowski and Polioudakis for FC 
#  Nowakowski for GE
#  Aldinger for Cer

##  Load Packages  --------------------------------------------------------------------
library(cowplot)
library(tidyverse)
library(Seurat)
library(Matrix) # to load Polioudakis dataset
library(clustifyr)


## Set variables  ---------------------------------------------------------------------
DATA_DIR <- '~/Desktop/single_cell/scRNAseq/batch2_CR5_200121/r_objects/final/'
FIG_DIR <- '~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/'
ALDINGER_DIR <- '~/Desktop/single_cell/public_datasets/aldinger2021_fetal_cer_prePreint/'
REGIONS <- c('cer', 'hip', 'pfc', 'tha', 'wge')
PCA_UPPER_DIM <- 17

## Load Data --------------------------------------------------------------------------
# Seurat objects  ----
for (REGION in REGIONS) { 
  
  seurat.obj <- readRDS(paste0(DATA_DIR, 'seurat.', REGION, '.final.rds'))
  
  # Add plot IDs - remove region labels
  seurat.obj$plotIDs <- gsub("^.*?-", "", seurat.obj$cellIDs)
  
  assign(paste0('seurat.', REGION), seurat.obj, .GlobalEnv)
  
  cluster_plot <- DimPlot(seurat.obj, reduction = "umap", group.by = 'plotIDs',
                          label = TRUE, label.size = 4,
                          pt.size = 0.1, repel = TRUE) + ggtitle(NULL) 
  
  cluster_plot_noLabel <- DimPlot(seurat.obj, reduction = "umap", group.by = 'plotIDs',
                                  pt.size = 0.1) + ggtitle(NULL) 
  
  assign(paste0(REGION, '_plot'), cluster_plot, .GlobalEnv)
  assign(paste0(REGION, '_noLabel_plot'), cluster_plot_noLabel, .GlobalEnv)
  
}

# Standardise colours
ExN_blue <- c('#9A6324', '#76B5C5', "#00BDD2", "#00B6EB", '#CEE5FD', 
              '#ABDBE3','#1E81B0', '#779CBA', '#B8D2EB')
InN_green <- c('#3CBB75FF', '#55C667FF', '#00FF00A5','#73D055FF', '#95D840FF', '#10A53DFF', '#006400')
RG_red <- c('#FAA0A0', '#FF5959', '#F75151', '#EF0029', '#D2042D')
cycPro_lavender <- c("#DCBEFF")
endo_brown <- '#9A6324'
ip_magenta <- "#D078FF"
OPC_yelow <- '#FFDE725FF'
MG_orange <- '#F58231'

# GE
ge_colours <- c('#DCBEFF', '#006400', '#55C667FF', '#00FF00A5', '#10A53DFF',
                 '#95D840FF', '#73D055FF',  '#3CBB75FF', '#FAA0A0', '#EF0029', 
                 '#D2042D')

ge_final_plot <- wge_plot +
  ggtitle('Ganglionic Eminence') +
  theme_bw() +
  NoLegend() +
  theme(legend.text=element_text(size = 12),
        legend.title = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 12)) +
  scale_color_manual(values = ge_colours)


# Cer
cer_colours <- c('#DCBEFF', '#9A6324', '#76B5C5', '#CEE5FD', '#00BDD2',  
                 '#00B6EB', '#ABDBE3', '#1E81B0', '#3CBB75FF', '#00FF00A5', 
                 '#006400', '#95D840FF', '#B7FFB7', '#10A53DFF', '#F58231', 
                 '#949494', '#CCCCCC', '#FDE725FF', '#FAA0A0','#EF0029', 
                 '#D2042D')

cer_final_plot <- cer_plot +
  ggtitle('Cerebellum') +
  theme_bw() +
  NoLegend() +
  theme(legend.text=element_text(size = 12),
        legend.title = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 12)) +
  scale_color_manual(values = cer_colours)

# FC
fc_colours <- c("#DCBEFF", '#9A6324', '#CEE5FD', '#76B5C5', '#ABDBE3', '#779CBA', 
                '#1E81B0', '#95D840FF', '#55C667FF', '#73D055FF', '#10A53DFF', 
                '#D078FF', '#F58231', '#CCCCCC', '#FDE725FF', '#FF5959', '#EF0029')

fc_final_plot <- pfc_plot +
  ggtitle('Frontal Cortex') +
  theme_bw() +
  NoLegend() +
  theme(legend.text=element_text(size = 12),
        legend.title = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 12)) +
  scale_color_manual(values = fc_colours)

# Hip
hip_colours <- c("#DCBEFF", '#9A6324', '#76B5C5', "#00BDD2", '#CEE5FD',  
                "#00B6EB", '#ABDBE3','#1E81B0', '#779CBA', '#73D055FF', 
                '#10A53DFF', '#95D840FF', '#F58231', "#949494", '#CCCCCC',  
                '#FDE725FF', '#FAA0A0', '#EF0029', '#D2042D')

hip_final_plot <- hip_plot +
  ggtitle('Hippocampus') +
  theme_bw() +
  NoLegend() +
  theme(legend.text=element_text(size = 12),
        legend.title = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 12)) +
  scale_color_manual(values = hip_colours)



# Tha
thalamus_colours <- c("#DCBEFF", '#9A6324', '#76B5C5', "#00BDD2", "#00B6EB", 
                      '#CEE5FD', '#ABDBE3','#1E81B0', '#779CBA', '#B8D2EB', 
                      '#CEE5FD', "#3CBB75FF", "#55C667FF", "#73D055FF",
                      '#95D840FF', "#D078FF", '#F58231', '#FDE725FF', '#FAA0A0',
                      '#FF5959', '#F75151', '#EF0029', '#D2042D')

tha_final_plot <- tha_plot +
  ggtitle('Thalamus') +
  theme_bw() +
  NoLegend() +
  theme(legend.text=element_text(size = 12),
        legend.title = element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 12)) +
  scale_color_manual(values = thalamus_colours)
  


# Plot
fig1_plot <- plot_grid(fc_final_plot, 
                       ge_final_plot,
                       hip_final_plot,
                       tha_final_plot,
                       cer_final_plot, 
                       labels = 'AUTO', label_size = 16)

# Save plot
tiff(paste0(FIG_DIR, "Fig_1.tiff"), height = 20, width = 30, units='cm',
     compression = "lzw", res = 300)
fig1_plot
dev.off()

##  Label Transfer  - Polioudakis  - FC only  -----------------------------------------
load("~/Desktop/single_cell/public_datasets/sc_dev_cortex_geschwind/raw_counts_mat.rdata")
gersch_matrix <-  as.matrix(raw_counts_mat)
gersch_meta <-  as.data.frame(read_csv("Desktop/single_cell/public_datasets/sc_dev_cortex_geschwind/cell_metadata.csv")) 

#gersch_meta$Cell <- paste0(gersch_meta$Cell, "-1")
gersch_meta <- column_to_rownames(gersch_meta, var = "Cell")

rna_gersch <- CreateSeuratObject(counts = gersch_matrix, meta.data = gersch_meta)

# Prepare gersch data for label transfer
rna_gersch <- NormalizeData(rna_gersch) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData()

# Set anchors based on main clusters
study.anchors <- FindTransferAnchors(reference = rna_gersch, query = seurat.pfc, 
                                     dims = 1:PCA_UPPER_DIM)
predictions <- TransferData(anchorset = study.anchors, refdata = rna_gersch$Cluster, 
                            dims = 1:PCA_UPPER_DIM)
seurat.pfc <- AddMetaData(seurat.pfc, metadata = predictions)
seurat.pfc$prediction.match <- seurat.pfc$predicted.id == seurat.pfc$seurat_clusters

fig_S1C <- DimPlot(seurat.pfc, pt.size = 0.1, reduction = "umap",
                                      label = TRUE, group.by = "predicted.id") +
  ggtitle(NULL)
  
  

##  ClustifyR - using reference dataset - FC and GE  ----------------------------------
load(file = paste0("Desktop/single_cell/public_datasets/Nowakowski/Nowakowski2018_cortex_dev_clustifyR.rda"))

## Run clustifyr  
for (REGION in c("pfc", "wge")) {
  
  seurat.obj <- get(paste0('seurat.', REGION))
  
  # Set varibles genes to consider - 500-1000 recommended 
  var_genes <- VariableFeatures(seurat.obj)[1:1000]

  res <- clustifyr::clustify(
    input = seurat.obj,       # a Seurat object
    ref_mat = ref_cortex_dev,    # matrix of RNA-seq expression data for each cell type
    cluster_col = "RNA_snn_res.0.5", # name of column in meta.data containing cell clusters
    obj_out = TRUE,# output Seurat object with cell type inserted as "type" column
    query_genes = var_genes
  )
  
  # Plots
  clustifyR_plot <- DimPlot(res, reduction = "umap", label = TRUE,
                            pt.size = 0.01, group.by = 'type') +
    ggtitle(NULL)
  
  # Assign plots
  if (REGION == 'pfc') {
    
    assign('fig_S1B', clustifyR_plot, envir = .GlobalEnv)
    
  } else {
    
    assign('fig_S2B', clustifyR_plot, envir = .GlobalEnv)
    
  }
  
}


##  ClustifyR - using gene set - Cer --------------------------------------------------
## Load and prep data  ----------------------------------------------------------------
aldinger_cer_markers <- read_tsv(paste0(ALDINGER_DIR, 
                                        "aldinger2019_supp_markers_table_S9.txt")) 

# matrixize_markers needs 'gene' and 'cluster' columns
colnames(aldinger_cer_markers) <- c("gene", "P value", "Ave LogFC", "Pct.1", 
                                    "Pct.2", "Adj P value", "Cluster", "cluster")

# Need to remove 10 'genes' which are labelled with a date instead of a gene name
row_remove <- "01-Mar|03-Mar|04-Mar|05-Sep|07-Sep|08-Sep|09-Sep|10-Sep|11-Mar|11-Sep"
aldinger_cer_markers <- aldinger_cer_markers[!grepl(row_remove, 
                                                    aldinger_cer_markers$gene),]

# Create gene list matrix
gene_matrix <- matrixize_markers(
  marker_df = aldinger_cer_markers,
  ranked = FALSE,
  remove_rp = TRUE,
  unique = FALSE
)

# Run clustify - extract seurat object
list_res <- clustify_lists(
  input = seurat.cer,            # matrix of normalized single-cell RNA-seq counts         # meta.data table containing cell clusters
  cluster_col = "cellIDs",    # name of column in meta.data containing cell clusters
  marker = gene_matrix,                   # list of known marker genes
  metric = "pct"                 # test to use for assigning cell types
)

# Plots
# Plot clusters
fig_S5B <- DimPlot(list_res, reduction = "umap", 
                   pt.size = 0.1, group.by = "type.clustify",
                   repel = TRUE) + ggtitle(NULL)


##  Create group plots ----------------------------------------------------------------
# S1
fig_S1 <- plot_grid(pfc_plot, fig_S1B, fig_S1C, pfc_noLabel_plot, 
                    labels = 'AUTO', label_size = 16)

# S2
fig_S2 <- plot_grid(wge_plot, fig_S2B, wge_noLabel_plot,
                    labels = 'AUTO', label_size = 16)

# S3
fig_S3 <- plot_grid(hip_plot, hip_noLabel_plot,
                    labels = 'AUTO', label_size = 16)

# S4
fig_S4 <- plot_grid(tha_plot, tha_noLabel_plot,
                    labels = 'AUTO', label_size = 16)

# S5
fig_S5 <- plot_grid(cer_plot, fig_S5B, cer_noLabel_plot,
                    labels = 'AUTO', label_size = 16)


# Save plots
# Fig S1 
tiff(paste0(FIG_DIR, "Fig_S1.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
fig_S1
dev.off()

# Fig S2
tiff(paste0(FIG_DIR, "Fig_S2.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
fig_S2
dev.off()

# Fig S3
tiff(paste0(FIG_DIR, "Fig_S3.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
fig_S3
dev.off()

# Fig S4
tiff(paste0(FIG_DIR, "Fig_S4.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
fig_S4
dev.off()

# Fig S5
tiff(paste0(FIG_DIR, "Fig_S5.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
fig_S5
dev.off()

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

