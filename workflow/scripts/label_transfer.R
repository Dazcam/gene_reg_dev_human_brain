#--------------------------------------------------------------------------------------
#
#    snRNAseq analysis QC 
#
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------
library(Matrix)
library(cowplot)
library(tidyverse)
library(PCAtools)
library(pheatmap)
library(Seurat)
library(clustifyr)

options(stringsAsFactors = FALSE)

##  Set variables  --------------------------------------------------------------------
DATA_DIR <- "Desktop/single_cell/scRNAseq/batch2_CR5_200121/"
REGION <- "WGE"

##  Load data  ------------------------------------------------------------------------
seurat.obj <- readRDS(paste0("~/sce.rna.", REGION, ".seurat.rds"))

##  Label Transfer  - Polioudakis  ----------------------------------------------------  
load("~/Desktop/single_cell/public_datasets/sc_dev_cortex_geschwind/raw_counts_mat.rdata")
gersch_matrix <-  as.matrix(raw_counts_mat)
gersch_meta <-  as.data.frame(read_csv("~/Desktop/single_cell/public_datasets/sc_dev_cortex_geschwind/cell_metadata.csv")) 

#gersch_meta$Cell <- paste0(gersch_meta$Cell, "-1")
gersch_meta <- column_to_rownames(gersch_meta, var = "Cell")

rna_gersch <- CreateSeuratObject(counts = gersch_matrix, meta.data = gersch_meta)

# Prepare gersch data for label transfer
rna_gersch <- NormalizeData(rna_gersch) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData()

# Set anchors based on main clusters
study.anchors <- FindTransferAnchors(reference = rna_gersch, query = seurat.obj, 
                                     dims = 1:PCA_UPPER_DIM)
predictions <- TransferData(anchorset = study.anchors, refdata = rna_gersch$Cluster, 
                            dims = 1:PCA_UPPER_DIM)
seurat.obj <- AddMetaData(seurat.obj, metadata = predictions)
seurat.obj$prediction.match <- seurat.obj$predicted.id == seurat.obj$seurat_clusters

polioudakis_main_umap_plot <- DimPlot(seurat.obj, pt.size = 0.3, reduction = "umap",
                                      label = TRUE, group.by = "predicted.id") +
  ggtitle("Polioudakis mappings - main")

##  ClustifyR  - Nowakowski  ----------------------------------------------------------
load(file = paste0("Desktop/single_cell/public_datasets/Nowakowski/Nowakowski2018_cortex_dev_clustifyR.rda"))

# Set varibles genes to consider - 500-1000 recommended 
var_genes <- VariableFeatures(seurat.obj)[1:1000]

## Run clustifyr  ---------------------------------------------------------------------
res <- clustifyr::clustify(
  input = seurat.obj,       # a Seurat object
  ref_mat = ref_cortex_dev,    # matrix of RNA-seq expression data for each cell type
  cluster_col = "RNA_snn_res.0.5", # name of column in meta.data containing cell clusters
  obj_out = TRUE,# output Seurat object with cell type inserted as "type" column
  query_genes = var_genes
)

# Plots
nowakowski_heatmap <- clustifyr::plot_cor_heatmap(cor_mat = res$r)

nowakowski_umap_plot <- DimPlot(res, reduction = "umap", label = TRUE,
                                pt.size = 0.01, group.by = 'type')
nowakowski_tsne_plot <- DimPlot(res, reduction = "tsne", label = TRUE,
                                pt.size = 0.01, group.by = 'type')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
