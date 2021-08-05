#--------------------------------------------------------------------------------------
#
#    snRNAseq analysis - Seurat
#
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------
library(scDblFinder)
library(SingleCellExperiment)
library(scater)
library(Matrix)
library(scran)
library(DropletUtils)
library(cowplot)
library(tidyverse)
library(knitr)
library(rmarkdown)
library(PCAtools)
library(pheatmap)
library(Seurat)
options(stringsAsFactors = FALSE)

## Parse arguments  -------------------------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("Read in scRNAseq QCed sce object and output info  ... \n")
p <- add_argument(p, "region", help = "No brain region specified")
p <- add_argument(p, "in_obj", help = "No snRNAseq QCed SCE object specified")
p <- add_argument(p, "out_seurat", help = "No seurat output RDS object specified")
p <- add_argument(p, "out_topDEgenes", help = "No top DE gene output object specified")
args <- parse_args(p)
print(args)

##  Define variables  -----------------------------------------------------------------
REGION <- args$region
IN_OBJ <- args$in_obj
OUT_SEURAT <- args$gwas
OUT_TOP_DE_GENES

VAR_FEAT  <- 2000       # No. variable features to use for FindVariableFeatures()
PCA_UPPER_DIM <- 17     # No. PCA dimensions to be considered for dim reduction
RESOLUTION <- 0.5       # Variable for cluster granularity

seurat_variables <- data.frame(Variable_Features=VAR_FEAT,
                               Principal_components=PCA_UPPER_DIM,
                               Resolution=RESOLUTION)

##  Load data  ------------------------------------------------------------------------
seurat.filt <- readRDS(IN_OBJ)


##  Convert to Seurat Object  ---------------------------------------------------------
seurat.filt <- as.Seurat(sce.filt, counts = "counts", data = NULL)


##  Normalisation and scaling  --------------------------------------------------------
seurat.filt <- NormalizeData(seurat.filt) %>% 
  
  FindVariableFeatures(selection.method = "vst", nfeatures = VAR_FEAT) 

seurat.filt <- ScaleData(seurat.filt)# doesn't work unless create new variable?

#SCTransform(seurat.filt, vars.to.regress = c("percent.mt", "percent.ribo_RPL", "percent.ribo_RPS"), verbose = FALSE)

# Identify the 10 most highly variable genes
top10_var_genes <- head(VariableFeatures(seurat.filt), 10)
var_feat_plot <- LabelPoints(plot = VariableFeaturePlot(seurat.filt), 
                             points = top10_var_genes, repel = TRUE)

##  Dimensionality reduction  and clustering  -----------------------------------------
seurat.filt <- RunPCA(seurat.filt, features = VariableFeatures(object = seurat.filt)) %>%
  
  FindNeighbors(dims = 1:PCA_UPPER_DIM) %>%
  FindClusters(resolution = RESOLUTION) %>%
  
  RunUMAP(dims = 1:PCA_UPPER_DIM) %>%
  RunTSNE() 

# Plot  #####
elbow.plot <- ElbowPlot(seurat.filt)
pca.plot <- DimPlot(seurat.filt, pt.size = 0.2, reduction = "pca", label = TRUE)
tsne.plot <- DimPlot(seurat.filt, pt.size = 0.2, reduction = "tsne", label = TRUE) 
umap.plot <- DimPlot(seurat.filt, pt.size = 0.2, reduction = "umap", label = TRUE) 

plot_grid(pca.plot, elbow.plot, tsne.plot, umap.plot)

print(seurat.filt[["pca"]], dims = 1:20, nfeatures = 10)
pca_loadings_ <- VizDimLoadings(seurat.filt, dims = 1:2, reduction = "pca", )
pca_heatmaps_1_6 <- DimHeatmap(seurat.filt, dims = 1:6, cells = 500, balanced = TRUE)
pca_heatmaps_7_12 <- DimHeatmap(seurat.filt, dims = 7:12, cells = 500, balanced = TRUE)
pca_heatmaps_13_18 <- DimHeatmap(seurat.filt, dims = 13:18, cells = 500, balanced = TRUE)
pca_heatmaps_19_24 <- DimHeatmap(seurat.filt, dims = 19:24, cells = 500, balanced = TRUE)

## Cell-cycle  -----------------------------------------------------------------------
seurat.filt_cell_cycle <- CellCycleScoring(seurat.filt, 
                                          s.features = cc.genes.updated.2019$s.genes,
                                          g2m.features = cc.genes.updated.2019$g2m.genes,
                                          set.ident = TRUE)


seurat.filt$S.Score <- seurat.filt_cell_cycle$S.Score
seurat.filt$G2M.Score <- seurat.filt_cell_cycle$G2M.Score
seurat.filt$Phase <- seurat.filt_cell_cycle$Phase

# Clustering QC plots
cell_cycle_plot <- DimPlot(seurat.filt, pt.size = 0.2, reduction = "umap", 
                           split.by = "Phase")
cluster_batch_plot <- DimPlot(seurat.filt, pt.size = 0.2, reduction = "umap", 
                              split.by = "batch", label = TRUE)
cluster_donor_plot <- DimPlot(seurat.filt, pt.size = 0.2, reduction = "umap", 
                              split.by = "donor", label = TRUE)

## Cluster annotation - Differential expression  --------------------------------------
# Get the top DE genes per cluster vs. all other clusters - 5-10min 
cat("Get the top DE genes per cluster vs. all other clusters ...\n")
seurat.filt.markers <- FindAllMarkers(seurat.filt, only.pos = TRUE, 
                                     min.pct = 0.25, logfc.threshold = 0.25) 

# Get top 30 genes per cluster
top30 <- seurat.filt.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 30, wt = avg_log2FC)


##  Save object  ----------------------------------------------------------------------
# Need to do this to carry over regional cluster IDs
seurat.filt <- StashIdent(seurat.filt, save.name = paste0(tolower(REGION),'_idents')
saveRDS(seurat.filt, OUT_SEURAT)
saveRDS(seurat.filt.markers, OUT_TOP_DE_GENES)


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

