#--------------------------------------------------------------------------------------
#
#    ArchR - Unconstrained integration
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# ArchR manual - https://www.archrproject.com/index.html
# ArchR GitHiub - https://github.com/GreenleafLab/ArchR
# Summarized Expriment - https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
# Harmony - Github -   https://github.com/immunogenomics/harmony
# Granges - https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html

## Info  ------------------------------------------------------------------------------

#  snATAC-seq - run unconstrained and contrained integration

#   Here we map the RNA-seq cluster IDs to our ATAC-seq clusters. It is a 2-step process.
#   First unconstrained integration is run which is a broad pass at mapping the RNA-seq 
#   cluster IDs to the ATAC-seq clusters. Then the constrained integration is run, this 
#   time by adding supervised cell groupings as IDed in the unconstrained analysis 
#   i.e. InNs (InN-1, InN-2) and ExNs (ExN-1, ExN-2) etc. are clumped. This enables more 
#   accurate cell mapping between the modalites.

##  Load Packages  --------------------------------------------------------------------
library(ArchR)
library(pheatmap)
library(tidyverse)
library(rmarkdown)
library(BSgenome.Hsapiens.UCSC.hg38) 
library(ComplexHeatmap)
library(clustree)
library(cowplot)
library(argparser)

## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("\nRead brain region and output directory for snATACseq QC ... \n")
p <- add_argument(p, "region", help = "No brain region specified")
p <- add_argument(p, "data_dir", help = "No input data directory specified")
p <- add_argument(p, "archR_out_dir", help = "No ArchR output directory specified")
p <- add_argument(p, "markdown_file", help = "No markdown file path specified")
p <- add_argument(p, "report_dir", help = "No report output directory specified")
p <- add_argument(p, "report_file", help = "No report filename specified")
args <- parse_args(p)
print(args)


##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
REGION <- args$region
DATA_DIR <- args$data_dir
OUT_DIR <- args$archR_out_dir
MARKDOWN_FILE <- args$markdown_file
REPORT_DIR <- args$report_dir
REPORT_FILE <- args$report_file

addArchRThreads(threads = 8) # Set Hawk to 32 cores so 0.75 of total
addArchRGenome("hg38")


##  Load ArchR project  -------------------------------------------------------------------
cat(paste0('\nLoading ArchR project for ', REGION, ' ... \n'))
archR <- loadArchRProject(path = OUT_DIR)


##  Integrating snATACseq with snRNAseq - Cptr 8  -----------------------------------------
# Load Seurat RNA data
if (REGION == 'Cer') {
  
  cat(paste0('\nLoading Seurat object for ', REGION, ' ... \n'))
  seurat.obj <- readRDS("R_objects/seurat.cer.final.rds")
  seurat.obj$cellIDs <- gsub('Cer-', '', seurat.obj$cellIDs)
  
} else if (REGION == 'FC') {
  
  cat(paste0('\nLoading Seurat object for ', REGION, ' ... \n'))
  seurat.obj <- readRDS("R_objects/seurat.pfc.final.rds")
  seurat.obj$cellIDs <- gsub('FC-', '', seurat.obj$cellIDs)
  
} else {
  
  cat(paste0('\nLoading Seurat object for ', REGION, ' ... \n'))
  seurat.obj <- readRDS("R_objects/seurat.wge.final.rds")
  seurat.obj$cellIDs <- gsub('GE-', '', seurat.obj$cellIDs)
  
}


#  Run unconstrained integration
cat('\nRunning unconstrained integration ... \n')
archR.2 <- addGeneIntegrationMatrix(
  ArchRProj = archR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony", # Note using Harmony here
  seRNA = seurat.obj,
  addToArrow = FALSE,
  groupRNA = "cellIDs",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

## Unconstrained integration - reporting  ---------------------------------------------
# Confusion matrix - unconstrained cell mappings 
cat('\nCreating tables and plots ... \n')
cM_geneExp <- as.matrix(confusionMatrix(archR.2$Clusters_harmony, archR.2$predictedGroup_Un))
clust_CM_geneExp <- pheatmap::pheatmap(
  mat = as.matrix(cM_geneExp), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f"
)

# Get df of top cellID matches from RNA for each ATAC cluster
preClust <- colnames(cM_geneExp)[apply(cM_geneExp, 1 , which.max)]
integration_df <- t(as.data.frame(cbind(preClust, rownames(cM_geneExp)))) #Assignments
rownames(integration_df) <- c("RNA", "ATAC")
colnames(integration_df) <- NULL 


# Plot RNA and ATAC UMAPs for comparison
clusters_UMAP_har <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                   name = "Clusters_harmony", embedding = "UMAPHarmony")
umap_seurat_plot <- DimPlot(seurat.obj, pt.size = 0.2, reduction = "umap", 
                            label = TRUE) + NoLegend()
integration_UMAP_plot <- plot_grid(umap_seurat_plot, clusters_UMAP_har)


# Prepare cell groupings for constrained integration
# Only cell-types in preClust need to be included 
cM_unconstrained <- as.matrix(confusionMatrix(archR.2$Clusters_harmony, archR.2$predictedGroup_Un))
preClust <- colnames(cM_unconstrained)[apply(cM_unconstrained, 1 , which.max)]
cM_unconstrained2 <- cbind(preClust, rownames(cM_unconstrained))
unique(unique(archR.2$predictedGroup_Un))


## Save ArchR project  ----------------------------------------------------------------
cat('\nSaving project ... \n')
saveArchRProject(ArchRProj = archR.2, 
                 outputDirectory = OUT_DIR, 
                 load = FALSE)


## Create markdown doc  ---------------------------------------------------------------
cat('\nCreating markdown report ... \n')
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)

cat('\nDONE.\n')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
