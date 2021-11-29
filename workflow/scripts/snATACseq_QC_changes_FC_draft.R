#--------------------------------------------------------------------------------------
#
#    ArchR - Bray 2021 scATACseq Frontal Cortex - altered fragments QC - 221121
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# ArchR manual - https://www.archrproject.com/index.html
# ArchR GitHiub - https://github.com/GreenleafLab/ArchR
# Summarized Expriment - https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
# Harmony - Github -   https://github.com/immunogenomics/harmony
# Granges - https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html

## Info  ------------------------------------------------------------------------------

#  Run the analysis up until cluster QC cell removal 

##  Load Packages  --------------------------------------------------------------------
library(ArchR)
library(pheatmap)
library(tidyverse)
library(rmarkdown)
library(BSgenome.Hsapiens.UCSC.hg38) 
library(ComplexHeatmap)
library(clustree)
library(cowplot)

##  Set variables  --------------------------------------------------------------------
WK_DIR <- "~/Desktop/single_cell/scATACseq/"

SAMPLES <- c("510_PFC_ATAC", "611_PFC_ATAC", "993_PFC_ATAC")
addArchRThreads(threads = 2)
addArchRGenome("hg38")
setwd(WK_DIR) # Required or saves all files to ~/

##  Load data - Cptr 1.5  -------------------------------------------------------------
ArrowFiles <- createArrowFiles(
  inputFiles = c(paste0(WK_DIR, SAMPLES[1], "/fragments.tsv.gz"),
                 paste0(WK_DIR, SAMPLES[2], "/fragments.tsv.gz"),
                 paste0(WK_DIR, SAMPLES[3], "/fragments.tsv.gz")),
  sampleNames = SAMPLES,
  minTSS = 4, # Dont set this too high because you can always increase later
  minFrags = 3160, # approx 10^3.5
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

##  Doublets  - Cptr 2  ---------------------------------------------------------------
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", # Refers to the embedding to use for nearest neighbor search 
  # with doublet projection.
  LSIMethod = 1
)

##  Create Arrow project  - Cptr 3  ---------------------------------------------------
fc.archR <- ArchRProject(ArrowFiles = ArrowFiles, 
                         outputDirectory = paste0(WK_DIR, "ArchR/fc.minFrag_3160"),
                         copyArrows = TRUE # This is recommened so that if you modify 
                         # the Arrow files you have an original copy for later usage.
)

##  Save and load Arrow project  - Cptr 3.5  ------------------------------------------
saveArchRProject(ArchRProj = fc.archR, 
                 outputDirectory = paste0(WK_DIR, "ArchR/fc.minFrag_3160"), 
                 load = FALSE)

# Load project
#fc.archR <- loadArchRProject(path = "~/Desktop/single_cell/scATACseq/ArchR/fc")

##  Add coldata  ----------------------------------------------------------------------
fc.archR$donor <- gsub("_PFC_ATAC","", fc.archR$Sample)

##  Inital ArchR QC -------------------------------------------------------------------
# ArchR does some QC when loading the files in so need to load the pre-QC info
# Pre-filter
for (SAMPLE in 1:length(SAMPLES)) {
  
  # Subset IDs
  sampleID <- substr(SAMPLES[SAMPLE], 1, 7)
  donorID <- substr(SAMPLES[SAMPLE], 1, 3)
  
  # Load Pre-filtered data
  preQC_df <- readRDS(paste0(WK_DIR, "QualityControl/", SAMPLES[SAMPLE], "/", 
                             SAMPLES[SAMPLE], "-Pre-Filter-Metadata.rds"))
  preQC_df$log10nFrags <- log10(preQC_df$nFrags)
  
  # TSS Plot
  preQC_tss_uFrag_plot <- ggPoint(
    x = preQC_df[,"log10nFrags"], 
    y = preQC_df[,"TSSEnrichment"], 
    title = SAMPLES[SAMPLE],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(preQC_df[,"log10nFrags"], probs = 0.99)),
    ylim = c(0, quantile(preQC_df[,"TSSEnrichment"], probs = 0.99))
  ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
  
  # Assign frag plots and counts_df
  assign(paste0("preQC_tss_uFrag_plot_", donorID), preQC_tss_uFrag_plot)
  assign(paste0("counts_df_", donorID), 
         data.frame("Sample" = sampleID,
                    "Cells_Pass_Filter" = sum(preQC_df$Keep),
                    "Cells_dropped" = sum(preQC_df$Keep == 0),
                    "Total_Frags" = sum(preQC_df$nFrags),
                    "Median_Frags" = median(preQC_df$nFrags[preQC_df$Keep==1]),
                    "Median_TSS_Enrichment" = median(preQC_df$TSSEnrichment[preQC_df$Keep==1])))
  
}

## Initial QC reporting  --------------------------------------------------------------
# Pre-filter tss-frag plot
preQC_tss_uFrag_plot <- ggAlignPlots(preQC_tss_uFrag_plot_510, preQC_tss_uFrag_plot_611, 
                                     preQC_tss_uFrag_plot_993, type = "h")
# Counts df
counts_df <- rbind(counts_df_510, counts_df_611, counts_df_993)

## PostQC
fc.archR.meta <- as.data.frame(getCellColData(fc.archR))
fc.archR.meta$log10nFrags <- log10(fc.archR.meta$nFrags)

tss_uFrag_plot <- ggPoint(
  x = fc.archR.meta[,"log10nFrags"], 
  y = fc.archR.meta[,"TSSEnrichment"], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(fc.archR.meta[,"log10nFrags"], probs = 0.99)),
  ylim = c(0, quantile(fc.archR.meta[,"TSSEnrichment"], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3.5, lty = "dashed")


fragSize_plot <- plotFragmentSizes(ArchRProj = fc.archR)
fragSize_plot 

tss_plot <- plotTSSEnrichment(ArchRProj = fc.archR)
tss_plot

tss_uFrag_plot

ridge_plot <- plotGroups(
  ArchRProj =  fc.archR, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)

## Filter doublets  -------------------------------------------------------------------
fc.archR.2 <- filterDoublets(fc.archR)

doublet_results <- readRDS(paste0("~/", WK_DIR, "/QualityControl/", SAMPLES[1], "/", 
                                  SAMPLES[1], "-Doublet-Summary.RDS"))
doublet_df <- cbind(as.data.frame(table(fc.archR$Sample)), as.data.frame(table(fc.archR.2$Sample)))
doublet_df[3] <- NULL
doublet_df$cells_removed <- 100 - doublet_df[3] / doublet_df[2] * 100
colnames(doublet_df) <- c("Sample", "Pre_DoubRem", "Post_DoubRem", "pc_cells_removed")
doublet_df

##  Dimensionality reduction  ---------------------------------------------------------
fc.archR.2 <- addIterativeLSI(
  ArchRProj = fc.archR.2,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

##  Clustering  -----------------------------------------------------------------------
fc.archR.2 <- addClusters(
  input = fc.archR.2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

##  Visualisation  --------------------------------------------------------------------
fc.archR.2 <- addUMAP(
  ArchRProj = fc.archR.2, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

## Clustering - reporting  ------------------------------------------------------------
# Cluster counts - after Iterative LSI based clustering
clusters_cnts <- as.data.frame(t(as.data.frame(as.vector((table(fc.archR.2$Clusters))))))
rownames(clusters_cnts) <- NULL
colnames(clusters_cnts) <- names(table(fc.archR.2$Clusters))

# Confusion matrix - cell counts per donor
cM_LSI <- confusionMatrix(paste0(fc.archR.2$Clusters), paste0(fc.archR.2$Sample))
clust_CM_LSI <- pheatmap::pheatmap(
  mat = as.matrix(cM_LSI), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f"
)
clust_CM_LSI

# Plot UMAP - for Integrated LSI clusters
clusters_UMAP <- plotEmbedding(ArchRProj = fc.archR.2, colorBy = "cellColData", 
                               name = "Clusters", embedding = "UMAP")
clusters_UMAP_BySample <- plotEmbedding(ArchRProj = fc.archR.2, colorBy = "cellColData", 
                                        name = "Sample", embedding = "UMAP")
cluster_plot <- ggAlignPlots(clusters_UMAP, clusters_UMAP_BySample, type = "h")

## Clusters QC - not required for FC  --------------------------------------------------
# # Remove clusters that have < 10 cells from 2 out of 3 donors 
# cells_to_keep <- rownames(donor_specific_cnts %>% 
#                             filter(Clusters_RNAmapped == "InN-7" | "MG"))
# 
# # Remove cells
# ge.arch
# ge.archR <- ge.archR[cells_to_keep, ]

## Save ArchR project 2  --------------------------------------------------------------
saveArchRProject(ArchRProj = fc.archR.2, 
                 outputDirectory = paste0(WK_DIR, "ArchR/fc.minFrag_3160"), 
                 load = FALSE)

## Batch effect correction - correcting for Sample based batch effects  ---------------
# Batch correct the LSI reduction using harmony save as new reduction named 'Harmony'
fc.archR.3 <- addHarmony(
  ArchRProj = fc.archR.2,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

# Re-cluster using batch corrected LSI reduction data as input - save as 'Clusters_harmony'
fc.archR.3 <- addClusters(
  input = fc.archR.3,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters_harmony",
  resolution = 0.8,
  force = TRUE
)

# Plot UMAP
fc.archR.3 <- addUMAP(
  ArchRProj = fc.archR.3,
  reducedDims = "Harmony",
  name = "UMAPHarmony",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

## Batch effects - reporting  -------------------------------------------------------------
# Cluster counts - after batch corrected Iterative LSI based clustering
clusters_cnts_harmony <- as.data.frame(t(as.data.frame(as.vector((table(fc.archR.3$Clusters_harmony))))))
rownames(clusters_cnts_harmony) <- NULL
colnames(clusters_cnts_harmony) <- names(table(fc.archR.3$Clusters_harmony))

# Confusion matrix - cell counts per donor
cM_harmony <- confusionMatrix(paste0(fc.archR.3$Clusters_harmony), 
                              paste0(fc.archR.3$Sample))
clust_CM_harmony <- pheatmap::pheatmap(
  mat = as.matrix(cM_harmony), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f"
)
clust_CM_harmony

# Plot UMAPs
clusters_UMAP_har <- plotEmbedding(ArchRProj = fc.archR.3, colorBy = "cellColData", 
                                   name = "Clusters_harmony", embedding = "UMAPHarmony")
clusters_UMAP_BySample_har <- plotEmbedding(ArchRProj = fc.archR.3, colorBy = "cellColData", 
                                            name = "Sample", embedding = "UMAPHarmony")
cluster_plot_har <- ggAlignPlots(clusters_UMAP_har, clusters_UMAP_BySample_har, type = "h")

# Confusion matrix to compare LSI based and batch corrected based clusters
cM_harmony_compare <- confusionMatrix(paste0(fc.archR.3$Clusters), 
                                      paste0(fc.archR.3$Clusters_harmony))
clust_CM_harmony_compare <- pheatmap::pheatmap(
  mat = as.matrix(cM_harmony_compare), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f"
)
clust_CM_harmony_compare

# Cluster tree to compare LSI based and batch corrected based clusters
clusttree_harmony_df <- as.data.frame(getCellColData(fc.archR.3,
                                                     select = c("Clusters", 
                                                                "Clusters_harmony")))
colnames(clusttree_harmony_df) <- c("K1", "K2")
clustTree_harmony_plot <- clustree(clusttree_harmony_df, prefix = "K", prop_filter = 0.01)


## NOTE: FROM THIS POINT ON THE HARMONY REDUCED DIM SHOULD BE USED
## THIS IS THE LSI REDUCTION CORRECTED FOR SAMPLE BATCH EFFECTS (USING HARMONY)

##  Integrating snATACseq with snRNAseq  - Cptr 8  ------------------------------------

#   Here we map the RNA-seq cluster IDs to our ATAC-seq clusters. It is a 2-step process.
#   First unconstrained integration is run which is a broad pass at mapping the RNA-seq 
#   cluster IDs to the ATAC-seq clusters. Then the constrained integration is run, this 
#   time by adding supervised cell groupings as IDed in the unconstrained analysis 
#   i.e. InNs (InN-1, InN-2) and ExNs (ExN-1, ExN-2) etc. are clumped. This enables more 
#   accurate cell mapping between the modalites.

# Load Seurat RNA data
seurat.fc <- readRDS("~/Desktop/single_cell/scRNAseq/batch2_CR5_200121/r_objects/final/seurat.pfc.final.rds")
Idents(seurat.fc) # Check cell IDs
seurat.fc$cellIDs <- gsub('FC-', '', seurat.fc$cellIDs)
DimPlot(seurat.fc, label = TRUE) + NoLegend()

#  Run unconstrained integration
fc.archR.3 <- addGeneIntegrationMatrix(
  ArchRProj = fc.archR.3, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony", # Note using Harmony here
  seRNA = seurat.fc,
  addToArrow = FALSE,
  groupRNA = "cellIDs",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

## Unconstrained integration - reporting  ---------------------------------------------
# Confusion matrix - unconstrained cell mappings 
cM_geneExp <- as.matrix(confusionMatrix(fc.archR.3$Clusters_harmony, fc.archR.3$predictedGroup_Un))
clust_CM_geneExp <- pheatmap::pheatmap(
  mat = as.matrix(cM_geneExp), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f"
)

# Get df of top matches cellID matches from RNA for each ATAC cluster
preClust <- colnames(cM_geneExp)[apply(cM_geneExp, 1 , which.max)]
integration_df <- t(as.data.frame(cbind(preClust, rownames(cM_geneExp)))) #Assignments
rownames(integration_df) <- c("RNA", "ATAC")
colnames(integration_df) <- NULL 

# Plot RNA and ATAC UMAPs for comparison
umap_seurat_plot <- DimPlot(seurat.fc, pt.size = 0.2, reduction = "umap", 
                            label = TRUE) + NoLegend()
clusters_UMAP_har_predGroup <- plotEmbedding(ArchRProj = fc.archR.3, colorBy = "cellColData", 
                                   name = "predictedGroup_Un", embedding = "UMAPHarmony")
integration_UMAP_plot <- plot_grid(clusters_UMAP_har, clusters_UMAP_har_predGroup)


# Prepare cell groupings for constrained integration
# Only cell-types in preClust need to be included 
cM_unconstrained <- as.matrix(confusionMatrix(fc.archR.3$Clusters_harmony, fc.archR.3$predictedGroup_Un))
preClust <- colnames(cM_unconstrained)[apply(cM_unconstrained, 1 , which.max)]
cM_unconstrained2 <- cbind(preClust, rownames(cM_unconstrained))
unique(unique(fc.archR.3$predictedGroup_Un))

# Group the unique RNA cluster IDs assigned to ATAC clusters into broad catagories
cExN <- paste0(c("ExN-1", "ExN-2", "ExN-3", "ExN-4", "ExN-5", "ExN-6"), collapse="|")
cInN <- paste0(c("InN-1", "InN-2", "InN-3", "InN-4"), collapse="|")
cRG <- paste0(c("RG-1", "RG-2"), collapse="|")
cMG <- "MG"
cOPC <- "OPC"
cIP <- "IP"
cCycPro <- "CycPro"
cNundef <- "N-undef"

# Pull out the ATAC cluster IDs that map to RNA clusters in unconstrained run
clustExN <- rownames(cM_unconstrained)[grep(cExN, preClust)]
clustInN <- rownames(cM_unconstrained)[grep(cInN, preClust)]
clustRG <- rownames(cM_unconstrained)[grep(cRG, preClust)]
clustMG <- rownames(cM_unconstrained)[grep(cMG, preClust)]
clustNundef<- rownames(cM_unconstrained)[grep(cNundef, preClust)]

# Get unique cells IDs in these categories from the seurat RNA object
rnaExN <- colnames(seurat.fc)[grep(cExN, seurat.fc$cellIDs)]
rnaInN <- colnames(seurat.fc)[grep(cInN, seurat.fc$cellIDs)]
rnaRG <- colnames(seurat.fc)[grep(cRG, seurat.fc$cellIDs)]
rnaMG <- colnames(seurat.fc)[grep(cMG, seurat.fc$cellIDs)]
rnaNundef <- colnames(seurat.fc)[grep(cNundef, seurat.fc$cellIDs)]

# Prepare group list for constrained run
groupList <- SimpleList(
  ExN = SimpleList(
    ATAC = fc.archR.3$cellNames[fc.archR.3$Clusters_harmony %in% clustExN],
    RNA = rnaExN
  ),
  InN = SimpleList(
    ATAC = fc.archR.3$cellNames[fc.archR.3$Clusters_harmony %in% clustInN],
    RNA = rnaInN
  ),    
  RG = SimpleList(
    ATAC = fc.archR.3$cellNames[fc.archR.3$Clusters_harmony %in% clustRG],
    RNA = rnaRG
  ), 
  MG = SimpleList(
    ATAC = fc.archR.3$cellNames[fc.archR.3$Clusters_harmony %in% clustMG],
    RNA = rnaMG
  ), 
  Nundef = SimpleList(
    ATAC = fc.archR.3$cellNames[fc.archR.3$Clusters_harmony %in% clustNundef],
    RNA = rnaNundef
  )
  
)

# Run constrained integration
fc.archR.3 <- addGeneIntegrationMatrix(
  ArchRProj = fc.archR.3, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = seurat.fc,
  addToArrow = FALSE, 
  groupList = groupList,
  groupRNA = "cellIDs",
  nameCell = "predictedCell_Co",
  nameGroup = "predictedGroup_Co",
  nameScore = "predictedScore_Co",
  dimsToUse = 1:26,
  k.score = 26,
  npcs = 26,
  dims = 1:26
  #num.cc = 26
)

## Constrained integration - reporting  -----------------------------------------------
# Compare constrained and unconstrained
pal <- paletteDiscrete(values = seurat.fc$cellIDs)
pal

unconstrained_plot <- plotEmbedding(
  fc.archR.3, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un", 
  pal = pal,
  
)

constrained_plot <- plotEmbedding(
  fc.archR.3, 
  colorBy = "cellColData", 
  name = "predictedGroup_Co", 
  pal = pal,
  
)

integration_compare_plot <- plot_grid(unconstrained_plot, 
                                      constrained_plot, labels = 'AUTO')

# Add integration matrix to ArchR project

#     Add the linked gene expression data to each of the Arrow files.
#     'GroupList' constrains integration and column names to nameCell, 
#     nameGroup, and nameScore for each metadata column we will add 
#     to cellColData. 

fc.archR.3 <- addGeneIntegrationMatrix(
  ArchRProj = fc.archR.3, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = seurat.fc,
  addToArrow = TRUE,
  force= TRUE,
  groupList = groupList,
  groupRNA = "cellIDs",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  dimsToUse = 1:26,
  k.score = 26,
  npcs = 26,
  dims = 1:26
)

# Check available matrices
getAvailableMatrices(fc.archR.3)

# Add imputed weights 
fc.archR.3 <- addImputeWeights(fc.archR.3)

# Get old ATAC-seq clusters IDs
cM_mappings <- confusionMatrix(fc.archR.3$Clusters_harmony, fc.archR.3$predictedGroup)
labelOld <- rownames(cM_mappings)
labelOld

# Identify the cell type from predictedGroup which best defines that each ATAC cluster
labelNew <- colnames(cM_mappings)[apply(cM_mappings, 1, which.max)]
labelNew
unique(labelNew)
atac_clustIDs <- rbind(labelOld, labelNew)

# Map new labels to clusters
fc.archR.3$Clusters_RNAmapped <- mapLabels(fc.archR.3$Clusters_harmony, 
                                           newLabels = labelNew, 
                                           oldLabels = labelOld)

## Final integration - reporting  -----------------------------------------------------
# Confusion matrix to compare old and new cluster mappings
cM_oldNew_compare <- confusionMatrix(paste0(fc.archR.3$Clusters_harmony), 
                                     paste0(fc.archR.3$Clusters_RNAmapped))
clust_CM_oldNew_compare <- pheatmap::pheatmap(
  mat = as.matrix(cM_oldNew_compare), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black", display_numbers = TRUE, number_format =  "%.0f"
)
clust_CM_oldNew_compare

# Cluster tree to compare old and new cluster mappings
clusttree_OldNew_mappings_df <- as.data.frame(getCellColData(fc.archR.3,
                                                             select = c("Clusters_harmony", 
                                                                        "Clusters_RNAmapped")))
colnames(clusttree_OldNew_mappings_df) <- c("K1", "K2")
clustTree_OldNew_mappings_plot <- clustree(clusttree_OldNew_mappings_df, 
                                           prefix = "K", prop_filter = 0.01) 

# Plot UMAPs
clust_UMAP_oldLabel <- plotEmbedding(fc.archR.3, colorBy = "cellColData", embedding = "UMAPHarmony", name = "Clusters_harmony")
clust_UMAP_newLabel <- plotEmbedding(fc.archR.3, colorBy = "cellColData", embedding = "UMAPHarmony", name = "Clusters_RNAmapped")
clust_UMAP_newOldLabel_compare <- plot_grid(clust_UMAP_oldLabel, clust_UMAP_newLabel)

# Donor specific counts after integration
rna_mapped_clusters_df <- as.data.frame(getCellColData(fc.archR.3, select = c("donor", "Clusters_RNAmapped")))
postRNAint_donor_cnts <- rna_mapped_clusters_df %>% dplyr::count(donor, Clusters_RNAmapped) %>%
  pivot_wider(names_from = Clusters_RNAmapped, values_from = n)

# Save archR project 3
saveArchRProject(ArchRProj = fc.archR.3, 
                 outputDirectory = paste0(WK_DIR, "ArchR/fc.minFrag_3160"), 
                 load = FALSE)

## Pseudo-bulk replicates - chtr 9  ---------------------------------------------------

#  The underlying assumption in this process is that the single cells 
#  that are being grouped together are sufficiently similar that we do 
#  not care to understand the differences between them. These cell 
#  groupings are almost always derived from individual clusters or 
#  supersets of clusters that correspond to known cell types. 

fc.archR.4 <- addGroupCoverages(ArchRProj = fc.archR.3, groupBy = "Clusters_RNAmapped")

##  Peak Calling  - cptr 10  ----------------------------------------------------------
# Set macs2 path - note that you need to set the default python env to python 2
pathToMacs2 <- "/Users/Darren/opt/anaconda3/envs/py2/bin/macs2"

# Call peaks
fc.archR.4 <- addReproduciblePeakSet(
  ArchRProj = fc.archR.4, 
  groupBy = "Clusters_RNAmapped", 
  pathToMacs2 = pathToMacs2
)

## Peak Calling - Reporting  ----------------------------------------------------------

# Print peak calling parameters
# Had to cobble code from ArchR repo to generate this - this is printed to screen
# During peak calling but only partially reproduced in log

coverageParams <- fc.archR.4@projectMetadata$GroupCoverages[["Clusters_RNAmapped"]]$Params
coverage_metadata <- fc.archR.4@projectMetadata$GroupCoverages[["Clusters_RNAmapped"]]$coverageMetadata
maxPeaks_default <- 150000
peaksPerCell_default <- 500

tableGroups <- table(getCellColData(fc.archR.4, "Clusters_RNAmapped", drop = TRUE))
peakCallParams_summary_df <- lapply(seq_along(coverageParams$cellGroups), function(y){
  x <- coverageParams$cellGroups[[y]]
  uniq <- unique(unlist(x))
  n <- lapply(x, length) %>% unlist %>% sum
  nmin <- lapply(x, length) %>% unlist %>% min
  nmax <- lapply(x, length) %>% unlist %>% max
  data.frame(
    Group=names(coverageParams$cellGroups)[y], 
    nCells=tableGroups[names(coverageParams$cellGroups)[y]], 
    nCellsUsed=length(uniq), 
    nReplicates=length(x), 
    nMin=nmin, 
    nMax=nmax, 
    maxPeaks = min(maxPeaks_default, length(uniq) * peaksPerCell_default)
  )
}) %>% Reduce("rbind",.)

# Get Granges peak set
fc.peaks <- getPeakSet(fc.archR.4)

# Extract peak matrix
fc.peaks.mat <- getMatrixFromProject(fc.archR.4, useMatrix = "PeakMatrix")

# Get peak matrix meta data
fc.peak_meta <- as.data.frame(colData(fc.peaks.mat))

# Get peak count per cell-type
peak_count_fc <- as.data.frame(fc.peaks$GroupReplicate)
peak_count_fc$Cell_type <- sub("\\..*", "", fc.peaks$GroupReplicate)
peak_count_fc$Donor <- str_extract(fc.peaks$GroupReplicate, "[0-9]{3}")
peak_count_fc[1] <- NULL

peak_count_fc <- peak_count_fc %>% 
  dplyr::count(Cell_type, Donor)

peak_count_fc %>% group_by(Cell_type)

peak_call_perDonor_summary_plot <- ggplot(peak_count_fc, 
                                          aes(fill=Donor, y=n, x=Cell_type)) + 
  geom_bar(position="stack", stat="identity") +
  viridis::scale_fill_viridis(discrete = T) +
  ggtitle("Peak calling summary - peaks per donor") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(title="Donor")) +
  xlab("") +
  ylab("No. of Peaks")


# Plot peak call summary
peak_call_summary <- metadata(fc.archR.4@peakSet)$PeakCallSummary
peak_call_summary_plot <- ggplot(peak_call_summary, 
                                 aes(fill=Var1, y=Freq, x=Group)) + 
  geom_bar(position="stack", stat="identity") +
  viridis::scale_fill_viridis(discrete = T) +
  ggtitle("Peak call summary") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(title="Annotation")) +
  xlab("") +
  ylab(expression("No. of Peaks"~(x10^{"3"})))

# Extract cell-specific peaks
fc_cell_types <- c('ExN.2', 'ExN.3', 'ExN.5', 'ExN.6',
                   'InN.1', 'InN.2', 'InN.3', 
                   'RG.1', 'RG.2', 'MG')

for (CELL_TYPE in fc_cell_types) {
  
  # Load reproducible peak set for each cell-type
  PEAKS <- readRDS(paste0(WK_DIR, 
                          'fc/PeakCalls/', CELL_TYPE, '-reproduciblePeaks.gr.rds'))
  
  # Convert to bed file
  PEAKS_DF <- data.frame(seqnames=seqnames(PEAKS),
                         starts=start(PEAKS)-1,
                         ends=end(PEAKS),
                         names=c(rep(".", length(PEAKS))),
                         scores=c(rep(".", length(PEAKS))),
                         strands=strand(PEAKS))
  
  # Write df to bed file - https://www.biostars.org/p/89341/ 
  write.table(PEAKS_DF, 
              file=paste0(WK_DIR, 'fc/PeakCalls/', CELL_TYPE, '.fc.bed'),
              quote=F, 
              sep="\t", 
              row.names=F, 
              col.names=F)
  
  # Assign Granges object 
  assign(paste0(CELL_TYPE, '_peaks'), PEAKS)
  
}

# Save ArchR project 4
saveArchRProject(ArchRProj = fc.archR.4, 
                 outputDirectory = paste0(WK_DIR, "ArchR/fc.minFrag_3160"), 
                 load = FALSE)

## sLDSC  -----------------------------------------------------------------------------

# This is run on the cluster just testing/reporting defaults here

# Load data
ldsc_results <- read_delim('~/Desktop/single_cell/scATACseq/ldsc/fc_partHerit_SCZ_results.tsv',
                           delim = "\t")

# Extract cell type info (between 2nd and 3rd underscore)
ldsc_results$Category = sapply(strsplit(ldsc_results$Category, "_"), function(x) x[3])


## Create markdown doc  ---------------------------------------------------------------
render("~/Desktop/single_cell/scATACseq/markdown/scATAC_ArchR_FC_fragsQC_3160.Rmd")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
