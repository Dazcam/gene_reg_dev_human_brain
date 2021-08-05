#--------------------------------------------------------------------------------------
#
#    snRNAseq analysis - QC 
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

## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("\nRead brain region and output directory for snRNAseq QC ... \n")
p <- add_argument(p, "region", help = "No brain region specified")
p <- add_argument(p, "in_dir", help = "No input data directory specified")
p <- add_argument(p, "out_obj", help = "No output object and path specified")
args <- parse_args(p)
print(args)

##  Define variables  -----------------------------------------------------------------
REGION <- args$region
DATA_DIR <- args$in_dir
OUT_OBJ <- args$out_obj
GWAS <- args$gwas

## Load and munge sample specific data  ----------------------------------------------
if (REGION == 'Cer') {
  
  SAMPLES <- c("510_Cer_B2", "611_Cer_B2", "993_Cer_B3")
  
  sce <- DropletUtils::read10xCounts(
    
    c( 
      paste0(DATA_DIR, "510_Cer_B2/filtered_feature_bc_matrix"), 
      paste0(DATA_DIR, "611_Cer_B2/filtered_feature_bc_matrix"), 
      paste0(DATA_DIR, "993_Cer_B3/filtered_feature_bc_matrix")
    ),
    sample.names = SAMPLES, type = "auto"
  )
  
  sce@colData$batch <- ifelse(grepl("B2", sce@colData$Sample), "B2", "B3")
  sce@colData$region <- substr(sce@colData$Sample, 5, 7)
  sce@colData$donor <- substr(sce@colData$Sample, 1, 3)
  
} else if (REGION == 'Hip') {
  
  SAMPLES <- c("510_Hip_B2", "611_Hip_B2", "993_Hip_B2")
  
  sce <- DropletUtils::read10xCounts(
    
    c( 
      paste0(DATA_DIR, "510_Hip_B2/filtered_feature_bc_matrix"), 
      paste0(DATA_DIR, "611_Hip_B2/filtered_feature_bc_matrix"), 
      paste0(DATA_DIR, "993_Hip_B2/filtered_feature_bc_matrix")
    ),
    sample.names = SAMPLES, type = "auto"
  )
  
  sce@colData$batch <- ifelse(grepl("B2", sce@colData$Sample), "B2", "B3")
  sce@colData$region <- substr(sce@colData$Sample, 5, 7)
  sce@colData$donor <- substr(sce@colData$Sample, 1, 3)
    
  
} else if (REGION == 'PFC') {
  
  SAMPLES <- c("510_PFC_B1", "611_PFC_B1", "993_PFC_B2")
  
  sce <- DropletUtils::read10xCounts(
    
    c( 
      paste0(DATA_DIR, "510_PFC_B1/filtered_feature_bc_matrix"), 
      paste0(DATA_DIR, "611_PFC_B1/filtered_feature_bc_matrix"), 
      paste0(DATA_DIR, "993_PFC_B2/filtered_feature_bc_matrix")
    ),
    sample.names = SAMPLES, type = "auto"
  )
  
  sce@colData$batch <- ifelse(grepl("B1", sce@colData$Sample), "B1", "B2")
  sce@colData$region <- substr(sce@colData$Sample, 5, 7)
  sce@colData$donor <- substr(sce@colData$Sample, 1, 3)
  
} else if (REGION == 'Tha') {
  
  SAMPLES <- c("510_Tha_B3", "611_Tha_B2", "993_Tha_B2")
  
  sce <- DropletUtils::read10xCounts(
    
    c( 
      paste0(DATA_DIR, "510_Tha_B3/filtered_feature_bc_matrix"), 
      paste0(DATA_DIR, "611_Tha_B2/filtered_feature_bc_matrix"), 
      paste0(DATA_DIR, "993_Tha_B2/filtered_feature_bc_matrix")
    ),
    sample.names = SAMPLES, type = "auto"
  )
  
  sce@colData$batch <- ifelse(grepl("B2", sce@colData$Sample), "B2", "B3")
  sce@colData$region <- substr(sce@colData$Sample, 5, 7)
  sce@colData$donor <- substr(sce@colData$Sample, 1, 3)
  
} else {
  
  SAMPLES <- c("510_WGE_B1", "611_WGE_B1", "993_WGE_B2")
  
  sce <- DropletUtils::read10xCounts(
    
    c( 
      paste0(DATA_DIR, "510_WGE_B1/filtered_feature_bc_matrix"), 
      paste0(DATA_DIR, "611_WGE_B1/filtered_feature_bc_matrix"), 
      paste0(DATA_DIR, "993_WGE_B2/filtered_feature_bc_matrix")
    ),
    sample.names = SAMPLES, type = "auto"
  )
  
  sce@colData$batch <- ifelse(grepl("B1", sce@colData$Sample), "B1", "B2")
  sce@colData$region <- substr(sce@colData$Sample, 5, 7)
  sce@colData$donor <- substr(sce@colData$Sample, 1, 3)
  
  }

## Adjust rownames to GeneIDs and colnames to cell barcodes
# Can probs omit batch info here to join samples - what happens with duplicated columns?
# sampleIDs <- substr(sce@colData$Sample, 1, 7)
# Check - length(unique(colnames(sce))) == dim(sce)[2]
colnames(sce) <- paste(substr(sce@colData$Sample, 1, 7), colData(sce)$Barcode, sep = "_")
rownames(sce) <- rowData(sce)$Symbol

# Check for duplicates - for now remove these
# See here for alternatives:
# http://biocworkshops2019.bioconductor.org.s3-website-us-east-1.amazonaws.com/page/OSCABioc2019__OSCABioc2019/
dim(sce)
duplicate_genes <- duplicated(rownames(sce))
duplicate_geneIDs <- rownames(sce)[duplicated(rownames(sce))]
sce <- sce[!duplicate_genes, ]
dim(sce)

##  Inital QC  ------------------------------------------------------------------------
# Add MT, Ribo, Exp and gene QCs
# Gene IDs to rownames
mt_genes <- grepl("^MT-", rowData(sce)$Symbol)
ribo_genes <- grepl("^RP[LS]", rowData(sce)$Symbol)
hb_genes <- rownames(sce)[grep("^HB[^(P)]", rownames(sce))]

sce <- addPerCellQC(sce, subsets=list(Mito=mt_genes, Ribo=ribo_genes, HB = hb_genes))
sce <- addPerFeatureQC(sce)
sce$log10sum <- log10(sce$sum)

# Pre-QC stats
sce_metadata <- as.data.frame(colData(sce))
tot_cells_preQC <- dim(sce)[2]
tot_genes_preQC <- dim(sce)[1]
sample_cells_preQC <- sce_metadata %>% group_by(Sample) %>% dplyr::count() 
region_cells_preQC <- sce_metadata %>% group_by(region) %>% dplyr::count() 

##  Initial plotting  -----------------------------------------------------------------
genes_plot_preQC <- plotColData(sce, y = "detected", x = "Sample", colour_by = "Sample") + 
  geom_hline(yintercept = 1000, color = "red") +
  geom_hline(yintercept = 5000, color = "red") +
  theme(axis.text.x = element_text(angle = 90))

reads_plot_preQC <- plotColData(sce, y = "log10sum", x = "Sample", colour_by = "Sample") +
  theme(axis.text.x = element_text(angle = 90))

mito_plot_preQC <- plotColData(sce, y = "subsets_Mito_percent", x = "Sample", colour_by = "Sample") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 5, color = "red") 

ribo_plot_preQC <- plotColData(sce, y = "subsets_Ribo_percent", x = "Sample", colour_by = "Sample") +
  theme(axis.text.x = element_text(angle = 90), legend.position="none") +
  geom_hline(yintercept = 10, color = "red") 

hb_plot_preQC <- plotColData(sce, y = "subsets_HB_percent", x = "Sample", colour_by = "Sample") +
  theme(axis.text.x = element_text(angle = 90))

reads_by_genes_preQC <- plotColData(sce, x = "total", y = "detected", colour_by = "Sample")

# Relative expression of each gene per cell
C.sce = counts(sce)
C.sce@x = C.sce@x/rep.int(colSums(C.sce), diff(C.sce@p))
most_expressed <- order(Matrix::rowSums(C.sce), decreasing = T)[50:1]
boxplot(as.matrix(t(C.sce[most_expressed, ])), cex = 0.05, 
        las = 1, xlab = "% total count per cell", cex.axis=0.5,
        col = (scales::hue_pal())(50)[50:1], horizontal = TRUE)

##  Set QC thresholds  ----------------------------------------------------------------
mito_thresh <- 5
genesInCells_thresh <- 3
gene_low_thresh <- 1000
gene_high_thresh <- 5000
ribo_thresh <- 10

sce$mito_discard <- sce$subsets_Mito_percent > mito_thresh
sce$ribo_discard <- sce$subsets_Ribo_percent > ribo_thresh
rowData(sce)$genesInXcells <- rowSums(counts(sce)) > genesInCells_thresh
sce$cellsWithXgenesLow <- sce$detected < gene_low_thresh
sce$cellsWithXgenesHigh <- sce$detected > gene_high_thresh

qc_thresholds <- data.frame(Mito=mito_thresh,
                            Ribo=ribo_thresh,
                            Gene_in_X_cells=genesInCells_thresh,
                            Gene_low_thresh=gene_low_thresh,
                            Gene_high_thresh=gene_high_thresh)


##  Detect Doublets  ------------------------------------------------------------------
set.seed(123)
sce <- scDblFinder(sce, samples=sce$Sample)
sce$doublet_logic <- ifelse(sce$scDblFinder.class == "doublet", TRUE, FALSE)
plotDoubletMap(sce)
table(c(sce$scDblFinder.class, sce$Sample))

sce_metadata <- as.data.frame(colData(sce))

doublet_stats_df <- sce_metadata %>%
  group_by(Sample) %>%
  dplyr::count(doublet_logic) %>%
  pivot_wider(names_from = doublet_logic, values_from = n) %>%
  dplyr::rename(Singlet = `FALSE`) %>% dplyr::rename(Doublet = `TRUE`) %>%
  mutate(Total=(Doublet+Singlet)) %>%
  mutate(Prop=(Doublet/Singlet)*100)

##  Calculate discarded numbers/proportions  ------------------------------------------
# Note this doesn't trim read counts at all - does Seurat?
discard <- sce$cellsWithXgenesLow | sce$cellsWithXgenesHigh | sce$mito_discard | sce$ribo_discard | sce$doublet_logic
sce$discard <- discard
cells_removed_df <- DataFrame(XgenesLow=sum(sce$cellsWithXgenesLow),
                              XgenesHigh=sum(sce$cellsWithXgenesHigh),
                              MitoProp=sum(sce$mito_discard),
                              RiboProp=sum(sce$ribo_discard),
                              Doublets=sum(sce$doublet_logic), 
                              Total=sum(discard),
                              PropAllCells=(sum(discard)/dim(sce)[2]*100))

sce_metadata <- as.data.frame(colData(sce))

discard_perSample_df <- sce_metadata %>% 
  group_by(Sample) %>%
  dplyr::count(discard) %>%
  pivot_wider(names_from = discard, values_from = n) %>%
  dplyr::rename(Keep = `FALSE`) %>% dplyr::rename(Discard = `TRUE`) %>%
  mutate(Pre_QC=(Keep+Discard)) %>%
  mutate(Prop_discard=(Discard/Pre_QC)*100) %>%
  relocate(Pre_QC, .after = Sample)

discard_perRegion_df <- sce_metadata %>% 
  group_by(region) %>%
  dplyr::count(discard) %>%
  pivot_wider(names_from = discard, values_from = n) %>%
  dplyr::rename(Keep = `FALSE`) %>% dplyr::rename(Discard = `TRUE`) %>%
  mutate(Pre_QC=(Keep+Discard)) %>%
  mutate(Prop_discard=(Discard/Pre_QC)*100)  %>%
  relocate(Pre_QC, .after = region)

##  Hard filters  ---------------------------------------------------------------------
# Cells
dim(sce)

sce.filt <- sce[rowData(sce)$genesInXcells, !sce$discard]

dim(sce.filt)

# Genes
sce.filt <- sce.filt[!grepl("MALAT1", rownames(sce.filt)), ]
sce.filt <- sce.filt[!grepl("^MT-", rownames(sce.filt)), ]

discard_total_df <- data.frame(Pre_QC=c(dim(sce)[1], dim(sce)[2]),
                               Post_QC=c(dim(sce.filt)[1], dim(sce.filt)[2]),
                               Prop_discard=c((dim(sce)[1]-dim(sce.filt)[1])/dim(sce)[1]*100, 
                                              (dim(sce)[2]-dim(sce.filt)[2])/dim(sce)[2]*100),
                               row.names = c("Genes", "Reads"))


##  Post-QC plots  --------------------------------------------------------------------
genes_plot_postQC <- plotColData(sce.filt, y = "detected", x = "Sample", colour_by = "Sample") + 
  theme(axis.text.x = element_text(angle = 90))
genes_QC_plot <- plot_grid(genes_plot_preQC + theme(legend.position="none"), 
                           genes_plot_postQC + theme(legend.position="none"), 
                           ncol = 1) 

reads_plot_postQC <- plotColData(sce.filt, y = "log10sum", x = "Sample", colour_by = "Sample") +
  theme(axis.text.x = element_text(angle = 90))
reads_QC_plot <- plot_grid(reads_plot_preQC + theme(legend.position="none"), 
                           reads_plot_postQC + theme(legend.position="none"),
                           ncol = 1)

mito_plot_postQC <- plotColData(sce.filt, y = "subsets_Mito_percent", 
                                x = "Sample", colour_by = "Sample") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 5, color = "red") 
mito_QC_plot <- plot_grid(mito_plot_preQC + theme(legend.position="none"),
                          mito_plot_postQC + theme(legend.position="none"), 
                          ncol = 1)

ribo_plot_postQC <- plotColData(sce.filt, y = "subsets_Ribo_percent", 
                                x = "Sample", colour_by = "Sample") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 10, color = "red") 

C.sce.filt = counts(sce.filt)
C.sce.filt@x = C.sce.filt@x/rep.int(colSums(C.sce.filt), diff(C.sce.filt@p))
most_expressed_postQC <- order(Matrix::rowSums(C.sce.filt), decreasing = T)[50:1]
boxplot(as.matrix(t(C.sce.filt[most_expressed_postQC, ])), cex = 0.05, 
        las = 1, xlab = "% total count per cell", cex.axis=0.5,
        col = (scales::hue_pal())(50)[50:1], horizontal = TRUE)

saveRDS(object = sce.filt, file = OUT_OBJ)


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
