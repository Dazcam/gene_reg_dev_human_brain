#--------------------------------------------------------------------------------------
#
#    snRNAseq GENIE3 analysis - laptop
#
#--------------------------------------------------------------------------------------

##  Resources  -------------------------------------------------------------------------

#  Paper - https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0012776
#  Github - https://github.com/vahuynh/GENIE3
#  Vingette - https://bioconductor.org/packages/devel/bioc/vignettes/GENIE3/inst/doc/GENIE3.html

## Requirements  -------------------------------------------------------------------------

# Required on Hawk before opening R
# module load libgit2/1.1.0
# module load R/4.0.3

## Initialise R library  --------------------------------------------------------------
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )

##  Load Packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(Seurat)
library(GENIE3)

##  Define variables  -----------------------------------------------------------------
DATA_DIR <- "/scratch/c.c1477909/R_objects/"
#REGIONS <- c("cer", "hip", "pfc", "tha", "wge")
REGION <- "pfc"

## Load Data --------------------------------------------------------------------------
cat(paste0('\nLoading ', REGION, ' seurat object ... \n'))
seurat.obj <- readRDS(paste0(DATA_DIR, 'seurat.', REGION, '.final.rds'))

set.seed(123) # For reproducibility of results

# Subset by cell-type
cat(paste0('\nSubsetting ', REGION, ' seurat object ... \n'))
FC_ExN_2 <- subset(seurat.pfc, idents = 'FC-ExN-2')

# Prepare count matrix
cat(paste0('\nCreating count matrix ', REGION, ' seurat object ... \n'))
counts <- as.matrix(FC_ExN_2@assays$RNA@counts)

# Run GENIE3 to create matrix of weights of putative regulatory links
cat(paste0('\nRunning GENIE3 ... \n'))
weightMat <- GENIE3(counts, nCores = 24)

# Get link list
cat(paste0('\nObtaining link list ... \n'))
linkList <- getLinkList(weightMat)

# Save data
cat(paste0('\nSaving data ... \n'))
write_tsv(weightMat, 
          paste0(DATA_DIR, "snRNAseq_GENIE3_", REGION, "_FC_ExN_2_weight_matrix.txt"), 
          escape = "double")
write_tsv(linkList, 
          paste0(DATA_DIR, "snRNAseq_GENIE3_", REGION, "_FC_ExN_2_link_list.txt"), 
          escape = "double")

cat('\nDone.')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
