#--------------------------------------------------------------------------------------
#
#    snRNAseq SCENIC analysis - laptop
#
#--------------------------------------------------------------------------------------

##  Resources  -------------------------------------------------------------------------

# Paper - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5937676/
# Github - https://github.com/aertslab/SCENIC

## Info  -------------------------------------------------------------------------

# Scenic servers were down so had to copy TF and DB feather files over from Hawk
# Still only using a single feather file (for testing) but 2 are required
# Authors recommend using pscenic as the R implemetaion is too slow

##  Load Packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(Seurat)
library(SCENIC)

##  Define variables  -----------------------------------------------------------------
DATA_DIR <- '~/Desktop/single_cell/scRNAseq/batch2_CR5_200121/r_objects/final/'
DB_DIR <- "~/Desktop/single_cell/scRNAseq/scenic"
#REGIONS <- c("cer", "hip", "pfc", "tha", "wge")
REGION <- "pfc"

## Load Data --------------------------------------------------------------------------
cat(paste0('\nLoading ', REGION, ' seurat object ... \n'))
seurat.obj <- readRDS(paste0(DATA_DIR, 'seurat.', REGION, '.final.rds'))

exprMat <- as.matrix(seurat.obj@assays$RNA@counts)
cellInfo <- data.frame(seuratCluster=Idents(seurat.obj))

## Initialize settings  --------------------------------------------------------------
# Server is down need to use 2 DBs files
# Need to rename default DB feather files
scenicOptions <- initializeScenic(org = "hgnc", 
                                  dbs = 'encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.feather', 
                                  dbDir = DB_DIR, nCores = 2)
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo

### Co-expression network - STEP 1
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)


#### ABORTED - Ran on laptop no progress after ~24 hrs #####

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
