# -------------------------------------------------------------------------------------
#
#    MAGMA Celltyping - Prepare gene lists for LDSC - 10 Quantiles
#
# -------------------------------------------------------------------------------------

##  Resources  -------------------------------------------------------------------------

    #  Paper - https://www.nature.com/articles/s41588-018-0129-5
    #  Github - https://github.com/NathanSkene/MAGMA_Celltyping

    # ECWE - for creating the ctd object
      # Github - https://github.com/neurogenomics/EWCE/

## Requirements  -------------------------------------------------------------------------

    # Required on Hawk before opening R
      # module load libgit2/1.1.0
      # module load R/4.0.3

    # GWAS:
      # SNP, CHR, BP as first three columns.
      # GWAS has one of: Z, OR, BETA, LOG_ODDS, SIGNED_SUMSTAT
      # GWAS has all of: SNP, CHR, BP, P, A1, A2

    # Net access:
      # Uses BiomaRt for annotations need to run SLURM jobs locally

## Initialise R library  --------------------------------------------------------------
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )

## Load packages  ---------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(Seurat)
library(devtools)
#library(EWCE) # using srun can't access Biomart
library(MAGMA.Celltyping) # Note the "." instead of "_" - using srun can't see magma executable
#library(One2One) # Not required for human data I think
library(argparser)
library(reshape2)
library(tidyverse)

## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("Read brain region for magma ... \n")
p <- add_argument(p, "region", help = "No brain region provided")
p <- add_argument(p, "ctd_obj", help = "No CTD obj provided")
args <- parse_args(p)
print(args)

##  Define variables  -----------------------------------------------------------------
REGION <- args$region
CTD_OBJ <- args$ctd_obj

MAGMA_DIR <- '/scratch/c.c1477909/magma_celltyping/'

##  Magma celltyping linear and top10% association analyses  ------------------------------
#   1-2 min on Hawk home  

##  Load ctd  -----------------------------------------------------------------------
# Load ctd data - Note this step is not in the vingette!! Loads data in as 'ctd'
cat(paste0("\n\n\n\n\nLoading ctd for ", REGION, " ... \n\n\n"))
load(CTD_OBJ)

##  Prepare quantile groups for celltypes  ------------------------------------------
ctd_10q <- prepare.quantile.groups(ctd, specificity_species = "human", numberOfBins = 10)

## Extract gene lists - 10 quantiles (0-10) for each cell type  ---------------------
# Get quantiles for brain region - alt from matrix to dataframe
ctd_10q_quantiles <- as.data.frame(ctd_10q[[1]]$specificity_quantiles)

# Get cell types for brain region
cell_types <- colnames(ctd_10q_quantiles)

# Get gene list - 1 list per quantile per cell type
for (CELL_TYPE in cell_types) {
  
  cat(paste0('\nExtracting genes for ', CELL_TYPE, ' ...\n'))
  
  for (QUANTILE in seq(0, 10)) {
    
    CELL_TYPE_DF <- ctd_10q_quantiles[CELL_TYPE] %>% 
      filter(.data[[CELL_TYPE]] == QUANTILE)
      
    GENE_LIST_HEAD <- paste(head(rownames(CELL_TYPE_DF)), collapse=", ")
    GENE_LIST_LENGTH <- length(rownames(CELL_TYPE_DF))
    
    cat(paste0('Q', QUANTILE, '. ', GENE_LIST_LENGTH, ' genes. 1st 6 genes:  ', 
               GENE_LIST_HEAD, '.\n'))
    
    write.table(rownames(CELL_TYPE_DF), paste0('/scratch/c.c1477909/magma_celltyping/ldsc_gene_lists/', 
                       CELL_TYPE, '_Q', QUANTILE, '_genes.tsv'), quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    
  }
  

}

# Create empty file to alert snakemake that job is complete
file.create(paste0('/scratch/c.c1477909/magma_celltyping/ldsc_gene_lists/',
                   REGION, '_complete.file'))



# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
