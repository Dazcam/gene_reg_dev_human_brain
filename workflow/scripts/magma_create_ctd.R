# -------------------------------------------------------------------------------------
#
#    MAGMA Celltyping - create CTD object
#
# -------------------------------------------------------------------------------------

##  Resources  -------------------------------------------------------------------------
   
    #  Paper - https://www.nature.com/articles/s41588-018-0129-5
    #  Github - https://github.com/NathanSkene/MAGMA_Celltyping

    # ECWE - for creating the ctd object
      # Github - https://github.com/neurogenomics/EWCE/
      # Vingette - https://nathanskene.github.io/EWCE/articles/EWCE.html

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
library(EWCE) # using srun can't access Biomart
library(MAGMA.Celltyping) # Note the "." instead of "_" - using srun can't see magma executable
library(One2One) # Not required for human data I think
library(argparser)
library(reshape2)

## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("Read brain region for magma ... \n")
p <- add_argument(p, "region", help = "No brain region provided")
p <- add_argument(p, "seurat_obj", help = "No Seurat obj provided")
args <- parse_args(p)
print(args)

##  Define variables  -----------------------------------------------------------------
REGION <- args$region
SEURAT_OBJ <- args$seurat_obj

MAGMA_DIR <- '/scratch/c.c1477909/magma_celltyping/'

GENOME_DIR <- paste0(MAGMA_DIR, 'g1000_eur')
CTD_DIR <- paste0(MAGMA_DIR, 'ctd_objects')

##  Load seurat objects  ------------------------------------------------------------
cat(paste0("\nLoading seurat object for ", toupper(REGION), " ... \n"))
seurat.obj <- readRDS(SEURAT_OBJ)

# Create raw count gene matrix .csv - needs to be cells x genes
cat("Creating count by gene matrix ...\n")
raw_counts <- as.matrix(seurat.obj@assays$RNA@counts)

# Create annotations - required “cell_id”, “level1class” and “level2class”
cat("Creating annotations ...\n")
annotations <- annotations <- as.data.frame(cbind(colnames(seurat.obj), 
                                                  seurat.obj[[c('cellIDs')]], 
                                                  seurat.obj[[c('cellIDs')]]))
colnames(annotations) <- c('cell_id', 'level1class', 'level2class')
rownames(annotations) <- NULL

##  Generate CTD object  --------------------------------------------------------------
  
  # Uses EWCE package - this creates an object called ctd 
  # Note: option here to scTransform data to correct for cell size - not done this
  # drop.uninformative.genes drops genes that do not show sig variance
  # between level 2 celltypes (based on ANOVA) - not necessary as we are only using
  # 1 level
  # generate.celltype.data calculates the cell specific averages for each gene

cat("Creating celltype data ... \n")
exp_DROPPED <- drop.uninformative.genes(exp = raw_counts, 
                                            level2annot = annotations$level2class)
annotLevels <- list(level1class = annotations$level2class, 
                    level2class = annotations$level2class)
ctd <- generate.celltype.data(exp = exp_DROPPED, 
                              annotLevels = annotLevels, 
                              groupName = REGION,
                              savePath = CTD_DIR)
print(ctd)
cat(paste0(toupper(REGION), " ctd ... created.\n"))

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
