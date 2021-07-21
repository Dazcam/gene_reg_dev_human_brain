# -------------------------------------------------------------------------------------
#
#    MAGMA Celltyping 3
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

## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("Read brain region for magma ... \n")
p <- add_argument(p, "region", help = "No brain region provided")
p <- add_argument(p, "ctd_obj", help = "No CTD obj provided")
p <- add_argument(p, "gwas", help = "No GWAS sumstats path provided")
args <- parse_args(p)
print(args)

##  Define variables  -----------------------------------------------------------------
REGION <- args$region
CTD_OBJ <- args$ctd_obj
GWAS <- args$gwas

MAGMA_DIR <- '/scratch/c.c1477909/magma_celltyping/'

GENOME_DIR <- paste0(MAGMA_DIR, 'g1000_eur/g1000_eur')
GWAS_PATH <- GWAS
CTD_DIR <- paste0(MAGMA_DIR, 'ctd_objects')
#MAGMA_OUTDIR <- paste0(MAGMA_DIR, 'MAGMA_output/PGC3_SCZ_hg38_magmaReady.tsv.10UP.1.5DOWN/')
#FIGS_OUTDIR <- paste0(MAGMA_DIR, 'MAGMA_Figures/*')  # * needed due to space in Figs dir name!
#LINEAR_OUTDIR <- paste0(MAGMA_OUTDIR, REGION, "_linear")
#TOP10_OUTDIR <- paste0(MAGMA_OUTDIR, REGION, "_top10")


##  Magma celltyping linear and top10% association analyses  ------------------------------
#   1-2 min on Hawk home  

##  Load ctd  -----------------------------------------------------------------------
# Load ctd data - Note this step is not in the vingette!! Loads data in as 'ctd'
cat(paste0("\n\n\n\n\nLoading ctd for ", REGION, " ... \n\n\n"))
load(CTD_OBJ)

##  Prepare quantile groups for celltypes  ------------------------------------------
ctd <- prepare.quantile.groups(ctd, specificity_species = "human", numberOfBins = 40)

# Examine how the quantile groups look
print(ctd[[1]]$specificity_quantiles[c("GFAP","DLG4"),])
print(table(ctd[[1]]$specificity_quantiles[, 1]))

##  Run the main cell type association analysis  ------------------------------------
##  Linear mode  --------------------------------------------------------------------
cat(paste0("\n\n\n\n\nRunning linear association analysis for ", REGION, " ...\n\n\n\n\n"))
ctAssocsLinear <- calculate_celltype_associations(ctd = ctd,
                                                  gwas_sumstats_path = GWAS_PATH,
                                                  genome_ref_path = GENOME_DIR,
                                                  specificity_species = "human",
                                                  analysis_name = paste0(REGION, "_linear"))
FigsLinear <- plot_celltype_associations(ctAssocs = ctAssocsLinear,
                                         ctd = ctd)

# Move files to region specific linear directory
#system(paste0("mv ", outDir, "*.MainRun.* ", LINEAR_OUTDIR))
#cat(paste0("Moving magma_celltyping output to linear assoc. folder\n", SAMPLE, " ... DONE"))

##  Top 10% mode - added specificity species param, not in vingette  ----------------
cat(paste0("\n\n\n\n\nRunning top10% association analysis for ", REGION, " ...\n\n\n\n\n"))
ctAssocsTop <- calculate_celltype_associations(ctd = ctd,
                                               gwas_sumstats_path = GWAS_PATH,
                                               genome_ref_path = GENOME_DIR,
                                               EnrichmentMode = "Top 10%",
                                               specificity_species = "human",
                                               analysis_name = paste0(REGION, "_top10"))
FigsTopDecile <- plot_celltype_associations(ctAssocs=ctAssocsTop,
                                            ctd = ctd)

# Move files to region specific top10 directory
#system(paste0("mv ", outDir, "*.MainRun.* ", sample_top10Dir))
#cat(paste0("/n/n/n/n/nMoving magma_celltyping output to top10 assoc. folder\n", SAMPLE, " ... "))


##  Plot linear together with the top decile mode  ----------------------------------
ctAssocMerged = merge_magma_results(ctAssoc1=ctAssocsLinear,
                                    ctAssoc2=ctAssocsTop)
FigsMerged = plot_celltype_associations(ctAssocs=ctAssocMerged,
                                        ctd=ctd)

# Move plots to region specific directory
#system(paste0('mv ', FIGS_OUTDIR, '/*.pdf ', FIGS_OUTDIR, '/', SAMPLE))
#cat(paste0("Moving figures to region specific directory ", SAMPLE, " ... DONE"))


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
