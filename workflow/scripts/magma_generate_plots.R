# -------------------------------------------------------------------------------------
#
#    MAGMA Celltyping plots and data-analysis
#
# -------------------------------------------------------------------------------------

## Initialise R library  --------------------------------------------------------------
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)
library(rmarkdown)

## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("Read magma input directory ... \n")
p <- add_argument(p, "magma_dir", help = "No markdown magma dir specified")
args <- parse_args(p)
print(args)


##  Initialise variables  -------------------------------------------------------------
REGIONS <- c('cer', 'hip', 'pfc', 'tha', 'wge')
GWAS <- c('ADHD', 'BPD', 'MDD', 'SCZ')
MAGMA_DIR <- args$magma_dir
NAME_BODY <- '_hg38_magma_ready.sumstats.tsv.10UP.1.5DOWN'
SUFFIX_LINEAR <- '_linear.gsa.out'
SUFFIX_TOP10 <- '_top10.gsa.out'


# Linear results
for (REGION in REGIONS) {
  
  for (DISORDER in GWAS) {
    
    ##  Load Data  ----------------------------------------------------------------------
    
    #   Need to skip the first 4 columns in datafile 
    
    linearData <- read.table(paste0(MAGMA_DIR, DISORDER, NAME_BODY, '/', DISORDER, NAME_BODY, 
                                     '.level1.', REGION, SUFFIX_LINEAR), header = FALSE)
    names(linearData) <- as.matrix(linearData[1, ])
    linearData <- linearData[-1, ]
    as.numeric(linearData$P)
    -log10(as.numeric(linearData$P))
    print(head(linearData))
    
    ##  Plot  ---------------------------------------------------------------------------
    linearPlot <- ggplot(data = linearData, aes(x = -log10(as.numeric(P)), y = VARIABLE)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_bw() +
      ggtitle(toupper(REGION)) +
      theme(plot.title = element_text(hjust = 0.5))
    
    assign(paste0(REGION, '_', DISORDER, '_magma_linear_plot'), linearPlot, envir = .GlobalEnv)
    assign(paste0(REGION, '_', DISORDER, '_magma_linear_data'), linearData, envir = .GlobalEnv)
    
  }
  
}


# Top 10% results
for (REGION in REGIONS) {
  
  for (DISORDER in GWAS) {
    
    ##  Load Data  ----------------------------------------------------------------------
    
    #   Need to skip the first 4 columns in datafile 
    
    top10Data <- read.table(paste0(MAGMA_DIR, DISORDER, NAME_BODY, '/', DISORDER, NAME_BODY, 
                                    '.level1.', REGION, SUFFIX_TOP10), header = FALSE)
    names(top10Data) <- as.matrix(top10Data[1, ])
    top10Data <- top10Data[-1, ]
    as.numeric(top10Data$P)
    -log10(as.numeric(top10Data$P))
    print(head(top10Data))
    
    ##  Plot  ---------------------------------------------------------------------------
    top10Plot <- ggplot(data = top10Data, aes(x = -log10(as.numeric(P)), y = VARIABLE)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_bw() +
      ggtitle(toupper(REGION)) +
      theme(plot.title = element_text(hjust = 0.5))
    
    assign(paste0(REGION, '_', DISORDER, '_magma_top10_plot'), top10Plot, envir = .GlobalEnv)
    assign(paste0(REGION, '_', DISORDER, '_magma_top10_data'), top10Data, envir = .GlobalEnv)
    
  }
  
}

# Create group plots
for (DISORDER in GWAS) {

  magma_linear_plot <- plot_grid(get(paste0('cer_', DISORDER, '_magma_linear_plot')),
                                 get(paste0('hip_', DISORDER, '_magma_linear_plot')), 
                                 get(paste0('pfc_', DISORDER, '_magma_linear_plot')),
                                 get(paste0('tha_', DISORDER, '_magma_linear_plot')),
                                 get(paste0('wge_', DISORDER, '_magma_linear_plot')))
  
  assign(paste0('all_regions_', DISORDER, '_magma_linear_plot'), magma_linear_plot, envir = .GlobalEnv)
  
  magma_top10_plot <- plot_grid(get(paste0('cer_', DISORDER, '_magma_top10_plot')),
                                 get(paste0('hip_', DISORDER, '_magma_top10_plot')), 
                                 get(paste0('pfc_', DISORDER, '_magma_top10_plot')),
                                 get(paste0('tha_', DISORDER, '_magma_top10_plot')),
                                 get(paste0('wge_', DISORDER, '_magma_top10_plot')))
  
  assign(paste0('all_regions_', DISORDER, '_magma_top10_plot'), magma_top10_plot, envir = .GlobalEnv)
                                        
}

## Rare variants  ---------------------------------------------------------------------
# Rare variants
library(dittoSeq)

rare_genes <- c('SETD1A', 'CUL1', 'XPO7', 'TRIO', 'CACNA1G',
                'SP4', 'GRIA3', 'GRIN2A', 'HERC1', 'RB1CC1',
                'HCN4', 'AKAP11', 'ZNF136', 'SRRM2', 'NR3C2',
                'ZMYM2', 'FAM120A', 'SLF2', 'KDM6B', 'DNM3')


for (REGION in REGIONS) {
  
  ROBJ_DIR <- '~/Desktop/single_cell/scRNAseq/batch2_CR5_200121/r_objects/final/'
  
  seurat.obj <- readRDS(paste0(ROBJ_DIR, 'seurat.', REGION, '.final.rds'))

  
  # Create heatmap
  # scaled.to.max normalisee all expression data to the max expression of each gene [0,1],
  # which is often useful for zero-enriched single-cell data.
  heatmap <- dittoHeatmap(seurat.obj, rare_genes,
               annot.by = c('cellIDs'), scaled.to.max = TRUE)
  dotplot <- dittoDotPlot(seurat.obj, vars = rare_genes, group.by = 'cellIDs')
  
  assign(paste0(REGION, '_heatmap'), heatmap, envir = .GlobalEnv)
  assign(paste0(REGION, '_dotplot'), dotplot, envir = .GlobalEnv)
  
  
  
}


# Render markdown report
render("/scratch/c.c1477909/markdown/magma_celltyping_plots.Rmd")

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
