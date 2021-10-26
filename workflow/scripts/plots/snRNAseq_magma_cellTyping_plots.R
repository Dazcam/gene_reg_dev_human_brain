# -------------------------------------------------------------------------------------
#
#    snRNAseq MAGMA Celltyping plots  - final
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures 2, S7, S9, S11 - laptop

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)

##  Initialise variables  -------------------------------------------------------------
REGIONS <- c('cer', 'hip', 'pfc', 'tha', 'wge')
GWAS <- c( 'BPD', 'MDD', 'SCZ', 'HEIGHT')
MAGMA_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/data_for_figures/"
PLOT_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/"
NAME_BODY <- '_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN'
SUFFIX_TOP10 <- '_top10.gsa.out'

# Top 10% results
cat('\nCreating Top 10% plots ... \n')
for (REGION in REGIONS) {
  
  for (DISORDER in GWAS) {
    
    ##  Load Data  ----------------------------------------------------------------------
    #   Need to skip the first 4 columns in datafile 
    top10Data <- read.table(paste0(MAGMA_DIR, DISORDER, NAME_BODY, 
                                   '.level1.', REGION, SUFFIX_TOP10), header = FALSE)
    names(top10Data) <- as.matrix(top10Data[1, ])
    top10Data <- top10Data[-1, ]
    as.numeric(top10Data$P)
    -log10(as.numeric(top10Data$P))
    print(head(top10Data))
    
    # Remove regional info from cell types on y axis
    top10Data$VARIABLE <- gsub("^.*?\\.", "", top10Data$VARIABLE)
    
  
    ##  Plot  ---------------------------------------------------------------------------
    # Update region names in plot titles to new ones for FC and GE
    if (REGION == "pfc") {
      
      TITLE <- "FC"
      
    } else if (REGION == "wge") {
      
      TITLE <- "GE"
      
    } else {
      
      TITLE <- REGION
    }
    
    top10Plot <- ggplot(data = top10Data, aes(x = -log10(as.numeric(P)), y = factor(VARIABLE, rev(levels(factor(VARIABLE)))))) +
      geom_bar(stat = "identity", fill = "steelblue", color = 'black') +
      geom_vline(xintercept=-log10(0.00054), linetype = "dashed", color = "black") +
      theme_bw() +
      ggtitle(toupper(TITLE)) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(colour = "#000000", size = 14),
            axis.title.y = element_text(colour = "#000000", size = 14),
            axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 12)) +
      xlab(expression('-log10(P)')) +
      ylab('Cell type')
    
    assign(paste0(REGION, '_', DISORDER, '_magma_top10_plot'), top10Plot, envir = .GlobalEnv)
    assign(paste0(REGION, '_', DISORDER, '_magma_top10_data'), top10Data, envir = .GlobalEnv)
    
  }
  
}

# Create group plots
cat('\nCreating group plots ... \n')
for (DISORDER in GWAS) {
  
  magma_top10_plot <- plot_grid(get(paste0('cer_', DISORDER, '_magma_top10_plot')),
                                get(paste0('hip_', DISORDER, '_magma_top10_plot')), 
                                get(paste0('pfc_', DISORDER, '_magma_top10_plot')),
                                get(paste0('tha_', DISORDER, '_magma_top10_plot')),
                                get(paste0('wge_', DISORDER, '_magma_top10_plot')))
  
  assign(paste0('all_regions_', DISORDER, '_magma_top10_plot'), magma_top10_plot, envir = .GlobalEnv)
  
}

# Save plots
# Fig 2 - SCZ
tiff(paste0(PLOT_DIR, "Fig_2.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_SCZ_magma_top10_plot
dev.off()

# Fig S7 - SCZ
tiff(paste0(PLOT_DIR, "Fig_S7.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_BPD_magma_top10_plot
dev.off()

# Fig S9 - SCZ
tiff(paste0(PLOT_DIR, "Fig_S9.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_MDD_magma_top10_plot
dev.off()

# Fig S11 - SCZ
tiff(paste0(PLOT_DIR, "Fig_S11.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_HEIGHT_magma_top10_plot
dev.off()

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
