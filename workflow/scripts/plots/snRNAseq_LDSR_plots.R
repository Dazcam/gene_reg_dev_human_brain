# -------------------------------------------------------------------------------------
#
#    snRNAseq sLSDC plots  - final
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures 3, S8, S10, S12 - laptop

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)

##  Initialise Categorys  -------------------------------------------------------------
REGIONS <- c("Cer", "FC", "GE", "Hipp", "Thal")
GWAS <- c( 'BPD', 'MDD', 'SCZ', 'HEIGHT')
DATA_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/data_for_figures/"
PLOT_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/"

## Load Data  -------------------------------------------------------------------------


for (DISORDER in GWAS) {
  
  for (REGION in REGIONS) {
    
    # Load data
    all_top10_df <- read_tsv(paste0(DATA_DIR, 'snRNAseq_LDSC_', DISORDER, '_top10pc.tsv'))
    subset_top10_df <- as.data.frame(filter(all_top10_df, grepl(REGION, Category)))
    
    # Remove regional info from cell types on y axis
    subset_top10_df$Category <- gsub("^.*?-", "", subset_top10_df$Category)
    
    ##  Plot  ---------------------------------------------------------------------------
    # Update region names in plot titles to new ones for FC and GE
    if (REGION == "FC") {
      
      top10Plot <- ggplot(data = subset_top10_df, aes(x = `Coefficient_z-score`, y = factor(Category, rev(levels(factor(Category)))))) +
        geom_bar(stat = "identity", fill = c("#DCBEFF", '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                                             '#CEE5FD', '#CEE5FD', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', 
                                             '#3CBB75FF', '#D078FF', '#F58231', '#CCCCCC', '#FDE725FF', 
                                             '#FF5959', '#FF5959'), color = 'black') +
#        geom_vline(xintercept=qnorm(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-qnorm(0.05), linetype = "dotted", color = "black") +
        geom_vline(xintercept=-qnorm(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=qnorm(0.05), linetype = "dotted", color = "black") +
        theme_bw() +
        ggtitle('Frontal Cortex') +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", size = 1),
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_text(colour = "#000000", size = 14),
              axis.title.y = element_text(colour = "#000000", size = 14),
              axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
              axis.text.y  = element_text(colour = "#000000", size = 12)) +
        xlab(expression('Coefficient Z score')) +
        ylab('Cell type') +
        xlim(-6, 6)
      
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_plot'), top10Plot, envir = .GlobalEnv)
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_data'), subset_top10_df, envir = .GlobalEnv)
      
    } else if (REGION == "Cer") {
      
      top10Plot <- ggplot(data = subset_top10_df, aes(x = `Coefficient_z-score`, y = factor(Category, rev(levels(factor(Category)))))) +
        geom_bar(stat = "identity", fill = c('#DCBEFF', '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                                             '#CEE5FD', '#CEE5FD', '#CEE5FD', '#3CBB75FF', '#3CBB75FF', 
                                             '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#F58231', 
                                             '#CCCCCC', '#CCCCCC', '#FDE725FF', '#FF5959', '#FF5959', 
                                             '#FF5959'), color = 'black') +
#        geom_vline(xintercept=qnorm(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-qnorm(0.05), linetype = "dotted", color = "black") +
#        geom_vline(xintercept=-qnorm(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=qnorm(0.05), linetype = "dotted", color = "black") +
        theme_bw() +
        ggtitle('Cerebellum') +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", size = 1),
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_text(colour = "#000000", size = 14),
              axis.title.y = element_text(colour = "#000000", size = 14),
              axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
              axis.text.y  = element_text(colour = "#000000", size = 12)) +
        xlab(expression('Coefficient Z score')) +
        ylab('Cell type') +
        xlim(-6, 6)
      
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_plot'), top10Plot, envir = .GlobalEnv)
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_data'), subset_top10_df, envir = .GlobalEnv)
      
    } else if (REGION == "Hipp") {
      
      top10Plot <- ggplot(data = subset_top10_df, aes(x = `Coefficient_z-score`, y = factor(Category, rev(levels(factor(Category)))))) +
        geom_bar(stat = "identity", fill = c("#DCBEFF", '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD',  
                                             '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#3CBB75FF', 
                                             '#3CBB75FF', '#3CBB75FF', '#F58231', '#CCCCCC', "#CCCCCC",   
                                             '#FDE725FF', '#FF5959', '#FF5959', '#FF5959'), color = 'black') +
#        geom_vline(xintercept=qnorm(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-qnorm(0.05), linetype = "dotted", color = "black") +
 #       geom_vline(xintercept=-qnorm(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=qnorm(0.05), linetype = "dotted", color = "black") +
        theme_bw() +
        ggtitle('Hippocampus') +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", size = 1),
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_text(colour = "#000000", size = 14),
              axis.title.y = element_text(colour = "#000000", size = 14),
              axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
              axis.text.y  = element_text(colour = "#000000", size = 12)) +
        xlab(expression('Coefficient Z score')) +
        ylab('Cell type') +
        xlim(-6, 6)
      
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_plot'), top10Plot, envir = .GlobalEnv)
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_data'), subset_top10_df, envir = .GlobalEnv)
      
    } else if (REGION == "GE") {
      
      top10Plot <- ggplot(data = subset_top10_df, aes(x = `Coefficient_z-score`, y = factor(Category, rev(levels(factor(Category)))))) +
        geom_bar(stat = "identity", fill = c('#DCBEFF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF',
                                             '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#FF5959', '#FF5959', 
                                             '#FF5959'), color = 'black') +
 #       geom_vline(xintercept=qnorm(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-qnorm(0.05), linetype = "dotted", color = "black") +
 #       geom_vline(xintercept=-qnorm(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=qnorm(0.05), linetype = "dotted", color = "black") +
        theme_bw() +
        ggtitle('Ganglionic Eminence') +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", size = 1),
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_text(colour = "#000000", size = 14),
              axis.title.y = element_text(colour = "#000000", size = 14),
              axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
              axis.text.y  = element_text(colour = "#000000", size = 12)) +
        xlab(expression('Coefficient Z score')) +
        ylab('Cell type') +
        xlim(-6, 6)
      
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_plot'), top10Plot, envir = .GlobalEnv)
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_data'), subset_top10_df, envir = .GlobalEnv)
      
      
    } else {
      
      top10Plot <- ggplot(data = subset_top10_df, aes(x = `Coefficient_z-score`, y = factor(Category, rev(levels(factor(Category)))))) +
        geom_bar(stat = "identity", fill = c("#DCBEFF", '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD',  
                                             '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD',
                                             '#CEE5FD', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF',
                                             "#D078FF", '#F58231', '#FDE725FF', '#FF5959',
                                             '#FF5959', '#FF5959', '#FF5959', '#FF5959'), color = 'black') +
 #       geom_vline(xintercept=qnorm(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-qnorm(0.05), linetype = "dotted", color = "black") +
 #       geom_vline(xintercept=-qnorm(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=qnorm(0.05), linetype = "dotted", color = "black") +
        theme_bw() +
        ggtitle('Thalamus') +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", size = 1),
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_text(colour = "#000000", size = 14),
              axis.title.y = element_text(colour = "#000000", size = 14),
              axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
              axis.text.y  = element_text(colour = "#000000", size = 12)) +
        xlab(expression('Coefficient Z score')) +
        ylab('Cell type') +
        xlim(-6, 6)
      
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_plot'), top10Plot, envir = .GlobalEnv)
      assign(x = paste0(REGION, '_', DISORDER, '_ldsc_top10_df'), value = subset_top10_df, envir = .GlobalEnv)
      
    }
    
  }
  
}

cat('\nCreating group plots ... \n')
for (DISORDER in GWAS) {
  
  ldsc_top10_plot <- plot_grid(get(paste0('FC_', DISORDER, '_ldsc_top10_plot')),
                               get(paste0('GE_', DISORDER, '_ldsc_top10_plot')),
                               get(paste0('Hipp_', DISORDER, '_ldsc_top10_plot')), 
                               get(paste0('Cer_', DISORDER, '_ldsc_top10_plot')),
                               get(paste0('Thal_', DISORDER, '_ldsc_top10_plot')))
  
  assign(paste0('all_regions_', DISORDER, '_ldsc_top10_plot'), ldsc_top10_plot, envir = .GlobalEnv)
  
}

# Save plots
# Fig 3 - SCZ
tiff(paste0(PLOT_DIR, "Fig_3.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_SCZ_ldsc_top10_plot
dev.off()

# Fig S8 - BPD
tiff(paste0(PLOT_DIR, "Fig_S8.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_BPD_ldsc_top10_plot
dev.off()

# Fig S10 - MDD
tiff(paste0(PLOT_DIR, "Fig_S10.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_MDD_ldsc_top10_plot
dev.off()

# Fig 5 - HEIGHT
tiff(paste0(PLOT_DIR, "Fig_5.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_HEIGHT_ldsc_top10_plot
dev.off()


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
