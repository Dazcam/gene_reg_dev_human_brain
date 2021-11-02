# -------------------------------------------------------------------------------------
#
#    snATACseq sLSDC condition plots  - final
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figure 6 - laptop

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)

##  Initialise variables  -------------------------------------------------------------
REGIONS <- c("fc", "ge")
CONDITIONS <- c("", ".vs.fc.ExN", ".vs.fc.RG", ".vs.ge.RG")
DATA_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/data_for_figures/"
PLOT_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/"

## Load Data  -------------------------------------------------------------------------

for (CONDITION in CONDITIONS) {
  
  df <- read_tsv(paste0(DATA_DIR, 'snATACseq_LDSC_summary_SCZ', CONDITION, '.tsv'))
  
  for (REGION in REGIONS) {
    
    df_region <- df %>%
      dplyr::filter(grepl(paste0(REGION, '.'), Category)) %>%
      dplyr::mutate(cell_type = gsub("^.*?\\.", "", Category))
    
    ldsc_plot <- ggplot(data = df_region, aes(x = `Coefficient_z-score`, y = factor(cell_type, rev(levels(factor(cell_type)))))) +
      geom_bar(stat = "identity", fill = "steelblue", color = 'black') +
      geom_vline(xintercept = qnorm(0.05), linetype = "dotted", color = "black") +
      geom_vline(xintercept = -qnorm(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      ggtitle(paste0(toupper(REGION), CONDITION)) +
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
      ylab('Cell type') + xlim(-2.5, 10)
    
    
    if (CONDITION == "" && REGION == 'fc') {
      
      ldsc_plot[["layers"]][[1]]$aes_params$fill <- c('#CEE5FD', '#3CBB75FF', '#F58231', '#CCCCCC', '#FF5959')
      ldsc_plot <- ldsc_plot +
        geom_vline(xintercept = qnorm(0.05/7), linetype = "dashed", color = "black") +
        geom_vline(xintercept = -qnorm(0.05/7), linetype = "dashed", color = "black")
      
    } else if (CONDITION == "" && REGION == 'ge') {
      
      ldsc_plot[["layers"]][[1]]$aes_params$fill <- c('#CEE5FD', '#FF5959') 
      ldsc_plot <- ldsc_plot +
        geom_vline(xintercept = qnorm(0.05/7), linetype = "dashed", color = "black") +
        geom_vline(xintercept = -qnorm(0.05/7), linetype = "dashed", color = "black")
      
    } else if (CONDITION == ".vs.fc.ExN" && REGION == 'fc') {
    
    ldsc_plot[["layers"]][[1]]$aes_params$fill <- c('#3CBB75FF', '#F58231', '#CCCCCC', '#FF5959')
    ldsc_plot <- ldsc_plot +
      geom_vline(xintercept = qnorm(0.05/6), linetype = "dashed", color = "black") +
      geom_vline(xintercept = -qnorm(0.05/6), linetype = "dashed", color = "black")
    
    } else if (CONDITION == ".vs.fc.RG" && REGION == 'fc') {
    
    ldsc_plot[["layers"]][[1]]$aes_params$fill <- c('#CEE5FD', '#3CBB75FF', '#F58231', '#CCCCCC')
    ldsc_plot <- ldsc_plot +
      geom_vline(xintercept = qnorm(0.05/6), linetype = "dashed", color = "black") +
      geom_vline(xintercept = -qnorm(0.05/6), linetype = "dashed", color = "black")
    
    } else if (CONDITION == ".vs.ge.RG" && REGION == 'fc') {
    
    ldsc_plot[["layers"]][[1]]$aes_params$fill <- c('#CEE5FD', '#3CBB75FF', '#F58231', '#CCCCCC', '#FF5959')
    ldsc_plot <- ldsc_plot +
      geom_vline(xintercept = qnorm(0.05/6), linetype = "dashed", color = "black") +
      geom_vline(xintercept = -qnorm(0.05/6), linetype = "dashed", color = "black")
    
    } else if (CONDITION == ".vs.fc.ExN" && REGION == 'ge') {
      
      ldsc_plot[["layers"]][[1]]$aes_params$fill <- c('#3CBB75FF', '#FF5959')
      ldsc_plot <- ldsc_plot +
        geom_vline(xintercept = qnorm(0.05/6), linetype = "dashed", color = "black") +
        geom_vline(xintercept = -qnorm(0.05/6), linetype = "dashed", color = "black")
      
    } else if (CONDITION == ".vs.fc.RG" && REGION == 'ge') {
      
      ldsc_plot[["layers"]][[1]]$aes_params$fill <- c('#3CBB75FF', '#FF5959')
      ldsc_plot <- ldsc_plot +
        geom_vline(xintercept = qnorm(0.05/6), linetype = "dashed", color = "black") +
        geom_vline(xintercept = -qnorm(0.05/6), linetype = "dashed", color = "black")
      
    } else if (CONDITION == ".vs.ge.RG" && REGION == 'ge') {
      
      ldsc_plot[["layers"]][[1]]$aes_params$fill <- '#3CBB75FF'
      ldsc_plot <- ldsc_plot +
        geom_vline(xintercept = qnorm(0.05/6), linetype = "dashed", color = "black") +
        geom_vline(xintercept = -qnorm(0.05/6), linetype = "dashed", color = "black")
      
    }
    
      assign(paste0(REGION, CONDITION, '_ldsc_plot'), ldsc_plot, envir = .GlobalEnv)
    
  }
  
}


  
ldsc_group_plot <- plot_grid(get('fc_ldsc_plot'),
                             get('ge_ldsc_plot'),
                             get('fc.vs.fc.ExN_ldsc_plot'),
                             get('ge.vs.fc.ExN_ldsc_plot'),
                             get('fc.vs.fc.RG_ldsc_plot'),
                             get('ge.vs.fc.RG_ldsc_plot'),
                             get('fc.vs.ge.RG_ldsc_plot'),
                             get('ge.vs.ge.RG_ldsc_plot'),
                             ncol = 2, labels = 'AUTO', label_size = 16) 

# Save plots
# Fig 6 
tiff(paste0(PLOT_DIR, "Fig_6.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
ldsc_group_plot
dev.off()



#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
