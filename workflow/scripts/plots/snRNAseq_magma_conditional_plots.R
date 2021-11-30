# -------------------------------------------------------------------------------------
#
#    snRNAseq MAGMA Celltyping plots conditional analyses  - final
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures 2, 4B, S7, S9, S11 - laptop

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)
library(janitor)

##  Initialise variables  -------------------------------------------------------------
MAGMA_DIR <- "~/Dropbox/BRAY_sc_analysis/snRNA-seq/magma_conditional_analyses/"
PLOT_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/"
SKENE_CELL_TYPES <- c("skene_CA1", "skene_InN", "skene_MSN", "skene_SS")

# Condition on FC-ExN-2
FC_ExN_2 <- read.table(paste0(MAGMA_DIR, 'magma_all_sig_and_skene_condition_FC_ExN_2.gsa.out'), header = FALSE) %>%
  row_to_names(row_number = 1) %>% 
  filter(!str_detect(VARIABLE, "FC_ExN_2|skene")) %>%
  mutate(VARIABLE = R.utils::capitalize(VARIABLE)) %>%
  mutate(VARIABLE = gsub("_", "-", VARIABLE))

FC_ExN_2_plot <- ggplot(data = FC_ExN_2, aes(x = -log10(as.numeric(P)), y = factor(VARIABLE, rev(levels(factor(VARIABLE)))))) +
  geom_bar(stat = "identity", fill = c('#CEE5FD', '#CEE5FD', '#CEE5FD', '#3CBB75FF', '#3CBB75FF', 
                                       '#3CBB75FF', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD'), 
           color = 'black') +
  geom_vline(xintercept=-log10(0.05/14), linetype = "dashed", color = "black") +
  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 12)) +
  xlab(expression(-log[10](P))) +
  ylab('Cell type') +
  xlim(0, 11.5)

# Condition on Skene cell types
for (CELL_TYPE in SKENE_CELL_TYPES) {
  
  ##  Load Data  ----------------------------------------------------------------------
  SKENE_DATA <- read.table(paste0(MAGMA_DIR, 'magma_all_sig_and_skene_condition_', CELL_TYPE, '.gsa.out'), header = FALSE) %>%
    row_to_names(row_number = 1) %>% 
    filter(!str_detect(VARIABLE, 'skene')) %>%
    #   mutate(VARIABLE = R.utils::capitalize(VARIABLE)) %>%
    mutate(VARIABLE = gsub("_", "-", VARIABLE)) 
  
  SKENE_PLOT <- ggplot(data = SKENE_DATA, aes(x = -log10(as.numeric(P)), y = factor(VARIABLE, rev(levels(factor(VARIABLE)))))) +
    geom_bar(stat = "identity", fill = c('#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#3CBB75FF', 
                                         '#3CBB75FF', '#3CBB75FF', '#CEE5FD', '#CEE5FD', '#CEE5FD',
                                         '#CEE5FD'), color = 'black') +
    geom_vline(xintercept=-log10(0.05/14), linetype = "dashed", color = "black") +
    geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
    theme_bw() +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 1),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(colour = "#000000", size = 14),
          axis.title.y = element_text(colour = "#000000", size = 14),
          axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
          axis.text.y  = element_text(colour = "#000000", size = 12)) +
    xlab(expression(-log[10](P))) +
    ylab('Cell type') +
    xlim(0, 11.5)
  
  assign(CELL_TYPE, SKENE_DATA, envir = .GlobalEnv)
  assign(paste0(CELL_TYPE, '_plot'), SKENE_PLOT, envir = .GlobalEnv)
  
}

# Save plots
# Fig 6B - SCZ
tiff(paste0(PLOT_DIR, "Figure_6B.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
FC_ExN_2_plot
dev.off()

# Fig 7B_E - SCZ
tiff(paste0(PLOT_DIR, "Figure_7B_E.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
plot_grid(skene_CA1_plot, skene_InN_plot, skene_MSN_plot, skene_CA1_plot, 
          labels = 'AUTO', label_size = 16)
dev.off()


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
