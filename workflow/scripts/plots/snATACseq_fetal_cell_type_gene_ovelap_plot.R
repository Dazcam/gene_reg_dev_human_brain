# -------------------------------------------------------------------------------------
#
#    snRNAseq fetal cell type gene overlap plot  - final
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures 4A - Desktop

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)


##  Initialise variables  -------------------------------------------------------------
REGIONS <- c("Cer", "FC", "GE", "Hipp", "Thal")
GWAS <- c( 'BPD', 'MDD', 'SCZ', 'HEIGHT')
DATA_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/data_for_figures/"
PLOT_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/"

##  Load data -------------------------------------------------------------------------
fetal_matrix <- readRDS(paste0(DATA_DIR, 'fetal_overlap_matrix.rds'))

##  Plot ------------------------------------------------------------------------------
fetal_plot <- reshape2::melt(fetal_matrix) %>%
  arrange(Var1) %>%
  group_by(Var1) %>%
  filter(row_number() >= which(Var1 == Var2)) %>%
  ggplot(aes(x = Var1, y = Var2, fill = 'white')) + 
  geom_tile(color = "black", size = 0.5, fill = '#DBF3FA') +
  geom_text(aes(label = value, size = 12)) +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = "#000000", size = 12),
        axis.text.y  = element_text(colour = "#000000", size = 12),
        legend.position = "none",
        panel.grid = element_blank()) +
  xlab("") + 
  ylab("") +
  coord_fixed() 

# Save plots
# Fig 4A - Fetal cell type gene overlap plot
tiff(paste0(PLOT_DIR, "Fig_4A.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
fetal_plot
dev.off()

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
