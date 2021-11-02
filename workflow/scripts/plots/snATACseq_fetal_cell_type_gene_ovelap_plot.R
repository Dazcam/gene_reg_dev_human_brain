# -------------------------------------------------------------------------------------
#
#    snRNAseq fetal cell type gene overlap plot  - final
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures 4A and 4B - Desktop

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
magma_matrix <- read_table(paste0(DATA_DIR, 'magma_fetal.vs.adult_summary.tsv'))

##  Plot ------------------------------------------------------------------------------
fig_4A <- reshape2::melt(fetal_matrix) %>%
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


fig_4B <- ggplot(data = magma_matrix, aes(x = -log10(as.numeric(P)), y = factor(VARIABLE, rev(levels(factor(VARIABLE)))))) +
  geom_bar(stat = "identity", fill = c('#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD',
                                       '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD',
                                       '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD',
                                       '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#3CBB75FF',
                                       '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF',
                                       '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF',
                                       '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#CEE5FD',
                                       '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD',
                                       '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD',
                                       '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD',
                                       '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD') , color = 'black') +
  geom_vline(xintercept=-log10(0.05/54), linetype = "dashed", color = "black") +
  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme_bw() +
#  ggtitle(toupper(TITLE)) +
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

plot_grid(fig_4A, fig_4B, labels = 'AUTO', label_size = 16)

# Save plots
# Fig 4A - Fetal cell type gene overlap plot
tiff(paste0(PLOT_DIR, "Fig_4A.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
fig_4A
dev.off()

# Fig 4B - Fetal cell type gene overlap plot
tiff(paste0(PLOT_DIR, "Fig_4B.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
fig_4B
dev.off()

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
