# -------------------------------------------------------------------------------------
#
#    snRNAseq fetal cell type gene overlap plot  - final
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures 4 A, B and ext_data_fig_1, 2 - Desktop

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
skene_matrix <- readRDS(paste0(DATA_DIR, 'skene_overlap_matrix.rds'))
magma_matrix <- read_table(paste0(DATA_DIR, 'magma_fetal.vs.adult_summary.tsv'))

##  Plot ------------------------------------------------------------------------------
# Figure 4 - Adult vs. fetal grid
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
  scale_y_discrete(limits = rev(levels(fetal_matrix_prep$Var2))) +
  xlab("") + 
  ylab("") +
  coord_fixed() 

fig_4B <- skene_matrix %>%
  mutate(percent = sprintf("%0.3f", percent)) %>%
  ggplot(aes(x=Var1, y=Var2, fill = 'white')) + 
  geom_tile(color = "black", size = 0.5, fill = '#DBF3FA') +
  geom_text(aes(label = paste(value, '\n', percent)), hjust = 'centre') +
  theme(axis.text.x = element_text(angle = -90, hjust = TRUE)) +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = "#000000", size = 12),
        axis.text.y  = element_text(colour = "#000000", size = 12),
        legend.position = "none",
        panel.grid = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank()) +
  scale_x_discrete(labels = as.factor(gsub("_", "-", as.vector(unique(skene_matrix$Var1))))) +
  scale_y_discrete(limits = rev(levels(skene_matrix$Var2)), labels = rev(gsub("_", "-", as.vector(unique(skene_matrix$Var2))))) +
  xlab("") + 
  ylab("") +
  coord_equal(ratio = 1) 



# Figure ext data fig 1 - Adult vs. fetal grid
magma_matrix_filt <- magma_matrix %>%
  filter(!grepl("_no_FC", VARIABLE))

# Fix y-axis labels
y_axis_labels <- rev(gsub("_", "-", as.vector(unique(magma_matrix_filt$VARIABLE))))
y_axis_labels <- gsub("-no-", " no ", y_axis_labels)

ext_data_fig_1 <- magma_matrix_filt %>%
  ggplot(aes(x = -log10(as.numeric(P)), y = factor(VARIABLE, rev(levels(factor(VARIABLE)))))) +
  geom_bar(stat = "identity", fill = c('#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD',
                                       '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD',
                                       '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD',
                                       '#CEE5FD', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF',
                                       '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF',
                                       '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#CEE5FD', '#CEE5FD',
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
  scale_y_discrete(labels = y_axis_labels) +
  xlab(expression(-log[10](P))) +
  ylab('Cell type')

# Figure 4B - fetal vs. fetal grid
magma_matrix_no_skene <- magma_matrix %>%
  filter(!grepl("skene", VARIABLE))

# Fix y-axis labels
y_axis_labels <- rev(gsub("_", "-", as.vector(unique(magma_matrix_no_skene$VARIABLE))))
y_axis_labels <- gsub("-no-", " no ", y_axis_labels)

ext_data_fig_2 <- magma_matrix_no_skene %>%
  ggplot(aes(x = -log10(as.numeric(P)), y = factor(VARIABLE, rev(levels(factor(VARIABLE)))))) +
  geom_bar(stat = "identity", fill = c('#CEE5FD', '#CEE5FD', '#CEE5FD', '#3CBB75FF', '#3CBB75FF', 
                                       '#3CBB75FF', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD'), 
           color = 'black') +
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
  scale_y_discrete(labels = y_axis_labels) +
  xlab(expression(-log[10](P))) +
  ylab('Cell type')

plot_grid(fig_4A, fig_4B )
# Save plots
# Fig 4A - Fetal cell type gene overlap plot
tiff(paste0(PLOT_DIR, "Fig_4A.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
fig_4A
dev.off()

tiff(paste0(PLOT_DIR, "Fig_4B.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
fig_4B
dev.off()

# Fig extended data figure 1
tiff(paste0(PLOT_DIR, "ext_data_fig_1.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
ext_data_fig_1
dev.off()

# Fig extended data figure 2
tiff(paste0(PLOT_DIR, "ext_data_fig_2.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
ext_data_fig_2
dev.off()

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
