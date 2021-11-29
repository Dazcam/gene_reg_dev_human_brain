# -------------------------------------------------------------------------------------
#
#    snATACseq - clusters, motif enrichment and footpring plots  - final
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures 5 - Laptop

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
for (PLOT in c('A', 'B', 'C', 'D')) {
  
  plot <- readRDS(paste0(DATA_DIR, 'Fig_5', PLOT, '.rds'))
  assign(paste0('fig_5', PLOT), plot)
  
}

# Old
# for (PLOT in c('E', 'F')) {
#   
#   plot_matrix <- readRDS(paste0(DATA_DIR, 'Fig_5', PLOT, '_matrix.rds'))
#   plot_matrix_melt <- reshape::melt(plot_matrix)
#   assign(paste0('fig_5', PLOT, '_matrix'), plot_matrix_melt)
#   
# }

# Using matrices pulled directly from full motifs object
for (REGION in c('FC', 'GE')) {
  
  motifs_matrix <- readRDS(paste0(DATA_DIR, REGION, '_motifs.rds'))
  motifs_pVal_matrix <- readRDS(paste0(DATA_DIR, REGION, '_motifs_mlogPs.rds'))
  assign(paste0(REGION, '_motifs_matrix'), motifs_matrix)
  assign(paste0(REGION, '_motifs_pVal_matrix'), motifs_pVal_matrix)       
  
}

for (PLOT in c('G', 'H')) {
  
  plot_grob <- readRDS(paste0(DATA_DIR, 'Fig_5', PLOT, '_grob.rds'))
  assign(paste0('fig_5', PLOT, '_grob'), plot_grob)
  
}

# Figs 5 A-D
cer_colours <- c('#DCBEFF', '#9A6324', '#76B5C5', '#CEE5FD', '#00BDD2',  
                 '#00B6EB', '#ABDBE3', '#1E81B0', '#3CBB75FF', '#00FF00A5', 
                 '#006400', '#95D840FF', '#B7FFB7', '#10A53DFF', '#F58231', 
                 '#949494', '#CCCCCC', '#FDE725FF', '#FAA0A0', '#EF0029', 
                 '#D2042D')

fc_colours <- c('#76B5C5', '#CEE5FD', '#00B6EB', '#1E81B0', '#3CBB75FF', 
                '#00FF00A5', '#10A53DFF', '#F58231', '#CCCCCC', '#FAA0A0',
                '#D2042D')

fig_5A <- fig_5A +
  ggtitle('Frontal Cortex') +
  scale_color_manual(values = fc_colours) +
  NoLegend()


# Can get rid of numbers on plot where are they????? This changes legend labels
# levels(fig_5A[["data"]][["color"]]) <- c('ExN-5', 'ExN-4', 'ExN-3', 'ExN-2', 'InN-3', 
#                                         'InN-2', 'InN-1', 'MG', 'N-undef', 'RG-2', 'RG-1')

ge_colours <- c('#3CBB75FF', '#00FF00A5', '#10A53DFF', '#B7FFB7', '#D2042D', 
                '#EF0029', '#FAA0A0')

fig_5B <- fig_5B +
  ggtitle('Ganglionic Eminence') +
  scale_color_manual(values = ge_colours) +
  NoLegend()

fc_broad_colours <- c('#76B5C5', '#3CBB75FF', '#F58231', '#CCCCCC', '#EF0029')

fig_5C <- fig_5C +
  ggtitle('Frontal Cortex') +
  scale_color_manual(values = fc_broad_colours) +
  NoLegend()

ge_broad_colours <- c('#3CBB75FF', '#EF0029')

fig_5D <- fig_5D +
  ggtitle('Ganglionic Eminence') +
  scale_color_manual(values = ge_broad_colours) +
  NoLegend()

# Motif enrichment plots  -------------------------------------------------------------
# Had an issue with N.undef and number of motifs being plotted. Issues due to low 
# p-vals for N.undef and the pVal cut of threshold for the motif plotting function
# Just pulled the code from the repo to generate these plots manually

for (REGION in c('FC', 'GE')) {
  
  # Set variables 
  cutOff = 20
  n = 7
  pMax = Inf
  binaryClusterRows = TRUE
  clusterCols = TRUE
  pal = paletteContinuous(set = "comet", n = 100)
  rastr = TRUE
  
  # Load pVal matrix
  mat <- get(paste0(REGION, '_motifs_pVal_matrix'))
  
  # Remove all text from motif names except motif IDs
  rownames(mat) <- gsub("_.*","\\1", rownames(mat))
  
  # Remove N.undef from FC
  if (REGION == 'FC') { mat <- mat[, c(1,2,3,5)] }
  
  # Keep top x motifs
  keep <- lapply(seq_len(ncol(mat)), function(x){
    idx <- head(order(mat[, x], decreasing = TRUE), n)
    rownames(mat)[idx[which(mat[idx, x] > cutOff)]]
  }) %>% unlist %>% unique
  mat <- mat[keep, , drop = FALSE]
  
  passMat <- lapply(seq_len(nrow(mat)), function(x){
    (mat[x, ] >= 0.9*max(mat[x, ])) * 1
  }) %>% Reduce("rbind", .) %>% data.frame
  colnames(passMat) <- colnames(mat)
  
  mat[mat > pMax] <- pMax
  
  mat <- ArchR:::.rowScale(as.matrix(mat), min = 0)
  
  # This adds maxes in brackets in rownames - not needed
  # if(pMax != 100){
  #   rownames(mat[[1]]) <- paste0(rownames(mat[[1]]), " (",round(mat$max),")")
  #   rownames(passMat) <- rownames(mat[[1]])
  # }
  
  mat2 <- mat[[1]] * 100
  
  assign(paste0(REGION, '_motifs_matix_filtered'), mat2)
  
}

fig_5E <- ggplot(reshape::melt(FC_motifs_matix_filtered), aes(x = X2, y = X1, fill = value)) +
  geom_tile(color = "black", size = 0.4) +
  scale_fill_gradient2(low = "#4575b4",
                       mid = "#fdfbba",
                       high = "#d73027") +
  theme(legend.position = "none",
        legend.title=element_blank(),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 11),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_y_discrete(expand = c(0, 0)) +
  xlab(NULL) +
  ylab(NULL) 

fig_5F <- ggplot(reshape::melt(GE_motifs_matix_filtered), aes(x = X2, y = X1, fill = value)) +
  geom_tile(color = "black", size = 0.4) +
  scale_fill_gradient2(low = "#4575b4",
                       mid = "#fdfbba",
                       high = "#d73027",
                       guide = guide_colorbar(frame.colour = "black",
                                              frame.linewidth = 0.5)) +
  theme(#legend.position = "none",
    legend.title=element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_text(colour = "#000000", size = 14),
    axis.title.y = element_text(colour = "#000000", size = 14),
    axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
    axis.text.y  = element_text(colour = "#000000", size = 12),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()) +
  scale_y_discrete(expand = c(0, 0)) +
  xlab(NULL) +
  ylab(NULL) 

# Footprinting plots  -------------------------------------------------------------
fig_5G <- ggplotify::as.ggplot(fig_5G_grob[['NEUROD2_73']])
fig_5H <- ggplotify::as.ggplot(fig_5H_grob[['DLX5_412']])


# Save plots
# Fig 5A-D - Fetal cell type gene overlap plot

fig_5AtoD <- plot_grid(fig_5A, fig_5B, fig_5C, fig_5D, labels = c('A', 'B', 'C', 'D'), 
                       label_size = 16, align = c("hv"))


fig_5EtoF <- plot_grid(fig_5E, fig_5F, labels = c('E', 'F'), 
                       label_size = 16)


fig_5GtoH <- plot_grid(fig_5G, fig_5H, labels = c('G', 'H'), 
                       label_size = 16, align = c("h"))


# Fig 5
tiff(paste0(PLOT_DIR, "Fig_5.tiff"), height = 60, width = 30, units='cm', 
     compression = "lzw", res = 300)
plot_grid(fig_5AtoD, fig_5EtoF, fig_5GtoH, ncol = 1, label_size = 16)
dev.off()

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
