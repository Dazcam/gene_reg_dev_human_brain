# -------------------------------------------------------------------------------------
#
#    snRNAseq MAGMA conditional analyses an gene overlap plots  - final
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures 6 and 7 - Desktop
#  https://stackoverflow.com/questions/55855426

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)
library(janitor)

##  Initialise variables  -------------------------------------------------------------
DATA_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/data_for_figures/"
MAGMA_DIR <- "~/Dropbox/BRAY_sc_analysis/snRNA-seq/magma_conditional_analyses/"
PLOT_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/"
SKENE_CELL_TYPES <- c("skene_CA1", "skene_InN", "skene_MSN", "skene_SS")

##  Load data -------------------------------------------------------------------------
fetal_matrix <- readRDS(paste0(DATA_DIR, 'fetal_overlap_SCZ_magma_P_0.05_matrix.rds'))
skene_matrix <- readRDS(paste0(DATA_DIR, 'skene_overlap_SCZ_magma_P_0.05_matrix.rds'))
magma_matrix <- read_table(paste0(DATA_DIR, 'magma_fetal.vs.adult_summary.tsv'))

##  Plot ------------------------------------------------------------------------------
# Figure 6A - Fetal cell type gene overlap plot
fig_6A_data <- reshape2::melt(fetal_matrix) %>%
  arrange(Var1) %>%
  group_by(Var1) %>%
  filter(row_number() >= which(Var1 == Var2)) %>%
  mutate(bold = case_when(Var1 == Var2 ~ TRUE,                  # Lines to omit
                          Var1 != Var2 ~ FALSE))                # if bold boxes
frames = fig_6A_data[fig_6A_data$bold, c("Var1", "Var2")]       # not
frames$Var1 = as.integer(frames$Var1)                           # required
frames$Var2 = rev(as.integer(frames$Var2))                      # and geom_rect

fig_6A <- fig_6A_data %>%  ggplot(aes(x = Var1, y = Var2, fill = 'white')) + 
  geom_tile(color = "black", size = 0.5, fill = '#DBF3FA') +
  geom_rect(data = frames, size = 1, fill = NA, colour = "black",
            aes(xmin = Var1 - 0.5, xmax = Var1 + 0.5, ymin = Var2 - 0.5, ymax = Var2 + 0.5)) +
  geom_text(aes(label = value, size = 12)) +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = "#000000", size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.text.y  = element_text(colour = "#000000", size = 12),
        legend.position = "none",
        panel.grid = element_blank()) +
  scale_y_discrete(limits = rev(levels(fig_6A_data$Var2))) +
  xlab("") + 
  ylab("") +
  coord_fixed() 

# Figure 7A - Adult vs. fetal gene overlap grid
fig_7A_data <- reshape2::melt(skene_matrix) %>%
  arrange(Var1) %>%
  group_by(Var1) %>%
  mutate(Var1 = R.utils::capitalize(Var1))

fig_7A <- fig_7A_data %>%  
  ggplot(aes(x=Var1, y=Var2, fill = 'white')) + 
  geom_tile(color = "black", size = 0.5, fill = '#DBF3FA') +
  geom_text(aes(label = value, size = 12)) +
  theme(axis.text.x = element_text(angle = -90, hjust = TRUE)) +
  theme_minimal() +
  theme(axis.text.x = element_text(colour = "#000000", size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.text.y  = element_text(colour = "#000000", size = 12),
        legend.position = "none",
        panel.grid = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank()) +
  scale_x_discrete(labels = as.factor(gsub("_", "-", as.vector(unique(fig_7A_data$Var1))))) +
  scale_y_discrete(limits = rev(levels(fig_7A_data$Var2)), labels = rev(gsub("_", "-", as.vector(unique(fig_7A_data$Var2))))) +
  xlab("") + 
  ylab("") +
  coord_equal(ratio = 1) 



# Conditional analysis bar chart conditioning fetal cell types on FC-ExN-2
FC_ExN_2 <- read.table(paste0(MAGMA_DIR, 'magma_all_sig_and_skene_condition_FC_ExN_2.gsa.out'), header = FALSE) %>%
  row_to_names(row_number = 1) %>% 
  filter(!str_detect(VARIABLE, "FC_ExN_2|skene")) %>%
  mutate(VARIABLE = R.utils::capitalize(VARIABLE)) %>%
  mutate(VARIABLE = gsub("_", "-", VARIABLE))

FC_ExN_2_plot <- ggplot(data = FC_ExN_2, aes(x = -log10(as.numeric(P)), y = factor(VARIABLE, rev(levels(factor(VARIABLE)))))) +
  geom_bar(stat = "identity", fill = c('#CEE5FD', '#CEE5FD', '#CEE5FD', '#3CBB75FF', '#3CBB75FF', 
                                       '#3CBB75FF', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD'), 
           color = 'black') +
  #  geom_vline(xintercept=-log10(0.05/14), linetype = "dashed", color = "black") +
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

for (CELL_TYPE in SKENE_CELL_TYPES) {
  
  ##  Load Data  ----------------------------------------------------------------------
  SKENE_DATA <- read.table(paste0(MAGMA_DIR, 'magma_all_sig_and_skene_condition_', CELL_TYPE, '.gsa.out'), header = FALSE) %>%
    row_to_names(row_number = 1) %>% 
    filter(!str_detect(VARIABLE, 'skene')) %>%
    mutate(VARIABLE = paste0(VARIABLE, " no ", CELL_TYPE)) %>%
    mutate(VARIABLE = gsub("_", "-", VARIABLE))
  
  assign(CELL_TYPE, SKENE_DATA)
  
}

fetal_vs_adult <- rbind(skene_CA1, skene_InN, skene_MSN, skene_SS) 
fetal_vs_adult$VARIABLE <- as.factor(fetal_vs_adult$VARIABLE)

fig_7B <- fetal_vs_adult %>%
  ggplot(aes(x = -log10(as.numeric(P)), y = factor(VARIABLE, rev(levels(factor(VARIABLE)))))) +
  geom_bar(stat = "identity", fill = c('#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#3CBB75FF',
                                       '#3CBB75FF', '#3CBB75FF', '#CEE5FD', '#CEE5FD', '#CEE5FD',
                                       '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                                       '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#CEE5FD', '#CEE5FD',
                                       '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD',
                                       '#CEE5FD', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#CEE5FD',
                                       '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD',
                                       '#CEE5FD', '#CEE5FD', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF',
                                       '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD'), color = 'black') +
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
  xlab(expression(-log[10](P))) +
  ylab('Cell type')

# # Conditional analysis bar charts conditioning Skene cell types
# for (CELL_TYPE in SKENE_CELL_TYPES) {
#   
#   ##  Load Data  ----------------------------------------------------------------------
#   SKENE_DATA <- read.table(paste0(MAGMA_DIR, 'magma_all_sig_and_skene_condition_', CELL_TYPE, '.gsa.out'), header = FALSE) %>%
#     row_to_names(row_number = 1) %>% 
#     filter(!str_detect(VARIABLE, 'skene')) %>%
#     #   mutate(VARIABLE = R.utils::capitalize(VARIABLE)) %>%
#     mutate(VARIABLE = gsub("_", "-", VARIABLE)) 
#   
#   SKENE_PLOT <- ggplot(data = SKENE_DATA, aes(x = -log10(as.numeric(P)), y = factor(VARIABLE, rev(levels(factor(VARIABLE)))))) +
#     geom_bar(stat = "identity", fill = c('#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#3CBB75FF', 
#                                          '#3CBB75FF', '#3CBB75FF', '#CEE5FD', '#CEE5FD', '#CEE5FD',
#                                          '#CEE5FD'), color = 'black') +
#     geom_vline(xintercept=-log10(0.05/14), linetype = "dashed", color = "black") +
#     geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
#     theme_bw() +
#     theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           panel.border = element_rect(colour = "black", size = 1),
#           plot.title = element_text(hjust = 0.5),
#           axis.title.x = element_text(colour = "#000000", size = 14),
#           axis.title.y = element_text(colour = "#000000", size = 14),
#           axis.text.x  = element_text(colour = "#000000", size = 12),
#           axis.text.y  = element_text(colour = "#000000", size = 12)) +
#     xlab(expression(-log[10](P))) +
#     ylab('Cell type') +
#     xlim(0, 11.5)
#   
#   assign(CELL_TYPE, SKENE_DATA, envir = .GlobalEnv)
#   assign(paste0(CELL_TYPE, '_plot'), SKENE_PLOT, envir = .GlobalEnv)
#   
# }
# # Fig 7B_E - SCZ
# tiff(paste0(PLOT_DIR, "Figure_7.tiff"), height = 40, width = 40, units='cm', 
#      compression = "lzw", res = 300)
# plot_grid(fig_7A, skene_CA1_plot, skene_InN_plot, skene_MSN_plot, skene_SS_plot, labels = 'AUTO', label_size = 16, align = 'h', axis = 'tb')
# dev.off()

# Save plots
# Fig 6 - SCZ
tiff(paste0(PLOT_DIR, "Figure_6.tiff"), height = 20, width = 40, units='cm', 
     compression = "lzw", res = 300)
plot_grid(fig_6A, FC_ExN_2_plot, align = 'h', labels = 'AUTO', label_size = 16, rel_heights = c(1, 0.1), axis = 'tb')
dev.off()

# Fig 7 - SCZ
tiff(paste0(PLOT_DIR, "Figure_7.tiff"), height = 20, width = 40, units='cm', 
     compression = "lzw", res = 300)
plot_grid(fig_7A, fig_7B, align = 'h', labels = 'AUTO', label_size = 16, rel_heights = c(1, 0.1), axis = 'tb')
dev.off()


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
