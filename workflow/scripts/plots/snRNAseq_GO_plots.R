# -------------------------------------------------------------------------------------
#
#    snRNAseq GO plots  - final
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures 6 and 7 - Desktop

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)
library(janitor)
library(readxl)

##  Info  -----------------------------------------------------------------------------

# 1. Table has 2 sets of col values - https://stackoverflow.com/questions/62613535
# 2. Y = cell-types, X = GO terms
# 3. Circle size = fold-change enrichment, colour FDR / P-value 
# Changing Fold enrichment dot radius:
# https://ggplot2-book.org/scale-other.html / https://stackoverflow.com/questions/56098080 

##  Initialise variables  -------------------------------------------------------------
DATA_DIR <- "Dropbox/BRAY_sc_analysis/snRNA-seq/GO/"
PLOT_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/"
CELL_TYPES <- c("FC_ExN_2", "FC_ExN_3", "FC_ExN_4", "FC_ExN_5", "FC_InN_1", 
                "GE_InN_1", "GE_InN_2", "Hipp_ExN_4", "Hipp_ExN_6", "Thal_ExN_5", 
                "Thal_ExN_9")


##  Load data  ------------------------------------------------------------------------
GO_DATA <- read_excel(paste0(DATA_DIR, "GOterms_to_plot_for_figure.xlsx"), skip = 2, col_names = F) %>%  
  mutate(across(everything(), ~replace_na(.x, 0))) # Get rid of NAs

colnames(GO_DATA) <- c("Term", "FC_ExN_2-FE", "FC_ExN_2-FDR", "FC_ExN_3-FE", "FC_ExN_3-FDR", 
  "FC_ExN_4-FE", "FC_ExN_4-FDR", "FC_ExN_5-FE", "FC_ExN_5-FDR", "FC_InN_1-FE", 
  "FC_InN_1-FDR", "GE_InN_1-FE", "GE_InN_1-FDR", "GE_InN_2-FE", "GE_InN_2-FDR", 
  "Hipp_ExN_4-FE", "Hipp_ExN_4-FDR", "Hipp_ExN_6-FE", "Hipp_ExN_6-FDR", "Thal_ExN_5-FE", 
  "Thal_ExN_5-FDR", "Thal_ExN_9-FE", "Thal_ExN_9-FDR")

# Create factor for y-axis order
ORDERED_LIST <- GO_DATA %>% select(Term) %>%  
#  mutate(Term = gsub(".*~", "", Term)) %>%  # Remove GO numbers
#  mutate(Term = R.utils::capitalize(Term)) %>%
  pull() %>%
  as_factor()
  
# Prep data
GO_DATA <- GO_DATA %>% 
  pivot_longer(-Term) %>% 
  separate(name, into = c("cell_type", "Score"), '-') %>% 
  pivot_wider(names_from = Score, values_from = value) %>%
  mutate_at('cell_type', str_replace_all, "_", "-") %>%
  # mutate(bio_term = gsub(".*~", "", Term)) %>%
  # mutate(bio_term = R.utils::capitalize(bio_term)) %>%
  rename(`Fold Enrichment` = FE) 

# Plot
GO_plot <- ggplot(data = GO_DATA, aes(y = factor(Term, level = rev(ORDERED_LIST)), x = cell_type, 
                                         color = -log10(FDR), size = `Fold Enrichment`)) +
  geom_point() +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        #      panel.grid.major = element_blank(), 
        #    panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 12),
        legend.text=element_text(size = 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme(axis.text.x = element_text(colour = "#000000", angle = 45, vjust = 1, hjust = 1)) +
  ylab("") + 
  xlab("") +
  scale_radius(limits = c(1, 5), range = c(1,6))

# Save plot
tiff(paste0(PLOT_DIR, "ALL_GO.tiff"), height = 30, width = 40, units='cm', 
     compression = "lzw", res = 300)
print(GO_plot)
dev.off()

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
