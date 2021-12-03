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

##  Info  -----------------------------------------------------------------------------

# 1. Extract fold-change enrichments and FDR for everything FDR < 0.05 
# 2. Y = cell-types, X = Processes
# 3. Circle size = fold-change enrichment, colour FDR / P-value 

##  Initialise variables  -------------------------------------------------------------
DATA_DIR <- "Dropbox/BRAY_sc_analysis/snRNA-seq/GO/"
PLOT_DIR <- "~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/"
CELL_TYPES <- c("FC_ExN_2", "FC_ExN_3", "FC_ExN_4", "FC_ExN_5", "FC_InN_1", 
                "GE_InN_1", "GE_InN_2", "Hipp_ExN_4", "Hipp_ExN_6", "Thal_ExN_5", 
                "Thal_ExN_9")

##  Load data -------------------------------------------------------------------------
for (CELL_TYPE in CELL_TYPES) {
  
  cat("\nImporting ", CELL_TYPE, " data\n")
  GO <- read_delim(paste0(DATA_DIR, CELL_TYPE, " top decile GO.txt"), 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE, progress = FALSE) %>%
    filter(FDR < 0.05) %>%
    top_n(-30)  # select the terms with the lowest FDR
    
  cat(paste0("\nNumber of terms with FDR < 0.05: ", nrow(GO), "\n")) 
    
  GO <- GO %>%  mutate(cell_type = rep(CELL_TYPE, nrow(GO))) %>%
    select(Term, FDR, `Fold Enrichment`, cell_type) 
    
  GO_plot <- ggplot(data = GO, aes(y = Term, x = cell_type, 
                          color = FDR, size = `Fold Enrichment`)) + 
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    theme_bw() + 
    ylab("") + 
    xlab("") + 
    ggtitle("GO enrichment analysis")
  
  assign(paste0(CELL_TYPE, '_GO'), GO)
  assign
  
  # Save plot
  tiff(paste0(PLOT_DIR, CELL_TYPE, "_GO.tiff"), height = 15, width = 30, units='cm', 
       compression = "lzw", res = 300)
  print(GO_plot)
  dev.off()
  
}



# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
