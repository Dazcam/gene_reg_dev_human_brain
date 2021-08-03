--------------------------------------------------------------------------------------
#
#    Create tables and plots sLDSC results - ATAC
#
#--------------------------------------------------------------------------------------

## Requirements  ----------------------------------------------------------------------

# Required on Hawk before opening R
# module load libgit2/1.1.0
# module load R/4.0.3

## Initialise R library  --------------------------------------------------------------
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )

## Load packages  ---------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)
library(rmarkdown)
library(argparser)

## Parse cell_type / set cell_type variable -------------------------------------------
cat('Parsing args ... \n')
p <- arg_parser("Read ldsc result table for snATAC-seq data ... \n")
p <- add_argument(p, "in_file", help = "No input file provided")
args <- parse_args(p)
print(args)

##  Define variables  -----------------------------------------------------------------
IN_FILE <- args$in_file


## Load Data  -------------------------------------------------------------------------
cat('Loading data ... \n')
ph_df <- read_delim('~/Desktop/prtHrt_ATACseq_SCZ_summary.tsv', "\t", 
                    escape_double = FALSE, trim_ws = TRUE, col_names = TRUE)

for (REGION in c('cer', 'fc', 'ge')) {

  ## Split into seprate dfs  ------------------------------------------------------------
  df_250bp <- ph_df %>% filter(grepl(REGION, Category)) %>%
    filter(grepl('250bp', Category)) %>%
    separate(Category, c("Cell Type", "Extra"), sep = '.250') %>%
    select(!Extra) %>%
    mutate(Zp = 2*pnorm(-abs(`Coefficient_z-score`))) %>%
    mutate(neglog10 = -log10(Zp))
  df_500bp <- ph_df %>% filter(grepl(REGION, Category)) %>%
    filter(grepl('500bp', Category)) %>%
    separate(Category, c("Cell Type", "Extra"), sep = '.500') %>%
    select(!Extra) %>%
    mutate(Zp = 2*pnorm(-abs(`Coefficient_z-score`))) %>%
    mutate(neglog10 = -log10(Zp))

  assign(paste0(REGION, '_250bp_df'), df_250bp)
  assign(paste0(REGION, '_500bp_df'), df_500bp)
  
}
  
## Plots  ---------------------------------------------------------------------------
## CER ------
# 250bp
cer_250bp_enrich_plot <- ggplot(cer_250bp_df, aes(x=`Cell Type`, y=Enrichment)) +
  geom_bar(stat = "identity", colour = "Black", fill = "#f8766d", width = 0.85) +
  scale_y_continuous(limits = c(-4, 6), expand = c(0.02, 0)) +
  theme_bw() + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), # Margin around plot
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
        axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5),
        axis.text.y  = element_text(colour="#000000", size=12)) +
  ylab("Enrichment") +
  geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                position=position_dodge(.9)) 

cer_250bp_pVal_plot <- ggplot(cer_250bp_df, aes(x=`Cell Type`, y=neglog10)) +
  geom_bar(stat = "identity", colour = "Black", fill = "#01b0f6", width = 0.85) +
  scale_y_continuous(limits = c(0, 1.8), expand = c(0.02, 0)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
        axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5),
        axis.text.y  = element_text(colour="#000000", size=12)) +
  ylab("-log10(Z-Score P)") 

cer_250bp_plot <- plot_grid(cer_250bp_enrich_plot, cer_250bp_pVal_plot, labels = c('A', 'B'), label_size = 16)

# 500bp
cer_500bp_enrich_plot <- ggplot(cer_500bp_df, aes(x=`Cell Type`, y=Enrichment)) +
  geom_bar(stat = "identity", colour = "Black", fill = "#f8766d", width = 0.85) +
  scale_y_continuous(limits = c(-2, 6), expand = c(0.02, 0)) +
  theme_bw() + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), # Margin around plot
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
        axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5),
        axis.text.y  = element_text(colour="#000000", size=12)) +
  ylab("Enrichment") +
  geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                position=position_dodge(.9)) 

cer_500bp_pVal_plot <- ggplot(cer_500bp_df, aes(x=`Cell Type`, y=neglog10)) +
  geom_bar(stat = "identity", colour = "Black", fill = "#01b0f6", width = 0.85) +
  scale_y_continuous(limits = c(0, 2.8), expand = c(0.02, 0)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
        axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5),
        axis.text.y  = element_text(colour="#000000", size=12)) +
  ylab("-log10(Z-Score P)") 

cer_500bp_plot <-plot_grid(cer_500bp_enrich_plot, cer_500bp_pVal_plot, labels = c('A', 'B'), label_size = 16)

## FC ------
# 250bp
fc_250bp_enrich_plot <- ggplot(fc_250bp_df, aes(x=`Cell Type`, y=Enrichment)) +
  geom_bar(stat = "identity", colour = "Black", fill = "#f8766d", width = 0.85) +
  scale_y_continuous(limits = c(-2, 8), expand = c(0.02, 0)) +
  theme_bw() + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), # Margin around plot
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
        axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5),
        axis.text.y  = element_text(colour="#000000", size=12)) +
  ylab("Enrichment") +
  geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                position=position_dodge(.9)) 

fc_250bp_pVal_plot <- ggplot(fc_250bp_df, aes(x=`Cell Type`, y=neglog10)) +
  geom_bar(stat = "identity", colour = "Black", fill = "#01b0f6", width = 0.85) +
  scale_y_continuous(limits = c(0, 4), expand = c(0.02, 0)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
        axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5),
        axis.text.y  = element_text(colour="#000000", size=12)) +
  ylab("-log10(Z-Score P)") 

fc_250bp_plot <-plot_grid(fc_250bp_enrich_plot, fc_250bp_pVal_plot, labels = c('A', 'B'), label_size = 16)

# 500bp
fc_500bp_enrich_plot <- ggplot(fc_500bp_df, aes(x=`Cell Type`, y=Enrichment)) +
  geom_bar(stat = "identity", colour = "Black", fill = "#f8766d", width = 0.85) +
  scale_y_continuous(limits = c(0, 8), expand = c(0.02, 0)) +
  theme_bw() + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), # Margin around plot
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
        axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5),
        axis.text.y  = element_text(colour="#000000", size=12)) +
  ylab("Enrichment") +
  geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                position=position_dodge(.9)) 

fc_500bp_pVal_plot <- ggplot(fc_500bp_df, aes(x=`Cell Type`, y=neglog10)) +
  geom_bar(stat = "identity", colour = "Black", fill = "#01b0f6", width = 0.85) +
  scale_y_continuous(limits = c(0, 6), expand = c(0.02, 0)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
        axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5),
        axis.text.y  = element_text(colour="#000000", size=12)) +
  ylab("-log10(Z-Score P)") 

fc_500bp_plot <- plot_grid(fc_500bp_enrich_plot, fc_500bp_pVal_plot, labels = c('A', 'B'), label_size = 16)


## GE ------
# 250bp
ge_250bp_enrich_plot <- ggplot(ge_250bp_df, aes(x=`Cell Type`, y=Enrichment)) +
  geom_bar(stat = "identity", colour = "Black", fill = "#f8766d", width = 0.85) +
  scale_y_continuous(limits = c(-2, 6), expand = c(0.02, 0)) +
  theme_bw() + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), # Margin around plot
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
        axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5),
        axis.text.y  = element_text(colour="#000000", size=12)) +
  ylab("Enrichment") +
  geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                position=position_dodge(.9)) 

ge_250bp_pVal_plot <- ggplot(ge_250bp_df, aes(x=`Cell Type`, y=neglog10)) +
  geom_bar(stat = "identity", colour = "Black", fill = "#01b0f6", width = 0.85) +
  scale_y_continuous(limits = c(0, 1.2), expand = c(0.02, 0)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
        axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5),
        axis.text.y  = element_text(colour="#000000", size=12)) +
  ylab("-log10(Z-Score P)") 

ge_250bp_plot <- plot_grid(ge_250bp_enrich_plot, ge_250bp_pVal_plot, labels = c('A', 'B'), label_size = 16)

# 500bp
ge_500bp_enrich_plot <- ggplot(ge_500bp_df, aes(x=`Cell Type`, y=Enrichment)) +
  geom_bar(stat = "identity", colour = "Black", fill = "#f8766d", width = 0.85) +
  scale_y_continuous(limits = c(0,4), expand = c(0.02, 0)) +
  theme_bw() + 
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), # Margin around plot
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
        axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5),
        axis.text.y  = element_text(colour="#000000", size=12)) +
  ylab("Enrichment") +
  geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error, ymax=Enrichment+Enrichment_std_error), width=.2,
                position=position_dodge(.9)) 

ge_500bp_pVal_plot <- ggplot(ge_500bp_df, aes(x=`Cell Type`, y=neglog10)) +
    geom_bar(stat = "identity", colour = "Black", fill = "#01b0f6", width = 0.85) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.02, 0)) +
    theme_bw() +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 0.8),
          axis.title.x = element_blank(),
          axis.title.y = element_text(colour="#000000", size=14, vjust = 3),
          axis.text.x  = element_text(colour="#000000", size=12, vjust = 0.5),
          axis.text.y  = element_text(colour="#000000", size=12)) +
    ylab("-log10(Z-Score P)") 
  
ge_500bp_plot <- plot_grid(ge_500bp_enrich_plot, ge_500bp_pVal_plot, labels = c('A', 'B'), label_size = 16)

# Generate markdown document
render('~/Desktop/single_cell/markdown/ldsc_report_atac.Rmd')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


