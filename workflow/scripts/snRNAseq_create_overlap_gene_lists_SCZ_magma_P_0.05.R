#--------------------------------------------------------------------------------------
#
#    snRNAseq - Fetal / Adult gene overlaps (SCZ MAGMA P < 0.05)
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  1. Determine number of SCZ MAGMA P < 0.05 in q10 of each fetal and adult cell type 
#  2. Get no. of gene overlaps between all fetal cell types (using gene lists created at stage 1)
#  3. Get no. of gene overlaps between all fetal and adult cell types (using gene lists created at stage 1)
#  4. Save lists of q10 / SCZ MAGMA < 0.05 genes for each cell type to run GO analysis

# Note: gene lists were created in snRNAseq_magma_conditonal_analysis_file_prep.R
# intersect matrix code from: snRNAseq_fetal_vs_adult_create_overlap_gene_lists.R

##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(biomaRt)

## Set variables  ---------------------------------------------------------------------
DATA_DIR <- '~/Desktop/single_cell/scRNAseq/batch2_CR5_200121/r_objects/final/'
MAGMA_DIR <- '~/Dropbox/BRAY_sc_analysis/snRNA-seq/magma_conditional_analyses/'
ENTREZ_DIR <- '~/Dropbox/BRAY_sc_analysis/snRNA-seq/top_decile_gene_lists/Entrez/'
MATRIX_DIR <- '~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/data_for_figures/'
CELL_TYPES <- c('FC_ExN_2', 'FC_ExN_3', 'FC_ExN_4', 'FC_ExN_5', 'FC_InN_1', 
                'GE_InN_1', 'GE_InN_2', 'Hipp_ExN_4', 'Hipp_ExN_6', 'Thal_ExN_5', 
                'Thal_ExN_9', 'skene_InN', 'skene_MSN', 'skene_CA1', 'skene_SS')
FETAL_CELL_TYPES <- c('FC_ExN_2', 'FC_ExN_3', 'FC_ExN_4', 'FC_ExN_5', 'FC_InN_1', 
                      'GE_InN_1', 'GE_InN_2', 'Hipp_ExN_4', 'Hipp_ExN_6', 'Thal_ExN_5', 
                      'Thal_ExN_9')


## Load Data --------------------------------------------------------------------------
# SCZ magma genes at P < 0.05
magma_df <- read_table(paste0(MAGMA_DIR, "SCZ_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.out")) %>%
  filter(P < 0.05)

cat('\nTotal genes in SCZ MAGMA: 18904')
cat('\nTotal genes in SCZ MAGMA P < 0.05: ', nrow(magma_df), '\n')

# Fetal and skene q10 cell type genes that overlap SCZ MAGMA P < 0.05 genes
for (CELL_TYPE in CELL_TYPES) {
  
  GENE_LIST <- read_tsv(paste0(ENTREZ_DIR, CELL_TYPE, '_entrez_gene_list.tsv'), show_col_types = FALSE) %>%
    pull(CELL_TYPE)
  
  GENE_OVERLAPS <- filter(magma_df, GENE %in% GENE_LIST) %>%
    pull(GENE)
  
  
  assign(paste0(CELL_TYPE, '_entrezID'), GENE_LIST)
  assign(paste0(CELL_TYPE, '_overlaps'), GENE_OVERLAPS)
  
  cat(paste0('\nGenes in ', CELL_TYPE, ' q10: ', length(GENE_LIST)))
  cat(paste0('\nGenes in ', CELL_TYPE, ' q10 overlapping SCZ MAGMA P < 0.05: ', length(GENE_OVERLAPS), '\n'))
  
}

# Calculate overlaps
## Create overlap grids for comaprisons across all cell-types  ------------------------
# Create overlap grids
cat('\nCreating overlap grids for comparisons across all cell-types  ... \n')
fetal_overlap_list <- list(FC_ExN_2_overlaps, FC_ExN_3_overlaps, FC_ExN_4_overlaps, 
                           FC_ExN_5_overlaps, FC_InN_1_overlaps, GE_InN_1_overlaps, 
                           GE_InN_2_overlaps, Hipp_ExN_4_overlaps, Hipp_ExN_6_overlaps,
                           Thal_ExN_5_overlaps, Thal_ExN_9_overlaps)
skene_overlap_list <- list(skene_SS_overlaps, skene_InN_overlaps, skene_MSN_overlaps, skene_CA1_overlaps)

## Intra-dataset comparisons - all fetal cells  ---------------------------------------
cat('Intra-dataset comparisons - all fetal cells  ... \n')

# Initialise vector
vect = vector()

cat('Creating pairwise intersections  ... \n')
# Get pairwise intersections - https://stackoverflow.com/questions/51697101
for (i in fetal_overlap_list) {
  for (j in fetal_overlap_list) {
    vect <- c(vect, length(intersect(i,j)))}}

# Create martix
fetal_matrix <- matrix(vect, nrow=length(fetal_overlap_list))

# Assign columns and rows names
colnames(fetal_matrix) <- c('FC-ExN-2', 'FC-ExN-3', 'FC-ExN-4', 'FC-ExN-5', 
                            'FC-InN-1', 'GE-InN-1', 'GE-InN-2', 'Hipp-ExN-4', 
                            'Hipp-ExN-6', 'Thal-ExN-5', 'Thal-ExN-9')
rownames(fetal_matrix) <- c('FC-ExN-2', 'FC-ExN-3', 'FC-ExN-4', 'FC-ExN-5', 
                            'FC-InN-1', 'GE-InN-1', 'GE-InN-2', 'Hipp-ExN-4', 
                            'Hipp-ExN-6', 'Thal-ExN-5', 'Thal-ExN-9')

## Inter-dataset comparisons - 4 Skene vs. all fetal cells  ---------------------------
cat('\nInter-dataset comparisons - 4 Skene vs. all fetal cells  ... \n')
# Get intersections all skene cell types in all fetal cell types

cat('Creating pairwise intersections  ... \n')
x <- 0
for (i in skene_overlap_list) {
  
  x <- x + 1   
  vect_skene <- vector()
  
  for (j in fetal_overlap_list) {
    vect_skene <- c(vect_skene, length(intersect(i,j)))
  }
  
  assign(paste0('vect_skene', x), vect_skene)
  
}

# Create DF of overlap counts
skene_matrix <- rbind(vect_skene1, vect_skene2, vect_skene3, vect_skene4)

# Assign columns and rows names
colnames(skene_matrix) <- c('FC_ExN_2', 'FC_ExN_3', 'FC_ExN_4', 'FC_ExN_5', 'FC_InN_1', 'GE_InN_1', 
                        'GE_InN_2', 'Hipp_ExN_4', 'Hipp_ExN_6', 'Thal_ExN_5', 'Thal_ExN_9')
rownames(skene_matrix) <- c('skene_SS', 'skene_InN', 'skene_MSN', 'skene_CA1')

##  Save files  -----------------------------------------------------------------------
saveRDS(fetal_matrix, paste0(MATRIX_DIR, 'fetal_overlap_SCZ_magma_P_0.05_matrix.rds'))
saveRDS(skene_matrix, paste0(MATRIX_DIR, 'skene_overlap_SCZ_magma_P_0.05_matrix.rds'))

for (CELL_TYPE in FETAL_CELL_TYPES) {
  
  GENE_LIST <- as.data.frame(get(paste0(CELL_TYPE, '_overlaps')))
  write_tsv(GENE_LIST, paste0(ENTREZ_DIR, CELL_TYPE, '_overlap_SCZ_magma_P_0.05_genes.tsv'), 
            col_names = FALSE)
  
}

x <- read_table(paste0(MAGMA_DIR, "SCZ_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.out"))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
