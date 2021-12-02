#--------------------------------------------------------------------------------------
#
#    snRNAseq - Extract all gene expressed in significant fetal cell types
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  1. Generate average expression values for each cluster
#  2. Remove genes that are not expressed from each of the 11 cell types
#  3. Convert cell type gene IDs to entrez to match magma out put
#  4. Intersect ~18K genes in SCZ MAGMA gene file with each cell type expressed genes
#     then rank by magma gene wide p-value (most negative p-value first)
#  5. Add ranks to the list - Can't use p-values as ranks as GSEA reverses the order
#  6. Save as a .rnk file - Col 1 = Gene IDs, Col 2 = ranks
#  6. Run GSE analyses on the pre-ranked gene list page on the GSEA website:
#     https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?_GSEAPreranked_Page


##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(biomaRt)

## Set variables  ---------------------------------------------------------------------
DATA_DIR <- '~/Desktop/single_cell/scRNAseq/batch2_CR5_200121/r_objects/final/'
MAGMA_DIR <- '~/Dropbox/BRAY_sc_analysis/snRNA-seq/magma_conditional_analyses/'
REGIONS <- c('cer', 'hip', 'pfc', 'tha', 'wge')
CELL_TYPES <- c('FC_ExN_2', 'FC_ExN_3', 'FC_ExN_4', 'FC_ExN_5', 'FC_InN_1', 
                'GE_InN_1', 'GE_InN_2', 'Hipp_ExN_4', 'Hipp_ExN_6', 'Thal_ExN_5', 
                'Thal_ExN_9')


## Load Data --------------------------------------------------------------------------
# Seurat objects and calulate average expression for each region
for (REGION in REGIONS) { 
  
  seurat.obj <- readRDS(paste0(DATA_DIR, 'seurat.', REGION, '.final.rds'))
  seurat.ave_exp <- AverageExpression(object = seurat.obj)
  
  assign(paste0('seurat.', REGION), seurat.obj, .GlobalEnv)
  assign(paste0(REGION, '_avExp'), seurat.ave_exp, .GlobalEnv)
  
}

magma_df <- read_table2("~/Dropbox/BRAY_sc_analysis/snRNA-seq/magma_conditional_analyses/SCZ_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.out")

## Extract gene list for each cell type   ---------------------------------------------
FC_ExN_2_gene_list <- as.data.frame(pfc_avExp$RNA) %>%
  dplyr::select(`FC-ExN-2`) %>%
  filter(`FC-ExN-2` > 0) %>%
  arrange(desc(`FC-ExN-2`))

FC_ExN_3_gene_list <- as.data.frame(pfc_avExp$RNA) %>%
  dplyr::select(`FC-ExN-3`) %>%
  filter(`FC-ExN-3` > 0) %>%
  arrange

FC_ExN_4_gene_list <- as.data.frame(pfc_avExp$RNA) %>%
  dplyr::select(`FC-ExN-4`) %>%
  filter(`FC-ExN-4` > 0) %>%
  arrange(desc(`FC-ExN-4`))

FC_ExN_5_gene_list <- as.data.frame(pfc_avExp$RNA) %>%
  dplyr::select(`FC-ExN-5`) %>%
  filter(`FC-ExN-5` > 0) %>%
  arrange(desc(`FC-ExN-5`))

FC_InN_1_gene_list <- as.data.frame(pfc_avExp$RNA) %>%
  dplyr::select(`FC-InN-1`) %>%
  filter(`FC-InN-1` > 0) %>%
  arrange(desc(`FC-InN-1`))

GE_InN_1_gene_list <- as.data.frame(wge_avExp$RNA) %>%
  dplyr::select(`GE-InN-1`) %>%
  filter(`GE-InN-1` > 0) %>%
  arrange(desc(`GE-InN-1`))

GE_InN_2_gene_list <- as.data.frame(wge_avExp$RNA) %>%
  dplyr::select(`GE-InN-2`) %>%
  filter(`GE-InN-2` > 0) %>%
  arrange(desc(`GE-InN-2`))

Hipp_ExN_4_gene_list <- as.data.frame(hip_avExp$RNA) %>%
  dplyr::select(`Hipp-ExN-4`) %>%
  filter(`Hipp-ExN-4` > 0) %>%
  arrange(desc(`Hipp-ExN-4`))

Hipp_ExN_6_gene_list <- as.data.frame(hip_avExp$RNA) %>%
  dplyr::select(`Hipp-ExN-6`) %>%
  filter(`Hipp-ExN-6` > 0) %>%
  arrange(desc(`Hipp-ExN-6`))

Thal_ExN_5_gene_list <- as.data.frame(tha_avExp$RNA) %>%
  dplyr::select(`Thal-ExN-5`) %>%
  filter(`Thal-ExN-5` > 0) %>%
  arrange(desc(`Thal-ExN-5`))

Thal_ExN_9_gene_list <- as.data.frame(tha_avExp$RNA) %>%
  dplyr::select(`Thal-ExN-9`) %>%
  filter(`Thal-ExN-9` > 0) %>%
  arrange(desc(`Thal-ExN-9`))

## Convert cell type gene IDs to Entrez   ---------------------------------------------
for (CELL_TYPE in CELL_TYPES) {
  
  GENE_FILE <- get(paste0(CELL_TYPE, '_gene_list'))
  GENE_LIST <- rownames(GENE_FILE)
  
  
  cat(paste0('\n\nRunning conversion for: ', CELL_TYPE, '\n'))
  cat(paste0('Total genes before converting IDs: ', nrow(GENE_LIST), '\n'))
  
  # Convert gene IDs 
  cat('Converting human gene IDs using BiomaRt ... \n')
  mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ENTREZ_GENES = getBM(attributes = c("hgnc_symbol", "entrezgene_id"), 
                       filters = "hgnc_symbol", 
                       values = GENE_LIST, 
                       bmHeader = T, 
                       mart = mart)
  ENTREZ_GENES_noNA <- ENTREZ_GENES %>% drop_na()
  ENTREZ_GENES_UNIQUE <- as.data.frame(unique(ENTREZ_GENES_noNA[, 2]))
  colnames(ENTREZ_GENES_UNIQUE) <- CELL_TYPE
  
  cat(paste0('Total genes after converting IDs: ', nrow(ENTREZ_GENES), '\n'))
  cat(paste0('Total genes after removing NAs: ', nrow(ENTREZ_GENES_noNA), '\n'))
  cat(paste0('Total genes after checking uniqueness (so final count): ', nrow(ENTREZ_GENES_UNIQUE), '\n'))
  
  assign(paste0(CELL_TYPE, '_entrezID'), ENTREZ_GENES_UNIQUE)
  
}

## Intersect with genes in SCZ MAGMA genewise p-value file  ---------------------------
for (CELL_TYPE in CELL_TYPES) {
  
  GENE_LIST <- get(paste0(CELL_TYPE, '_entrezID'))
  
  cat(paste0('\n\nRunning overlap analysis for: ', CELL_TYPE, '\n'))
  cat(paste0('Total genes before overlap: ', nrow(GENE_LIST), '\n'))
  
  GENE_OVERLAP <- magma_df %>% 
    filter(GENE %in% GENE_LIST[,1]) %>%
    arrange(P) %>%
    dplyr::select(GENE) # Can't use p-values as ranks as GSEA reverses the order
  
  # Add ranks
  GENE_OVERLAP$ranks <- seq(1, nrow(GENE_OVERLAP), 1)
  
  cat(paste0('Total genes after overlapping gene sets: ', nrow(GENE_OVERLAP), '\n'))
  
  assign(paste0(CELL_TYPE, '_overlap'), GENE_OVERLAP)
  
  write.table(GENE_OVERLAP, paste0(MAGMA_DIR, CELL_TYPE, '_all_genesExp_in_magmaSCZ.rnk'),
              quote = FALSE, row.names = FALSE,  sep = '\t')
  
}


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
