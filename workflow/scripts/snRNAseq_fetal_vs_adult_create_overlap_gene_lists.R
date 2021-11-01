#--------------------------------------------------------------------------------------
#
#    Fetal vs. Adult - Find overlapping genes & create condtional gene lists 
#    
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(ggVennDiagram)
library(readxl)
library(biomaRt)
library(cowplot)
library(org.Hs.eg.db)

# -------------------------------------------------------------------------------------
#                 
#              SECTION 1 - Find overlapping genes
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  1. Get genes in the top decile of each cell type
#  2. Change Skene mouse gene IDs to human using BiomaRt
#  3. Get Magma gene-wide P-values 
#  4. Take top 1000 SCZ genes in terms of Magma gene-wide P-values
#  5. Convert SCZ gene IDs from Entrez to HGNC 
#  5. Find out how many top 1K SCZ genes are in the top decile of each cell type 
#  6. Create a Venn Diagram to check overlaps between various cell types
#  7. Creat overlap plots for multiple comparisons  
#  https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/
#  7. https://github.com/gaospecial/ggVennDiagram

## Inputs  ----------------------------------------------------------------------------

# 1. Q10 gene lists for all fetal cell types
# 2. Skene specificity scores xslx for there 4 sig. cell types
# 3. Magma SCZ gene-wide p-values

## Intra-dataset gene lists  ----------------------------------------------------------
cat('\n\n\n|-------------------------------------------|\n')
cat('\n               SECTION 1 \n')
cat('\n|-------------------------------------------|\n\n')

##  Set variables  --------------------------------------------------------------------
cat('\nDefining variables for section 1 ... \n')
DATA_DIR <- '~/Desktop/single_nuclei/snRNAseq/tests/'
DISORDER <- 'SCZ'
CELL_TYPES <- c('FC-ExN-2', 'FC-ExN-3', 'FC-ExN-4', 'FC-ExN-5', 'FC-InN-1', 
                'GE-InN-1', 'GE-InN-2', 'Hipp-ExN-4', 'Hipp-ExN-6', 'Thal-ExN-5',
                'Thal-ExN-9')
ALL_CELL_TYPES <- c('FC-ExN-2', 'FC-ExN-3', 'FC-ExN-4', 'FC-ExN-5', 'FC-InN-1', 
                    'GE-InN-1', 'GE-InN-2', 'Hipp-ExN-4', 'Hipp-ExN-6', 'skene_InN', 
                    'skene_MSN', 'skene_CA1', 'skene_SS', 'Thal-ExN-5',
                    'Thal-ExN-9')
GO_OUTDIR <- '~/Dropbox/BRAY_sc_analysis/files_for_paper/RNA/fetal_vs_adult/GO_analysis/'

## Load and prep data  ----------------------------------------------------------------
# Gene lists from Q10s for our significant cell types in snRNAseq data
cat('Load the fetal brain cell type Q10 gene lists ... \n')
for (CELL_TYPE in CELL_TYPES) {

  GENE_LIST <- read_csv(paste0(DATA_DIR, CELL_TYPE, "_Q10_genes.tsv"), col_names = FALSE)
  CELL_TYPE <- gsub("-", "_", CELL_TYPE)
  assign(CELL_TYPE, GENE_LIST)
                        
}

# Load Skene gene specificity data
cat('Load the Skene cell specificity scores ... \n')
skene_genes <- read_xlsx('~/Desktop/single_nuclei/snRNAseq/tests/skene_gene_specificity_scores.xlsx')
colnames(skene_genes) <- c('gene', 'skene_InN', 'skene_MSN', 'skene_CA1', 'skene_SS')

# Extract top 10% most sig genes in each Skene cell-type and convert gene IDs to human
for (CELL_TYPE in c('skene_InN', 'skene_MSN', 'skene_CA1', 'skene_SS')) {
  
  # Pull out Q10 vector
  cat(paste0('Extracting genes that have highest expression specificity in: ', CELL_TYPE, '\n'))
  GENES <- skene_genes %>%
  dplyr::select(gene, .data[[CELL_TYPE]]) %>%
  arrange(desc(.data[[CELL_TYPE]])) %>%
  top_frac(.10) %>%
  pull(gene) 
  
  cat(paste0('Total genes before converting IDs: ', length(GENES), '\n'))

  # Convert mouse gene IDs to human
  cat('Convert gene IDs from mouse to human using BiomaRt ... \n')
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), 
                   filters = "mgi_symbol", 
                   values = GENES, 
                   mart = mouse, 
                   attributesL = c("hgnc_symbol"), 
                   martL = human, 
                   uniqueRows=T)
  humanx <- as.data.frame(unique(genesV2[, 2]))
  colnames(humanx) <- "X1"
  
  cat(paste0('Total genes after converting IDs: ', dim(humanx)[1], '\n'))
  
  assign(CELL_TYPE, humanx)

}

# Magma gene wide p-values
cat('Load Magma SCZ gene-wide p-values ... \n')
gene_wide_df <- read_table2(paste0(DATA_DIR, DISORDER, 
                                   "_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.out"))
top_1k_genes_df <- gene_wide_df %>%
  arrange(P) 
top_1k_genes_df <- top_1k_genes_df[1:1000,]
top_1k_genes_entrez <- as.vector(as.character(top_1k_genes_df$GENE))

# Convert Magma ENTREZ IDs to HGNC 
cat('Convert SCZ gene IDs from Entrez to HGNC ... \n')
my.symbols <- c(top_1k_genes_entrez)
hs <- org.Hs.eg.db
gene_IDs <- select(hs,
                   keys = my.symbols,
                   columns = c("ENTREZID", "SYMBOL"),
                   keytype = "ENTREZID")
top_1k_genes_hgnc <- gene_IDs$SYMBOL


## Get highest expressing genes in cell type that overlap top GWAS genes  -------------
for (CELL_TYPE in ALL_CELL_TYPES) {
  
  cat(paste0('\nExtracting genes in top ', DISORDER, 
             ' genes that have highest expression specificity in: ', CELL_TYPE, '\n'))
  CELL_TYPE <- gsub("-", "_", CELL_TYPE)
  TOP_10_GENES <- as.vector(get(CELL_TYPE))
  OVERLAPPING_GENES <- unique(grep(paste(top_1k_genes_hgnc, collapse="|"), 
                                   as.vector(TOP_10_GENES$X1), value = TRUE))
  OVERLAPPING_GENES_COLLAPSE <- paste(OVERLAPPING_GENES, collapse="|")
  assign(paste0(CELL_TYPE, '_overlaps'), OVERLAPPING_GENES)
  
  cat(paste0('Total overlapping genes: ', length(OVERLAPPING_GENES), '\n'))
  cat(paste0('\n', OVERLAPPING_GENES_COLLAPSE, '\n'))
  
  # Save overlap genes to file for GO analysis
  write_tsv(base::as.data.frame(OVERLAPPING_GENES),
            paste0(GO_OUTDIR, CELL_TYPE, '_overlaps_Q10_vs_magmaTop1000.tsv'),
            col_names = FALSE)
  
}

## Create Venn Diagrams of overlaps of top GWAS genes across cell-types  --------------
# Intra-dataset comparisons  ----------------------------------------------------------
# Fig4a - FC-ExN-2, FC-ExN-3, FC-ExN-4, FC-ExN-5, FC-InN-1
cat('\nCreating Venn gene lists ... \n')
gene_list_4a <- list(
  FC_ExN_2 = FC_ExN_2_overlaps, 
  FC_ExN_3 = FC_ExN_3_overlaps, 
  FC_ExN_4 = FC_ExN_4_overlaps,
  FC_ExN_5 = FC_ExN_5_overlaps,
  FC_InN_1 = FC_InN_1_overlaps
)

## Inter-dataset comparisons  ---------------------------------------------------------
# Fig4b - FC-ExN-(2-5) with Skene's pyramidal SS
gene_list_4b <- list(
  FC_ExN_2 = FC_ExN_2_overlaps, 
  FC_ExN_3 = FC_ExN_3_overlaps, 
  FC_ExN_4 = FC_ExN_4_overlaps,
  FC_ExN_5 = FC_ExN_5_overlaps,
  Skene_SS = skene_SS_overlaps
)

# Fig4c - FC-InN-1, GE-InN-1, GE-InN-2, skene_InN
gene_list_4c <- list(
  FC_InN_1 = FC_InN_1_overlaps, 
  GE_InN_1 = GE_InN_1_overlaps, 
  GE_InN_2 = GE_InN_2_overlaps,
  Skene_InN = skene_InN_overlaps
)

# Fig4d - FC-InN-1, GE-InN-1, GE-InN-2, skene_MSN
gene_list_4d <- list(
  FC_InN_1 = FC_InN_1_overlaps, 
  GE_InN_1 = GE_InN_1_overlaps, 
  GE_InN_2 = GE_InN_2_overlaps,
  Skene_MSN = skene_MSN_overlaps
)

# Fig4d - Hipp-ExN4 and Hipp-ExN-6, skene_CA1
gene_list_4e <- list(
  Hipp_ExN_4 = Hipp_ExN_4_overlaps, 
  Hipp_ExN_6 = Hipp_ExN_4_overlaps, 
  Skene_CA1 = skene_CA1_overlaps
)

# Plot
cat('Plotting Venn diagrams ... \n')
fig_4a <- ggVennDiagram(gene_list_4a, label_alpha = 0, edge_size = 0.5) +
  ggplot2::scale_fill_gradient(low="#FEE2E8", high = "#FF708E")
fig_4b <- ggVennDiagram(gene_list_4b, label_alpha = 0, edge_size = 0.5) +
  ggplot2::scale_fill_gradient(low="#FEE2E8", high = "#FF708E")
fig_4c <- ggVennDiagram(gene_list_4c, label_alpha = 0, edge_size = 0.5) +
  ggplot2::scale_fill_gradient(low="#FEE2E8", high = "#FF708E") +
  fig_4d <- ggVennDiagram(gene_list_4d, label_alpha = 0, edge_size = 0.5) +
  ggplot2::scale_fill_gradient(low="#FEE2E8", high = "#FF708E")
fig_4e <- ggVennDiagram(gene_list_4e, label_alpha = 0, edge_size = 0.5) +
  ggplot2::scale_fill_gradient(low="#FEE2E8", high = "#FF708E") 

# Intersections
cat('Producing intersection lists ... \n')
fig_4a_intersection <- process_region_data(Venn(gene_list_4a))
fig_4b_intersection <- process_region_data(Venn(gene_list_4b))
fig_4c_intersection <- process_region_data(Venn(gene_list_4c))
fig_4d_intersection <- process_region_data(Venn(gene_list_4d))
fig_4e_intersection <- process_region_data(Venn(gene_list_4e))

# Data info
process_data(Venn(gene_list_4e))

# Save plots
#for (PLOT in c('a', 'b', 'c', 'd', 'e')) {
#  
#  FIG <- get(paste0('fig_4', PLOT))
#  tiff(paste0(DATA_DIR, 'fig_4', PLOT, '.tiff'), height = 30, width = 20, units='cm', 
#       compression = 'lzw', res = 300)
#  print(FIG)
#  dev.off()
  
#}

# Output tables
#top_GWAS_genes <- cbind(top_1k_genes_df, top_1k_genes_hgnc)

#write_tsv(top_GWAS_genes, paste0(DATA_DIR, 'top_', DISORDER, '_genes.tsv'))

#for (CELL_TYPE in CELL_TYPES) {
  
#  cat(paste0('\nWriting tables for ', DISORDER, '\n'))
#  CELL_TYPE <- gsub("-", "_", CELL_TYPE)
#  GENES <- as.data.frame(get(paste0(CELL_TYPE, '_overlaps')))
#  write_tsv(GENES, paste0(DATA_DIR, CELL_TYPE, '_overlapping_', 
#                          DISORDER, '_genes.tsv'), col_names = FALSE)


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
skene_df <- rbind(vect_skene1, vect_skene2, vect_skene3, vect_skene4)

# Assign columns and rows names
colnames(skene_df) <- c('FC_ExN_2', 'FC_ExN_3', 'FC_ExN_4', 'FC_ExN_5', 'FC_InN_1', 'GE_InN_1', 
                        'GE_InN_2', 'Hipp_ExN_4', 'Hipp_ExN_6', 'Thal_ExN_5', 'Thal_ExN_9')
rownames(skene_df) <- c('skene_SS', 'skene_InN', 'skene_MSN', 'skene_CA1')

cat('Calculate percentages for Skene vs. fetal comparisons  ... \n')
# Initialise vector to add percentages to grid
prcnt <- vector()

# Calculate percentages
skene_df <- reshape2::melt(skene_df) %>%
  arrange(Var1) %>%
  group_by(Var1) 
  
for (i in 1:length(skene_df$Var2)) {
    
  GENE_NUMBER <- length(get(paste0(as.character(skene_df$Var2[i]), '_overlaps')))
  prcnt <- c(prcnt, (skene_df$value[i] / GENE_NUMBER))
    
}

cat('Creating plots  ... \n')
# Create plots
fetal_plot <- reshape2::melt(fetal_matrix) %>%
  arrange(Var1) %>%
  group_by(Var1) %>%
  filter(row_number() >= which(Var1 == Var2)) %>%
  ggplot(aes(x = Var1, y = Var2, fill = 'white')) + 
  geom_tile(color = "black", size = 0.5, fill = '#DBF3FA') +
  geom_text(aes(label = value, size = 12)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, colour = "#000000", size = 12, hjust = TRUE),
        axis.text.y  = element_text(colour = "#000000", size = 12),
        legend.position = "none",
        panel.grid = element_blank()) +
  xlab("") + 
  ylab("") +
  coord_fixed() 

saveRDS(fetal_matrix, '~/Dropbox/BRAY_sc_analysis/files_for_paper/figures/data_for_figures/fetal_overlap_matrix.rds')

skene_plot <- skene_df %>%
  add_column(percent = prcnt) %>%
  mutate(percent = sprintf("%0.3f", percent)) %>%
  ggplot(aes(x=Var1, y=Var2, fill = 'white')) + 
  geom_tile(color = "black") +
  geom_text(aes(label = paste(value, '\n', percent))) +
  theme(axis.text.x = element_text(angle = -90, hjust = TRUE)) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(panel.grid = element_blank()) +
  xlab("") + 
  ylab("") +
  coord_flip() 

cat('Writing plots  ... \n')
# Write plots
tiff(paste0(DATA_DIR, 'fetal_intersect_plot.tiff'), height = 30, width = 30, units='cm', 
     compression = 'lzw', res = 300)
print(fetal_plot)
dev.off()

tiff(paste0(DATA_DIR, 'skene_intersect_plot.tiff'), height = 30, width = 30, units='cm', 
     compression = 'lzw', res = 300)
print(skene_plot)
dev.off()
  
cat('Done.\n')

# -------------------------------------------------------------------------------------
#                 
#              SECTION 2 - Create conditional gene lists for Magma
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  1. Intra-dataset - Remove overlapping FC-ExN-2 genes from other fetal cell
#                     Extract 10 gene lists                      
#  2. Intra-dataset - Individually remove overlapping Skene cell type genes from all 
#                   - fetal cell types: Extract 44 (4 x 11) gene lists to test on Magma


## Intra-dataset gene lists  ----------------------------------------------------------
cat('\n\n\n|-------------------------------------------|\n')
cat('\n               SECTION 2 \n')
cat('\n|-------------------------------------------|\n\n')

cat('\nRemove overlapping FC-ExN-2 genes from other fetal cell type gene lists ... \n')
cat('\nDefining variables for Intra-dataset gene lists  ... \n')
CELL_TYPES_INTRA <- c('FC_ExN_3', 'FC_ExN_4', 'FC_ExN_5', 'FC_InN_1', 'GE_InN_1', 
                      'GE_InN_2', 'Hipp_ExN_4', 'Hipp_ExN_6', 'Thal_ExN_5', 'Thal_ExN_9')

for (CELL_TYPE in CELL_TYPES_INTRA) {
  
  CELL_TYPE_DF <- get(CELL_TYPE)
    
  cat(paste0('\nProducing gene list for ', CELL_TYPE, ' ... \n'))
  
  gene_list <- list(
    FC_ExN_2 = FC_ExN_2_overlaps, 
    OTHER = get(paste0(CELL_TYPE, '_overlaps')))
  
  GENE_LIST_PLOT <- ggVennDiagram(gene_list, label_alpha = 0, edge_size = 0.5) +
    ggplot2::scale_fill_gradient(low="#FEE2E8", high = "#FF708E")
                
  # Intersections
  intersection_genes <- unlist(process_region_data(Venn(gene_list))[[3,3]])
  cat(paste0('Genes intersecting ', CELL_TYPE, ' and FC-ExN-2: ', length(intersection_genes), '\n'))
  cat(paste0('The genes are: \n\n'))
  print(intersection_genes)
  
  # Get cell type genes and remove overlapping FC-ExN-2 genes from gene list
  cat(paste0('\nGenes in ', CELL_TYPE, ' before removing intersection: ', length(CELL_TYPE_DF$X1), '\n'))
  
  cell_type_genes <- CELL_TYPE_DF %>% # https://stackoverflow.com/questions/26813667
    filter(!str_detect(X1, paste0(paste0('^', paste0(intersection_genes, collapse="$|^")), '$')))
  
  cat(paste0('Genes in ', CELL_TYPE, ' after removing intersection: ', length(cell_type_genes$X1), '\n'))
  cat(paste0('Have ', length(intersection_genes), ' genes been removed from ', CELL_TYPE, ' gene list: ', 
             length(CELL_TYPE_DF$X1) - length(cell_type_genes$X1) == length(intersection_genes), '\n'))
  
  # Convert to Entrez IDs for MAGMA
  cat(paste0('Converting HGNC IDs to Entrez for MAGMA ... \n'))
  my.genes <- as.vector(cell_type_genes$X1)
  gene_IDs <- select(hs,
                     keys = my.genes,
                     columns = c("ENTREZID", "SYMBOL"),
                     keytype = "SYMBOL")
  genes_entrez <- gene_IDs$ENTREZID
  cat(paste0('Gene count after conversion: ', length(genes_entrez)  , ' ... \n'))
  cat(paste0('NA count after conversion: ', sum(is.na(genes_entrez))  , ' ... \n'))
  
  
  #cell_type_genes[is.na(genes_entrez),]
  
  # Remove NAs
  genes_entrez_noNA <- genes_entrez[!is.na(genes_entrez)]
  cat(paste0('Duplicate count after removing NAs: ', sum(duplicated(genes_entrez_noNA))  , ' ... \n'))
  
  # Final gene count
  cat(paste0('Final gene count: ', length(genes_entrez_noNA)  , ' ... \n'))
  
  # Append name to vector for final gene list for Magma
  genes_entrez_noNA <- as.data.frame(genes_entrez_noNA)
  genes_entrez_noNA$gene_set <- rep(paste0(CELL_TYPE, '_no_FC_ExN_2'), length(genes_entrez_noNA))
  
  
  # Assign gene list
  assign(paste0(CELL_TYPE, '_no_FC_ExN_2'), genes_entrez_noNA)
  
  # Write gene list
  cat('Writing gene list ...\n')
  write_tsv(genes_entrez_noNA, paste0(DATA_DIR, CELL_TYPE, '_no_FC_ExN_2.tsv'), col_names = FALSE)
  
}

cat('Done.\n')
  
## Inter-dataset gene lists  ----------------------------------------------------------
cat('\nIndividually remove overlapping Skene cell type genes from all fetal cell type gene lists ... \n')

cat('\n---- Defining variables for Intra-dataset gene lists ----\n')
CELL_TYPES_INTER <- c('FC_ExN_2', 'FC_ExN_3', 'FC_ExN_4', 'FC_ExN_5', 'FC_InN_1', 'GE_InN_1', 
                      'GE_InN_2', 'Hipp_ExN_4', 'Hipp_ExN_6', 'Thal_ExN_5', 'Thal_ExN_9')
CELL_TYPES_SKENE <- c('skene_SS', 'skene_InN', 'skene_MSN', 'skene_CA1')
  
for (SKENE_CELL_TYPE in CELL_TYPES_SKENE) {  
  
  cat(paste0('\n--------  ', SKENE_CELL_TYPE, '  --------\n'))
  intersections_vector <- vector()
  
  for (CELL_TYPE in CELL_TYPES_INTER) {
    
    CELL_TYPE_DF <- get(CELL_TYPE)
    
    cat(paste0('\nProducing gene list for ', CELL_TYPE, ' ... \n'))
    
    gene_list <- list(
      FC_ExN_2 = get(paste0(CELL_TYPE, '_overlaps')),
      SKENE = get(paste0(SKENE_CELL_TYPE, '_overlaps')))
    
    GENE_LIST_PLOT <- ggVennDiagram(gene_list, label_alpha = 0, edge_size = 0.5) +
      ggplot2::scale_fill_gradient(low="#FEE2E8", high = "#FF708E")
    
    # Intersections
    intersection_genes <- unlist(process_region_data(Venn(gene_list))[[3,3]])
    cat(paste0('Genes intersecting ', SKENE_CELL_TYPE, ' and ', CELL_TYPE,': ', length(intersection_genes), '\n'))
    cat(paste0('The genes are: \n\n'))
    print(intersection_genes)
    intersections_vector <- c(intersections_vector, length(intersection_genes))
    
    # Get cell type genes and remove overlapping Skene genes from gene list
    cat(paste0('\nGenes in ', CELL_TYPE, ' before removing intersection: ', length(CELL_TYPE_DF$X1), '\n'))
    
    cell_type_genes <- CELL_TYPE_DF %>%
      filter(!str_detect(X1, paste0(paste0('^', paste0(intersection_genes, collapse="$|^")), '$')))
    
    cat(paste0('Genes in ', CELL_TYPE, ' after removing intersection: ', length(cell_type_genes$X1), '\n'))
    cat(paste0('Have ', length(intersection_genes), ' genes been removed from ', CELL_TYPE, ' gene list: ', 
               length(CELL_TYPE_DF$X1) - length(cell_type_genes$X1) == length(intersection_genes), '\n'))
    
    # Convert to Entrez IDs for MAGMA
    cat(paste0('Converting HGNC IDs to Entrez for MAGMA ... \n'))
    my.genes <- as.vector(cell_type_genes$X1)
    gene_IDs <- select(hs,
                       keys = my.genes,
                       columns = c("ENTREZID", "SYMBOL"),
                       keytype = "SYMBOL")
    genes_entrez <- gene_IDs$ENTREZID
    cat(paste0('Gene count after conversion: ', length(genes_entrez)  , ' ... \n'))
    cat(paste0('NA count after conversion: ', sum(is.na(genes_entrez))  , ' ... \n'))
    
    
    #cell_type_genes[is.na(genes_entrez),]
    
    # Remove NAs
    genes_entrez_noNA <- genes_entrez[!is.na(genes_entrez)]
    cat(paste0('Duplicate count after removing NAs: ', sum(duplicated(genes_entrez_noNA))  , ' ... \n'))
    
    # Final gene count
    cat(paste0('Final gene count: ', length(genes_entrez_noNA)  , ' ... \n'))
    
    # Compare to using Skene's magma celltypng hgnc to entrez conversion list
    # Taken from:
    #   https://github.com/neurogenomics/MAGMA_Celltyping/blob/f931877dedcde103723f14242d5242fdca7b3af6/R/map_specificity_to_entrez.r#L19
    #   all_hgnc_wtEntrez<- MAGMA.Celltyping::all_hgnc_wtEntrez
    #all_hgnc_wt_entrez <- read_delim("Desktop/single_nuclei/snRNAseq/final_analyses/all_hgnc_wt_entrez.tsv", 
     #                                delim = "\t", escape_double = FALSE, 
     #                                trim_ws = TRUE)
    
    #colnames(all_hgnc_wt_entrez)[1] = "human.symbol"
    humanSymsPresent = as.character(all_hgnc_wt_entrez$human.symbol[all_hgnc_wt_entrez$human.symbol %in% cell_type_genes$X1])
    cat(paste0('Gene conversion count after using Skene conversion list: ', length(humanSymsPresent)  , ' ... \n'))
    
    conversion_check <- cbind(paste0(CELL_TYPE, '_no_', SKENE_CELL_TYPE), length(cell_type_genes$X1), sum(is.na(genes_entrez)), 
          length(genes_entrez_noNA), length(humanSymsPresent))
    
    if (exists("conversion_check_df")) {
      
      conversion_check_df <- rbind(conversion_check_df, conversion_check)
      
    } else {
      
      conversion_check_df <- conversion_check 
      colnames(conversion_check_df) <- c('cell_type', 'initial_gene_cnt', 'NAs_BiomaRt', 'final_cnt', 'Skene_cnt')
      
    }
    
    # Assign gene list
    assign(paste0(CELL_TYPE, '_no_', SKENE_CELL_TYPE, '_hgnc'), cell_type_genes)
    assign(paste0(CELL_TYPE, '_no_', SKENE_CELL_TYPE, '_entrez'), genes_entrez_noNA)
    
    # Append name to vector for final gene list for Magma
    genes_entrez_noNA <- as.data.frame(genes_entrez_noNA)
    genes_entrez_noNA$gene_set <- rep(paste0(CELL_TYPE, '_no_', SKENE_CELL_TYPE), length(genes_entrez_noNA))
    
    
    # Write gene list
    cat('Writing gene list ...\n')
    write_tsv(genes_entrez_noNA, paste0(DATA_DIR, CELL_TYPE, '_no_', SKENE_CELL_TYPE, '.tsv'), col_names = FALSE)
    
  }
  
  # Assign the intersections vector for the skene cell types
  assign(paste0(SKENE_CELL_TYPE, '_intersections'), intersections_vector)
  
  cat(paste0('List of ', SKENE_CELL_TYPE, '  intersection counts: \n\n'))
  print(intersections_vector)
  
}

# Need to convert gene lists back to Entrez for MAGMA
cat('Convert gene IDs for each cell group from HGNC to Entrez ... \n')




#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


