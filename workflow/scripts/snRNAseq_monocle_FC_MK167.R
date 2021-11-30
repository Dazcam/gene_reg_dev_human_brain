#--------------------------------------------------------------------------------------
#
#    snRNAseq monocle analysis MK167 - laptop
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Running monocle analysis using MK167 as the root gene 
#  Using the seurat wrapper code verbatim:
#  http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html

## Set variables  ---------------------------------------------------------------------
DATA_DIR <- '~/Desktop/single_cell/scRNAseq/batch2_CR5_200121/r_objects/final/'

## Load Data  -------------------------------------------------------------------------
seurat.pfc <- readRDS(paste0(DATA_DIR, 'seurat.pfc.final.rds'))

## Load Data  -------------------------------------------------------------------------
cds <- as.cell_data_set(seurat.pfc)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

## Run pseudotime analysis  -----------------------------------------------------------
# Choose partition
integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)

# Use MK167 as root
max.avp <- which.max(unlist(FetchData(integrated.sub, "MKI67")))
max.avp <- colnames(integrated.sub)[max.avp]
cds <- order_cells(cds, root_cells = max.avp)
plot_cells(cds, color_cells_by = "cellIDs", label_leaves = FALSE, 
           label_branch_points = FALSE, label_groups_by_cluster = TRUE, 
           show_trajectory_graph = FALSE)

# Save cds object
saveRDS(cds, paste0(DATA_DIR, "pfc_monocle_MK167.rds"))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
