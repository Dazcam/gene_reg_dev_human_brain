#--------------------------------------------------------------------------------------
#
#    ArchR - Motif and peak to gene linkage analysis analysis
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# ArchR manual - https://www.archrproject.com/index.html
# ArchR GitHiub - https://github.com/GreenleafLab/ArchR
# Summarized Expriment - https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
# Harmony - Github -   https://github.com/immunogenomics/harmony
# Granges - https://kasperdanielhansen.github.io/genbioconductor/html/GenomicRanges_GRanges_Usage.html

## Info  ------------------------------------------------------------------------------

#  snATAC-seq - Running additonal ArchR analyses on peaks

##  Load Packages  --------------------------------------------------------------------
library(ArchR)
library(tidyverse)
library(rmarkdown)
library(chromVARmotifs)
library(ComplexHeatmap)
library(argparser)
library(pheatmap)
library(cowplot)
library(ggplotify) # required to convert grob to ggplot object - footprinting


## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("\nRead markdown and report file and directory info ... \n")
p <- add_argument(p, "markdown_file", help = "No markdown file path specified")
p <- add_argument(p, "report_dir", help = "No report output directory specified")
p <- add_argument(p, "report_file", help = "No report filename specified")
args <- parse_args(p)
print(args)

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
OUT_DIR <- "results/reports/ATAC/ARCHR/"
MARKDOWN_FILE <- args$markdown_file
REPORT_DIR <- args$report_dir
REPORT_FILE <- args$report_file
addArchRThreads(threads = 12) # Set Hawk to 32 cores so 0.75 of total
addArchRGenome("hg38")

for (REGION in c("FC", "GE")) {
  
  ##  Load ArchR project  -------------------------------------------------------------------
  cat(paste0('\nLoading ArchR project for ', REGION, ' ... \n'))
  archR <- loadArchRProject(paste0(OUT_DIR, REGION))
  
  cat(paste0('\nAdding peak matrix ... \n'))
  archR <- addPeakMatrix(archR) 
  
  ##  Create data for figs 5A-D  --------------------------------------------------------
  cat(paste0('\nCreating rds files for figs 5A-D ... \n'))
  if (REGION == 'FC') {
  
    fig_5A <- plotEmbedding(archR, colorBy = "cellColData", name = "Clusters_RNAmapped") +
      ggtitle(NULL) +
      theme_bw() +
      theme(legend.text=element_text(size = 12),
            legend.title = element_blank(),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(colour = "#000000", size = 14),
            axis.title.y = element_text(colour = "#000000", size = 14),
            axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 12)) 
    
    fig_5C <- plotEmbedding(archR, colorBy = "cellColData", name = "Clusters_RNAmapped_broad") +
      ggtitle(NULL) +
      theme_bw() +
      theme(legend.text=element_text(size = 12),
            legend.title = element_blank(),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(colour = "#000000", size = 14),
            axis.title.y = element_text(colour = "#000000", size = 14),
            axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 12)) 
    
    saveRDS(fig_5A, 'results/figures/Fig_5A.rds') 
    saveRDS(fig_5C, 'results/figures/Fig_5C.rds') 
  
  } else {
    
    fig_5B <- plotEmbedding(archR, colorBy = "cellColData", name = "Clusters_RNAmapped") +
      ggtitle(NULL) +
      theme_bw() +
      theme(legend.text=element_text(size = 12),
            legend.title = element_blank(),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(colour = "#000000", size = 14),
            axis.title.y = element_text(colour = "#000000", size = 14),
            axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 12)) 
    
    fig_5D <- plotEmbedding(archR, colorBy = "cellColData", name = "Clusters_RNAmapped_broad") +
      ggtitle(NULL) +
      theme_bw() +
      theme(legend.text=element_text(size = 12),
            legend.title = element_blank(),
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(colour = "#000000", size = 14),
            axis.title.y = element_text(colour = "#000000", size = 14),
            axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 12)) 
    
    saveRDS(fig_5B, 'results/figures/Fig_5B.rds') 
    saveRDS(fig_5D, 'results/figures/Fig_5D.rds') 
    
  }

  ##  IDing Marker peaks that are unique to individual cluster - chptr 11  --------------
  # scRNA labels
  table(archR$Clusters_RNAmapped_broad)
  
  cat(paste0('\nAdding marker peaks for each cell type ... \n'))
  # Add marker peaks - returns summarisedExperiment with 6 assays
  markersPeaks <- getMarkerFeatures(
    ArchRProj = archR, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters_RNAmapped_broad",
    bias = c("TSSEnrichment", "log10(nFrags)"), # Correction fof diffs in data quality
    testMethod = "wilcoxon"
  )
  markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
  markerList
  
  ## Motif enrichment  -  Chptr 12  -----------------------------------------------------
  # Add motif anns to archR project - creates binary matrix for presence/absence of motif in peak
  cat(paste0('\nAdding motif annotations  ... \n'))
  archR <- addMotifAnnotations(ArchRProj = archR, motifSet = "cisbp", name = "Motif", force = TRUE)
  
  ##   --------  Motif enrichment in marker peaks  ----------------- 
  cat(paste0('\nAssessing motif enrichment ... \n'))
  enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = archR,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  enrichMotifs
  
  # Export enrichMotifs and motif p-val table used for plotting
  saveRDS(enrichMotifs, paste0('results/figures/', REGION, '_motifs.rds'))
  saveRDS(assays(enrichMotifs)[["mlog10Padj"]], paste0('results/figures/', REGION, '_motifs_mlogPs.rds'))

  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 15, transpose = TRUE)
  assign(paste0("motif_heatmap_", REGION), heatmapEM)
  
  ## Create motif heatmap RDS files for figures 5E-F  -----------------------------------
  cat(paste0('\nCreating rds files for figs 5E-F ... \n'))
  # Return heatmap matrix
  heatmap_matrix <- plotEnrichHeatmap(enrichMotifs, n = 15, transpose = TRUE,
                                      returnMatrix = TRUE)
  
  # Remove all extra info from motif names
  new_cols <- as.data.frame(str_split(colnames(heatmap_matrix), 
                                      fixed("_"), simplify = TRUE)) %>%
    pull(V1)
  colnames(heatmap_matrix) <- new_cols
    
  # Note that Complex heatmap also has a function called heatmap for object conversion!!
  motif_heatmap_plot <- as.ggplot(pheatmap::pheatmap(heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                           cellwidth = 12, cellheight = 15, fontsize_row = 12,
                           fontsize_col = 12)) 
  
  # Save motif heatmap RDS files
  if (REGION == 'FC') {
  
    saveRDS(motif_heatmap_plot, 'results/figures/Fig_5E.rds') 
    saveRDS(heatmap_matrix, 'results/figures/Fig_5E_matrix.rds')    

  } else {
    
    saveRDS(motif_heatmap_plot, 'results/figures/Fig_5F.rds') 
    saveRDS(heatmap_matrix, 'results/figures/Fig_5F_matrix.rds')    

  }


  ## Motif Footprinting  -  Chptr 14  ---------------------------------------------------
  cat(paste0('\nRunning footprinting ... \n'))
  if (REGION == 'FC') {

  
    cat(paste0('\nGenerating motif footprints for ', REGION, ' ... \n'))
    motifPositions <- getPositions(archR)
    motifPositions
  
  cat(paste0('\nCreating rds files for figs 5G-H ... \n'))
  for (MOTIF in c('NEUROD2_73', 'DLX5_412')) {
      
      motifs <- MOTIF
      markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
      markerMotifs
      
      # Compute footprints
      cat(paste0(MOTIF, ' ... \n'))
      seFoot <- getFootprints(
        ArchRProj = archR, 
        positions = motifPositions[markerMotifs], 
        groupBy = "Clusters_RNAmapped_broad"
      )
      
      # Plot footprint - plot = FALSE required to get grob object
      cat(paste0('Plotting ... \n'))
      footprint_grob <-  plotFootprints(
          seFoot = seFoot,
          ArchRProj = archR, 
          normMethod = "Subtract",
          plotName = "Footprints_subtract_bias",
          addDOC = FALSE,
          smoothWindow = 5,
          plot = FALSE
          
        )
    
      # Convert to ggplot object 
      cat(paste0('Converting grob to ggplot object ... \n'))
      footprint_plot <- ggplotify::as.ggplot(footprint_grob[[MOTIF]])
      
      if (MOTIF == 'NEUROD2_73') {
        
        saveRDS(footprint_plot, 'results/figures/Fig_5G.rds') 
        saveRDS(footprint_grob, 'results/figures/Fig_5G_grob.rds')        

      } else {
        
        saveRDS(footprint_plot, 'results/figures/Fig_5H.rds') 
        saveRDS(footprint_grob, 'results/figures/Fig_5H_grob.rds')
      }
    
    }

  }

  ## Peak to gene links  -  Chptr 15.3  -------------------------------------------------
  cat(paste0("Starting peak to gene linkage analysis for ", REGION, " ... \n"))
  cat("Adding Peak to gene links ... \n")
  archR <- addPeak2GeneLinks(
  ArchRProj = archR,
  reducedDims = "Harmony"
  )

  cat("\nRetrieving Peak to gene links ... \n")
  p2g <- getPeak2GeneLinks(
  ArchRProj = archR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
  )

  cat("\nCreating Peak to gene links dataframes ... \n")
  p2g_df <- base::as.data.frame(p2g)
  geneIDs_df <- Repitools::annoGR2DF(metadata(p2g)$geneSet)
  peakIDs_df <- Repitools::annoGR2DF(metadata(p2g)$peakSet)
  
  markerGenes  <- c("CALN1", "FUT9", "DLX1", "DLX2")

  p <- plotBrowserTrack(
    ArchRProj = archR, 
    groupBy = "Clusters_RNAmapped_broad", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(archR)
)

  # Plot p2g links
  cat(paste0("Plotting p2g links for ", REGION, " ... \n"))
  tiff(paste0("results/figures/CALN1_peak2gene", REGION, ".tiff"), height = 10, width = 10, units='cm', 
       compression = "lzw", res = 300)
  grid::grid.newpage()
  grid::grid.draw(p$CALN1)
  dev.off()

  tiff(paste0("results/figures/FUT9_peak2gene", REGION, ".tiff"), height = 10, width = 10, units='cm', 
       compression = "lzw", res = 300)
  grid::grid.newpage()    
  grid::grid.draw(p$FUT9)
  dev.off()

  tiff(paste0("results/figures/DLX1_peak2gene", REGION, ".tiff"), height = 10, width = 10, units='cm', 
       compression = "lzw", res = 300)
  grid::grid.newpage()
  grid::grid.draw(p$DLX1)
  dev.off()

  tiff(paste0("results/figures/DLX2_peak2gene", REGION, ".tiff"), height = 10, width = 10, units='cm', 
       compression = "lzw", res = 300)
  grid::grid.newpage()
  grid::grid.draw(p$DLX2)
  dev.off()
  

  cat("\nObtaining peak start/stop coordinates and gene IDs for links ... \n")

  # Need to remove dfs between regions or region specific entries will be appended
  if (exists("peak2gene_df")) { rm(peak2gene_df) }
  if (exists("peak2gene_final_df")) { rm(peak2gene_final_df) }
  
  for (LOOP in 1:nrow(p2g_df)) {
    
    gene_index <- p2g_df[LOOP, 2]
    peak_index <- p2g_df[LOOP, 1]
    
    gene_TSS <- geneIDs_df[gene_index, ]$start
    gene_ID <- geneIDs_df[gene_index, ]$name
    peak_start <- peakIDs_df[peak_index, ]$start
    peak_end <- peakIDs_df[peak_index, ]$end
    chr <- as.vector(peakIDs_df[peak_index, ]$chr)
    
    if (exists("peak2gene_df")) {
      
      peak2gene_row <- cbind(chr, peak_start, peak_end, gene_TSS, gene_ID)
      peak2gene_df <- rbind(peak2gene_df, peak2gene_row)
      
    } else {
      
      peak2gene_df <- cbind(chr, peak_start, peak_end, gene_TSS, gene_ID)
      
    }
    
    
  }
  
  # Join peak start/stop coordinates and gene IDs to original table
  peak2gene_final_df <- cbind(peak2gene_df, p2g_df)
  
  
  cat("\nWriting table to file ... \n")
  write_tsv(peak2gene_final_df, paste0('results/peak2gene_table_', REGION, '.tsv'))
  
}  




## Create markdown doc  ---------------------------------------------------------------
cat('\nCreating markdown report ... \n')
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)

cat('\nDONE.\n')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
