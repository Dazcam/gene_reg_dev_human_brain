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
  
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 15, transpose = TRUE)
  assign(paste0("motif_heatmap_", REGION), heatmapEM)

  ## Motif Footprinting  -  Chptr 14  ---------------------------------------------------
  if (REGION == 'FC') {

  
    cat(paste0('\nGenerating motif footprints for ', REGION, ' ... \n'))
    motifPositions <- getPositions(archR)
    motifPositions

    # Subset the GRangesList to a few TFs
    motifs <- c('MEF2C', 'NEUROD4', 'NEUROD2', 'ASCL1', 'ASCL2', 'DLX5', 'PRRX1')
    markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
    #markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
    markerMotifs

    # Compute footprints - Note: Any user-defined feature set can be footprinted
    seFoot <- getFootprints(
    ArchRProj = archR, 
    positions = motifPositions[markerMotifs], 
    groupBy = "Clusters_RNAmapped_broad"
    )

    # Plot footprints - Note: multiple normalisation methods that can be used
    plotFootprints(
      seFoot = seFoot,
      ArchRProj = archR, 
      normMethod = "Subtract",
      plotName = "Footprints-Subtract-Bias",
      addDOC = FALSE,
      smoothWindow = 5
    )
  
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

