#--------------------------------------------------------------------------------------
#
#    ArchR - Motif and Co-accessibility analysis
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
  
  heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
  assign(paste0("motif_heatmap_", REGION), heatmapEM)

  ## Motif Footprinting  -  Chptr 14  ---------------------------------------------------
  if (REGION == 'FC') {

  
    cat(paste0('\nGenerating motif footprints for ', REGION, ' ... \n'))
    motifPositions <- getPositions(archR)
    motifPositions

    # Subset the GRangesList to a few TFs
    motifs <- c("NEUROD2", "NEUROD6", "ASCL1", "TCF12")
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

}

## Create markdown doc  ---------------------------------------------------------------
cat('\nCreating markdown report ... \n')
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)

cat('\nDONE.\n')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
