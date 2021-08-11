# -------------------------------------------------------------------------------------
#
#    QC and processing snATAC-seq data 
#
# -------------------------------------------------------------------------------------

# Not tested

# ---------  SET SMK PARMAS  ----------

configfile: "config/config.yaml"

# ----------  SET VARIABLES  ----------

SCRATCH =   "/scratch/c.c1477909/"
ATAC_DATA_DIR = SCRATCH + "results/snATACseq_CR-atac_1.2.0"
MARKDOWN_DIR = SCRATCH + "markdown/"
ROBJ_DIR = SCRATCH + "results/R_objects/"
LOG_DIR = SCRATCH + "logs/"


# -------------  RULES  ---------------

rule all:
    input: 
        expand(MARKDOWN_DIR + "snATACseq_QC.html"),
#        expand(ROBJ_DIR + "rna.{REGION}.DEgenes.rds", REGION = config['REGIONS'])
          
          
rule snATAC_seq_QC:
    # region in inout may not be needed here
    input:  data_dir = dir(ATAC_DATA_DIR),
            region = config['ATAC_REGIONS']
    output: MARKDOWN_DIR + "snATACseq_QC.html" 
    log:    LOG_DIR + "data_processing/snATAC_{REGION}_QC.log" 
    shell: 
            """
            
            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/snATACseq_QC.R {wildcards.REGION} {input.data_dir} {output} 2> {log}
            """
        
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
