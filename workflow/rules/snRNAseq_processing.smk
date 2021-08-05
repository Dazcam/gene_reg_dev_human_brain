# -------------------------------------------------------------------------------------
#
#    QC and processing snRNA-seq data 
#
# -------------------------------------------------------------------------------------

# Not tested

# ---------  SET SMK PARMAS  ----------

configfile: "config/config.yaml"

# ----------  SET VARIABLES  ----------

SCRATCH =   "/scratch/c.c1477909/"
RNA_DATA_DIR = SCRATCH + "results/cell_ranger_RNA"
ROBJ_DIR = SCRATCH + "results/R_objects/"
LOG_DIR = SCRATCH + "logs/"


# -------------  RULES  ---------------

rule all:
    input: 
        expand(ROBJ_DIR + "rna.{REGION}.seurat.rds", REGION = config['REGIONS']),
        expand(ROBJ_DIR + "rna.{REGION}.DEgenes.rds", REGION = config['REGIONS'])
          
          
rule snRNA_seq_QC:
    # region in inout may not be needed here
    input:  data_dir = dir(RNA_DATA_DIR),
            region = config['REGIONS']
    output: ROBJ_DIR + "rna.{REGION}.sce.QC.rds" 
    log:    LOG_DIR + "data_processing/snRNA_{REGION}_QC.log" 
    shell: 
            """
            
            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/snRNAseq_seurat.R {wildcards.REGION} {input.data_dir} {output} 2> {log}
            """
        
rule snRNA_seq_Seurat:
    input:  ROBJ_DIR + "sce.rna.{REGION}.QC.rds"
    output: seurat_obj = ROBJ_DIR + "sce.rna.{REGION}.seurat.rds" 
            topDEgenes = ROBJ_DIR + "~/sce.rna.{REGION}.DEgenes.rds"
    log:    LOG_DIR + "data_processing/snRNA_{REGION}_seurat.log" 
    shell: 
            """
            
            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/snRNAseq_QC.R {wildcards.REGION} {input} {output.seurat_obj} {output.topDEgenes} 2> {log}
            """
        
        
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
