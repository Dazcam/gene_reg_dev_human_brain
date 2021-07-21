# -------------------------------------------------------------------------------------
#
#
#    Script for MAGMA cell typing analysis
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARMAS  ----------

configfile: "config.yaml"
localrules: map_genes_to_snps


# ----------  SET VARIABLES  ----------

SCRATCH =   "/scratch/c.c1477909/"
ROBJ_DIR = SCRATCH + "R_objects/"
MAGMA_DIR = SCRATCH + "magma_celltyping/"
LOG_DIR = SCRATCH + "logs/"
GWAS_DIR =  MAGMA_DIR + "GWAS_for_magma/"
MAGMA_OUTDIR = GWAS_DIR + "MAGMA_Files/"
GENE_LIST_OUTDIR = MAGMA_DIR + "ldsc_gene_lists/"

# -------------  RULES  ---------------

rule all:
  #  input: expand(MAGMA_DIR + "ctd_objects/CellTypeData_{REGION}.rda", REGION = config["REGIONS"])
  #  input:  expand(MAGMA_OUTDIR + "{GWAS}_hg38_magma_ready.sumstats.tsv.10UP.1.5DOWN/{GWAS}_hg38_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.raw", GWAS = config["SUMSTATS_MAGMA"]) 
    input: 
        expand(MAGMA_OUTDIR + "{GWAS}_hg38_magma_ready.sumstats.tsv.10UP.1.5DOWN/{GWAS}_hg38_magma_ready.sumstats.tsv.10UP.1.5DOWN.level2.{REGION}_top10.gsa.out", REGION = config["REGIONS"], GWAS = config["SUMSTATS"]),
        expand(GENE_LIST_OUTDIR + "{REGION}_complete.file", REGION = config["REGIONS"])
          
rule create_ctd:
    input:  ROBJ_DIR + "seurat.{REGION}.final.rds"
    output: MAGMA_DIR + "ctd_objects/CellTypeData_{REGION}.rda"
    log:   LOG_DIR + "magma/create_ctd_{REGION}.log" 
    shell: 
            """
            
            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/magma_create_ctd.R {wildcards.REGION} {input} 2> {log}

            """

rule map_genes_to_snps:
    # Requires net access to run
    input:  GWAS_DIR + "{GWAS}_hg38_magma_ready.sumstats.tsv"
    output: MAGMA_OUTDIR + "{GWAS}_hg38_magma_ready.sumstats.tsv.10UP.1.5DOWN/{GWAS}_hg38_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.raw"
    log:   LOG_DIR + "magma/map_genes_to_snps_{GWAS}.log"
    shell:
            """
       
            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/magma_map_genes2snps.R {input} 2> {log}
              
            """

rule magma_analysis:
    input:  ctd_obj = MAGMA_DIR + "ctd_objects/CellTypeData_{REGION}.rda",
            gene_file = MAGMA_OUTDIR + "{GWAS}_hg38_magma_ready.sumstats.tsv.10UP.1.5DOWN/{GWAS}_hg38_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.raw", 
            gwas_file = GWAS_DIR + "{GWAS}_hg38_magma_ready.sumstats.tsv" 
    output: MAGMA_OUTDIR + "{GWAS}_hg38_magma_ready.sumstats.tsv.10UP.1.5DOWN/{GWAS}_hg38_magma_ready.sumstats.tsv.10UP.1.5DOWN.level2.{REGION}_top10.gsa.out"
    log:    LOG_DIR + "magma/magma_analysis_{REGION}_{GWAS}.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/magma_analysis.R {wildcards.REGION} {input.ctd_obj} {input.gwas_file} 2> {log}

            """

rule magma_create_LDSC_geneLists:
    input:  ctd_obj = MAGMA_DIR + "ctd_objects/CellTypeData_{REGION}.rda",
    output: GENE_LIST_OUTDIR + "{REGION}_complete.file"
    log:    LOG_DIR + "magma/magma_ldsc_gene_lists_{REGION}.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/magma_create_LDSC_geneLists.R {wildcards.REGION} {input.ctd_obj} 2> {log}

            """


