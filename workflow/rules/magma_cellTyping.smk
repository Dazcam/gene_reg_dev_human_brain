# -------------------------------------------------------------------------------------
#
#    Script for MAGMA cell typing analysis
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARMAS  ----------

configfile: "config/config.yaml"
localrules: get_genes_in_MHC, create_ctd, map_genes_to_snps


# ----------  SET VARIABLES  ----------

SCRATCH =   "/scratch/c.c1477909/"
ROBJ_DIR = SCRATCH + "R_objects/"
MAGMA_DIR = SCRATCH + "magma_celltyping/"
LOG_DIR = SCRATCH + "logs/"
GWAS_DIR =  MAGMA_DIR + "GWAS_for_magma/"
MAGMA_OUTDIR = GWAS_DIR + "MAGMA_Files/"
GENE_LIST_OUTDIR = MAGMA_DIR + "ldsc_gene_lists/"
MARKDOWN_DIR = SCRATCH + "markdown/"
RESULTS_DIR = SCRATCH + "results/"

# -------------  RULES  ---------------

rule all:
    input:
        RESULTS_DIR + "snRNAseq_magma_generate_plots.html", 
        expand(GENE_LIST_OUTDIR + "{REGION}_complete.file", REGION = config["RNA_REGIONS"])

rule get_genes_in_MHC:
    input:  ROBJ_DIR + "seurat.{REGION}.final.rds"
    output: ROBJ_DIR + "{REGION}_MHC_overlapping_genes.rds"
    params: outdir = ROBJ_DIR
    log:    LOG_DIR + "magma/get_genes_in_MHC_{REGION}.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/snRNAseq_magma_get_MHC_overlap_genes.R {wildcards.REGION} {input} {params.outdir} 2> {log}          

            """

rule create_ctd:
    input:  seurat_obj = ROBJ_DIR + "seurat.{REGION}.final.rds",
            MHC_gene_list = ROBJ_DIR + "{REGION}_MHC_overlapping_genes.rds"
    output: MAGMA_DIR + "ctd_objects/CellTypeData_{REGION}.rda"
    log:    LOG_DIR + "magma/create_ctd_{REGION}.log" 
    shell: 
            """
            
            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/snRNAseq_magma_create_ctd.R {wildcards.REGION} {input.seurat_obj} {input.MHC_gene_list} 2> {log}

            """

rule map_genes_to_snps:
    # Requires net access to run
    input:  GWAS_DIR + "{GWAS}_hg19_magma_ready.sumstats.tsv"
    output: MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN/{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.raw"
    log:   LOG_DIR + "magma/map_snps2genes_hg19_{GWAS}.log"
    shell:
            """
       
            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/snRNAseq_magma_map_snps2genes.R {input} 2> {log}
              
            """

rule magma_analysis:
    input:  ctd_obj = MAGMA_DIR + "ctd_objects/CellTypeData_{REGION}.rda",
            gene_file = MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN/{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.raw", 
            gwas_file = GWAS_DIR + "{GWAS}_hg19_magma_ready.sumstats.tsv" 
    output: MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN/{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.level2.{REGION}_top10.gsa.out"
    log:    LOG_DIR + "magma/magma_analysis_hg19_{REGION}_{GWAS}.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/snRNAseq_magma_analysis.R {wildcards.REGION} {input.ctd_obj} {input.gwas_file} 2> {log}

            """

rule magma_generate_plots:
    input:  magma = expand(MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN/{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.level2.{REGION}_top10.gsa.out", REGION = config["RNA_REGIONS"], GWAS = config["SUMSTATS"]), 
            markdown = MARKDOWN_DIR + "snRNAseq_magma_generate_plots.Rmd"     
    output: RESULTS_DIR + "snRNAseq_magma_generate_plots.html"
    params: magma_dir = MAGMA_OUTDIR,
            report_dir = RESULTS_DIR,
            report_file = "snRNAseq_magma_generate_plots.html",
    log:    LOG_DIR + "magma/magma_generate_plots.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/snRNAseq_magma_generate_plots.R {params.magma_dir} {input.markdown} {params.report_dir} {params.report_file}  2> {log}

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
            scripts/R/snRNAseq_magma_create_LDSC_geneLists.R {wildcards.REGION} {input.ctd_obj} 2> {log}

            """


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
