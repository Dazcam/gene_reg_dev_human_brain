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
ATAC_DATA_DIR = SCRATCH + "snATACseq_CR-atac_1.2.0/"
MARKDOWN_DIR = SCRATCH + "markdown/"
ROBJ_DIR = SCRATCH + "results/R_objects/"
LOG_DIR = SCRATCH + "logs/"
RESULTS_DIR = SCRATCH + "results/"

# -------------  RULES  ---------------

rule all:
    input:
	expand(RESULTS_DIR + "snATACseq_remove_batch_effects_{REGION}.html", REGION = config['ATAC_REGIONS']),


rule snATAC_seq_QC:
    input:  markdown = MARKDOWN_DIR + "snATACseq_QC.Rmd"
    output: RESULTS_DIR + "snATACseq_QC_{REGION}.html"
    params: data_dir = ATAC_DATA_DIR,
            archR_out_dir = RESULTS_DIR + "ARCHR/{REGION}",
            report_dir = RESULTS_DIR,
            report_file = "snATACseq_QC_{REGION}.html"
    log:    LOG_DIR + "data_processing/snATAC_QC_{REGION}.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/snATACseq_QC.R {wildcards.REGION} {params.data_dir} \
            {params.archR_out_dir} {input.markdown} {params.report_dir} {params.report_file} 2> {log}
            """
            
rule snATAC_remove_batch_effects:
    input:  markdown = MARKDOWN_DIR + "snATACseq_remove_batch_effects.Rmd",
            qc_html = RESULTS_DIR + "snATACseq_QC_{REGION}.html" # Needed for rule order
    output: RESULTS_DIR + "snATACseq_remove_batch_effects_{REGION}.html"
    params: data_dir = ATAC_DATA_DIR,
            archR_out_dir = RESULTS_DIR + "ARCHR/{REGION}",
            report_dir = RESULTS_DIR,
            report_file = "snATACseq_remove_batch_effects_{REGION}.html"
    log:    LOG_DIR + "data_processing/snATAC_remove_batch_effects_{REGION}.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/snATACseq_remove_batch_effects.R {wildcards.REGION} {params.data_dir} \
            {params.archR_out_dir} {input.markdown} {params.report_dir} {params.report_file} 2> {log}
            """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
