# -------------------------------------------------------------------------------------
#
#
#    Script for creating snRNA-seq data LDSC plots
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "config/config.yaml"
report: "report/workflow.rst" # Experimented with this but not currently using

# ----------  SET VARIABLES  ----------
SCRATCH = "/scratch/c.c1477909/"
LDSC_DIR = SCRATCH + "ldsc/"
REF_DIR = LDSC_DIR + "reference_files/"
RESULTS_DIR = SCRATCH + "results/"
MARKDOWN_DIR = SCRATCH + "markdown/"
MAGMA_DIR = SCRATCH + "magma_celltyping/"
ANN_DIR = RESULTS_DIR + "LDSR_annotation_files/RNA/"
PART_HERIT_DIR = RESULTS_DIR + "LDSR_part_herit_files/RNA/"
GWAS_DIR = LDSC_DIR + "GWAS_for_ldsc/"
QUANTILES_DIR = MAGMA_DIR + "ldsc_gene_lists/"


# -------------  RULES  ---------------

rule all:
    input:
#         expand(PART_HERIT_DIR + "snRNAseq_LDSC_{CELL_TYPE}_{GWAS}_summary.tsv", CELL_TYPE = config["RNA_CELL_TYPES"], GWAS = config["LDSC_GWAS"]),
#         expand(PART_HERIT_DIR + "snRNAseq_LDSC_{GWAS}_top10pc.tsv", GWAS = config["LDSC_GWAS"])
#          expand(PART_HERIT_DIR + "snRNAseq_LDSC_{CELL_TYPE}_{GWAS}.rds", CELL_TYPE = config["RNA_CELL_TYPES"], GWAS = config["LDSC_GWAS"])          
          RESULTS_DIR + "snRNAseq_LDSC_report.html"

rule create_partHerit_summary:
    input:   expand(PART_HERIT_DIR + "snRNAseq_LDSC_{CELL_TYPE}_Q{QUANTILE}_{GWAS}.results", CELL_TYPE = config["RNA_CELL_TYPES"], QUANTILE = range(0,11), GWAS = config["LDSC_GWAS"])
    output:  PART_HERIT_DIR + "snRNAseq_LDSC_{CELL_TYPE}_{GWAS}_summary.tsv"
    message: "Creating summary file for {wildcards.CELL_TYPE} and {wildcards.GWAS} GWAS"
    params:  dir = PART_HERIT_DIR
    log:     "logs/LDSR/snRNAseq.{CELL_TYPE}.{GWAS}_partHerit.summary.log"
    shell:
        """
        
        head -1 {params.dir}snRNAseq_LDSC_Cer-RG-1_Q1_SCZ.results > {output} 
        grep L2_1 {params.dir}snRNAseq_LDSC_{wildcards.CELL_TYPE}_Q*_{wildcards.GWAS}.results >> {output}
        
        """


rule create_top_decile_tables:
    input:   expand(PART_HERIT_DIR + "snRNAseq_LDSC_{CELL_TYPE}_{GWAS}_summary.tsv", CELL_TYPE = config["RNA_CELL_TYPES"], GWAS = config["LDSC_GWAS"])
    output:  PART_HERIT_DIR + "snRNAseq_LDSC_{GWAS}_top10pc.tsv"
    message: "Creating LDSC top decile tables for {wildcards.GWAS} GWAS"
    params:  dir = PART_HERIT_DIR
    log:     "logs/LDSR/snRNAseq.{GWAS}_partHerit.top10pc_summary.log"
    shell:
             """
             
             head -1 {params.dir}snRNAseq_LDSC_Cer-OPC_Q10_BPD.results > {output}

             for file in `ls {params.dir}*Q10_{wildcards.GWAS}*`; do

             CELL_TYPE=$(echo ${{file}} | cut -d'_' -f6) 
             tail -1 ${{file}} >> {output}
             sed -i "s/L2_1/${{CELL_TYPE}}/g" {output}
             sed -i '/Total time elapsed/d' {output}

             done 

             """
        
rule create_ldsc_plots:
    input:   PART_HERIT_DIR + "snRNAseq_LDSC_{CELL_TYPE}_{GWAS}_summary.tsv"
    output:  rds = PART_HERIT_DIR + "snRNAseq_LDSC_{CELL_TYPE}_{GWAS}.rds",
             plot = PART_HERIT_DIR + "snRNAseq_LDSC_{CELL_TYPE}_{GWAS}.png"
    message: "Creating ldsc plots for {wildcards.CELL_TYPE} and {wildcards.GWAS} GWAS"
    log:     "logs/LDSR/snRNAseq.{CELL_TYPE}.{GWAS}_partHerit.plots.log"
    shell:
             """
             
             export R_LIBS_USER=/scratch/c.c1477909/R/library
             module load libgit2/1.1.0
             /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
             scripts/R/snRNAseq_LDSC_create_plots.R {wildcards.CELL_TYPE} {input} {output.rds} {output.plot} 2> {log}

             """

rule ldsc_mrkdwn_report_rna:
    # R produces 5 plots but only tracking the final plot here 
    input:   all = expand(PART_HERIT_DIR + "snRNAseq_LDSC_{CELL_TYPE}_{GWAS}.rds", CELL_TYPE = config["RNA_CELL_TYPES"], GWAS = config["LDSC_GWAS"]),
             top10 = expand(PART_HERIT_DIR + "snRNAseq_LDSC_{GWAS}_top10pc.tsv", GWAS = config["LDSC_GWAS"])    
    output:  RESULTS_DIR + "snRNAseq_LDSC_report.html"
    params:  markdown_file = MARKDOWN_DIR + "snRNAseq_LDSC_report.Rmd",
             out_dir = PART_HERIT_DIR,
             mrkdn_out_dir = RESULTS_DIR
    message: "Creating ldsc group plots for all regions and GWAS"
    log:     "logs/LDSR/snRNAseq.AllRegions.All.GWAS_partHerit.group.plots.log"
    shell:
             """
             export R_LIBS_USER=/scratch/c.c1477909/R/library
             module load libgit2/1.1.0
             module load pandoc/2.7.3
             /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
             scripts/R/snRNAseq_LDSC_report.R {params.markdown_file} {params.out_dir} {params.mrkdn_out_dir} 2> {log}
             
             """
             
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
