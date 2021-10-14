# -------------------------------------------------------------------------------------
#
#
#    Script for running comparing fetal and adult cell types
#
#
# -------------------------------------------------------------------------------------


# ---------  SET SMK PARAMS  ----------
configfile: "config/config.yaml"

# ----------  SET VARIABLES  ----------
SCRATCH = "/scratch/c.c1477909/"
REF_DIR = SCRATCH + "ldsc/reference_files/"
RESULTS_DIR = SCRATCH + "results/fetal_vs_adult/"
ANN_DIR = RESULTS_DIR + "LDSR_annotation_files/ATAC/"
PART_HERIT_DIR = RESULTS_DIR + "LDSR_part_herit_files/ATAC/"
GWAS_DIR = SCRATCH + "ldsc/GWAS_for_ldsc/"
MARKDOWN_DIR = SCRATCH + "markdown/"
LIFT_OVER_DIR = SCRATCH + "liftover/"
INPUT_DIR = SCRATCH + "input_files"
MAGMA_DIR = SCRATCH + "magma_celltyping/GWAS_for_magma/MAGMA_Files/SCZ_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN/"

# -------------  RULES  ---------------

rule all:
    input:
	expand(RESULTS_DIR + "{CELL_TYPE}.gsa.out", CELL_TYPE = config['FETAL_VS_ADULT_CELL_TYPES'])


rule fetal_vs_adult_magma:
    input:   gene_list = RESULTS_DIR + "{CELL_TYPE}.tsv",
             scz_magma = MAGMA_DIR + "SCZ_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.raw"
    output:  RESULTS_DIR + "{CELL_TYPE}.gsa.out"
    message: "Running magma for {wildcards.CELL_TYPE} to hg38"
    log:     SCRATCH + "logs/fetal_vs_adult/snRNAseq.magma.{CELL_TYPE}.log"
    shell:
             """

             magma --gene-results {input.scz_magma} --set-annot {input.gene_list} col=1,2 --out {output}

             """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
