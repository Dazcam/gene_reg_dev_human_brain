# -------------------------------------------------------------------------------------
#
#
#    Script for running LDSC on snATAC-seq data
#
#
# -------------------------------------------------------------------------------------


# ---------  SET SMK PARAMS  ----------
configfile: "config/config.yaml"

# ----------  SET VARIABLES  ----------
SCRATCH = "/scratch/c.c1477909/"
REF_DIR = SCRATCH + "ldsc/reference_files/"
RESULTS_DIR = SCRATCH + "results/"
ANN_DIR = RESULTS_DIR + "LDSR_annotation_files/ATAC/"
PART_HERIT_DIR = RESULTS_DIR + "LDSR_part_herit_files/ATAC/conditional/"
BED_DIR = RESULTS_DIR + "bed_files_for_LDSC_ATAC/"
GWAS_DIR = SCRATCH + "ldsc/GWAS_for_ldsc/"
MARKDOWN_DIR = SCRATCH + "markdown/"
LIFT_OVER_DIR = SCRATCH + "liftover/"
REPORT_DIR = RESULTS_DIR + "reports/ATAC/"

# -------------  RULES  ---------------
rule all:
    input:
        expand(PART_HERIT_DIR + "snATACseq_LDSC_{CELL_TYPE}.vs.fc.ExN_{GWAS}.results", CELL_TYPE = config['ATAC_CELL_TYPES_COND'], GWAS = config['LDSC_GWAS']), 
        expand(RESULTS_DIR + "snATACseq_LDSC_summary_{GWAS}.vs.fc.ExN.tsv", GWAS = config['LDSC_GWAS']),
        REPORT_DIR + "snATACseq_LDSC_conditional_report.html"
 
rule ldsr:
    input:   annot = ANN_DIR + "snATACseq.{CELL_TYPE}.{CHR}.annot.gz",
             bfile_folder = REF_DIR + "1000G_EUR_Phase3_plink",
             snps_folder = REF_DIR + "hapmap3_snps"
    output:  ANN_DIR + "snATACseq.{CELL_TYPE}.{CHR}.l2.ldscore.gz"
    conda:   "envs/ldsc.yml"
    params:  bfile = REF_DIR + "1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores = ANN_DIR + "snATACseq.{CELL_TYPE}.{CHR}",
             snps = REF_DIR + "w_hm3.snplist_rsIds"
    message: "Running LDSR Phase3 {wildcards.CHR}" 
    log:     "logs/LDSR/snATACseq.{CELL_TYPE}.{CHR}_ldsc.log"
    shell:
             "python ldsc/ldsc.py --l2 --bfile {params.bfile} --ld-wind-cm 1 "
             "--annot {input.annot} --out {params.ldscores} --print-snps {params.snps} 2> {log}"

rule partitioned_heritability_conditional:
    input:   SUMSTATS = GWAS_DIR + "{GWAS}_hg19_ldsc_ready.sumstats.gz",
             LDSR = expand(rules.ldsr.output, CELL_TYPE = config['ATAC_CELL_TYPES_COND'], CHR = range(1,23))
    output:  PART_HERIT_DIR + "snATACseq_LDSC_{CELL_TYPE}.vs.fc.ExN_{GWAS}.results"
    conda:   SCRATCH + "envs/ldsc.yml"
    params:  weights = REF_DIR + "weights_hm3_no_hla/weights.",
             baseline = REF_DIR + "baselineLD_v2.2_1000G_Phase3/baselineLD.",
             frqfile = REF_DIR + "1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = ANN_DIR + "snATACseq.{CELL_TYPE}.",
             out_file = PART_HERIT_DIR + "snATACseq_LDSC_{CELL_TYPE}.vs.fc.ExN_{GWAS}",
             condition = ANN_DIR + "snATACseq.fc.ExN.",
    message: "Running Prt Hrt with {wildcards.CELL_TYPE} and {wildcards.GWAS} GWAS"
    log:     "logs/LDSR/snATACseq.{CELL_TYPE}.vs.fc.ExN.{GWAS}_partHerit.log"
    shell:
             "python ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.condition},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule create_partHerit_summary_conditional:
    # Requires list of snATACseq cell types in atac_celltypes_conditional.tsv 
    input:   expand(PART_HERIT_DIR + "snATACseq_LDSC_{CELL_TYPE}.vs.fc.ExN_{GWAS}.results", CELL_TYPE = config["ATAC_CELL_TYPES_COND"], GWAS = config['LDSC_GWAS'])
    output:  RESULTS_DIR + "snATACseq_LDSC_summary_{GWAS}.vs.fc.ExN.tsv"
    message: "Creating summary file for {wildcards.GWAS} GWAS .vs.fc.ExN"
    params:  ph_dir = PART_HERIT_DIR,
             results_dir = RESULTS_DIR,
             cell_types = SCRATCH + "config/atac_celltypes_conditional.tsv"
    log:     "logs/LDSR/snATACseq.{GWAS}_partHerit.vs.fc.ExN.summary.log"
    shell:
             """

             head -1 {params.ph_dir}snATACseq_LDSC_fc.MG.vs.fc.ExN_MDD.results > {output}
             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_2 {params.ph_dir}snATACseq_LDSC_"$Line".vs.fc.ExN_{wildcards.GWAS}.results | sed "s/L2_2/$Line/g" >> {output}
             done

             """

rule ldsc_mrkdwn_report:
    input:   summary = expand(RESULTS_DIR + "snATACseq_LDSC_summary_{GWAS}.vs.fc.ExN.tsv", GWAS = config['LDSC_GWAS']),
             markdown = MARKDOWN_DIR + "snATACseq_LDSC_conditional_report.Rmd"
    output:  REPORT_DIR + "snATACseq_LDSC_conditional_report.html"
    message: "Creating LDSC conditional Rmarkdown report"
    params:  out_dir = REPORT_DIR
    log:     "logs/LDSR/snATACseq.LDSC.conditional.report.log"
    shell:
             """
             
             export R_LIBS_USER=/scratch/c.c1477909/R/library
             module load libgit2/1.1.0
             module load pandoc/2.7.3 
             /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
             scripts/R/snATACseq_LDSC_conditional_report.R {input.markdown} {params.out_dir} 2> {log}

             """
    
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
