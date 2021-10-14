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
PART_HERIT_DIR = RESULTS_DIR + "LDSR_part_herit_files/ATAC/"
BED_DIR = RESULTS_DIR + "bed_files_for_LDSC_ATAC/"
GWAS_DIR = SCRATCH + "ldsc/GWAS_for_ldsc/"
MARKDOWN_DIR = SCRATCH + "markdown/"
LIFT_OVER_DIR = SCRATCH + "liftover/"

# -------------  RULES  ---------------
rule all:
    input:
#        expand(RESULTS_DIR + "snATACseq_LDSC_summary_{GWAS}.tsv", GWAS = config['LDSC_GWAS']),
        expand(RESULTS_DIR + "snATACseq_LDSC_report.html") 
#         expand(PART_HERIT_DIR + "snATACseq_LDSC_{CELL_TYPE}_{GWAS}.results", CELL_TYPE = config['ATAC_CELL_TYPES'], GWAS = config['LDSC_GWAS'])

rule annot2bed:
    input:   folder = REF_DIR + "baselineLD_v2.2_1000G_Phase3"
    params:  file = REF_DIR + "baselineLD_v2.2_1000G_Phase3/baselineLD.{CHR}.annot.gz"
    output:  ANN_DIR + "baseline.{CHR}_no_head.bed"
    message: "Preparing annotation files for {wildcards.CHR}"
    shell:
             "zcat {params.file} | tail -n +2 | awk -v OFS=\"\t\" '{{print \"chr\"$1, $2-1, $2, $3, $4}}' "
             "| sort -k1,1 -k2,2n > {output}"

rule lift_over:
    input:   mybed = BED_DIR + "{CELL_TYPE}.hg38.bed",
             chain_file = LIFT_OVER_DIR + "hg38ToHg19.over.chain.gz"
    output:  mybed = BED_DIR + "{CELL_TYPE}.hg19.bed"
    message: "Lifting {input.mybed} to hg38"
    log:     SCRATCH + "logs/LDSR/scATACseq.{CELL_TYPE}_liftOver.log"
    params:  unlift = BED_DIR + "{CELL_TYPE}_hg38_unlifted.bed"
    shell:
             """

             ./liftover/liftOver {input.mybed} {input.chain_file} {output} {params.unlift} 2> {log}

             """

rule intersect_mybed:
    input:   annot = rules.annot2bed.output,
             mybed = BED_DIR + "{CELL_TYPE}.hg19.bed"
    output:  ANN_DIR + "snATACseq.{CELL_TYPE}.{CHR}.annot.gz"
    params:  out = ANN_DIR + "snATACseq.{CELL_TYPE}.{CHR}.annot"
    message: "Creating LDSR annotation files for {wildcards.CHR}"
    shell:
             "module load bedtools; " 
             "echo -e \"CHR\tBP\tSNP\tCM\tANN\" > {params.out}; "
             "bedtools intersect -a {input.annot} -b {input.mybed} -c | "
             "sed 's/^chr//g' | awk -v OFS=\"\t\" '{{print $1, $2, $4, $5, $6}}' >> {params.out}; "
             "gzip {params.out}"

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

rule partitioned_heritability:
    input:   SUMSTATS = GWAS_DIR + "{GWAS}_hg19_ldsc_ready.sumstats.gz",
             LDSR = expand(rules.ldsr.output, CELL_TYPE = config['ATAC_CELL_TYPES'], CHR = range(1,23))
    output:  PART_HERIT_DIR + "snATACseq_LDSC_{CELL_TYPE}_{GWAS}.results"
    conda:   SCRATCH + "envs/ldsc.yml"
    params:  weights = REF_DIR + "weights_hm3_no_hla/weights.",
             baseline = REF_DIR + "baselineLD_v2.2_1000G_Phase3/baselineLD.",
             frqfile = REF_DIR + "1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = ANN_DIR + "snATACseq.{CELL_TYPE}.",
             out_file = PART_HERIT_DIR + "snATACseq_LDSC_{CELL_TYPE}_{GWAS}"
    message: "Running Prt Hrt with {wildcards.CELL_TYPE} and {wildcards.GWAS} GWAS"
    log:     "logs/LDSR/snATACseq.{CELL_TYPE}.{GWAS}_partHerit.log"
    shell:
             "python ldsc/ldsc.py --h2 {input.SUMSTATS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule create_partHerit_summary:
    # Requires list of snATACseq cell types in atac_celltypes.tsv 
    input:   expand(PART_HERIT_DIR + "snATACseq_LDSC_{CELL_TYPE}_{GWAS}.results", CELL_TYPE = config["ATAC_CELL_TYPES"], GWAS = config['LDSC_GWAS'])
    output:  RESULTS_DIR + "snATACseq_LDSC_summary_{GWAS}.tsv"
    message: "Creating summary file for {wildcards.GWAS} GWAS"
    params:  ph_dir = PART_HERIT_DIR,
             results_dir = RESULTS_DIR,
             cell_types = SCRATCH + "config/atac_celltypes.tsv"
    log:     "logs/LDSR/snATACseq.{GWAS}_partHerit.summary.log"
    shell:
             """

             head -1 {params.ph_dir}snATACseq_LDSC_cer.ExN.1_BPD.results > {output}

             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_1 {params.ph_dir}snATACseq_LDSC_"$Line"_{wildcards.GWAS}.results | sed "s/L2_1/$Line/g" >> {output}
             done

             """

rule ldsc_mrkdwn_report:
    input:   summary = RESULTS_DIR + "snATACseq_LDSC_summary_SCZ.tsv",
             markdown = MARKDOWN_DIR + "snATACseq_LDSC_report.Rmd"
    output:  RESULTS_DIR + "snATACseq_LDSC_report.html"
    message: "Creating LDSC Rmarkdown report"
    params:  out_dir = RESULTS_DIR
    log:     "logs/LDSR/snATACseq.LDSC.report.log"
    shell:
             """
             
            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3 
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/snATACseq_LDSC_report.R {input.markdown} {params.out_dir} 2> {log}

             """
    
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
