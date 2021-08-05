# -------------------------------------------------------------------------------------
#
#
#    Script for running LDSC on snATAC-seq data - on SCZ GWAS only at this point
#
#
# -------------------------------------------------------------------------------------

# Need rule for markdown

# ---------  SET SMK PARAMS  ----------
configfile: "config/config.yaml"

# ----------  SET VARIABLES  ----------
SCRATCH = "/scratch/c.c1477909/"
REF_DIR = SCRATCH + "ldsc/reference_files/"
ANN_DIR = SCRATCH + "ldsc/LDSR_annotation_files/ATAC/"
BED_DIR = SCRATCH + "bed_files/"
GWAS_DIR = SCRATCH + "ldsc/GWAS_for_ldsc/"
PART_HERIT_DIR = SCRATCH + "ldsc/part_herit_files/ATAC/"
RESULTS_DIR = SCRATCH + "results/"
MARKDOWN_DIR = SCRATCH + "markdown/"

# -------------  RULES  ---------------
rule all:
    input:
#        expand(PART_HERIT_DIR + "prtHrt_snATACseq_{CELL_TYPE}_SCZ.results", CELL_TYPE = config['ATAC_CELL_TYPES'])
        expand(RESULTS_DIR + "prtHrt_ATACseq_SCZ_summary.tsv")

rule annot2bed:
    input:   folder = REF_DIR + "baselineLD_v2.2_1000G_Phase3"
    params:  file = REF_DIR + "baselineLD_v2.2_1000G_Phase3/baselineLD.{CHR}.annot.gz"
    output:  ANN_DIR + "baseline.{CHR}_no_head.bed"
    message: "Preparing annotation files for {wildcards.CHR}"
    shell:
             "zcat {params.file} | tail -n +2 | awk -v OFS=\"\t\" '{{print \"chr\"$1, $2-1, $2, $3, $4}}' "
             "| sort -k1,1 -k2,2n > {output}"

rule intersect_mybed:
    input:   annot = rules.annot2bed.output,
             mybed = BED_DIR + "{CELL_TYPE}.bed"
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
    input:   GWAS = GWAS_DIR + "SCZ_hg38_ldsc_ready.sumstats.gz",
             #GWAS = lambda wildcards: config['SUMSTATS'][wildcards.GWASsumstat],
             LDSR = expand(rules.ldsr.output, CELL_TYPE = config['ATAC_CELL_TYPES'], CHR = range(1,23))
    output:  PART_HERIT_DIR + "prtHrt_snATACseq_{CELL_TYPE}_SCZ.results"
    conda:   SCRATCH + "envs/ldsc.yml"
    params:  weights = REF_DIR + "weights_hm3_no_hla/weights.",
             baseline = REF_DIR + "baselineLD_v2.2_1000G_Phase3/baselineLD.",
             frqfile = REF_DIR + "1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = ANN_DIR + "snATACseq.{CELL_TYPE}.",
             out_file = PART_HERIT_DIR + "prtHrt_snATACseq_{CELL_TYPE}_SCZ"
    message: "Running Prt Hrt with {wildcards.CELL_TYPE} and SCZ GWAS"
    log:     "logs/LDSR/snATACseq.{CELL_TYPE}.SCZ_partHerit.log"
    shell:
             "python ldsc/ldsc.py --h2 {input.GWAS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule create_partHerit_summary:
    # Requires list of snATACseq cell types in atac_celltypes.tsv 
    input:   expand(PART_HERIT_DIR + "prtHrt_snATACseq_{CELL_TYPE}_SCZ.results", CELL_TYPE = config["ATAC_CELL_TYPES"])
    output:  RESULTS_DIR + "prtHrt_ATACseq_SCZ_summary.tsv"
    message: "Creating summary file for SCZ GWAS"
    params:  ph_dir = PART_HERIT_DIR,
             results_dir = RESULTS_DIR,
             cell_types = SCRATCH + "config/atac_celltypes.tsv"
    log:     "logs/LDSR/snATACseq.SCZ_partHerit.summary.log"
    shell:
             """

             head -1 {params.ph_dir}prtHrt_snATACseq_InN.1.500bpExt.FDR.0.01.ge_SCZ.results > {output}

             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_1 {params.ph_dir}prtHrt_snATACseq_"$Line"_SCZ.results | sed "s/L2_1/$Line/g" >> {output}
             done

             """

rule ldsc_mrkdwn_report:
    input:   summary = RESULTS_DIR + "prtHrt_ATACseq_SCZ_summary.tsv",
             markdown = MARKDOWN_DIR + "snATACseq_LDSC_report.Rmd"
    output:  RESULTS_DIR + "snATACseq_LDSC_report.html"
    message: "Creating LDSC Rmarkdown reportfor SCZ GWAS"
    params:  out_dir = RESULTS_DIR
             cell_types = config["ATAC_CELL_TYPES"]
    log:     "logs/LDSR/snATACseq.SCZ_partHerit.markdown.log"
    shell:
             """
             
            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/R/snATACseq_LDSC_report.R {input.summary} {input.markdown} {params.out_dir} 2> {log}

             """
    
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
