# -------------------------------------------------------------------------------------
#
#
#    Script for running LDSC on snATAC-seq data - on SCZ GWAS only at this point
#
#
# -------------------------------------------------------------------------------------

# Need rule for markown

# ---------  SET SMK PARAMS  ----------
configfile: "config/config.yaml"

# ----------  SET VARIABLES  ----------
SCRATCH = "/scratch/c.c1477909/"
LDSC_DIR = SCRATCH + "ldsc/"
MAGMA_DIR = SCRATCH + "magma_celltyping/"
REF_DIR = LDSC_DIR + "reference_files/"
ANN_DIR = LDSC_DIR + "LDSR_annotation_files/ATAC/"
GWAS_DIR = LDSC_DIR + "GWAS_for_ldsc/"
QUANTILES_DIR = MAGMA_DIR + "ldsc_gene_lists/"
PART_HERIT_DIR = LDSC_DIR + "part_herit_files/ATAC/"

# -------------  RULES  ---------------configfile: "config.yaml"
rule all:
    input:
#        expand(PART_HERIT_DIR + "prtHrt_snATACseq_{CELL_TYPE}_SCZ.results", CELL_TYPE = config['CELL_TYPES'])
        expand(RESULTS_DIR + "prtHrt_ATACseq_SCZ_summary.tsv")

rule annot2bed:
    input:   folder = REF_DIR + "baselineLD_v2.2_1000G_Phase3"
    params:  file = REF_DIR + "baselineLD_v2.2_1000G_Phase3/baselineLD.{CHR}.annot.gz"
    output:  ANN_DIR + "baseline.{CHR}_no_head.bed"
    message: "Preparing annotation files for {wildcards.CHR}"
    shell:
        "zcat {params.file} | tail -n +2 |awk -v OFS=\"\t\" '{{print \"chr\"$1, $2-1, $2, $3, $4}}' "
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
             LDSR = expand(rules.ldsr.output, CELL_TYPE = config['CELL_TYPES'], CHR = range(1,23))
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
    input:   expand(PART_HERIT_DIR + "prtHrt_snATACseq_{CELL_TYPE}_SCZ.results", CELL_TYPE = config["ATAC_CELL_TYPES"])
    output:  RESULTS_DIR + "prtHrt_ATACseq_SCZ_summary.tsv"
    message: "Creating summary file for SCZ GWAS"
    params:  dir = PART_HERIT_DIR,
             cell_types = config["ATAC_CELL_TYPES"]
    log:     "logs/LDSR/snRNAseq.SCZ_partHerit.summary.log"
    shell:
	         """

	         head -1 {params.dir}prtHrt_snATACseq_InN.1.500bpExt.FDR.0.01.ge_SCZ.results > {output}

             File=atac_celltypes.tsv
             Lines=$(cat $File)
             for Line in $Lines
             do
	         grep L2_1 /scratch/c.c1477909/ldsc/part_herit_files/prtHrt_snATACseq_"$Line"_SCZ.results | sed "s/L2_1/$Line/g" >> /scratch/c.c1477909/results/prtHrt_ATACseq_SCZ_summary.tsv
             done

	         """
    
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
