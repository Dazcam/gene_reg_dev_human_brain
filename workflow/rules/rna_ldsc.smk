
# -------------------------------------------------------------------------------------
#
#
#    Script for processing snATAC-seq FASTQ files with Cell Ranger (10X)
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "config.yaml"
report: "report/workflow.rst" # Experimented with this but not currently using

# ----------  SET VARIABLES  ----------
SCRATCH = "/scratch/c.c1477909/"
LDSC_DIR = SCRATCH + "ldsc/"
MAGMA_DIR = SCRATCH + "magma_celltyping/"
REF_DIR = LDSC_DIR + "reference_files/"
ANN_DIR = LDSC_DIR + "LDSR_annotation_files/RNA/"
GWAS_DIR = LDSC_DIR + "GWAS_for_ldsc/"
QUANTILES_DIR = MAGMA_DIR + "ldsc_gene_lists/"
PART_HERIT_DIR = LDSC_DIR + "part_herit_files/RNA/"

# -------------  RULES  ---------------
rule all:
    input:
#        expand(ANN_DIR + "snRNAseq.test.Q{QUANTILE}.{CHR}.annot.gz", QUANTILE = range(1,11), CHR = range(1,23))
#        expand(ANN_DIR + "snRNAseq.test.Q{QUANTILE}.{CHR}.l2.ldscore.gz", QUANTILE = range(1,11), CHR = range(1,23))
#        expand(PART_HERIT_DIR + "prtHrt_snRNAseq_Q{QUANTILE}_SCZ.results", QUANTILE = range(1,11))
#         expand(PART_HERIT_DIR + "prtHrt_snRNAseq_{CELL_TYPE}_Q{QUANTILE}_SCZ.results", CELL_TYPE = config["RNA_CELL_TYPES"], QUANTILE = range(0,11)), 
#         expand(PART_HERIT_DIR + "prtHrt_snRNAseq_{CELL_TYPE}_SCZ_summary.tsv", CELL_TYPE = config["RNA_CELL_TYPES"])
#         expand(PART_HERIT_DIR + "prtHrt_snRNAseq_{CELL_TYPE}_SCZ_plot.png", CELL_TYPE = config["RNA_CELL_TYPES"]) 
        expand(PART_HERIT_DIR + "Thal_ldsc_RNA_group_plot_lst.rds")

rule make_annot:
    input:   gene_set = QUANTILES_DIR + "{CELL_TYPE}_Q{QUANTILE}_genes.tsv",
             gene_coord = LDSC_DIR + "ldsc_gene_coords.tsv",
             bim_file = REF_DIR + "1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}.bim" 
    output:  ANN_DIR + "snRNAseq.{CELL_TYPE}.Q{QUANTILE}.{CHR}.annot.gz"
    conda:   "envs/ldsc.yml"    
    message: "Creating annotation files for snRNAseq: {wildcards.CELL_TYPE} Quantile {wildcards.QUANTILE}, Chr {wildcards.CHR} "
    log:     "logs/ldsc/make_annot.snRNAseq.{CELL_TYPE}.Q{QUANTILE}.Chr{CHR}.log"
    shell:
        """
        
        python ldsc/make_annot.py \
        --gene-set-file {input.gene_set} \
        --gene-coord-file {input.gene_coord} \
        --windowsize 100000 \
        --bimfile {input.bim_file} \
        --annot-file {output} 2> {log} \
        
        """

rule ldsr:
    input:   annot = ANN_DIR + "snRNAseq.{CELL_TYPE}.Q{QUANTILE}.{CHR}.annot.gz",
             bfile_folder = REF_DIR + "1000G_EUR_Phase3_plink",
             snps_folder = REF_DIR + "hapmap3_snps"
    output:  ANN_DIR + "snRNAseq.{CELL_TYPE}.Q{QUANTILE}.{CHR}.l2.ldscore.gz"
    conda:   "envs/ldsc.yml"
    params:  bfile = REF_DIR + "1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores = ANN_DIR + "snRNAseq.{CELL_TYPE}.Q{QUANTILE}.{CHR}",
             snps = REF_DIR + "w_hm3.snplist_rsIds"
    message: "Running LDSR Phase 3 for {wildcards.CELL_TYPE} Quantile {wildcards.QUANTILE} CHR {wildcards.CHR}" 
    log:     "logs/LDSR/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.Chr{CHR}_ldsc.log"
    shell:
        "python ldsc/ldsc.py --thin-annot --l2 --bfile {params.bfile} --ld-wind-cm 1 "
        "--annot {input.annot} --out {params.ldscores} --print-snps {params.snps} 2> {log}"

rule partitioned_heritability:
    input:   GWAS = GWAS_DIR + "SCZ_hg38_ldsc_ready.sumstats.gz",
             #GWAS = lambda wildcards: config['SUMSTATS'][wildcards.GWASsumstat],
             LDSR = expand(rules.ldsr.output, CELL_TYPE = config["RNA_CELL_TYPES"], QUANTILE = range(0,11), CHR = range(1,23))
    output:  PART_HERIT_DIR + "prtHrt_snRNAseq_{CELL_TYPE}_Q{QUANTILE}_SCZ.results"
    conda:   SCRATCH + "envs/ldsc.yml"
    params:  weights = REF_DIR + "weights_hm3_no_hla/weights.",
             baseline = REF_DIR + "baselineLD_v2.2_1000G_Phase3/baselineLD.",
             frqfile = REF_DIR + "1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = ANN_DIR + "snRNAseq.{CELL_TYPE}.Q{QUANTILE}.",
             out_file = PART_HERIT_DIR + "prtHrt_snRNAseq_{CELL_TYPE}_Q{QUANTILE}_SCZ"
    message: "Running Prt Hrt with {wildcards.CELL_TYPE} Q{wildcards.QUANTILE} and SCZ GWAS"
    log:     "logs/LDSR/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.SCZ_partHerit.log"
    shell:
        "python ldsc/ldsc.py --h2 {input.GWAS} --w-ld-chr {params.weights} "
        "--ref-ld-chr {params.baseline},{params.LD_anns} --overlap-annot "
        "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"


rule create_partHerit_summary:
    input:   expand(PART_HERIT_DIR + "prtHrt_snRNAseq_{CELL_TYPE}_Q{QUANTILE}_SCZ.results", CELL_TYPE = config["RNA_CELL_TYPES"], QUANTILE = range(0,11))
    output:  PART_HERIT_DIR + "prtHrt_snRNAseq_{CELL_TYPE}_SCZ_summary.tsv"
    message: "Creating summary file for {wildcards.CELL_TYPE} and SCZ GWAS"
    params:  dir = PART_HERIT_DIR
    log:     "logs/LDSR/snRNAseq.{CELL_TYPE}.SCZ_partHerit.summary.log"
    shell:
        """
        
        head -1 {params.dir}prtHrt_snRNAseq_Cer-RG-1_Q1_SCZ.results > {output} 
        grep L2_1 {params.dir}prtHrt_snRNAseq_{wildcards.CELL_TYPE}_Q*_SCZ.results >> {output}

        """

        # Last two rules were run in seprate  test script need to test whole thing is working
        
rule create_ldsc_plots:
    input:   PART_HERIT_DIR + "prtHrt_snRNAseq_{CELL_TYPE}_SCZ_summary.tsv"
    output:  rds = PART_HERIT_DIR + "prtHrt_snRNAseq_{CELL_TYPE}_SCZ.rds"
             plot = PART_HERIT_DIR + "prtHrt_snRNAseq_{CELL_TYPE}_SCZ.png"
    message: "Creating ldsc plots for {wildcards.CELL_TYPE} and SCZ GWAS"
    log:     "logs/LDSR/snRNAseq.{CELL_TYPE}.SCZ_partHerit.plots.log"
    shell:
             """
             
             export R_LIBS_USER=/scratch/c.c1477909/R/library
             module load libgit2/1.1.0
             /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
             scripts/R/scRNAseq_LDSC_create_plots.R {wildcards.CELL_TYPE} {input} {output.rds} {output.plot} 2> {log}

             """

rule create_ldsc_group_plots:
    # R produces 5 plots but only tracking the final plot here 
    input:   expand(PART_HERIT_DIR + "prtHrt_snRNAseq_{CELL_TYPE}_SCZ.rds", CELL_TYPE = config["RNA_CELL_TYPES"])
    output:  PART_HERIT_DIR + "Thal_ldsc_RNA_group_plot_lst.rds"
    params:  out_dir = PART_HERIT_DIR 
    message: "Creating ldsc group plots for all regions and SCZ GWAS"
    log:     "logs/LDSR/snRNAseq.AllRegions.SCZ_partHerit.group.plots.log"
    shell:
             """

             export R_LIBS_USER=/scratch/c.c1477909/R/library
             module load libgit2/1.1.0
             /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
             scripts/R/scRNAseq_LDSC_create_group_plots.R {params.out_dir}  2> {log}
             
             """
        
        
        
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
