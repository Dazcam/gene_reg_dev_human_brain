# -------------------------------------------------------------------------------------
#
#    Script to standarise and munge GWAS sumstats for MAGMA and LDSC analysis
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARMAS  ----------
configfile: "config.yaml"


# ----------  SET VARIABLES  ----------
SCRATCH =  "/scratch/c.c1477909/"
GWAS_DIR = SCRATCH + "GWAS_sumstats/"
MAGMA_GWAS_DIR =  SCRATCH + "magma_celltyping/GWAS_for_magma/"
LDSC_GWAS_DIR = SCRATCH + "ldsc/GWAS_for_ldsc/"


# -------------  RULES  --------------
rule all:
    input:
        expand(MAGMA_GWAS_DIR + "{GWAS}_hg19_magma_ready.sumstats.tsv", GWAS=config['SUMSTATS']),
        expand(MAGMA_GWAS_DIR + "{GWAS}_hg38_magma_ready.sumstats.tsv", GWAS=config['SUMSTATS']),
#        expand(LDSC_GWAS_DIR + "{GWAS}_hg19_ldsc_ready.sumstats.gz", GWAS=config['SUMSTATS']) 
        expand(LDSC_GWAS_DIR + "{GWAS}_hg38_ldsc_ready.sumstats.gz", GWAS=config['SUMSTATS']) 

rule standardise_sumstats:
    # Standardises sumstats: SNP, CHR. BP, PVAL, A1, A2 + additional GWAS dependant cols
    input:   lambda wildcards: config['SUMSTATS'][wildcards.GWAS]
    output:  GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_sumstats.tsv"
    message: "Formatting {input}"
    log:     SCRATCH + "logs/munge_sumstats/{GWAS}_standardise_sumstats.log"
    shell:
             """ 
       
             python python_convert/sumstats.py csv \
             --sumstats {input} \
       	     --out {output} --force --auto --head 5 \
             --log {log}

             """

rule add_z_score:
    # Adds z-scores to GWAS sumstats lacking (SCZ and BPD)
    input:   GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_sumstats.tsv"
    output:  GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_withZ_sumstats.tsv"
    message: "Adding Z score to {input}"
    log:     SCRATCH + "logs/munge_sumstats/{GWAS}_addZscore.log"
    shell:
             """

             python python_convert/sumstats.py zscore \
             --sumstats {input} \
             --out {output} --force \
             --log {log} \
             --a1-inc

             """

rule final_standardisation_edits:
    # SCZ - removing SNP with 13 cols, adding N, sort cols to match rest of sumstats files
    # BPD - Adding N
    # Format all GWAS to SNP, CHR, BP, PVAL, A1, A2, N, Z, OR, BETA/INFO SE
    input:   GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_withZ_sumstats.tsv"
    output:  GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_withZv2_sumstats.tsv"
    message: "Adding Z score to {input}"
    log:     SCRATCH + "logs/munge_sumstats/{GWAS}_SCZ_edits.log"
    run:
             if "SCZ" in wildcards.GWAS:
             
                 shell("""
 
                 sed -i '/rs148878475/d' {input}
                 awk '{{s=(NR==1)?"N":"161405";$0=$0 OFS s}}1' {input} |\
                 awk '{{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$12"\t"$11"\t"$7"\t"$9"\t"$8}}' >\
                 {output} 2> {log} 

                 """)

             elif "BPD" in wildcards.GWAS:

                 shell("""
 
                 awk '{{s=(NR==1)?"N":"413466";$0=$0 OFS s}}1' {input} |\
                 awk '{{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$12"\t"$11"\t"$7"\t"$9"\t"$8}}' >\
                 {output}

                 """)

             else:

                 shell("cp {input} {output}")


rule sumstats_to_bed:
     # Convert GWAS sumstats to bed format to run liftover to hg38
    input:   GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_withZv2_sumstats.tsv"
    output:  GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_sumstats.bed"
    message: "Creating {input} bed file"
    log:     SCRATCH + "logs/munge_sumstats/{GWAS}_create_bed.log"
    shell:
             """

             awk '{{print "chr"$2"\t"$3"\t"($3 + 1)"\t"$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}}' {input} |\
             sed -e '1s/^/#/' > {output} 2> {log}

             """


rule lift_over:
    # Lift GWAS data from hg19 to hg38
    input:   GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_sumstats.bed"
    output:  GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg38_sumstats.bed"
    message: "Lifting {input} to hg38"
    log:     SCRATCH + "logs/munge_sumstats/{GWAS}_liftOver.log"
    params:  chain = GWAS_DIR + "hg19ToHg38.over.chain.gz",
             unlift = GWAS_DIR + "{wildcards.GWAS}_hg38_unlifted.bed"
    shell:
             """

             ./liftover/liftOver {input} {params.chain} {output} {params.unlift} -bedPlus=4 2> {log}

             """


rule bed_to_sumstats:
    # Convert hg38 bed file to sumstats format
    input:   bed = GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg38_sumstats.bed",
             header = GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_withZv2_sumstats.tsv"
    output:  GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg38_sumstats.tsv"
    message: "Revert {input} to sumstats"
    log:     SCRATCH + "logs/munge_sumstats/{GWAS}_revert2sumstats.log"
    shell:
             """
             
             cut --complement -d$'\\t' -f3 {input.bed} |\
             sed -e 's/^chr//g' |\
             awk '{{print $3"\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}}' |\
             sed -e "1s/^/$(head -n1 {input.header})\\n/" > {output} 2> {log} 

             """

rule sumstats_for_magma:
    # MAGMA needs PVAL to be labeled P
    input:   hg19 = GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg19_withZv2_sumstats.tsv",
             hg38 = GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg38_sumstats.tsv"
    output:  hg19 = MAGMA_GWAS_DIR + "{GWAS}_hg19_magma_ready.sumstats.tsv",
             hg38 = MAGMA_GWAS_DIR + "{GWAS}_hg38_magma_ready.sumstats.tsv"
    message: "Munging sumstats for {input} for magma compatibility"
    log:     SCRATCH + "logs/munge_sumstats/{GWAS}_sumstats_for_magma.log"
    shell:   
             """

             sed 's/PVAL/P/g' {input.hg19} > {output.hg19}
             sed 's/PVAL/P/g' {input.hg38} > {output.hg38}

             """


rule sumstats_for_ldsc_hg38:
    input:   snps = SCRATCH + "ldsc/reference_files/w_hm3.snplist",
             gwas = GWAS_DIR + "GWAS_sumstats_standardised/{GWAS}_hg38_sumstats.tsv"
    output:  LDSC_GWAS_DIR + "{GWAS}_hg38_ldsc_ready.sumstats.gz"
    conda:   "envs/ldsc.yml"
    message: "Munging sumstats for {input.gwas} for ldsc compatibility"
    params:  out = SCRATCH + GWAS_DIR + "GWAS_for_ldsc/{GWAS}_hg38_ldsc_ready"
    log:     "logs/munge_sumstats/{GWAS}_sumstats_for_ldsc.log"
    shell:
        """

        python ldsc/munge_sumstats.py --sumstats {input.gwas} \
        --merge-alleles {input.snps} \
        --out {params.out} \
        --a1-inc 2> {log}

        """


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------


