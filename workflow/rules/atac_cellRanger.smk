# -------------------------------------------------------------------------------------
#
#
#    Script for processing snATAC-seq FASTQ files with Cell Ranger (10X)
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "config_ATAC.yaml"


# ----------  SET VARIABLES  ----------
shell("cellranger-atac -V;")

# -------------  RULES  ---------------

rule all:
    input:
       expand("{SAMPLE}.stamp", SAMPLE=config["SAMPLES"]),
#       "scRNA_CR_output.csv"

rule CR_cnt_ATAC:
    # Diminishing returns > 128Gs
    output: "{SAMPLE}.stamp" # Need .stamp and touch {output} in shell command as both SM and CR want to mkdir
    params: fastq_dir=config["FASTQ_DIR"],
            reference=config["REFERENCE_ATAC"]
    log:    "logs/{SAMPLE}.log"
    shell:
            """
            cellranger-atac count --id={wildcards.SAMPLE} \
            --fastqs={params.FASTQ_DIR} \
            --sample={wildcards.SAMPLE} \
            --reference={params.REFERENCE} \
            --localcores=32 \
            --localmem=128 2> {log}
            
            touch {output}
            """

rule CR_summary_doc:
    input:  csv = lambda wildcards: "{wildcards.SAMPLE}/outs/metrics_summary.csv", 
            stamp = lambda wildcards: "{wildcards.SAMPLE}.stamp"
    output: "scRNA_CR_output.csv"
    shell:
            """
            awk 'NR == 1' {input.csv} | sed -e 's/^/Sample,/' > {output}

            while read sample; do

              awk 'NR == 2' {input.csv} | sed -e "s/^/${sample},/" >> {output}
 
            done <file_list
            """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
