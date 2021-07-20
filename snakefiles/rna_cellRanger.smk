# -------------------------------------------------------------------------------------
#
#
#    Script for processing snRNA-seq FASTQ files with Cell Ranger (10X)
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "config.yaml"


# ----------  SET VARIABLES  ----------
shell("cellranger -V;")


# -------------  RULES  ---------------
rule all:
    input:
       expand("{SAMPLE}.stamp", SAMPLE=config["SAMPLES_RNA"]),
#       "scRNA_CR_output.csv"

rule CR_cnt_RNA:
    output: "{SAMPLE}.stamp" # Need .stamp and touch {output} in shell command as both SM and CR want to mkdir
    params: fastq_dir=config["FASTQ_DIR"],
            transcriptome=config["REFERENCE_RNA"]
    log:    "logs/{SAMPLE}.log"
    shell:
            """
            cellranger count --id={wildcards.SAMPLE} \
            --fastqs={params.fastq_dir} \
            --sample={wildcards.SAMPLE} \
            --transcriptome={params.transcriptome} \
            --include-introns \
            --chemistry=SC3Pv3 \
            --localcores=32 \
            --localmem=230 2> {log}
            
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
