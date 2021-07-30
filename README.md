# Gene regulation in the developing human brain

Code repository for Gene regulation in the devloping human brain 

***

## Samples


## Alignment

snRNA-seq data aligned to Hg38 genome build
snATAC-seq data aligned to Hg38 genome build

***

## Data availability

## Pipeline

This snakemake pipeline

[Distribution and reproducibility](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)

## Data processing

1. Cell Ranger - process FASTQ files
2. RNA-seq - Seurat
3. ATAC-seq - ArchR

4. MAGMA - Reference data - 1000 genomes (Hg19)
5. sLDSC - refernce data - 1000 geneomes (Hg19)

6. Cluster stability (RNA) - SVM - [scRNAseq Benchmark](https://github.com/tabdelaal/scRNAseq_Benchmark/tree/snakemake_and_docker)

# GWAS

