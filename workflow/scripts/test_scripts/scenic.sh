#!/bin/bash --login
#SBATCH --job-name=scenic
#SBATCH --output=scenic.out.%J
#SBATCH --error=scenic.err.%J
#SBATCH --ntasks=24 
#SBATCH --mem-per-cpu=10G
#SBATCH --qos=maxjobs500
#SBATCH --account=scw1641

echo -e "\nLoading singularity\n"

module load singularity

echo -e "Run GRNboost2 ... \n"

singularity run aertslab-pyscenic-0.10.0.sif \
  pyscenic grn \
    --num_workers 24 \
    -o FC_raw_counts.adjacencies.tsv \
    FC_raw_counts.tsv \
    hs_hgnc_tfs.txt

echo -e "Get Regulons ... \n"

singularity run aertslab-pyscenic-0.10.0.sif \
  pyscenic ctx \
    FC_raw_counts.adjacencies.tsv \
    encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.feather \
    --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname FC_raw_counts.tsv \
    --mode "dask_multiprocessing" \
    --output regulons.csv \
    --num_workers 24

echo -e "Run AUC ... \n"

singularity run aertslab-pyscenic-0.10.0.sif \
  pyscenic aucell \
    FC_raw_counts.tsv \
    regulons.csv \
    -o auc_mtx.csv \
    --num_workers 24

echo "Done."
