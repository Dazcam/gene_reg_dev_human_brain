#!/bin/bash --login
#SBATCH --job-name=scenic
#SBATCH --output=scenic_20G.out.%J
#SBATCH --error=scenic_20G.err.%J
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=10G
#SBATCH --qos=maxjobs500
#SBATCH --account=scw1641

echo -e "\nLoading singularity\n"

module load singularity

echo -e "Run GRNboost2 ... \n"

singularity run aertslab-pyscenic-0.10.0.sif \
  pyscenic grn \
    --num_workers 2 \
    -o FC_raw_counts.adjacencies.tsv \
    FC_raw_counts.txt \
    hs_hgnc_tfs.txt

echo -e "Get Regulons ... \n"

singularity run aertslab-pyscenic-0.10.0.sif \
  pyscenic ctx \
    FC_raw_counts.adjacencies.tsv \
    encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.feather \
    encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__500bp_up_and_100bp_down_tss.max.feather \
    --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname FC_raw_counts.txt \
    --mode "dask_multiprocessing" \
    --output regulons.csv \
    --num_workers 2

echo -e "Run AUC ... \n"

singularity run aertslab-pyscenic-0.10.0.sif \
  pyscenic aucell \
    FC_raw_counts.txt \
    regulons.csv \
    -o auc_mtx.csv \
    --num_workers 2


echo "Done."
