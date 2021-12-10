#!/bin/bash --login
#SBATCH --job-name=scenic
#SBATCH --output=CLI_scenic.out.%J
#SBATCH --error=CLI_scenic.err.%J
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=10G
#SBATCH --qos=maxjobs500
#SBATCH --account=scw1641

echo -e "\nLoading pyscenic environment ... \n"

source activate pyscenic

echo -e "Run GRNboost2 ... \n"

pyscenic grn \
  --num_workers 20 \
  -o FC_raw_counts.adjacencies_CLI.tsv \
  FC_raw_counts.tsv \
  hs_hgnc_tfs.txt

echo -e "Get Regulons ... \n"

pyscenic ctx \
  FC_raw_counts.adjacencies_CLI.tsv \
  encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.feather \
  encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__500bp_up_and_100bp_down_tss.max.feather \
  --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname FC_raw_counts.tsv \
  --mode "dask_multiprocessing" \
  --output regulons_CLI.csv \
  --num_workers 20

echo -e "Run AUC ... \n"

singularity run aertslab-pyscenic-0.10.0.sif \
  pyscenic aucell \
    FC_raw_counts.tsv \
    regulons_CLI.csv \
    -o auc_mtx.csv \
    --num_workers 20


echo "Done."
