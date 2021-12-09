# run this for scenic.py
exec(open('scenic_env.py').read())


# Docs: https://buildmedia.readthedocs.org/media/pdf/pyscenic/stable/pyscenic.pdf

import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.distributed import Client
from dask.jobqueue import SLURMCluster # Needed to add this for parallelisation see here: http://jobqueue.dask.org/en/latest/
#import time # Maybe not needed? https://github.com/dask/distributed/issues/2941

from dask.diagnostics import ProgressBar 

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase 
from pyscenic.utils import modules_from_adjacencies, load_motifs 
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import seaborn as sns

DATABASES_GLOB = "SCENIC/encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.feather"

MOTIF_ANNOTATIONS_FNAME = "SCENIC/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
MM_TFS_FNAME = "SCENIC/hs_hgnc_tfs.txt"
SC_EXP_FNAME = "SCENIC/FC_raw_counts.tsv"
REGULONS_FNAME = "SCENIC/regulons.p"
MOTIFS_FNAME = "SCENIC/motifs.csv"


if __name__ == '__main__':
    local_cluster = SLURMCluster(n_workers=4, project="scw1641", cores=4, memory= "4G")
    custom_client = Client(local_cluster)
    local_cluster.adapt(minimum=1)
#    time.sleep(10) # Added for time function so may not be needed
#    [future] = Client.scatter([1.0], broadcast=True)

    # Load data
    ex_matrix = pd.read_csv(SC_EXP_FNAME, sep='\t', header=0, index_col=0).T 
    ex_matrix.shape

    tf_names = load_tf_names(MM_TFS_FNAME)

    db_fnames = glob.glob(DATABASES_GLOB) 
    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames] 
    dbs

    # Run GRNboost2
    adjacencies = grnboost2(ex_matrix, tf_names=tf_names, verbose=True, client_or_address=custom_client)
    adjacencies.to_csv('FC_network.tsv', sep='\t', header=False, index=False)

    # Create modules
    modules = list(modules_from_adjacencies(adjacencies, ex_matrix))

    # Calculate a list of enriched motifs and the corresponding target genes for all modules.
    with ProgressBar():
        df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

    # Create regulons from this table of enriched motifs.
    regulons = df2regulons(df)

    # Save the enriched motifs and the discovered regulons to disk.
    df.to_csv(MOTIFS_FNAME)
    with open(REGULONS_FNAME, "wb") as f:
        pickle.dump(regulons, f)

    # The clusters can be leveraged via the dask framework:
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME, client_or_address=SCHEDULER)

    # close the Client and LocalCluster after use
    custom_client.close()
    local_cluster.close()
