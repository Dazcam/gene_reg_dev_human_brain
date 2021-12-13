# run this for pyscenic.py to skip to line 59
# exec(open('pyscenic_env.py').read())
# Run in scratch folder


# Docs: https://buildmedia.readthedocs.org/media/pdf/pyscenic/stable/pyscenic.pdf

import os
import glob
import pickle
import pandas as pd
import numpy as np

from dask.distributed import Client
from dask_jobqueue import SLURMCluster # Needed to add this for parallelisation see here: http://jobqueue.dask.org/en/latest/
#import time # Maybe not needed? https://github.com/dask/distributed/issues/2941

from dask.diagnostics import ProgressBar 

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase 
from pyscenic.utils import modules_from_adjacencies, load_motifs 
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import seaborn as sns

SCENIC_DIR = "SCENIC/"
DATABASES_GLOB = os.path.join(SCENIC_DIR, "encode*.feather") 

MOTIF_ANNOTATIONS_FNAME = os.path.join(SCENIC_DIR, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
MM_TFS_FNAME = os.path.join(SCENIC_DIR, "hs_hgnc_tfs.txt")
SC_EXP_FNAME = os.path.join(SCENIC_DIR, "FC_raw_counts.tsv")
REGULONS_FNAME = os.path.join(SCENIC_DIR, "regulons.p")
MOTIFS_FNAME = os.path.join(SCENIC_DIR, "motifs.csv")


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

