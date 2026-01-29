"""
pySCENIC Analysis Pipeline for GitHub Submission.
Original analysis code for inferring Gene Regulatory Networks (GRN) and calculating regulon activity (AUCell).
"""

import os
import glob
import pickle
import pandas as pd
import numpy as np

# Distributed computing imports
from dask.diagnostics import ProgressBar
from dask.distributed import Client

# Arboreto (GRN inference) imports
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

# pySCENIC imports
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

import seaborn as sns

def run_analysis():
    # =========================================================================
    # 1. Configuration & Path Definitions
    # =========================================================================
    DATA_FOLDER = "tmpdir"
    RESOURCES_FOLDER = "resources"
    DATABASE_FOLDER = "tmpdir/cisTarget_databases"
    
    # Path patterns
    DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "mm9-*.mc9nr.genes_vs_motifs.rankings.feather")
    MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v9-nr.mgi-m0.001-o0.0.tbl")
    MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'mm_mgi_tfs.txt')
    
    # Output paths
    REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.p")
    MOTIFS_FNAME = os.path.join(DATA_FOLDER, "motifs.csv")
    OUTPUT_AUCELL_FNAME = "tmpdir/GSE231728_RAW/aucell_scores_t.csv"
    
    # Input Data
    ex_matrix_path = os.path.join(RESOURCES_FOLDER, "data_fibroblasts.csv")

    # Dask Configuration
    # Note: Define scheduler IP if connecting to an existing cluster.
    SCHEDULER = "123.122.8.24:8786" 

    # =========================================================================
    # 2. Initialize Resources & Data
    # =========================================================================
    
    # Load Ranking Databases
    db_fnames = glob.glob(DATABASES_GLOB)
    def name(fname):
        return os.path.splitext(os.path.basename(fname))[0]
    
    dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
    # Print loaded databases for verification
    print("Loaded databases:", dbs)

    # Load Expression Matrix
    # Assuming input CSV is Genes (rows) x Cells (cols). 
    # Transposing to Cells (rows) x Genes (cols) for pySCENIC.
    exprMat = pd.read_csv(ex_matrix_path, sep=',', header=0, index_col=0).T
    print("Expression Matrix Shape:", exprMat.shape)

    # Load TF list
    tf_names = load_tf_names(MM_TFS_FNAME)

    # =========================================================================
    # 3. Step I: Inference of Co-expression Modules (GRNBoost2)
    # =========================================================================
    
    # Initialize Dask Client
    # Note: Client() without arguments starts a local cluster. 
    # If SCHEDULER variable was intended, use Client(SCHEDULER).
    client = Client() 
    print("Dask Client Info:", client)

    # Run GRNBoost2
    print("Running GRNBoost2...")
    adjacencies = grnboost2(expression_data=exprMat, tf_names=tf_names, verbose=True)

    # Generate Modules
    modules = list(modules_from_adjacencies(adjacencies, exprMat, rho_mask_dropouts=True))
    print(f"Inferred {len(modules)} co-expression modules.")

    # =========================================================================
    # 4. Step II: Prune Modules for Targets (cisTarget)
    # =========================================================================
    print("Running cisTarget pruning...")
    with ProgressBar():
        df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

    # Create regulons from the enriched motifs table
    regulons = df2regulons(df)
    
    # =========================================================================
    # 5. Save & Checkpoint Results
    # =========================================================================
    # Save motifs to CSV
    df.to_csv(MOTIFS_FNAME)
    
    # Save regulons to pickle
    with open(REGULONS_FNAME, "wb") as f:
        pickle.dump(regulons, f)

    # Reloading specifically as per original workflow logic
    df_reloaded = load_motifs(MOTIFS_FNAME)
    with open(REGULONS_FNAME, "rb") as f:
        regulons_reloaded = pickle.load(f)
    
    print("Motifs dataframe head:")
    print(df_reloaded.head())

    # =========================================================================
    # 6. Step III: Cellular Regulon Enrichment (AUCell)
    # =========================================================================
    print("Running AUCell...")
    # Calculate AUC matrix
    auc_mtx = aucell(exprMat, regulons_reloaded, num_workers=4)
    print("AUCell Matrix preview:")
    print(auc_mtx.head())

    # Transpose: Rows = Regulons, Columns = Cells
    auc_mtx_t = auc_mtx.T
    print("Transposed AUCell Matrix preview:")
    print(auc_mtx_t.head())

    # Save final results
    # Ensure directory exists
    os.makedirs(os.path.dirname(OUTPUT_AUCELL_FNAME), exist_ok=True)
    
    auc_mtx_t.to_csv(
        OUTPUT_AUCELL_FNAME,
        index=True
    )
    print(f"Analysis complete. Results saved to {OUTPUT_AUCELL_FNAME}")

if __name__ == "__main__":
    run_analysis()