#!/usr/bin/env python3

"""
Download the combined BAL dataset

Usage:
    dataset-combined.py [options]
    
Options:
    -h --help           Show this screen.
    -o --overwrite      Overwrite existing file.
"""

import scanpy as sc
import pandas as pd

# Global parameters
file_input_adata = "data/combined.h5ad"
file_input_obs = "docs/combined_coldata.csv"
file_output_adata_random = "data/combined-random.h5ad"
file_output_adata = "data/combined-qc.h5ad"

# Read files
adata = sc.read(file_input_adata)
obs = pd.read_csv(file_input_obs)

# Add sample metadata (obs)
obs.index = obs["Unnamed: 0"].values
del obs["Unnamed: 0"]
obs["sample"] = obs["sample"].astype("category")
adata.obs = obs

# Write dataset with random cells
index = adata.obs.random
ad_random = adata[index, :]
sc.pp.filter_genes(ad_random, min_cells=3)
ad_random.write(file_output_adata_random)

# Write dataset with high quality cells
key = "qc.loess"
index = adata.obs[key].isin(["good"])
adata = adata[index, :]
sc.pp.filter_genes(adata, min_cells=3)
adata.write(file_output_adata)
