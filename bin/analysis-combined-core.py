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

# Global parameters
file = "data/combined.h5ad"
cells = "Something R outputs"
genes = "Something R outputs"

# Read files
adata = sc.read(file)

# Subset adata
adata = adata[adata.obs_names.isin(cells),:]
adata = adata[:,adata.var_names.isin(genes)]

# Continue with Tutorial pipeline...