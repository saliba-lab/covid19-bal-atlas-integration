#!/usr/bin/env python3

"""
Run scVI

Usage:
    method-scVI.py [options]
    
Options:
    -h --help           Show this screen.
    -o --overwrite      Overwrite existing file.
"""

import scanpy as sc

# Global parameters
file_input_adata = "data/combined-random.h5ad"

# Explore
sc.pl.embedding(adata, basis = "X_mde", color = "libsize")
