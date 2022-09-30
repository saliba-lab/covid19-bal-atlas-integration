#!/usr/bin/env python3

"""
Download the combined BAL dataset

Usage:
    dataset-combined.py [options]
    
Options:
    -h --help           Show this screen.
    -o --overwrite      Overwrite existing file.
"""

import os
import os.path
import tempfile
import scanpy as sc
from docopt import docopt

# Global settings
url = "https://nubes.helmholtz-berlin.de/s/A7XQ3yZYGMnFR3H"
url = os.path.join(url, "download")
out_file = "data/combined.h5ad"

# Check file presence / overwrite
args = docopt(__doc__)
if (os.path.exists(out_file) and not args["--overwrite"]):
    print(f"'{out_file}' already exists. Exiting")
    exit(2)

print(f"Downloading matrix from '{url}'")
temp_dir = tempfile.TemporaryDirectory()
temp_file = os.path.join(temp_dir.name, "matrix.h5")
call = "curl" + " " + url + " -o " + temp_file
os.system(call)

print("Create AnnData object")
adata = sc.read_10x_h5(temp_file)
adata.var_names_make_unique()

temp_dir.cleanup()

print(f"Saving AnnData to '{out_file}'")
adata.write_h5ad(out_file)