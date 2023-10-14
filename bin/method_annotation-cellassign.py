#!/usr/bin/env python3

"""
Feature selection

Usage:
    method-annotation-cellassign.py [options] <file>
    
Options:
    -h --help               Show this screen.
"""

def run_cellassign(adata, mrk):
    """
    Run the scanpy core analysis

    Parameters
    ----------
    adata : AnnData

    Returns
    -------
    CellAssign model.

    """
    import scvi
    
    # Filter data
    bdata = adata[:, mrk.index].copy()
    # del adata
    
    # Use raw counts
    bdata.X = bdata.layers["counts"]
    
    # Set up object
    scvi.external.CellAssign.setup_anndata(bdata, 
                                           size_factor_key="size_factor")
    
    # Train model
    model = scvi.external.CellAssign(bdata, mrk)
    model.train()
    
    return(model)


def main():
    """The main script function"""
    from docopt import docopt
    import scanpy as sc
    import pandas as pd
    import numpy as np

    # Variables
    args = docopt(__doc__)
    
    in_file = args["<file>"]
    in_markers = args["<markers>"]
    out_file = in_file
    
    print(f"Reading dataset from '{in_file}'")
    ds = sc.read(in_file)

    print(f"Reading markers from '{in_markers}'")
    markers = pd.read_excel(in_markers, "level_3")
    mrk = pd.DataFrame()
    mrk["type"] = markers["type"]
    mrk.index = markers["gene"].values
    mrk = pd.get_dummies(mrk)
    
    if "size_factor" in ds.obs:
        print("Size factor already present.")
    else:
        print("Adding size factor...")
        lib_size = ds.obs["libsize"].values
        ds.obs["size_factor"]= lib_size / np.mean(lib_size)
    
    print("Running CellAssign...")
    label = run_cellassign(ds, mrk)

    # Write data
    print(f"Writing output to '{out_file}'")
    ds.write_h5ad(out_file)

    print("Done.")

if __name__ == "__main__":
    main()