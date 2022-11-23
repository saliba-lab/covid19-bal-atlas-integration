#!/usr/bin/env python3

"""
Feature selection

Usage:
    method-hvfeatures.py [options] <file>
    
Options:
    -h --help               Show this screen.
"""

def select_features(adata):
    """
    Run the scanpy core analysis

    Parameters
    ----------
    adata : AnnData

    Returns
    -------
    AnnData.

    """
    import scanpy as sc
    
    adata.obs["batch"] = adata.obs["patient"].astype("str")
    adata.X = adata.layers["cp10k"].copy()
    adata.uns["log1p"] = {'base': None}
    
    # Add keys to uns
    adata.uns["hvg_keys"] = {"keys" : []}
    
    # Seurat
    for i in [2000,5000]:
        sc.pp.highly_variable_genes(adata, batch_key="batch", n_top_genes=i)
        key = "seurat" + str(i)
        adata.var[key] = adata.var.highly_variable
        adata.uns["hvg_keys"]["keys"].append(key)
    
    return(adata)


def main():
    """The main script function"""
    from docopt import docopt
    import scanpy as sc

    # Variables
    args = docopt(__doc__)
    
    in_file = args["<file>"]
    out_file = in_file
    
    print(f"Reading dataset from '{in_file}'")
    ds = sc.read(in_file)
    
    print("Selecting variable features...")
    ds = select_features(ds)

    # Write data
    print(f"Writing output to '{out_file}'")
    ds.write_h5ad(out_file)

    print("Done.")

if __name__ == "__main__":
    main()