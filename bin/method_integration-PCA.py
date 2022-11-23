#!/usr/bin/env python3

"""
Integration by PCA

Usage:
    method_integration-PCA.py [options] <file>
    
Options:
    -h --help           Show this screen.
"""

def run_PCA(adata, dims, hvg_method):
    """
    Run PCA

    Parameters
    ----------
    adata : AnnData
    dims : Number of dimensions

    Returns
    -------
    AnnData.

    """
    import scanpy as sc

    # Subset to HVGs
    adata = adata[:, adata.var[hvg_method]]

    # Run PCA
    sc.tl.pca(adata, n_comps = dims)

    return(adata.obsm["X_pca"])


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

    print("Running PCA...")
    layer = ["cp10k"]
    dims = [10,30,50]
    hvg_method = ds.uns["hvg_keys"]["keys"]
    for l in layer:
        for n in dims:
            for hvg in hvg_method:
                print(f"{n} dimensions from {hvg} genes ({l})")
                ds.X = ds.layers[l]
                arch = "n" + str(n)
                key = "cp10k" + "_" + hvg + "_" + "PCA" + "_" + arch
                ds.obsm[key] = run_PCA(ds, n, hvg)

    # Write data
    print(f"Writing output to '{out_file}'")
    ds.write_h5ad(out_file)


if __name__ == "__main__":
    main()