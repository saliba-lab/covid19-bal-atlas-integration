#!/usr/bin/env python3

"""
Normalization

Usage:
    method_normalization.py [options] <file>
    
Options:
    -h --help               Show this screen.
"""

def create_layers_from_X(adata, layers):
    """
    Create layers

    Parameters
    ----------
    adata : AnnData

    Returns
    -------
    AnnData

    """
    for i in layers:
        print(f"Creating layer '{i}'")
        adata.layers[i] = adata.X.copy()
    
    return(adata)


def normalize_data(adata):
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
    
    # Store counts
    create_layers_from_X(adata, ["counts"])
    
    # CP10K
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    create_layers_from_X(adata, ["cp10k"])
    
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

    print("Running normalization...")
    ds = normalize_data(ds)

    # Write data
    print(f"Writing output to '{out_file}'")
    ds.write_h5ad(out_file)

    print("Done.")

if __name__ == "__main__":
    main()