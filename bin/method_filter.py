#!/usr/bin/env python3

"""
Filter AnnData

Usage:
    method_filter.py --file-out=<path> [options] <file>
    
Options:
    -h --help               Show this screen.
    -o --file-out=<path>    Output file (AnnData)
    --random                Whether to select 200 random cells per sample
    --subset                Whether to select subset of samples
    --cohortA               Whether to select cohort A
"""


def subset_by_qc_labels(adata, label_key, label_value):
    """
    Subset AnnData using quality labels
    
    Parameters
    ----------
    adata : AnnData
    label_key : str
        Column key for 'obs' used as label.
    label_value : str
        Value to select in 'label_key'.

    Returns
    -------
    AnnData

    """
    
    # Create index
    ind = adata.obs[label_key].isin([label_value])
    
    # Remove cells
    adata = adata[ind, :]
    
    return(adata)


def subset_random(adata, n_cells, batch_key=None):
    """
    Subset AnnData by random

    Parameters
    ----------
    adata : AnnData
    n_cells : int
        Number of cells to sample (per batch).
    batch_key : str
        Column key for 'obs' used as batch.

    Returns
    -------
    AnnData

    """
    import numpy as np
    import pandas as pd

    if (batch_key in adata.obs_keys()):
        print(f"Selecting {n_cells} cells from each {batch_key}")
        batches = adata.obs[batch_key]
    else:
        print(f"No column '{batch_key}' in obs.")
        print(f"Selecting {n_cells} cells total.")
        batches = pd.array(np.repeat(1, len(adata.obs)))

    # Sample cells from batches
    cells = []
    for i in batches.unique():
        pop = adata.obs.index[batches.isin([i])]
        sel = np.random.choice(pop, n_cells, replace=False)
        cells = np.append(cells, sel)
    
    # Subset based on cells
    adata = adata[cells, :]

    return(adata)


def subset_sample(adata):
    """
    Subset AnnData by sample
    """       
    # Select samples
    samples = [
        "C19-CB-0062",
        "C19-CB-0082",
        "C19-CB-0474",
        "C19-CB-0660",
        "C19-CB-0694",
        "X3"
        ]
    index = adata.obs.patient.isin(samples)
    sum(index)

    # Subset
    adata = adata[index, :]

    return(adata)

def subset_cohortA(adata):
    """
    Subset to cohort A
    """
    index = adata.obs.cohort.isin(["A"])
    sum(index)    
    
    adata = adata[index, :]
    
    return(adata)
    

def main():
    """The main script function"""
    from docopt import docopt
    import scanpy as sc
    import pandas as pd

    # Variables
    args = docopt(__doc__)
    
    in_file = args["<file>"]
    in_file_csv = in_file.replace(".h5ad", "_qc_colData.csv")
    in_random = args["--random"]
    in_subset = args["--subset"]
    in_cohortA = args["--cohortA"]
    
    out_file = args["--file-out"]
    
    print(f"Reading dataset from '{in_file}'")
    ds = sc.read(in_file)
    csv = pd.read_csv(in_file_csv, index_col = 0)
    ds.obs = csv
    
    print("Subsetting data")
    ds = subset_by_qc_labels(ds, "qc_linear", "good")
    
    if in_random:
        ds = subset_random(ds, 200, "sample")
        
    if in_subset:
        ds = subset_sample(ds)
        
    if in_cohortA:
        ds = subset_cohortA(ds)

    print("Dataset after subsetting:")
    print(ds)

    # Write data
    print(f"Writing output to '{out_file}'")
    ds.write_h5ad(out_file)

    print("Done.")

if __name__ == "__main__":
    main()