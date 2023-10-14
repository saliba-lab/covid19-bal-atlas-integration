#!/usr/bin/env python3

"""
Download the combined BAL dataset

Usage:
    dataset-combined.py [options] <file>
    
Options:
    -h --help                   Show this screen.
    -o --overwrite              Overwrite existing file.
"""


def get_combined():
    """
    Get the combined BAL dataset

    Returns
    -------
    AnnData containing the combined BAL dataset.

    """
    import os.path
    import scanpy as sc
    import pandas as pd
    
    # Variables ---------------------------------------------------------------
    h5file = "data/BCB/filtered_feature_bc_matrix.h5"
    url = "https://nubes.helmholtz-berlin.de/s/yPaeYMtAqmrbJYD"
    url = os.path.join(url, "download")
    
    aggr_csv = "data/BCB/aggregation.csv"
    aggr_url = "https://nubes.helmholtz-berlin.de/s/Z3KMkrYnirwfCqb"
    aggr_url = os.path.join(aggr_url, "download")

    # Download data -----------------------------------------------------------
    
    # Counts
    if os.path.exists(h5file):
        print(f"File '{h5file}' already exists.")
    else:
        print(f"Downloading matrix from '{url}'...")
        os.system("curl" + " " + url + " -o " + h5file)
        
    # Library info
    if os.path.exists(aggr_csv):
        print(f"File '{aggr_csv}' already exists.")
    else:
        print(f"Downloading aggregation.csv from '{aggr_url}'...")
        os.system("curl" + " " + aggr_url + " -o " + aggr_csv)

    # Create AnnData ----------------------------------------------------------
    print("Creating AnnData object")
    adata = sc.read_10x_h5(h5file, gex_only=False)
    adata.var_names_make_unique()

    print("Adding sample identifiers in 'obs'")
    aggr = pd.read_csv(aggr_csv)
    aggr.index = aggr.index+1
    v = []
    for i in adata.obs.index:
        v.append(int(i.split("-")[1]))
    
    adata.obs["sample_id"] = aggr["sample_id"][v].values
    
    return adata


def main():
    """The main script function"""
    from docopt import docopt
    import os.path
    import os
    import scanpy as sc

    args = docopt(__doc__)
    
    flag_overwrite = args["--overwrite"]
    out_file = args["<file>"]
    dname=os.path.dirname(out_file)
    out_file_viral = dname + "/alt_" + os.path.basename(out_file)
    
    # Checks
    fpath = os.path.dirname(out_file)
    if not(os.path.exists(fpath)):
        print(f"Creating directory '{fpath}'")
        os.makedirs(fpath)
        
    if (os.path.exists(out_file) and not flag_overwrite):
        print(f"'{out_file}' already exists. Exiting")
        exit(2)
        
    print("Read dataset:")
    ds = get_combined()
    print(ds)
    
    # Separate assays
    print("Moving other counts in separate AnnData...")
    viral_ids = [s for s in ds.var.gene_ids if "ENSG" not in s]
    viral_ind = ds.var.gene_ids.isin(viral_ids)
    viral_genes = ds.var.index[viral_ind]
    vdata = sc.AnnData(ds.X[:, viral_ind])
    vdata.var.index = viral_genes
    vdata.var["gene_ids"] = viral_ids
    vdata.obs = ds.obs
    
    print("Retaining only human gene counts in X...")
    human_ids = [s for s in ds.var.gene_ids if "ENSG" in s]
    human_ind = ds.var.gene_ids.isin(human_ids)
    ds = ds[:, human_ind]
    
    print(f"Writing output to '{out_file}...'")
    ds.write_h5ad(out_file)
    vdata.write_h5ad(out_file_viral)


if __name__ == "__main__":
    main()
