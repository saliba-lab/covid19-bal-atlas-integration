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
    import tempfile
    import scanpy as sc
    
    # Variables
    url = "https://nubes.helmholtz-berlin.de/s/A7XQ3yZYGMnFR3H"
    url = os.path.join(url, "download")

    # Download
    print(f"Downloading matrix from '{url}'...")
    temp_dir = tempfile.TemporaryDirectory()
    temp_file = os.path.join(temp_dir.name, "matrix.h5")
    os.system("curl" + " " + url + " -o " + temp_file)

    # Object creation
    print("Creating AnnData object")
    adata = sc.read_10x_h5(temp_file)
    adata.var_names_make_unique()

    temp_dir.cleanup()

    print("Adding sample identifiers in 'obs'")
    v = []
    for i in adata.obs.index:
        v.append(int(i.split("-")[1]))
    adata.obs["sample_id"] = v
    adata.obs["sample_id"] = adata.obs["sample_id"].astype("string")

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
    out_file_viral = os.path.dirname(out_file) + "/SCoV2_" + os.path.basename(out_file)
    
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
    
    # Separate human & viral counts
    print("Moving viral counts in separate layer")
    viral_ids = [s for s in ds.var.gene_ids if "ENSSASG" in s]
    viral_ind = ds.var.gene_ids.isin(viral_ids)
    viral_genes = ds.var.index[viral_ind]
    vdata = sc.AnnData(ds.X[:, viral_ind])
    vdata.var.index = viral_genes
    vdata.var["gene_ids"] = viral_ids
    vdata.obs = ds.obs
    
    print("Removing viral counts from X")
    human_ids = [s for s in ds.var.gene_ids if "ENSG" in s]
    human_ind = ds.var.gene_ids.isin(human_ids)
    ds = ds[:, human_ind]
    
    print(f"Writing output to '{out_file}'")
    ds.write_h5ad(out_file)
    vdata.write_h5ad(out_file_viral)


if __name__ == "__main__":
    main()
