#!/usr/bin/env python3

"""
Download the combined BAL dataset

Usage:
    dataset-combined.py --file-out=<path> [options]
    
Options:
    -h --help                   Show this screen.
    -f --file-out=<path>        Output file
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
    import pandas as pd
    
    # Variables
    url = "https://nubes.helmholtz-berlin.de/s/A7XQ3yZYGMnFR3H"
    url = os.path.join(url, "download")

    # Download
    print(f"Downloading matrix from '{url}'")
    temp_dir = tempfile.TemporaryDirectory()
    temp_file = os.path.join(temp_dir.name, "matrix.h5")
    call = "curl" + " " + url + " -o " + temp_file
    os.system(call)

    # Object creation
    print("Create AnnData object")
    adata = sc.read_10x_h5(temp_file)
    adata.var_names_make_unique()

    temp_dir.cleanup()

    # Data wrangling
    print("Move viral counts in separate layer")
    viral_ids = [s for s in adata.var.gene_ids if "ENSSASG" in s]
    viral_ind = adata.var.gene_ids.isin(viral_ids)
    viral_genes = adata.var.index[viral_ind]
    vmat = adata.X[:, viral_ind].todense()
    vmat = pd.DataFrame(vmat)
    vmat.columns = viral_genes
    vmat.index = adata.obs.index
    adata.obsm["SCoV2_counts"] = vmat

    return adata


def main():
    """The main script function"""
    from docopt import docopt
    import os.path

    args = docopt(__doc__)
    
    flag_overwrite = args["--overwrite"]
    file_out = args["--file-out"]
    
    # Checks
    if (os.path.exists(file_out) and not flag_overwrite):
        print(f"'{file_out}' already exists. Exiting")
        exit(2)
        
    print("Read dataset:")
    output = get_combined()
    print(output)
    
    print(f"Writing output to '{file_out}'")
    output.write_h5ad(file_out)


if __name__ == "__main__":
    main()