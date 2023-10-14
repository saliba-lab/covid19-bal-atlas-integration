#!/usr/bin/env python3

"""
Integration by scVI

Usage:
    method_integration-scVI.py [options] <file>

Options:
    -h --help           Show this screen.
    -g --gpu		Use GPU for training of the VAE
"""

def run_scVI(adata, n_latent, n_layers, n_hidden, hvg_method, use_gpu):
    """
    Run scVI

    Parameters
    ----------
    input : AnnData

    Returns
    -------
    AnnData.

    """
    import scvi

    # Subset to HVGs
    adata = adata[:, adata.var[hvg_method]]

    # Batch correction with scVI
    adata = adata.copy()
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
    vae = scvi.model.SCVI(adata, n_latent=n_latent, n_layers=n_layers, 
                          n_hidden=n_hidden, gene_likelihood="nb")
    vae.train(use_gpu=use_gpu)
    adata.obsm["X_emb"] = vae.get_latent_representation()

    return(adata.obsm["X_emb"])


def main():
    """The main script function"""
    from docopt import docopt
    import scanpy as sc

    # Variables
    args = docopt(__doc__)

    in_file = args["<file>"]
    gpu_flag = args["--gpu"]

    out_file = in_file

    print(f"Reading dataset from '{in_file}'")
    ds = sc.read(in_file)

    print("Running scVI...")
    layer = "counts"
    ds.X = ds.layers[layer]
    n_latent = [10,30,50]
    n_layers = [1]
    n_hidden = [128]
    hvg_method = ds.uns["hvg_keys"]["keys"]
    for n in n_latent:
        for l in n_layers:
            for h in n_hidden:
                for hvg in hvg_method:
                    print(f"{n} dimensions from {l} layer(s) with {h} nodes from {hvg} features")
                    arch = "n" + str(n) + "l" + str(l) + "h" + str(h)
                    key = layer + "_" + hvg + "_" + "scVI" + "_" + arch
                    ds.obsm[key] = run_scVI(ds, n, l, h, hvg, gpu_flag)

    # Write data
    print(f"Writing output to '{out_file}'")
    ds.write_h5ad(out_file)


if __name__ == "__main__":
    main()
