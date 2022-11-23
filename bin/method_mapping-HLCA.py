#!/usr/bin/env python3

"""
Mapping by scArches

Usage:
    method_mapping-HLCA.py [options] <file>
    
Options:
    -h --help           Show this screen.
"""

def map_HLCA(query, ref):
    """
    Run scArches

    Parameters
    ----------
    input : AnnData

    Returns
    -------
    AnnData.

    """
    import tempfile
    import os.path
    import scvi
    import pynndescent
    import scarches as sca
    import scanpy as sc
    import numba
    import numpy as np
    
    # Prepare query
    query = query.copy()
    query.X = query.layers["counts"]
    query.obs["dataset"] = query.obs.patient
    query.obs["scanvi_label"] = "unlabeled"
    query.var["gene_name"] = query.var.index.tolist()
    query.var.index = query.var["gene_ids"].values
    
    # Get reference model
    url = "https://zenodo.org/record/6337966/files/HLCA_reference_model.zip"
    temp_dir = tempfile.TemporaryDirectory()
    temp_file = os.path.join(temp_dir.name, "HLCA_reference_model.zip")
    os.system("curl" + " " + url + " -o " + temp_file)
    os.system("unzip" + " " + temp_file + " -d " + temp_dir.name)
    ref_model = os.path.join(temp_dir.name, "HLCA_reference_model")
    scvi.model.SCANVI.convert_legacy_save(ref_model, ref_model, True, )
    
    # Build NN index
    X_train = ref.obsm["X_scanvi_emb"]
    ref_nn_index = pynndescent.NNDescent(X_train)
    ref_nn_index.prepare()    
    
    # Prepare scArches object
    scvi.model.SCANVI.prepare_query_anndata(query, ref_model)
    sca.models.SCANVI.load_query_data(query, ref_model)
    scvi.model.SCANVI.view_setup_args(ref_model)
    query_model = scvi.model.SCANVI.load_query_data(query, ref_model)

    # scArches settings
    surgery_epochs = 500
    train_kwargs_surgery = {
        "early_stopping": True,
        "early_stopping_monitor": "elbo_train",
        "early_stopping_patience": 10,
        "early_stopping_min_delta": 0.001,
        "plan_kwargs": {"weight_decay": 0.0},
    }
    
    # Train model
    query_model.train(max_epochs=surgery_epochs, **train_kwargs_surgery)
    
    # Cleanup
    temp_dir.cleanup()

    # Store latent space
    query_emb = sc.AnnData(query_model.get_latent_representation())
    query_emb.obs_names = query.obs_names
    
    ref_neighbors, ref_distances = ref_nn_index.query(query_emb.X)    
    
    # convert distances to affinities
    stds = np.std(ref_distances, axis=1)
    stds = (2.0 / stds) ** 2
    stds = stds.reshape(-1, 1)
    ref_distances_tilda = np.exp(-np.true_divide(ref_distances, stds))
    weights = ref_distances_tilda / np.sum(
        ref_distances_tilda, axis=1, keepdims=True
        )

    @numba.njit
    def weighted_prediction(weights, ref_cats):
        """Get highest weight category."""
        N = len(weights)
        predictions = np.zeros((N,), dtype=ref_cats.dtype)
        uncertainty = np.zeros((N,))
        for i in range(N):
            obs_weights = weights[i]
            obs_cats = ref_cats[i]
            best_prob = 0
            for c in np.unique(obs_cats):
                cand_prob = np.sum(obs_weights[obs_cats == c])
                if cand_prob > best_prob:
                    best_prob = cand_prob
                    predictions[i] = c
                    uncertainty[i] = max(1 - best_prob, 0)
                    
        return predictions, uncertainty
    
    # for each annotation level, get prediction and uncertainty
    label_keys = [f"ann_level_{i}" for i in range(1, 6)] + ["ann_finest_level"]
    for l in label_keys:
        ref_cats = ref.obs[l].cat.codes.to_numpy()[ref_neighbors]
        p, u = weighted_prediction(weights, ref_cats)
        p = np.asarray(ref.obs[l].cat.categories)[p]
        query_emb.obs[l + "_pred"], query_emb.obs[l + "_uncertainty"] = p, u
        
    return(query_emb)


def main():
    """The main script function"""
    from docopt import docopt
    import scanpy as sc

    # Variables
    args = docopt(__doc__)
    
    in_file = args["<file>"]
    
    in_ref = "data/HLCA/core.h5ad"
    out_file = in_file
    
    print("Reading datasets...")
    ds = sc.read(in_file)
    ref = sc.read(in_ref)

    print("Running scArches...")
    emb = map_HLCA(ds, ref)
    ds.obsm["hlca"] = emb.X
    for i in emb.obs.columns:
        print(f"Transferring {i}")
        ds.obs[i] = emb.obs[i]
    
    # Write data
    print(f"Writing output to '{out_file}'")
    ds.write_h5ad(out_file)


if __name__ == "__main__":
    main()