#!/usr/bin/env python3

"""
Compute integration metrics using scIB

Usage:
    method_integration-metrics.py [options] <file>
    
Options:
    -h --help       Show this screen.
    -l --labels     CSV with cell labels
"""

def calc_metrics_labelfree(adata):
    """Calculate integration metrics using scIB"""
    import scanpy as sc
    import scib
    import pandas as pd
    
    print("Create control object")
    key = adata.uns["control"].split("_")
    ctrl = sc.AnnData(obs = adata.obs, var=adata.var, X=adata.layers[key[0]])
    ctrl = ctrl[:, ctrl.var[key[1]]]
    
    print("Create index of embeddings to compute metrics for...")
    index = adata.obsm_keys()
    cols = ["batch_PCR", "iLISI"]
    df = pd.DataFrame(index=cols)
    
    for i in index:
        print(f"Computing metrics for {i}")
        v=[]
        test = sc.AnnData(obs=adata.obs, var=adata.var)
        test.obsm["X_emb"] = adata.obsm[i]
        
        print("batch PCR")
        v.append(scib.metrics.pcr_comparison(ctrl,test,covariate="batch",
                                             embed="X_emb",
                                             n_comps=int(key[3][1:3])
                                             ))
        
        print("iLISI")
        sc.pp.neighbors(test, use_rep="X_emb")
        v.append(scib.metrics.ilisi_graph(test, "batch"))
        df[i]=v
    
    return df.T


def calc_metrics_label(adata, label):
    """Calculate integration metrics using scIB"""
    import scanpy as sc
    import scib
    import pandas as pd
    
    print("Create control object")
    key = adata.uns["control"].split("_")
    ctrl = sc.AnnData(obs = adata.obs, var=adata.var, X=adata.layers[key[0]])
    ctrl = ctrl[:, ctrl.var[key[1]]]
    
    print("Create index of embeddings to compute metrics for...")
    index = adata.obsm_keys()
    cols = ["ASW", "isolated_labels", "NMI", "ARI", "cLISI", 
            "graph_connectivity", "batch_ASW", "kBET"]
    df = pd.DataFrame(index=cols)
    
    for i in index:
        print(f"Computing metrics for {i}")
        v=[]
        test = sc.AnnData(obs=adata.obs, var=adata.var)
        test.obsm["X_emb"] = adata.obsm[i]
        
        print("ASW")
        v.append(scib.metrics.silhouette(test, label, embed="X_emb"))

        print("isolated_labels")
        v.append(scib.metrics.isolated_labels(test, label, "batch", "X_emb"))

        print("NMI")
        v.append(scib.metrics.nmi(test,label,"batch"))
        sc.pp.neighbors(test, use_rep="X_emb")
        nmi={}
        for j in range(1, 21):
            sc.tl.louvain(test, resolution=j)
            nmi[str(j)] = scib.metrics.nmi(test, label, "louvain")
        sc.tl.louvain(test, resolution=int(max(nmi))/10)
        
        max(nmi.values())
        np.where(list(nmi.values()) == max(nmi.values()))
        
        # Add metric output
        df[i] = v
    
    return df.T
    

    # ASW batch
    print("Calculating ASW for 'batch'")
    v = []
    for i in df.index:
        v.append(scib.metrics.silhouette_batch(adata, "batch", 
                                               group_key=label, embed=i))
    df["ASW_batch"]=v
    
    # ASW label
    print("Calculating silhouette width...")
    v = []
    for i in df.index:
        v.append(scib.metrics.silhouette(adata, group_key=label, embed=i))
    df["ASW_label"]=v
    
    # Isolated label ASW/F1
    #print("Calculating isolated label score...")
    #v = []
    #for i in df.index:
    #    v.append(scib.metrics.isolated_labels(adata, "celltype", "batch", i))
    #df["isolated_labels"]=v
    
    return df


def main():
    """The main script function"""
    from docopt import docopt
    import scanpy as sc
    import pandas as pd

    # Variables
    args = docopt(__doc__)
    
    in_file = args["<file>"]
    in_file_labels = args["--labels"]
    out_file = in_file
    
    print(f"Reading dataset from '{in_file}'")
    ds = sc.read(in_file)
    if not "metrics" in ds.uns_keys():
        ds.uns["metrics"] = {}
    
    print("Select control embedding")
    ds.uns["control"] = "cp10k_seurat5000_PCA_n30"

    print("Calculating label-free metrics")
    ds.uns["metrics"]["labelfree"] = calc_metrics_labelfree(ds)
    
    labels = pd.read_csv(in_file_labels, index_col = 0)
    for i in labels.columns:
        ds.obs[i] = labels[i]
        print(f"Calculating integration metrics with label '{i}'")
        ds.uns["metrics"][i] = calc_metrics_label(ds, label = i)

    # Write data
    print(f"Writing output to '{out_file}'")
    ds.write_h5ad(out_file)
    

if __name__ == "__main__":
    main()