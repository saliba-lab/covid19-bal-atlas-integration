# Setup the combined BAL samples as a Anndata object

# Import packages --------------------------------------------------------------
import scanpy as sc

# Global variables -------------------------------------------------------------
fromfile = "data/combined/filtered_feature_bc_matrix.h5"
tofile = "data/combined/anndata.h5ad"

# Load data --------------------------------------------------------------------
ad = sc.read_10x_h5(fromfile)
ad.var_names_make_unique()

# Save data
ad.write_h5ad(tofile)
