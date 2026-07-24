
import scanpy as sc
import numpy as np
### VIASH START
par = {
  'input': 'resources_test/common/pancreas/dataset.h5ad',
  'input_layer': 'log_cp10k',
  'output': 'dataset.h5ad',
  'obsm_embedding': 'X_pca',
  'varm_loadings': 'pca_loadings',
  'uns_variance': 'pca_variance',
  'num_components': 25
}
### VIASH END

print(">> Load data", flush=True)
adata = sc.read(par['input'])

print(">> Look for layer", flush=True)
is_empty = (par['input_layer'] and par['input_layer'] not in adata.layers) or (0 in adata.shape)

if not is_empty:
    layer = adata.X if not par['input_layer'] else adata.layers[par['input_layer']]

    print(">> Run PCA", flush=True)
    X_pca, loadings, variance, variance_ratio = sc.tl.pca(
        layer, 
        n_comps=par["num_components"], 
        return_info=True
    )

    print(">> Storing output", flush=True)
    adata.obsm[par["obsm_embedding"]] = X_pca
    adata.varm[par["varm_loadings"]] = loadings.T
    adata.uns[par["uns_variance"]] = {
        "variance": variance, 
        "variance_ratio": variance_ratio
    }
else:
    nc = 1 if par["num_components"] is None else par["num_components"]
    adata.obsm[par["obsm_embedding"]] = np.zeros((adata.shape[0],nc))

print(">> Writing data", flush=True)
adata.write_h5ad(par['output'], compression=par["output_compression"])
