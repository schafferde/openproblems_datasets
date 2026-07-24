from typing import Any, Callable, Dict, Tuple
import scanpy as sc
import scipy
import os
import time

## VIASH START
par = {
    "file_name": "pancreas.h5ad",
    "dataset_id": "pancreas",
    "obs_cell_type": "cell_type",
    "obs_batch": "tech",
    "obs_tissue": "tissue",
    "layer_counts": "counts",
    "output": "test_data.h5ad",
}
meta = {
    "resources_dir": "src/datasets/loaders/openproblems_v1/"
}
## VIASH END

# fetch dataset
file_name = par["file_name"]
print("Looking for file:", file_name, flush=True)
while not os.path.isfile(file_name):
    time.sleep(30)
time.sleep(15)

adata = sc.read_h5ad(file_name) #TODO: did the old dataset loading functions do anything?

#This was an alternative for setting hardcoded uns metadata values. See below for argument-based
"""
# override values one by one because adata.uns and
# metadata are two different classes.
for key, value in dataset_fun.metadata.items():
    print(f"Setting .uns['{key}']", flush=True)
    adata.uns[key] = value
"""
print("Setting .obs['cell_type']", flush=True)
if par["obs_cell_type"] and par["obs_cell_type"] != "cell_type":
    if par["obs_cell_type"] in adata.obs:
        adata.obs["cell_type"] = adata.obs[par["obs_cell_type"]]
    else:
        print(f"Warning: key '{par['obs_cell_type']}' could not be found in adata.obs.", flush=True)

print("Setting .obs['batch']", flush=True)
if par["obs_batch"] and par["obs_batch"] != "batch":
    if par["obs_batch"] in adata.obs:
        adata.obs["batch"] = adata.obs[par["obs_batch"]]
    else:
        print(f"Warning: key '{par['obs_batch']}' could not be found in adata.obs.", flush=True)

print("Setting .obs['tissue']", flush=True)
if par["obs_tissue"] and par["obs_tissue"] != "tissue":
    if par["obs_tissue"] in adata.obs:
        adata.obs["tissue"] = adata.obs[par["obs_tissue"]]
    else:
        print(f"Warning: key '{par['obs_tissue']}' could not be found in adata.obs.", flush=True)

if par["layer_counts"] and par["layer_counts"] in adata.layers:
    print(f"Temporarily moving .layers['{par['layer_counts']}'] to .X", flush=True)
    adata.X = adata.layers[par["layer_counts"]]
    del adata.layers[par["layer_counts"]]

#Empty X can arrise from datasets that contain only pre-processed obsm embedding fields
is_empty = (adata.X is None) or (0 in adata.X.shape)

if is_empty:
    print("This AnnData does not contain any expression values")
    adata.layers["counts"] = adata.X
    print(f"However, it contains the following obsm fields:{adata.obsm_keys()}")
else:
    if par["sparse"] and not scipy.sparse.issparse(adata.X):
        print("Make counts sparse", flush=True)
        adata.X = scipy.sparse.csr_matrix(adata.X)

    print("Removing empty genes", flush=True)
    sc.pp.filter_genes(adata, min_cells=1)

    print("Removing empty cells", flush=True)
    sc.pp.filter_cells(adata, min_counts=2)

    print("Moving .X to .layers['counts']", flush=True)
    adata.layers["counts"] = adata.X
    del adata.X

print("Add metadata to uns", flush=True)
metadata_fields = [
    "dataset_id", "dataset_name", "dataset_url", "dataset_reference",
    "dataset_summary", "dataset_description", "dataset_organism"
]
uns_metadata = {
    id: par[id]
    for id in metadata_fields
    if id in par
}
adata.uns.update(uns_metadata)

print("Setting .var['feature_name']", flush=True)
if not is_empty:
    if par["var_feature_name"] == "index":
        adata.var["feature_name"] = adata.var.index
    else:
        if par["var_feature_name"] in adata.var:
            adata.var["feature_name"] = adata.var[par["var_feature_name"]]
            del adata.var[par["var_feature_name"]]
        else:
            print(f"Warning: key '{par['var_feature_name']}' could not be found in adata.var.", flush=True)

    print("Setting .var['feature_id']", flush=True)

    if par["var_feature_id"] == "index":
        adata.var["feature_id"] = adata.var.index
    else:
        if par["var_feature_id"] in adata.var:
            adata.var["feature_id"] = adata.var[par["var_feature_id"]]
            del adata.var[par["var_feature_id"]]
        else:
            adata.var["feature_id"] = ""
            print(f"Warning: key '{par['var_feature_id']}' could not be found in adata.var.", flush=True)
else:
    adata.var["feature_id"] = None
    adata.var["feature_name"] = None

print("Writing adata to file", flush=True)
adata.write_h5ad(par["output"], compression=par["output_compression"])
