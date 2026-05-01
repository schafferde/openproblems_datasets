from scipy.io import mmread
from scipy.sparse import vstack, csr_matrix
import pandas as pd
import anndata as ad
import numpy as np

### DROPSEQ

droplet_obs_full = pd.read_csv("annotations_droplet.csv")

#Problem: some cell ontology classses are na
#free_annotation values of missing cell_ontology_class
#chondrocyte-like                                                                         383
#unknown                                                                                  244
#club cells, neuroendocrine cells, alveolar epithelial type 1 cells, and unknown cells     45
#Plan to add back chondrocyte-like
droplet_obs_full.loc[droplet_obs_full.free_annotation == 'chondrocyte-like', 'cell_ontology_class'] = droplet_obs_full.free_annotation
droplet_obs = droplet_obs_full[['cell', 'cell_ontology_class', 'channel', 'mouse.id', 'tissue']]
droplet_obs = droplet_obs[~droplet_obs.cell_ontology_class.isna()]
#Left with 55348 cells, filtering out 408

#Some samples contain the exact # of cells/barcodes. Some contain a few more, presumably ones that were filtered at a later point
#Some contain 737280 barcodes, which is the max # for 10x v2
grouped_obs = [(name, group) for name, group in droplet_obs.groupby(['tissue', 'channel'])]
matrices = []
genes = None
for (tissue, channel), sample_obs in grouped_obs:
    dir = f"droplet/{tissue}-{channel}"
    matrix = mmread(f"{dir}/matrix.mtx").transpose().tocsr()
    barcodes_filtered = [x.rsplit("_", 1)[1] for x in sample_obs.cell.tolist()]
    with open(f"{dir}/barcodes.tsv") as f:
        barcodes_matrix = [x.strip().split("-")[0] for x in f.readlines()]
    matrix_indices = [barcodes_matrix.index(b) for b in barcodes_filtered]
    matrix_filtered = matrix[matrix_indices]
    if genes is None:
        with open(f"{dir}/genes.tsv") as f:
            genes = [x.split()[0] for x in f.readlines()]
    matrices.append(matrix_filtered)

reordered_obs =  pd.concat([df for _, df in grouped_obs], ignore_index=True)
merged_matrix = vstack(matrices)
droplet_adata = ad.AnnData(X=merged_matrix, obs=reordered_obs, var={'gene_name': genes})
droplet_adata.var.index = genes
droplet_adata.write_h5ad("Muris_dropseq_dataset.h5ad")


###FACS
facs_obs_full = pd.read_csv("annotations_facs.csv")
facs_obs = facs_obs_full[['cell', 'cell_ontology_class', 'mouse.id', 'tissue']]
facs_obs = facs_obs[~facs_obs.cell_ontology_class.isna()]
#Left with 44779 cells, after filtering 170

grouped_obs = [(name, group) for name, group in facs_obs.groupby('tissue')]
matrices = []
genes = None
for tissue, sample_obs in grouped_obs:
    print(tissue)
    file = f"FACS/{tissue}-counts.csv"
    with open(file, 'r') as f:
        matrix_cells = f.readline().strip().split(',')[1:]
        # We load everything as strings first to handle the row names in the first column
        raw_data = np.genfromtxt(f, delimiter=',', dtype=str, encoding='utf-8')
    matrix_genes = [x.replace('"', '') for x in raw_data[:, 0].tolist()]
    data_array = raw_data[:, 1:].astype(np.int64)
    matrix = csr_matrix(data_array.transpose()) #Now, cells are rows
    del data_array
    filtered_cells = sample_obs.cell.tolist()
    matrix_cells = [x.replace('"', '') for x in matrix_cells]
    matrix_indices = [matrix_cells.index(c) for c in filtered_cells]
    matrix_filtered = matrix[matrix_indices]
    if genes is None:
        genes = matrix_genes
    else:
        if genes != matrix_genes:
            print("mismatched genes for", tissue)
    matrices.append(matrix_filtered)

reordered_obs =  pd.concat([df for _, df in grouped_obs], ignore_index=True)
merged_matrix = vstack(matrices)
facs_adata = ad.AnnData(X=merged_matrix, obs=reordered_obs, var={'gene_name': genes})
facs_adata.var.index = genes
facs_adata.write_h5ad("Muris_facs_dataset.h5ad")

import scanpy as sc
droplet_adata = sc.read_h5ad("Muris_dropseq_dataset.h5ad")
facs_adata = sc.read_h5ad("Muris_facs_dataset.h5ad")

#Observation: gene sets are the same, but FACS's is sorted and Droplet's is not
droplet_adata = droplet_adata[:, sorted(droplet_adata.var_names)].copy()
assert droplet_adata.var_names.tolist() == facs_adata.var_names.tolist()
facs_adata.obs["channel"] = "-"
#Need to harmonize mouse names - using style from Droplet
def change_mouse(id):
    if "/" in id:
        id1, id2 = id.split("/")
        second_mouse = "/" + id2.split("_")[1]
        id = id1
    else:
        second_mouse = ""
    tokens = id.split("_")
    new_id = tokens[0] + "-" + tokens[2] + "-" + tokens[1]
    new_id += second_mouse
    return new_id
#Somewhere, the gene name field is lost, but it is retained in the index.
facs_adata.obs["mouse.id"] = [change_mouse(m) for m in facs_adata.obs["mouse.id"].tolist()]
droplet_adata.obs_names = droplet_adata.obs.cell
facs_adata.obs_names = facs_adata.obs.cell
merged_adata = ad.concat([droplet_adata, facs_adata], label="assay", keys=["Droplet", "FACS"])
print(merged_adata.shape)
merged_adata.obs["batch"] = merged_adata.obs.assay.str.cat(merged_adata.obs["mouse.id"], sep="+")
merged_adata.write_h5ad("Muris_merged_dataset.h5ad")
#We have 100027 cells

#There are a couple of small batches (n=51, n=61). I decided to trim these
#This leaves 99915 cells
batch_counts = merged_adata.obs['batch'].value_counts()
drop_batches = batch_counts[batch_counts < 100].index
merged_adata_trim = merged_adata[~merged_adata.obs['batch'].isin(drop_batches)].copy()
print(merged_adata_trim.shape)
merged_adata_trim.write_h5ad("Muris_merged_dataset_trim.h5ad")


"""
#Could also consider trimming small cell types, e.g. 13 cell types with counts < 50 (22-49), totalling 478 cells
#This would 99549 cells (not considering trimming batches)
cell_type_counts = merged_adata.obs['cell_ontology_class'].value_counts()
drop_cell_types = cell_type_counts[cell_type_counts < 50].index
merged_adata_trim = merged_adata[~merged_adata.obs['cell_ontology_class'].isin(drop_cell_types)].copy()
print(merged_adata_trim.shape)
merged_adata_trim.write_h5ad("Muris_merged_dataset_trim.h5ad")
"""









