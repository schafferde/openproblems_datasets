import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite
import scanpy as sc


adata = sc.read_h5ad("Adult_Peaks.h5ad") #Pre-gzipped
adata.obs.index.name = "cell_label"

#Extract metadata as is
df = adata.obs[["cell_type"]].copy()
df["batch"] = adata.obs["batch_key"]
df.to_csv("Kanemaru_metadata.tsv", sep="\t")
"""
peaks = adata.var.index.values
with open("Kanemaru_peaks.tsv", "w") as f:
    for peak in peaks:
        print(peak.replace("_", "-", 1), file=f) #Convert chr:start_end to chr:start-end

#Extract ATAC data
atac_matrix = adata.X
sparse_matrix = sp.csr_matrix(atac_matrix.T)
print("Matrix shape:", sparse_matrix.shape)
mmwrite("Kanemaru_counts.mtx", sparse_matrix)
print("Wrote data")
"""
