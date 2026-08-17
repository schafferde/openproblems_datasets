import scanpy as sc
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import diags

#From https://cells.ucsc.edu/?ds=retina-atac
#Via https://explore.data.humancellatlas.org/projects/9c20a245-f2c0-43ae-82c9-2232ec6b594f
adata = sc.read_h5ad("retina-atac.peaks.h5ad")

metadata = adata.obs[["cell_type", "tissue", "sample"]].copy()
metadata.index.name = "cell_label"
metadata["batch"] = adata.obs["donor_uuid"].values
metadata.to_csv("Liang_metadata.tsv", sep="\t")

peaks = adata.var.index.values
with open("Liang_peaks.tsv", "w") as f:
    f.write("\n".join(peaks))

#Reverse normalization, which was to 10000 counts followed by log1p
#It is not clear why the nCount_RNA field was used for counts/cell,
#but emperically treating it as the factor gives integer counts
matrix = adata.X
matrix.data = np.expm1(matrix.data)
factors = adata.obs["nCount_RNA"].to_numpy() / 10000
F = diags(factors)
raw_matrix = F @ matrix
raw_matrix.data = np.round(raw_matrix.data).astype(np.int32)
mmwrite("Liang_counts.mtx", raw_matrix.T)
