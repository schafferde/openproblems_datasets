import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite
#import gzip
import scanpy as sc


adata = sc.read_h5ad("GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad") #Pre-gzipped
adata.obs.index.name = "cell_label"

#Extract metadata as is
df = adata.obs[["batch", "cell_type"]]
df.to_csv("Burkhardt_metadata.tsv", sep="\t")
peaks = adata.var.index[adata.var.feature_types=="ATAC"].values
with open("Burkhardt_peaks.tsv", "w") as f:
    for peak in peaks:
        print(peak.replace("-", ":", 1), file=f) #Convert chr-start-end to chr:start-end

#Extract ATAC data
atac_matrix = adata.layers["counts"][:,adata.var.feature_types=="ATAC"]
sparse_matrix = sp.csr_matrix(atac_matrix.T)
print("Matrix shape:", sparse_matrix.shape)
#Write out the standard MatrixMarket file
#with gzip.open("Burkhardt_counts.mtx.gz", "wb") as f:
#    mmwrite(f, sparse_matrix)
mmwrite("Burkhardt_counts.mtx", sparse_matrix)
print("Wrote data")

#Preserve RNA as small h5ad
adata = adata[:,adata.var.feature_types=="GEX"].copy()
adata.X = adata.layers["counts"]
del adata.layers["counts"]
adata.var.index.name = "gene_name"
adata.write_h5ad("Burkhardt_RNA_counts.h5ad")

