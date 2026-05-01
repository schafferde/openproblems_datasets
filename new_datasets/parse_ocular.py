from scipy.io import mmread
import pandas as pd
import anndata as ad
matrix_coo = mmread('raw_count.mtx')
matrix_csr = matrix_coo.transpose().tocsr()
obs_data = obs_data = pd.read_csv('SCP_meta.tsv', index_col=0, skiprows=[1], sep='\t')
with open('gene_names.tsv') as f:
    genes = [x.strip() for x in f.readlines()]
adata = ad.AnnData(X=matrix_csr, obs=obs_data, var={'gene_name': genes})
adata.var.index = genes
adata.write_h5ad("Ocular_dataset.h5ad")
