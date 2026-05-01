from scipy.io import mmread
import pandas as pd
import anndata as ad

matrix_coo = mmread('GSE126954_gene_by_cell_count_matrix.txt')
matrix_csr = matrix_coo.transpose().tocsr()
obs_data = pd.read_csv("GSE126954_cell_annotation.csv", index_col=0)
var_data = pd.read_csv("GSE126954_gene_annotation.csv", index_col=0)
adata = ad.AnnData(X=matrix_csr, obs=obs_data, var=var_data)
print(adata.shape)
#Initial 89701
#QC filtering
adata = adata[adata.obs['passed_initial_QC_or_later_whitelisted'], :] 
#Now 86024
print(adata.shape)

#Now, apply cell type filtering
#Read in mapping of possible labels (fine or coarse) to new labels
#See notes file for details
df = pd.read_csv('celegans_cell_type_map.csv')
df = df[['Cell type from Table S2', 'Assigned new label']]
df_clean = df.dropna()
label_dict = dict(zip(df_clean['Cell type from Table S2'], df_clean['Assigned new label']))
new_label_set = set(df_clean['Assigned new label'])

filter = []
new_cell_types = []
skipped = []

#General idea is:
#For mature cells with top level cell type labels (all those listed in new_label_set), keep those
#We still look these up, but they map to themselves, except for HMC, where we merge 3 coarse labels into one
#For other cells, map fine cell labels to coarse labels as specified
#For late progenitor cells, keep if they have one of the permitted coarse cell type labels already
#(this is done implicitly since coarse labels are considered first)
#Or, if cutting off "*_parent" or "Parent_of_*" from one of the labels (usually the fine label)
# matches one of the allowed coarse labels, or mapped fine labels, we then keep the cell
#As a result, we exclude all early progenitors
#Disallowed labels include: G2_and_W_blasts, ABarpaaa_lineage (too early); AUA, T, XXX, (vague, see notes)

for coarse_label, fine_label in adata.obs[["cell.type", "cell.subtype"]].itertuples(index=False, name=None):
    skip = False
    if not pd.isna(coarse_label): 
        query = coarse_label.replace("Parent_of_", "").replace("_parent", "")
        if query in new_label_set: #Easiest case
            new_cell_types.append(query)
            filter.append(True)
            continue            
        elif coarse_label in label_dict: #Case where we choose to rename a coarse type - only two in practice
            new_cell_types.append(label_dict[coarse_label])
            filter.append(True)
            continue
        else:
            skip = True
    if not pd.isna(fine_label):
        query = fine_label.replace("parent_of_", "").replace("_parent", "")
        if query in label_dict:
            new_cell_types.append(label_dict[query])
            filter.append(True)
            continue
        skipped.append(fine_label)
    if skip:
        skipped.append(coarse_label)
    filter.append(False)

print(set(skipped))
adata = adata[filter].copy()
print(adata.shape)
adata.obs['new_cell_type'] = new_cell_types

adata.write_h5ad("Celegans_embryo_dataset_relabel.h5ad")
#Left with 62178 cells