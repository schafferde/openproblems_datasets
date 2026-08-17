import pandas as pd
"""
import scipy.sparse as sp
from scipy.io import mmwrite
import gzip
#Counts matrix has row=peaks and columns=cells, with no labels
dense_df = pd.read_csv('GSE162170_atac_counts.tsv', sep='\t', header=None, dtype=int) #Change to gz
sparse_matrix = sp.csr_matrix(dense_df.values)
print("Matrix shape:", sparse_matrix.shape)
#Write out the standard MatrixMarket file
with gzip.open("Trevino_counts.mtx.gz", "wb") as f:
    mmwrite(f, sparse_matrix)
print("Wrote data")
"""
#For memory reasons, use to_mm.sh instead
import subprocess
result = subprocess.run(["./tsv_to_mm.sh", "GSE162170_atac_counts.tsv.gz", "Trevino_counts.mtx"], capture_output=True, text=True)
print(result.stdout)

#Now gather brief metadata
#Map cluster labels to cell types. These correspond to the labeles in Figure 1F, using text from Table S1E

df_label = pd.read_csv("trevino_cluster_map.tsv", sep="\t")

label_map = dict(zip(df_label["Cluster ID"], df_label["Name (long)"].apply(lambda x: x[:-2] if (x[-1].isdigit() and x[-2] == " ") else x)))


df = pd.read_csv('GSE162170_atac_cell_metadata.txt', sep="\t", header=0, index_col='Cell.ID')
df_new = df[["Sample.ID", "Iterative.LSI.Clusters"]].rename(columns={"Sample.ID": "batch", "Iterative.LSI.Clusters": "cluster"})
df_new["cell_type"] = df_new["cluster"].map(label_map)
print("Metadata shape:", df_new.shape)
df_new.to_csv("Trevino_metadata.tsv", sep="\t")


#And, reformat peaks

inFile = open("GSE162170_atac_consensus_peaks.bed")
outFile = open("Trevino_peaks.tsv", "w")
for line in inFile:
    tokens = line.split()
    outFile.write(f"{tokens[0]}:{tokens[1]}-{tokens[2]}\n")
inFile.close()
outFile.close()


