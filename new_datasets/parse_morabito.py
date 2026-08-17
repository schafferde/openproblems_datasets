import h5py
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.io import mmwrite

# 1. Load the metadata and retain order
metadata_path = "GSE174367_snATAC-seq_cell_meta.csv"
meta_df = pd.read_csv(metadata_path)

# Ensure barcodes are strings
filtered_barcodes = meta_df["Barcode"].astype(str).tolist()

# 2. Read the H5 file structure
h5_path = "GSE174367_snATAC-seq_filtered_peak_bc_matrix.h5"

with h5py.File(h5_path, "r") as f:
    # Read cell barcodes and feature names (decoding bytes to str if needed)
    h5_barcodes = f["/matrix/barcodes"][:].astype(str)
    feature_names = f["/matrix/features/name"][:].astype(str)

    # Reconstruct the sparse CSC matrix (Features x Cells)
    data = f["/matrix/data"][:]
    indices = f["/matrix/indices"][:]
    indptr = f["/matrix/indptr"][:]
    shape = f["/matrix/shape"][:]

    matrix = csc_matrix((data, indices, indptr), shape=shape)

# 3. Create a lookup mapping H5 cell barcodes to matrix column indices
barcode_to_idx = {bc: idx for idx, bc in enumerate(h5_barcodes)}

# Keep only metadata rows whose barcodes actually exist in the H5 matrix
valid_mask = meta_df["Barcode"].isin(barcode_to_idx)
meta_df_filtered = meta_df[valid_mask].copy()

# Get the ordered column indices matching the filtered CSV sequence
target_indices = [barcode_to_idx[bc] for bc in meta_df_filtered["Barcode"]]

# Slice and reorder columns in one operation
filtered_matrix = matrix[:, target_indices]

# 4. Save Outputs
# Save Matrix Market file (.mtx)
mmwrite("Morabito_counts.mtx", filtered_matrix)

# Save Feature names (one per line)
with open("Morabito_peaks.tsv", "w") as f:
    for name in feature_names:
        f.write(f"{name}\n")


meta_df_small = pd.DataFrame({"batch":meta_df["Sample.ID"].values, "cell_type":meta_df["Cell.Type"].values}, index=meta_df["Barcode"].values)
print(meta_df_small.cell_type.value_counts())

meta_df_small.to_csv("Morabito_metadata.tsv", sep="\t", index=True)

print(
    f"Successfully exported {filtered_matrix.shape[1]} cells and {filtered_matrix.shape[0]} features."
)
