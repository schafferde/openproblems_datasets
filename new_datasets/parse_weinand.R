library(Matrix)
mat <- readRDS("dataset1/dataset1_peak_matrix_nonzero.rds")
writeMM(mat, "Weinand_counts.mtx")
writeLines(rownames(mat), "Weinand_peaks.tsv")

#Metadata
df = readRDS("dataset1/dataset1_meta.rds")
df = transform(df, cell_label=rownames(df))
df_small = df[, c("cell_label", "sample", "phenotype")]
colnames(df_small) = c("cell_label", "batch", "cell_type")
df_ordered <- df_small[match(colnames(mat), df$cell_label), ]

write.table(df_ordered, "Weinand_metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

mat2 <- readRDS("dataset1/dataset1_gene_matrix.rds")
# Check that all column names match between the two matrices
stopifnot(all(colnames(mat) %in% colnames(mat2)))

# Subset/reorder mat2 columns to match mat1's column order
mat2_ordered <- mat2[, colnames(mat)]

# Verify column names are 100% identical and in order
identical(colnames(mat), colnames(mat2_ordered))

library(anndata)

var_df <- data.frame(row.names = rownames(mat2_ordered))
rownames(df_ordered) <- df_ordered$cell_label
adata <- AnnData(
  X = t(mat2_ordered),
  obs = df_ordered,
  var = var_df
)

write_h5ad(adata, "Weinand_RNA_counts.h5ad")


