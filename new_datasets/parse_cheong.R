library(Matrix)


mat <- readRDS("GSE196987_ATAC_ReadCounts_Final_Multiome.rds")
writeMM(mat, "Cheong_counts.mtx")
writeLines(sub("-", ";", rownames(mat)), "Cheong_peaks.tsv")

df = read.table("GSE196987_ATAC_Metadata_Final_Multiome.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
df = transform(df, cell_label=rownames(df))
df_small = df[, c("cell_label", "MarkerAnnotations", "sample")]
colnames(df_small) = c("cell_label", "cell_type", "batch")
df_ordered <- df_small[match(colnames(mat), df$cell_label), ]

write.table(df_ordered, "Cheong_metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

mat2 <- readRDS("GSE196988_RNA_ReadCounts_Final_Multiome.rds")
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

write_h5ad(adata, "Cheong_RNA_counts.h5ad")


