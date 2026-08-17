library(SummarizedExperiment)
library(Matrix)

#See https://github.com/GreenleafLab/MPAL-Single-Cell-2019

data = readRDS("scATAC-Healthy-Hematopoiesis-191120.rds")

df = colData(data)
df = transform(df, CellLabel=rownames(df))
df = df[, c("CellLabel", "Group", "BioClassification")]
df$Cluster <- as.integer(sub("_.*", "", df$BioClassification))
map = read.delim("granja_cluster_map.tsv", header = TRUE)$scATAC
df$cell_type = map[df$Cluster]
filter = df$cell_type != "Unknown"

df = df[filter, ]
df$batch = df$Group
df = df[, c("CellLabel", "batch", "cell_type", "BioClassification")]
write.table(df, "Granja_metadata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

writeLines(sub("_", "-", sub("_", ":", rownames(data))), "Granja_peaks.tsv")

sparse_counts <- as(assay(data, "counts"), "dgTMatrix")
filter_counts = sparse_counts[, filter]
writeMM(filter_counts, file = "Granja_counts.mtx")
