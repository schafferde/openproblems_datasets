library(Matrix)
library(anndata)
library(dplyr)
library(readr)
library(stringr)
library(Seurat)

# 1. Read metadata (row.names = 1 handles the offset header automatically)
meta_raw <- read.table("GSE149683_Metadata_of_high_quality_cells_after_doublet_filtering.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

# 2. Clean, filter, and modify cell types
meta_clean <- meta_raw %>%
  select(cell, tissue, cell_type, donor_id, batch) %>%
  mutate(batch = paste(donor_id, batch, sep = "_")) %>%
  select(cell, tissue, cell_type, batch) %>%
  
  # --- Filter out "Unknown" ---
  filter(!str_detect(cell_type, "Unknown")) %>%
  
  # --- Clean up cell_type labels ---
  mutate(
    # Remove trailing '?' (e.g., "T-cell?" -> "T-cell")
    cell_type = str_remove(cell_type, "\\?$"),
    
    # Specific substring replacements
    cell_type = str_replace_all(cell_type, "Syncytiotrophoblast$", "Syncytiotrophoblasts"),
    cell_type = str_replace_all(cell_type, "Lymphoid and Myeloid", "Lymphoid/Myeloid"),
    cell_type = str_replace_all(cell_type, "interneurons", "neurons")
  ) %>%
  
  # Ensure no duplicate cell barcodes
  filter(!duplicated(cell))

rownames(meta_clean) <- meta_clean$cell

# 3. Split by tissue (as before)
meta_by_tissue <- split(meta_clean, meta_clean$tissue)

# 1. Setup output directories for chunks
dir.create("atac_chunks", showWarnings = FALSE)
dir.create("rna_chunks", showWarnings = FALSE)
dir.create("meta_chunks", showWarnings = FALSE)

# Vector of your Seurat object file paths
seurat_files <- list.files(path = ".", 
                           pattern = "_filtered\\.seurat\\.RDS$", 
                           full.names = TRUE)

peak_list <- list()
rna_list  <- list()
meta_ordered_list <- list()

for (file in seurat_files) {
  file_name <- basename(file)
  this_tissue <- sub("^GSM\\d{7}_(.*?)_filtered\\.seurat\\.RDS$", "\\1", file_name, ignore.case = TRUE)
  
  tissue_meta <- meta_by_tissue[[this_tissue]]
  if (is.null(tissue_meta)) next
  
  message(sprintf("[PROCESSING] Tissue: '%s' (%s)", this_tissue, file_name))
  
  obj <- readRDS(file)
  obj <- UpdateSeuratObject(obj)
  
  peak_mat <- GetAssayData(obj, assay = "peaks", layer = "counts")
  valid_cells <- intersect(colnames(peak_mat), tissue_meta$cell)
  
  if (length(valid_cells) == 0) {
    rm(obj, peak_mat)
    gc()
    next
  }
  
  # Align metadata and matrices to exact matching cell order
  tissue_meta_sub <- tissue_meta[match(valid_cells, tissue_meta$cell), ]
  peak_sub <- peak_mat[, valid_cells]
  
  # --- Save Metadata Chunk --- - this can actually be ignored since it's in the h5ad files as well
  saveRDS(tissue_meta_sub, file.path("meta_chunks", paste0(this_tissue, "_meta.rds")))
  
  # --- Save ATAC Peak Chunk as standard Matrix Market files ---
  tissue_atac_dir <- file.path("atac_chunks", this_tissue)
  dir.create(tissue_atac_dir, showWarnings = FALSE)
  
  writeMM(peak_sub, file.path(tissue_atac_dir, "matrix.mtx"))
  writeLines(colnames(peak_sub), file.path(tissue_atac_dir, "barcodes.tsv"))
  writeLines(rownames(peak_sub), file.path(tissue_atac_dir, "features.tsv"))
  #Apparently this is actually gene accessibility, not paired RNAseq
  # --- Save RNA Chunk (.h5ad) ---
  if ("RNA" %in% names(obj@assays)) {
    rna_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
    rna_sub <- rna_mat[, valid_cells]
    
    adata_chunk <- AnnData(
      X = t(rna_sub),
      obs = tissue_meta_sub,
      var = data.frame(row.names = rownames(rna_sub))
    )
    write_h5ad(adata_chunk, file.path("rna_chunks", paste0(this_tissue, "_rna.h5ad")))
    rm(rna_mat, rna_sub, adata_chunk)
  }
  
  rm(obj, peak_mat, peak_sub, tissue_meta_sub)
  gc()
}

message("\nAll R chunks successfully written to disk!")
quit(status = 0)

# 1. Reassemble global metadata (now cleanly ordered by tissue)
meta_final <- do.call(rbind, meta_ordered_list)

# 2. Combine Peak matrices & Export
combined_peaks <- do.call(cbind, peak_list)
stopifnot(identical(colnames(combined_peaks), meta_final$cell))

writeMM(combined_peaks, "Domcke_counts.mtx")
writeLines(rownames(combined_peaks), "Domcke_peaks.tsv")
write_tsv(meta_final, "Domcke_metadata.tsv")

# Free peak memory before building AnnData
rm(combined_peaks, peak_list)
gc()

# 3. Combine RNA matrices & Export to AnnData
library(anndata)

combined_rna <- do.call(cbind, rna_list)
stopifnot(identical(colnames(combined_rna), meta_final$cell))

adata_rna <- AnnData(
  X = t(combined_rna),
  obs = meta_final,
  var = data.frame(row.names = rownames(combined_rna))
)

write_h5ad(adata_rna, "Domcke_RNA_counts.h5ad")
