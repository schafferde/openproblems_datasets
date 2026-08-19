# Signac / Seurat: emits two embeddings from one load of the dataset:
#   signac       - raw LSI (TF-IDF + SVD, drop comp 1). No batch correction.
#   seurat_rlsi  - reciprocal-LSI anchor integration of the LSI across samples (ATAC-standard).
# Base assay/TF-IDF/LSI computed once. Integration wrapped in tryCatch so a failure doesn't
# lose the signac embedding already written.
suppressMessages({library(Signac); library(Seurat)})
.self <- normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]))
source(file.path(dirname(.self), "..", "common", "manifest.R"))

name <- commandArgs(trailingOnly = TRUE)[1]
ds <- load_dataset(name)
f <- rownames(ds$counts)                       # chr_start_end -> GRanges (chr may contain "_")
gr <- GenomicRanges::GRanges(sub("_[0-9]+_[0-9]+$", "", f),
                             IRanges::IRanges(as.numeric(sub(".*_([0-9]+)_[0-9]+$", "\\1", f)) + 1,
                                              as.numeric(sub(".*_([0-9]+)$", "\\1", f))))
obj <- CreateSeuratObject(CreateChromatinAssay(counts = ds$counts, ranges = gr), assay = "peaks",
                          meta.data = `rownames<-`(ds$obs, ds$obs$cell_id))
obj <- RunTFIDF(obj); obj <- FindTopFeatures(obj, min.cutoff = "q0"); obj <- RunSVD(obj, n = 101)

# --- signac: raw LSI, drop comp 1 (depth) ---
emb <- Embeddings(obj, "lsi")[, 2:101]
write_embedding(name, "signac", rownames(emb), emb)

# per-sample LSI, used only to find the reciprocal-LSI anchors (dims 2:30 per the Signac vignette)
parts <- SplitObject(obj, split.by = ds$batch_key)
parts <- lapply(parts, function(x) RunSVD(FindTopFeatures(RunTFIDF(x), min.cutoff = "q0"), n = 50))

# --- seurat_rlsi: integrate Signac's LSI (obj[["lsi"]]) across samples via rLSI anchors ---
# Output name encodes the k.weight actually used: seurat_rlsi_kw<k>.
tryCatch({
  afeat <- head(VariableFeatures(obj), 25000)              # cap anchor features (top-accessible) for tractability
  a <- FindIntegrationAnchors(object.list = parts, anchor.features = afeat,
                              reduction = "rlsi", dims = 2:30)
  # k.weight: requested (default 100) vs the sparsest per-dataset anchor count.
  # Default: silently reduce k.weight to the min anchor count (name records what was used).
  # BR_RLSI_STRICT=1 -> instead skip the dataset when it can't support the requested k.weight.
  kw_req <- as.integer(Sys.getenv("BR_RLSI_KWEIGHT", "100"))
  min_anch <- min(table(c(slot(a, "anchors")$dataset1, slot(a, "anchors")$dataset2)))
  strict <- Sys.getenv("BR_RLSI_STRICT", "0") == "1"
  if (strict && min_anch < kw_req) {
    message("skipped seurat_rlsi: min anchors ", min_anch, " < requested k.weight ", kw_req)
  } else {
    kw <- min(kw_req, min_anch)
    repeat {
      ri <- tryCatch(IntegrateEmbeddings(anchorset = a, reductions = obj[["lsi"]],  # input = the signac LSI output
                                         new.reduction.name = "integrated_lsi",
                                         dims.to.integrate = 2:101, k.weight = kw),
                     error = function(e) if (!strict && grepl("k.weight", conditionMessage(e)) && kw > 2) NULL else stop(e))
      if (!is.null(ri) || kw <= 2) break
      kw <- max(2, kw %/% 2L)                              # still too high somewhere: halve and retry
    }
    e <- Embeddings(ri, "integrated_lsi")
    write_embedding(name, paste0("seurat_rlsi_kw", kw), rownames(e), e) # prior kw runs kept (distinct names)
  }
}, error = function(err) message("seurat_rlsi failed: ", conditionMessage(err)))
