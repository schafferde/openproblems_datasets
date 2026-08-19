# Shared I/O for R methods: manifest + dataset loading + embedding output.
# The matrix is read straight from the collaborator's read-only .mtx (one source of truth,
# no derived copy). Only the two small sidecars -- features.tsv (peak ids normalized to
# chr_start_end) and obs.tsv (cell_id/barcode/sample/label) -- are built by
# src/common/manifest.py (build_cache); R just reads them.
suppressMessages({library(yaml); library(Matrix)})

BR_ROOT <- Sys.getenv("BR_ATAC_ROOT", unset = getwd())

load_manifest <- function(name) {
  p <- if (file.exists(name)) name else file.path(BR_ROOT, "config/datasets", paste0(name, ".yaml"))
  yaml::read_yaml(p)
}

# Returns list(counts = features x cells dgCMatrix, obs, batch_key, genome).
load_dataset <- function(name) {
  m <- load_manifest(name); out <- file.path(BR_ROOT, "data/prep", m$name)
  if (!file.exists(file.path(out, "obs.tsv")))
    stop("sidecars missing; run: python src/common/manifest.py prep ", name)
  src <- file.path(m$data_dir, m$matrix)                  # READ-ONLY source, never written
  # gzfile() reads plain .mtx transparently, so the same call covers both .mtx and .mtx.gz
  mat <- as(readMM(gzfile(src)), "CsparseMatrix")
  obs <- read.delim(file.path(out, "obs.tsv"), colClasses = "character")
  feats <- readLines(file.path(out, "features.tsv"))
  stopifnot(nrow(mat) == length(feats), ncol(mat) == nrow(obs))   # sidecars must match source
  rownames(mat) <- feats; colnames(mat) <- obs$cell_id
  list(counts = mat, obs = obs, batch_key = m$batch_key, genome = m$genome)
}

# Standard output: outputs/<dataset>/<method>/embedding.tsv.gz (barcode, dim1..dimN).
write_embedding <- function(name, method, cell_ids, emb) {
  m <- load_manifest(name); out <- file.path(BR_ROOT, "outputs", m$name, method)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  df <- data.frame(barcode = cell_ids, as.matrix(emb), check.names = FALSE)
  colnames(df)[-1] <- paste0("dim", seq_len(ncol(emb)))
  gz <- gzfile(file.path(out, "embedding.tsv.gz"), "w")
  write.table(df, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
  message(sprintf("[%s] %s: %d cells x %d dims -> %s", method, m$name, nrow(df), ncol(emb), out))
}
