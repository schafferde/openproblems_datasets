"""Shared I/O: manifest loading, sidecar prep, dataset loading, embedding output.

The datasets come pre-standardized by the collaborator as a read-only triplet
    <data_dir>/<Name>_counts.mtx      MatrixMarket, peaks x cells
    <data_dir>/<Name>_peaks.tsv       one peak id per row of the matrix
    <data_dir>/<Name>_metadata.tsv    one row per column of the matrix (order == column order)

We deliberately do NOT re-materialize the matrix. `build_cache` writes only two small
sidecars under data/prep/<name>/, and `load_dataset` reads the source .mtx directly:

  * one source of truth -- upstream re-exports (Granja was fixed mid-project, Kanemaru is
    still being iterated) are never shadowed by a stale local copy;
  * we do not pre-convert the .mtx. A derived dataset MAY drop a CSR <matrix>.npz sibling
    (large geosketch subsets, to skip a slow billion-nnz mmread); load_dataset prefers it
    ONLY when it is at least as new as the .mtx (mtime guard), else falls back to mmread.

The sidecars are what adapt the sources to the method scripts, and neither touches the matrix:
    features.tsv  peak ids normalized "chr:start-end" | "chr;start-end" -> "chr_start_end"
                  (run_signac.R regex-parses chr_start_end; run_pycistopic.py rsplit("_", 2))
    obs.tsv       cell_id, barcode, sample, label   (sample <- batch, label <- cell_type)

Usage: python manifest.py prep <dataset> [--force]
"""
import os
import re
import sys
import json
import yaml
import numpy as np
import pandas as pd
from scipy import sparse, io as sio

ROOT = os.environ.get("BR_ATAC_ROOT") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

PEAK_RE = re.compile(r"^(.+)[:;](\d+)-(\d+)$")     # chr may itself contain ':' / '_' / '-'


def load_manifest(name):
    p = name if os.path.exists(name) else os.path.join(ROOT, "config/datasets", f"{name}.yaml")
    with open(p) as fh:
        return yaml.safe_load(fh)


def cache_dir(name):
    return os.path.join(ROOT, "data/prep", load_manifest(name)["name"])


def source_paths(m):
    """(matrix, peaks, metadata) absolute paths in the collaborator's READ-ONLY tree."""
    d = m["data_dir"]
    return (os.path.join(d, m["matrix"]), os.path.join(d, m["peaks"]),
            os.path.join(d, m["metadata"]))


def _mtx_header(path):
    """(n_rows, n_cols, nnz) from the MatrixMarket banner, without reading the body."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            r, c, nnz = line.split()
            return int(r), int(c), int(nnz)
    raise ValueError(f"no dimension line in {path}")


def _stamp(path):
    st = os.stat(path)
    return {"path": path, "size": st.st_size, "mtime": int(st.st_mtime)}


def normalize_peaks(raw):
    """chr:start-end | chr;start-end -> chr_start_end. Raises on anything unparseable."""
    out, bad = [], []
    for p in raw:
        mt = PEAK_RE.match(p)
        if mt:
            out.append(f"{mt.group(1)}_{mt.group(2)}_{mt.group(3)}")
        else:
            bad.append(p)
            if len(bad) > 5:
                break
    if bad:
        raise ValueError(f"{len(bad)}+ unparseable peak ids, e.g. {bad[:5]}")
    return out


def build_cache(name, force=False):
    """Write the two sidecars for <name>. Cheap: never reads the matrix body."""
    m = load_manifest(name)
    out = os.path.join(ROOT, "data/prep", m["name"])
    mtx, peaks_f, meta_f = source_paths(m)
    stamp_f = os.path.join(out, "source.json")

    stamp = {"matrix": _stamp(mtx), "peaks": _stamp(peaks_f), "metadata": _stamp(meta_f)}
    if os.path.exists(stamp_f) and not force:
        with open(stamp_f) as fh:
            if json.load(fh) == stamp:                     # sources unchanged -> sidecars valid
                return out
        print(f"[prep] {m['name']}: source files changed since last prep -- rebuilding")

    n_peaks, n_cells, nnz = _mtx_header(mtx)

    raw = [l.rstrip("\n") for l in open(peaks_f) if l.strip()]
    if len(raw) != n_peaks:
        raise ValueError(f"{m['name']}: peaks.tsv has {len(raw)} rows, matrix has {n_peaks}")
    feats = normalize_peaks(raw)

    md = pd.read_csv(meta_f, sep="\t", dtype=str, keep_default_na=False)
    if len(md) != n_cells:
        raise ValueError(f"{m['name']}: metadata has {len(md)} rows, matrix has {n_cells} cols")
    for col in ("batch", "cell_type"):                     # named, since column ORDER varies
        if col not in md.columns:
            raise ValueError(f"{m['name']}: metadata lacks '{col}' (has {list(md.columns)})")
    bc = md.iloc[:, 0].astype(str)                         # barcode col is unnamed in some sets
    if not bc.is_unique:
        raise ValueError(f"{m['name']}: barcodes in column 0 are not unique")

    obs = pd.DataFrame({"cell_id": bc, "barcode": bc,
                        "sample": md["batch"].astype(str), "label": md["cell_type"].astype(str)})

    os.makedirs(out, exist_ok=True)
    pd.Series(feats).to_csv(os.path.join(out, "features.tsv"), index=False, header=False)
    obs.to_csv(os.path.join(out, "obs.tsv"), sep="\t", index=False)
    with open(stamp_f, "w") as fh:
        json.dump(stamp, fh, indent=2)
    print(f"[prep] {m['name']}: {n_peaks} peaks x {n_cells} cells (nnz={nnz}), "
          f"{obs['sample'].nunique()} batches / {obs['label'].nunique()} labels -> {out}")
    return out


def load_dataset(name):
    """AnnData (cells x features) with obs cell_id/barcode/sample/label. Reads the source .mtx."""
    import anndata as ad
    m = load_manifest(name)
    out = build_cache(name)                                # cheap; also validates the sidecars
    mtx, _, _ = source_paths(m)
    npz = os.path.splitext(mtx)[0] + ".npz"                # optional CSR fast-load sibling (peaks x cells)
    use_npz = os.path.exists(npz) and os.path.getmtime(npz) >= os.path.getmtime(mtx)  # never a stale copy
    mat = (sparse.load_npz(npz) if use_npz else sio.mmread(mtx)).tocsr()
    feats = pd.read_csv(os.path.join(out, "features.tsv"), header=None)[0].astype(str).tolist()
    obs = pd.read_csv(os.path.join(out, "obs.tsv"), sep="\t", dtype=str,
                      keep_default_na=False).set_index("cell_id")
    assert mat.shape == (len(feats), len(obs)), (mat.shape, len(feats), len(obs))
    A = ad.AnnData(X=mat.T.tocsr().astype(np.float32), obs=obs)
    A.var_names = feats
    return A


def write_embedding(name, method, cell_ids, emb):
    """Standard method output: outputs/<dataset>/<method>/embedding.tsv.gz (barcode, dim1..dimN)."""
    m = load_manifest(name)
    out = os.path.join(ROOT, "outputs", m["name"], method)
    os.makedirs(out, exist_ok=True)
    df = pd.DataFrame(np.asarray(emb), columns=[f"dim{i + 1}" for i in range(np.shape(emb)[1])])
    df.insert(0, "barcode", list(cell_ids))
    df.to_csv(os.path.join(out, "embedding.tsv.gz"), sep="\t", index=False)
    print(f"[{method}] {m['name']}: {df.shape[0]} cells x {df.shape[1] - 1} dims -> {out}")
    return out


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "prep":
        for d in [a for a in sys.argv[2:] if not a.startswith("--")]:
            build_cache(d, force="--force" in sys.argv)
    else:
        sys.exit("usage: python manifest.py prep <dataset> [<dataset> ...] [--force]")
