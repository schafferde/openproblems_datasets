"""Build a per-dataset evaluation AnnData:
  .X                  = paired raw RNA counts (cells x genes), or empty (n_obs x 0) if none.
  .obsm['X_<method>'] = each method embedding, row-aligned to .obs by barcode.
  .obs                = the prep sidecar obs (cell_id index, barcode/sample/label), + 'has_rna'.

Paired RNA exists for Burkhardt / Weinand / Cheong as ready-made <Name>_RNA_counts.h5ad in
the collaborator's read-only tree; those were written with obs_names == the ATAC cell ids,
so the join is a straight barcode lookup. The ATAC-only datasets (Trevino / Granja / Morabito /
Kanemaru / Liang / Domcke_gs174k -- Domcke's "RNA" is ATAC-derived gene activity, not a
transcriptome) are absent from RNA_H5AD and get an empty .X.

Note: the ATAC matrix is intentionally NOT copied into .X or .layers -- it is large and
already available via manifest.load_dataset() straight from the source .mtx.

Usage: python build_anndata.py <dataset> [<dataset> ...]
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.manifest import load_manifest, build_cache  # noqa: E402

ROOT = os.environ.get("BR_ATAC_ROOT") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(ROOT, "outputs", "anndata")

RNA_H5AD = {                                   # dataset -> file in its (read-only) data_dir
    "Burkhardt": "Burkhardt_RNA_counts.h5ad",
    "Weinand":   "Weinand_RNA_counts.h5ad",
    "Cheong":    "Cheong_RNA_counts.h5ad",
}                                              # ATAC-only datasets (incl. Domcke, gene-activity
#   only -- no true scRNA-seq) are absent here and get an empty (n_obs, 0) .X.


def _rna(ds):
    """Paired RNA as (csr cells x genes, var, index of join keys). Raw counts in .X."""
    m = load_manifest(ds)
    a = ad.read_h5ad(os.path.join(m["data_dir"], RNA_H5AD[ds]))
    X = a.X.tocsr() if sparse.issparse(a.X) else sparse.csr_matrix(a.X)
    return X.astype(np.float32), a.var.copy(), pd.Index(a.obs_names.astype(str))


def build(ds):
    prep = build_cache(ds)                     # cheap; validates sidecars against the source
    obs = pd.read_csv(os.path.join(prep, "obs.tsv"), sep="\t", dtype=str,
                      keep_default_na=False).set_index("cell_id")   # "" not NaN -> h5ad-writable
    n = len(obs)

    obsm, uns = {}, {}
    mdir = os.path.join(ROOT, "outputs", ds)
    for meth in sorted(os.listdir(mdir)) if os.path.isdir(mdir) else []:
        ef = os.path.join(mdir, meth, "embedding.tsv.gz")
        if not os.path.isfile(ef):
            continue
        e = pd.read_csv(ef, sep="\t").set_index("barcode").reindex(obs.index)
        miss = int(e.isnull().any(axis=1).sum())
        if miss:
            print(f"  WARN {meth}: {miss}/{n} cells missing embedding (filled NaN)")
        obsm["X_" + meth] = e.values.astype(np.float32)
        vf = os.path.join(mdir, meth, "variance.tsv")   # scanpy-pca-style variance, if the method emits it
        if os.path.isfile(vf):
            v = pd.read_csv(vf, sep="\t")
            uns["X_" + meth] = {c: v[c].to_numpy() for c in
                                ("variance", "variance_ratio", "singular_value") if c in v}

    if ds in RNA_H5AD:
        M, var, keys = _rna(ds)
        pos = keys.get_indexer(obs.index)       # -1 where that cell has no RNA
        has = pos >= 0
        X = M[np.where(has, pos, 0)].multiply(has[:, None]).tocsr().astype(np.float32)
        obs["has_rna"] = has
        print(f"  RNA .X: {has.sum()}/{n} cells matched ({100 * has.mean():.1f}%), {M.shape[1]} genes")
    else:
        X, var = sparse.csr_matrix((n, 0), dtype=np.float32), pd.DataFrame()

    A = ad.AnnData(X=X, obs=obs, var=var)
    for k, v in obsm.items():
        A.obsm[k] = v
    for k, v in uns.items():
        A.uns[k] = v
    os.makedirs(OUT, exist_ok=True)
    A.write_h5ad(os.path.join(OUT, ds + ".h5ad"))
    print(f"[anndata] {ds}: X {A.n_obs}x{A.n_vars}, obsm={sorted(A.obsm)} -> outputs/anndata/{ds}.h5ad")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("usage: python build_anndata.py <dataset> [<dataset> ...]")
    for d in sys.argv[1:]:
        print(f"=== {d} ===")
        build(d)
