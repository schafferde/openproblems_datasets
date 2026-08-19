"""Raw-data LSI baseline -> outputs/<ds>/TF-IDF_SVD/ (build_anndata puts it in .obsm['X_TF-IDF_SVD']).

TF-IDF on the raw peak counts, then TruncatedSVD to 101 components, keep components 2:101 (100
dims; the 1st is dropped as it tracks sequencing depth). No batch correction, no feature
selection. TF-IDF uses the log-IDF variant idf = log(1 + N/(df+1)); the TF-IDF nonzeros are
then log1p(tf * idf * 1e4) before SVD.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.manifest import build_cache, load_dataset, write_embedding  # noqa: E402


def main(name, dims=101):                                  # 101 comps -> drop 1 -> 100-dim embedding
    build_cache(name)
    A = load_dataset(name)                                 # cells x peaks
    import numpy as np
    from scipy import sparse
    from sklearn.decomposition import TruncatedSVD

    X = A.X.tocsr().astype(np.float32)                      # cells x peaks
    n_cells = X.shape[0]
    cell_sums = np.asarray(X.sum(1)).ravel()               # per-cell depth (TF denom)
    peak_sums = np.asarray(X.sum(0)).ravel()               # per-peak prevalence (IDF denom)
    idf = np.log(1.0 + n_cells / (peak_sums + 1.0)).astype(np.float32)   # +1 guards empty peaks
    # row-scale by 1/depth (TF), col-scale by idf; max(.,1) guards empty cells. .astype keeps f32.
    tfidf = (sparse.diags(1.0 / np.maximum(cell_sums, 1.0)).astype(np.float32) @ X) \
        @ sparse.diags(idf).astype(np.float32)
    tfidf = tfidf.tocsr()
    tfidf.data = np.log1p(tfidf.data * 1e4).astype(np.float32)   # log on nonzeros -> sparsity kept

    svd = TruncatedSVD(n_components=dims, random_state=0, algorithm="randomized", n_iter=7)
    emb = svd.fit_transform(tfidf)[:, 1:]                  # cells x (dims-1); drop depth comp 1
    out = write_embedding(name, "TF-IDF_SVD", A.obs_names, emb)

    # scanpy-PCA-style explained-variance sidecar (slice [1:] to match the dropped depth comp).
    # build_anndata reads this into A.uns['X_TF-IDF_SVD'] = {variance, variance_ratio, singular_value}.
    import pandas as pd
    pd.DataFrame({
        "dim": np.arange(1, emb.shape[1] + 1),
        "variance": svd.explained_variance_[1:],
        "variance_ratio": svd.explained_variance_ratio_[1:],
        "singular_value": svd.singular_values_[1:],
    }).to_csv(os.path.join(out, "variance.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main(sys.argv[1])
