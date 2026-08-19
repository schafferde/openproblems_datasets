"""peakVI (scvi-tools): batch-corrected latent embedding from a peak/window count matrix."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.manifest import load_manifest, build_cache, load_dataset, write_embedding  # noqa: E402


def main(name):
    m = load_manifest(name)
    build_cache(name)
    A = load_dataset(name)
    import scanpy as sc
    sc.pp.filter_genes(A, min_cells=max(1, int(0.01 * A.n_obs)))  # drop peaks in <1% of cells
    A.X.data[:] = (A.X.data > 0).astype(A.X.dtype)          # PEAKVI expects binary accessibility
    import scvi
    scvi.model.PEAKVI.setup_anndata(A, batch_key=m["batch_key"])
    model = scvi.model.PEAKVI(A, n_latent=100)              # 100-dim latent
    # orig: model = scvi.model.PEAKVI(A)                    # default n_latent auto (smaller)
    ep = os.environ.get("BR_MAX_EPOCHS")                    # optional cap (CPU/smoke); default = scvi heuristic
    model.train(max_epochs=int(ep) if ep else None)
    write_embedding(name, "peakvi", A.obs_names, model.get_latent_representation())


if __name__ == "__main__":
    main(sys.argv[1])
