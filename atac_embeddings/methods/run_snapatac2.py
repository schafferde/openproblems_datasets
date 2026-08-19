"""snapATAC2: spectral (Laplacian eigenmap) embedding from a count matrix. No batch correction."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.manifest import load_manifest, build_cache, load_dataset, write_embedding  # noqa: E402


def main(name):
    build_cache(name)
    A = load_dataset(name)
    import snapatac2 as snap
    snap.pp.select_features(A)                              # most-accessible features -> var['selected']
    nc = int(os.environ.get("BR_N_COMPS", 100))            # request 100 dims (spectral may return fewer)
    snap.tl.spectral(A, n_comps=nc)                        # spectral embedding (Laplacian eigenmap, cosine)
    write_embedding(name, "snapatac2", A.obs_names, A.obsm["X_spectral"])   # no batch correction here


if __name__ == "__main__":
    main(sys.argv[1])
