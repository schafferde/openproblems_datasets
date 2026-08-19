"""pycisTopic: LDA topic model -> cell-topic matrix (the embedding). No batch correction here.

Matrix mode (fragment files not required). pycisTopic's API shifts across versions;
adjust the create/model calls if they differ in the installed build.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.manifest import build_cache, load_dataset, write_embedding, cache_dir  # noqa: E402


def main(name, n_topics=(100,)):                            # 100 topics -> 100-dim embedding (heavy on wide data)
    # orig: n_topics=(10, 20, 30)  # candidate counts selected by evaluate_models
    build_cache(name)
    A = load_dataset(name)
    import matplotlib
    matplotlib.use("Agg")                                   # headless: evaluate_models may plot
    from pycisTopic.cistopic_class import create_cistopic_object
    from pycisTopic.lda_models import run_cgs_models_mallet, evaluate_models

    regions = [f"{c}:{s}-{e}" for c, s, e in (r.rsplit("_", 2) for r in A.var_names)]  # chr_start_end -> chr:start-end (chr may hold "_")
    counts = A.X.T.tocsr()                                  # regions x cells
    counts.data = counts.data.astype("int32")              # LDA (CGS) requires integer counts
    obj = create_cistopic_object(counts, cell_names=list(A.obs_names), region_names=regions)
    os.environ.setdefault("_JAVA_OPTIONS", "-Xmx64g")       # Java heap; mallet wrapper hardcodes -Xmx1g, _JAVA_OPTIONS overrides
    n_iter = int(os.environ.get("BR_LDA_ITER", 150))        # optional cap (smoke); default 150
    n_cpu = int(os.environ.get("BR_MALLET_CPU", 32))        # Mallet is multithreaded (CPU, not GPU)
    tmp = os.path.join(cache_dir(name), "mallet_tmp")       # per-dataset (Mallet corpus.txt is else /tmp, collides)
    os.makedirs(tmp, exist_ok=True)
    models = run_cgs_models_mallet(obj, n_topics=list(n_topics), n_cpu=n_cpu, n_iter=n_iter,
                                   random_state=0, tmp_path=tmp)
    model = evaluate_models(models, return_model=True) if len(models) > 1 else models[0]  # select only if >1 candidate
    obj.add_LDA_model(model)

    ct = obj.selected_model.cell_topic.T                    # cells x topics = the embedding (index=cell_id)
    ids = [c.rsplit("___", 1)[0] for c in ct.index]         # strip pycisTopic's "___cisTopic" tag so barcodes match other methods
    write_embedding(name, "pycistopic", ids, ct.to_numpy())


if __name__ == "__main__":
    main(sys.argv[1])
