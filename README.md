# Data Processing for Batch Integration with BatchRefiner
## Description

This branch of this fork of `openprobems/datasets` contains code and materials used to add three additional scRNA-seq and nine additional scATAC-seq benchmarking datasets to the Openproblems `task_batch_integration` pipeline for use with BatchRefiner.

Standalone package - [BatchRefiner](https://github.com/schafferde/BatchRefiner)

Main reproducibility repository - [https://github.com/schafferde/task_batch_integration/tree/batchrefiner_reproducibility](https://github.com/schafferde/task_batch_integration/tree/batchrefiner_reproducibility)

Schäffer, D. E, Kang, H., Aksu, E. D., Edelman, D., Berger, B.: BatchRefiner: fast, significant improvement in batch integration of single-cell embeddings with ensemble refinement. *In preparation*

## Data
``new_datasets`` contains fourteen scripts and three text files used to parse datasets into `.h5ad` (scRNA-seq) or MatrixMarket-like (scATAC-seq) format from raw counts:
- The script `parse_ocular.py` parses the Ocular Atlas scRNA-seq dataset. Raw data consist of one `.mtx` file and two `.tsv` files downloaded from [Single Cell Portal accession SCP2310](https://singlecell.broadinstitute.org/single_cell/study/SCP2310/) (account required for download). 
- The script `parse_muris.py` parses the Tabula Muris scRNA-seq dataset. Raw data consist of two `.zip` archives and four `.csv` files downloaded from Figshare, [dataset 5968960 for droplet](https://doi.org/10.6084/m9.figshare.5968960) and [dataset 5829687 for FACS](https://doi.org/10.6084/m9.figshare.5829687). Each method is associated with one archive and two metadata files. The archives must first be extracted. 
- The script `parse_celegans.py` parses the C. elegans Embryo scRNA-seq dataset. Raw data consist of one `.txt` matrix file and two `.csv` files downloaded from [NCBI GEO accession GSE126954](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126954).
- The file `celegans_cell_type_map.csv` contains a human- and machine-readable mapping of cell type labels for the C. elegans dataset into a new, consistent form for differentiated cells. The mapping logic is described in a comment in `parse_celegans.py`. These labels were based on the notes in the supplementary information of the original study.
- The script `add_emsembl_robust.py` was used to add Ensembl IDs to many gene names in the Ocular Atlas and Tabula Muris datasets. Ensembl IDs are expected by the OpenProblems pipeline because a few downstream methods expect them. However, none of the methods we used expect Ensembl IDs, so use of this script is not strictly necessary. 
- The script `parse_burkhardt.py` parses the Burkhardt scATAC-seq dataset. Raw data consist of one `.h5ad` file downloaded from [NCBI GEO accession GSE19412](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19412).
- The script `parse_cheong.R` parses the Cheong scATAC-seq dataset. Raw data consist of two `.rds` data files and one `.txt` metadata file downloaded from [NCBI GEO accession GSE196987](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196987).
- The script `parse_domcke.R` parses the Domcke scATAC-seq dataset. Raw data consist of fifteen `.rds` data files and one `.txt` metadata file downloaded from [NCBI GEO accession GSE149683](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149683).
- The script `parse_granja.R` parses the Granja scATAC-seq dataset. Raw data consist of one `.rds` file downloaded from [github.com/GreenleafLab/MPAL-Single-Cell-2019](https://github.com/GreenleafLab/MPAL-Single-Cell-2019).
- The file `granja_cluster_map.tsv` contains a human- and machine-readable mapping of cluster labels to cell type labels for the Granja dataset. This mapping is based on the supplementary information of the original study.
- The script `parse_kanemaru.py` parses the Kanemaru scATAC-seq dataset. Raw data consist of one `.h5ad` file downloaded from [www.heartcellatlas.org/](https://www.heartcellatlas.org/).
- The script `parse_liang.py` parses the Liang scATAC-seq dataset, including reversing the normalization that was applied. Raw data consist of one `.h5ad` file downloaded from [cells.ucsc.edu/?ds=retina-atac](https://cells.ucsc.edu/?ds=retina-atac).
- The script `parse_morabito.py` parses the Morabito scATAC-seq dataset. Raw data consist of one `.h5` data file and one `.csv` metadata file downloaded from [NCBI GEO accession GSE174367](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367).
- The script `parse_trevino.py` parses the Trevino scATAC-seq dataset. Raw data consist of one `.tsv.gz` data file, and one `.txt` metadata file, and one `.bed` peak list downloaded from [NCBI GEO accession GSE162170](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162170).
- The script `tsv_to_mm.sh` provies a more memory-efficient conversion of a plain-text TSV matrix to compressed (MatrixMarket) format, and is used for the Trevino dataset.
- The file `trevino_cluster_map.tsv` contains a human- and machine-readable mapping of cluster labels to cell type labels for the Trevino dataset. This mapping is based on the supplementary information of the original study.
- The script `parse_weinand.R` parses the Granja scATAC-seq dataset. Raw data consist of the `snATAC_benchmark_inputs_dataset1.zip` archive downloaded (and extracted) from [Zenodo record 15073137](https://zenodo.org/records/15073137).


### Dataset References for New Datasets
Monavarfeshani, A., Yan, W., Pappas, C., *et al*: Transcriptomic analysis of the ocular posterior segment completes a cell atlas of the human eye. *Proc. Natl. Acad. Sci. U. S. A.* **120**(34), e2306153120 (Aug 2023)

Tabula Muris Consortium: Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris. *Nature* **562**(7727), 367-372 (Oct 2018)

Packer, J. S., Zhu, Q., Huynh, C., *et al*: A lineage-resolved molecular atlas of C. elegans embryogenesis at single-cell resolution. *Science* **365**(6459), eaax1971 (Sep 2019)

Luecken, M., Burkhardt, D., Cannoodt, R., *et al*: A sandbox for prediction and integration of DNA, RNA, and proteins in single cells. In: Vanschoren, J., Yeung, S. (eds.) *Proceedings of the Neural Information Processing Systems Track on Datasets and Benchmarks*. **vol. 1** (2021)

Cheong, J.G., Ravishankar, A., Sharma, S., *et al*: Epigenetic memory of coronavirus infection in innate immune cells and their progenitors. *Cell* **186**(18), 3882-3902.e24 (Aug 2023)

Domcke, S., Hill, A.J., Daza, R.M., *et al*: A human cell atlas of fetal chromatin accessibility. *Science* **370**(6518), eaba7612 (Nov 2020)

Granja, J.M., Klemm, S., McGinnis, L.M., *et al*: Single-cell multiomic analysis identifies regulatory programs in mixed-phenotype acute leukemia. *Nat. Biotechnol.* **37**(12), 1458-1465 (Dec 2019)

Kanemaru, K., Cranley, J., Muraro, D.,*et al*: Spatially resolved multiomics of human cardiac niches. *Nature* **619**(7971), 801-810 (Jul 2023)

Liang, Q., Cheng, X., Wang, J., *et al*: A multi-omics atlas of the human retina at single-cell resolution. *Cell Genom.* **3**(6), 100298 (Jun 2023)

Morabito, S., Miyoshi, E., Michael, N., *et al*: Single-nucleus chromatin accessibility and transcriptomic characterization of Alzheimer’s disease. *Nat. Genet.* **53**(8), 1143-1155 (Aug 2021)

Trevino, A.E., Müller, F., Andersen, J., *et al*: Chromatin and gene-regulatory dynamics of the developinghuman cerebral cortex at single-cell resolution. *Cell* **184**(19), 5053-5069.e23 (Sep 2021)

Weinand, K., Langan, E.M., Curtis, M., Raychaudhuri, S.: Defining effective strategies to integrate multi-sample single-nucleus ATAC-seq datasets via a multimodal-guided approach. *bioRxiv* p. 2025.04.02.646871 (Apr 2025)

##  scATAC-seq Embeddings
The directory `atac_embeddings` contains the scripts we used to generate baseline embeddings of the nine scATAC-seq datasets. 

##  Modifications for OpenProblems Pipeline
- We added a new proprocessing workflow and corresponding data loader to process local scRNA-seq datasets in `h5ad` format into OpenProblems common format. This workflow also preserves `obsm` fields and permits empty data fields, thereby supporting scATAC-seq datasets with pre-computed embeddings. 
  - These can be found at `src/workflows/scrnaseq/process_local_h5ad/` and `src/loaders/scrnaseq/local_h5ad`
  - Because of limitations in read permissions when running via docker, they require each `h5ad` file to be manually copied into the loader's working directory when running. 
- The configuration files `celegans_config.yaml`, `muris_config.yaml`, and `ocular_config.yaml` contain the parameters used to process our three added scRNA-seq datasets through this pipeline. The nine confirgration files in `atac_configs` contain the parameters used to process the nine scATAC-seq datasets with pre-computed embeddings.
- The script `./run_local.sh` runs the data processing pipeline on datasets whose configs are described by the `input_states` field, subject to the input `h5ad` files being copied to the working directory of each loader process as mentioned above.
  - The script `./build_all_docker_containers.sh` builds just the docker containers needed for processing local scRNA-seq (or scATAC-seq) datasets.

---
## The original README from the OpenProblems repository follows below.

# openproblems datasets

This repository contains dataset loaders and processing workflows.

- <a href="#common-datasets" id="toc-common-datasets">Common datasets</a>
  - <a href="#pipeline-topology" id="toc-pipeline-topology">Pipeline
    topology</a>
  - <a href="#file-format-api" id="toc-file-format-api">File format API</a>
    - <a href="#datasetpcahvg"
      id="toc-datasetpcahvg"><code>Dataset+Pca+Hvg</code></a>
    - <a href="#normalized-dataset"
      id="toc-normalized-dataset"><code>Normalized Dataset</code></a>
    - <a href="#datasetpca" id="toc-datasetpca"><code>Dataset+Pca</code></a>
    - <a href="#raw-dataset" id="toc-raw-dataset"><code>Raw Dataset</code></a>
  - <a href="#component-api" id="toc-component-api">Component API</a>
    - <a href="#dataset-loader"
      id="toc-dataset-loader"><code>Dataset Loader</code></a>
    - <a href="#normalization"
      id="toc-normalization"><code>Normalization</code></a>
    - <a href="#processor-hvg"
      id="toc-processor-hvg"><code>Processor Hvg</code></a>
    - <a href="#processor-pca"
      id="toc-processor-pca"><code>Processor Pca</code></a>

## Pipeline topology

``` mermaid
%%| column: screen-inset-shaded
flowchart LR
  file_dataset(Dataset+Pca+Hvg)
  file_normalized(Normalized Dataset)
  file_pca(Dataset+Pca)
  file_raw(Raw Dataset)
  comp_dataset_loader[/Dataset Loader/]
  comp_normalization[/Normalization/]
  comp_processor_hvg[/Processor Hvg/]
  comp_processor_pca[/Processor Pca/]
  file_raw---comp_normalization
  file_pca---comp_processor_hvg
  file_normalized---comp_processor_pca
  comp_dataset_loader-->file_raw
  comp_normalization-->file_normalized
  comp_processor_hvg-->file_dataset
  comp_processor_pca-->file_pca
```

## File format API

### `Dataset+Pca+Hvg`

A normalised data with a PCA embedding and HVG selection

Used in:

- [processor hvg](#processor%20hvg): output (as output)

Slots:

| struct | name             | type    | description                                                             |
|:-------|:-----------------|:--------|:------------------------------------------------------------------------|
| layers | counts           | integer | Raw counts                                                              |
| layers | normalized       | double  | Normalised expression values                                            |
| obs    | celltype         | string  | Cell type information                                                   |
| obs    | batch            | string  | Batch information                                                       |
| obs    | tissue           | string  | Tissue information                                                      |
| obs    | size_factors     | double  | The size factors created by the normalisation method, if any.           |
| var    | hvg              | boolean | Whether or not the feature is considered to be a ‘highly variable gene’ |
| var    | hvg_score        | integer | A ranking of the features by hvg.                                       |
| obsm   | X_pca            | double  | The resulting PCA embedding.                                            |
| varm   | pca_loadings     | double  | The PCA loadings matrix.                                                |
| uns    | dataset_id       | string  | A unique identifier for the dataset                                     |
| uns    | normalization_id | string  | Which normalization was used                                            |
| uns    | pca_variance     | double  | The PCA variance objects.                                               |

Example:

    AnnData object
     obs: 'celltype', 'batch', 'tissue', 'size_factors'
     var: 'hvg', 'hvg_score'
     uns: 'dataset_id', 'normalization_id', 'pca_variance'
     obsm: 'X_pca'
     varm: 'pca_loadings'
     layers: 'counts', 'normalized'

### `Normalized Dataset`

A normalized dataset

Used in:

- [normalization](#normalization): output (as output)
- [processor pca](#processor%20pca): input (as input)

Slots:

| struct | name             | type    | description                                                   |
|:-------|:-----------------|:--------|:--------------------------------------------------------------|
| layers | counts           | integer | Raw counts                                                    |
| layers | normalized       | double  | Normalised expression values                                  |
| obs    | celltype         | string  | Cell type information                                         |
| obs    | batch            | string  | Batch information                                             |
| obs    | tissue           | string  | Tissue information                                            |
| obs    | size_factors     | double  | The size factors created by the normalisation method, if any. |
| uns    | dataset_id       | string  | A unique identifier for the dataset                           |
| uns    | normalization_id | string  | Which normalization was used                                  |

Example:

    AnnData object
     obs: 'celltype', 'batch', 'tissue', 'size_factors'
     uns: 'dataset_id', 'normalization_id'
     layers: 'counts', 'normalized'

### `Dataset+Pca`

A normalised data with a PCA embedding

Used in:

- [processor hvg](#processor%20hvg): input (as input)
- [processor pca](#processor%20pca): output (as output)

Slots:

| struct | name             | type    | description                                                   |
|:-------|:-----------------|:--------|:--------------------------------------------------------------|
| layers | counts           | integer | Raw counts                                                    |
| layers | normalized       | double  | Normalised expression values                                  |
| obs    | celltype         | string  | Cell type information                                         |
| obs    | batch            | string  | Batch information                                             |
| obs    | tissue           | string  | Tissue information                                            |
| obs    | size_factors     | double  | The size factors created by the normalisation method, if any. |
| obsm   | X_pca            | double  | The resulting PCA embedding.                                  |
| varm   | pca_loadings     | double  | The PCA loadings matrix.                                      |
| uns    | dataset_id       | string  | A unique identifier for the dataset                           |
| uns    | normalization_id | string  | Which normalization was used                                  |
| uns    | pca_variance     | double  | The PCA variance objects.                                     |

Example:

    AnnData object
     obs: 'celltype', 'batch', 'tissue', 'size_factors'
     uns: 'dataset_id', 'normalization_id', 'pca_variance'
     obsm: 'X_pca'
     varm: 'pca_loadings'
     layers: 'counts', 'normalized'

### `Raw Dataset`

An unprocessed dataset as output by a dataset loader.

Used in:

- [dataset loader](#dataset%20loader): output (as output)
- [normalization](#normalization): input (as input)

Slots:

| struct | name       | type    | description                         |
|:-------|:-----------|:--------|:------------------------------------|
| layers | counts     | integer | Raw counts                          |
| obs    | celltype   | string  | Cell type information               |
| obs    | batch      | string  | Batch information                   |
| obs    | tissue     | string  | Tissue information                  |
| uns    | dataset_id | string  | A unique identifier for the dataset |

Example:

    AnnData object
     obs: 'celltype', 'batch', 'tissue'
     uns: 'dataset_id'
     layers: 'counts'

## Component API

### `Dataset Loader`

Arguments:

| Name       | Type                          | Direction | Description                                           |
|:-----------|:------------------------------|:----------|:------------------------------------------------------|
| `--output` | [Raw Dataset](#Raw%20dataset) | output    | An unprocessed dataset as output by a dataset loader. |

### `Normalization`

Arguments:

| Name                 | Type                                        | Direction | Description                                                  |
|:---------------------|:--------------------------------------------|:----------|:-------------------------------------------------------------|
| `--input`            | [Raw Dataset](#Raw%20dataset)               | input     | An unprocessed dataset as output by a dataset loader.        |
| `--output`           | [Normalized Dataset](#Normalized%20dataset) | output    | A normalized dataset                                         |
| `--layer_output`     | `string`                                    | input     | The name of the layer in which to store the normalized data. |
| `--obs_size_factors` | `string`                                    | input     | In which .obs slot to store the size factors (if any).       |

### `Processor Hvg`

Arguments:

| Name              | Type                                | Direction | Description                                                                |
|:------------------|:------------------------------------|:----------|:---------------------------------------------------------------------------|
| `--input`         | [Dataset+Pca](#Dataset+PCA)         | input     | A normalised data with a PCA embedding                                     |
| `--layer_input`   | `string`                            | input     | Which layer to use as input for the PCA.                                   |
| `--output`        | [Dataset+Pca+Hvg](#Dataset+PCA+HVG) | output    | A normalised data with a PCA embedding and HVG selection                   |
| `--var_hvg`       | `string`                            | input     | In which .var slot to store whether a feature is considered to be hvg.     |
| `--var_hvg_score` | `string`                            | input     | In which .var slot to store whether a ranking of the features by variance. |
| `--num_features`  | `integer`                           | input     | The number of HVG to select                                                |

### `Processor Pca`

Arguments:

| Name               | Type                                        | Direction | Description                                                                                                          |
|:-------------------|:--------------------------------------------|:----------|:---------------------------------------------------------------------------------------------------------------------|
| `--input`          | [Normalized Dataset](#Normalized%20dataset) | input     | A normalized dataset                                                                                                 |
| `--layer_input`    | `string`                                    | input     | Which layer to use as input for the PCA.                                                                             |
| `--output`         | [Dataset+Pca](#Dataset+PCA)                 | output    | A normalised data with a PCA embedding                                                                               |
| `--obsm_embedding` | `string`                                    | input     | In which .obsm slot to store the resulting embedding.                                                                |
| `--varm_loadings`  | `string`                                    | input     | In which .varm slot to store the resulting loadings matrix.                                                          |
| `--uns_variance`   | `string`                                    | input     | In which .uns slot to store the resulting variance objects.                                                          |
| `--num_components` | `integer`                                   | input     | Number of principal components to compute. Defaults to 50, or 1 - minimum dimension size of selected representation. |
