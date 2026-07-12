# Gene-centric machine learning for osteogenesis-associated genes

Analysis code for the study in which a **gene-centric** machine-learning framework is used to
identify genes associated with osteoblast differentiation, with **ENPEP** as the top candidate
under simulated microgravity.

In this framework the usual roles are transposed: **genes are the instances to be classified**
(labelled *positive* = osteoblast differentiation / *negative* = adipogenesis, from Gene Ontology),
and the **expression profile of each gene across 223 biological samples (34 datasets) is the feature
vector**. Because every gene is measured across the same set of samples, sample-level batch effects
are confined to the feature dimension and do not require cross-sample correction (see the manuscript
Methods/Discussion).

---

## Repository layout

```
gene-centric-osteogenesis-ML/
├── data/                          # inputs (see data/README.md)
├── 01_original_four_models/       # exploratory 4-model analysis
│   ├── knn_optimise_k.R           #   scan k for k-NN over repeated resampling
│   ├── knn_k116.R                 #   k-NN fixed at the selected k = 116
│   ├── svm.R                      #   SVM training (kernlab::ksvm, RBF kernel)
│   ├── svm_roc_analysis.R         #   ROC/threshold analysis of the SVM results
│   ├── logistic_cubist.R          #   logistic-regression (Cubist) model
│   └── random_forest.R            #   random-forest model
├── 02_model_benchmarking/
│   └── caret_15_models.R          # expanded benchmarking of 15 models (caret)
├── 03_best_model_selection/
│   └── select_best_models.R       # pick the best cycle/model by AUC; ROC plots
├── 04_enpep_prediction/
│   └── enpep_prediction.R         # per-model ENPEP prediction  →  Fig 3j
├── 05_figures/
│   ├── fig2_performance_comparison.R   # model-performance figures  →  Fig 2
│   ├── figS4_pca_batch.R               # PCA coloured by dataset     →  Fig S4
│   └── figS5_heldout_violin.R          # held-out P(OB) violin       →  Fig S5
└── outputs/                       # generated results/figures (git-ignored)
```

The numeric prefixes indicate the order in which the scripts are intended to run.

---

## Data

The scripts expect two input files in `data/` (details in [`data/README.md`](data/README.md)):

| File | Description |
|------|-------------|
| `log2_dataset_frame.csv` | Merged expression matrix, genes × 223 samples, log2-transformed. **Not distributed here** (≈26 MB, derived entirely from public datasets); rebuild from the accessions listed in the manuscript Supplementary Table, or obtain from the authors. |
| `Label_genes.csv` | Two columns, `Symbol` and `Label` (1 = positive/osteoblast differentiation, 0 = negative/adipogenesis). |

All expression data originate from publicly available repositories (GEO / ArrayExpress); the complete
list of accession numbers is provided in the manuscript Supplementary Table. Only a log2 transformation
was applied; no batch correction was performed (see manuscript).

---

## Requirements

- **R ≥ 4.2**
- Core: `caret`, `pROC`, `randomForest`, `kernlab`, `class`, `Cubist`, `dplyr`, `stringr`, `tidyr`
- Figures: `ggplot2`, `scales`, `pals` (optional palette), `openxlsx`
- Individual `caret` model back-ends (installed on demand by `02_model_benchmarking/caret_15_models.R`):
  `gbm`, `e1071`, `xgboost`, `earth`, `klaR`, `pls`, `ada`, `ranger`, `Rborist`, `wsrf`, `rotationForest`, …

```r
install.packages(c("caret","pROC","randomForest","kernlab","class","Cubist",
                   "dplyr","stringr","tidyr","ggplot2","scales","pals","openxlsx"))
```

---

## How to reproduce

Run the scripts from the **repository root** (each reads `data/…` and writes to `outputs/…`; adjust
paths if your layout differs).

1. **Original four-model analysis** — `01_original_four_models/`
   Exploratory feasibility test that the gene-centric matrix is learnable. Each classifier is
   evaluated by repeated random 80/20 gene-level hold-out; per-cycle predictions (including the
   held-out probabilities and the ENPEP prediction) are written to `OB_Model_result_frame_<model>.csv`.
   - `knn_optimise_k.R` scans *k* and selects **k = 116**; `knn_k116.R` runs the fixed-*k* model.
   - `svm.R`, `logistic_cubist.R`, `random_forest.R` train the SVM, LR and RF models.
   - `svm_roc_analysis.R` produces the ROC/threshold analysis from the SVM result table.

2. **15-model benchmarking** — `02_model_benchmarking/caret_15_models.R`
   Re-implements the pipeline in `caret` and evaluates 15 models under identical 80/20 hold-out
   partitions (1,000 repeats), writing one `OB_Model_result_frame_<model>.csv` per model and saving
   fitted models to `cycles/`.

3. **Best-model selection** — `03_best_model_selection/select_best_models.R`
   Selects, per model, the best cycle by AUC and draws the corresponding ROC curves.

4. **ENPEP prediction** — `04_enpep_prediction/enpep_prediction.R`
   Loads each model's best fit, predicts ENPEP, and produces the AUC-vs-P(ENPEP=OB) scatter (**Fig 3j**).

5. **Figures** — `05_figures/`
   - `fig2_performance_comparison.R` → performance comparison (**Fig 2**).
   - `figS4_pca_batch.R` → PCA of the merged matrix coloured by dataset of origin (**Fig S4**);
     requires a `sample_dataset_map` file mapping each GSM sample to its dataset.
   - `figS5_heldout_violin.R` → violin of out-of-sample P(OB) for positive vs negative genes,
     parsed directly from the `OB_Model_result_frame_*.csv` tables (**Fig S5**).

---

## Notes

- **Held-out evaluation.** All reported metrics come from repeated random 80/20 gene-level hold-out
  resampling; each gene is scored only in cycles where it is held out of training, so figures such as
  Fig S5 reflect out-of-sample behaviour rather than training fit.
- **SVM.** `svm.R` is the training loop (`kernlab::ksvm`, RBF kernel `rbfdot`, `C = 1`), structurally
  identical to the other three original models; `svm_roc_analysis.R` operates on its output table.
- Absolute `setwd()` paths from the original working environment have been removed; set the working
  directory to the repository root before running.

---

## Citation

If you use this code, please cite the associated manuscript (details to be added on publication).

## License

Released under the MIT License (see `LICENSE`).
