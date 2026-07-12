# Data

Place the following input files in this folder before running the scripts.

## `log2_dataset_frame.csv` 
Merged expression matrix used as the feature space.
- Rows: genes (first column `gene`/`Symbol`), Columns: 223 biological samples (GSM/accession IDs).
- Values: log2-transformed expression; no batch correction applied.

## `Label_genes.csv`  (include this)
Gene-level labels used for training.
- Columns: `Symbol`, `Label`.
- `Label = 1` → positive (osteoblast differentiation); `Label = 0` → negative (adipogenesis).
- Labels are derived from Gene Ontology (161 genes total).

## `sample_dataset_map.xlsx`  (only for figS4_pca_batch.R)
Maps each sample to its dataset of origin, for colouring the PCA.
- Columns: `Sample` (GSM ID, matching the matrix column names), `Dataset` (e.g. GSE1367).
- Build from the manuscript Supplementary Table.
