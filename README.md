# FastPotts

FastPotts is an R package for spatial transcriptomics segmentation. It includes
the original Potts-model conditional random field workflow over transcript
graphs, plus a continuous spatial basis segmentation model for inferring
cell-type fields from transcript locations and gene signatures.

Version 0.1 adds the spatial basis segmentation workflow, including 2D
triangular and 3D BCC interpolation bases, mesh-normalized spatial
regularization, basis-point purity penalties, and optional conservative
signature refinement.

## Installation

Install from GitHub using `devtools`:

```r
install.packages("devtools")
devtools::install_github("danieljstein/FastPotts")
```

## Basic usage

```r
library(FastPotts)

# Transcript-graph CRF workflow
crf_fit <- run_crf(transcripts_df, cell_signatures)

# Continuous spatial basis workflow
basis_fit <- spatial_basis_segmentation(
  transcripts_df,
  cell_signatures,
  s = 2,
  lambda = 1,
  regularization = "bounded"
)
```

## Main functions

- `run_crf()` — run CRF-style inference on a transcript-neighbor graph
- `spatial_basis_segmentation()` — infer a continuous spatial cell-type field
  using 2D or 3D barycentric basis functions
- `tri_barycentric()` — compute 2D triangular-lattice barycentric weights
- `bcc_barycentric()` — compute 3D BCC-lattice barycentric weights
- `hilbert_index()` — compute Hilbert ordering for spatial coordinates
- `build_potts_lbp_graph()` — convert sparse adjacency to LBP graph format
- `potts_lbp_parallel_cpp()` — parallel loopy belief propagation backend

## Spatial basis segmentation

`spatial_basis_segmentation()` models each transcript's cell type as a spatial
softmax field interpolated from nearby lattice basis points. The observed gene
identity is scored against the supplied cell-type signatures, and the returned
`marginals` matrix gives posterior transcript probabilities over cell types.

By default, `spatial_basis_segmentation()` uses `basis = "3d"` for a BCC
lattice over `x`/`y`/`z` and `regularization = "bounded"` for boundary-aware
smoothing. Use `basis = "2d"` for a triangular lattice over `x`/`y`. The mesh
size `s` controls the Delaunay cell scale. The spatial regularization parameter
`lambda` is normalized to an average transcript likelihood and an average
basis-edge slope penalty, making it more comparable across mesh sizes and
transcript counts.

Available regularizers are:

- `regularization = "quadratic"` for smooth fields
- `regularization = "huber"` for boundary-preserving linear tails
- `regularization = "bounded"` for Potts-like bounded edge penalties

Basis-point purity can be encouraged with `purity = "entropy"` or
`purity = "gini"`, and signatures can be conservatively refined with
`refine_signatures = TRUE`.

## License

GPL-3 (see `LICENSE.md`).
