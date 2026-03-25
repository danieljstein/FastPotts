# FastPotts

FastPotts is an R package for running loopy belief propagation in a Potts model, with a focus on spatial transcriptomics segmentation workflows.

## Installation

Install from GitHub using `devtools`:

```r
install.packages("devtools")
devtools::install_github("danieljstein/FastPotts")
```

## Basic usage

```r
library(FastPotts)

# Main high-level workflow
# result <- run_crf(transcripts_df, cell_signatures)
```

## Main functions

- `run_crf()` — run CRF-style inference for transcript labeling
- `hilbert_index()` — compute Hilbert ordering for spatial coordinates
- `build_potts_lbp_graph()` — convert sparse adjacency to LBP graph format
- `potts_lbp_parallel_cpp()` — parallel loopy belief propagation backend

## License

MIT (see LICENSE).
