# Transcript Graph Watershed Partitioning

This note describes the prototype algorithm implemented by
`partition_transcripts_watershed()`.

The goal is to partition transcripts into putative individual cells after a
cell-type posterior model, such as `spatial_basis_segmentation()`, has inferred
posterior probabilities for each transcript.

## Observed Data

For each transcript `i`, we assume:

- spatial location `x_i`
- posterior cell-type probabilities `q_i`
- optionally, an existing prior segmentation label `cell_id_i`

The prior segmentation may contain a value such as `"UNASSIGNED"` for
transcripts that were not assigned to a cell by the original segmentation.

The posterior vector is:

```text
q_i[k] = P(cell type = k | transcript i)
```

The algorithm treats these posteriors as soft marks on the spatial point cloud.
It does not require a hard cell-type label for each transcript.

## Neighborhood Graph

The first step builds a spatial transcript graph.

Each transcript is connected to its `n_neighbors` nearest spatial neighbors,
optionally discarding edges longer than `max_distance`.

The resulting undirected edge list contains:

```text
i, j, distance_ij
```

For posterior-aware diagnostics, each edge also receives a Jensen-Shannon
divergence:

```text
JS(q_i, q_j)
```

This gives a local graph on which density ascent and boundary diagnostics can
be computed without rasterizing the full spatial domain.

## Local Density

Each transcript receives a scalar density value `rho_i`.

The prototype supports two density modes.

### Total Density

In total-density mode, every neighboring transcript contributes equally:

```text
rho_i = 1 + sum_j K(distance_ij / h)
```

where `h` is `density_bandwidth` and the current kernel is Gaussian:

```text
K(r) = exp(-0.5 * r^2)
```

### Type-Weighted Density

In type-weighted mode, density is estimated separately for each posterior type:

```text
rho_i[k] = q_i[k] + sum_j q_j[k] K(distance_ij / h)
```

The scalar density used for graph ascent is then:

```text
rho_i = sum_k q_i[k] rho_i[k]
```

This favors density produced by nearby transcripts with compatible posterior
cell-type probabilities, while still preserving uncertainty.

## Graph Density Ascent

Each transcript chooses a parent among its higher-density graph neighbors.

For a directed candidate step from transcript `i` to neighbor `j`, the prototype
uses:

```text
score(i -> j)
  = rho_j - rho_i
    - distance_weight * distance_ij
    - posterior_weight * JS(q_i, q_j)
```

Only neighbors with `rho_j > rho_i` are eligible. If no higher-density neighbor
exists, transcript `i` is a local density mode.

Following parent pointers from every transcript produces a partition into
initial watershed basins:

```text
transcript -> ... -> density mode
```

Each basin is a candidate cell before merge correction.

## Basin Saddles

Adjacent basins are detected from transcript graph edges whose endpoints belong
to different basins.

For two adjacent basins `A` and `B`, the saddle density is:

```text
saddle(A, B)
  = max over boundary edges i-j of min(rho_i, rho_j)
```

The saddle ratio is:

```text
saddle_ratio(A, B)
  = saddle(A, B) / min(mode_density_A, mode_density_B)
```

High saddle ratios mean two basins meet through a relatively high-density
connection and may be over-split pieces of the same cell. Low saddle ratios
mean the basins are separated by a stronger density trough.

The basin adjacency table also records:

- number of boundary edges
- mean and minimum boundary-edge distance
- mean edge-level posterior divergence
- posterior divergence between basin-average posterior profiles

## Basin Merging

The first prototype applies a simple agglomerative merge rule.

Adjacent basins are considered in decreasing saddle-ratio order. Two basins are
merged when:

```text
saddle_ratio >= saddle_ratio_threshold
and basin_JS <= posterior_js_threshold
```

Small basins are also merged into compatible neighbors when possible:

```text
basin_size < min_transcripts_per_cell
and basin_JS <= posterior_js_threshold
```

This merge rule is intentionally simple. It is meant to expose interpretable
diagnostics first, not to be the final stopping criterion.

Future versions can replace this rule with an objective that combines:

- density prominence
- posterior compatibility
- expected transcript count per cell
- prior cell-count calibration
- shape or compactness penalties

## Prior Cell-Count Calibration

The function `estimate_prior_cell_type_counts()` uses an existing `cell_id`
column as a weak calibration signal.

The key idea is to estimate the total number of cells of each type by comparing
posterior-weighted transcript mass in all transcripts to posterior-weighted
transcript mass in reliable prior cells.

Prior cells with too few transcripts are filtered out:

```text
n_c >= min_transcripts_per_prior_cell
```

For a reliable prior cell `c`, the posterior-weighted type fraction is:

```text
f_c[k] = sum_{i in c} q_i[k] / sum_l sum_{i in c} q_i[l]
```

The assigned number of cells of type `k` is:

```text
N_assigned[k] = sum_c f_c[k]
```

The posterior transcript mass assigned to reliable prior cells is:

```text
M_assigned[k] = sum_{i in reliable prior cells} q_i[k]
```

The total posterior transcript mass is:

```text
M_total[k] = sum_i q_i[k]
```

The extrapolated total cell count is:

```text
N_hat[k] = N_assigned[k] * M_total[k] / M_assigned[k]
```

If `M_assigned[k] = 0`, the count is not estimable from the reliable prior
cells. In that case, the implementation returns `NA` for `N_hat[k]` and marks
the type with `estimable = FALSE`.

This estimate should be interpreted as a calibration prior, not a hard truth.
It assumes reliable assigned cells are approximately representative of all
transcripts of the same type. That assumption may fail when unassigned
transcripts are enriched for boundaries, low-quality regions, rare cell types,
or segmentation failures.

## Returned Diagnostics

The prototype returns:

- `transcripts_df`: input transcript table with `watershed_basin`,
  `watershed_cell`, and `watershed_density`
- `cell_assignment`: final merged cell assignment per transcript
- `initial_basin`: pre-merge density-ascent basin per transcript
- `parent`: density-ascent parent pointer per transcript
- `density`: scalar density per transcript
- `edges`: transcript graph edges with posterior divergence
- `basin_adjacency`: boundary/saddle diagnostics between initial basins
- `merge_history`: accepted basin merges
- `cells`: summary of final cells
- `initial_basins`: summary of initial basins
- `prior_cell_type_counts`: optional prior-cell calibration table
- `parameters`: effective algorithm parameters

## Suggested Initial Workflow

Run `spatial_basis_segmentation()` first:

```r
fit <- spatial_basis_segmentation(
    transcripts_df,
    cell_signatures,
    basis = "3d",
    s = 2
)
```

Then run the watershed prototype:

```r
partition <- partition_transcripts_watershed(
    transcripts_df = fit$transcripts_df,
    posterior = fit$marginals,
    n_neighbors = 20,
    density_mode = "type_weighted"
)
```

Inspect:

```r
partition$cells
partition$basin_adjacency
partition$merge_history
partition$prior_cell_type_counts
```

The most important early diagnostic is whether the initial basins are
over-splitting cells, under-splitting adjacent cells, or both. The saddle table
is intended to make those errors visible enough to tune or replace the merge
criterion.
