# Continuous Spatial Basis Segmentation Model

This note describes the inference model implemented by
`spatial_basis_segmentation()`.

## Observed Data

For each transcript `i`, we observe:

- location `x_i`
- gene identity `g_i`

The model assumes a fixed cell-type signature matrix

```text
theta[g, k] = P(gene = g | cell type = k)
```

with `G` genes and `K` cell types.

In mathematical notation:

$$
\theta_{gk} = P(g_i = g \mid z_i = k),
$$

where $z_i \in \{1,\dots,K\}$ is the latent cell type of transcript $i$.

## Spatial Basis

The spatial domain is represented by a lattice basis:

- `basis = "2d"`: 2D triangular lattice
- `basis = "3d"`: 3D body-centered cubic lattice

The older values `basis = "tri"` and `basis = "bcc"` are accepted as aliases,
but the public API documents the dimensional names because they are clearer for
segmentation users.

The user-facing mesh size `s` is defined as the long Delaunay-cell edge scale:

- in 2D, triangular Delaunay cells are equilateral triangles with side length
  `s`
- in 3D, BCC Delaunay tetrahedra have two opposite long edges of length `s`
  and four shorter edges of length `sqrt(3) * s / 2`

With this convention, increasing or decreasing `s` has a comparable
interpretation in the 2D and 3D APIs.

Each transcript location has a small number of active basis functions:

- 3 active basis functions for the 2D triangular lattice
- 4 active basis functions for the 3D BCC lattice

Let

```text
phi_im
```

be the barycentric weight of basis point `m` at transcript location `x_i`.
Most `phi_im` are zero.

We write this as:

$$
\phi_{im} = \phi_m(x_i),
$$

with $\phi_{im} = 0$ for all but the locally active basis points.

The implemented code stores these sparse interpolation values through:

- `basis_id[i, a]`: the active basis point index
- `basis_weight[i, a]`: the corresponding barycentric weight

where `a = 1, ..., 3` for the triangular lattice and `a = 1, ..., 4` for the
BCC lattice.

## Spatial Latent Field

Each basis point `m` has a latent logit value for each cell type:

```text
w[k, m]
```

The spatial logit for transcript `i` and cell type `k` is interpolated from
nearby basis points:

```text
f[i, k] = sum_m phi_im * w[k, m]
```

Equivalently:

$$
f_{ik}
= f_k(x_i)
= \sum_m \phi_{im} w_{km}.
$$

Equivalently, using only the active basis functions:

```text
f[i, k] = sum_a basis_weight[i, a] * w[k, basis_id[i, a]]
```

That sparse form is:

$$
f_{ik}
= \sum_{a \in \mathcal A(i)}
\phi_{i,a} \, w_{k,m_{i,a}},
$$

where $\mathcal A(i)$ is the active basis set for transcript $i$.

The spatial prior over cell types is:

```text
p[i, k] = P(z_i = k | x_i, W)
        = softmax_k(f[i, 1], ..., f[i, K])
```

where `z_i` is the latent cell type of transcript `i`.

In equation form:

$$
p_{ik}
= P(z_i = k \mid x_i, W)
= \frac{\exp(f_{ik})}{\sum_{\ell=1}^K \exp(f_{i\ell})}.
$$

## Identifiability

The softmax is invariant to adding the same constant to every cell-type logit
at a location. To remove this redundancy, the implementation uses a reference
class parameterization:

```text
w[K, m] = 0
```

for every basis point `m`. Only the first `K - 1` cell-type fields are
optimized.

Thus:

$$
w_{Km} = 0
\quad \text{for all basis points } m.
$$

## Gene Emission Model

Given the latent cell type, the transcript gene is sampled from the fixed
signature:

```text
P(g_i = g | z_i = k) = theta[g, k]
```

Equivalently:

$$
P(g_i = g \mid z_i = k) = \theta_{gk}.
$$

The marginal likelihood for transcript `i` is:

```text
P(g_i | x_i, W, theta)
  = sum_k theta[g_i, k] * p[i, k]
```

In math:

$$
P(g_i \mid x_i, W, \Theta)
= \sum_{k=1}^K \theta_{g_i k} p_{ik}.
$$

Using logits, this is evaluated as:

```text
log P(g_i | x_i, W, theta)
  = logsumexp_k(log theta[g_i, k] + f[i, k])
    - logsumexp_k(f[i, k])
```

Equivalently:

$$
\log P(g_i \mid x_i, W, \Theta)
=
\operatorname{logsumexp}_{k}
\left(
\log \theta_{g_i k} + f_{ik}
\right)
-
\operatorname{logsumexp}_{k}
\left(
f_{ik}
\right).
$$

The second term appears because `p[i, k]` is the softmax-normalized spatial
prior.

## Transcript Posterior Probabilities

After fitting the spatial field, transcript posterior probabilities are:

```text
q[i, k] = P(z_i = k | g_i, x_i, W, theta)
        = theta[g_i, k] * p[i, k]
          / sum_l theta[g_i, l] * p[i, l]
```

In equation form:

$$
q_{ik}
= P(z_i = k \mid g_i, x_i, W, \Theta)
=
\frac{
\theta_{g_i k} p_{ik}
}{
\sum_{\ell=1}^K \theta_{g_i \ell} p_{i\ell}
}.
$$

Equivalently:

```text
q[i, k] = softmax_k(log theta[g_i, k] + f[i, k])
```

That is:

$$
q_{ik}
=
\frac{
\exp(\log \theta_{g_i k} + f_{ik})
}{
\sum_{\ell=1}^K
\exp(\log \theta_{g_i \ell} + f_{i\ell})
}.
$$

These posterior probabilities are returned as `marginals`.

## Optional Signature Refinement

The reference signatures can optionally be updated from the current posterior
assignments. This is disabled by default because the reference signatures are
often the main anchor preventing neighboring or abundant cell types from
absorbing rare ones.

Let $c_{gk}$ be the posterior-weighted transcript count for gene $g$ and cell
type $k$:

$$
c_{gk} = \sum_{i:g_i=g} q_{ik}.
$$

For each cell type, the update uses a Dirichlet posterior mean anchored to the
input reference signature $\theta^{(0)}_{\cdot k}$:

$$
\tilde{\theta}_{gk}
=
\frac{c_{gk} + \alpha \theta^{(0)}_{gk}}
{\sum_{g'} c_{g'k} + \alpha}.
$$

The prior strength $\alpha$ is controlled by `signature_prior_strength`. If it
is `NULL`, the default is $\alpha = G$, the number of genes, so the reference
contributes approximately one pseudo-transcript per gene. This is conservative:
large cell-type-specific transcript counts can move the signature, but weak or
rare assignments remain close to the reference.

The update is damped before refitting the spatial field:

$$
\theta^{(t+1)}_{\cdot k}
=
(1-\rho)\theta^{(t)}_{\cdot k}
+
\rho\tilde{\theta}_{\cdot k},
$$

where $\rho$ is `signature_update_rate`. The default is $\rho = 0.25$, meaning
each update moves one quarter of the way toward the current posterior estimate.
After the update, columns are floored and normalized.

The alternating algorithm is:

1. Fit the spatial field with the current signatures.
2. Compute posterior transcript assignments.
3. Update signatures from posterior-weighted gene counts.
4. Warm-start another spatial field fit from the previous basis weights.

This repeats `signature_update_iters` times when `refine_signatures = TRUE`.
The returned object includes `cell_signatures_initial`, final
`cell_signatures`, `signature_history`, `signature_update_history`, and
`optim_history`.

## MAP Objective

The implemented estimator minimizes the negative log posterior:

```text
L(W) = mean_i[-log P(g_i | x_i, W, theta)] + R(W)
```

where the first term is the mean negative transcript log likelihood over the
$N$ retained transcripts.

More explicitly:

$$
\mathcal L(W)
=
\frac{1}{N}
\sum_i
\left[
-\log
\left[
\sum_{k=1}^K \theta_{g_i k} p_{ik}
\right]
\right]
+
R(W).
$$

## Spatial Smoothing

The regularizer smooths neighboring lattice basis points:

```text
R(W) = lambda * mean_edges sum_k rho((w[k, m] - w[k, n]) / d[m, n])
```

where `E` is the lattice neighbor graph, $d_{mn}$ is the physical distance
between neighboring basis points, and `lambda >= 0` controls smoothing. The
implementation supports three choices for `rho`.

This convention makes `lambda` a per-transcript, per-edge spatial-slope
penalty. For the BCC basis, axis edges have distance $s$ and center-to-corner
edges have distance $\sqrt{3}s/2$; for the triangular basis, all neighbor edges
have distance $s$.

### Quadratic

The default is quadratic smoothing:

$$
R(W)
=
\frac{\lambda}{|E|}
\sum_{(m,n)\in E}
\sum_{k=1}^{K-1}
\frac{1}{2}
\left(
\frac{w_{km} - w_{kn}}{d_{mn}}
\right)^2.
$$

This is diffusive and strongly penalizes large spatial slopes.

### Huber

Huber smoothing is quadratic for small spatial logit slopes and linear for
larger slopes:

$$
\rho_\delta(r)
=
\begin{cases}
\frac{1}{2}r^2, & |r| \le \delta, \\
\delta \left(|r| - \frac{1}{2}\delta\right), & |r| > \delta.
\end{cases}
$$

The transition scale `delta` is measured in logit units per spatial unit. The
default `delta = 1` treats slopes below roughly one logit unit per spatial unit
as smooth variation and larger slopes more like boundaries.

### Bounded

The bounded option is a Potts-like smooth approximation:

$$
\rho_\sigma(r)
=
\sigma^2
\left[
1 - \exp\left(-\frac{r^2}{2\sigma^2}\right)
\right].
$$

The saturation scale `sigma` is measured in logit units per spatial unit. The
small-slope behavior is $\rho_\sigma(r) \approx r^2 / 2$, matching quadratic
smoothing with the same `lambda`. The asymptotic per-component penalty is
$\sigma^2$.

The sum is over $K-1$ optimized classes because the final class is the
reference class with $w_{Km}=0$. These penalties are currently applied
component-wise to each optimized cell-type logit.

For the 2D triangular lattice, neighbors are the six adjacent triangular
lattice points, represented by three undirected offset directions:

```text
(1, 0), (0, 1), (1, -1)
```

For the 3D BCC lattice, neighbors include:

- same-parity axis neighbors at offsets `(+-2, 0, 0)`, `(0, +-2, 0)`,
  `(0, 0, +-2)`
- opposite-parity body-diagonal neighbors at offsets `(+-1, +-1, +-1)`

The implementation only includes neighbor edges where both endpoint basis
points are present in the local basis set touched by the transcripts.

## Basis-Point Purity

The spatial smoothing penalty controls how much neighboring basis points agree.
A separate purity penalty can encourage each individual basis point to place
most of its probability mass on one or a few cell types.

Let

$$
\pi_{mk}
=
\frac{\exp(w_{km})}{\sum_{\ell=1}^K \exp(w_{\ell m})}
$$

be the cell-type prior at basis point $m$, with the reference-class logit
$w_{Km}=0$.

The implementation supports:

- `purity = "none"`: no basis-point purity penalty
- `purity = "entropy"`: penalize high entropy
- `purity = "gini"`: penalize Gini impurity

The entropy penalty is:

$$
R_{\text{entropy}}(W)
=
\frac{\alpha}{M}
\sum_m
\left[
-
\sum_{k=1}^K
\pi_{mk}\log \pi_{mk}
\right].
$$

The Gini impurity penalty is:

$$
R_{\text{gini}}(W)
=
\frac{\alpha}{M}
\sum_m
\left[
1 -
\sum_{k=1}^K
\pi_{mk}^2
\right].
$$

Here $\alpha$ is `purity_lambda` and $M$ is the number of active basis points.
Both penalties are minimized when each basis point is close to a one-hot
cell-type prior. Gini is bounded and often a gentler first choice; entropy is
sharper near the simplex corners.

## Gradient Derivation

Let

```text
p[i, k] = softmax_k(f[i, .])
q[i, k] = softmax_k(log theta[g_i, k] + f[i, .])
```

For one transcript, the negative log likelihood is:

```text
ell_i = -logsumexp_k(log theta[g_i, k] + f[i, k])
        + logsumexp_k(f[i, k])
```

Equivalently:

$$
\ell_i
=
-
\operatorname{logsumexp}_{k}
\left(
\log \theta_{g_i k} + f_{ik}
\right)
+
\operatorname{logsumexp}_{k}
\left(
f_{ik}
\right).
$$

The derivative with respect to the spatial logit is:

```text
d ell_i / d f[i, k] = p[i, k] - q[i, k]
```

That is:

$$
\frac{\partial \ell_i}{\partial f_{ik}}
=
p_{ik} - q_{ik}.
$$

By the chain rule,

```text
d L / d w[k, m]
  = (1 / N) * sum_i phi_im * (p[i, k] - q[i, k])
    + d R / d w[k, m]
```

So:

$$
\frac{\partial \mathcal L}{\partial w_{km}}
=
\frac{1}{N}
\sum_i
\phi_{im}
\left(
p_{ik} - q_{ik}
\right)
+
\frac{\partial R}{\partial w_{km}}.
$$

Because each transcript touches only 3 or 4 basis functions, the likelihood
gradient is sparse.

For the quadratic smoothing term:

```text
d R / d w[k, m]
  = (lambda / |E|) * sum_(n adjacent to m)
      (w[k, m] - w[k, n]) / d[m, n]^2
```

In equation form:

$$
\frac{\partial R}{\partial w_{km}}
=
\frac{\lambda}{|E|}
\sum_{n : (m,n)\in E}
\frac{
w_{km} - w_{kn}
}{
d_{mn}^2
}.
$$

The C++ backend computes the objective and analytic gradient together.

For Huber smoothing, the derivative of the scalar penalty is:

$$
\rho_\delta'(r)
=
\begin{cases}
r, & |r| \le \delta, \\
\delta \operatorname{sign}(r), & |r| > \delta.
\end{cases}
$$

For bounded smoothing, the derivative is:

$$
\rho_\sigma'(r)
=
\exp\left(-\frac{r^2}{2\sigma^2}\right)
r.
$$

For a general purity penalty $h(\pi_m)$, the softmax chain rule gives:

$$
\frac{\partial h}{\partial w_{km}}
=
\pi_{mk}
\left[
\frac{\partial h}{\partial \pi_{mk}}
-
\sum_{\ell=1}^K
\pi_{m\ell}
\frac{\partial h}{\partial \pi_{m\ell}}
\right].
$$

For entropy,

$$
\frac{\partial h}{\partial \pi_{mk}}
=
-(\log \pi_{mk} + 1).
$$

For Gini impurity,

$$
\frac{\partial h}{\partial \pi_{mk}}
=
-2\pi_{mk}.
$$

## Optimization

The implementation uses `stats::optim(method = "L-BFGS-B")`.

L-BFGS-B uses objective values and gradients. It does not compute the full
Hessian. Instead, it builds a limited-memory quasi-Newton approximation to
curvature from recent parameter and gradient changes.

This is a good default for the current model because:

- the objective is smooth
- gradients are deterministic full-batch gradients
- the analytic gradient is available
- full Hessians would be too large for realistic numbers of basis points

## Returned Quantities

`spatial_basis_segmentation()` returns:

- `marginals`: posterior transcript probabilities `q[i, k]`
- `spatial_prior`: fitted spatial priors `p[i, k]`
- `logits`: interpolated transcript logits `f[i, k]`
- `basis_weights`: fitted basis coefficients `w[k, m]`
- `basis_points`: spatial coordinates of lattice basis points
- `basis_lattice`: integer lattice coordinates
- `basis_edges`: neighboring basis point graph with physical edge distances
- `transcripts_df`: filtered transcript data with MAP labels
- `optim`: the optimizer result
- `optim_history`: optimizer results from each spatial field fit
- `cell_signatures_initial`: input signatures after filtering and normalization
- `cell_signatures`: final signatures used for the returned posterior
- `signature_history`: signatures after each refinement step
- `signature_update_history`: effective counts and max change per update

## Current Limitations

The current implementation is intentionally a first MAP estimator:

- Signature refinement is optional and conservative; poor initial signatures
  can still bias the posterior updates.
- Quadratic regularization encourages smooth fields and may blur sharp cell
  boundaries. Huber and bounded regularization are available to reduce this,
  but introduce additional logit-scale parameters.
- Only basis points touched by at least one transcript are included.
- The final cell type is used as the reference class.

Natural next extensions include:

- robust or total-variation-like smoothing to preserve sharper boundaries
- minibatch Adam for very large datasets, followed by L-BFGS-B polishing
- explicit background/noise components
- priors or penalties that encourage sparse cell-type occupancy per basis point
- GPU acceleration

Other open questions and directions:

- Can we go beyond discrete cell type signatures to a more continuous view of cell state?
  - Either with probabilistic cell embeddings or logistic normal models of expression within each cell type?
  - If learned from the data, how do we avoid learning segmentation contamination initially and getting stuck without improving cell type purity?
  - Could these methods plug into expanding cell atlas foundation models?
- What is the best way to reliably / automatically initialize the cell type signatures?
  - Can you support finer cell states, or will that lead to over-splitting of individual cells due to ambiguity? Would appropriate spatial regularization prevent this?
  - How can we deal with noise / spatial contamination in the initial cell type signatures? We want to avoid getting stuck at a signature (even after refinement) that includes this contamination
  - Subcellular compartmentalization: Are nuclear transcripts leading to over-segmentation of some cells? Would it be worth defining a separate nuclear signature that can be added to other cell type signatures?
- What is the best way to choose the mesh size and the regularization? Is there sufficient sharing of information across basis points to allow for small mesh sizes?
  - Why were the myeloid cells being called B/plasma/pDC in some runs with small 3D mesh size? Similarly, small fragments of different kinds of epithelial cells when lambda = 0? Would different cell type signatures help?
  - How are we controlling the rate of change / diffusion allowed at boundaries?
- How do we go from cell types to individual cells?
  - One proposal: HDBSCAN-style partitioning of transcripts based on density / transcript NN-graph, only allowing contacts to the same cell type (how would this be done with continuous embeddings rather than discrete signatures?)
  - Could try to cut so that the cell counts match the number when running Cellpose with low filtering and some post-QC
- Can we build similar multiscale models?
  - From subcellular up to full tissues

Datasets of interest:

- 3D spatial transcriptomics segmentation
