///|
/// splatter - single-cell RNA-seq count data simulation.
/// Inspired by the Bioconductor splatter package (Zappia et al., 2017,
/// Genome Biology). Implements the splat model: gamma-distributed gene
/// base means, log-normal library size factors, trended BCV dispersion,
/// negative-binomial counts, logistic dropout, group-wise differential
/// expression, batch effects and pseudotime paths. All randomness comes
/// from a deterministic LCG stream so simulations are reproducible.

///|
/// Simulation parameters. Defaults mirror splatter's splatEstimate /
/// splatSimulate defaults.
pub struct SplatParams {
  /// Number of genes to simulate.
  n_genes : Int
  /// Gamma shape for base gene means (mean.shape).
  mean_shape : Double
  /// Gamma rate for base gene means (mean.rate).
  mean_rate : Double
  /// Log-normal location for library sizes (lib.loc).
  lib_loc : Double
  /// Log-normal scale for library sizes (lib.scale).
  lib_scale : Double
  /// Common biological coefficient of variation (bcv.common).
  bcv_common : Double
  /// Cells per group; an empty array means one group with all cells.
  group_cells : Array[Int]
  /// Fraction of genes that are differentially expressed per group (de.prob).
  de_prob : Double
  /// Fraction of DE genes down-regulated (de.downProb).
  de_down_prob : Double
  /// Location of log-normal DE factor (de.facLoc).
  de_fac_loc : Double
  /// Scale of log-normal DE factor (de.facScale).
  de_fac_scale : Double
  /// Cells per batch; an empty array means a single batch.
  batch_cells : Array[Int]
  /// Location of log-normal batch factor (batch.facLoc).
  batch_fac_loc : Double
  /// Scale of log-normal batch factor (batch.facScale).
  batch_fac_scale : Double
  /// Dropout mid point on log10(mean) (dropout.mid); used when
  /// enable_dropout is true.
  dropout_mid : Double
  /// Dropout logistic slope on log10(mean) (dropout.shape).
  dropout_shape : Double
  /// Whether to apply dropout.
  enable_dropout : Bool
  /// Simulate the first group as a pseudotime path.
  enable_path : Bool
  /// Path length (number of cells along the path when enable_path).
  n_path_steps : Int
  /// Scale of log-normal end-state shift for path genes (path.skew style).
  path_shift_scale : Double
  /// Fraction of genes that change along the path.
  path_probability : Double
}

///|
/// Default parameters: 100 genes, 100 cells, splatter defaults otherwise.
pub fn SplatParams::default() -> SplatParams {
  SplatParams::{
    n_genes: 100,
    mean_shape: 0.6,
    mean_rate: 0.3,
    lib_loc: 11.0,
    lib_scale: 0.2,
    bcv_common: 0.1,
    group_cells: [100],
    de_prob: 0.1,
    de_down_prob: 0.5,
    de_fac_loc: 0.1,
    de_fac_scale: 0.4,
    batch_cells: [],
    batch_fac_loc: 0.1,
    batch_fac_scale: 0.5,
    dropout_mid: 2.0,
    dropout_shape: -1.0,
    enable_dropout: false,
    enable_path: false,
    n_path_steps: 10,
    path_shift_scale: 0.1,
    path_probability: 0.5,
  }
}

///|
/// Build parameters with an explicit gene/cell count; optional arguments
/// override the splatter defaults.
pub fn SplatParams::create(
  n_genes : Int,
  n_cells : Int,
  mean_shape? : Double = 0.6,
  mean_rate? : Double = 0.3,
  lib_loc? : Double = 11.0,
  lib_scale? : Double = 0.2,
  bcv_common? : Double = 0.1,
  group_cells? : Array[Int] = [],
  de_prob? : Double = 0.1,
  de_down_prob? : Double = 0.5,
  de_fac_loc? : Double = 0.1,
  de_fac_scale? : Double = 0.4,
  batch_cells? : Array[Int] = [],
  batch_fac_loc? : Double = 0.1,
  batch_fac_scale? : Double = 0.5,
  dropout_mid? : Double = 2.0,
  dropout_shape? : Double = -1.0,
  enable_dropout? : Bool = false,
  enable_path? : Bool = false,
  n_path_steps? : Int = 10,
  path_shift_scale? : Double = 0.1,
  path_probability? : Double = 0.5,
) -> SplatParams {
  let effective_groups = if group_cells.length() > 0 {
    group_cells
  } else {
    [n_cells]
  }
  SplatParams::{
    n_genes,
    mean_shape,
    mean_rate,
    lib_loc,
    lib_scale,
    bcv_common,
    group_cells: effective_groups,
    de_prob,
    de_down_prob,
    de_fac_loc,
    de_fac_scale,
    batch_cells,
    batch_fac_loc,
    batch_fac_scale,
    dropout_mid,
    dropout_shape,
    enable_dropout,
    enable_path,
    n_path_steps,
    path_shift_scale,
    path_probability,
  }
}

///|
/// Simulated dataset.
pub struct SplatResult {
  /// Counts matrix, rows = genes, columns = cells.
  counts : Array[Array[Int]]
  /// Base gene means (gamma draws, before DE/batch effects).
  gene_base_means : Array[Double]
  /// Library size factor per cell (multiplicative).
  cell_lib_factors : Array[Double]
  /// Group index of each cell.
  cell_group : Array[Int]
  /// Batch index of each cell.
  cell_batch : Array[Int]
  /// Pseudotime progress of each cell in [0, 1] (0 for non-path cells).
  cell_pseudotime : Array[Double]
  /// Number of (gene, cell) dropout zero events.
  n_dropout_events : Int
  /// Final per-group gene log-means: log_means[group][gene].
  group_log_means : Array[Array[Double]]
}

///|
/// Deterministic LCG (Park-Miller style, as used elsewhere in the
/// codebase). Returns (next_state, uniform in (0, 1)).
pub fn splat_lcg(state : Int) -> (Int, Double) {
  let normalized = (state - 1) % 2147483646 + 1
  let quotient = normalized / 44488
  let remainder = normalized % 44488
  let candidate = 48271 * remainder - 3399 * quotient
  let next = if candidate > 0 { candidate } else { candidate + 2147483647 }
  (next, next.to_double() / 2147483647.0)
}

///|
/// Standard normal via Box-Muller.
pub fn splat_normal(state : Int) -> (Int, Double) {
  let (state1, raw1) = splat_lcg(state)
  let (state2, raw2) = splat_lcg(state1)
  let first = raw1.max(1.0e-15)
  let value = (-2.0 * @math.ln(first)).sqrt() *
    @math.cos(6.283185307179586 * raw2)
  (state2, value)
}

///|
/// Gamma sample (Marsaglia-Tsang), shape > 0.
pub fn splat_gamma(shape : Double, state : Int) -> (Int, Double) {
  if shape < 1.0 {
    let (state1, uniform) = splat_lcg(state)
    let (state2, sampled) = splat_gamma(shape + 1.0, state1)
    return (state2, sampled * @math.exp(@math.ln(uniform.max(1.0e-15)) / shape))
  }
  let d = shape - 1.0 / 3.0
  let c = 1.0 / (9.0 * d).sqrt()
  let mut current = state
  while true {
    let (state1, normal) = splat_normal(current)
    current = state1
    let base = 1.0 + c * normal
    if base > 0.0 {
      let volume = base * base * base
      let (state2, uniform) = splat_lcg(current)
      current = state2
      if uniform < 1.0 - 0.0331 * normal * normal * normal * normal {
        return (current, d * volume)
      }
      if @math.ln(uniform.max(1.0e-15)) <
        0.5 * normal * normal + d * (1.0 - volume + @math.ln(volume)) {
        return (current, d * volume)
      }
    }
  }
  (current, shape)
}

///|
/// Poisson sample. Knuth's method for small lambda, transformed-rejection
/// style normal approximation for large lambda.
pub fn splat_poisson(lambda : Double, state : Int) -> (Int, Int) {
  if lambda <= 0.0 {
    return (state, 0)
  }
  if lambda > 30.0 {
    let (state1, normal) = splat_normal(state)
    let value = (lambda + normal * lambda.sqrt()).floor().max(0.0)
    return (state1, value.to_int())
  }
  let l = @math.exp(-lambda)
  let mut k = 0
  let mut p = 1.0
  let mut current = state
  while true {
    let (state1, uniform) = splat_lcg(current)
    current = state1
    p = p * uniform
    if p <= l {
      return (current, k)
    }
    k = k + 1
  }
  (current, 0)
}

///|
/// Log-normal sample: exp(normal(loc, scale)).
fn splat_lognormal(loc : Double, scale : Double, state : Int) -> (Int, Double) {
  let (state1, normal) = splat_normal(state)
  (state1, @math.exp(loc + scale * normal))
}

///|
/// Logistic dropout probability for a log10 mean.
fn splat_dropout_probability(
  log10_mean : Double,
  mid : Double,
  shape : Double,
) -> Double {
  1.0 / (1.0 + @math.exp(shape * (log10_mean - mid)))
}

///|
/// Draw DE / batch multiplicative factors: `n` selections, each factor is
/// a log-normal(loc, scale) applied as an up or down shift.
fn splat_draw_factors(
  n : Int,
  down_probability : Double,
  loc : Double,
  scale : Double,
  state : Int,
) -> (Int, Array[Double]) {
  let mut current = state
  let factors : Array[Double] = Array::new(capacity=n)
  for _i in 0..<n {
    let (s1, direction) = splat_lcg(current)
    current = s1
    let (s2, magnitude) = splat_lognormal(loc, scale, current)
    current = s2
    factors.push(
      if direction < down_probability {
        1.0 / magnitude
      } else {
        magnitude
      },
    )
  }
  (current, factors)
}

///|
/// Select gene indices below a uniform probability threshold.
fn splat_select_genes(
  n_genes : Int,
  probability : Double,
  state : Int,
) -> (Int, Array[Int]) {
  let mut current = state
  let selected : Array[Int] = []
  for g in 0..<n_genes {
    let (s1, uniform) = splat_lcg(current)
    current = s1
    if uniform < probability {
      selected.push(g)
    }
  }
  (current, selected)
}

///|
/// Expand per-group / per-batch cell-count arrays into per-cell indices.
fn splat_expand_cells(counts_per_unit : Array[Int]) -> Array[Int] {
  let out : Array[Int] = []
  for unit in 0..<counts_per_unit.length() {
    for _c in 0..<counts_per_unit[unit] {
      out.push(unit)
    }
  }
  out
}

///|
/// Run the splat simulation.
pub fn splat_simulate(params : SplatParams, seed : Int) -> SplatResult {
  let mut state = if seed <= 0 { 1 } else { seed }
  let n_genes = params.n_genes
  // 1. base gene means
  let base_means : Array[Double] = Array::new(capacity=n_genes)
  for _g in 0..<n_genes {
    let (s1, mean) = splat_gamma(params.mean_shape, state)
    state = s1
    base_means.push(mean / params.mean_rate)
  }
  // 2. cell layout: groups, batches
  let group_cells = if params.group_cells.length() > 0 {
    params.group_cells
  } else {
    [1]
  }
  let cell_group = splat_expand_cells(group_cells)
  let n_cells = cell_group.length()
  let batch_cells = if params.batch_cells.length() > 0 {
    params.batch_cells
  } else {
    [n_cells]
  }
  // align batches to the total number of cells by trimming/extending
  let effective_batch_cells : Array[Int] = Array::new()
  let mut assigned = 0
  for c in batch_cells {
    if assigned >= n_cells {
      break
    }
    let take = if assigned + c > n_cells { n_cells - assigned } else { c }
    effective_batch_cells.push(take)
    assigned = assigned + take
  }
  if assigned < n_cells {
    effective_batch_cells.push(n_cells - assigned)
  }
  let cell_batch = splat_expand_cells(effective_batch_cells)
  // 3. library size factors
  let lib_sizes : Array[Double] = Array::new(capacity=n_cells)
  for _c in 0..<n_cells {
    let (s1, size) = splat_lognormal(params.lib_loc, params.lib_scale, state)
    state = s1
    lib_sizes.push(size)
  }
  let lib_geomean = {
    let mut total = 0.0
    for v in lib_sizes {
      total = total + @math.ln(v.max(1.0e-15))
    }
    @math.exp(total / n_cells.to_double())
  }
  let lib_factors : Array[Double] = Array::new(capacity=n_cells)
  for v in lib_sizes {
    lib_factors.push(v / lib_geomean)
  }
  // 4. group log-means: base means + optional path + DE shifts
  let n_groups = group_cells.length()
  let group_log_means : Array[Array[Double]] = Array::new(capacity=n_groups)
  for group in 0..<n_groups {
    let log_means : Array[Double] = Array::new(capacity=n_genes)
    let is_path_group = params.enable_path && group == 0
    let (state1, path_genes) = if is_path_group {
      splat_select_genes(n_genes, params.path_probability, state)
    } else {
      (state, [])
    }
    state = state1
    let (state2, path_shifts) = if is_path_group {
      splat_draw_factors(
        path_genes.length(),
        0.5,
        0.0,
        params.path_shift_scale,
        state,
      )
    } else {
      (state, [])
    }
    state = state2
    let (state3, de_genes) = splat_select_genes(n_genes, params.de_prob, state)
    state = state3
    let (state4, de_factors) = splat_draw_factors(
      de_genes.length(),
      params.de_down_prob,
      params.de_fac_loc,
      params.de_fac_scale,
      state,
    )
    state = state4
    for g in 0..<n_genes {
      let mut log_mean = @math.ln(base_means[g].max(1.0e-15))
      if is_path_group {
        let mut in_path = false
        let mut path_idx = -1
        for pi in 0..<path_genes.length() {
          if path_genes[pi] == g {
            in_path = true
            path_idx = pi
            break
          }
        }
        if in_path {
          log_mean = log_mean + path_shifts[path_idx]
        }
      }
      let mut de_idx = -1
      for di in 0..<de_genes.length() {
        if de_genes[di] == g {
          de_idx = di
          break
        }
      }
      if de_idx >= 0 {
        log_mean = log_mean + @math.ln(de_factors[de_idx])
      }
      log_means.push(log_mean)
    }
    group_log_means.push(log_means)
  }
  // 5. batch factors (per batch, applied to all genes)
  let n_batches = effective_batch_cells.length()
  let batch_factors : Array[Array[Double]] = Array::new(capacity=n_batches)
  for _b in 0..<n_batches {
    let factors : Array[Double] = Array::new(capacity=n_genes)
    for _g in 0..<n_genes {
      let (s1, factor) = splat_lognormal(
        params.batch_fac_loc,
        params.batch_fac_scale,
        state,
      )
      state = s1
      factors.push(factor)
    }
    batch_factors.push(factors)
  }
  // 6. counts: NB(mu = exp(log_mean) * batch_factor * lib_factor, dispersion)
  //    dispersion = bcv_common + 1/sqrt(mean)
  let counts : Array[Array[Int]] = Array::new(capacity=n_genes)
  for g in 0..<n_genes {
    let row : Array[Int] = Array::new(capacity=n_cells)
    for c in 0..<n_cells {
      let group = cell_group[c]
      let batch = cell_batch[c]
      let mut log_mean = group_log_means[group][g]
      if params.enable_path && group == 0 {
        // pseudotime progress along the path
        let step = params.n_path_steps.max(1)
        let progress = (c % step).to_double() / step.to_double()
        log_mean = log_mean +
          progress * @math.ln(base_means[g].max(1.0e-15) * 10.0)
      }
      let mu = @math.exp(log_mean) * batch_factors[batch][g] * lib_factors[c]
      let dispersion = params.bcv_common + 1.0 / mu.max(1.0e-15).sqrt()
      // NB = Poisson(gamma mixture)
      let (s1, gamma_scale) = splat_gamma(1.0 / dispersion.max(1.0e-9), state)
      state = s1
      let lambda = mu * gamma_scale * dispersion.max(1.0e-9)
      let (s2, count) = splat_poisson(lambda, state)
      state = s2
      row.push(count)
    }
    counts.push(row)
  }
  // 7. dropout
  let mut n_dropouts = 0
  if params.enable_dropout {
    for g in 0..<n_genes {
      for c in 0..<n_cells {
        let group = cell_group[c]
        let log10_mean = @math.log10(
          @math.exp(group_log_means[group][g]).max(1.0e-15),
        )
        let probability = splat_dropout_probability(
          log10_mean,
          params.dropout_mid,
          params.dropout_shape,
        )
        let (s1, uniform) = splat_lcg(state)
        state = s1
        if uniform < probability {
          counts[g][c] = 0
          n_dropouts = n_dropouts + 1
        }
      }
    }
  }
  let pseudotime : Array[Double] = Array::new(capacity=n_cells)
  for c in 0..<n_cells {
    if params.enable_path && cell_group[c] == 0 {
      let step = params.n_path_steps.max(1)
      pseudotime.push((c % step).to_double() / step.to_double())
    } else {
      pseudotime.push(0.0)
    }
  }
  SplatResult::{
    counts,
    gene_base_means: base_means,
    cell_lib_factors: lib_factors,
    cell_group,
    cell_batch,
    cell_pseudotime: pseudotime,
    n_dropout_events: n_dropouts,
    group_log_means,
  }
}
