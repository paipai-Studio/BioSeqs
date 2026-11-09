///|
/// scanpy `pp.highly_variable_genes(flavor='seurat_v3')` — the original
/// Seurat v3 variance-stabilizing-transform path (Stuart et al. 2019),
/// implemented from first principles: clip counts, normalize to the median
/// depth, `asinh` transform, loess-fit `log10(variance)` against
/// `log10(mean)`, and rank genes by the ratio of observed to expected
/// (stabilized) variance.

///|
/// One-dimensional local-linear (loess-style) regression: for every query
/// point the `k = ceil(span * n)` nearest x-neighbours are weighted by the
/// tricube kernel and a weighted linear model is fitted; returns the fitted
/// values at the input x positions.  `xs` must contain at least 2 points.
fn scanpy_loess_1d(
  xs : Array[Double],
  ys : Array[Double],
  span~ : Double = 0.3,
) -> Array[Double] {
  let n = xs.length()
  let out = Array::make(n, 0.0)
  if n == 0 {
    return out
  }
  if n == 1 {
    out[0] = ys[0]
    return out
  }
  let order = Array::make(n, 0)
  for t = 0; t < n; t = t + 1 {
    order[t] = t
  }
  order.sort_by(fn(a, b) {
    if xs[a] < xs[b] {
      -1
    } else if xs[a] > xs[b] {
      1
    } else if a < b {
      -1
    } else if a > b {
      1
    } else {
      0
    }
  })
  let k = {
    let raw = (span * n.to_double()).to_int()
    let kk = if raw < 2 { 2 } else { raw }
    if kk > n {
      n
    } else {
      kk
    }
  }
  for qi = 0; qi < n; qi = qi + 1 {
    // find the window of k nearest x values around xs[qi] in sorted order
    let pos = {
      let mut p = 0
      for t = 0; t < n; t = t + 1 {
        if order[t] == qi {
          p = t
          break
        }
      }
      p
    }
    let mut lo = pos - k / 2
    if lo < 0 {
      lo = 0
    }
    if lo + k > n {
      lo = n - k
    }
    let x0 = xs[order[lo]]
    let x1 = xs[order[lo + k - 1]]
    let width = if x1 - x0 > 1.0e-12 { x1 - x0 } else { 1.0 }
    // weighted linear fit sum_j w_j (a + b x_j) ~ y_j
    let mut sw = 0.0
    let mut swx = 0.0
    let mut swy = 0.0
    let mut swxx = 0.0
    let mut swxy = 0.0
    for t = lo; t < lo + k; t = t + 1 {
      let j = order[t]
      let u = ((xs[j] - xs[qi]) / width).abs()
      let u1 = 1.0 - u * u * u
      let w = if u1 > 0.0 { u1 * u1 * u1 } else { 0.0 }
      if w <= 0.0 {
        continue
      }
      sw = sw + w
      swx = swx + w * xs[j]
      swy = swy + w * ys[j]
      swxx = swxx + w * xs[j] * xs[j]
      swxy = swxy + w * xs[j] * ys[j]
    }
    let det = sw * swxx - swx * swx
    if det.abs() < 1.0e-12 {
      out[qi] = swy / (if sw > 0.0 { sw } else { 1.0 })
    } else {
      let a = (swy * swxx - swx * swxy) / det
      let b = (sw * swxy - swx * swy) / det
      out[qi] = a + b * xs[qi]
    }
  }
  out
}

///|
/// Result of the seurat_v3 VST HVG selection over the genes that passed
/// `min_cells` (`kept_genes` indexes into the input columns).
pub struct ScanpyHvgSeuratV3 {
  highly_variable : Array[Bool]
  means : Array[Double]
  variances : Array[Double]
  residual_variance : Array[Double]
  kept_genes : Array[Int32]
}

///|
/// Seurat-v3 highly variable genes via the original VST route:
///
/// 1. keep genes expressed in at least `min_cells` cells,
/// 2. clip counts at `sqrt(mean total counts per cell)`,
/// 3. normalize each cell to the median count depth and `asinh`-transform,
/// 4. loess-fit `log10(variance)` on `log10(mean)` (span 0.3),
/// 5. score genes by `variance / 10^fit` (stabilized variance ratio) and mark
///    the `n_top_genes` highest as highly variable.
pub fn scanpy_hvg_seurat_v3(
  x : Array[Array[Double]],
  n_top_genes~ : Int32 = 2000,
  min_cells~ : Int32 = 3,
  span~ : Double = 0.3,
) -> ScanpyHvgSeuratV3 {
  let n_cells = x.length()
  let n_genes = if n_cells > 0 { x[0].length() } else { 0 }
  if n_cells == 0 || n_genes == 0 {
    return ScanpyHvgSeuratV3::{
      highly_variable: [],
      means: [],
      variances: [],
      residual_variance: [],
      kept_genes: [],
    }
  }
  // gene filter + count depth
  let totals : Array[Double] = Array(capacity=n_cells)
  for i = 0; i < n_cells; i = i + 1 {
    let mut s = 0.0
    for j = 0; j < n_genes; j = j + 1 {
      s = s + x[i][j]
    }
    totals.push(s)
  }
  let kept : Array[Int32] = Array(capacity=n_genes)
  for j = 0; j < n_genes; j = j + 1 {
    let mut cnt = 0
    for i = 0; i < n_cells; i = i + 1 {
      if x[i][j] > 0.0 {
        cnt = cnt + 1
      }
    }
    if cnt >= min_cells {
      kept.push(j)
    }
  }
  let nk = kept.length()
  if nk == 0 {
    return ScanpyHvgSeuratV3::{
      highly_variable: [],
      means: [],
      variances: [],
      residual_variance: [],
      kept_genes: [],
    }
  }
  let mut total_sum = 0.0
  for t in totals {
    total_sum = total_sum + t
  }
  let clip_max = (total_sum / n_cells.to_double()).sqrt()
  // median count depth
  let sorted_totals = totals.copy()
  sorted_totals.sort()
  let median_depth = if n_cells % 2 == 1 {
    sorted_totals[n_cells / 2]
  } else {
    (sorted_totals[n_cells / 2 - 1] + sorted_totals[n_cells / 2]) / 2.0
  }
  // clipped / normalized / asinh matrix over kept genes
  let means = Array::make(nk, 0.0)
  let vars = Array::make(nk, 0.0)
  for t = 0; t < nk; t = t + 1 {
    let j = kept[t]
    let mut s = 0.0
    for i = 0; i < n_cells; i = i + 1 {
      let c = if x[i][j] > clip_max { clip_max } else { x[i][j] }
      let norm = if totals[i] > 0.0 {
        c / totals[i] * median_depth
      } else {
        0.0
      }
      s = s + asinh_transform(norm)
    }
    let mean = s / n_cells.to_double()
    means[t] = mean
    let mut v = 0.0
    for i = 0; i < n_cells; i = i + 1 {
      let c = if x[i][j] > clip_max { clip_max } else { x[i][j] }
      let norm = if totals[i] > 0.0 {
        c / totals[i] * median_depth
      } else {
        0.0
      }
      let v2 = asinh_transform(norm)
      let d = v2 - mean
      v = v + d * d
    }
    vars[t] = v / (n_cells - 1).to_double()
  }
  // loess fit of log10(var) ~ log10(mean) over non-constant genes
  let log_mean : Array[Double] = Array(capacity=nk)
  let log_var : Array[Double] = Array(capacity=nk)
  let fit_idx : Array[Int32] = Array(capacity=nk)
  for t = 0; t < nk; t = t + 1 {
    if vars[t] > 0.0 && means[t] > 0.0 {
      fit_idx.push(t)
      log_mean.push(@math.log10(means[t]))
      log_var.push(@math.log10(vars[t]))
    }
  }
  let fitted = scanpy_loess_1d(log_mean, log_var, span~)
  let residual_variance = Array::make(nk, 0.0)
  for q = 0; q < fit_idx.length(); q = q + 1 {
    let t = fit_idx[q]
    let expected = @math.pow(10.0, fitted[q])
    let ratio = if expected > 1.0e-12 { vars[t] / expected } else { 0.0 }
    residual_variance[t] = if ratio > 1.0e12 { 1.0e12 } else { ratio }
  }
  // top n by stabilized variance
  let order = Array::make(nk, 0)
  for t = 0; t < nk; t = t + 1 {
    order[t] = t
  }
  order.sort_by(fn(a, b) {
    if residual_variance[a] > residual_variance[b] {
      -1
    } else if residual_variance[a] < residual_variance[b] {
      1
    } else if a < b {
      -1
    } else if a > b {
      1
    } else {
      0
    }
  })
  let n_select = if n_top_genes < nk { n_top_genes } else { nk }
  let highly_variable = Array::make(nk, false)
  for t = 0; t < n_select; t = t + 1 {
    highly_variable[order[t]] = true
  }
  ScanpyHvgSeuratV3::{
    highly_variable,
    means,
    variances: vars,
    residual_variance,
    kept_genes: kept,
  }
}

///|
/// `asinh` via the log formula (`ln(x + sqrt(x^2 + 1))`) — the repo has no
/// direct `@math.asinh`.
fn asinh_transform(x : Double) -> Double {
  if x == 0.0 {
    return 0.0
  }
  @math.ln(x + (x * x + 1.0).sqrt())
}
