///|
/// scanpy `preprocessing` 核心：`normalize_total`（每细胞文库归一化）、`log1p`、
/// `scale`（每基因 z-score + 裁剪）、`highly_variable_genes`（Seurat flavor 的
/// 高变基因检测）。语义对齐 scanpy 源码（`_normalization.py` / `_simple.py` /
/// `_scale.py` / `_highly_variable_genes.py`），`x` 均为 n_obs × n_vars（行=细胞、
/// 列=基因）。

///|
fn col_mean(x : Array[Array[Double]], j : Int) -> Double {
  let n = x.length()
  if n == 0 {
    return 0.0
  }
  let mut s = 0.0
  let mut i = 0
  while i < n {
    s = s + x[i][j]
    i = i + 1
  }
  s / n.to_double()
}

///|
fn col_sample_var(x : Array[Array[Double]], j : Int, mean : Double) -> Double {
  let n = x.length()
  if n <= 1 {
    return 0.0
  }
  let mut s = 0.0
  let mut i = 0
  while i < n {
    let d = x[i][j] - mean
    s = s + d * d
    i = i + 1
  }
  s / (n - 1).to_double()
}

///|
/// Normalize each cell to `target_sum` total counts (`pp.normalize_total`):
/// `x[i][j] = x[i][j] * target_sum / sum(x[i])`. Cells with a zero total are left
/// at zero.
pub fn scanpy_normalize_total(
  x : Array[Array[Double]],
  target_sum : Double,
) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array::new()
  for row in x {
    let mut s = 0.0
    for v in row {
      s = s + v
    }
    let nor : Array[Double] = Array::new()
    if s == 0.0 {
      for v in row {
        nor.push(0.0)
      }
    } else {
      let f = target_sum / s
      for v in row {
        nor.push(v * f)
      }
    }
    out.push(nor)
  }
  out
}

///|
/// Natural-log transform `x -> ln(1 + x)` element-wise (`pp.log1p`).
pub fn scanpy_log1p(x : Array[Array[Double]]) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array::new()
  for row in x {
    let r : Array[Double] = Array::new()
    for v in row {
      r.push(@math.ln(1.0 + v))
    }
    out.push(r)
  }
  out
}

///|
/// Zero-center and unit-variance scale each gene across cells (`pp.scale`,
/// ddof=1). Genes with zero variance are set to 0. When `max_value` is given the
/// result is clipped to `[-max_value, max_value]` (matching `pp.scale(max_value=)`).
pub fn scanpy_scale(
  x : Array[Array[Double]],
  max_value : Double?,
) -> Array[Array[Double]] {
  let n_genes = if x.length() == 0 { 0 } else { x[0].length() }
  let out : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < x.length() {
    out.push(Array::make(n_genes, 0.0))
    i = i + 1
  }
  for j = 0; j < n_genes; j = j + 1 {
    let mean = col_mean(x, j)
    let variance = col_sample_var(x, j, mean)
    let std = if variance > 0.0 { variance.sqrt() } else { 0.0 }
    let mut r = 0
    while r < x.length() {
      let z = if std == 0.0 {
        0.0
      } else {
        (x[r][j] - mean) / std
      }
      let mut v = z
      match max_value {
        None => ()
        Some(mx) => {
          if v > mx {
            v = mx
          } else if v < -mx {
            v = -mx
          }
        }
      }
      out[r][j] = v
      r = r + 1
    }
  }
  out
}

///|
/// Seurat flavor `highly_variable_genes`.
pub struct HvgResult {
  means : Array[Double]
  dispersions : Array[Double]
  dispersions_norm : Array[Double]
  highly_variable : Array[Bool]
}

///|
fn dmin(a : Array[Double]) -> Double {
  let mut m = a[0]
  for v in a {
    if v < m {
      m = v
    }
  }
  m
}

///|
fn dmax(a : Array[Double]) -> Double {
  let mut m = a[0]
  for v in a {
    if v > m {
      m = v
    }
  }
  m
}

///|
fn is_nan(v : Double) -> Bool {
  v != v
}

///|
/// Identify highly variable genes using the Seurat `vst` heuristic (`flavor=
/// "seurat"`): per gene compute mean and sample variance across cells, dispersion
/// = var/mean, then log both dispersion and (log1p'd) mean; bin genes into
/// `n_bins` equal-width bins on the log1p mean, z-score the log dispersion within
/// each bin, and flag genes whose log1p mean falls in `(min_mean, max_mean)` and
/// whose binned z-score exceeds `min_disp`. Matches scanpy defaults
/// (n_bins=20, min_disp=0.5, min_mean=0.0125, max_mean=3).
pub fn scanpy_highly_variable_genes(
  x : Array[Array[Double]],
  n_bins : Int,
  min_mean : Double,
  max_mean : Double,
  min_disp : Double,
) -> HvgResult {
  let n_genes = if x.length() == 0 { 0 } else { x[0].length() }
  let means : Array[Double] = Array::new()
  let dispersions : Array[Double] = Array::new()
  let mut j = 0
  while j < n_genes {
    let raw_mean = col_mean(x, j)
    let variance = col_sample_var(x, j, raw_mean)
    let m_adj = if raw_mean == 0.0 { 1.0e-12 } else { raw_mean }
    let d = variance / m_adj
    means.push(@math.ln(1.0 + m_adj))
    if d == 0.0 {
      dispersions.push(0.0 / 0.0)
    } else {
      dispersions.push(@math.ln(d))
    }
    j = j + 1
  }
  let dispersions_norm : Array[Double] = Array::make(n_genes, 0.0)
  if n_genes > 0 {
    let mn = dmin(means)
    let mx = dmax(means)
    let width = if mx == mn { 1.0 } else { (mx - mn) / n_bins.to_double() }
    // per-bin accumulation of valid log-dispersions
    let bin_count : Array[Int] = Array::make(n_bins, 0)
    let bin_sum : Array[Double] = Array::make(n_bins, 0.0)
    let bin_sumsq : Array[Double] = Array::make(n_bins, 0.0)
    let gene_bin : Array[Int] = Array::make(n_genes, 0)
    let mut g = 0
    while g < n_genes {
      let b = if width == 1.0 {
        0
      } else {
        let t = ((means[g] - mn) / width).to_int()
        if t >= n_bins { n_bins - 1 } else if t < 0 { 0 } else { t }
      }
      gene_bin[g] = b
      if not(is_nan(dispersions[g])) {
        bin_count[b] = bin_count[b] + 1
        bin_sum[b] = bin_sum[b] + dispersions[g]
        bin_sumsq[b] = bin_sumsq[b] + dispersions[g] * dispersions[g]
      }
      g = g + 1
    }
    // z-score within bin
    let mut k = 0
    while k < n_genes {
      let b = gene_bin[k]
      if bin_count[b] >= 2 {
        let avg = bin_sum[b] / bin_count[b].to_double()
        let dev = ((bin_sumsq[b] - bin_count[b].to_double() * avg * avg) /
          (bin_count[b] - 1).to_double()).sqrt()
        if dev > 0.0 && not(is_nan(dispersions[k])) {
          dispersions_norm[k] = (dispersions[k] - avg) / dev
        }
      } else if not(is_nan(dispersions[k])) {
        // single valid gene in bin -> normalized dispersion of 1
        dispersions_norm[k] = 1.0
      } else {
        // constant gene (zero dispersion) -> NaN, excluded from HV
        dispersions_norm[k] = 0.0 / 0.0
      }
      k = k + 1
    }
  }
  let highly_variable : Array[Bool] = Array::make(n_genes, false)
  let mut h = 0
  while h < n_genes {
    if means[h] > min_mean && means[h] < max_mean && dispersions_norm[h] > min_disp {
      highly_variable[h] = true
    }
    h = h + 1
  }
  HvgResult::{ means, dispersions, dispersions_norm, highly_variable }
}