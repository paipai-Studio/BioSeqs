/// skbio_stats_validate.mbt — Bootstrap, adonis, DissimilarityMatrix
/// Items 15-17: skbio_bootstrap, skbio_adonis, skbio_dissim_matrix

///|
/// Bootstrap resampling result.
pub struct SkbioBootstrapResult {
  sample_indices : Array[Int32] // indices of resampled data
  statistic : Double // statistic computed on this resample
}

///|
/// Bootstrap resampling: generate n_resamples resampled datasets
/// and compute a statistic on each.
/// data: the original data values
/// stat_fn: function that takes a subsample and returns a statistic
pub fn skbio_bootstrap(
  data : Array[Double],
  n_resamples~ : Int32 = 1000,
  seed~ : Int32 = 42,
) -> Array[SkbioBootstrapResult] {
  let results : Array[SkbioBootstrapResult] = Array(capacity=n_resamples)
  let n = data.length()
  if n == 0 {
    return results
  }
  let mut rng_state = seed
  for r = 0; r < n_resamples; r = r + 1 {
    let indices : Array[Int32] = Array(capacity=n)
    let mut sum = 0.0
    for i = 0; i < n; i = i + 1 {
      // Simple LCG random
      rng_state = (rng_state * 1103515245 + 12345) & 0x7FFFFFFF
      let idx = rng_state % n
      indices.push(idx)
      sum = sum + data[idx]
    }
    let stat = sum / n.to_double()
    results.push(SkbioBootstrapResult::{
      sample_indices: indices,
      statistic: stat,
    })
  }
  results
}

///|
/// Compute bootstrap confidence interval from bootstrap results.
pub fn skbio_bootstrap_ci(
  results : Array[SkbioBootstrapResult],
  confidence~ : Double = 0.95,
) -> (Double, Double) {
  if results.length() == 0 {
    return (0.0, 0.0)
  }
  let stats : Array[Double] = Array(capacity=results.length())
  for i = 0; i < results.length(); i = i + 1 {
    stats.push(results[i].statistic)
  }
  stats.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  let alpha = 1.0 - confidence
  let lo_idx = (alpha / 2.0 * stats.length().to_double()).floor().to_int()
  let hi_idx = ((1.0 - alpha / 2.0) * stats.length().to_double())
    .floor()
    .to_int()
  let lo = if lo_idx < stats.length() { stats[lo_idx] } else { 0.0 }
  let hi = if hi_idx < stats.length() { stats[hi_idx] } else { 0.0 }
  (lo, hi)
}

///|
/// Adonis result (PERMANOVA with covariates).
pub struct SkbioAdonisResult {
  term : String
  df : Int
  ss : Double // sum of squares
  r_squared : Double // R^2 (effect size)
  f_stat : Double // F statistic
  p_value : Double
}

///|
/// Adonis: PERMANOVA with covariate decomposition (顺序检验, vegan::adonis2 风格).
///
/// Decomposes distance matrix variation into components explained by each
/// covariate, in **sequential** (order-dependent) order — the order you pass
/// `groupings` is the order they enter the model.
///
/// **Mathematics:** 先对距离矩阵做 Gower double-centering 得到 B，然后对每
/// 个分类变量做 dummy coding，用投影矩阵 Pₖ = Xₖ(Xₖ'Xₖ)⁻¹Xₖ' 计算增量 SS：
/// ΔSSₖ = trace(Pₖ B Pₖ) − trace(Pₖ₋₁ B Pₖ₋₁)。F 统计量 = (ΔSSₖ/dfₖ) /
/// (SS_resid/df_resid)。
///
/// `n_perm > 0` 时返回真正的 permutation p-value（在 strata 内 Fisher-Yates
/// shuffle 第一个分组变量，跨 strata 不动），`n_perm = 0` 时 p-value 为 NaN。
///
/// distance_matrix: n x n 距离矩阵
/// groupings: list of (name, assignments) for each covariate (顺序检验)
/// n_perm: permutation 次数 (默认 999)
/// strata: optional blocking / strata 变量（长度 n），在每个 strata 内独立 shuffle
/// seed: permutation 的 RNG 种子
pub fn skbio_adonis(
  distance_matrix : Array[Array[Double]],
  groupings : Array[(String, Array[Int32])],
  n_perm~ : Int32 = 999,
  strata~ : Array[Int32]? = None,
  seed~ : Int32 = 42,
) -> Array[SkbioAdonisResult] {
  let n = distance_matrix.length()
  let results : Array[SkbioAdonisResult] = Array(capacity=0)
  if n == 0 {
    return results
  }
  if groupings.length() == 0 {
    // 退化情况：空模型（残差 = 总 SS）
    let b = skbio_center_distance_matrix(distance_matrix)
    let ss_total = skbd_trace(b)
    results.push(SkbioAdonisResult::{
      term: "Residual".to_string(),
      df: n - 1,
      ss: ss_total,
      r_squared: 1.0,
      f_stat: 0.0 / 0.0,
      p_value: 0.0 / 0.0,
    })
    return results
  }

  // Gower center → B 矩阵
  let b = skbio_center_distance_matrix(distance_matrix)
  let ss_total = skbd_trace(b)

  // 累积设计矩阵
  let mut cum_x : Array[Array[Double]] = []
  let mut ss_cum = 0.0
  let mut df_cum : Int32 = 0

  for g = 0; g < groupings.length(); g = g + 1 {
    let (name, groups) = groupings[g]
    let new_x = skbd_dummy(groups)
    let new_df = new_x[0].length()
    let cur_x = if cum_x.length() > 0 { skbd_hstack(cum_x, new_x) } else { new_x }
    let cur_cols = cur_x[0].length()
    let ss_full = skbd_proj_trace(cur_x, b)
    let delta_ss = ss_full - ss_cum
    ss_cum = ss_full
    df_cum = cur_cols

    let df_resid = n - cur_cols
    let ss_resid = ss_total - ss_cum
    let obs_f = if df_resid > 0 && new_df > 0 {
      (delta_ss / new_df.to_double()) / (ss_resid / df_resid.to_double())
    } else {
      0.0 / 0.0
    }

    // permutation
    let mut count_better = 0
    let mut rng = seed
    let first_grp = groupings[0].1
    for _p = 0; _p < n_perm; _p = _p + 1 {
      let shuffled = skbd_stratified_shuffle(first_grp, strata, rng)
      rng = skbd_a_lcg(rng)
      // 重建前 g+1 个变量的累积 X（替换第一个为 shuffled）
      let mut perm_cum : Array[Array[Double]] = []
      for k = 0; k <= g; k = k + 1 {
        let p_groups = if k == 0 { shuffled } else { groupings[k].1 }
        let pk = skbd_dummy(p_groups)
        perm_cum = if perm_cum.length() > 0 { skbd_hstack(perm_cum, pk) } else { pk }
      }
      let mut perm_prev : Array[Array[Double]] = []
      for k = 0; k < g; k = k + 1 {
        let p_groups = if k == 0 { shuffled } else { groupings[k].1 }
        let pk = skbd_dummy(p_groups)
        perm_prev = if perm_prev.length() > 0 { skbd_hstack(perm_prev, pk) } else { pk }
      }
      let perm_ss_full = skbd_proj_trace(perm_cum, b)
      let perm_ss_prev = skbd_proj_trace(perm_prev, b)
      let perm_delta = perm_ss_full - perm_ss_prev
      let perm_df_res = n - cur_cols
      let perm_ss_res = ss_total - perm_ss_full
      let perm_f = if perm_df_res > 0 && new_df > 0 {
        (perm_delta / new_df.to_double()) / (perm_ss_res / perm_df_res.to_double())
      } else {
        0.0 / 0.0
      }
      if perm_f >= obs_f {
        count_better = count_better + 1
      }
    }

    let p_value = if n_perm == 0 {
      0.0 / 0.0
    } else {
      (count_better + 1).to_double() / (n_perm + 1).to_double()
    }
    let r_sq = if ss_total > 0.0 { delta_ss / ss_total } else { 0.0 }
    results.push(SkbioAdonisResult::{
      term: name,
      df: new_df,
      ss: delta_ss,
      r_squared: r_sq,
      f_stat: obs_f,
      p_value,
    })
    cum_x = cur_x
  }

  // Residual 行
  let df_res_total = n - df_cum
  let ss_res_total = ss_total - ss_cum
  results.push(SkbioAdonisResult::{
    term: "Residual".to_string(),
    df: if df_res_total > 0 { df_res_total } else { 0 },
    ss: ss_res_total,
    r_squared: if ss_total > 0.0 { ss_res_total / ss_total } else { 0.0 },
    f_stat: 0.0 / 0.0,
    p_value: 1.0,
  })
  results
}

// ---------------------------------------------------------------------------
// ADONIS 内部工具（与 skbio_ordination / skbio_distance 风格一致）
// ---------------------------------------------------------------------------

///|
fn skbd_trace(m : Array[Array[Double]]) -> Double {
  let mut s = 0.0
  for i = 0; i < m.length(); i = i + 1 {
    s = s + m[i][i]
  }
  s
}

///|
/// categorical grouping → dummy design matrix（levels−1 columns，首个 level 为 baseline）。
fn skbd_dummy(grouping : Array[Int32]) -> Array[Array[Double]] {
  let n = grouping.length()
  let levels : Array[Int32] = Array(capacity=0)
  for g in grouping {
    let mut dup = false
    for lv in levels {
      if lv == g { dup = true; break }
    }
    if !dup { levels.push(g) }
  }
  let p = if levels.length() > 1 { levels.length() - 1 } else { 0 }
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(p, 0.0)
    for k = 0; k < p; k = k + 1 {
      if grouping[i] == levels[k + 1] {
        row[k] = 1.0
      }
    }
    out.push(row)
  }
  out
}

///|
/// 水平堆叠两个设计矩阵（n 相同）。
fn skbd_hstack(a : Array[Array[Double]], b : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = a.length()
  let pa = if n > 0 { a[0].length() } else { 0 }
  let pb = if b.length() > 0 { b[0].length() } else { 0 }
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for j = 0; j < pa; j = j + 1 { row.push(a[i][j]) }
    for j = 0; j < pb; j = j + 1 { row.push(b[i][j]) }
    out.push(row)
  }
  out
}

///|
/// SS = trace(P B) where P = X(X'X)^-1 X' 是投影矩阵。
fn skbd_proj_trace(
  x : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Double {
  let n = b.length()
  let p = if x.length() > 0 { x[0].length() } else { 0 }
  if p == 0 { return 0.0 }
  let xt = skb_transpose(x)
  let xtx_inv = linalg_inverse(skb_mul(xt, x))
  let p_mat = skb_mul(skb_mul(x, xtx_inv), xt)
  let pb = skb_mul(p_mat, b)
  let mut trace_v = 0.0
  for i = 0; i < n; i = i + 1 { trace_v = trace_v + pb[i][i] }
  trace_v
}

///|
fn skbd_a_lcg(seed : Int32) -> Int32 {
  (seed * 1103515245 + 12345) & 0x7FFFFFFF
}

///|
/// 按 strata 分层 Fisher–Yates shuffle；strata=None 时退化为全局 shuffle。
fn skbd_stratified_shuffle(
  arr : Array[Int32],
  strata : Array[Int32]?,
  seed : Int32,
) -> Array[Int32] {
  match strata {
    None => {
      let out = arr.copy()
      let mut s = seed
      let n = out.length()
      let mut i = n - 1
      while i > 0 {
        s = skbd_a_lcg(s)
        let j = s % (i + 1)
        let t = out[i]; out[i] = out[j]; out[j] = t
        i = i - 1
      }
      out
    }
    Some(str) => {
      let out = arr.copy()
      let mut s = seed
      let seen : Array[Int32] = Array(capacity=0)
      for v in str {
        let mut dup = false
        for lv in seen {
          if lv == v { dup = true; break }
        }
        if !dup { seen.push(v) }
      }
      for lv in seen {
        let members : Array[Int32] = Array(capacity=0)
        for i = 0; i < str.length(); i = i + 1 {
          if str[i] == lv { members.push(i) }
        }
        let mut vi = members.length() - 1
        while vi > 0 {
          s = skbd_a_lcg(s)
          let vj = s % (vi + 1)
          let vt = members[vi]; members[vi] = members[vj]; members[vj] = vt
          vi = vi - 1
        }
        let original_vals : Array[Int32] = Array(capacity=0)
        for m in members { original_vals.push(arr[m]) }
        let mut oi = original_vals.length() - 1
        while oi > 0 {
          s = skbd_a_lcg(s)
          let oj = s % (oi + 1)
          let ot = original_vals[oi]
          original_vals[oi] = original_vals[oj]
          original_vals[oj] = ot
          oi = oi - 1
        }
        for idx = 0; idx < members.length(); idx = idx + 1 {
          out[members[idx]] = original_vals[idx]
        }
      }
      out
    }
  }
}

///|
/// DissimilarityMatrix: structured distance matrix.
pub struct SkbioDissimMatrix {
  ids : Array[String] // sample identifiers
  data : Array[Array[Double]] // n x n distance matrix
}

///|
/// Create a DissimilarityMatrix from ids and data.
pub fn SkbioDissimMatrix::new(
  ids : Array[String],
  data : Array[Array[Double]],
) -> SkbioDissimMatrix {
  SkbioDissimMatrix::{ ids, data, }
}

///|
/// Get the size (number of samples) of the matrix.
pub fn SkbioDissimMatrix::size(self : SkbioDissimMatrix) -> Int32 {
  self.ids.length()
}

///|
/// Get distance between two samples by ID.
pub fn SkbioDissimMatrix::get(
  self : SkbioDissimMatrix,
  id1 : String,
  id2 : String,
) -> Double {
  let mut idx1 = -1
  let mut idx2 = -1
  for i = 0; i < self.ids.length(); i = i + 1 {
    if self.ids[i] == id1 {
      idx1 = i
    }
    if self.ids[i] == id2 {
      idx2 = i
    }
  }
  if idx1 < 0 || idx2 < 0 {
    return -1.0
  }
  self.data[idx1][idx2]
}

///|
/// Check if the matrix is symmetric.
pub fn SkbioDissimMatrix::is_symmetric(self : SkbioDissimMatrix) -> Bool {
  let n = self.ids.length()
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      let diff = self.data[i][j] - self.data[j][i]
      if diff > 0.0000000001 || diff < -0.0000000001 {
        return false
      }
    }
  }
  true
}

///|
/// Get condensed (upper triangle) form of the matrix.
pub fn SkbioDissimMatrix::condensed(self : SkbioDissimMatrix) -> Array[Double] {
  let n = self.ids.length()
  let result : Array[Double] = Array(capacity=n * (n - 1) / 2)
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      result.push(self.data[i][j])
    }
  }
  result
}

///|
/// Get mean distance within a group of sample IDs.
pub fn SkbioDissimMatrix::within(
  self : SkbioDissimMatrix,
  group_ids : Array[String],
) -> Double {
  let mut sum = 0.0
  let mut count = 0.0
  for i = 0; i < group_ids.length(); i = i + 1 {
    for j = i + 1; j < group_ids.length(); j = j + 1 {
      let d = self.get(group_ids[i], group_ids[j])
      if d >= 0.0 {
        sum = sum + d
        count = count + 1.0
      }
    }
  }
  if count == 0.0 {
    0.0
  } else {
    sum / count
  }
}

///|
/// Get mean distance between two groups of sample IDs.
pub fn SkbioDissimMatrix::between(
  self : SkbioDissimMatrix,
  group1 : Array[String],
  group2 : Array[String],
) -> Double {
  let mut sum = 0.0
  let mut count = 0.0
  for i = 0; i < group1.length(); i = i + 1 {
    for j = 0; j < group2.length(); j = j + 1 {
      let d = self.get(group1[i], group2[j])
      if d >= 0.0 {
        sum = sum + d
        count = count + 1.0
      }
    }
  }
  if count == 0.0 {
    0.0
  } else {
    sum / count
  }
}
