///|
/// scanpy `tools` 核心（第一批）：`rank_genes_groups` 的 t-test 差异表达 —
/// 对每个基因，比较 `target_label` 与 `ref_label` 两组细胞的 Welch t 检验，
/// 返回 t 统计量、p 值（经项目 `stat_pvalue_t` 的正态近似，大样本下准确）、
/// 两组均值与 log2(伪计数+1) 倍率变化。语义对齐 scanpy
/// `_rank_genes_groups.py` 的 t-test flavor（Welch，不等方差）。

///|
pub struct RankGenesResult {
  t_stat : Double
  p_value : Double
  mean_ref : Double
  mean_target : Double
  logfoldchange : Double
}

///|
/// Mean and sample variance (ddof=1) of column `j` over the given row indices.
fn sub_mean_var(x : Array[Array[Double]], j : Int, idx : Array[Int]) -> (Double, Double) {
  let n = idx.length()
  if n == 0 {
    return (0.0, 0.0)
  }
  let mut s = 0.0
  for r in idx {
    s = s + x[r][j]
  }
  let mean = s / n.to_double()
  if n == 1 {
    return (mean, 0.0)
  }
  let mut ss = 0.0
  for r in idx {
    let d = x[r][j] - mean
    ss = ss + d * d
  }
  (mean, ss / (n - 1).to_double())
}

///|
/// Welch two-sample t-test p-value for fractional df.
fn welch_p(t : Double, df : Double) -> Double {
  if df < 1.0 {
    1.0
  } else {
    stat_pvalue_t(t, df.to_int())
  }
}

///|
/// Differential expression of every gene between two cell groups using Welch's
/// t-test (`rank_genes_groups(method="t-test")`). `groups` aligns with `x` rows;
/// every row not in either group is ignored. Returns one result per gene (per
/// `x` column) in column order.
pub fn scanpy_rank_genes_groups(
  x : Array[Array[Double]],
  groups : Array[String],
  ref_label : String,
  target_label : String,
) -> Array[RankGenesResult] {
  // collect row indices for the two groups
  let ref_idx : Array[Int] = Array(capacity=0)
  let target_idx : Array[Int] = Array(capacity=0)
  for i = 0; i < groups.length(); i = i + 1 {
    if groups[i] == ref_label {
      ref_idx.push(i)
    } else if groups[i] == target_label {
      target_idx.push(i)
    }
  }
  let n_genes = if x.length() == 0 { 0 } else { x[0].length() }
  let n_ref = ref_idx.length()
  let n_target = target_idx.length()
  let out : Array[RankGenesResult] = Array(capacity=0)
  let mut j = 0
  while j < n_genes {
    let (mean_ref, var_ref) = sub_mean_var(x, j, ref_idx)
    let (mean_target, var_target) = sub_mean_var(x, j, target_idx)
    let denom = if n_ref > 0 && n_target > 0 {
      var_ref / n_ref.to_double() + var_target / n_target.to_double()
    } else {
      0.0
    }
    let (t, p) = if denom == 0.0 {
      (0.0, 1.0)
    } else {
      let tt = (mean_target - mean_ref) / denom.sqrt()
      let df = if n_ref > 1 && n_target > 1 {
        let a = var_ref / n_ref.to_double()
        let b = var_target / n_target.to_double()
        (a + b) * (a + b) /
          (a * a / (n_ref - 1).to_double() + b * b / (n_target - 1).to_double())
      } else {
        1.0
      }
      (tt, welch_p(tt, df))
    }
    let lfc = @math.ln((mean_target + 1.0) / (mean_ref + 1.0)) / @math.ln(2.0)
    out.push(
      RankGenesResult::{ t_stat: t, p_value: p, mean_ref, mean_target, logfoldchange: lfc },
    )
    j = j + 1
  }
  out
}