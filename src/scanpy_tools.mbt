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

// ============================================================================
// tsne / dendrogram
// ============================================================================

///|
/// t-SNE embedding of `x` (n_obs × n_vars) into `dims` dimensions
/// (`tl.tsne`), a thin wrapper over the repo's Barnes–Hut t-SNE.
pub fn scanpy_tsne(
  x : Array[Array[Double]],
  perplexity~ : Double = 30.0,
  dims~ : Int = 2,
) -> Array[Array[Double]] {
  let config = TsneConfig::new_custom(perplexity, 1000, dims, 42)
  let result = tsne(x, config)
  result.embedding
}

///|
/// Hierarchical clustering of genes (`tl.dendrogram`): complete linkage over
/// the Euclidean distances between gene expression vectors, returned as a
/// merge/height/order tree.
pub fn scanpy_dendrogram(x : Array[Array[Double]]) -> HclustResult {
  let n_cells = x.length()
  let n_genes = if n_cells > 0 { x[0].length() } else { 0 }
  let dist : Array[Array[Double]] = Array::make(n_genes, [])
  for j = 0; j < n_genes; j = j + 1 {
    dist[j] = Array::make(n_genes, 0.0)
  }
  for a = 0; a < n_genes; a = a + 1 {
    for b = a + 1; b < n_genes; b = b + 1 {
      let mut d2 = 0.0
      for i = 0; i < n_cells; i = i + 1 {
        let diff = x[i][a] - x[i][b]
        d2 = d2 + diff * diff
      }
      let d = d2.sqrt()
      dist[a][b] = d
      dist[b][a] = d
    }
  }
  hclust_complete(dist)
}
