///|
/// scanpy `tl.rank_genes_groups` with a `method=` dispatcher.  Unifies the
/// three differential-expression flavours scanpy offers into one entry point
/// returning a common result shape: `t-test` (Welch), `wilcoxon`
/// (rank-sum with tie correction) and `logreg` (L2-regularized logistic
/// regression, scores = absolute coefficients; no p-values are defined for
/// this method and they are reported as 1.0).  Results are sorted by score in
/// descending order (most group-upregulated first), optionally truncated to
/// `n_top`, and p-values get Benjamini-Hochberg correction.

///|
/// Unified differential-expression result across methods.
pub struct ScanpyRankGenesMethodResult {
  method_ : String
  gene_indices : Array[Int32] // sorted by score descending
  scores : Array[Double]
  logfoldchanges : Array[Double]
  pvals : Array[Double]
  pvals_adj : Array[Double]
}

///|
/// Benjamini-Hochberg FDR correction.  Returns adjusted p-values in the
/// original order (values are clipped to 1.0).
fn scanpy_bh_adjust(pvals : Array[Double]) -> Array[Double] {
  let n = pvals.length()
  let out = Array::make(n, 1.0)
  if n == 0 {
    return out
  }
  let order = Array::make(n, 0)
  for t = 0; t < n; t = t + 1 {
    order[t] = t
  }
  order.sort_by(fn(a, b) {
    let pa = pvals[a]
    let pb = pvals[b]
    if pa < pb {
      -1
    } else if pa > pb {
      1
    } else {
      0
    }
  })
  let mut prev = 1.0
  let mut t = n - 1
  while t >= 0 {
    let idx = order[t]
    let rank = t + 1
    let mut q = pvals[idx] * n.to_double() / rank.to_double()
    if q > prev {
      q = prev
    }
    prev = q
    out[idx] = if q > 1.0 { 1.0 } else { q }
    t = t - 1
  }
  out
}

///|
/// Mean of column `j` over the given rows.
fn scanpy_rgm_mean_var(rows : Array[Array[Double]], j : Int32) -> Double {
  let n = rows.length()
  if n == 0 {
    return 0.0
  }
  let mut s = 0.0
  for r in rows {
    s = s + r[j]
  }
  s / n.to_double()
}

///|
/// Extract the rows of `x` whose group label equals `label`.
fn scanpy_rgm_pick_rows(
  x : Array[Array[Double]],
  groups : Array[String],
  label : String,
) -> Array[Array[Double]] {
  let rows : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < groups.length(); i = i + 1 {
    if groups[i] == label && i < x.length() {
      rows.push(x[i])
    }
  }
  rows
}

///|
/// L2-regularized binary logistic regression trained on the standardized
/// columns of the stacked group/rest matrix; returns the absolute coefficient
/// of every feature (gene) as its score.  Full-batch gradient descent with a
/// fixed learning rate; the class balance is left unweighted, matching
/// scanpy's `method='logreg'` default of an unweighted `LogisticRegression`.
fn scanpy_logreg_abs_coefs(
  x_group : Array[Array[Double]],
  x_rest : Array[Array[Double]],
  max_iter~ : Int32 = 200,
  l2~ : Double = 0.01,
  lr~ : Double = 0.1,
) -> Array[Double] {
  let n1 = x_group.length()
  let n0 = x_rest.length()
  let n = n1 + n0
  let g = if n > 0 && n1 + n0 > 0 {
    if n1 > 0 {
      x_group[0].length()
    } else {
      x_rest[0].length()
    }
  } else {
    0
  }
  let scores = Array::make(g, 0.0)
  if n == 0 || g == 0 {
    return scores
  }
  // stack and standardize columns
  let xs : Array[Array[Double]] = Array(capacity=n)
  for r in x_group {
    xs.push(r)
  }
  for r in x_rest {
    xs.push(r)
  }
  let means = Array::make(g, 0.0)
  let stds = Array::make(g, 1.0)
  for j = 0; j < g; j = j + 1 {
    let mut s = 0.0
    for i = 0; i < n; i = i + 1 {
      s = s + xs[i][j]
    }
    means[j] = s / n.to_double()
    let mut v = 0.0
    for i = 0; i < n; i = i + 1 {
      let d = xs[i][j] - means[j]
      v = v + d * d
    }
    let sd = (v / (n - 1).to_double()).sqrt()
    stds[j] = if sd > 1.0e-12 { sd } else { 1.0 }
  }
  let w = Array::make(g, 0.0)
  let mut b = 0.0
  let n_f = n.to_double()
  for _iter = 0; _iter < max_iter; _iter = _iter + 1 {
    let grad_w = Array::make(g, 0.0)
    let mut grad_b = 0.0
    for i = 0; i < n; i = i + 1 {
      let y = if i < n1 { 1.0 } else { 0.0 }
      let mut z = b
      for j = 0; j < g; j = j + 1 {
        z = z + w[j] * (xs[i][j] - means[j]) / stds[j]
      }
      // numerically stable sigmoid
      let p = if z >= 30.0 {
        1.0
      } else if z <= -30.0 {
        0.0
      } else {
        1.0 / (1.0 + @math.exp(-z))
      }
      let e = p - y
      grad_b = grad_b + e
      for j = 0; j < g; j = j + 1 {
        grad_w[j] = grad_w[j] + e * (xs[i][j] - means[j]) / stds[j]
      }
    }
    for j = 0; j < g; j = j + 1 {
      w[j] = w[j] - lr * (grad_w[j] / n_f + l2 * w[j])
    }
    b = b - lr * (grad_b / n_f)
  }
  for j = 0; j < g; j = j + 1 {
    scores[j] = if w[j] < 0.0 { -w[j] } else { w[j] }
  }
  scores
}

///|
/// Rank genes between two cell groups with a selectable method
/// (`rank_genes_groups(method=...)`):
///
/// - `"t-test"` — Welch's t-test (score = t statistic, two-sided p-value)
/// - `"wilcoxon"` — rank-sum test with tie correction (score = z statistic,
///   two-sided normal-approximation p-value)
/// - `"logreg"` — L2 logistic regression (score = |coefficient|; p-values are
///   not defined and reported as 1.0)
///
/// Log-fold-changes are `log2((mean_target + 1) / (mean_ref + 1))` throughout.
/// Rows are sorted by score in descending order (group-upregulated first) and
/// truncated to `n_top` when positive.  Adjusted p-values use
/// Benjamini-Hochberg over all genes.
pub fn scanpy_rank_genes_groups_method(
  x : Array[Array[Double]],
  groups : Array[String],
  ref_label : String,
  target_label : String,
  method_~ : String = "t-test",
  n_top~ : Int32 = 0,
) -> ScanpyRankGenesMethodResult {
  let n_genes = if x.length() > 0 { x[0].length() } else { 0 }
  let x_group = scanpy_rgm_pick_rows(x, groups, target_label)
  let x_rest = scanpy_rgm_pick_rows(x, groups, ref_label)
  let scores = Array::make(n_genes, 0.0)
  let pvals = Array::make(n_genes, 1.0)
  let lfcs = Array::make(n_genes, 0.0)
  let resolved = match method_ {
    "wilcoxon" => "wilcoxon"
    "logreg" => "logreg"
    _ => "t-test"
  }
  if resolved == "t-test" {
    let per_gene = scanpy_rank_genes_groups(x, groups, ref_label, target_label)
    for j = 0; j < n_genes; j = j + 1 {
      scores[j] = per_gene[j].t_stat
      pvals[j] = per_gene[j].p_value
      lfcs[j] = per_gene[j].logfoldchange
    }
  } else if resolved == "wilcoxon" {
    let per_gene = scanpy_rank_genes_wilcoxon_tie_correct(x_group, x_rest)
    let mean_g = Array::make(n_genes, 0.0)
    let mean_r = Array::make(n_genes, 0.0)
    for j = 0; j < n_genes; j = j + 1 {
      mean_g[j] = scanpy_rgm_mean_var(x_group, j)
      mean_r[j] = scanpy_rgm_mean_var(x_rest, j)
    }
    for kv in per_gene {
      let j = kv.0
      scores[j] = kv.1
      pvals[j] = normal_pvalue_two_sided(kv.1)
      lfcs[j] = @math.ln((mean_g[j] + 1.0) / (mean_r[j] + 1.0)) / @math.ln(2.0)
    }
  } else {
    let coefs = scanpy_logreg_abs_coefs(x_group, x_rest)
    let mean_g = Array::make(n_genes, 0.0)
    let mean_r = Array::make(n_genes, 0.0)
    for j = 0; j < n_genes; j = j + 1 {
      mean_g[j] = scanpy_rgm_mean_var(x_group, j)
      mean_r[j] = scanpy_rgm_mean_var(x_rest, j)
      scores[j] = coefs[j]
      lfcs[j] = @math.ln((mean_g[j] + 1.0) / (mean_r[j] + 1.0)) / @math.ln(2.0)
    }
  }
  let adj = scanpy_bh_adjust(pvals)
  // order by score descending, then by gene index ascending (stable)
  let order = Array::make(n_genes, 0)
  for j = 0; j < n_genes; j = j + 1 {
    order[j] = j
  }
  order.sort_by(fn(a, b) {
    if scores[a] > scores[b] {
      -1
    } else if scores[a] < scores[b] {
      1
    } else if a < b {
      -1
    } else if a > b {
      1
    } else {
      0
    }
  })
  let take = if n_top > 0 && n_top < n_genes { n_top } else { n_genes }
  let gene_indices : Array[Int32] = Array(capacity=take)
  let out_scores : Array[Double] = Array(capacity=take)
  let out_lfcs : Array[Double] = Array(capacity=take)
  let out_pvals : Array[Double] = Array(capacity=take)
  let out_adj : Array[Double] = Array(capacity=take)
  for t = 0; t < take; t = t + 1 {
    let j = order[t]
    gene_indices.push(j)
    out_scores.push(scores[j])
    out_lfcs.push(lfcs[j])
    out_pvals.push(pvals[j])
    out_adj.push(adj[j])
  }
  ScanpyRankGenesMethodResult::{
    method_: resolved,
    gene_indices,
    scores: out_scores,
    logfoldchanges: out_lfcs,
    pvals: out_pvals,
    pvals_adj: out_adj,
  }
}
