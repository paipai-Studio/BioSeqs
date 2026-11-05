/// Scanpy second-round: Advanced differential expression
/// Inspired by scanpy.tl.rank_genes_groups extensions

///|
/// Extended rank genes result with additional methods and stats.
pub struct RankGenesExtResult {
  gene_indices : Array[Int32]
  scores : Array[Double]
  logfoldchanges : Array[Double]
  pvals : Array[Double]
  pvals_adj : Array[Double]
  pts_group : Array[Double] // fraction expressing in group
  pts_rest : Array[Double] // fraction expressing in reference
  method_ : String
}

///|
/// Rank genes using t-test with overestimated variance.
/// This method (from Seurat) adds the mean of the group to the variance,
/// which helps with low-expressing genes.
pub fn scanpy_rank_genes_ttest_overestim(
  x_group : Array[Array[Double]],
  x_rest : Array[Array[Double]],
  gene_names : Array[String],
) -> RankGenesExtResult {
  let n_genes = gene_names.length()
  let n_group = x_group.length()
  let n_rest = x_rest.length()
  let gene_indices = Array(capacity=n_genes)
  let scores = Array(capacity=n_genes)
  let logfoldchanges = Array(capacity=n_genes)
  let pvals = Array(capacity=n_genes)
  let pvals_adj = Array(capacity=n_genes)
  let pts_group = Array(capacity=n_genes)
  let pts_rest = Array(capacity=n_genes)
  for g in 0..<n_genes {
    // Group stats
    let mut sum_g = 0.0
    let mut sum_g2 = 0.0
    let mut n_expr_g = 0
    for i in 0..<n_group {
      let v = x_group[i][g]
      sum_g = sum_g + v
      sum_g2 = sum_g2 + v * v
      if v > 0.0 {
        n_expr_g = n_expr_g + 1
      }
    }
    let mean_g = if n_group > 0 { sum_g / n_group.to_double() } else { 0.0 }
    let var_g = if n_group > 1 {
      (sum_g2 - sum_g * sum_g / n_group.to_double()) / (n_group - 1).to_double()
    } else {
      0.0
    }
    // Rest stats
    let mut sum_r = 0.0
    let mut sum_r2 = 0.0
    let mut n_expr_r = 0
    for i in 0..<n_rest {
      let v = x_rest[i][g]
      sum_r = sum_r + v
      sum_r2 = sum_r2 + v * v
      if v > 0.0 {
        n_expr_r = n_expr_r + 1
      }
    }
    let mean_r = if n_rest > 0 { sum_r / n_rest.to_double() } else { 0.0 }
    let var_r = if n_rest > 1 {
      (sum_r2 - sum_r * sum_r / n_rest.to_double()) / (n_rest - 1).to_double()
    } else {
      0.0
    }
    // Overestimated variance: add mean to variance
    let var_g_oe = var_g + mean_g
    let var_r_oe = var_r + mean_r
    // t-statistic
    let se = @math.pow(
      var_g_oe / n_group.to_double() + var_r_oe / n_rest.to_double(),
      0.5,
    )
    let t = if se > 0.0 { (mean_g - mean_r) / se } else { 0.0 }
    // Log fold change
    let lfc = if mean_r > 0.0 && mean_g > 0.0 {
      @math.log2(mean_g) - @math.log2(mean_r)
    } else if mean_g > 0.0 {
      @math.log2(mean_g)
    } else {
      0.0
    }
    // Approximate p-value from t (simplified normal approximation)
    let pval = normal_pvalue_two_sided(t)
    gene_indices.push(g)
    scores.push(t)
    logfoldchanges.push(lfc)
    pvals.push(pval)
    pvals_adj.push(pval * n_genes.to_double()) // Bonferroni
    pts_group.push(
      if n_group > 0 {
        n_expr_g.to_double() / n_group.to_double()
      } else {
        0.0
      },
    )
    pts_rest.push(
      if n_rest > 0 {
        n_expr_r.to_double() / n_rest.to_double()
      } else {
        0.0
      },
    )
  }
  RankGenesExtResult::{
    gene_indices,
    scores,
    logfoldchanges,
    pvals,
    pvals_adj,
    pts_group,
    pts_rest,
    method_: "t-test_overestim_var",
  }
}

///|
/// Filter rank genes by minimum fold change and adjusted p-value.
pub fn scanpy_filter_rank_genes_by_fc(
  result : RankGenesExtResult,
  min_fold_change~ : Double = 1.0,
  max_pval_adj~ : Double = 0.05,
) -> Array[Int32] {
  let filtered = Array(capacity=0)
  for i in 0..<result.gene_indices.length() {
    let fc = if result.logfoldchanges[i] < 0.0 {
      -result.logfoldchanges[i]
    } else {
      result.logfoldchanges[i]
    }
    if fc >= min_fold_change && result.pvals_adj[i] <= max_pval_adj {
      filtered.push(result.gene_indices[i])
    }
  }
  filtered
}

///|
/// Prepare dotplot data from rank genes result.
/// Returns (gene_index, mean_expr_group, pct_expr_group, mean_expr_rest, pct_expr_rest)
pub fn scanpy_rank_genes_dotplot_data(
  result : RankGenesExtResult,
) -> Array[(Int32, Double, Double, Double, Double)] {
  let data = Array(capacity=result.gene_indices.length())
  for i in 0..<result.gene_indices.length() {
    // Use scores as proxy for mean expression difference
    let mean_group = if result.pts_group[i] > 0.0 {
      result.scores[i] * result.pts_group[i]
    } else {
      0.0
    }
    let mean_rest = if result.pts_rest[i] > 0.0 {
      result.scores[i] * result.pts_rest[i]
    } else {
      0.0
    }
    data.push(
      (
        result.gene_indices[i],
        mean_group,
        result.pts_group[i],
        mean_rest,
        result.pts_rest[i],
      ),
    )
  }
  data
}

///|
/// Rank genes with tie correction (Wilcoxon rank-sum with tie handling).
pub fn scanpy_rank_genes_wilcoxon_tie_correct(
  x_group : Array[Array[Double]],
  x_rest : Array[Array[Double]],
) -> Array[(Int32, Double)] {
  let n_genes = if x_group.length() > 0 { x_group[0].length() } else { 0 }
  let n_group = x_group.length()
  let n_rest = x_rest.length()
  let n_total = n_group + n_rest
  let result = Array(capacity=n_genes)
  for g in 0..<n_genes {
    // Combine all values with group labels
    let combined = Array(capacity=n_total)
    for i in 0..<n_group {
      combined.push((x_group[i][g], 0))
    }
    for i in 0..<n_rest {
      combined.push((x_rest[i][g], 1))
    }
    // Sort by value
    combined.sort_by(fn(a, b) {
      if a.0 < b.0 {
        -1
      } else if a.0 > b.0 {
        1
      } else {
        0
      }
    })
    // Assign ranks with tie correction
    let ranks = Array(capacity=n_total)
    let mut i = 0
    while i < n_total {
      let mut j = i + 1
      while j < n_total && combined[j].0 == combined[i].0 {
        j = j + 1
      }
      // Ties from i to j-1, assign average rank
      let avg_rank = (i + j + 1).to_double() / 2.0
      for _k in i..<j {
        ranks.push(avg_rank)
      }
      i = j
    }
    // Sum of ranks for group 0
    let mut sum_ranks_g = 0.0
    for k in 0..<n_total {
      if combined[k].1 == 0 {
        sum_ranks_g = sum_ranks_g + ranks[k]
      }
    }
    // Mann-Whitney U
    let u = sum_ranks_g - n_group.to_double() * (n_group + 1).to_double() / 2.0
    // Z-score with tie correction
    let mut tie_term = 0.0
    let mut ti = 0
    let mut k = 0
    while k < n_total {
      let mut j = k + 1
      while j < n_total && combined[j].0 == combined[k].0 {
        j = j + 1
      }
      ti = j - k
      tie_term = tie_term + ti.to_double() * (ti * ti - 1).to_double()
      k = j
    }
    let n_f = n_total.to_double()
    let ng_f = n_group.to_double()
    let nr_f = n_rest.to_double()
    let sigma = @math.pow(
      ng_f * nr_f / 12.0 * (n_f + 1.0 - tie_term / (n_f * (n_f - 1.0))),
      0.5,
    )
    let z = if sigma > 0.0 { (u - ng_f * nr_f / 2.0) / sigma } else { 0.0 }
    result.push((g, z))
  }
  result
}

///|
/// Helper: two-sided p-value from normal approximation.
fn normal_pvalue_two_sided(z : Double) -> Double {
  // Simplified: use erfc approximation
  // p = 2 * (1 - Phi(|z|)) ≈ 2 * 0.5 * erfc(|z|/sqrt(2))
  let abs_z = if z < 0.0 { -z } else { z }
  // erfc approximation (Abramowitz & Stegun 7.1.26)
  let t = 1.0 / (1.0 + 0.3275911 * abs_z / @math.pow(2.0, 0.5))
  let z2 = abs_z / @math.pow(2.0, 0.5)
  let erfc_val = t *
    (
      0.254829592 +
      t *
      (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))
    ) *
    @math.pow(2.718281828459045, -z2 * z2)
  let p = erfc_val
  if p > 1.0 {
    1.0
  } else if p < 0.0 {
    0.0
  } else {
    p
  }
}
