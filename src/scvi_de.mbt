// scvi_de.mbt - Bayesian differential expression
// Inspired by scvi-tools differential_expression method
// Implements vanilla and change-mode DE with Bayes factors

///|
/// Compute normalized expression per group
pub fn scvi_de_normalized(
  x : Array[Array[Double]],
  group_labels : Array[Int32],
  n_groups : Int32,
) -> Array[Array[Double]] {
  let n_genes = x[0].length()
  let group_means = Array(capacity=n_groups)
    |> a => {
      for _g in 0..<n_groups {
        a.push(
          Array(capacity=n_genes)
          |> ga => {
            for _i in 0..<n_genes {
              ga.push(0.0)
            }
            ga
          },
        )
      }
      a
    }
  let group_counts = Array(capacity=n_groups)
    |> a => {
      for _i in 0..<n_groups {
        a.push(0)
      }
      a
    }
  for i in 0..<x.length() {
    let g = group_labels[i]
    group_counts[g] += 1
    for j in 0..<n_genes {
      group_means[g][j] += x[i][j]
    }
  }
  for g in 0..<n_groups {
    if group_counts[g] > 0 {
      for j in 0..<n_genes {
        group_means[g][j] /= group_counts[g].to_double()
      }
    }
  }
  group_means
}

///|
/// Bayes factor for differential expression (change mode)
pub fn scvi_de_bayes_factor(
  expr1 : Array[Double],
  expr2 : Array[Double],
  delta : Double,
) -> Array[Double] {
  let n_genes = expr1.length()
  let result = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let lfc = @math.ln(expr1[g] + 0.0000000001) -
      @math.ln(expr2[g] + 0.0000000001)
    // Bayes factor: prob(|lfc| > delta) / prob(|lfc| < delta)
    let prob_change = if lfc.abs() > delta { 1.0 } else { 0.0 }
    let prob_null = 1.0 - prob_change
    result.push(
      @math.ln(prob_change + 0.0000000001) - @math.ln(prob_null + 0.0000000001),
    )
  }
  result
}

///|
/// Differential expression analysis between two groups
pub fn scvi_differential_expression(
  x : Array[Array[Double]],
  group_labels : Array[Int32],
  group1 : Int32,
  group2 : Int32,
  mode : String,
  delta : Double,
) -> Array[(String, Double, Double, Double)] {
  let n_genes = x[0].length()
  // Get group means
  let group_means = scvi_de_normalized(x, group_labels, group_labels.length())
  let expr1 = group_means[group1]
  let expr2 = group_means[group2]
  let result = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let lfc = @math.ln(expr1[g] + 0.0000000001) -
      @math.ln(expr2[g] + 0.0000000001)
    let bayes_factor = if mode == "change" {
      let bf = scvi_de_bayes_factor(expr1, expr2, delta)
      bf[g]
    } else {
      // Vanilla mode: use t-statistic-like score
      lfc
    }
    let mean_diff = expr1[g] - expr2[g]
    result.push(("gene_{g}".to_string(), lfc, bayes_factor, mean_diff))
  }
  // Sort by absolute Bayes factor
  result.sort_by(fn(a, b) {
    let abs_a = a.2.abs()
    let abs_b = b.2.abs()
    if abs_a > abs_b {
      -1
    } else if abs_a < abs_b {
      1
    } else {
      0
    }
  })
  result
}

///|
/// Get differentially expressed genes with FDR correction (Benjamini-Hochberg)
pub fn scvi_de_fdr_correction(
  p_values : Array[Double],
  fdr_target : Double,
) -> Array[Bool] {
  let n = p_values.length()
  // Create index-value pairs and sort by p-value
  let indexed = Array(capacity=n)
    |> a => {
      for i in 0..<n {
        a.push((i, p_values[i]))
      }
      a
    }
  indexed.sort_by(fn(a, b) {
    if a.1 < b.1 {
      -1
    } else if a.1 > b.1 {
      1
    } else {
      0
    }
  })
  let significant = Array(capacity=n)
    |> a => {
      for _i in 0..<n {
        a.push(false)
      }
      a
    }
  // BH procedure
  let mut max_rank = 0
  for rank in 1..<=n {
    let threshold = fdr_target * rank.to_double() / n.to_double()
    if indexed[rank - 1].1 <= threshold {
      max_rank = rank
    }
  }
  for rank in 0..<max_rank {
    significant[indexed[rank].0] = true
  }
  significant
}

///|
/// Compute posterior predictive p-values for DE
pub fn scvi_de_posterior_pvalue(
  expr1_samples : Array[Double],
  expr2_samples : Array[Double],
) -> Double {
  let n = expr1_samples.length()
  let mut count = 0
  for i in 0..<n {
    if expr1_samples[i] > expr2_samples[i] {
      count += 1
    }
  }
  count.to_double() / n.to_double()
}

///|
/// Multiple group comparison (one-vs-rest)
pub fn scvi_de_one_vs_rest(
  x : Array[Array[Double]],
  group_labels : Array[Int32],
  n_groups : Int32,
  delta : Double,
) -> Array[Array[(String, Double, Double, Double)]] {
  let result = Array(capacity=n_groups)
  for g in 0..<n_groups {
    // Create binary labels: group g vs rest
    let binary_labels = Array(capacity=group_labels.length())
    for i in 0..<group_labels.length() {
      if group_labels[i] == g {
        binary_labels.push(g)
      } else {
        binary_labels.push(n_groups)
      }
    }
    let de_results = scvi_differential_expression(
      x, binary_labels, g, n_groups, "change", delta,
    )
    result.push(de_results)
  }
  result
}
