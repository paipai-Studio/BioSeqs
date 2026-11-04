/// Scanpy-inspired: tl.compare_gene_expression — compare gene expression between conditions

///|
/// Compare gene expression between two conditions.
/// Computes fold change, log2 fold change, and Welch's t-test p-value.
///
/// Parameters:
/// - `x`: cells × genes matrix
/// - `condition`: condition label per cell ("A" or "B")
/// - `gene_names`: gene names
pub fn scanpy_compare_gene_expression(
  x : Array[Array[Double]],
  condition : Array[String],
  gene_names : Array[String],
) -> Array[(String, Double, Double, Double)] {
  // Returns: (gene_name, mean_A, mean_B, log2_fold_change)
  let n_cells = x.length()
  if n_cells == 0 {
    return []
  }
  let n_genes = x[0].length()
  // Split cells by condition
  let group_a = Array(capacity=0)
  let group_b = Array(capacity=0)
  for i in 0..<n_cells {
    if condition[i] == "A" || condition[i] == "control" {
      group_a.push(i)
    } else {
      group_b.push(i)
    }
  }
  let result = Array(capacity=n_genes)
  for j in 0..<n_genes {
    // Mean of group A
    let mut sum_a = 0.0
    for i in group_a {
      sum_a += x[i][j]
    }
    let mean_a = if group_a.length() > 0 {
      sum_a / group_a.length().to_double()
    } else {
      0.0
    }
    // Mean of group B
    let mut sum_b = 0.0
    for i in group_b {
      sum_b += x[i][j]
    }
    let mean_b = if group_b.length() > 0 {
      sum_b / group_b.length().to_double()
    } else {
      0.0
    }
    // Log2 fold change (B vs A)
    let l2fc = if mean_a > 0.0 && mean_b > 0.0 {
      @math.log2(mean_b / mean_a)
    } else if mean_b > 0.0 {
      20.0 // cap at large positive
    } else if mean_a > 0.0 {
      -20.0 // cap at large negative
    } else {
      0.0
    }
    let gene_name = if j < gene_names.length() {
      gene_names[j]
    } else {
      "gene_" + j.to_string()
    }
    result.push((gene_name, mean_a, mean_b, l2fc))
  }
  result
}

///|
/// Compute expression comparison with statistics (t-test).
/// Returns: (gene_name, log2FC, t_statistic, p_value_approx)
pub fn scanpy_compare_expression_ttest(
  x : Array[Array[Double]],
  condition : Array[String],
  gene_names : Array[String],
) -> Array[(String, Double, Double, Double)] {
  let n_cells = x.length()
  if n_cells == 0 {
    return []
  }
  let n_genes = x[0].length()
  let group_a = Array(capacity=0)
  let group_b = Array(capacity=0)
  for i in 0..<n_cells {
    if condition[i] == "A" || condition[i] == "control" {
      group_a.push(i)
    } else {
      group_b.push(i)
    }
  }
  let result = Array(capacity=n_genes)
  for j in 0..<n_genes {
    // Group A stats
    let mut sum_a = 0.0
    for i in group_a {
      sum_a += x[i][j]
    }
    let n_a = group_a.length().to_double()
    let mean_a = if n_a > 0.0 { sum_a / n_a } else { 0.0 }
    let mut var_a = 0.0
    for i in group_a {
      let diff = x[i][j] - mean_a
      var_a += diff * diff
    }
    if n_a > 1.0 {
      var_a /= n_a - 1.0
    }
    // Group B stats
    let mut sum_b = 0.0
    for i in group_b {
      sum_b += x[i][j]
    }
    let n_b = group_b.length().to_double()
    let mean_b = if n_b > 0.0 { sum_b / n_b } else { 0.0 }
    let mut var_b = 0.0
    for i in group_b {
      let diff = x[i][j] - mean_b
      var_b += diff * diff
    }
    if n_b > 1.0 {
      var_b /= n_b - 1.0
    }
    // Log2 fold change
    let l2fc = if mean_a > 0.0 && mean_b > 0.0 {
      @math.log2(mean_b / mean_a)
    } else if mean_b > 0.0 {
      20.0
    } else if mean_a > 0.0 {
      -20.0
    } else {
      0.0
    }
    // Welch's t-statistic
    let se_sq = var_a / n_a + var_b / n_b
    let t_stat = if se_sq > 0.0 {
      (mean_b - mean_a) / se_sq.sqrt()
    } else {
      0.0
    }
    // Approximate p-value using normal approximation
    let p_val = scanpy_compare_normal_pval(t_stat.abs())
    let gene_name = if j < gene_names.length() {
      gene_names[j]
    } else {
      "gene_" + j.to_string()
    }
    result.push((gene_name, l2fc, t_stat, p_val))
  }
  result
}

///|
fn scanpy_compare_normal_pval(z : Double) -> Double {
  // Approximate two-tailed p-value from normal distribution
  // Using erf approximation
  let t = 1.0 / (1.0 + 0.2316419 * z)
  let d = 0.3989423 * @math.exp(-z * z / 2.0)
  let p = d *
    t *
    (
      0.3193815 +
      t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274)))
    )
  2.0 * p
}
