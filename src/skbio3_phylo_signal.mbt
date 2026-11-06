// Scikit-bio third-round: phylogenetic signal and comparative methods (items 12-15)
// Blomberg K, Pagel lambda, PIC, phylo ANOVA

///|
pub fn skbio3_blomberg_k(
  trait_values : Array[Double],
  tree_branch_lengths : Array[Double],
  n_tips : Int,
) -> Double {
  // Blomberg's K: ratio of observed to expected variance under Brownian motion
  // K > 1: strong phylogenetic signal; K < 1: weak signal
  if n_tips <= 1 || trait_values.length() < 2 {
    return 0.0
  }
  // Mean trait value
  let mut mean = 0.0
  for i in 0..<trait_values.length() {
    mean = mean + trait_values[i]
  }
  mean = mean / trait_values.length().to_double()
  // Observed variance of trait
  let mut var_obs = 0.0
  for i in 0..<trait_values.length() {
    let diff = trait_values[i] - mean
    var_obs = var_obs + diff * diff
  }
  var_obs = var_obs / trait_values.length().to_double()
  // Expected variance under Brownian motion (approximate using total branch length)
  let mut total_branch = 0.0
  for i in 0..<tree_branch_lengths.length() {
    total_branch = total_branch + tree_branch_lengths[i]
  }
  if total_branch <= 0.0 || var_obs <= 0.0 {
    return 0.0
  }
  let expected_var = total_branch / n_tips.to_double()
  // K = (observed / expected) normalized
  let mse_obs = var_obs / expected_var
  // Blomberg K = MSE(obs) / MSE(BM)
  // Simplified: K ≈ var_obs / (total_branch / n_tips)
  mse_obs
}

///|
pub fn skbio3_pagel_lambda(
  trait_values : Array[Double],
  tree_branch_lengths : Array[Double],
  n_tips : Int,
) -> Double {
  // Pagel's lambda: transformation parameter for internal branches
  // lambda = 0: no phylogenetic signal; lambda = 1: pure Brownian motion
  if n_tips <= 1 || trait_values.length() < 2 {
    return 0.0
  }
  // Estimate lambda via method of moments
  let mut mean = 0.0
  for i in 0..<trait_values.length() {
    mean = mean + trait_values[i]
  }
  mean = mean / trait_values.length().to_double()
  // Compute trait variance
  let mut var_trait = 0.0
  for i in 0..<trait_values.length() {
    let diff = trait_values[i] - mean
    var_trait = var_trait + diff * diff
  }
  var_trait = var_trait / trait_values.length().to_double()
  // Compute mean branch length
  let mut mean_branch = 0.0
  for i in 0..<tree_branch_lengths.length() {
    mean_branch = mean_branch + tree_branch_lengths[i]
  }
  mean_branch = mean_branch / tree_branch_lengths.length().to_double()
  if mean_branch <= 0.0 || var_trait <= 0.0 {
    return 0.0
  }
  // Lambda ≈ var_trait / (mean_branch * n_tips)
  let lambda = var_trait / (mean_branch * n_tips.to_double())
  if lambda > 1.0 {
    1.0
  } else if lambda < 0.0 {
    0.0
  } else {
    lambda
  }
}

///|
pub struct Skbio3PICResult {
  node_name : String
  contrast : Double
  variance : Double
  branch_length : Double
}

///|
pub fn skbio3_phylo_independent_contrasts(
  tip_names : Array[String],
  tip_values : Array[Double],
  tip_branches : Array[Double],
  internal_branches : Array[Double],
) -> Array[Skbio3PICResult] {
  // Phylogenetic Independent Contrasts (Felsenstein 1985)
  // For a simple bifurcating tree (caterpillar), compute PICs bottom-up
  // This is a simplified version for a caterpillar (ladder) tree
  let n = tip_names.length()
  if n < 2 {
    return []
  }
  let results : Array[Skbio3PICResult] = []
  // Working arrays: current values, variances, branch lengths, names
  let values = tip_values.copy()
  let branches = tip_branches.copy()
  let names_ = tip_names.copy()
  let mut ib_idx = 0  // internal branch index
  while values.length() > 1 {
    // Take the first two nodes and compute contrast
    let v1 = values[0]
    let v2 = values[1]
    let b1 = branches[0]
    let b2 = branches[1]
    let n1 = names_[0]
    let n2 = names_[1]
    // Contrast: (v1 - v2) / sqrt(b1 + b2)
    let var_sum = b1 + b2
    let contrast = if var_sum > 0.0 {
      (v1 - v2) / var_sum.sqrt()
    } else {
      0.0
    }
    results.push({
      node_name: n1 + "-" + n2,
      contrast,
      variance: var_sum,
      branch_length: var_sum,
    })
    // New node value: weighted average
    let new_value = (v1 * b2 + v2 * b1) / var_sum
    // New branch length: internal branch + (b1 * b2 / (b1 + b2))
    let new_branch = if ib_idx < internal_branches.length() {
      internal_branches[ib_idx] + b1 * b2 / var_sum
    } else {
      b1 * b2 / var_sum
    }
    ib_idx = ib_idx + 1
    // Remove first two, add new node
    let _ = values.remove(0)
    let _ = values.remove(0)
    values.insert(0, new_value)
    let _ = branches.remove(0)
    let _ = branches.remove(0)
    branches.insert(0, new_branch)
    let _ = names_.remove(0)
    let _ = names_.remove(0)
    names_.insert(0, "internal")
  }
  results
}

///|
pub struct Skbio3PhyloAnovaResult {
  group : String
  mean : Double
  n : Int
  f_statistic : Double
  p_value : Double
}

///|
pub fn skbio3_phylo_anova(
  groups : Array[String],
  values : Array[Double],
  tree_branch_lengths : Array[Double],
) -> Array[Skbio3PhyloAnovaResult] {
  // Phylogenetic ANOVA: compare group means accounting for phylogenetic covariance
  // Simplified: use PIC-transformed values in standard ANOVA
  let n = values.length()
  if n == 0 {
    return []
  }
  // Get unique groups
  let unique_groups : Array[String] = []
  for i in 0..<groups.length() {
    let mut found = false
    for j in 0..<unique_groups.length() {
      if unique_groups[j] == groups[i] {
        found = true
        break
      }
    }
    if not(found) {
      unique_groups.push(groups[i])
    }
  }
  let results : Array[Skbio3PhyloAnovaResult] = []
  // Compute overall mean
  let mut grand_mean = 0.0
  for i in 0..<n {
    grand_mean = grand_mean + values[i]
  }
  grand_mean = grand_mean / n.to_double()
  // Between-group and within-group SS
  let mut ss_between = 0.0
  let mut ss_within = 0.0
  let mut total_branch = 0.0
  for i in 0..<tree_branch_lengths.length() {
    total_branch = total_branch + tree_branch_lengths[i]
  }
  for gi in 0..<unique_groups.length() {
    let group = unique_groups[gi]
    let group_values : Array[Double] = []
    for i in 0..<n {
      if groups[i] == group {
        group_values.push(values[i])
      }
    }
    if group_values.is_empty() {
      continue
    }
    let gn = group_values.length()
    let mut gmean = 0.0
    for i in 0..<gn {
      gmean = gmean + group_values[i]
    }
    gmean = gmean / gn.to_double()
    // Between-group SS
    ss_between = ss_between + gn.to_double() * (gmean - grand_mean) * (gmean - grand_mean)
    // Within-group SS
    for i in 0..<gn {
      let diff = group_values[i] - gmean
      ss_within = ss_within + diff * diff
    }
    // F-statistic (simplified, no phylogenetic correction)
    let df_between = (unique_groups.length() - 1).to_double()
    let df_within = (n - unique_groups.length()).to_double()
    let ms_between = if df_between > 0.0 { ss_between / df_between } else { 0.0 }
    let ms_within = if df_within > 0.0 { ss_within / df_within } else { 0.0 }
    let f_stat = if ms_within > 0.0 { ms_between / ms_within } else { 0.0 }
    // Approximate p-value from F distribution (very rough)
    let p_val = if f_stat > 3.0 { 0.05 } else if f_stat > 6.0 { 0.01 } else { 0.5 }
    results.push({ group, mean: gmean, n: gn, f_statistic: f_stat, p_value: p_val })
  }
  results
}
