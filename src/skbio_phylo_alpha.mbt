/// Phylogenetic alpha diversity metrics.
///
/// scikit-bio equivalent: `diversity.alpha` phylogenetic metrics
///
/// Beyond Faith's PD, these metrics incorporate phylogenetic structure:
/// - PE: Phylogenetic Entropy (Shannon-like with branch lengths)
/// - PSV: Phylogenetic Species Variability
/// - PSE: Phylogenetic Species Evenness
/// - PSR: Phylogenetic Species Richness

///|
/// Phylogenetic entropy (PE).
/// PE = -sum(p_i * ln(p_i) * branch_length_i)
/// where p_i is the relative abundance of species i.
///
/// Requires: species abundances and corresponding branch lengths
/// (distance from root to each species tip).
pub fn skbio_phylo_entropy(
  abundances : Array[Double],
  branch_lengths : Array[Double],
) -> Double {
  let n = abundances.length()
  if n == 0 || n != branch_lengths.length() {
    return 0.0
  }
  let mut total = 0.0
  for a in abundances {
    total = total + a
  }
  if total <= 0.0 {
    return 0.0
  }
  let mut pe = 0.0
  for i = 0; i < n; i = i + 1 {
    let p = abundances[i] / total
    if p > 0.0 {
      pe = pe - p * @math.ln(p) * branch_lengths[i]
    }
  }
  pe
}

///|
/// Phylogenetic Species Variability (PSV).
/// PSV = 1 - (sum(sum(p_i * p_j * d_ij)) / mean_branch_length)
/// where d_ij is the phylogenetic distance between species i and j.
/// PSV ranges from 0 (species are phylogenetically identical) to 1
/// (species are maximally distantly related).
///
/// Requires: relative abundances and pairwise phylogenetic distance matrix.
pub fn skbio_phylo_psv(
  abundances : Array[Double],
  dist_matrix : Array[Array[Double]],
) -> Double {
  let n = abundances.length()
  if n == 0 || n != dist_matrix.length() {
    return 0.0
  }
  let mut total = 0.0
  for a in abundances {
    total = total + a
  }
  if total <= 0.0 {
    return 0.0
  }
  // Compute weighted average pairwise distance
  let mut weighted_dist = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      let pi = abundances[i] / total
      let pj = abundances[j] / total
      weighted_dist = weighted_dist + pi * pj * dist_matrix[i][j]
    }
  }
  // PSV = 1 - weighted_dist / mean_dist
  // Mean distance of the tree (sum of branch lengths / number of tips)
  let mut mean_dist = 0.0
  let mut count = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      mean_dist = mean_dist + dist_matrix[i][j]
      count = count + 1.0
    }
  }
  if count > 0.0 {
    mean_dist = mean_dist / count
  }
  if mean_dist > 0.0 {
    1.0 - weighted_dist / (2.0 * mean_dist)
  } else {
    1.0
  }
}

///|
/// Phylogenetic Species Evenness (PSE).
/// PSE = PSV * (species_evenness)
/// Combines phylogenetic structure with species evenness.
pub fn skbio_phylo_pse(
  abundances : Array[Double],
  dist_matrix : Array[Array[Double]],
) -> Double {
  let n = abundances.length()
  if n == 0 {
    return 0.0
  }
  let psv = skbio_phylo_psv(abundances, dist_matrix)
  // Compute species evenness (Simpson's evenness)
  let mut total = 0.0
  for a in abundances {
    total = total + a
  }
  if total <= 0.0 {
    return 0.0
  }
  let mut sum_sq = 0.0
  for a in abundances {
    let p = a / total
    sum_sq = sum_sq + p * p
  }
  let simpson_e = if sum_sq > 0.0 { 1.0 / sum_sq / n.to_double() } else { 0.0 }
  psv * simpson_e
}

///|
/// Phylogenetic Species Richness (PSR).
/// PSR = PSV * species_richness
pub fn skbio_phylo_psr(
  abundances : Array[Double],
  dist_matrix : Array[Array[Double]],
) -> Double {
  let psv = skbio_phylo_psv(abundances, dist_matrix)
  let mut richness = 0
  for a in abundances {
    if a > 0.0 {
      richness = richness + 1
    }
  }
  psv * richness.to_double()
}

///|
/// Faith's Phylogenetic Diversity (PD) - sum of branch lengths for present species.
/// This is a re-implementation for standalone use (also exists in skbio_tree.mbt).
pub fn skbio_phylo_pd(
  present : Array[Bool],
  branch_lengths : Array[Double],
) -> Double {
  let n = present.length()
  if n != branch_lengths.length() {
    return 0.0
  }
  let mut pd = 0.0
  for i = 0; i < n; i = i + 1 {
    if present[i] {
      pd = pd + branch_lengths[i]
    }
  }
  pd
}

///|
/// Phylogenetic Diversity Variance (PDV).
/// Measures the variance of pairwise phylogenetic distances among present species.
pub fn skbio_phylo_pdv(
  abundances : Array[Double],
  dist_matrix : Array[Array[Double]],
) -> Double {
  let n = abundances.length()
  if n == 0 {
    return 0.0
  }
  let mut total = 0.0
  for a in abundances {
    total = total + a
  }
  if total <= 0.0 {
    return 0.0
  }
  // Collect weighted pairwise distances
  let dists : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let pi = abundances[i] / total
      let pj = abundances[j] / total
      if pi > 0.0 && pj > 0.0 {
        dists.push(dist_matrix[i][j])
      }
    }
  }
  if dists.length() == 0 {
    return 0.0
  }
  let mut mean = 0.0
  for d in dists {
    mean = mean + d
  }
  mean = mean / dists.length().to_double()
  let mut variance = 0.0
  for d in dists {
    variance = variance + (d - mean) * (d - mean)
  }
  variance / dists.length().to_double()
}
