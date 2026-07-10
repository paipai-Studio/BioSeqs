///|
/// scikit-bio Diversity module - Alpha and beta diversity metrics.

///|
/// Helper: sum of all counts.
fn sum_counts(counts : Array[Double]) -> Double {
  let mut total = 0.0
  for c in counts {
    total = total + c
  }
  total
}

///|
/// Helper: number of observed features (non-zero counts).
fn observed_count(counts : Array[Double]) -> Int {
  let mut n = 0
  for c in counts {
    if c > 0.0 {
      n = n + 1
    }
  }
  n
}

///|
/// Number of species observed (richness).
pub fn observed_otus(counts : Array[Double]) -> Int {
  observed_count(counts)
}

///|
/// Alias for observed_otus.
pub fn sobs(counts : Array[Double]) -> Int {
  observed_count(counts)
}

///|
/// Alias for observed_otus.
pub fn observed_features(counts : Array[Double]) -> Int {
  observed_count(counts)
}

///|
/// Shannon-Wiener diversity index.
/// H = -sum(p_i * @math.ln(p_i))
pub fn shannon(counts : Array[Double]) -> Double {
  let total = sum_counts(counts)
  if total == 0.0 {
    return 0.0
  }
  let mut h = 0.0
  for c in counts {
    if c > 0.0 {
      let p = c / total
      h = h - p * @math.ln(p)
    }
  }
  h
}

///|
/// Shannon-Wiener diversity index with custom base.
pub fn shannon_base(counts : Array[Double], base : Double) -> Double {
  let h = shannon(counts)
  if base == 2.718281828459045 {
    h
  } else {
    h / @math.ln(base)
  }
}

///|
/// Simpson's diversity index.
/// D = 1 - sum(p_i^2)
pub fn simpson(counts : Array[Double]) -> Double {
  let total = sum_counts(counts)
  if total == 0.0 {
    return 0.0
  }
  let mut sum_p2 = 0.0
  for c in counts {
    let p = c / total
    sum_p2 = sum_p2 + p * p
  }
  1.0 - sum_p2
}

///|
/// Simpson's dominance index.
/// D = sum(p_i^2)
pub fn simpson_d(counts : Array[Double]) -> Double {
  let total = sum_counts(counts)
  if total == 0.0 {
    return 0.0
  }
  let mut sum_p2 = 0.0
  for c in counts {
    let p = c / total
    sum_p2 = sum_p2 + p * p
  }
  sum_p2
}

///|
/// Inverse Simpson index.
/// D = 1 / sum(p_i^2)
pub fn inv_simpson(counts : Array[Double]) -> Double {
  let total = sum_counts(counts)
  if total == 0.0 {
    return 0.0
  }
  let mut sum_p2 = 0.0
  for c in counts {
    let p = c / total
    sum_p2 = sum_p2 + p * p
  }
  if sum_p2 == 0.0 {
    0.0
  } else {
    1.0 / sum_p2
  }
}

///|
/// Dominance index (same as simpson_d).
pub fn dominance(counts : Array[Double]) -> Double {
  simpson_d(counts)
}

///|
/// Gini-Simpson index (same as simpson).
pub fn gini_simpson(counts : Array[Double]) -> Double {
  simpson(counts)
}

///|
/// Berger-Parker dominance index.
/// d = N_max / N_total
pub fn berger_parker_d(counts : Array[Double]) -> Double {
  let total = sum_counts(counts)
  if total == 0.0 {
    return 0.0
  }
  let mut max_c = 0.0
  for c in counts {
    if c > max_c {
      max_c = c
    }
  }
  max_c / total
}

///|
/// Alias for berger_parker_d.
pub fn berger_parker(counts : Array[Double]) -> Double {
  berger_parker_d(counts)
}

///|
/// Pielou's evenness index.
/// J = H' / @math.ln(S)
pub fn pielou_e(counts : Array[Double]) -> Double {
  let s = observed_count(counts).to_double()
  if s <= 1.0 {
    return 0.0
  }
  let h = shannon(counts)
  h / @math.ln(s)
}

///|
/// Simpson's evenness index.
/// E_D = (1/D) / S
pub fn simpson_e(counts : Array[Double]) -> Double {
  let s = observed_count(counts).to_double()
  if s == 0.0 {
    return 0.0
  }
  let inv_d = inv_simpson(counts)
  inv_d / s
}

///|
/// Heip's evenness index.
/// E_Heip = (e^H' - 1) / (S - 1)
pub fn heip_e(counts : Array[Double]) -> Double {
  let s = observed_count(counts).to_double()
  if s <= 1.0 {
    return 1.0
  }
  let h = shannon(counts)
  (@math.exp(h) - 1.0) / (s - 1.0)
}

///|
/// Chao1 richness estimator (bias-corrected, matches scikit-bio default).
/// S_chao1 = S_obs + F1*(F1-1) / (2*(F2+1))
pub fn chao1(counts : Array[Double]) -> Double {
  let s_obs = observed_count(counts).to_double()
  let mut f1 = 0.0
  let mut f2 = 0.0
  for c in counts {
    if c == 1.0 {
      f1 = f1 + 1.0
    } else if c == 2.0 {
      f2 = f2 + 1.0
    }
  }
  // Bias-corrected Chao1 (matches scikit-bio default)
  s_obs + f1 * (f1 - 1.0) / (2.0 * (f2 + 1.0))
}

///|
/// Classic Chao1 richness estimator.
/// S_chao1 = S_obs + F1^2 / (2*F2)
pub fn chao1_classic(counts : Array[Double]) -> Double {
  let s_obs = observed_count(counts).to_double()
  let mut f1 = 0.0
  let mut f2 = 0.0
  for c in counts {
    if c == 1.0 {
      f1 = f1 + 1.0
    } else if c == 2.0 {
      f2 = f2 + 1.0
    }
  }
  if f2 > 0.0 {
    s_obs + f1 * f1 / (2.0 * f2)
  } else {
    s_obs + f1 * (f1 - 1.0) / (2.0 * (f2 + 1.0))
  }
}

///|
/// ACE richness estimator.
/// Merges two traversals into one by accumulating sum_c2 in the first pass,
/// then deriving sum_p2 = sum_c2 / n_rare^2 without a second pass.
pub fn ace(counts : Array[Double], rare_threshold : Int) -> Double {
  let mut s_rare = 0.0
  let mut s_abund = 0.0
  let mut n_rare = 0.0
  let mut f1 = 0.0
  let mut sum_c2 = 0.0
  let threshold = rare_threshold.to_double()

  for c in counts {
    if c > 0.0 {
      if c <= threshold {
        s_rare = s_rare + 1.0
        n_rare = n_rare + c
        sum_c2 = sum_c2 + c * c
        if c == 1.0 {
          f1 = f1 + 1.0
        }
      } else {
        s_abund = s_abund + 1.0
      }
    }
  }

  if n_rare == 0.0 {
    return s_abund
  }

  let c_ace = 1.0 - f1 / n_rare

  if c_ace == 0.0 {
    return s_rare + s_abund
  }

  // sum_p2 = sum((c/n_rare)^2) = sum_c2 / n_rare^2
  let sum_p2 = sum_c2 / (n_rare * n_rare)

  let gamma = if c_ace > 0.0 {
    let s_rare_val = s_rare / c_ace
    s_rare_val * sum_p2 / c_ace - 1.0
  } else {
    0.0
  }

  s_abund + s_rare / c_ace + f1 / c_ace * gamma.max(0.0)
}

///|
/// ACE richness estimator with default threshold of 10.
pub fn ace_default(counts : Array[Double]) -> Double {
  ace(counts, 10)
}

///|
/// Fisher's alpha diversity index.
/// Solves S = alpha * @math.ln(1 + N/alpha) iteratively.
pub fn fisher_alpha(counts : Array[Double]) -> Double {
  let s = observed_count(counts).to_double()
  let n = sum_counts(counts)
  if n == 0.0 || s == 0.0 {
    return 0.0
  }
  // Iterative solution using Newton's method
  let mut alpha = s // initial guess
  let mut i = 0
  while i < 50 {
    let x = n / (n + alpha)
    let f = alpha * -@math.ln(1.0 - x) - s
    let df = -@math.ln(1.0 - x) - n / (n + alpha)
    if df == 0.0 {
      break
    }
    let step = f / df
    alpha = alpha - step
    if alpha <= 0.0 {
      alpha = 0.001
    }
    if step.abs() < 0.0001 {
      break
    }
    i = i + 1
  }
  alpha
}

///|
/// Good's coverage index.
/// C = 1 - (F1 / N)
pub fn goods_coverage(counts : Array[Double]) -> Double {
  let n = sum_counts(counts)
  if n == 0.0 {
    return 0.0
  }
  let mut f1 = 0.0
  for c in counts {
    if c == 1.0 {
      f1 = f1 + 1.0
    }
  }
  1.0 - f1 / n
}

///|
/// Margalef's richness index.
/// D_Mg = (S - 1) / @math.ln(N)
pub fn margalef(counts : Array[Double]) -> Double {
  let s = observed_count(counts).to_double()
  let n = sum_counts(counts)
  if n <= 1.0 {
    return 0.0
  }
  (s - 1.0) / @math.ln(n)
}

///|
/// Menhinick's richness index.
/// D_Mn = S / sqrt(N)
pub fn menhinick(counts : Array[Double]) -> Double {
  let s = observed_count(counts).to_double()
  let n = sum_counts(counts)
  if n == 0.0 {
    return 0.0
  }
  s / n.sqrt()
}

///|
/// McIntosh diversity index.
pub fn mcintosh_d(counts : Array[Double]) -> Double {
  let n = sum_counts(counts)
  let s = observed_count(counts).to_double()
  if n == 0.0 || s <= 1.0 {
    return 0.0
  }
  let mut sum_sq = 0.0
  for c in counts {
    sum_sq = sum_sq + c * c
  }
  let u = sum_sq.sqrt()
  (n - u) / (n - n / s.sqrt())
}

///|
/// McIntosh evenness index.
pub fn mcintosh_e(counts : Array[Double]) -> Double {
  let n = sum_counts(counts)
  let s = observed_count(counts).to_double()
  if n == 0.0 || s <= 1.0 {
    return 1.0
  }
  let mut sum_sq = 0.0
  for c in counts {
    sum_sq = sum_sq + c * c
  }
  let u = sum_sq.sqrt()
  (n - u) / (n - n / s.sqrt())
}

///|
/// Brillouin diversity index.
/// HB = (1/N) * @math.ln(N! / prod(n_i!))
pub fn brillouin_d(counts : Array[Double]) -> Double {
  let n = sum_counts(counts)
  if n == 0.0 {
    return 0.0
  }
  let ln_n_fact = log_factorial(n)
  let mut sum_ln_fact = 0.0
  for c in counts {
    if c > 0.0 {
      sum_ln_fact = sum_ln_fact + log_factorial(c)
    }
  }
  (ln_n_fact - sum_ln_fact) / n
}

///|
/// Natural log factorial (Stirling's approximation).
fn log_factorial(n : Double) -> Double {
  if n < 2.0 {
    return 0.0
  }
  n * @math.ln(n) - n + 0.5 * @math.ln(2.0 * 3.141592653589793 * n)
}

///|
/// Robbins estimator.
pub fn robbins(counts : Array[Double]) -> Double {
  let n = sum_counts(counts)
  if n == 0.0 {
    return 1.0
  }
  let mut f1 = 0.0
  for c in counts {
    if c == 1.0 {
      f1 = f1 + 1.0
    }
  }
  f1 / (n + 1.0)
}

///|
/// Number of singletons (species with count == 1).
pub fn singles(counts : Array[Double]) -> Int {
  let mut n = 0
  for c in counts {
    if c == 1.0 {
      n = n + 1
    }
  }
  n
}

///|
/// Number of doubletons (species with count == 2).
pub fn doubles(counts : Array[Double]) -> Int {
  let mut n = 0
  for c in counts {
    if c == 2.0 {
      n = n + 1
    }
  }
  n
}

// ============================================================
// Beta diversity
// ============================================================

///|
/// Bray-Curtis dissimilarity.
pub fn braycurtis(counts1 : Array[Double], counts2 : Array[Double]) -> Double {
  let len = if counts1.length() < counts2.length() {
    counts1.length()
  } else {
    counts2.length()
  }
  let mut sum_diff = 0.0
  let mut sum_total = 0.0
  for i in 0..<len {
    let diff = (counts1[i] - counts2[i]).abs()
    sum_diff = sum_diff + diff
    sum_total = sum_total + counts1[i] + counts2[i]
  }
  if counts1.length() > len {
    for i in len..<counts1.length() {
      sum_diff = sum_diff + counts1[i].abs()
      sum_total = sum_total + counts1[i]
    }
  }
  if counts2.length() > len {
    for i in len..<counts2.length() {
      sum_diff = sum_diff + counts2[i].abs()
      sum_total = sum_total + counts2[i]
    }
  }
  if sum_total == 0.0 {
    0.0
  } else {
    sum_diff / sum_total
  }
}

///|
/// Jaccard dissimilarity (presence/absence).
pub fn jaccard(counts1 : Array[Double], counts2 : Array[Double]) -> Double {
  let len = if counts1.length() < counts2.length() {
    counts1.length()
  } else {
    counts2.length()
  }
  let mut intersection = 0.0
  let mut union = 0.0
  for i in 0..<len {
    let p1 = counts1[i] > 0.0
    let p2 = counts2[i] > 0.0
    if p1 && p2 {
      intersection = intersection + 1.0
      union = union + 1.0
    } else if p1 || p2 {
      union = union + 1.0
    }
  }
  if counts1.length() > len {
    for i in len..<counts1.length() {
      if counts1[i] > 0.0 {
        union = union + 1.0
      }
    }
  }
  if counts2.length() > len {
    for i in len..<counts2.length() {
      if counts2[i] > 0.0 {
        union = union + 1.0
      }
    }
  }
  if union == 0.0 {
    0.0
  } else {
    1.0 - intersection / union
  }
}

///|
/// Euclidean distance.
pub fn euclidean(counts1 : Array[Double], counts2 : Array[Double]) -> Double {
  let len = if counts1.length() < counts2.length() {
    counts1.length()
  } else {
    counts2.length()
  }
  let mut sum_sq = 0.0
  for i in 0..<len {
    let diff = counts1[i] - counts2[i]
    sum_sq = sum_sq + diff * diff
  }
  if counts1.length() > len {
    for i in len..<counts1.length() {
      sum_sq = sum_sq + counts1[i] * counts1[i]
    }
  }
  if counts2.length() > len {
    for i in len..<counts2.length() {
      sum_sq = sum_sq + counts2[i] * counts2[i]
    }
  }
  sum_sq.sqrt()
}

///|
/// Manhattan / cityblock distance.
pub fn manhattan(counts1 : Array[Double], counts2 : Array[Double]) -> Double {
  let len = if counts1.length() < counts2.length() {
    counts1.length()
  } else {
    counts2.length()
  }
  let mut sum = 0.0
  for i in 0..<len {
    sum = sum + (counts1[i] - counts2[i]).abs()
  }
  if counts1.length() > len {
    for i in len..<counts1.length() {
      sum = sum + counts1[i].abs()
    }
  }
  if counts2.length() > len {
    for i in len..<counts2.length() {
      sum = sum + counts2[i].abs()
    }
  }
  sum
}

///|
/// Chebyshev distance.
pub fn chebyshev(counts1 : Array[Double], counts2 : Array[Double]) -> Double {
  let len = if counts1.length() < counts2.length() {
    counts1.length()
  } else {
    counts2.length()
  }
  let mut max_diff = 0.0
  for i in 0..<len {
    let diff = (counts1[i] - counts2[i]).abs()
    if diff > max_diff {
      max_diff = diff
    }
  }
  if counts1.length() > len {
    for i in len..<counts1.length() {
      if counts1[i].abs() > max_diff {
        max_diff = counts1[i].abs()
      }
    }
  }
  if counts2.length() > len {
    for i in len..<counts2.length() {
      if counts2[i].abs() > max_diff {
        max_diff = counts2[i].abs()
      }
    }
  }
  max_diff
}

///|
/// Cosine dissimilarity = 1 - cosine_similarity.
pub fn cosine(counts1 : Array[Double], counts2 : Array[Double]) -> Double {
  let len = if counts1.length() < counts2.length() {
    counts1.length()
  } else {
    counts2.length()
  }
  let mut dot = 0.0
  let mut norm1 = 0.0
  let mut norm2 = 0.0
  for i in 0..<len {
    dot = dot + counts1[i] * counts2[i]
    norm1 = norm1 + counts1[i] * counts1[i]
    norm2 = norm2 + counts2[i] * counts2[i]
  }
  if counts1.length() > len {
    for i in len..<counts1.length() {
      norm1 = norm1 + counts1[i] * counts1[i]
    }
  }
  if counts2.length() > len {
    for i in len..<counts2.length() {
      norm2 = norm2 + counts2[i] * counts2[i]
    }
  }
  let n1 = norm1.sqrt()
  let n2 = norm2.sqrt()
  if n1 == 0.0 || n2 == 0.0 {
    0.0
  } else {
    1.0 - dot / (n1 * n2)
  }
}

///|
/// Cosine distance (alias for cosine).
pub fn cosine_distance(
  counts1 : Array[Double],
  counts2 : Array[Double],
) -> Double {
  cosine(counts1, counts2)
}

///|
/// Canberra distance.
pub fn canberra(counts1 : Array[Double], counts2 : Array[Double]) -> Double {
  let len = if counts1.length() < counts2.length() {
    counts1.length()
  } else {
    counts2.length()
  }
  let mut sum = 0.0
  for i in 0..<len {
    let abs_diff = (counts1[i] - counts2[i]).abs()
    let abs_sum = counts1[i].abs() + counts2[i].abs()
    if abs_sum > 0.0 {
      sum = sum + abs_diff / abs_sum
    }
  }
  if counts1.length() > len {
    for i in len..<counts1.length() {
      if counts1[i].abs() > 0.0 {
        sum = sum + 1.0
      }
    }
  }
  if counts2.length() > len {
    for i in len..<counts2.length() {
      if counts2[i].abs() > 0.0 {
        sum = sum + 1.0
      }
    }
  }
  sum
}

// ============================================================
// Phylogenetic diversity
// ============================================================

///|
/// Faith's Phylogenetic Diversity (PD).
/// Uses HashMap for O(1) taxon lookup at tip nodes.
pub fn faith_pd(
  counts : Array[Double],
  taxa : Array[String],
  tree : TreeNode,
) -> Double {
  let present_set : @hashmap.HashMap[String, Bool] = @hashmap.HashMap::new()
  for i in 0..<counts.length() {
    if counts[i] > 0.0 && i < taxa.length() {
      present_set.set(taxa[i], true)
    }
  }

  if present_set.is_empty() {
    return 0.0
  }

  compute_pd(tree, present_set)
}

///|
fn compute_pd(
  tree : TreeNode,
  present : @hashmap.HashMap[String, Bool],
) -> Double {
  let (pd, _) = pd_helper(tree, present)
  pd
}

///|
fn pd_helper(
  node : TreeNode,
  present : @hashmap.HashMap[String, Bool],
) -> (Double, Bool) {
  let mut total_pd = 0.0
  let mut has_present = false

  if node.is_tip() {
    match node.name {
      Some(name) => has_present = present.contains(name)
      None => ()
    }
    (total_pd, has_present)
  } else {
    for child in node.children {
      let (child_pd, child_has) = pd_helper(child, present)
      total_pd = total_pd + child_pd
      if child_has {
        has_present = true
        total_pd = total_pd + child.length.unwrap_or(0.0)
      }
    }
    (total_pd, has_present)
  }
}

// ============================================================
// Diversity profile (Hill numbers)
// ============================================================

///|
/// Hill number of order q.
pub fn hill(counts : Array[Double], q : Double) -> Double {
  let total = sum_counts(counts)
  if total == 0.0 {
    return 0.0
  }

  if q == 1.0 {
    @math.exp(shannon(counts))
  } else if q == 0.0 {
    observed_count(counts).to_double()
  } else {
    let mut sum_pq = 0.0
    for c in counts {
      if c > 0.0 {
        let p = c / total
        sum_pq = sum_pq + @math.pow(p, q)
      }
    }
    @math.pow(sum_pq, 1.0 / (1.0 - q))
  }
}

///|
/// Renyi entropy of order q.
pub fn renyi(counts : Array[Double], q : Double) -> Double {
  let total = sum_counts(counts)
  if total == 0.0 {
    return 0.0
  }

  if q == 1.0 {
    shannon(counts)
  } else if q == 0.0 {
    @math.ln(observed_count(counts).to_double())
  } else {
    let mut sum_pq = 0.0
    for c in counts {
      if c > 0.0 {
        let p = c / total
        sum_pq = sum_pq + @math.pow(p, q)
      }
    }
    @math.ln(sum_pq) / (1.0 - q)
  }
}

///|
/// Tsallis entropy of order q.
pub fn tsallis(counts : Array[Double], q : Double) -> Double {
  let total = sum_counts(counts)
  if total == 0.0 {
    return 0.0
  }

  if q == 1.0 {
    shannon(counts)
  } else {
    let mut sum_pq = 0.0
    for c in counts {
      if c > 0.0 {
        let p = c / total
        sum_pq = sum_pq + @math.pow(p, q)
      }
    }
    (1.0 - sum_pq) / (q - 1.0)
  }
}
