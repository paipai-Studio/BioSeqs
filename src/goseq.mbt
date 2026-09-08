///|
/// goseq - Gene Ontology enrichment analysis with length bias correction.
/// Inspired by Bioconductor goseq package (Young et al., 2010, Genome Biology).
/// Uses the Wallenius non-central hypergeometric distribution to correct for
/// gene length bias in GO enrichment analysis.

///|
/// Result of goseq analysis for a single GO term.
pub struct GoseqResult {
  go_id : String
  term_name : String
  ontology : String
  p_value : Double
  adjusted_p_value : Double
  overlap : Int
  num_de_in_term : Int
  num_in_term : Int
  num_de : Int
  term_genes : Array[String]
}

///|
pub fn GoseqResult::new(
  go_id : String,
  term_name : String,
  ontology : String,
  p_value : Double,
  adjusted_p_value : Double,
  overlap : Int,
  num_de_in_term : Int,
  num_in_term : Int,
  num_de : Int,
  term_genes : Array[String],
) -> GoseqResult {
  GoseqResult::{
    go_id, term_name, ontology, p_value, adjusted_p_value, overlap,
    num_de_in_term, num_in_term, num_de, term_genes,
  }
}

///|
/// Stirling's approximation for lgamma(x)
fn goseq_lgamma(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x < 7.0 {
    let mut shift = 0.0
    let mut xx = x
    while xx < 7.0 {
      shift = shift - @math.ln(xx)
      xx = xx + 1.0
    }
    shift + goseq_lgamma_stirling(xx)
  } else {
    goseq_lgamma_stirling(x)
  }
}

///|
fn goseq_lgamma_stirling(x : Double) -> Double {
  let ln_2pi = 0.9189385332046727
  let inv_12x = 1.0 / (12.0 * x)
  let inv_360x3 = 1.0 / (360.0 * x * x * x)
  ln_2pi + (x - 0.5) * @math.ln(x) - x + inv_12x - inv_360x3
}

///|
fn goseq_log_choose(n : Int, k : Int) -> Double {
  if k < 0 || k > n {
    return 0.0
  }
  if k == 0 || k == n {
    return 0.0
  }
  let nk = n - k
  goseq_lgamma(n.to_double() + 1.0) - goseq_lgamma(k.to_double() + 1.0) -
    goseq_lgamma(nk.to_double() + 1.0)
}

///|
/// Compute the Wallenius non-central hypergeometric p-value P(X >= k).
/// Uses the Manly (1978) approximation with numerical stability via log-sum-exp.
///
/// Parameters:
///   k - observed overlap between DE genes and GO term genes
///   m - number of genes in the GO term
///   n - number of DE genes
///   n_total - total number of genes in the universe
  ///   odds - odds ratio for the GO term (from PWF)
  pub fn wallenius_pvalue(k : Int, m : Int, n : Int, n_total : Int, odds : Double) -> Double {
    if odds <= 0.0 || odds == 1.0 {
      // Fall back to standard hypergeometric when odds=1 (no bias)
      return hypergeo_pvalue(k, m, n, n_total)
    }

    let upper = if m < n { m } else { n }
    if k > upper {
      return 1.0
    }
    if k <= 0 {
      return 1.0
    }

    // Compute log numerator and denominator for the Manly approximation:
    // P(X >= k) = Σ_{i=k}^{upper} C(m,i) * C(n_total-m,n-i) * odds^i
    //           / Σ_{j=0}^{upper} C(m,j) * C(n_total-m,n-j) * odds^j
    //
    // Use log-space for numerical stability.

    // First pass: compute all log terms
    let log_terms : Array[Double] = Array::make(upper + 1, 0.0)
    let mut i = 0
    while i <= upper {
      let log_choose_mi = goseq_log_choose(m, i)
      let log_choose_nm_ni = goseq_log_choose(n_total - m, n - i)
      let log_odds_term = i.to_double() * @math.ln(odds)
      log_terms[i] = log_choose_mi + log_choose_nm_ni + log_odds_term
      i = i + 1
    }

    // Find max log term for numerical stability (log-sum-exp trick)
    let mut log_max = -1.0e300
    i = 0
    while i <= upper {
      if log_terms[i] > log_max {
        log_max = log_terms[i]
      }
      i = i + 1
    }

    // Compute denominator: sum of exp(log_term - log_max) for all terms
    let mut denom_sum = 0.0
    i = 0
    while i <= upper {
      denom_sum = denom_sum + @math.exp(log_terms[i] - log_max)
      i = i + 1
    }
    let log_denom = log_max + @math.ln(denom_sum)

    // Compute numerator: sum of exp(log_term - log_denom) for k..upper
    let mut numer_sum = 0.0
    i = k
    while i <= upper {
      numer_sum = numer_sum + @math.exp(log_terms[i] - log_denom)
      i = i + 1
    }

    // P = numerator / denominator = numer_sum / 1.0 (since denominator = exp(log_denom))
    let pval = numer_sum

    if pval > 1.0 {
      1.0
    } else {
      pval
    }
}

///|
/// Standard hypergeometric p-value P(X >= k) (no bias correction).
  pub fn hypergeo_pvalue(k : Int, m : Int, n : Int, n_total : Int) -> Double {
    let upper = if m < n { m } else { n }
    if k > upper {
      return 1.0
    }
    if k <= 0 {
      return 1.0
    }

    let log_denom = goseq_log_choose(n_total, n)

    // Find max log term for numerical stability
    let mut log_max = -1.0e300
    let mut i = k
    while i <= upper {
      let log_term = goseq_log_choose(m, i) + goseq_log_choose(n_total - m, n - i) - log_denom
      if log_term > log_max {
        log_max = log_term
      }
      i = i + 1
    }

    let mut sum = 0.0
    i = k
    while i <= upper {
      let log_term = goseq_log_choose(m, i) + goseq_log_choose(n_total - m, n - i) - log_denom
      sum = sum + @math.exp(log_term - log_max)
      i = i + 1
    }

    let pval = @math.exp(log_max) * sum
  if pval > 1.0 {
    1.0
  } else {
    pval
  }
}

///|
/// Calculate the Probability Weighting Function (PWF) using LOWESS.
/// The PWF estimates the probability of a gene being differentially expressed
/// (DE) as a function of its transcript length.
///
/// Parameters:
///   gene_lengths - array of gene lengths (bp)
///   de_status - array of binary indicators (1 = DE, 0 = not DE)
///   f - LOWESS smoothing parameter (fraction of points to use)
///
/// Returns the fitted probabilities (PWF values) for each gene.
pub fn goseq_nullp(
  gene_lengths : Array[Double],
  de_status : Array[Int],
  f? : Double = 0.6667,
) -> Array[Double] {
  let n = gene_lengths.length()
  if n == 0 {
    return Array::new()
  }

  // Sort gene lengths and DE status together
  let sorted_pairs : Array[(Double, Int)] = Array::new()
  let mut i = 0
  while i < n {
    sorted_pairs.push((gene_lengths[i], de_status[i]))
    i = i + 1
  }

  let mut j = 0
  while j < n - 1 {
    let mut k = j + 1
    while k < n {
      if sorted_pairs[j].0 > sorted_pairs[k].0 {
        let tmp = sorted_pairs[j]
        sorted_pairs[j] = sorted_pairs[k]
        sorted_pairs[k] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }

  let sorted_lengths : Array[Double] = Array::make(n, 0.0)
  let sorted_de : Array[Double] = Array::make(n, 0.0)
  i = 0
  while i < n {
    sorted_lengths[i] = sorted_pairs[i].0
    sorted_de[i] = sorted_pairs[i].1.to_double()
    i = i + 1
  }

  // Apply LOWESS to smooth the binary DE status against log(length)
  let log_lengths : Array[Double] = Array::make(n, 0.0)
  i = 0
  while i < n {
    log_lengths[i] = @math.ln(sorted_lengths[i].max(1.0))
    i = i + 1
  }

  let smoothed = lowess(log_lengths, sorted_de, f=f, n_steps=2, delta=0.0)

  // Clamp to [0, 1]
  let result : Array[Double] = Array::make(n, 0.0)
  i = 0
  while i < n {
    let val = smoothed[i]
    result[i] = if val < 0.0 {
      0.0
    } else if val > 1.0 {
      1.0
    } else {
      val
    }
    i = i + 1
  }

  result
}

///|
/// Calculate the odds ratio for a GO term based on the PWF.
/// The odds ratio represents the bias for genes in this GO category.
///
/// Returns the Wallenius odds ratio ω.
fn goseq_calculate_odds(
  term_gene_indices : Array[Int],
  all_gene_indices : Array[Int],
  pwf : Array[Double],
) -> Double {
  let n_term = term_gene_indices.length()
  if n_term == 0 {
    return 1.0
  }

  // Calculate mean odds for genes in the term
  // odds_i = pwf_i / (1 - pwf_i + epsilon)
  let mut term_odds_sum = 0.0
  let mut term_count = 0
  let mut i = 0
  while i < n_term {
    let idx = term_gene_indices[i]
    let pwf_val = pwf[idx]
    let odds = pwf_val / (1.0 - pwf_val + 1.0e-10).max(1.0e-10)
    if !odds.is_nan() && odds.abs() <= 1.0e300 {
      term_odds_sum = term_odds_sum + odds
      term_count = term_count + 1
    }
    i = i + 1
  }

  // Calculate mean odds for genes not in the term
  let in_term_set : Map[Int, Bool] = Map([], capacity=n_term)
  i = 0
  while i < n_term {
    in_term_set.set(term_gene_indices[i], true)
    i = i + 1
  }

  let mut bg_odds_sum = 0.0
  let mut bg_count = 0
  i = 0
  while i < all_gene_indices.length() {
    if !in_term_set.contains(all_gene_indices[i]) {
      let idx = all_gene_indices[i]
      let pwf_val = pwf[idx]
      let odds = pwf_val / (1.0 - pwf_val + 1.0e-10).max(1.0e-10)
      if !odds.is_nan() && odds.abs() <= 1.0e300 {
        bg_odds_sum = bg_odds_sum + odds
        bg_count = bg_count + 1
      }
    }
    i = i + 1
  }

  let mean_term_odds = if term_count > 0 {
    term_odds_sum / term_count.to_double()
  } else {
    1.0
  }

  let mean_bg_odds = if bg_count > 0 {
    bg_odds_sum / bg_count.to_double()
  } else {
    1.0
  }

  // Wallenius odds ratio
  if mean_bg_odds > 0.0 {
    mean_term_odds / mean_bg_odds
  } else {
    1.0
  }
}

///|
/// Benjamini-Hochberg FDR correction.
fn goseq_bh_correction(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return Array::new()
  }

  let indices : Array[Int] = Array::new()
  let mut i = 0
  while i < n {
    indices.push(i)
    i = i + 1
  }

  let sorted_with_idx : Array[(Double, Int)] = Array::new()
  i = 0
  while i < n {
    sorted_with_idx.push((p_values[i], indices[i]))
    i = i + 1
  }

  let mut j = 0
  while j < n - 1 {
    let mut k = j + 1
    while k < n {
      if sorted_with_idx[j].0 > sorted_with_idx[k].0 {
        let tmp = sorted_with_idx[j]
        sorted_with_idx[j] = sorted_with_idx[k]
        sorted_with_idx[k] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }

  let corrected_sorted : Array[Double] = Array::make(n, 0.0)
  i = 0
  while i < n {
    let rank = (i + 1).to_double()
    let p = sorted_with_idx[i].0
    let corrected = p * n.to_double() / rank
    corrected_sorted[i] = if corrected > 1.0 { 1.0 } else { corrected }
    i = i + 1
  }

  for i = n - 2; i >= 0; i = i - 1 {
    if corrected_sorted[i] > corrected_sorted[i + 1] {
      corrected_sorted[i] = corrected_sorted[i + 1]
    }
  }

  let result : Array[Double] = Array::make(n, 0.0)
  i = 0
  while i < n {
    let (_, original_idx) = sorted_with_idx[i]
    result[original_idx] = corrected_sorted[i]
    i = i + 1
  }

  result
}

///|
/// Gene Ontology structure for goseq.
pub struct GoseqGeneSet {
  go_id : String
  term_name : String
  ontology : String
  gene_ids : Array[String]
}

///|
pub fn GoseqGeneSet::new(
  go_id : String,
  term_name : String,
  ontology : String,
  gene_ids : Array[String],
) -> GoseqGeneSet {
  GoseqGeneSet::{ go_id, term_name, ontology, gene_ids }
}

///|
/// Perform GO enrichment analysis with length bias correction.
/// This is the main goseq function, analogous to the goseq() function
/// in the Bioconductor goseq package.
///
/// Parameters:
///   gene_names - array of all gene names/IDs (order must match gene_lengths and de_status)
///   gene_lengths - array of gene lengths in bp
///   de_status - array of binary indicators (1 = DE, 0 = not DE)
///   go_sets - array of GoseqGeneSet representing GO categories
///   f - LOWESS smoothing parameter for PWF estimation
///
/// Returns sorted array of GoseqResult (by p-value, ascending).
pub fn goseq(
  gene_names : Array[String],
  gene_lengths : Array[Double],
  de_status : Array[Int],
  go_sets : Array[GoseqGeneSet],
  f? : Double = 0.6667,
) -> Array[GoseqResult] {
  let n = gene_names.length()
  if n == 0 {
    return Array::new()
  }

  // Step 1: Calculate PWF
  let pwf = goseq_nullp(gene_lengths, de_status, f~)

  // Step 2: Build gene name -> index map
  let gene_to_idx : Map[String, Int] = Map([], capacity=n)
  let mut i = 0
  while i < n {
    gene_to_idx.set(gene_names[i], i)
    i = i + 1
  }

  // Step 3: Count DE genes
  let mut de_count = 0
  i = 0
  while i < n {
    if de_status[i] == 1 {
      de_count = de_count + 1
    }
    i = i + 1
  }

  // Step 4: For each GO term, compute Wallenius p-value
  let results : Array[GoseqResult] = Array::new()
  let p_values : Array[Double] = Array::new()

  let all_indices : Array[Int] = Array::new()
  i = 0
  while i < n {
    all_indices.push(i)
    i = i + 1
  }

  let mut gs = 0
  while gs < go_sets.length() {
    let go_set = go_sets[gs]

    // Find indices of genes in this GO term that are in our gene list
    let term_gene_indices : Array[Int] = Array::new()
    let term_genes_in_list : Array[String] = Array::new()
    let mut gi = 0
    while gi < go_set.gene_ids.length() {
      let gid = go_set.gene_ids[gi]
      match gene_to_idx.get(gid) {
        Some(idx) => {
          term_gene_indices.push(idx)
          term_genes_in_list.push(gid)
        }
        None => ()
      }
      gi = gi + 1
    }

    let m = term_gene_indices.length()
    if m == 0 {
      gs = gs + 1
      continue
    }

    // Count DE genes in this term
    let mut k = 0
    gi = 0
    while gi < term_gene_indices.length() {
      let idx = term_gene_indices[gi]
      if de_status[idx] == 1 {
        k = k + 1
      }
      gi = gi + 1
    }

    if k == 0 {
      gs = gs + 1
      continue
    }

    // Calculate Wallenius odds ratio for this term
    let odds = goseq_calculate_odds(term_gene_indices, all_indices, pwf)

    // Compute Wallenius p-value
    let p_value = wallenius_pvalue(k, m, de_count, n, odds)
    p_values.push(p_value)

    results.push(
      GoseqResult::new(
        go_set.go_id,
        go_set.term_name,
        go_set.ontology,
        p_value,
        0.0,
        k,
        k,
        m,
        de_count,
        term_genes_in_list,
      ),
    )

    gs = gs + 1
  }

  // Step 5: BH correction
  let adjusted = goseq_bh_correction(p_values)

  let mut ri = 0
  while ri < results.length() {
    results[ri] = GoseqResult::{
      go_id: results[ri].go_id,
      term_name: results[ri].term_name,
      ontology: results[ri].ontology,
      p_value: results[ri].p_value,
      adjusted_p_value: adjusted[ri],
      overlap: results[ri].overlap,
      num_de_in_term: results[ri].num_de_in_term,
      num_in_term: results[ri].num_in_term,
      num_de: results[ri].num_de,
      term_genes: results[ri].term_genes,
    }
    ri = ri + 1
  }

  // Sort by p-value ascending
  let mut a = 0
  while a < results.length() - 1 {
    let mut b = a + 1
    while b < results.length() {
      if results[a].p_value > results[b].p_value {
        let tmp = results[a]
        results[a] = results[b]
        results[b] = tmp
      }
      b = b + 1
    }
    a = a + 1
  }

  results
}

///|
/// Filter goseq results by adjusted p-value threshold and minimum overlap.
pub fn goseq_filter_results(
  results : Array[GoseqResult],
  p_adjust_cutoff : Double,
  min_overlap : Int,
) -> Array[GoseqResult] {
  let filtered : Array[GoseqResult] = Array::new()
  let mut i = 0
  while i < results.length() {
    let r = results[i]
    if r.adjusted_p_value <= p_adjust_cutoff && r.overlap >= min_overlap {
      filtered.push(r)
    }
    i = i + 1
  }
  filtered
}

