///|
/// DOSE - Disease Ontology Semantic and Enrichment analysis.
/// Ported from Bioconductor DOSE package (Yu et al., 2020).
///
/// Core functionality:
///
/// - Fisher's exact test on a 2x2 contingency table (two-sided /
///   greater / less) with Haldane–Anscombe corrected odds ratio.
/// - Hypergeometric enrichment test (one-sided "greater", equivalent
///   to Fisher exact enrichment).
/// - Benjamini–Hochberg FDR correction and Bonferroni correction.
/// - Fold enrichment and odds ratio for each term.
/// - Batch enrichment: ``enrich_disease_fisher`` and
///   ``enrich_disease_hyper`` over a disease term database.
///
/// For the 2x2 contingency table the cells are labeled
///
///   | a = overlap  | b = query − overlap  |
///   | c = term − a | d = bg − term − b    |
///
/// with row margins a+b (query size) and c+d (term size),
/// column margins a+c (overlap in term) and b+d (overlap out of term),
/// and total a+b+c+d = background size.

// ===========================================================================
// Private mathematical helpers
// ===========================================================================

///|
/// Log-factorial ln(n!) computed by summation of natural logs.
/// Returns 0.0 for n ≤ 1; no overflow for large n because the logs stay small.
fn dose_log_factorial(n : Int) -> Double {
  if n <= 1 {
    return 0.0
  }
  let mut result = 0.0
  let mut i = 2
  while i <= n {
    result = result + @math.ln(i.to_double())
    i = i + 1
  }
  result
}

///|
/// Log binomial coefficient ln C(n, k) = ln(n!) − ln(k!) − ln((n−k)!).
/// Returns a very negative sentinel (−1e300) when k is out of range.
fn dose_log_binom(n : Int, k : Int) -> Double {
  if k < 0 || k > n {
    return -1.0e300
  }
  dose_log_factorial(n) - dose_log_factorial(k) - dose_log_factorial(n - k)
}

// ===========================================================================
// Fisher's exact test
// ===========================================================================

///|
/// Result of a Fisher's exact test on a 2x2 contingency table.
pub struct FisherExactResult {
  p_value : Double
  odds_ratio : Double
}

///|
/// Two-sided Fisher's exact test on the 2x2 table
///
///   | a | b |
///   | c | d |
///
/// Returns a ``FisherExactResult`` whose ``p_value`` is computed over
/// the hypergeometric distribution of the (1,1) cell summed over all
/// tables whose probability does not exceed that of the observed table.
///
/// The odds ratio uses Haldane–Anscombe correction (add 0.5 to every
/// cell) whenever any cell is zero, matching the convention in
/// Bioconductor DOSE and R's ``fisher.test``.
pub fn fisher_exact_two_sided(
  a : Int,
  b : Int,
  c : Int,
  d : Int,
) -> FisherExactResult {
  dose_fisher_impl(a, b, c, d, "two-sided")
}

///|
/// One-sided Fisher's exact test for enrichment ("greater" tail, P(X ≥ a)).
pub fn fisher_exact_greater(
  a : Int,
  b : Int,
  c : Int,
  d : Int,
) -> FisherExactResult {
  dose_fisher_impl(a, b, c, d, "greater")
}

///|
/// One-sided Fisher's exact test for depletion ("less" tail, P(X ≤ a)).
pub fn fisher_exact_less(
  a : Int,
  b : Int,
  c : Int,
  d : Int,
) -> FisherExactResult {
  dose_fisher_impl(a, b, c, d, "less")
}

///|
/// Internal implementation shared by the three public Fisher wrappers.
fn dose_fisher_impl(
  a : Int,
  b : Int,
  c : Int,
  d : Int,
  test_type : String,
) -> FisherExactResult {
  let n_total = a + b + c + d
  let row1 = a + b
  let row2 = c + d
  let col1 = a + c

  // Odds ratio with Haldane–Anscombe correction for zero cells.
  let odds_ratio = if a == 0 || b == 0 || c == 0 || d == 0 {
    let a2 = a.to_double() + 0.5
    let b2 = b.to_double() + 0.5
    let c2 = c.to_double() + 0.5
    let d2 = d.to_double() + 0.5
    a2 * d2 / (b2 * c2)
  } else {
    a.to_double() * d.to_double() / (b.to_double() * c.to_double())
  }

  // Possible range for the (1,1) cell under fixed margins.
  let k_min = if col1 - row2 > 0 { col1 - row2 } else { 0 }
  let k_max = if col1 < row1 { col1 } else { row1 }

  let log_denom = dose_log_binom(n_total, col1)
  let log_p_observed = dose_log_binom(row1, a) +
    dose_log_binom(row2, c) -
    log_denom
  let p_observed = @math.exp(log_p_observed)

  let mut p_value = 0.0
  let mut k = k_min
  while k <= k_max {
    let log_p_k = dose_log_binom(row1, k) +
      dose_log_binom(row2, col1 - k) -
      log_denom
    let p_k = @math.exp(log_p_k)
    match test_type {
      "greater" => if k >= a { p_value = p_value + p_k }
      "less" => if k <= a { p_value = p_value + p_k }
      _ => if p_k <= p_observed + 1.0e-15 { p_value = p_value + p_k }
    }
    k = k + 1
  }

  if p_value > 1.0 {
    p_value = 1.0
  }
  if p_value < 0.0 {
    p_value = 0.0
  }

  FisherExactResult::{ p_value, odds_ratio }
}

// ===========================================================================
// Hypergeometric enrichment test (Fisher exact greater, explicit API)
// ===========================================================================

///|
/// Hypergeometric enrichment test.
///
/// Parameters:
///   - ``overlap`` (k): number of query genes present in the term
///   - ``query_size`` (m): size of the query gene list
///   - ``term_size`` (msize): size of the disease / GO term
///   - ``background_size`` (n_bg): size of the background gene universe
///
/// Returns the enrichment p-value P(X ≥ k) under a Hypergeometric
/// (n_bg, term_size, query_size) distribution. This is numerically
/// identical to ``fisher_exact_greater`` on the 2x2 table above.
pub fn dose_hypergeometric_test(
  overlap : Int,
  query_size : Int,
  term_size : Int,
  background_size : Int,
) -> Double {
  let k = overlap
  let m = query_size
  let ms = term_size
  let nb = background_size
  let a = k
  let b = m - k
  let c = ms - k
  let d = nb - ms - (m - k)
  fisher_exact_greater(a, b, c, d).p_value
}

// ===========================================================================
// Multiple testing corrections (namespaced with dose_ to avoid collisions)
// ===========================================================================

///|
/// Benjamini–Hochberg FDR correction.
///
/// For p-values sorted in ascending order p_(1) ≤ … ≤ p_(m), the
/// BH-adjusted value at rank i is
///
///   q_(i) = min_{j≥i} (p_(j) * m / j)
///
/// so the output array has the same size as the input array and preserves
/// the original order (ranked element ``q[i]`` corresponds to ``p[i]``).
pub fn dose_bh_fdr(p_values : Array[Double]) -> Array[Double] {
  let m = p_values.length()
  if m == 0 {
    return []
  }

  // Build indexed array, sort by p-value ascending.
  let indexed : Array[Array[Double]] = Array::make(m, [0.0, 0.0])
  let mut i = 0
  while i < m {
    indexed[i] = [p_values[i], i.to_double()]
    i = i + 1
  }
  indexed.sort_by(fn(a : Array[Double], b : Array[Double]) -> Int {
    if a[0] < b[0] {
      -1
    } else if a[0] > b[0] {
      1
    } else {
      0
    }
  })

  let result : Array[Double] = Array::make(m, 0.0)
  let mut min_so_far = 1.0e300
  // Backward pass to enforce monotonicity of q-values.
  let mut idx = m - 1
  while idx >= 0 {
    let rank = (idx + 1).to_double()
    let raw_q = indexed[idx][0] * m.to_double() / rank
    let q = if raw_q < min_so_far { raw_q } else { min_so_far }
    min_so_far = q
    let orig = indexed[idx][1].to_int()
    result[orig] = if q > 1.0 { 1.0 } else { q }
    idx = idx - 1
  }
  result
}

///|
/// Bonferroni correction: ``adj[i] = min(p[i] * m, 1.0)`` where ``m``
/// is the number of tests.
pub fn dose_bonferroni(p_values : Array[Double]) -> Array[Double] {
  let m = p_values.length()
  if m == 0 {
    return []
  }
  let mm = m.to_double()
  let result : Array[Double] = Array::make(m, 0.0)
  let mut i = 0
  while i < m {
    let adj = p_values[i] * mm
    result[i] = if adj > 1.0 { 1.0 } else { adj }
    i = i + 1
  }
  result
}

// ===========================================================================
// Fold enrichment / odds ratio helpers
// ===========================================================================

///|
/// Fold enrichment (also called "enrichment ratio").
///
///   FE = (overlap / query_size) / (term_size / background_size)
///      = k * n_bg / (m * msize)
///
/// Returns 0.0 when ``overlap == 0``.  Callers should ensure that
/// ``query_size``, ``term_size``, and ``background_size`` are
/// strictly positive; otherwise the result is ``@double.not_a_number``.
pub fn dose_fold_enrichment(
  overlap : Int,
  query_size : Int,
  term_size : Int,
  background_size : Int,
) -> Double {
  if overlap <= 0 {
    return 0.0
  }
  let k = overlap.to_double()
  let m = query_size.to_double()
  let ms = term_size.to_double()
  let nb = background_size.to_double()
  if m <= 0.0 || ms <= 0.0 || nb <= 0.0 {
    return 0.0 / 0.0
  }
  k * nb / (m * ms)
}

///|
/// Odds ratio of the 2x2 table
///
///   | overlap | query − overlap |
///   | term − overlap | bg − term − (query − overlap) |
///
/// With Haldane–Anscombe 0.5-cell correction whenever any cell is zero,
/// matching the convention in Bioconductor DOSE.
pub fn dose_odds_ratio(
  overlap : Int,
  query_size : Int,
  term_size : Int,
  background_size : Int,
) -> Double {
  let a = overlap
  let b = query_size - overlap
  let c = term_size - overlap
  let d = background_size - term_size - (query_size - overlap)
  dose_fisher_impl(a, b, c, d, "two-sided").odds_ratio
}

// ===========================================================================
// Enrichment result types
// ===========================================================================

///|
/// One row of an enrichment analysis result.
pub struct DOSETerm {
  term_id : String
  term_name : String
  gene_count : Int
  genes : Array[String]
  overlap_count : Int
  p_value : Double
  adjusted_p_value : Double
  fold_enrichment : Double
  odds_ratio : Double
}

///|
/// Full enrichment result (one row per term, sorted by adjusted p).
pub struct DOSEResult {
  terms : Array[DOSETerm]
  query_count : Int
  background_count : Int
  n_terms : Int
}

// ===========================================================================
// Batch enrichment
// ===========================================================================

///|
/// Perform enrichment analysis using Fisher's exact test for every term
/// in the supplied disease / gene-set database.
///
/// Parameters:
///   - ``gene_list``: query gene symbols
///   - ``gene_database``: term_id → Array[gene symbols]
///   - ``background_genes``: background gene universe
///   - ``p_adjust_method``: "BH" (default) or "bonferroni"
///
/// The result is sorted by adjusted p-value (ascending) with ties broken
/// by fold enrichment descending.
pub fn enrich_disease_fisher(
  gene_list : Array[String],
  gene_database : Map[String, Array[String]],
  background_genes : Array[String],
  p_adjust_method? : String = "BH",
) -> DOSEResult {
  let n_query = gene_list.length()
  let n_bg = background_genes.length()

  // Build a set of background genes for fast look-up.
  let bg_set : Map[String, Bool] = Map([])
  let mut i = 0
  while i < n_bg {
    bg_set[background_genes[i]] = true
    i = i + 1
  }

  let term_ids = gene_database.keys().collect()
  let n_terms = term_ids.length()

  // Collect overlaps and store raw p-values per term.
  let rows : Array[Array[Double]] = Array::make(n_terms, [0.0, 0.0, 0.0, 0.0, 0.0])
  let overlap_genes : Array[Array[String]] = Array::make(n_terms, Array::new())

  i = 0
  while i < n_terms {
    let term_id = term_ids[i]
    let term_genes = gene_database[term_id]
    let n_term = term_genes.length()

    // Count overlap, restricted to genes present in the background universe.
    let overlap : Array[String] = Array::new()
    let mut j = 0
    while j < gene_list.length() {
      let g = gene_list[j]
      if bg_set.contains(g) {
        let mut found = false
        let mut k = 0
        while k < n_term {
          if g == term_genes[k] {
            found = true
            break
          }
          k = k + 1
        }
        if found {
          overlap.push(g)
        }
      }
      j = j + 1
    }

    let k_overlap = overlap.length()
    let m_eff = n_query
    let ms_eff = n_term
    let nb_eff = n_bg

    let a = k_overlap
    let b = m_eff - k_overlap
    let c = ms_eff - k_overlap
    let d = nb_eff - ms_eff - (m_eff - k_overlap)
    let fr = dose_fisher_impl(a, b, c, d, "greater")

    let fe = if k_overlap == 0 { 0.0 } else {
      k_overlap.to_double() * nb_eff.to_double() /
        (m_eff.to_double() * ms_eff.to_double())
    }
    let or_val = if a == 0 || b == 0 || c == 0 || d == 0 {
      let a2 = a.to_double() + 0.5
      let b2 = b.to_double() + 0.5
      let c2 = c.to_double() + 0.5
      let d2 = d.to_double() + 0.5
      a2 * d2 / (b2 * c2)
    } else {
      a.to_double() * d.to_double() / (b.to_double() * c.to_double())
    }

    rows[i] = [fr.p_value, fe, or_val, a.to_double(), n_term.to_double()]
    overlap_genes[i] = overlap
    i = i + 1
  }

  // Multiple testing correction.
  let raw_pvals : Array[Double] = Array::make(n_terms, 0.0)
  i = 0
  while i < n_terms {
    raw_pvals[i] = rows[i][0]
    i = i + 1
  }
  let adj_pvals = match p_adjust_method {
    "bonferroni" => dose_bonferroni(raw_pvals)
    _ => dose_bh_fdr(raw_pvals)
  }

  // Build term list and sort by adjusted p ascending (ties → FE desc).
  let term_objs : Array[DOSETerm] = Array::new()
  i = 0
  while i < n_terms {
    let term_id = term_ids[i]
    let n_genes_term = rows[i][4].to_int()
    let overlap_count = rows[i][3].to_int()
    term_objs.push(DOSETerm::{
      term_id,
      term_name: term_id,
      gene_count: n_genes_term,
      genes: overlap_genes[i],
      overlap_count,
      p_value: rows[i][0],
      adjusted_p_value: adj_pvals[i],
      fold_enrichment: rows[i][1],
      odds_ratio: rows[i][2],
    })
    i = i + 1
  }

  term_objs.sort_by(fn(x : DOSETerm, y : DOSETerm) -> Int {
    if x.adjusted_p_value < y.adjusted_p_value {
      -1
    } else if x.adjusted_p_value > y.adjusted_p_value {
      1
    } else {
      if x.fold_enrichment > y.fold_enrichment {
        -1
      } else if x.fold_enrichment < y.fold_enrichment {
        1
      } else {
        0
      }
    }
  })

  DOSEResult::{
    terms: term_objs,
    query_count: n_query,
    background_count: n_bg,
    n_terms,
  }
}

///|
/// Perform enrichment analysis using the hypergeometric test
/// (numerically identical to Fisher's exact test greater-tail) and
/// Benjamini–Hochberg FDR correction.
pub fn enrich_disease_hyper(
  gene_list : Array[String],
  gene_database : Map[String, Array[String]],
  background_genes : Array[String],
) -> DOSEResult {
  enrich_disease_fisher(gene_list, gene_database, background_genes)
}

// ===========================================================================
// DOSEResult accessors and filtering
// ===========================================================================

///|
/// Keep only terms whose adjusted p-value is below the given cutoff.
pub fn dose_filter_by_adjp(result : DOSEResult, cutoff : Double) -> DOSEResult {
  let filtered : Array[DOSETerm] = Array::new()
  let mut i = 0
  while i < result.terms.length() {
    if result.terms[i].adjusted_p_value < cutoff {
      filtered.push(result.terms[i])
    }
    i = i + 1
  }
  DOSEResult::{
    terms: filtered,
    query_count: result.query_count,
    background_count: result.background_count,
    n_terms: result.n_terms,
  }
}

///|
/// Keep at most ``n`` top terms (lowest adjusted p-value).
pub fn dose_top_terms(result : DOSEResult, n : Int) -> DOSEResult {
  let top : Array[DOSETerm] = Array::new()
  let mut i = 0
  while i < result.terms.length() && i < n {
    top.push(result.terms[i])
    i = i + 1
  }
  DOSEResult::{
    terms: top,
    query_count: result.query_count,
    background_count: result.background_count,
    n_terms: result.n_terms,
  }
}

///|
/// Pretty-printable summary of the enrichment result.
pub fn dose_summarize(result : DOSEResult) -> String {
  let mut s = "DOSE Enrichment Summary:\n"
  s = s + "Query genes: " + result.query_count.to_string() + "\n"
  s = s + "Background genes: " + result.background_count.to_string() + "\n"
  s = s + "Terms tested: " + result.n_terms.to_string() + "\n"
  let mut i = 0
  while i < result.terms.length() && i < 5 {
    let term = result.terms[i]
    s = s +
      "  " +
      (i + 1).to_string() +
      ". " +
      term.term_id +
      " (overlap=" +
      term.overlap_count.to_string() +
      ", adj_p=" +
      term.adjusted_p_value.to_string() +
      ")\n"
    i = i + 1
  }
  s
}

///|
/// Term identifiers in the current sorted order.
pub fn dose_term_ids(result : DOSEResult) -> Array[String] {
  let ids : Array[String] = Array::new()
  let mut i = 0
  while i < result.terms.length() {
    ids.push(result.terms[i].term_id)
    i = i + 1
  }
  ids
}

///|
/// Adjusted p-values, matching the current term order.
pub fn dose_adjusted_pvalues(result : DOSEResult) -> Array[Double] {
  let vals : Array[Double] = Array::new()
  let mut i = 0
  while i < result.terms.length() {
    vals.push(result.terms[i].adjusted_p_value)
    i = i + 1
  }
  vals
}
