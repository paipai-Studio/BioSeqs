///|
/// GO Enrichment Analysis module, inspired by Bioconductor GOstats/clusterProfiler.
/// Provides functions for Gene Ontology enrichment analysis using hypergeometric test.

///|
pub struct GOTerm {
  go_id : String
  term_name : String
  namespace : String
  genes : Array[String]
}

///|
pub fn GOTerm::new(
  go_id : String,
  term_name : String,
  namespace : String,
  genes : Array[String],
) -> GOTerm {
  GOTerm::{ go_id, term_name, namespace, genes }
}

///|
pub struct GOEnrichmentResult {
  go_id : String
  term_name : String
  namespace : String
  p_value : Double
  mut adjusted_p_value : Double
  odds_ratio : Double
  gene_ratio : String
  bg_ratio : String
  count : Int
  bg_count : Int
  genes : Array[String]
}

///|
pub fn GOEnrichmentResult::new(
  go_id : String,
  term_name : String,
  namespace : String,
  p_value : Double,
  adjusted_p_value : Double,
  odds_ratio : Double,
  gene_ratio : String,
  bg_ratio : String,
  count : Int,
  bg_count : Int,
  genes : Array[String],
) -> GOEnrichmentResult {
  GOEnrichmentResult::{
    go_id,
    term_name,
    namespace,
    p_value,
    adjusted_p_value,
    odds_ratio,
    gene_ratio,
    bg_ratio,
    count,
    bg_count,
    genes,
  }
}

///|
fn log_factorial(n : Int) -> Double {
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
fn log_choose(n : Int, k : Int) -> Double {
  if k < 0 || k > n {
    return 0.0
  }
  if k == 0 || k == n {
    return 0.0
  }
  log_factorial(n) - log_factorial(k) - log_factorial(n - k)
}

///|
pub fn hypergeometric_test(x : Int, m : Int, n : Int, k : Int) -> Double {
  let max_x = if m < k { m } else { k }
  let min_x = if n + k - m > 0 { n + k - m } else { 0 }

  let mut sum = 0.0
  let mut i = min_x
  while i <= max_x {
    let log_prob = log_choose(m, i) +
      log_choose(n, k - i) -
      log_choose(m + n, k)
    sum = sum + @math.exp(log_prob)
    i = i + 1
  }

  let mut p_value = 0.0
  let mut j = x
  while j <= max_x {
    let log_prob = log_choose(m, j) +
      log_choose(n, k - j) -
      log_choose(m + n, k)
    p_value = p_value + @math.exp(log_prob)
    j = j + 1
  }

  if p_value > 1.0 {
    1.0
  } else {
    p_value
  }
}

///|
pub fn fishers_exact_test(a : Int, b : Int, c : Int, d : Int) -> Double {
  let n = a + b + c + d
  let log_prob = log_choose(a + b, a) +
    log_choose(c + d, c) +
    log_choose(a + c, a) +
    log_choose(b + d, b) -
    log_choose(n, a) -
    log_choose(n, b) -
    log_choose(n, c)
  @math.exp(log_prob)
}

///|
pub fn bonferroni_correction(p_values : Array[Double]) -> Array[Double] {
  let m = p_values.length().to_double()
  let corrected : Array[Double] = Array::new()
  let mut i = 0
  while i < p_values.length() {
    let corrected_val = p_values[i] * m
    corrected.push(if corrected_val > 1.0 { 1.0 } else { corrected_val })
    i = i + 1
  }
  corrected
}

///|
pub fn bh_correction(p_values : Array[Double]) -> Array[Double] {
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

  let corrected_sorted : Array[Double] = Array::new()
  i = 0
  while i < n {
    let rank = (i + 1).to_double()
    let p = sorted_with_idx[i].0
    let corrected = p * n.to_double() / rank
    corrected_sorted.push(if corrected > 1.0 { 1.0 } else { corrected })
    i = i + 1
  }

  for i = n - 2; i >= 0; i = i - 1 {
    if corrected_sorted[i] > corrected_sorted[i + 1] {
      corrected_sorted[i] = corrected_sorted[i + 1]
    }
  }

  let result : Array[Double] = Array::new()
  i = 0
  while i < n {
    result.push(0.0)
    i = i + 1
  }

  i = 0
  while i < n {
    let (_, original_idx) = sorted_with_idx[i]
    result[original_idx] = corrected_sorted[i]
    i = i + 1
  }

  result
}

///|
pub fn go_enrich(
  genes : Array[String],
  go_terms : Array[GOTerm],
  universe_size : Int,
) -> Array[GOEnrichmentResult] {
  let gene_set = Array::new()
  let mut i = 0
  while i < genes.length() {
    gene_set.push(genes[i])
    i = i + 1
  }

  let results : Array[GOEnrichmentResult] = Array::new()
  let p_values : Array[Double] = Array::new()

  i = 0
  while i < go_terms.length() {
    let go_term = go_terms[i]

    let mut count = 0
    let enriched_genes : Array[String] = Array::new()
    let mut j = 0
    while j < gene_set.length() {
      let mut k = 0
      while k < go_term.genes.length() {
        if gene_set[j] == go_term.genes[k] {
          count = count + 1
          enriched_genes.push(gene_set[j])
          break
        }
        k = k + 1
      }
      j = j + 1
    }

    if count == 0 {
      i = i + 1
      continue
    }

    let m = go_term.genes.length()
    let n = universe_size - m
    let k = gene_set.length()
    let x = count

    let p_value = hypergeometric_test(x, m, n, k)
    p_values.push(p_value)

    let gene_ratio = count.to_string() + "/" + k.to_string()
    let bg_ratio = m.to_string() + "/" + universe_size.to_string()

    let odds_ratio = if x > 0 && m - x > 0 && k - x > 0 && n - k + x > 0 {
      x.to_double() *
      (n - k + x).to_double() /
      ((k - x).to_double() * (m - x).to_double())
    } else {
      0.0
    }

    results.push(
      GOEnrichmentResult::new(
        go_term.go_id,
        go_term.term_name,
        go_term.namespace,
        p_value,
        0.0,
        odds_ratio,
        gene_ratio,
        bg_ratio,
        count,
        m,
        enriched_genes,
      ),
    )

    i = i + 1
  }

  let adjusted_p_values = bh_correction(p_values)

  i = 0
  while i < results.length() {
    results[i].adjusted_p_value = adjusted_p_values[i]
    i = i + 1
  }

  let mut j = 0
  while j < results.length() - 1 {
    let mut k = j + 1
    while k < results.length() {
      if results[j].adjusted_p_value > results[k].adjusted_p_value {
        let tmp = results[j]
        results[j] = results[k]
        results[k] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }

  results
}

///|
pub fn filter_enrichment(
  results : Array[GOEnrichmentResult],
  p_adjust_cutoff : Double,
  min_count : Int,
) -> Array[GOEnrichmentResult] {
  let filtered : Array[GOEnrichmentResult] = Array::new()
  let mut i = 0
  while i < results.length() {
    let r = results[i]
    if r.adjusted_p_value <= p_adjust_cutoff && r.count >= min_count {
      filtered.push(r)
    }
    i = i + 1
  }
  filtered
}

///|
pub fn go_terms_by_namespace(
  go_terms : Array[GOTerm],
  namespace : String,
) -> Array[GOTerm] {
  let filtered : Array[GOTerm] = Array::new()
  let mut i = 0
  while i < go_terms.length() {
    if go_terms[i].namespace == namespace {
      filtered.push(go_terms[i])
    }
    i = i + 1
  }
  filtered
}
