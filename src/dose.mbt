///|
pub struct DOSETerm {
  term_id : String
  term_name : String
  gene_count : Int
  genes : Array[String]
  p_value : Double
  adjusted_p_value : Double
  enrichment_score : Double
}

///|
pub struct DOSEResult {
  terms : Array[DOSETerm]
  background_count : Int
  query_count : Int
}

///|
pub fn enrich_disease(
  gene_list : Array[String],
  disease_database : Map[String, Array[String]],
  background_genes : Array[String],
) -> DOSEResult {
  let terms : Array[DOSETerm] = []

  let db_keys = disease_database.keys().collect()
  let mut i = 0
  while i < db_keys.length() {
    let term_id = db_keys[i]
    let term_genes = disease_database[term_id]
    let overlap = count_overlap(gene_list, term_genes)
    let query_size = gene_list.length()
    let term_size = term_genes.length()
    let background_size = background_genes.length()

    let p_value = hypergeometric_test(
      overlap,
      term_size,
      background_size - term_size,
      query_size,
    )
    let enrichment_score = if term_size > 0 {
      overlap.to_double() / term_size.to_double()
    } else {
      0.0
    }

    if overlap > 0 {
      terms.push(DOSETerm::{
        term_id,
        term_name: term_id,
        gene_count: overlap,
        genes: find_overlap(gene_list, term_genes),
        p_value,
        adjusted_p_value: p_value * disease_database.size().to_double(),
        enrichment_score,
      })
    }
    i = i + 1
  }

  let n = terms.length()
  let mut i = 0
  while i < n - 1 {
    let mut j = i + 1
    while j < n {
      if terms[i].p_value > terms[j].p_value {
        let tmp = terms[i]
        terms[i] = terms[j]
        terms[j] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }

  DOSEResult::{
    terms,
    background_count: background_genes.length(),
    query_count: gene_list.length(),
  }
}

///|
fn count_overlap(a : Array[String], b : Array[String]) -> Int {
  let mut count = 0
  let mut i = 0
  while i < a.length() {
    let gene = a[i]
    let mut j = 0
    while j < b.length() {
      if gene == b[j] {
        count = count + 1
      }
      j = j + 1
    }
    i = i + 1
  }
  count
}

///|
fn find_overlap(a : Array[String], b : Array[String]) -> Array[String] {
  let overlap : Array[String] = []
  let mut i = 0
  while i < a.length() {
    let gene = a[i]
    let mut j = 0
    while j < b.length() {
      if gene == b[j] {
        overlap.push(gene)
      }
      j = j + 1
    }
    i = i + 1
  }
  overlap
}

///|
pub fn filter_results(result : DOSEResult, p_cutoff : Double) -> DOSEResult {
  let filtered : Array[DOSETerm] = []
  let mut i = 0
  while i < result.terms.length() {
    let term = result.terms[i]
    if term.adjusted_p_value < p_cutoff {
      filtered.push(term)
    }
    i = i + 1
  }
  DOSEResult::{
    terms: filtered,
    background_count: result.background_count,
    query_count: result.query_count,
  }
}

///|
pub fn get_top_terms(result : DOSEResult, n : Int) -> Array[DOSETerm] {
  let top : Array[DOSETerm] = []
  let mut i = 0
  while i < result.terms.length() && i < n {
    top.push(result.terms[i])
    i = i + 1
  }
  top
}

///|
pub fn summarize_result(result : DOSEResult) -> String {
  let mut summary = "DOSE Enrichment Analysis Summary:\n"
  summary = summary + "- Query genes: " + result.query_count.to_string() + "\n"
  summary = summary +
    "- Background genes: " +
    result.background_count.to_string() +
    "\n"
  summary = summary +
    "- Significant terms: " +
    result.terms.length().to_string() +
    "\n"

  let mut i = 0
  while i < result.terms.length() && i < 5 {
    let term = result.terms[i]
    summary = summary + "  "
    summary = summary + (i + 1).to_string()
    summary = summary + ". "
    summary = summary + term.term_id
    summary = summary + " (p="
    summary = summary + term.p_value.to_string()
    summary = summary + ")\n"
    i = i + 1
  }

  summary
}
