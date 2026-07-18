///|
pub struct ClusterProfilerResult {
  terms : Array[ClusterProfilerTerm]
  query_count : Int
  background_count : Int
  organism : String
  keytype : String
}

///|
pub struct ClusterProfilerTerm {
  term_id : String
  term_name : String
  gene_count : Int
  genes : Array[String]
  p_value : Double
  adjusted_p_value : Double
  enrichment_score : Double
  gene_ratio : Double
  bg_ratio : Double
}

///|
pub fn cp_enrich_generic(
  gene_list : Array[String],
  term_database : Map[String, (String, Array[String])],
  background_genes : Array[String],
  organism : String,
  keytype : String,
) -> ClusterProfilerResult {
  let terms : Array[ClusterProfilerTerm] = []
  let n = background_genes.length().to_double()
  let k = gene_list.length().to_double()

  let db_keys = term_database.keys().collect()
  let mut i = 0
  while i < db_keys.length() {
    let term_id = db_keys[i]
    let (term_name, term_genes) = term_database[term_id]

    let mut overlap_count = 0
    let overlap_genes : Array[String] = []

    let mut j = 0
    while j < gene_list.length() {
      let gene = gene_list[j]
      let mut found = false
      let mut m = 0
      while m < term_genes.length() {
        if term_genes[m] == gene {
          found = true
          overlap_count = overlap_count + 1
          overlap_genes.push(gene)
        }
        m = m + 1
      }
      j = j + 1
    }

    if overlap_count > 0 {
      let m_size = term_genes.length().to_double()
      let x = overlap_count.to_double()

      let p_value = cp_hypergeometric_test(x, m_size, k, n)
      let adjusted_p_value = p_value * db_keys.length().to_double()

      let gene_ratio = x / k
      let bg_ratio = m_size / n
      let enrichment_score = if bg_ratio > 0.0 {
        gene_ratio / bg_ratio
      } else {
        1.0
      }

      terms.push(ClusterProfilerTerm::{
        term_id,
        term_name,
        gene_count: overlap_count,
        genes: overlap_genes,
        p_value,
        adjusted_p_value,
        enrichment_score,
        gene_ratio,
        bg_ratio,
      })
    }

    i = i + 1
  }

  let n_terms = terms.length()
  let mut idx = 0
  while idx < n_terms - 1 {
    let mut j = idx + 1
    while j < n_terms {
      if terms[idx].p_value > terms[j].p_value {
        let tmp = terms[idx]
        terms[idx] = terms[j]
        terms[j] = tmp
      }
      j = j + 1
    }
    idx = idx + 1
  }

  ClusterProfilerResult::{
    terms,
    query_count: gene_list.length(),
    background_count: background_genes.length(),
    organism,
    keytype,
  }
}

///|
fn cp_hypergeometric_test(
  x : Double,
  m_size : Double,
  k_size : Double,
  n_size : Double,
) -> Double {
  let mut p = 0.0
  let mut i = x
  while i <= m_size && i <= k_size {
    let term = cp_combination(m_size, i) *
      cp_combination(n_size - m_size, k_size - i) /
      cp_combination(n_size, k_size)
    p = p + term
    i = i + 1.0
  }
  if p > 1.0 {
    1.0
  } else {
    p
  }
}

///|
fn cp_combination(n : Double, k : Double) -> Double {
  if k < 0.0 || k > n {
    0.0
  } else {
    let mut result = 1.0
    let mut i = 0.0
    while i < k {
      result = result * (n - i) / (i + 1.0)
      i = i + 1.0
    }
    result
  }
}

///|
pub fn cp_enrich_go(
  gene_list : Array[String],
  go_database : Map[String, (String, Array[String], String)],
  background_genes : Array[String],
  organism : String,
  ontology : String,
) -> ClusterProfilerResult {
  let filtered_db = Map([], capacity=go_database.length())

  let db_keys = go_database.keys().collect()
  let mut i = 0
  while i < db_keys.length() {
    let key = db_keys[i]
    let (name, genes, ont) = go_database[key]
    if ontology == "" || ont == ontology {
      filtered_db[key] = (name, genes)
    }
    i = i + 1
  }

  cp_enrich_generic(gene_list, filtered_db, background_genes, organism, "GO")
}

///|
pub fn cp_enrich_kegg(
  gene_list : Array[String],
  kegg_database : Map[String, (String, Array[String])],
  background_genes : Array[String],
  organism : String,
) -> ClusterProfilerResult {
  cp_enrich_generic(
    gene_list, kegg_database, background_genes, organism, "KEGG",
  )
}

///|
pub fn cp_filter_results(
  result : ClusterProfilerResult,
  cutoff : Double,
) -> ClusterProfilerResult {
  let filtered : Array[ClusterProfilerTerm] = []

  let mut i = 0
  while i < result.terms.length() {
    let term = result.terms[i]
    if term.adjusted_p_value < cutoff {
      filtered.push(term)
    }
    i = i + 1
  }

  ClusterProfilerResult::{
    terms: filtered,
    query_count: result.query_count,
    background_count: result.background_count,
    organism: result.organism,
    keytype: result.keytype,
  }
}

///|
pub fn cp_get_top_terms(
  result : ClusterProfilerResult,
  n : Int,
) -> Array[ClusterProfilerTerm] {
  let top : Array[ClusterProfilerTerm] = []
  let mut i = 0
  while i < result.terms.length() && i < n {
    top.push(result.terms[i])
    i = i + 1
  }
  top
}

///|
pub fn cp_dotplot(result : ClusterProfilerResult, n : Int) -> String {
  let mut plot = "Enrichment Dot Plot:\n"
  plot = plot + "Term\tGene Ratio\tP-value\tGene Count\n"

  let top = cp_get_top_terms(result, n)
  let mut i = 0
  while i < top.length() {
    let term = top[i]
    plot = plot + term.term_name
    plot = plot + "\t"
    plot = plot + term.gene_ratio.to_string()
    plot = plot + "\t"
    plot = plot + term.p_value.to_string()
    plot = plot + "\t"
    plot = plot + term.gene_count.to_string()
    plot = plot + "\n"
    i = i + 1
  }

  plot
}

///|
pub fn cp_enrich_map(result : ClusterProfilerResult) -> String {
  let mut map = "Enrichment Map:\n"
  map = map + "Number of terms: " + result.terms.length().to_string() + "\n"
  map = map + "Average gene ratio: "

  let mut total_ratio = 0.0
  let mut i = 0
  while i < result.terms.length() {
    total_ratio = total_ratio + result.terms[i].gene_ratio
    i = i + 1
  }

  if result.terms.length() > 0 {
    map = map + (total_ratio / result.terms.length().to_double()).to_string()
  } else {
    map = map + "0"
  }
  map = map + "\n"
  map = map + "Significant terms (adjusted p < 0.05): "

  let mut sig_count = 0
  i = 0
  while i < result.terms.length() {
    if result.terms[i].adjusted_p_value < 0.05 {
      sig_count = sig_count + 1
    }
    i = i + 1
  }

  map = map + sig_count.to_string() + "\n"
  map
}

///|
pub fn cp_pairwise_termsim(
  result : ClusterProfilerResult,
) -> Map[(String, String), Double] {
  let similarity = Map(
    [],
    capacity=result.terms.length() * result.terms.length(),
  )

  let mut i = 0
  while i < result.terms.length() {
    let mut j = i + 1
    while j < result.terms.length() {
      let term1 = result.terms[i]
      let term2 = result.terms[j]

      let mut common = 0
      let mut k = 0
      while k < term1.genes.length() {
        let mut m = 0
        while m < term2.genes.length() {
          if term1.genes[k] == term2.genes[m] {
            common = common + 1
          }
          m = m + 1
        }
        k = k + 1
      }

      let union = term1.genes.length() + term2.genes.length() - common
      let sim = if union > 0 {
        common.to_double() / union.to_double()
      } else {
        0.0
      }

      similarity[(term1.term_id, term2.term_id)] = sim
      j = j + 1
    }
    i = i + 1
  }

  similarity
}

///|
pub fn cp_simplify(
  result : ClusterProfilerResult,
  cutoff : Double,
) -> ClusterProfilerResult {
  let similarities = cp_pairwise_termsim(result)

  let to_remove : Array[String] = []

  let keys = similarities.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let (id1, id2) = keys[i]
    if similarities[(id1, id2)] > cutoff {
      let mut idx1 = -1
      let mut idx2 = -1
      let mut j = 0
      while j < result.terms.length() {
        if result.terms[j].term_id == id1 {
          idx1 = j
        }
        if result.terms[j].term_id == id2 {
          idx2 = j
        }
        j = j + 1
      }

      if idx1 >= 0 && idx2 >= 0 {
        if result.terms[idx1].p_value > result.terms[idx2].p_value {
          to_remove.push(id1)
        } else {
          to_remove.push(id2)
        }
      }
    }
    i = i + 1
  }

  let filtered : Array[ClusterProfilerTerm] = []
  let mut j = 0
  while j < result.terms.length() {
    let term = result.terms[j]
    let mut keep = true
    let mut k = 0
    while k < to_remove.length() {
      if to_remove[k] == term.term_id {
        keep = false
      }
      k = k + 1
    }
    if keep {
      filtered.push(term)
    }
    j = j + 1
  }

  ClusterProfilerResult::{
    terms: filtered,
    query_count: result.query_count,
    background_count: result.background_count,
    organism: result.organism,
    keytype: result.keytype,
  }
}

///|
pub fn cp_cnetplot(result : ClusterProfilerResult, n : Int) -> String {
  let mut plot = "Category-Network Plot:\n"
  plot = plot + "Nodes: Genes and Terms\n"
  plot = plot + "Edges: Gene-Term associations\n\n"

  let top = cp_get_top_terms(result, n)
  let genes_seen : Map[String, Unit] = Map([], capacity=100)

  let mut i = 0
  while i < top.length() {
    let term = top[i]
    plot = plot + "Term: " + term.term_name + "\n"
    plot = plot + "  Genes: "

    let mut j = 0
    while j < term.genes.length() {
      let gene = term.genes[j]
      genes_seen[gene] = ()
      if j > 0 {
        plot = plot + ", "
      }
      plot = plot + gene
      j = j + 1
    }
    plot = plot + "\n\n"
    i = i + 1
  }

  plot = plot + "Total unique genes: " + genes_seen.length().to_string() + "\n"
  plot
}
