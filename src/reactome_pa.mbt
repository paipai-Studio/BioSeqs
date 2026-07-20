///|
pub struct ReactomePathway {
  pathway_id : String
  pathway_name : String
  gene_count : Int
  genes : Array[String]
  p_value : Double
  adjusted_p_value : Double
  enrichment_score : Double
  top_level_term : String
}

///|
pub struct ReactomeResult {
  pathways : Array[ReactomePathway]
  background_count : Int
  query_count : Int
}

///|
pub fn enrich_reactome(
  gene_list : Array[String],
  reactome_database : Map[String, (String, Array[String], String)],
  background_genes : Array[String],
) -> ReactomeResult {
  let pathways : Array[ReactomePathway] = []

  let db_keys = reactome_database.keys().collect()
  let mut i = 0
  while i < db_keys.length() {
    let pathway_id = db_keys[i]
    let (pathway_name, pathway_genes, top_level) = reactome_database[pathway_id]
    let overlap = reactome_count_overlap(gene_list, pathway_genes)
    let query_size = gene_list.length()
    let pathway_size = pathway_genes.length()
    let background_size = background_genes.length()

    let p_value = hypergeometric_test(
      overlap,
      pathway_size,
      background_size - pathway_size,
      query_size,
    )
    let enrichment_score = if pathway_size > 0 {
      overlap.to_double() / pathway_size.to_double()
    } else {
      0.0
    }

    if overlap > 0 {
      pathways.push(ReactomePathway::{
        pathway_id,
        pathway_name,
        gene_count: overlap,
        genes: reactome_find_overlap(gene_list, pathway_genes),
        p_value,
        adjusted_p_value: p_value * reactome_database.size().to_double(),
        enrichment_score,
        top_level_term: top_level,
      })
    }
    i = i + 1
  }

  let n = pathways.length()
  let mut i = 0
  while i < n - 1 {
    let mut j = i + 1
    while j < n {
      if pathways[i].p_value > pathways[j].p_value {
        let tmp = pathways[i]
        pathways[i] = pathways[j]
        pathways[j] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }

  ReactomeResult::{
    pathways,
    background_count: background_genes.length(),
    query_count: gene_list.length(),
  }
}

///|
fn reactome_count_overlap(a : Array[String], b : Array[String]) -> Int {
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
fn reactome_find_overlap(a : Array[String], b : Array[String]) -> Array[String] {
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
pub fn reactome_filter_results(
  result : ReactomeResult,
  p_cutoff : Double,
) -> ReactomeResult {
  let filtered : Array[ReactomePathway] = []
  let mut i = 0
  while i < result.pathways.length() {
    let pathway = result.pathways[i]
    if pathway.adjusted_p_value < p_cutoff {
      filtered.push(pathway)
    }
    i = i + 1
  }
  ReactomeResult::{
    pathways: filtered,
    background_count: result.background_count,
    query_count: result.query_count,
  }
}

///|
pub fn get_top_pathways(
  result : ReactomeResult,
  n : Int,
) -> Array[ReactomePathway] {
  let top : Array[ReactomePathway] = []
  let mut i = 0
  while i < result.pathways.length() && i < n {
    top.push(result.pathways[i])
    i = i + 1
  }
  top
}

///|
pub fn group_by_top_level(
  result : ReactomeResult,
) -> Map[String, Array[ReactomePathway]] {
  let groups = Map([], capacity=16)

  let mut i = 0
  while i < result.pathways.length() {
    let pathway = result.pathways[i]
    let current = if groups.contains(pathway.top_level_term) {
      groups[pathway.top_level_term]
    } else {
      []
    }
    current.push(pathway)
    groups[pathway.top_level_term] = current
    i = i + 1
  }

  groups
}

///|
pub fn reactome_summarize_result(result : ReactomeResult) -> String {
  let mut summary = "Reactome Pathway Analysis Summary:\n"
  summary = summary + "- Query genes: " + result.query_count.to_string() + "\n"
  summary = summary +
    "- Background genes: " +
    result.background_count.to_string() +
    "\n"
  summary = summary +
    "- Significant pathways: " +
    result.pathways.length().to_string() +
    "\n"

  let mut i = 0
  while i < result.pathways.length() && i < 5 {
    let pathway = result.pathways[i]
    summary = summary + "  "
    summary = summary + (i + 1).to_string()
    summary = summary + ". "
    summary = summary + pathway.pathway_name
    summary = summary + " (p="
    summary = summary + pathway.p_value.to_string()
    summary = summary + ")\n"
    i = i + 1
  }

  summary
}
