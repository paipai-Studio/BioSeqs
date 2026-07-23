///|
/// topGO - Topology-based GO enrichment analysis.

///|
pub struct TopGOTerm {
  go_id : String
  name : String
  namespace : String
  level : Int
  parents : Array[String]
  children : Array[String]
  gene_count : Int
}

///|
pub fn TopGOTerm::new(go_id : String, name : String, namespace : String, level : Int) -> TopGOTerm {
  TopGOTerm::{
    go_id,
    name,
    namespace,
    level,
    parents: Array::new(),
    children: Array::new(),
    gene_count: 0
  }
}

///|
pub struct TopGOGraph {
  terms : Map[String, TopGOTerm]
  root : String
}

///|
pub fn TopGOGraph::new(root : String) -> TopGOGraph {
  TopGOGraph::{
    terms: Map([], capacity=1000),
    root
  }
}

///|
pub fn TopGOGraph::add_term(self : TopGOGraph, term : TopGOTerm) -> TopGOGraph {
  let terms = self.terms.copy()
  terms.set(term.go_id, term)
  TopGOGraph::{ terms: terms, root: self.root }
}

///|
pub fn TopGOGraph::add_edge(self : TopGOGraph, parent_id : String, child_id : String) -> TopGOGraph {
  let terms = self.terms.copy()
  
  if terms.contains(parent_id) {
    let parent = terms.get(parent_id).unwrap()
    let children = parent.children.copy()
    children.push(child_id)
    terms.set(parent_id, TopGOTerm::{
      go_id: parent.go_id,
      name: parent.name,
      namespace: parent.namespace,
      level: parent.level,
      parents: parent.parents,
      children: children,
      gene_count: parent.gene_count
    })
  }
  
  if terms.contains(child_id) {
    let child = terms.get(child_id).unwrap()
    let parents = child.parents.copy()
    parents.push(parent_id)
    terms.set(child_id, TopGOTerm::{
      go_id: child.go_id,
      name: child.name,
      namespace: child.namespace,
      level: child.level,
      parents: parents,
      children: child.children,
      gene_count: child.gene_count
    })
  }
  
  TopGOGraph::{ terms: terms, root: self.root }
}

///|
pub struct TopGOEnrichmentResult {
  go_id : String
  term_name : String
  namespace : String
  level : Int
  p_value : Double
  adjusted_p_value : Double
  enrichment_score : Double
  gene_count : Int
  total_genes : Int
  algorithm : String
}

///|
pub fn TopGOEnrichmentResult::new(go_id : String, term_name : String, namespace : String, level : Int, p_value : Double, adjusted_p_value : Double, enrichment_score : Double, gene_count : Int, total_genes : Int, algorithm : String) -> TopGOEnrichmentResult {
  TopGOEnrichmentResult::{ go_id, term_name, namespace, level, p_value, adjusted_p_value, enrichment_score, gene_count, total_genes, algorithm }
}

///|
pub fn topgo_fisher_exact(count : Int, total : Int, bg_count : Int, bg_total : Int) -> Double {
  let mut p = 1.0
  let max_x = if count < bg_count { count } else { bg_count }
  let mut x = count
  
  while x >= 0 && x <= max_x && p > 0.0 {
    let numerator = topgo_factorial(bg_count).to_double() * topgo_factorial(bg_total - bg_count).to_double() * topgo_factorial(total).to_double() * topgo_factorial(bg_total - total).to_double()
    let denominator = topgo_factorial(x).to_double() * topgo_factorial(bg_count - x).to_double() * topgo_factorial(total - x).to_double() * topgo_factorial(bg_total - bg_count - total + x).to_double()
    if denominator > 0.0 {
      p = p + numerator / denominator
    }
    x = x - 1
  }
  
  if p > 1.0 { 1.0 } else { p }
}

///|
fn topgo_factorial(n : Int) -> Int {
  if n <= 1 { 1 } else { n * topgo_factorial(n - 1) }
}

///|
pub fn topgo_elim_algorithm(graph : TopGOGraph, genes_of_interest : Array[String], all_genes : Array[String], namespace : String) -> Array[TopGOEnrichmentResult] {
  let results : Array[TopGOEnrichmentResult] = Array::new()
  let visited : Map[String, Bool] = Map([], capacity=1000)
  
  let queue : Array[String] = Array::new()
  queue.push(graph.root)
  
  let mut idx = 0
  while idx < queue.length() {
    let node_id = queue[idx]
    
    if !visited.contains(node_id) {
      visited.set(node_id, true)
      
      if graph.terms.contains(node_id) {
        let term = graph.terms.get(node_id).unwrap()
        
        if namespace == "" || term.namespace == namespace {
          let gene_count = term.gene_count
          let total_genes = all_genes.length()
          let bg_count = gene_count
          let bg_total = total_genes
          
          let p_value = topgo_fisher_exact(gene_count, genes_of_interest.length(), bg_count, bg_total)
          
          let children = term.children.copy()
          let mut j = 0
          while j < children.length() {
            queue.push(children[j])
            j = j + 1
          }
          
          let result = TopGOEnrichmentResult::new(
            term.go_id,
            term.name,
            term.namespace,
            term.level,
            p_value,
            p_value * 0.1,
            1.0 - p_value,
            gene_count,
            total_genes,
            "elim"
          )
          results.push(result)
        }
      }
    }
    
    idx = idx + 1
  }
  
  results
}

///|
pub fn topgo_weight01_algorithm(graph : TopGOGraph, genes_of_interest : Array[String], all_genes : Array[String], namespace : String) -> Array[TopGOEnrichmentResult] {
  let results : Array[TopGOEnrichmentResult] = Array::new()
  
  for key in graph.terms.keys() {
    let term = graph.terms.get(key).unwrap()
    
    if namespace == "" || term.namespace == namespace {
      let gene_count = term.gene_count
      let total_genes = all_genes.length()
      let bg_count = gene_count
      let bg_total = total_genes
      
      let p_value = topgo_fisher_exact(gene_count, genes_of_interest.length(), bg_count, bg_total)
      
      let weight = if term.children.length() > 0 { 0.1 } else { 1.0 }
      let adjusted_p = p_value * weight
      
      let result = TopGOEnrichmentResult::new(
        term.go_id,
        term.name,
        term.namespace,
        term.level,
        p_value,
        adjusted_p,
        1.0 - p_value,
        gene_count,
        total_genes,
        "weight01"
      )
      results.push(result)
    }
  }
  
  results
}

///|
pub fn create_example_topgo_graph() -> TopGOGraph {
  let mut graph = TopGOGraph::new("GO:0008150")
  
  let root = TopGOTerm::new("GO:0008150", "biological_process", "BP", 0)
  graph = graph.add_term(root)
  
  let term1 = TopGOTerm::new("GO:0009987", "cellular_process", "BP", 1)
  let term1 = TopGOTerm::{ go_id: term1.go_id, name: term1.name, namespace: term1.namespace, level: term1.level, parents: term1.parents, children: term1.children, gene_count: 50 }
  graph = graph.add_term(term1)
  graph = graph.add_edge("GO:0008150", "GO:0009987")
  
  let term2 = TopGOTerm::new("GO:0050789", "regulation of biological process", "BP", 1)
  let term2 = TopGOTerm::{ go_id: term2.go_id, name: term2.name, namespace: term2.namespace, level: term2.level, parents: term2.parents, children: term2.children, gene_count: 30 }
  graph = graph.add_term(term2)
  graph = graph.add_edge("GO:0008150", "GO:0050789")
  
  let term3 = TopGOTerm::new("GO:0006810", "transport", "BP", 2)
  let term3 = TopGOTerm::{ go_id: term3.go_id, name: term3.name, namespace: term3.namespace, level: term3.level, parents: term3.parents, children: term3.children, gene_count: 20 }
  graph = graph.add_term(term3)
  graph = graph.add_edge("GO:0009987", "GO:0006810")
  
  graph
}