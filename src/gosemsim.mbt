///|
/// GOSemSim - GO semantic similarity measures.
/// Inspired by the GOSemSim R package (Yu et al., 2010).
/// Implements Resnik, Lin, Rel, Jiang, and Wang similarity measures
/// based on Information Content (IC) of GO terms in the GO DAG.

///|
/// Semantic similarity measure types.
pub enum SimMeasure {
  Resnik
  Lin
  Rel
  Jiang
  Wang
}

///|
/// Construct Resnik measure.
pub fn resnik_measure() -> SimMeasure { Resnik }

///|
/// Construct Lin measure.
pub fn lin_measure() -> SimMeasure { Lin }

///|
/// Construct Rel measure.
pub fn rel_measure() -> SimMeasure { Rel }

///|
/// Construct Jiang measure.
pub fn jiang_measure() -> SimMeasure { Jiang }

///|
/// Construct Wang measure.
pub fn wang_measure() -> SimMeasure { Wang }

///|
/// A node in the GO graph.
pub struct GOTermNode {
  pub go_id : String
  pub term_name : String
  pub ic : Double
  pub parents : Array[String]
  pub children : Array[String]
}

///|
/// Get GO term ID.
pub fn GOTermNode::get_go_id(self : GOTermNode) -> String {
  self.go_id
}

///|
/// Get term name.
pub fn GOTermNode::get_term_name(self : GOTermNode) -> String {
  self.term_name
}

///|
/// Get information content.
pub fn GOTermNode::get_ic(self : GOTermNode) -> Double {
  self.ic
}

///|
/// Get parent IDs.
pub fn GOTermNode::get_parents(self : GOTermNode) -> Array[String] {
  self.parents
}

///|
/// Get child IDs.
pub fn GOTermNode::get_children(self : GOTermNode) -> Array[String] {
  self.children
}

///|
/// Create a GO term node.
pub fn GOTermNode::new(
  go_id : String,
  ic? : Double = 0.0,
  parents? : Array[String] = [],
  term_name? : String = "",
) -> GOTermNode {
  GOTermNode::{ go_id, term_name, ic, parents, children: [] }
}

///|
/// Gene Ontology Directed Acyclic Graph.
pub struct GOGraph {
  pub root_id : String
  terms : Map[String, GOTermNode]
}

///|
/// Get root term ID.
pub fn GOGraph::get_root_id(self : GOGraph) -> String {
  self.root_id
}

///|
/// Create an empty GO graph.
pub fn GOGraph::new() -> GOGraph {
  GOGraph::{ root_id: "", terms: Map([], capacity=64) }
}

///|
/// Set the root term ID.
pub fn GOGraph::add_root(self : GOGraph, root_id : String) -> GOGraph {
  GOGraph::{ root_id, terms: self.terms }
}

///|
/// Add a term to the GO graph and update parent references.
pub fn GOGraph::add_term(self : GOGraph, node : GOTermNode) -> GOGraph {
  self.terms.set(node.go_id, node)
  // Update parent's children list
  match self.terms.get(node.go_id) {
    Some(_) => {
      for parent_id in node.parents {
        match self.terms.get(parent_id) {
          Some(parent) => {
            let new_children : Array[String] = Array::new()
            for c in parent.children {
              new_children.push(c)
            }
            new_children.push(node.go_id)
            let updated = GOTermNode::{
              go_id: parent.go_id,
              term_name: parent.term_name,
              ic: parent.ic,
              parents: parent.parents,
              children: new_children,
            }
            self.terms.set(parent_id, updated)
          }
          None => ()
        }
      }
    }
    None => ()
  }
  GOGraph::{ root_id: self.root_id, terms: self.terms }
}

///|
/// Get a term by ID.
pub fn GOGraph::get_term(self : GOGraph, go_id : String) -> GOTermNode? {
  match self.terms.get(go_id) {
    Some(node) => Some(node)
    None => None
  }
}

///|
/// Compute semantic similarity between two GO terms.
pub fn gosemsim_term_sim(
  graph : GOGraph,
  id1 : String,
  id2 : String,
  measure : SimMeasure,
) -> Double {
  match (graph.get_term(id1), graph.get_term(id2)) {
    (Some(t1), Some(t2)) => {
      let (mica_id, mica_ic) = gosemsim_mica_ic(graph, id1, id2)
      match measure {
        Resnik => mica_ic
        Lin => {
          if t1.ic + t2.ic > 1.0e-10 {
            2.0 * mica_ic / (t1.ic + t2.ic)
          } else {
            0.0
          }
        }
        Rel => {
          if t1.ic + t2.ic > 1.0e-10 {
            let lin = 2.0 * mica_ic / (t1.ic + t2.ic)
            let rel = lin * (1.0 - mica_ic / (if t1.ic > t2.ic { t1.ic } else { t2.ic }))
            rel
          } else {
            0.0
          }
        }
        Jiang => {
          let common_ancestors = gosemsim_ancestors(graph, id1)
          let d1 = gosemsim_distance_to(graph, id1, common_ancestors)
          let d2 = gosemsim_distance_to(graph, id2, common_ancestors)
          let d_mica = gosemsim_distance_to(graph, mica_id, common_ancestors)
          let combined = d1 + d2 - 2.0 * d_mica
          if combined > 1.0e-10 {
            1.0 / combined
          } else {
            1.0
          }
        }
        Wang => {
          let s1 = gosemsim_wang_score(graph, id1)
          let s2 = gosemsim_wang_score(graph, id2)
          let s_common = gosemsim_wang_common_score(graph, id1, id2)
          if s1 + s2 > 1.0e-10 {
            s_common / (s1 + s2 - s_common)
          } else {
            0.0
          }
        }
        _ => 0.0
      }
    }
    _ => 0.0
  }
}

///|
/// Compute semantic similarity between two gene products (sets of GO terms).
pub fn gosemsim_gen_sim(
  graph : GOGraph,
  go_ids1 : Array[String],
  go_ids2 : Array[String],
  measure : SimMeasure,
  combine? : String = "max",
) -> Double {
  let sim_matrix : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < go_ids1.length() {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < go_ids2.length() {
      row.push(gosemsim_term_sim(graph, go_ids1[i], go_ids2[j], measure))
      j = j + 1
    }
    sim_matrix.push(row)
    i = i + 1
  }
  
  match combine {
    "max" => gosemsim_matrix_max(sim_matrix)
    "avg" => gosemsim_matrix_avg(sim_matrix)
    "rcmax" => gosemsim_matrix_rcmax(sim_matrix)
    _ => gosemsim_matrix_max(sim_matrix)
  }
}

fn gosemsim_matrix_max(matrix : Array[Array[Double]]) -> Double {
  let mut max_val = 0.0
  for row in matrix {
    for v in row {
      if v > max_val { max_val = v }
    }
  }
  max_val
}

fn gosemsim_matrix_avg(matrix : Array[Array[Double]]) -> Double {
  let mut sum = 0.0
  let mut count = 0
  for row in matrix {
    for v in row {
      sum = sum + v
      count = count + 1
    }
  }
  if count > 0 { sum / count.to_double() } else { 0.0 }
}

fn gosemsim_matrix_rcmax(matrix : Array[Array[Double]]) -> Double {
  let n_rows = matrix.length()
  let n_cols = if n_rows > 0 { matrix[0].length() } else { 0 }
  
  // Row-wise max and column-wise max, average of the larger set
  let mut row_max_sum = 0.0
  let mut i = 0
  while i < n_rows {
    let mut row_max = 0.0
    let mut j = 0
    while j < n_cols {
      if matrix[i][j] > row_max { row_max = matrix[i][j] }
      j = j + 1
    }
    row_max_sum = row_max_sum + row_max
    i = i + 1
  }
  
  let mut col_max_sum = 0.0
  let mut j = 0
  while j < n_cols {
    let mut col_max = 0.0
    let mut i2 = 0
    while i2 < n_rows {
      if matrix[i2][j] > col_max { col_max = matrix[i2][j] }
      i2 = i2 + 1
    }
    col_max_sum = col_max_sum + col_max
    j = j + 1
  }
  
  let n = if n_rows > n_cols { n_rows } else { n_cols }
  if n > 0 {
    (row_max_sum + col_max_sum) / (2.0 * n.to_double())
  } else {
    0.0
  }
}

///|
/// Find the Most Informative Common Ancestor (MICA) and its IC.
fn gosemsim_mica_ic(graph : GOGraph, id1 : String, id2 : String) -> (String, Double) {
  let anc1 = gosemsim_ancestors(graph, id1)
  let anc2 = gosemsim_ancestors(graph, id2)
  
  // Find common ancestors
  let seen : Map[String, Bool] = Map([], capacity=anc1.length())
  for a in anc1 {
    seen.set(a, true)
  }
  
  let mut best_id = ""
  let mut best_ic = -1.0
  
  for a in anc2 {
    if seen.get_or_default(a, false) {
      match graph.get_term(a) {
        Some(term) => {
          if term.ic > best_ic {
            best_ic = term.ic
            best_id = a
          }
        }
        None => ()
      }
    }
  }
  
  if best_ic < 0.0 {
    (graph.root_id, 0.0)
  } else {
    (best_id, best_ic)
  }
}

///|
/// Get all ancestors of a term (including itself).
fn gosemsim_ancestors(graph : GOGraph, go_id : String) -> Array[String] {
  let ancestors : Array[String] = Array::new()
  let visited : Map[String, Bool] = Map([], capacity=32)
  gosemsim_ancestors_recursive(graph, go_id, ancestors, visited)
  ancestors
}

fn gosemsim_ancestors_recursive(
  graph : GOGraph,
  go_id : String,
  ancestors : Array[String],
  visited : Map[String, Bool],
) -> Unit {
  if visited.get_or_default(go_id, false) {
    return
  }
  visited.set(go_id, true)
  ancestors.push(go_id)

  match graph.get_term(go_id) {
    Some(term) => {
      for parent in term.parents {
        gosemsim_ancestors_recursive(graph, parent, ancestors, visited)
      }
    }
    None => ()
  }
}

///|
/// Compute path distance to a common ancestor set.
fn gosemsim_distance_to(graph : GOGraph, go_id : String, ancestors : Array[String]) -> Double {
  // Simple BFS to find shortest path to any ancestor
  if ancestors.length() == 0 {
    return 0.0
  }
  
  let target_set : Map[String, Bool] = Map([], capacity=ancestors.length())
  for a in ancestors {
    target_set.set(a, true)
  }
  
  let mut current : Array[(String, Double)] = Array::new()
  let visited : Map[String, Bool] = Map([], capacity=64)
  current.push((go_id, 0.0))
  visited.set(go_id, true)
  
  let mut found = false
  while current.length() > 0 {
    let (node_id, dist) = current[0]
    let rest : Array[(String, Double)] = Array::new()
    let mut i = 1
    while i < current.length() {
      rest.push(current[i])
      i = i + 1
    }
    current = rest
    
    if target_set.get_or_default(node_id, false) {
      found = true
      return dist
    }
    
    match graph.get_term(node_id) {
      Some(term) => {
        for parent in term.parents {
          if !visited.get_or_default(parent, false) {
            visited.set(parent, true)
            current.push((parent, dist + 1.0))
          }
        }
      }
      None => ()
    }
  }
  
  if found { 0.0 } else { 100.0 }
}

///|
/// Compute Wang's similarity score components.
fn gosemsim_wang_score(graph : GOGraph, go_id : String) -> Double {
  let ancestors = gosemsim_ancestors(graph, go_id)
  let mut score = 0.0
  for anc_id in ancestors {
    match graph.get_term(anc_id) {
      Some(term) => {
        // Score based on depth from root
        score = score + 1.0
        let mut parent_count = term.parents.length().to_double()
        if parent_count < 1.0 { parent_count = 1.0 }
        score = score + 1.0 / parent_count
      }
      None => ()
    }
  }
  score
}

///|
/// Compute Wang common score (overlap of ancestor sets).
fn gosemsim_wang_common_score(graph : GOGraph, id1 : String, id2 : String) -> Double {
  let anc1 = gosemsim_ancestors(graph, id1)
  let anc2 = gosemsim_ancestors(graph, id2)
  
  let set1 : Map[String, Bool] = Map([], capacity=anc1.length())
  for a in anc1 {
    set1.set(a, true)
  }
  
  let mut score = 0.0
  for anc_id in anc2 {
    if set1.get_or_default(anc_id, false) {
      score = score + 1.0
      match graph.get_term(anc_id) {
        Some(term) => {
          let mut parent_count = term.parents.length().to_double()
          if parent_count < 1.0 { parent_count = 1.0 }
          score = score + 1.0 / parent_count
        }
        None => ()
      }
    }
  }
  score
}

///|
/// Build a small example GO graph for testing.
pub fn gosemsim_example_graph() -> GOGraph {
  let g = GOGraph::new()
  let g1 = g.add_root("GO:0008150")
  let g2 = g1.add_term(GOTermNode::new("GO:0008150", ic=0.0, parents=[], term_name="biological_process"))
  let g3 = g2.add_term(GOTermNode::new("GO:0050896", ic=0.3, parents=["GO:0008150"], term_name="response to stimulus"))
  let g4 = g3.add_term(GOTermNode::new("GO:0009987", ic=0.5, parents=["GO:0050896"], term_name="cellular process"))
  let g5 = g4.add_term(GOTermNode::new("GO:0006810", ic=0.7, parents=["GO:0009987"], term_name="transport"))
  let g6 = g3.add_term(GOTermNode::new("GO:0023046", ic=0.4, parents=["GO:0050896"], term_name="signaling"))
  let g7 = g6.add_term(GOTermNode::new("GO:0007165", ic=0.6, parents=["GO:0023046"], term_name="signal transduction"))
  let g8 = g2.add_term(GOTermNode::new("GO:0009058", ic=0.2, parents=["GO:0008150"], term_name="biosynthetic process"))
  let g9 = g8.add_term(GOTermNode::new("GO:0010467", ic=0.35, parents=["GO:0009058"], term_name="gene expression"))
  g9
}
