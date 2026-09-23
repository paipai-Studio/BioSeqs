///|
/// GOSemSim - GO semantic similarity measures.
/// Inspired by the GOSemSim R package (Yu et al., 2010; Bioconductor).
/// Implements Resnik, Lin, Rel (Schlicker), Jiang, and Wang similarity measures
/// based on Information Content (IC) of GO terms in the GO DAG.
///
/// Five measure formulas (canonical from GOSemSim R source):
///   Resnik = IC(MICA)
///   Lin    = 2*IC(MICA) / (IC(id1) + IC(id2))
///   Jiang  = 1 - min(1, IC(id1) + IC(id2) - 2*IC(MICA))
///   Rel    = Lin * (1 - exp(-IC_orig(MICA)))
///   Wang   = sum of SV on common nodes / sum of all SV (with is_a weight=0.8)
///
/// Reference: https://rdrr.io/bioc/GOSemSim/src/R/ICMethods.R
///            https://rdrr.io/bioc/GOSemSim/src/R/WangMethod.R

///|
/// Semantic similarity measure types.
pub enum SimMeasure {
  Resnik
  Lin
  Rel
  Jiang
  Wang
} derive(Eq, Debug)
pub extend SimMeasure with Eq::{equal, not_equal}
pub extend SimMeasure with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Construct Resnik measure.
pub fn resnik_measure() -> SimMeasure {
  Resnik
}

///|
/// Construct Lin measure.
pub fn lin_measure() -> SimMeasure {
  Lin
}

///|
/// Construct Rel (Schlicker) measure.
pub fn rel_measure() -> SimMeasure {
  Rel
}

///|
/// Construct Jiang measure.
pub fn jiang_measure() -> SimMeasure {
  Jiang
}

///|
/// Construct Wang measure (graph-based, no IC needed).
pub fn wang_measure() -> SimMeasure {
  Wang
}

// ──────────────────────────────────────────────────────────────────────
//  GO DAG data structures
// ──────────────────────────────────────────────────────────────────────

///|
/// A node in the GO DAG.
pub struct GOTermNode {
  pub go_id : String
  pub term_name : String
  pub ic : Double
  pub parents : Array[String]
  pub children : Array[String]
} derive(Debug)
pub extend GOTermNode with @moonbitlang/core/debug.Debug::{to_repr}

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
  GOTermNode::{ go_id, term_name, ic, parents, children: [], }
}

///|
/// Gene Ontology Directed Acyclic Graph.
pub struct GOGraph {
  pub root_id : String
  terms : Map[String, GOTermNode]
} derive(Debug)
pub extend GOGraph with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Get root term ID.
pub fn GOGraph::get_root_id(self : GOGraph) -> String {
  self.root_id
}

///|
/// Create an empty GO graph.
pub fn GOGraph::new() -> GOGraph {
  GOGraph::{ root_id: "", terms: Map([], capacity=64), }
}

///|
/// Set the root term ID.
pub fn GOGraph::add_root(self : GOGraph, root_id : String) -> GOGraph {
  GOGraph::{ root_id, terms: self.terms, }
}

///|
/// Add a term to the GO graph and update parent references' children.
pub fn GOGraph::add_term(self : GOGraph, node : GOTermNode) -> GOGraph {
  // Update parents' children lists first (before setting our node)
  let terms : Map[String, GOTermNode] = self.terms
  for parent_id in node.parents {
    match terms.get(parent_id) {
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
        terms.set(parent_id, updated)
      }
      None => ()
    }
  }
  terms.set(node.go_id, node)
  GOGraph::{ root_id: self.root_id, terms, }
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
/// Get all term IDs in the graph.
pub fn GOGraph::term_ids(self : GOGraph) -> Array[String] {
  self.terms.keys().collect()
}

///|
/// Get ancestors of a term (including itself) — iterative BFS.
pub fn GOGraph::ancestors(self : GOGraph, go_id : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let visited : Map[String, Bool] = Map([], capacity=32)
  let mut queue : Array[String] = Array::new()
  queue.push(go_id)
  visited.set(go_id, true)
  while queue.length() > 0 {
    let cur = queue[0]
    let rest : Array[String] = Array::new()
    let mut i = 1
    while i < queue.length() {
      rest.push(queue[i])
      i = i + 1
    }
    queue = rest
    result.push(cur)
    match self.terms.get(cur) {
      Some(term) =>
        for p in term.parents {
          if !visited.get_or_default(p, false) {
            visited.set(p, true)
            queue.push(p)
          }
        }
      None => ()
    }
  }
  result
}

// ──────────────────────────────────────────────────────────────────────
//  IC-based measures
// ──────────────────────────────────────────────────────────────────────

///|
/// Compute semantic similarity between two GO terms.
/// Five canonical measures from GOSemSim R package.
///
/// Reference formulas (from ICMethods.R commented canonical code):
///   Resnik = IC(MICA)
///   Lin    = 2*IC(MICA) / (IC(id1) + IC(id2))
///   Jiang  = 1 - min(1, IC(id1) + IC(id2) - 2*IC(MICA))
///   Rel    = Lin * (1 - exp(-IC_orig(MICA)))
///   Wang   = graph-based, see WangMethod.R SV propagation
///
/// All IC-based measures use ORIGINAL IC values (not normalized by max_IC).
pub fn gosemsim_term_sim(
  graph : GOGraph,
  id1 : String,
  id2 : String,
  measure : SimMeasure,
) -> Double {
  // Wang measure is graph-based (no IC needed), handles missing terms below
  if measure == Wang {
    return gosemsim_wang_sim(graph, id1, id2)
  }

  match (graph.get_term(id1), graph.get_term(id2)) {
    (Some(t1), Some(t2)) => {
      let (_, mica_ic) = gosemsim_mica(graph, id1, id2)
      let ic1 = t1.ic
      let ic2 = t2.ic
      match measure {
        Resnik => mica_ic
        Lin =>
          if ic1 + ic2 > 1.0e-15 {
            2.0 * mica_ic / (ic1 + ic2)
          } else {
            0.0
          }
        Jiang => {
          let delta = ic1 + ic2 - 2.0 * mica_ic
          if delta > 1.0 {
            0.0
          } else {
            1.0 - delta
          }
        }
        Rel => {
          let lin_v =
            if ic1 + ic2 > 1.0e-15 {
              2.0 * mica_ic / (ic1 + ic2)
            } else {
              0.0
            }
          lin_v * (1.0 - @math.exp(0.0 - mica_ic))
        }
        _ => 0.0
      }
    }
    _ => 0.0
  }
}

// ──────────────────────────────────────────────────────────────────────
//  MICA (Most Informative Common Ancestor)
// ──────────────────────────────────────────────────────────────────────

///|
/// Find the Most Informative Common Ancestor (MICA):
/// the common ancestor of id1 and id2 with the highest IC.
/// Returns (mica_id, mica_ic_orig).
/// Falls back to root with IC=0.0 if no common ancestor found.
fn gosemsim_mica(graph : GOGraph, id1 : String, id2 : String) -> (String, Double) {
  let anc1 = graph.ancestors(id1)
  let anc2 = graph.ancestors(id2)

  // Build lookup for anc1
  let seen : Map[String, Bool] = Map([], capacity=anc1.length())
  for a in anc1 {
    seen.set(a, true)
  }

  let mut best_id = ""
  let mut best_ic = -1.0

  for a in anc2 {
    if seen.get_or_default(a, false) {
      match graph.get_term(a) {
        Some(term) =>
          if term.ic > best_ic {
            best_ic = term.ic
            best_id = a
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

// ──────────────────────────────────────────────────────────────────────
//  Wang measure — SV propagation (canonical from WangMethod.R)
// ──────────────────────────────────────────────────────────────────────

///|
/// Wang semantic similarity — pure graph-based, no IC needed.
///
/// Reference: https://rdrr.io/bioc/GOSemSim/src/R/WangMethod.R
///
/// SV (Semantic Value) propagation:
///   1. Initial: SV[term] = 1.0
///   2. Upward propagation: for each parent p of current nodes,
///      SV[p] = weight(is_a) * SV[child], with weight = 0.8 for is_a edges
///   3. topNode SV is forced to 0 (if not DO ontology — for GO this applies)
///   4. Similarity = sum(SV on common nodes for both terms) / sum(all SV for both terms)
pub fn gosemsim_wang_sim(graph : GOGraph, id1 : String, id2 : String) -> Double {
  if id1 == id2 {
    return 1.0
  }
  match (graph.get_term(id1), graph.get_term(id2)) {
    (Some(_), Some(_)) => ()
    _ => return 0.0
  }

  let sv1 = gosemsim_compute_sv(graph, id1)
  let sv2 = gosemsim_compute_sv(graph, id2)

  // Find common nodes (node IDs present in both SV maps)
  let inter_sum = gosemsim_inter_sum(sv1, sv2)
  let total_sum = gosemsim_total_sum(sv1) + gosemsim_total_sum(sv2)

  if total_sum < 1.0e-15 {
    0.0
  } else {
    inter_sum / total_sum
  }
}

///|
/// Compute SV (Semantic Value) vector for a term via upward weight propagation.
///
/// Algorithm (from WangMethod.R getSV):
///   sv[term] = 1.0  (initial)
///   BFS upward: for each child in frontier, each parent gets
///     sv[parent] = weight(relationship) * sv[child]
///   weight for is_a = 0.8 (default GO relationship weight)
///   topNode SV forced to 0.0
///
/// Returns Map[node_id -> sv_value]
fn gosemsim_compute_sv(graph : GOGraph, start_id : String) -> Map[String, Double] {
  let sv : Map[String, Double] = Map([], capacity=16)
  sv.set(start_id, 1.0)

  let root_id = graph.root_id
  let weight : Double = 0.8 // is_a relationship weight

  // BFS upward
  let mut frontier : Array[String] = Array::new()
  frontier.push(start_id)
  let visited : Map[String, Bool] = Map([], capacity=32)
  visited.set(start_id, true)

  while frontier.length() > 0 {
    let cur = frontier[0]
    let rest : Array[String] = Array::new()
    let mut i = 1
    while i < frontier.length() {
      rest.push(frontier[i])
      i = i + 1
    }
    frontier = rest

    let cur_sv = sv.get_or_default(cur, 0.0)
    match graph.get_term(cur) {
      Some(term) =>
        for p in term.parents {
          if !visited.get_or_default(p, false) {
            visited.set(p, true)
            let parent_sv = weight * cur_sv
            sv.set(p, parent_sv)
            frontier.push(p)
          }
        }
      None => ()
    }
  }

  // Force topNode SV = 0 (per GOSemSim R implementation for non-DO ontologies)
  if root_id.length() > 0 {
    sv.set(root_id, 0.0)
  }

  sv
}

///|
/// Sum of SV values for common nodes across both maps.
fn gosemsim_inter_sum(sv_a : Map[String, Double], sv_b : Map[String, Double]) -> Double {
  let mut total = 0.0
  for k in sv_a.keys().collect() {
    let va = sv_a.get_or_default(k, 0.0)
    match sv_b.get(k) {
      Some(vb) => total = total + va + vb
      None => ()
    }
  }
  total
}

///|
/// Sum of all SV values in one map.
fn gosemsim_total_sum(sv : Map[String, Double]) -> Double {
  let mut total = 0.0
  for v in sv.values().to_array() {
    total = total + v
  }
  total
}

// ──────────────────────────────────────────────────────────────────────
//  Gene-level similarity (combine sets of GO terms)
// ──────────────────────────────────────────────────────────────────────

///|
/// Compute semantic similarity between two gene products
/// (each represented by a set of GO terms).
///
/// Combine methods (from GOSemSim R combineScores):
///   "max"     : maximum of term-pair similarities
///   "avg"     : arithmetic mean of all term-pair similarities
///   "rcmax"   : max(row-wise mean, col-wise mean)  (Bioconductor-style)
///   "BMA"     : (sum(row_max) + sum(col_max)) / (n_rows + n_cols)
///   "rcmax.avg": same as BMA
///
/// Default: "max"
pub fn gosemsim_gen_sim(
  graph : GOGraph,
  go_ids1 : Array[String],
  go_ids2 : Array[String],
  measure : SimMeasure,
  combine? : String = "max",
) -> Double {
  let sim_matrix = gosemsim_build_sim_matrix(graph, go_ids1, go_ids2, measure)
  match combine {
    "max" => gosemsim_matrix_max(sim_matrix)
    "avg" => gosemsim_matrix_avg(sim_matrix)
    "rcmax" => gosemsim_matrix_rcmax(sim_matrix)
    "BMA" => gosemsim_matrix_bma(sim_matrix)
    "rcmax.avg" => gosemsim_matrix_bma(sim_matrix)
    _ => gosemsim_matrix_max(sim_matrix)
  }
}

///|
/// Build similarity matrix between two GO term sets.
fn gosemsim_build_sim_matrix(
  graph : GOGraph,
  go_ids1 : Array[String],
  go_ids2 : Array[String],
  measure : SimMeasure,
) -> Array[Array[Double]] {
  let matrix : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < go_ids1.length() {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < go_ids2.length() {
      row.push(gosemsim_term_sim(graph, go_ids1[i], go_ids2[j], measure))
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  matrix
}

///|
/// Max combine: maximum value in matrix.
fn gosemsim_matrix_max(matrix : Array[Array[Double]]) -> Double {
  let mut best = 0.0
  for row in matrix {
    for v in row {
      if v > best {
        best = v
      }
    }
  }
  best
}

///|
/// Avg combine: arithmetic mean of all entries.
fn gosemsim_matrix_avg(matrix : Array[Array[Double]]) -> Double {
  let mut sum = 0.0
  let mut count = 0
  for row in matrix {
    for v in row {
      sum = sum + v
      count = count + 1
    }
  }
  if count > 0 {
    sum / count.to_double()
  } else {
    0.0
  }
}

///|
/// rcmax combine: max(row-wise max mean, col-wise max mean).
/// Reference: GOSemSim R combineMethods.R
///   rowScore = mean(apply(SimScores, 1, max))
///   colScore = mean(apply(SimScores, 2, max))
///   result = max(rowScore, colScore)
fn gosemsim_matrix_rcmax(matrix : Array[Array[Double]]) -> Double {
  let n_rows = matrix.length()
  let n_cols = if n_rows > 0 { matrix[0].length() } else { 0 }
  if n_rows == 0 || n_cols == 0 {
    return 0.0
  }

  // Row-wise max values
  let mut row_max_sum = 0.0
  let mut i = 0
  while i < n_rows {
    let mut row_max = 0.0
    let mut j = 0
    while j < n_cols {
      if matrix[i][j] > row_max {
        row_max = matrix[i][j]
      }
      j = j + 1
    }
    row_max_sum = row_max_sum + row_max
    i = i + 1
  }
  let row_score = row_max_sum / n_rows.to_double()

  // Column-wise max values
  let mut col_max_sum = 0.0
  let mut j = 0
  while j < n_cols {
    let mut col_max = 0.0
    let mut i2 = 0
    while i2 < n_rows {
      if matrix[i2][j] > col_max {
        col_max = matrix[i2][j]
      }
      i2 = i2 + 1
    }
    col_max_sum = col_max_sum + col_max
    j = j + 1
  }
  let col_score = col_max_sum / n_cols.to_double()

  if row_score > col_score {
    row_score
  } else {
    col_score
  }
}

///|
/// BMA / rcmax.avg combine: (sum(row_max) + sum(col_max)) / (n_rows + n_cols).
/// Reference: GOSemSim R combineMethods.R
///   result = ( sum(apply(SimScores, 1, max)) + sum(apply(SimScores, 2, max)) ) / sum(dim(SimScores))
fn gosemsim_matrix_bma(matrix : Array[Array[Double]]) -> Double {
  let n_rows = matrix.length()
  let n_cols = if n_rows > 0 { matrix[0].length() } else { 0 }
  if n_rows == 0 || n_cols == 0 {
    return 0.0
  }

  let mut row_max_sum = 0.0
  let mut i = 0
  while i < n_rows {
    let mut row_max = 0.0
    let mut j = 0
    while j < n_cols {
      if matrix[i][j] > row_max {
        row_max = matrix[i][j]
      }
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
      if matrix[i2][j] > col_max {
        col_max = matrix[i2][j]
      }
      i2 = i2 + 1
    }
    col_max_sum = col_max_sum + col_max
    j = j + 1
  }

  (row_max_sum + col_max_sum) / (n_rows + n_cols).to_double()
}

// ──────────────────────────────────────────────────────────────────────
//  Example graph (must match MoonBit Python frozen-value generator)
// ──────────────────────────────────────────────────────────────────────

///|
/// Build a small example GO DAG for testing.
///
/// Topology (all is_a edges, root=GO:0008150 biological_process):
///
///          GO:0008150 (ic=0.0)
///           /              \
///   GO:0050896(0.3)     GO:0009058(0.2)
///    /      \               |
/// GO:0009987 GO:0023046   GO:0010467(0.35)
/// (0.5)       (0.4)
///   |           |
/// GO:0006810  GO:0007165
/// (0.7)        (0.6)
///
/// Used by both MoonBit tests and Python frozen-value generator
/// (/tmp/gosemsim_frozen.py) to ensure exact correspondence.
pub fn gosemsim_example_graph() -> GOGraph {
  let g = GOGraph::new()
  let g1 = g.add_root("GO:0008150")
  let g2 = g1.add_term(
    GOTermNode::new(
      "GO:0008150",
      ic=0.0,
      parents=[],
      term_name="biological_process",
    ),
  )
  let g3 = g2.add_term(
    GOTermNode::new(
      "GO:0050896",
      ic=0.3,
      parents=["GO:0008150"],
      term_name="response_to_stimulus",
    ),
  )
  let g4 = g3.add_term(
    GOTermNode::new(
      "GO:0009987",
      ic=0.5,
      parents=["GO:0050896"],
      term_name="cellular_process",
    ),
  )
  let g5 = g4.add_term(
    GOTermNode::new(
      "GO:0006810",
      ic=0.7,
      parents=["GO:0009987"],
      term_name="transport",
    ),
  )
  let g6 = g3.add_term(
    GOTermNode::new(
      "GO:0023046",
      ic=0.4,
      parents=["GO:0050896"],
      term_name="signaling",
    ),
  )
  let g7 = g6.add_term(
    GOTermNode::new(
      "GO:0007165",
      ic=0.6,
      parents=["GO:0023046"],
      term_name="signal_transduction",
    ),
  )
  let g8 = g2.add_term(
    GOTermNode::new(
      "GO:0009058",
      ic=0.2,
      parents=["GO:0008150"],
      term_name="biosynthetic_process",
    ),
  )
  let g9 = g8.add_term(
    GOTermNode::new(
      "GO:0010467",
      ic=0.35,
      parents=["GO:0009058"],
      term_name="gene_expression",
    ),
  )
  g9
}
