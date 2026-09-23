///|
/// Bio.Pathway.Rep - Graph data structures ported from Biopython 1.88.
///
/// Two variants are provided:
/// - **PathwayGraph**      — a simple directed graph where each edge label is
///   unique.  ``add_edge(from, to, label)`` silently ignores subsequent
///   additions with the same label.
/// - **PathwayMultiGraph** — a directed multigraph that keeps parallel edges
///   distinct; duplicate labels are allowed.
///
/// Both structures expose the same inspection API (``nodes``,
/// ``labels``, ``edges``, ``children``, ``parents``, ``child_edges``,
/// ``parent_edges``).  Breadth-first and depth-first traversals are
/// provided as top-level functions so that they work with any graph.
///
/// Naming note: the pre-existing ``graph.mbt`` module already exposes a
/// weighted ``Graph`` / ``Node`` / ``Edge`` family, so the types here
/// are prefixed with ``Pathway_`` to avoid collision.

// ===========================================================================
// PathwayGraph  (simple directed graph with labelled edges)
// ===========================================================================

///|
/// A simple directed graph whose edges carry unique string labels.
///
/// The canonical Biopython API:
///
/// ```
/// g = PathwayGraph::new()
/// g.add_node("A"); g.add_node("B"); g.add_edge("A", "B", "e1")
/// g.nodes()        // -> ["A", "B"]
/// g.labels()       // -> ["e1"]
/// g.edges("e1")    // -> [("A", "B")]
/// g.children("A")  // -> ["e1"]          (edge labels leaving A)
/// g.parents("B")   // -> ["e1"]          (edge labels arriving at B)
/// g.child_edges("A")   // -> [("B", "e1")]
/// g.parent_edges("B")  // -> [("A", "e1")]
/// ```
pub struct PathwayGraph {
  nodes : Map[String, Bool]
  /// label -> (source, target)
  edge_map : Map[String, (String, String)]
}

///|
/// Create an empty PathwayGraph.
pub fn PathwayGraph::new() -> PathwayGraph {
  PathwayGraph::{ nodes: Map([]), edge_map: Map([]) }
}

///|
/// Register a node (no-op if already present).
pub fn PathwayGraph::add_node(self : PathwayGraph, label : String) -> Unit {
  self.nodes[label] = true
}

///|
/// Add a directed edge ``from -> to`` identified by ``label``.
/// Both endpoints are auto-registered as nodes.  Duplicate edge
/// labels are silently ignored (PathwayGraph invariant).
pub fn PathwayGraph::add_edge(self : PathwayGraph, from_node : String, to_node : String, label : String) -> Unit {
  self.add_node(from_node)
  self.add_node(to_node)
  if !self.edge_map.contains(label) {
    self.edge_map[label] = (from_node, to_node)
  }
}

///|
/// The set of node labels.
pub fn PathwayGraph::nodes(self : PathwayGraph) -> Array[String] {
  self.nodes.keys().collect()
}

///|
/// The set of edge labels.
pub fn PathwayGraph::labels(self : PathwayGraph) -> Array[String] {
  self.edge_map.keys().collect()
}

///|
/// Returns the single ``(source, target)`` pair registered under
/// ``label``.
pub fn PathwayGraph::edges(self : PathwayGraph, label : String) -> (String, String) {
  self.edge_map[label]
}

///|
/// All edge labels that leave ``node``.
pub fn PathwayGraph::children(self : PathwayGraph, node : String) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  let keys = self.edge_map.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let lbl = keys[i]
    let (from, _to) = self.edge_map[lbl]
    if from == node {
      result.push(lbl)
    }
    i = i + 1
  }
  result
}

///|
/// All edge labels that arrive at ``node``.
pub fn PathwayGraph::parents(self : PathwayGraph, node : String) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  let keys = self.edge_map.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let lbl = keys[i]
    let (_from, to) = self.edge_map[lbl]
    if to == node {
      result.push(lbl)
    }
    i = i + 1
  }
  result
}

///|
/// ``(target, edge_label)`` pairs for every edge leaving ``node``.
pub fn PathwayGraph::child_edges(self : PathwayGraph, node : String) -> Array[(String, String)] {
  let result : Array[(String, String)] = Array(capacity=0)
  let keys = self.edge_map.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let lbl = keys[i]
    let (from, to) = self.edge_map[lbl]
    if from == node {
      result.push((to, lbl))
    }
    i = i + 1
  }
  result
}

///|
/// ``(source, edge_label)`` pairs for every edge arriving at ``node``.
pub fn PathwayGraph::parent_edges(self : PathwayGraph, node : String) -> Array[(String, String)] {
  let result : Array[(String, String)] = Array(capacity=0)
  let keys = self.edge_map.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let lbl = keys[i]
    let (from, to) = self.edge_map[lbl]
    if to == node {
      result.push((from, lbl))
    }
    i = i + 1
  }
  result
}

// ===========================================================================
// PathwayMultiGraph  (directed multigraph — parallel edges allowed)
// ===========================================================================

///|
/// A directed multigraph that keeps parallel edges distinct.
///
/// Unlike ``PathwayGraph``, multiple edges with the same label (or
/// different labels between the same pair of nodes) are preserved.
pub struct PathwayMultiGraph {
  nodes : Map[String, Bool]
  /// label -> list of (source, target)  — multiple entries allowed
  edge_map : Map[String, Array[(String, String)]]
}

///|
/// Create an empty PathwayMultiGraph.
pub fn PathwayMultiGraph::new() -> PathwayMultiGraph {
  PathwayMultiGraph::{ nodes: Map([]), edge_map: Map([]) }
}

///|
/// Register a node (no-op if already present).
pub fn PathwayMultiGraph::add_node(self : PathwayMultiGraph, label : String) -> Unit {
  self.nodes[label] = true
}

fn _pmg_get_edge_list(m : Map[String, Array[(String, String)]], label : String) -> Array[(String, String)] {
  match m.get(label) {
    Some(arr) => arr
    None => []
  }
}

///|
/// Add a directed edge ``from -> to`` labelled by ``label``.
/// Duplicate labels are **allowed** (PathwayMultiGraph invariant).
pub fn PathwayMultiGraph::add_edge(self : PathwayMultiGraph, from_node : String, to_node : String, label : String) -> Unit {
  self.add_node(from_node)
  self.add_node(to_node)
  let list = _pmg_get_edge_list(self.edge_map, label)
  list.push((from_node, to_node))
  self.edge_map[label] = list
}

///|
/// The set of node labels.
pub fn PathwayMultiGraph::nodes(self : PathwayMultiGraph) -> Array[String] {
  self.nodes.keys().collect()
}

///|
/// The set of **distinct** edge labels.
pub fn PathwayMultiGraph::labels(self : PathwayMultiGraph) -> Array[String] {
  self.edge_map.keys().collect()
}

///|
/// All ``(source, target)`` pairs registered under ``label``.
/// Returns an Array because PathwayMultiGraph allows parallel edges
/// with the same label.
pub fn PathwayMultiGraph::edges(self : PathwayMultiGraph, label : String) -> Array[(String, String)] {
  _pmg_get_edge_list(self.edge_map, label)
}

///|
/// All distinct edge labels that leave ``node``.
pub fn PathwayMultiGraph::children(self : PathwayMultiGraph, node : String) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  let seen : Map[String, Bool] = Map([])
  let keys = self.edge_map.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let lbl = keys[i]
    let list = self.edge_map[lbl]
    let mut j = 0
    while j < list.length() {
      let (from, _to) = list[j]
      if from == node && !seen.contains(lbl) {
        seen[lbl] = true
        result.push(lbl)
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// All distinct edge labels that arrive at ``node``.
pub fn PathwayMultiGraph::parents(self : PathwayMultiGraph, node : String) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  let seen : Map[String, Bool] = Map([])
  let keys = self.edge_map.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let lbl = keys[i]
    let list = self.edge_map[lbl]
    let mut j = 0
    while j < list.length() {
      let (_from, to) = list[j]
      if to == node && !seen.contains(lbl) {
        seen[lbl] = true
        result.push(lbl)
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// ``(target, edge_label)`` pairs for every edge leaving ``node``.
/// For parallel edges this returns **one entry per physical edge**.
pub fn PathwayMultiGraph::child_edges(self : PathwayMultiGraph, node : String) -> Array[(String, String)] {
  let result : Array[(String, String)] = Array(capacity=0)
  let keys = self.edge_map.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let lbl = keys[i]
    let list = self.edge_map[lbl]
    let mut j = 0
    while j < list.length() {
      let (from, to) = list[j]
      if from == node {
        result.push((to, lbl))
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// ``(source, edge_label)`` pairs for every edge arriving at ``node``.
pub fn PathwayMultiGraph::parent_edges(self : PathwayMultiGraph, node : String) -> Array[(String, String)] {
  let result : Array[(String, String)] = Array(capacity=0)
  let keys = self.edge_map.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let lbl = keys[i]
    let list = self.edge_map[lbl]
    let mut j = 0
    while j < list.length() {
      let (from, to) = list[j]
      if to == node {
        result.push((from, lbl))
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

// ===========================================================================
// Adjacency helpers (used by BFS / DFS)
// ===========================================================================

fn _pg_neighbors(g : PathwayGraph, node : String) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  let keys = g.edge_map.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let (from, to) = g.edge_map[keys[i]]
    if from == node {
      result.push(to)
    }
    i = i + 1
  }
  result
}

fn _pmg_neighbors(g : PathwayMultiGraph, node : String) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  let keys = g.edge_map.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let list = g.edge_map[keys[i]]
    let mut j = 0
    while j < list.length() {
      let (from, to) = list[j]
      if from == node {
        result.push(to)
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

// ===========================================================================
// Breadth-first search  (Biopython: bf_search)
// ===========================================================================

///|
/// Breadth-first traversal of a ``PathwayMultiGraph`` rooted at
/// ``root``.  Returns the list of **node labels** in discovery order
/// (matching ``Bio.Pathway.Rep.MultiGraph.bf_search``).
pub fn bf_search(g : PathwayMultiGraph, root? : String) -> Array[String] {
  let all_nodes = g.nodes()
  if all_nodes.length() == 0 {
    return []
  }
  let start = match root {
    Some(r) => r
    None => all_nodes[0]
  }
  let seen : Map[String, Bool] = Map([])
  let result : Array[String] = Array(capacity=0)
  let queue : Array[String] = Array(capacity=0)

  seen[start] = true
  result.push(start)
  queue.push(start)

  let mut head = 0
  while head < queue.length() {
    let current = queue[head]
    head = head + 1
    let neighbors = _pmg_neighbors(g, current)
    let mut i = 0
    while i < neighbors.length() {
      let nb = neighbors[i]
      if !seen.contains(nb) {
        seen[nb] = true
        result.push(nb)
        queue.push(nb)
      }
      i = i + 1
    }
  }
  result
}

// ===========================================================================
// Depth-first search  (Biopython: df_search)
// ===========================================================================

///|
/// Depth-first traversal of a ``PathwayMultiGraph`` rooted at
/// ``root``.  Returns node labels in pre-order discovery order
/// (matching ``Bio.Pathway.Rep.MultiGraph.df_search``).
pub fn df_search(g : PathwayMultiGraph, root? : String) -> Array[String] {
  let all_nodes = g.nodes()
  if all_nodes.length() == 0 {
    return []
  }
  let start = match root {
    Some(r) => r
    None => all_nodes[0]
  }
  let seen : Map[String, Bool] = Map([])
  let result : Array[String] = Array(capacity=0)
  let stack : Array[String] = Array(capacity=0)

  seen[start] = true
  result.push(start)
  stack.push(start)

  while stack.length() > 0 {
    let current_opt = stack.pop()
    let current = match current_opt {
      Some(s) => s
      None => ""
    }
    let neighbors = _pmg_neighbors(g, current)
    let mut i = 0
    // Reverse so iteration order matches recursive pre-order DFS.
    while i < neighbors.length() {
      let nb = neighbors[neighbors.length() - 1 - i]
      if !seen.contains(nb) {
        seen[nb] = true
        result.push(nb)
        stack.push(nb)
      }
      i = i + 1
    }
  }
  result
}

///|
/// Overloaded BFS on ``PathwayGraph``.
pub fn bf_search_graph(g : PathwayGraph, root? : String) -> Array[String] {
  let all_nodes = g.nodes()
  if all_nodes.length() == 0 {
    return []
  }
  let start = match root {
    Some(r) => r
    None => all_nodes[0]
  }
  let seen : Map[String, Bool] = Map([])
  let result : Array[String] = Array(capacity=0)
  let queue : Array[String] = Array(capacity=0)

  seen[start] = true
  result.push(start)
  queue.push(start)

  let mut head = 0
  while head < queue.length() {
    let current = queue[head]
    head = head + 1
    let neighbors = _pg_neighbors(g, current)
    let mut i = 0
    while i < neighbors.length() {
      let nb = neighbors[i]
      if !seen.contains(nb) {
        seen[nb] = true
        result.push(nb)
        queue.push(nb)
      }
      i = i + 1
    }
  }
  result
}

///|
/// Overloaded DFS on ``PathwayGraph``.
pub fn df_search_graph(g : PathwayGraph, root? : String) -> Array[String] {
  let all_nodes = g.nodes()
  if all_nodes.length() == 0 {
    return []
  }
  let start = match root {
    Some(r) => r
    None => all_nodes[0]
  }
  let seen : Map[String, Bool] = Map([])
  let result : Array[String] = Array(capacity=0)
  let stack : Array[String] = Array(capacity=0)

  seen[start] = true
  result.push(start)
  stack.push(start)

  while stack.length() > 0 {
    let current_opt = stack.pop()
    let current = match current_opt {
      Some(s) => s
      None => ""
    }
    let neighbors = _pg_neighbors(g, current)
    let mut i = 0
    while i < neighbors.length() {
      let nb = neighbors[neighbors.length() - 1 - i]
      if !seen.contains(nb) {
        seen[nb] = true
        result.push(nb)
        stack.push(nb)
      }
      i = i + 1
    }
  }
  result
}
