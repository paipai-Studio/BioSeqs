///|
/// De Bruijn Graph for sequence assembly.
///
/// References:
/// - Pevzner PA, Tang H, Waterman MS (2001). "An Eulerian path approach to
///   DNA fragment assembly." Proc Natl Acad Sci U S A. 98(17):9748-53.
/// - Compeau P, Pevzner P, Tesler G (2011). "How to apply de Bruijn graphs to
///   genome assembly." Nature Biotechnology. 29(11):987-991.

// ---------------------------------------------------------------------------
// De Bruijn Graph using adjacency list representation
// ---------------------------------------------------------------------------

///|
/// De Bruijn graph using adjacency list representation.
///
/// Nodes represent (k-1)-mers.
/// Edges represent k-mers: edge from "ACG" to "CGT" represents k-mer "ACGT".
pub struct DeBruijnGraph {
  k : Int
  edges : Map[String, Array[String]]
  node_info : Map[String, (Int, Int)]
  mut edge_count : Int
}

///|
pub fn DeBruijnGraph::new(k : Int) -> DeBruijnGraph {
  DeBruijnGraph::{
    k,
    edges: Map([], capacity=0),
    node_info: Map([], capacity=0),
    edge_count: 0,
  }
}

///|
pub fn DeBruijnGraph::k(self : DeBruijnGraph) -> Int {
  self.k
}

///|
pub fn DeBruijnGraph::node_count(self : DeBruijnGraph) -> Int {
  self.node_info.size()
}

///|
pub fn DeBruijnGraph::edge_count(self : DeBruijnGraph) -> Int {
  self.edge_count
}

///|
pub fn DeBruijnGraph::add_kmer(self : DeBruijnGraph, kmer : String) -> Unit {
  if kmer.length() < self.k {
    return
  }
  let prefix = kmer[0:self.k - 1].to_owned()
  let suffix = kmer[1:self.k].to_owned()
  match self.node_info.get(prefix) {
    None => self.node_info[prefix] = (0, 0)
    Some(_) => ()
  }
  match self.node_info.get(suffix) {
    None => self.node_info[suffix] = (0, 0)
    Some(_) => ()
  }
  match self.edges.get(prefix) {
    None => {
      self.edges[prefix] = [suffix]
      match self.node_info.get(prefix) {
        Some((in_d, out_d)) => self.node_info[prefix] = (in_d, out_d + 1)
        None => ()
      }
      match self.node_info.get(suffix) {
        Some((in_d, out_d)) => self.node_info[suffix] = (in_d + 1, out_d)
        None => ()
      }
      self.edge_count = self.edge_count + 1
    }
    Some(list) => {
      let mut exists = false
      for i = 0; i < list.length(); i = i + 1 {
        if list[i] == suffix {
          exists = true
          break
        }
      }
      if not(exists) {
        let new_list = list
        new_list.push(suffix)
        self.edges[prefix] = new_list
        match self.node_info.get(prefix) {
          Some((in_d, out_d)) => self.node_info[prefix] = (in_d, out_d + 1)
          None => ()
        }
        match self.node_info.get(suffix) {
          Some((in_d, out_d)) => self.node_info[suffix] = (in_d + 1, out_d)
          None => ()
        }
        self.edge_count = self.edge_count + 1
      }
    }
  }
}

///|
pub fn DeBruijnGraph::add_sequence(self : DeBruijnGraph, seq : String) -> Unit {
  let n = seq.length()
  if n < self.k {
    return
  }
  let total = n - self.k + 1
  for i = 0; i < total; i = i + 1 {
    let kmer = seq[i:i + self.k].to_owned()
    self.add_kmer(kmer)
  }
}

///|
pub fn DeBruijnGraph::add_sequences(
  self : DeBruijnGraph,
  seqs : Array[String],
) -> Unit {
  for i = 0; i < seqs.length(); i = i + 1 {
    self.add_sequence(seqs[i])
  }
}

///|
pub fn DeBruijnGraph::contains_node(
  self : DeBruijnGraph,
  kmer : String,
) -> Bool {
  match self.node_info.get(kmer) {
    Some(_) => true
    None => false
  }
}

///|
pub fn DeBruijnGraph::in_degree(self : DeBruijnGraph, kmer : String) -> Int {
  match self.node_info.get(kmer) {
    Some((in_d, _)) => in_d
    None => 0
  }
}

///|
pub fn DeBruijnGraph::out_degree(self : DeBruijnGraph, kmer : String) -> Int {
  match self.node_info.get(kmer) {
    Some((_, out_d)) => out_d
    None => 0
  }
}

///|
pub fn DeBruijnGraph::outgoing(
  self : DeBruijnGraph,
  kmer : String,
) -> Array[String] {
  match self.edges.get(kmer) {
    Some(list) => list
    None => []
  }
}

///|
pub fn DeBruijnGraph::all_nodes(self : DeBruijnGraph) -> Array[String] {
  self.node_info.keys().collect()
}

///|
pub fn DeBruijnGraph::start_nodes(self : DeBruijnGraph) -> Array[String] {
  let result : Array[String] = Array::new()
  let keys = self.node_info.keys().collect()
  for i = 0; i < keys.length(); i = i + 1 {
    match self.node_info.get(keys[i]) {
      Some((in_d, _)) => if in_d == 0 { result.push(keys[i]) }
      None => ()
    }
  }
  result
}

///|
pub fn DeBruijnGraph::end_nodes(self : DeBruijnGraph) -> Array[String] {
  let result : Array[String] = Array::new()
  let keys = self.node_info.keys().collect()
  for i = 0; i < keys.length(); i = i + 1 {
    match self.node_info.get(keys[i]) {
      Some((_, out_d)) => if out_d == 0 { result.push(keys[i]) }
      None => ()
    }
  }
  result
}

///|
pub fn DeBruijnGraph::branch_nodes(self : DeBruijnGraph) -> Array[String] {
  let result : Array[String] = Array::new()
  let keys = self.node_info.keys().collect()
  for i = 0; i < keys.length(); i = i + 1 {
    match self.node_info.get(keys[i]) {
      Some((_, out_d)) => if out_d > 1 { result.push(keys[i]) }
      None => ()
    }
  }
  result
}

///|
pub fn DeBruijnGraph::hierholzer(
  self : DeBruijnGraph,
  start : String,
) -> Array[String] {
  let current_path : Array[String] = Array::new()
  let circuit : Array[String] = Array::new()
  let available : Map[String, Array[String]] = Map([], capacity=0)
  let edge_keys = self.edges.keys().collect()
  for i = 0; i < edge_keys.length(); i = i + 1 {
    match self.edges.get(edge_keys[i]) {
      Some(list) => available[edge_keys[i]] = list.clone()
      None => ()
    }
  }
  current_path.push(start)
  while current_path.length() > 0 {
    let current = current_path[current_path.length() - 1]
    match available.get(current) {
      Some(list) =>
        if list.length() > 0 {
          let next_node = list[list.length() - 1]
          let new_list = list[0:list.length() - 1].to_owned()
          available[current] = new_list
          current_path.push(next_node)
        } else {
          let popped = current_path.pop()
          match popped {
            Some(v) => circuit.push(v)
            None => ()
          }
        }
      None => {
        let popped = current_path.pop()
        match popped {
          Some(v) => circuit.push(v)
          None => ()
        }
      }
    }
  }
  let result : Array[String] = Array::new()
  for i = circuit.length() - 1; i >= 0; i = i - 1 {
    result.push(circuit[i])
  }
  result
}

///|
pub fn DeBruijnGraph::find_eulerian_path(self : DeBruijnGraph) -> String? {
  let mut start_candidate : String? = None
  let mut end_candidate : String? = None
  let keys = self.node_info.keys().collect()
  for i = 0; i < keys.length(); i = i + 1 {
    match self.node_info.get(keys[i]) {
      Some((in_d, out_d)) => {
        let diff = out_d - in_d
        if diff == 1 {
          match start_candidate {
            None => start_candidate = Some(keys[i])
            Some(_) => return None
          }
        } else if diff == -1 {
          match end_candidate {
            None => end_candidate = Some(keys[i])
            Some(_) => return None
          }
        } else if diff != 0 {
          return None
        }
      }
      None => ()
    }
  }
  let start = match start_candidate {
    Some(s) => s
    None => {
      let mut found : String? = None
      let keys2 = self.node_info.keys().collect()
      for j = 0; j < keys2.length(); j = j + 1 {
        match self.node_info.get(keys2[j]) {
          Some((_, out_d)) =>
            if out_d > 0 {
              found = Some(keys2[j])
              break
            }
          None => ()
        }
      }
      match found {
        Some(s) => s
        None => return None
      }
    }
  }
  let path = self.hierholzer(start)
  if path.length() == 0 {
    return None
  }
  let result = StringBuilder::new()
  result.write_string(path[0])
  for i = 1; i < path.length(); i = i + 1 {
    let node = path[i]
    if node.length() > 0 {
      result.write_char(node.unsafe_get(node.length() - 1).unsafe_to_char())
    }
  }
  Some(result.to_string())
}

///|
pub fn DeBruijnGraph::to_string(self : DeBruijnGraph) -> String {
  let buf = StringBuilder::new()
  buf.write_string("De Bruijn Graph (k=")
  buf.write_string(self.k.to_string())
  buf.write_string("):\n")
  buf.write_string("  Nodes: ")
  buf.write_string(self.node_count().to_string())
  buf.write_string("\n")
  buf.write_string("  Edges: ")
  buf.write_string(self.edge_count().to_string())
  buf.write_string("\n")
  buf.write_string("  Start nodes: ")
  let starts = self.start_nodes()
  for i = 0; i < starts.length(); i = i + 1 {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(starts[i])
  }
  buf.write_string("\n")
  buf.write_string("  End nodes: ")
  let ends = self.end_nodes()
  for i = 0; i < ends.length(); i = i + 1 {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(ends[i])
  }
  buf.write_string("\n")
  buf.write_string("  Branch nodes: ")
  let branches = self.branch_nodes()
  for i = 0; i < branches.length(); i = i + 1 {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(branches[i])
  }
  buf.write_string("\n")
  buf.write_string("  Adjacency list:\n")
  let edge_keys = self.edges.keys().collect()
  for i = 0; i < edge_keys.length(); i = i + 1 {
    buf.write_string("    ")
    buf.write_string(edge_keys[i])
    buf.write_string(" -> ")
    match self.edges.get(edge_keys[i]) {
      Some(list) =>
        for j = 0; j < list.length(); j = j + 1 {
          if j > 0 {
            buf.write_string(", ")
          }
          buf.write_string(list[j])
        }
      None => ()
    }
    buf.write_char('\n')
  }
  buf.to_string()
}

///|
pub fn assemble_sequences(sequences : Array[String], k : Int) -> String? {
  let graph = DeBruijnGraph::new(k)
  graph.add_sequences(sequences)
  graph.find_eulerian_path()
}

///|
pub fn assemble_single(seq : String, k : Int) -> String? {
  let graph = DeBruijnGraph::new(k)
  graph.add_sequence(seq)
  graph.find_eulerian_path()
}
