///|
/// Overlap-Layout-Consensus (OLC) algorithm for sequence assembly.
///
/// References:
/// - Myers EW, Miller W (1988). "Optimal alignments in linear space."
///   Comput Appl Biosci. 4(1):11-17.
/// - Pevzner PA (2000). "Computational Molecular Biology: An Algorithmic Approach."
///   MIT Press.

// ---------------------------------------------------------------------------
// Helper Functions
// ---------------------------------------------------------------------------

///|
fn sa_copy_array(arr : Array[Int]) -> Array[Int] {
  let result : Array[Int] = Array::new()
  for i = 0; i < arr.length(); i = i + 1 {
    result.push(arr[i])
  }
  result
}

// ---------------------------------------------------------------------------
// Overlap Detection
// ---------------------------------------------------------------------------

///|
pub struct Overlap {
  from : Int
  to : Int
  overlap_len : Int
}

///|
pub fn Overlap::new(from : Int, to : Int, overlap_len : Int) -> Overlap {
  Overlap::{ from, to, overlap_len }
}

///|
pub fn Overlap::from(self : Overlap) -> Int {
  self.from
}

///|
pub fn Overlap::to(self : Overlap) -> Int {
  self.to
}

///|
pub fn Overlap::overlap_len(self : Overlap) -> Int {
  self.overlap_len
}

///|
pub fn compute_overlap(a : String, b : String, min_overlap : Int) -> Int {
  let a_len = a.length()
  let b_len = b.length()
  let max_possible = if a_len < b_len { a_len } else { b_len }
  if max_possible < min_overlap {
    return 0
  }
  for l = max_possible; l >= min_overlap; l = l - 1 {
    let suffix = a[a_len - l:a_len].to_owned()
    let prefix = b[0:l].to_owned()
    if suffix == prefix {
      return l
    }
  }
  0
}

///|
pub fn compute_all_overlaps(
  reads : Array[String],
  min_overlap : Int,
) -> Array[Overlap] {
  let n = reads.length()
  let overlaps : Array[Overlap] = Array::new()
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      if i == j {
        continue
      }
      let overlap_len = compute_overlap(reads[i], reads[j], min_overlap)
      if overlap_len >= min_overlap {
        overlaps.push(Overlap::new(i, j, overlap_len))
      }
    }
  }
  overlaps
}

///|
pub fn compute_overlap_graph(
  reads : Array[String],
  min_overlap : Int,
) -> Map[Int, Array[(Int, Int)]] {
  let n = reads.length()
  let graph : Map[Int, Array[(Int, Int)]] = Map([], capacity=n)
  for i = 0; i < n; i = i + 1 {
    let edges : Array[(Int, Int)] = Array::new()
    for j = 0; j < n; j = j + 1 {
      if i == j {
        continue
      }
      let overlap_len = compute_overlap(reads[i], reads[j], min_overlap)
      if overlap_len >= min_overlap {
        edges.push((j, overlap_len))
      }
    }
    if edges.length() > 0 {
      graph[i] = edges
    }
  }
  graph
}

// ---------------------------------------------------------------------------
// Layout: Find Hamiltonian Path in Overlap Graph
// ---------------------------------------------------------------------------

///|
pub fn find_layout_path(reads : Array[String], min_overlap : Int) -> Array[Int] {
  let n = reads.length()
  if n == 0 {
    return Array::new()
  }
  if n == 1 {
    let path : Array[Int] = Array::new()
    path.push(0)
    return path
  }
  let graph = compute_overlap_graph(reads, min_overlap)
  find_hamiltonian_path(graph, n)
}

///|
fn find_hamiltonian_path(
  graph : Map[Int, Array[(Int, Int)]],
  n : Int,
) -> Array[Int] {
  let visited : Array[Bool] = Array::new()
  for i = 0; i < n; i = i + 1 {
    visited.push(false)
  }
  for start = 0; start < n; start = start + 1 {
    let path : Array[Int] = Array::new()
    let (found, result_path) = dfs_hamiltonian(graph, start, visited, path, n)
    if found {
      return result_path
    }
  }
  let path : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    path.push(i)
  }
  path
}

///|
fn dfs_hamiltonian(
  graph : Map[Int, Array[(Int, Int)]],
  current : Int,
  visited : Array[Bool],
  path : Array[Int],
  n : Int,
) -> (Bool, Array[Int]) {
  visited[current] = true
  let new_path = sa_copy_array(path)
  new_path.push(current)
  if new_path.length() == n {
    visited[current] = false
    return (true, new_path)
  }
  match graph.get(current) {
    Some(neighbors) =>
      for i = 0; i < neighbors.length(); i = i + 1 {
        let (next, _) = neighbors[i]
        if !visited[next] {
          let (found, result_path) = dfs_hamiltonian(
            graph, next, visited, new_path, n,
          )
          if found {
            visited[current] = false
            return (true, result_path)
          }
        }
      }
    None => ()
  }
  visited[current] = false
  (false, Array::new())
}

// ---------------------------------------------------------------------------
// Consensus Generation
// ---------------------------------------------------------------------------

///|
pub fn generate_consensus(reads : Array[String], path : Array[Int]) -> String {
  if path.length() == 0 {
    return ""
  }
  let consensus = StringBuilder::new()
  consensus.write_string(reads[path[0]])
  for i = 1; i < path.length(); i = i + 1 {
    let prev_idx = path[i - 1]
    let curr_idx = path[i]
    let overlap_len = compute_overlap(reads[prev_idx], reads[curr_idx], 1)
    let curr_str = reads[curr_idx]
    if overlap_len < curr_str.length() {
      consensus.write_string(curr_str[overlap_len:curr_str.length()].to_owned())
    }
  }
  consensus.to_string()
}

///|
pub fn olc_assemble(reads : Array[String], min_overlap : Int) -> String {
  let path = find_layout_path(reads, min_overlap)
  generate_consensus(reads, path)
}

// ---------------------------------------------------------------------------
// Overlap Graph Visualization
// ---------------------------------------------------------------------------

///|
pub fn overlap_graph_to_dot(reads : Array[String], min_overlap : Int) -> String {
  let graph = compute_overlap_graph(reads, min_overlap)
  let buf = StringBuilder::new()
  buf.write_string("digraph OverlapGraph {\n")
  buf.write_string("  node [shape=ellipse];\n")
  for i = 0; i < reads.length(); i = i + 1 {
    buf.write_string("  node")
    buf.write_string(i.to_string())
    buf.write_string(" [label=\"")
    let label = if reads[i].length() > 10 {
      reads[i][0:10].to_owned() + "..."
    } else {
      reads[i]
    }
    buf.write_string(label)
    buf.write_string("\"];\n")
  }
  let keys = graph.keys().collect()
  for i = 0; i < keys.length(); i = i + 1 {
    let from = keys[i]
    match graph.get(from) {
      Some(edges) =>
        for j = 0; j < edges.length(); j = j + 1 {
          let (to, overlap_len) = edges[j]
          buf.write_string("  node")
          buf.write_string(from.to_string())
          buf.write_string(" -> node")
          buf.write_string(to.to_string())
          buf.write_string(" [label=\"")
          buf.write_string(overlap_len.to_string())
          buf.write_string("\"];\n")
        }
      None => ()
    }
  }
  buf.write_string("}\n")
  buf.to_string()
}

// ---------------------------------------------------------------------------
// Error Correction
// ---------------------------------------------------------------------------

///|
pub fn correct_errors(
  reads : Array[String],
  min_overlap : Int,
) -> Array[String] {
  let n = reads.length()
  let corrected : Array[String] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let read = reads[i]
    let corrected_read = correct_read(reads, i, min_overlap)
    corrected.push(corrected_read)
  }
  corrected
}

///|
fn correct_read(reads : Array[String], idx : Int, min_overlap : Int) -> String {
  let n = reads.length()
  let read = reads[idx]
  let overlaps : Array[(String, Int)] = Array::new()
  for i = 0; i < n; i = i + 1 {
    if i == idx {
      continue
    }
    let forward_overlap = compute_overlap(read, reads[i], min_overlap)
    let backward_overlap = compute_overlap(reads[i], read, min_overlap)
    if forward_overlap >= min_overlap {
      overlaps.push((reads[i], forward_overlap))
    }
    if backward_overlap >= min_overlap {
      overlaps.push((reads[i], -backward_overlap))
    }
  }
  if overlaps.length() == 0 {
    return read
  }
  read
}
