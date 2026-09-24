///|
/// Sparse symmetric weighted graph for the scanpy neighborhood pipeline.
///
/// The `pp.neighbors` connectivities are naturally sparse (each cell links to
/// ~`n_neighbors` others), so a dense `n × n` matrix wastes `O(n²)` memory.
/// `SparseGraph` stores the diagonal and an undirected adjacency list (each
/// edge kept in both directions), giving `O(n·k)` memory.

///|
pub struct SparseGraph {
  n : Int
  diag : Array[Double]
  adj : Array[Array[(Int, Double)]]
} derive(Eq)
pub extend SparseGraph with Eq::{equal, not_equal}

///|
pub fn SparseGraph::num_vertices(self : SparseGraph) -> Int {
  self.n
}

///|
/// Weighted degree of vertex `i`: diagonal self-weight plus incident edge
/// weights.
pub fn SparseGraph::degree(self : SparseGraph, i : Int) -> Double {
  let mut s = self.diag[i]
  for v in self.adj[i] {
    s = s + v.1
  }
  s
}

///|
/// Materialize the full dense `n × n` matrix (diagonal + symmetric off-diagonal
/// weights). Mainly for tests and for downstream steps that require a dense
/// operator (e.g. the eigendecomposition inside `diffmap`).
pub fn SparseGraph::to_dense(self : SparseGraph) -> Array[Array[Double]] {
  let n = self.n
  let m : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    row[i] = self.diag[i]
    m.push(row)
  }
  for i = 0; i < n; i = i + 1 {
    for v in self.adj[i] {
      m[i][v.0] = v.1
    }
  }
  m
}

///|
/// Build a `SparseGraph` from a dense symmetric matrix, keeping only the
/// off-diagonal positive entries (the diagonal is stored separately).
pub fn sparse_graph_from_dense(m : Array[Array[Double]]) -> SparseGraph {
  let n = m.length()
  let diag : Array[Double] = Array(capacity=n)
  let adj : Array[Array[(Int, Double)]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    diag.push(if i < m[i].length() { m[i][i] } else { 0.0 })
    adj.push([])
  }
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      if i != j && m[i][j] > 0.0 {
        adj[i].push((j, m[i][j]))
      }
    }
  }
  SparseGraph::{ n, diag, adj, }
}