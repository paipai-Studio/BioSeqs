///|
/// Louvain community detection, full multilevel scheme (Blondel et al. 2008).
/// The first pass moves individual vertices between communities until the
/// modularity gain is exhausted; the communities are then aggregated into
/// super-vertices (self-loops = twice the internal edge weight, preserving
/// total degree) and the local moving is repeated on the coarsened graph.
/// Iterating aggregation/moving until no merge is possible is what
/// distinguishes this from the single-pass `scanpy_louvain`.

///|
/// Result of the multilevel Louvain: the flattened partition (dense labels
/// `0..k-1`), the number of communities, the number of aggregation levels
/// performed and the final weighted modularity (with the `resolution`
/// multiplier).
pub struct ScanpyLouvainMultiLevel {
  partition : Array[Int32]
  n_communities : Int32
  n_levels : Int32
  modularity : Double
}

///|
/// Weighted modularity of a partition on an undirected weighted graph:
/// `Q = sum_c [ k_in_c / (2m) - resolution * (Sigma_tot_c / (2m))^2 ]` where
/// `k_in_c` counts twice the internal edge weight (self-loops included once)
/// and `m` is the total edge weight.
pub fn scanpy_modularity(
  graph : SparseGraph,
  partition : Array[Int32],
  resolution~ : Double = 1.0,
) -> Double {
  let n = graph.n
  let deg = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let mut d = if i < graph.diag.length() { graph.diag[i] } else { 0.0 }
    for nb in graph.adj[i] {
      d = d + nb.1
    }
    deg[i] = d
  }
  let mut m2 = 0.0
  for d in deg {
    m2 = m2 + d
  }
  if m2 <= 0.0 {
    return 0.0
  }
  // community sums: sigma_tot (degree) and k_in (2 x internal weight)
  let sigma_tot : Map[Int32, Double] = Map([], capacity=16)
  let k_in : Map[Int32, Double] = Map([], capacity=16)
  for i = 0; i < n; i = i + 1 {
    let c = partition[i]
    match sigma_tot.get(c) {
      None => sigma_tot.set(c, deg[i])
      Some(x) => sigma_tot.set(c, x + deg[i])
    }
    match k_in.get(c) {
      None =>
        k_in.set(c, if i < graph.diag.length() { graph.diag[i] } else { 0.0 })
      Some(x) => if i < graph.diag.length() { k_in.set(c, x + graph.diag[i]) }
    }
  }
  for i = 0; i < n; i = i + 1 {
    let ci = partition[i]
    for nb in graph.adj[i] {
      if partition[nb.0] == ci {
        // each undirected edge is stored in both directions: accumulate once
        // per direction, halve at the end via the factor below (x2 -> /2)
        match k_in.get(ci) {
          None => k_in.set(ci, nb.1)
          Some(x) => k_in.set(ci, x + nb.1)
        }
      }
    }
  }
  let mut q = 0.0
  for kv in k_in {
    let c = kv.0
    // k_in already counts both directed copies (= 2 * internal weight), which
    // is exactly the numerator of L_c/m once divided by 2m = m2.
    let kin = kv.1
    let stot = match sigma_tot.get(c) {
      None => 0.0
      Some(x) => x
    }
    q = q + kin / m2 - resolution * (stot / m2) * (stot / m2)
  }
  q
}

///|
/// One full local-moving pass on the graph `(n, diag, adj)`: repeatedly move
/// each vertex to the neighbouring community with the largest positive
/// modularity gain until a sweep changes nothing.  The gain of joining a
/// community `c` after removing the vertex from its own is
/// `k_uc - resolution * k_u * Sigma_c / (2m)` (all edge weights, `2m` the sum
/// of degrees); staying behind as a singleton is the zero-gain baseline.  That
/// formulation is monotone, so the sweep cannot oscillate between two
/// communities.  Returns the community id per vertex (not necessarily dense).
fn louvain_ml_local_moving(
  n : Int32,
  diag : Array[Double],
  adj : Array[Array[(Int32, Double)]],
  resolution : Double,
) -> Array[Int32] {
  let deg = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    deg[i] = diag[i]
    for nb in adj[i] {
      deg[i] = deg[i] + nb.1
    }
  }
  let mut m2 = 0.0
  for d in deg {
    m2 = m2 + d
  }
  if m2 <= 0.0 {
    return Array::make(n, 0)
  }
  let comm = Array::make(n, 0)
  let stot : Map[Int32, Double] = Map([], capacity=n)
  for i = 0; i < n; i = i + 1 {
    comm[i] = i
    stot.set(i, deg[i])
  }
  let mut changed = true
  let mut sweeps = 0
  while changed && sweeps < 200 {
    changed = false
    sweeps = sweeps + 1
    for u = 0; u < n; u = u + 1 {
      let k_u = deg[u]
      let cu = comm[u]
      // temporarily remove u from its community
      let stot_cu = match stot.get(cu) {
        None => 0.0
        Some(x) => x
      }
      stot.set(cu, stot_cu - k_u)
      // weight from u to every neighbouring community (self-loop belongs to
      // whichever community u ends up in)
      let kin : Map[Int32, Double] = Map([], capacity=16)
      for nb in adj[u] {
        let cv = comm[nb.0]
        match kin.get(cv) {
          None => kin.set(cv, nb.1)
          Some(x) => kin.set(cv, x + nb.1)
        }
      }
      if diag[u] > 0.0 {
        match kin.get(cu) {
          None => kin.set(cu, diag[u])
          Some(x) => kin.set(cu, x + diag[u])
        }
      }
      // best target: stay as a singleton (gain 0) or join a neighbour
      let mut best_c = cu
      let mut best_gain = 0.0
      for kv in kin {
        let stot_c = match stot.get(kv.0) {
          None => 0.0
          Some(x) => x
        }
        let gain = kv.1 - resolution * k_u * stot_c / m2
        if gain > best_gain + 1.0e-12 {
          best_gain = gain
          best_c = kv.0
        }
      }
      let stot_b = match stot.get(best_c) {
        None => 0.0
        Some(x) => x
      }
      stot.set(best_c, stot_b + k_u)
      if best_c != cu {
        comm[u] = best_c
        changed = true
      }
    }
  }
  comm
}

///|
/// Aggregate a graph into its communities: the new self-loop weight of a
/// community is twice its internal edge weight plus the sum of its members'
/// self-loops (degree-preserving), and the new edge weight between two
/// communities is the sum of all inter-community edge weights.
fn louvain_ml_aggregate(
  n : Int32,
  diag : Array[Double],
  adj : Array[Array[(Int32, Double)]],
  comm : Array[Int32],
  n_comm : Int32,
) -> (Array[Double], Array[Array[(Int32, Double)]]) {
  let new_diag = Array::make(n_comm, 0.0)
  let edges : Map[Int32, Double] = Map([], capacity=n)
  for i = 0; i < n; i = i + 1 {
    let ci = comm[i]
    new_diag[ci] = new_diag[ci] + diag[i]
    for nb in adj[i] {
      let cj = comm[nb.0]
      if cj == ci {
        // internal edge appears twice (both directions) -> 2 x weight
        new_diag[ci] = new_diag[ci] + nb.1
      } else {
        let key = ci * n_comm + cj
        match edges.get(key) {
          None => edges.set(key, nb.1)
          Some(x) => edges.set(key, x + nb.1)
        }
      }
    }
  }
  let new_adj : Array[Array[(Int32, Double)]] = Array(capacity=n_comm)
  for c = 0; c < n_comm; c = c + 1 {
    new_adj.push([])
  }
  for e in edges {
    let a = e.0 / n_comm
    let b = e.0 % n_comm
    new_adj[a].push((b, e.1))
  }
  (new_diag, new_adj)
}

///|
/// Relabel community ids to dense `0..k-1` in first-appearance order.
fn louvain_ml_dense_relabel(comm : Array[Int32]) -> (Array[Int32], Int32) {
  let n = comm.length()
  let mapping : Map[Int32, Int32] = Map([], capacity=16)
  let mut nc = 0
  let out = Array::make(n, 0)
  for u = 0; u < n; u = u + 1 {
    match mapping.get(comm[u]) {
      None => {
        mapping.set(comm[u], nc)
        out[u] = nc
        nc = nc + 1
      }
      Some(c) => out[u] = c
    }
  }
  (out, nc)
}

///|
/// Multilevel Louvain community detection (`tl.louvain` full scheme).  Local
/// moving and community aggregation are alternated until an aggregation step
/// produces no merge, yielding the classic Blondel et al. (2008) hierarchy
/// flattened back to the original vertices.  `resolution` scales the
/// modularity null term (higher values produce more, smaller communities).
pub fn scanpy_louvain_multilevel(
  graph : SparseGraph,
  resolution~ : Double = 1.0,
) -> ScanpyLouvainMultiLevel {
  let n = graph.n
  if n == 0 {
    return ScanpyLouvainMultiLevel::{
      partition: [],
      n_communities: 0,
      n_levels: 0,
      modularity: 0.0,
    }
  }
  let identity = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    identity[i] = i
  }
  let mut partition = identity
  let mut cur_diag = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    if i < graph.diag.length() {
      cur_diag[i] = graph.diag[i]
    }
  }
  let mut cur_adj : Array[Array[(Int32, Double)]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[(Int32, Double)] = Array(capacity=0)
    for nb in graph.adj[i] {
      row.push((nb.0, nb.1))
    }
    cur_adj.push(row)
  }
  let mut cur_n = n
  let mut levels = 0
  for _lvl = 0; _lvl < 64; _lvl = _lvl + 1 {
    let comm = louvain_ml_local_moving(cur_n, cur_diag, cur_adj, resolution)
    let (dense, n_comm) = louvain_ml_dense_relabel(comm)
    if n_comm >= cur_n {
      break
    }
    // map the running partition through this level's dense communities
    let mapped = Array::make(n, 0)
    for i = 0; i < n; i = i + 1 {
      mapped[i] = dense[partition[i]]
    }
    partition = mapped
    let (nd, na) = louvain_ml_aggregate(cur_n, cur_diag, cur_adj, dense, n_comm)
    cur_n = n_comm
    cur_diag = nd
    cur_adj = na
    levels = levels + 1
  }
  let (final_partition, n_comm) = louvain_ml_dense_relabel(partition)
  let q = scanpy_modularity(graph, final_partition, resolution~)
  ScanpyLouvainMultiLevel::{
    partition: final_partition,
    n_communities: n_comm,
    n_levels: levels,
    modularity: q,
  }
}
