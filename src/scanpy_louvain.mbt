///|
/// scanpy `tl.louvain` — Louvain modularity-based community detection on the
/// sparse `pp.neighbors` graph. Implements the local-moving pass (one level)
/// with the standard modularity gain of Blondel et al. (2008); the result is
/// relabelled to dense 0..k−1.

///|
pub fn scanpy_louvain(graph : SparseGraph) -> Array[Int] {
  let n = graph.n
  if n == 0 {
    return []
  }
  // weighted degrees (no self-loops: `adj` stores only off-diagonal edges)
  let deg = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    for nb in graph.adj[i] {
      deg[i] = deg[i] + nb.1
    }
  }
  let mut m = 0.0
  for d in deg {
    m = m + d
  }
  m = m / 2.0
  if m <= 0.0 {
    return Array::make(n, 0)
  }

  let comm = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    comm[i] = i
  }
  let stot : Map[Int, Double] = Map([], capacity=32)
  for i = 0; i < n; i = i + 1 {
    stot.set(i, deg[i])
  }
  let mut next_id = n

  let mut changed = true
  while changed {
    changed = false
    for u = 0; u < n; u = u + 1 {
      let k_u = deg[u]
      let cu = comm[u]
      // k_u,in per neighbour community
      let kin : Map[Int, Double] = Map([], capacity=16)
      for nb in graph.adj[u] {
        let cv = comm[nb.0]
        match kin.get(cv) {
          None => kin.set(cv, nb.1)
          Some(x) => kin.set(cv, x + nb.1)
        }
      }
      let kin_a = match kin.get(cu) {
        None => 0.0
        Some(x) => x
      }
      let stot_cu = match stot.get(cu) {
        None => 0.0
        Some(x) => x
      }
      let mut best_gain = 0.0
      let mut best_c = cu
      for c in kin {
        let key = c.0
        if key == cu {
          continue
        }
        let stot_c = match stot.get(key) {
          None => 0.0
          Some(x) => x
        }
        let gain = (c.1 - kin_a) / m -
          k_u * (stot_c - stot_cu + k_u) / (2.0 * m * m)
        if gain > best_gain {
          best_gain = gain
          best_c = key
        }
      }
      // leave current community and become a singleton
      let gain_single = (0.0 - kin_a) / m +
        k_u * (stot_cu - k_u) / (2.0 * m * m)
      if gain_single > best_gain {
        best_gain = gain_single
        best_c = -1
      }

      if best_c != cu {
        stot.set(cu, stot_cu - k_u)
        let target = if best_c == -1 {
          let t = next_id
          next_id = next_id + 1
          t
        } else {
          best_c
        }
        match stot.get(target) {
          None => stot.set(target, k_u)
          Some(x) => stot.set(target, x + k_u)
        }
        comm[u] = target
        changed = true
      }
    }
  }

  // relabel to dense 0..k−1
  let mapping : Map[Int, Int] = Map([], capacity=32)
  let mut nc = 0
  for u = 0; u < n; u = u + 1 {
    match mapping.get(comm[u]) {
      None => {
        mapping.set(comm[u], nc)
        nc = nc + 1
      }
      Some(_) => ()
    }
  }
  let out = Array::make(n, 0)
  for u = 0; u < n; u = u + 1 {
    out[u] = mapping.get(comm[u]).unwrap()
  }
  out
}

///|
/// Leiden-style connectivity refinement: run Louvain local moving, then split
/// any disconnected community (in the induced subgraph) into its connected
/// components, guaranteeing every returned community is connected.
pub fn scanpy_leiden(graph : SparseGraph) -> Array[Int] {
  let n = graph.n
  let base = scanpy_louvain(graph)
  let out = Array::make(n, 0)
  let visited = Array::make(n, false)
  let mut next_label = 0
  for seed = 0; seed < n; seed = seed + 1 {
    if visited[seed] {
      continue
    }
    let comm = base[seed]
    let frontier : Array[Int] = Array(capacity=0)
    frontier.push(seed)
    visited[seed] = true
    let mut head = 0
    while head < frontier.length() {
      let u = frontier[head]
      head = head + 1
      out[u] = next_label
      for nb in graph.adj[u] {
        let v = nb.0
        if !(visited[v]) && base[v] == comm {
          visited[v] = true
          frontier.push(v)
        }
      }
    }
    next_label = next_label + 1
  }
  out
}