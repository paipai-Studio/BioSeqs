///|
/// scanpy `tl.draw_graph` layout backends (`layout='fr' | 'fa2' | 'kk'`).
/// The existing `scanpy_draw_graph` implements the Fruchterman-Reingold force
/// model on a dense matrix; this file adds the `SparseGraph`-native variants:
/// a weighted FR / ForceAtlas2-style layout (`fr` / `fa2`, attraction along
/// graph edges with edge weights, repulsion over all pairs) and a
/// Kamada-Kawai stress layout (`kk`) driven by all-pairs shortest paths.

///|
/// All-pairs weighted shortest paths (plain Dijkstra, O(n^2) per source).
/// Unreachable pairs carry `1.0e30`.
fn scanpy_graph_shortest_paths(graph : SparseGraph) -> Array[Array[Double]] {
  let n = graph.n
  let out : Array[Array[Double]] = Array(capacity=n)
  for src = 0; src < n; src = src + 1 {
    let dist = Array::make(n, 1.0e30)
    let done = Array::make(n, false)
    dist[src] = 0.0
    for _t = 0; _t < n; _t = _t + 1 {
      // pick the unfinished vertex with minimal distance
      let mut u = -1
      let mut best = 1.0e30
      for v = 0; v < n; v = v + 1 {
        if !done[v] && dist[v] < best {
          best = dist[v]
          u = v
        }
      }
      if u < 0 {
        break
      }
      done[u] = true
      if dist[u] >= 1.0e30 {
        continue
      }
      for nb in graph.adj[u] {
        let nd = dist[u] + nb.1
        if nd < dist[nb.0] {
          dist[nb.0] = nd
        }
      }
    }
    out.push(dist)
  }
  out
}

///|
/// Deterministic LCG step (same parameters as `scanpy_draw_graph`).
fn scanpy_layout_rand(state : Int32) -> Int32 {
  (state * 1103515245 + 12345) & 0x7FFFFFFF
}

///|
/// Force-directed layout (Fruchterman-Reingold with ForceAtlas2-style
/// attraction scaling): O(n^2) pairwise repulsion `k^2 / d` and edge
/// attraction `w * d^2 / k`, temperature cooling over the iterations.
fn scanpy_draw_graph_fr(
  graph : SparseGraph,
  n_iter : Int32,
  seed : Int32,
  init_pos : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = graph.n
  if n == 0 {
    return []
  }
  let dim = 2
  let mut rng = if seed == 0 { 1 } else { seed }
  let pos : Array[Array[Double]] = Array(capacity=n)
  if init_pos.length() == n {
    for i = 0; i < n; i = i + 1 {
      let row : Array[Double] = Array(capacity=0)
      for d = 0; d < dim; d = d + 1 {
        row.push(if d < init_pos[i].length() { init_pos[i][d] } else { 0.0 })
      }
      pos.push(row)
    }
  } else {
    for _i = 0; _i < n; _i = _i + 1 {
      let row : Array[Double] = Array(capacity=0)
      for _d = 0; _d < dim; _d = _d + 1 {
        rng = scanpy_layout_rand(rng)
        row.push(rng.to_double() / 1073741823.5)
      }
      pos.push(row)
    }
  }
  let k = (1.0 / n.to_double()).sqrt()
  let k_sq = k * k
  for iter = 0; iter < n_iter; iter = iter + 1 {
    let temperature = (1.0 - iter.to_double() / n_iter.to_double()).max(0.01)
    let disp : Array[Array[Double]] = Array(capacity=n)
    for _i = 0; _i < n; _i = _i + 1 {
      disp.push(Array::make(dim, 0.0))
    }
    // repulsion (all pairs)
    for i = 0; i < n; i = i + 1 {
      for j = i + 1; j < n; j = j + 1 {
        let mut dist_sq = 0.0
        for d = 0; d < dim; d = d + 1 {
          let diff = pos[i][d] - pos[j][d]
          dist_sq = dist_sq + diff * diff
        }
        let dist = if dist_sq < 1.0e-4 { 1.0e-2 } else { dist_sq.sqrt() }
        let force = k_sq / dist
        for d = 0; d < dim; d = d + 1 {
          let delta = pos[i][d] - pos[j][d]
          disp[i][d] = disp[i][d] + delta / dist * force
          disp[j][d] = disp[j][d] - delta / dist * force
        }
      }
    }
    // attraction (edges only)
    for i = 0; i < n; i = i + 1 {
      for nb in graph.adj[i] {
        let j = nb.0
        if j <= i {
          continue
        }
        let mut dist_sq = 0.0
        for d = 0; d < dim; d = d + 1 {
          let diff = pos[j][d] - pos[i][d]
          dist_sq = dist_sq + diff * diff
        }
        let dist = if dist_sq < 1.0e-4 { 1.0e-2 } else { dist_sq.sqrt() }
        let force = dist / k * nb.1
        for d = 0; d < dim; d = d + 1 {
          let delta = pos[j][d] - pos[i][d]
          disp[i][d] = disp[i][d] + delta / dist * force
          disp[j][d] = disp[j][d] - delta / dist * force
        }
      }
    }
    // apply with cooling
    for i = 0; i < n; i = i + 1 {
      let mut mag = 0.0
      for d = 0; d < dim; d = d + 1 {
        mag = mag + disp[i][d] * disp[i][d]
      }
      mag = mag.sqrt()
      if mag > 1.0e-12 {
        let limited = if mag > temperature { temperature / mag } else { 1.0 }
        for d = 0; d < dim; d = d + 1 {
          pos[i][d] = pos[i][d] + disp[i][d] / mag * limited
        }
      }
    }
  }
  pos
}

///|
/// Kamada-Kawai stress layout: target distances are the weighted shortest
/// paths (unreachable pairs use twice the graph diameter), minimized by
/// gradient descent on `E = sum_{i<j} (d_ij - L_ij)^2 / L_ij`.
fn scanpy_draw_graph_kk(
  graph : SparseGraph,
  n_iter : Int32,
  seed : Int32,
  init_pos : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = graph.n
  if n == 0 {
    return []
  }
  let dim = 2
  let sp = scanpy_graph_shortest_paths(graph)
  let mut max_sp = 1.0
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      if sp[i][j] < 1.0e29 && sp[i][j] > max_sp {
        max_sp = sp[i][j]
      }
    }
  }
  let far = 2.0 * max_sp
  let mut rng = if seed == 0 { 2 } else { seed }
  let pos : Array[Array[Double]] = Array(capacity=n)
  if init_pos.length() == n {
    for i = 0; i < n; i = i + 1 {
      let row : Array[Double] = Array(capacity=0)
      for d = 0; d < dim; d = d + 1 {
        row.push(if d < init_pos[i].length() { init_pos[i][d] } else { 0.0 })
      }
      pos.push(row)
    }
  } else {
    for _i = 0; _i < n; _i = _i + 1 {
      let row : Array[Double] = Array(capacity=0)
      for _d = 0; _d < dim; _d = _d + 1 {
        rng = scanpy_layout_rand(rng)
        row.push(rng.to_double() / 536870911.5 - 1.0)
      }
      pos.push(row)
    }
  }
  for _iter = 0; _iter < n_iter; _iter = _iter + 1 {
    let grad : Array[Array[Double]] = Array(capacity=n)
    for _i = 0; _i < n; _i = _i + 1 {
      grad.push(Array::make(dim, 0.0))
    }
    for i = 0; i < n; i = i + 1 {
      for j = i + 1; j < n; j = j + 1 {
        let target = if sp[i][j] >= 1.0e29 { far } else { sp[i][j] }
        let mut dist_sq = 0.0
        for d = 0; d < dim; d = d + 1 {
          let diff = pos[i][d] - pos[j][d]
          dist_sq = dist_sq + diff * diff
        }
        let delta = if dist_sq < 1.0e-8 { 1.0e-4 } else { dist_sq.sqrt() }
        // dE/dpos_i = 2 (delta - target) (pos_i - pos_j) / (delta * target^2)
        let coef = 2.0 * (delta - target) / (delta * target * target)
        for d = 0; d < dim; d = d + 1 {
          let diff = pos[i][d] - pos[j][d]
          grad[i][d] = grad[i][d] + coef * diff
          grad[j][d] = grad[j][d] - coef * diff
        }
      }
    }
    // normalized step to keep the iteration numerically stable
    let mut max_grad = 1.0e-12
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < dim; d = d + 1 {
        let g = grad[i][d].abs()
        if g > max_grad {
          max_grad = g
        }
      }
    }
    let step = max_sp * 0.05 / max_grad
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < dim; d = d + 1 {
        pos[i][d] = pos[i][d] - step * grad[i][d]
      }
    }
  }
  pos
}

///|
/// SparseGraph-native graph layout dispatch (`tl.draw_graph(layout=...)`):
/// `"fr"` / `"fa2"` select the force-directed layout, `"kk"` the
/// Kamada-Kawai stress layout.  Returns an `n x 2` coordinate matrix.
pub fn scanpy_draw_graph_sparse(
  graph : SparseGraph,
  layout~ : String = "fr",
  n_iter~ : Int32 = 500,
  seed~ : Int32 = 0,
  init_pos~ : Array[Array[Double]] = [],
) -> Array[Array[Double]] {
  match layout {
    "kk" => scanpy_draw_graph_kk(graph, n_iter, seed, init_pos)
    _ => scanpy_draw_graph_fr(graph, n_iter, seed, init_pos)
  }
}
