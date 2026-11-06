/// Scanpy second-round: RNA velocity analysis
/// Inspired by scanpy.tl.velocity, velocity_graph, velocity_confidence, etc.

///|
/// RNA velocity result: estimated velocity for each cell-gene pair.
pub struct VelocityResult {
  velocity : Array[Array[Double]] // n_cells x n_genes velocity matrix
  confidence : Array[Double] // per-cell velocity confidence
  pseudotime : Array[Double] // velocity pseudotime per cell
}

///|
/// Compute RNA velocity from spliced and unspliced counts.
/// velocity = unspliced / (spliced + unspliced) * gamma - spliced * alpha
/// Simplified: velocity = unspliced - gamma * spliced
pub fn scanpy_velocity(
  spliced : Array[Array[Double]],
  unspliced : Array[Array[Double]],
  gamma~ : Double = 1.0,
) -> Array[Array[Double]] {
  let n_cells = spliced.length()
  let n_genes = if n_cells > 0 { spliced[0].length() } else { 0 }
  let velocity = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let row = Array(capacity=n_genes)
    for j in 0..<n_genes {
      let s = spliced[i][j]
      let u = if i < unspliced.length() && j < unspliced[i].length() {
        unspliced[i][j]
      } else {
        0.0
      }
      // velocity = unspliced - gamma * spliced
      row.push(u - gamma * s)
    }
    velocity.push(row)
  }
  velocity
}

///|
/// Build velocity graph: transition probabilities between cells.
/// For each cell, compute cosine similarity of velocity vector to
/// displacement vectors of neighboring cells.
pub fn scanpy_velocity_graph(
  velocity : Array[Array[Double]],
  neighbors : Array[Array[Double]],
  n_neighbors~ : Int32 = 30,
) -> Array[Array[Double]] {
  let n = velocity.length()
  let graph = Array(capacity=n)
  for i in 0..<n {
    let row = Array(capacity=n)
    for j in 0..<n {
      if neighbors[i][j] > 0.0 && i != j {
        // Cosine similarity between velocity[i] and displacement (j - i)
        let mut dot = 0.0
        let mut norm_v = 0.0
        let mut norm_d = 0.0
        for g in 0..<velocity[i].length() {
          let v = velocity[i][g]
          let d = if g < velocity[j].length() {
            velocity[j][g] - velocity[i][g]
          } else {
            0.0
          }
          dot = dot + v * d
          norm_v = norm_v + v * v
          norm_d = norm_d + d * d
        }
        if norm_v > 0.0 && norm_d > 0.0 {
          row.push(dot / (@math.pow(norm_v, 0.5) * @math.pow(norm_d, 0.5)))
        } else {
          row.push(0.0)
        }
      } else {
        row.push(0.0)
      }
    }
    // Keep only top n_neighbors
    let indices = Array(capacity=n)
    for j in 0..<n {
      indices.push(j)
    }
    indices.sort_by(fn(a, b) {
      if row[a] > row[b] {
        -1
      } else if row[a] < row[b] {
        1
      } else {
        0
      }
    })
    let k = if n_neighbors < n { n_neighbors } else { n }
    let keep = Set([])
    for kk in 0..<k {
      if kk < indices.length() {
        keep.add(indices[kk])
      }
    }
    for j in 0..<n {
      if not(keep.contains(j)) {
        row[j] = 0.0
      }
    }
    graph.push(row)
  }
  graph
}

///|
/// Compute velocity confidence per cell.
/// Confidence = mean of positive transition probabilities.
pub fn scanpy_velocity_confidence(
  velocity_graph : Array[Array[Double]],
) -> Array[Double] {
  let n = velocity_graph.length()
  let result = Array(capacity=n)
  for i in 0..<n {
    let mut sum = 0.0
    let mut count = 0
    for j in 0..<n {
      if velocity_graph[i][j] > 0.0 {
        sum = sum + velocity_graph[i][j]
        count = count + 1
      }
    }
    if count > 0 {
      result.push(sum / count.to_double())
    } else {
      result.push(0.0)
    }
  }
  result
}

///|
/// Compute velocity pseudotime via random walk on velocity graph.
/// Starting from root cell, propagate through graph.
pub fn scanpy_velocity_pseudotime(
  velocity_graph : Array[Array[Double]],
  root_cell~ : Int32 = 0,
) -> Array[Double] {
  let n = velocity_graph.length()
  let pseudotime = Array(capacity=n)
  let visited = Array(capacity=n)
  for _i in 0..<n {
    pseudotime.push(-1.0)
    visited.push(false)
  }
  // BFS from root
  let queue = Array(capacity=n)
  queue.push(root_cell)
  pseudotime[root_cell] = 0.0
  visited[root_cell] = true
  let mut head = 0
  while head < queue.length() {
    let curr = queue[head]
    head = head + 1
    for j in 0..<n {
      if velocity_graph[curr][j] > 0.0 && not(visited[j]) {
        visited[j] = true
        pseudotime[j] = pseudotime[curr] + 1.0
        queue.push(j)
      }
    }
  }
  // Normalize to [0, 1]
  let mut max_t = 0.0
  for i in 0..<n {
    if pseudotime[i] > max_t {
      max_t = pseudotime[i]
    }
  }
  if max_t > 0.0 {
    for i in 0..<n {
      pseudotime[i] = pseudotime[i] / max_t
    }
  }
  pseudotime
}

///|
/// Cluster cells based on velocity direction.
/// Uses k-means on velocity vectors.
pub fn scanpy_velocity_clusters(
  velocity : Array[Array[Double]],
  n_clusters~ : Int32 = 5,
  max_iter~ : Int32 = 100,
) -> Array[Int32] {
  let n = velocity.length()
  let n_genes = if n > 0 { velocity[0].length() } else { 0 }
  if n == 0 || n_genes == 0 {
    return []
  }
  // Initialize centroids (pick evenly spaced cells)
  let centroids = Array(capacity=n_clusters)
  for k in 0..<n_clusters {
    let idx = k * n / n_clusters
    let c = Array(capacity=n_genes)
    for g in 0..<n_genes {
      c.push(velocity[idx][g])
    }
    centroids.push(c)
  }
  let labels = Array(capacity=n)
  for _i in 0..<n {
    labels.push(0)
  }
  for _iter in 0..<max_iter {
    let mut changed = false
    for i in 0..<n {
      let mut best_k = 0
      let mut best_dist = 999999999999999.0
      for k in 0..<n_clusters {
        let mut dist = 0.0
        for g in 0..<n_genes {
          let diff = velocity[i][g] - centroids[k][g]
          dist = dist + diff * diff
        }
        if dist < best_dist {
          best_dist = dist
          best_k = k
        }
      }
      if labels[i] != best_k {
        changed = true
        labels[i] = best_k
      }
    }
    // Update centroids
    for k in 0..<n_clusters {
      for g in 0..<n_genes {
        let mut sum = 0.0
        let mut count = 0
        for i in 0..<n {
          if labels[i] == k {
            sum = sum + velocity[i][g]
            count = count + 1
          }
        }
        if count > 0 {
          centroids[k][g] = sum / count.to_double()
        }
      }
    }
    if not(changed) {
      break
    }
  }
  labels
}

///|
/// Recover dynamics parameters from velocity.
/// Estimates alpha (splicing rate), beta (degradation rate), gamma.
pub fn scanpy_recover_dynamics(
  spliced : Array[Array[Double]],
  unspliced : Array[Array[Double]],
  velocity : Array[Array[Double]],
) -> Array[(Double, Double, Double)] {
  let n_genes = if spliced.length() > 0 { spliced[0].length() } else { 0 }
  let params = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let mut sum_s = 0.0
    let mut sum_u = 0.0
    let mut sum_v = 0.0
    let mut count = 0
    for i in 0..<spliced.length() {
      sum_s = sum_s + spliced[i][g]
      if i < unspliced.length() {
        sum_u = sum_u + unspliced[i][g]
      }
      if i < velocity.length() {
        sum_v = sum_v + velocity[i][g]
      }
      count = count + 1
    }
    if count > 0 {
      let mean_s = sum_s / count.to_double()
      let mean_u = sum_u / count.to_double()
      let mean_v = sum_v / count.to_double()
      // alpha: transcription rate (estimated from unspliced production)
      let alpha = if mean_u > 0.0 { mean_u } else { 0.0 }
      // beta: splicing rate (estimated from s/u ratio)
      let beta = if mean_u > 0.0 && mean_s > 0.0 {
        mean_s / mean_u
      } else {
        1.0
      }
      // gamma: degradation rate (estimated from velocity/s ratio)
      let gamma = if mean_s > 0.0 { (mean_u - mean_v) / mean_s } else { 1.0 }
      params.push((alpha, beta, gamma))
    } else {
      params.push((0.0, 1.0, 1.0))
    }
  }
  params
}

///|
pub fn _velocity_result_dummy() -> Unit {
  let _ = VelocityResult::{ velocity: [], confidence: [], pseudotime: [], }
}
