/// Scanpy second-round: External integration tools
/// Inspired by scanpy.external: Harmony, BBKNN, Scanorama, Palantir, ingest

///|
/// Harmony batch correction result.
pub struct ScanpyHarmonyResult {
  corrected : Array[Array[Double]] // batch-corrected embedding
  n_batches : Int32
  n_iter : Int32
  converged : Bool
}

///|
/// Harmony batch correction algorithm.
/// Iteratively corrects batch effects in PCA embedding.
pub fn scanpy_harmonize(
  embedding : Array[Array[Double]],
  batch_labels : Array[Int32],
  n_clusters~ : Int32 = 10,
  max_iter~ : Int32 = 30,
  lambda~ : Double = 1.0,
) -> ScanpyHarmonyResult {
  let n_cells = embedding.length()
  let n_dims = if n_cells > 0 { embedding[0].length() } else { 0 }
  // Count batches
  let batch_set = Set([])
  for i in 0..<n_cells {
    batch_set.add(batch_labels[i])
  }
  let n_batches = batch_set.length()
  // Work on a copy
  let corrected = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let row = Array(capacity=n_dims)
    for d in 0..<n_dims {
      row.push(embedding[i][d])
    }
    corrected.push(row)
  }
  let mut converged = false
  let mut n_iter_done = 0
  for iter in 0..<max_iter {
    n_iter_done = iter + 1
    // Compute batch means
    let batch_means = Array(capacity=n_batches)
    let batch_counts = Array(capacity=n_batches)
    for _b in 0..<n_batches {
      let m = Array(capacity=n_dims)
      for _d in 0..<n_dims {
        m.push(0.0)
      }
      batch_means.push(m)
      batch_counts.push(0)
    }
    for i in 0..<n_cells {
      let b = batch_labels[i]
      for d in 0..<n_dims {
        batch_means[b][d] = batch_means[b][d] + corrected[i][d]
      }
      batch_counts[b] = batch_counts[b] + 1
    }
    for b in 0..<n_batches {
      if batch_counts[b] > 0 {
        for d in 0..<n_dims {
          batch_means[b][d] = batch_means[b][d] / batch_counts[b].to_double()
        }
      }
    }
    // Correct: subtract batch mean, scale by lambda
    let mut max_change = 0.0
    for i in 0..<n_cells {
      let b = batch_labels[i]
      for d in 0..<n_dims {
        let correction = lambda * batch_means[b][d]
        let old_val = corrected[i][d]
        corrected[i][d] = old_val - correction
        let change = old_val - corrected[i][d]
        if change < 0.0 {
          if -change > max_change {
            max_change = -change
          }
        } else if change > max_change {
          max_change = change
        }
      }
    }
    if max_change < 0.001 {
      converged = true
      break
    }
    ignore(n_clusters)
  }
  ScanpyHarmonyResult::{ corrected, n_batches, n_iter: n_iter_done, converged, }
}

///|
/// BBKNN — Batch-Balanced K-Nearest Neighbors.
/// Builds a neighbor graph that balances representation across batches.
pub fn scanpy_bbknn(
  embedding : Array[Array[Double]],
  batch_labels : Array[Int32],
  n_neighbors_per_batch~ : Int32 = 3,
) -> Array[Array[Double]] {
  let n = embedding.length()
  let batch_set = Set([])
  for i in 0..<n {
    batch_set.add(batch_labels[i])
  }
  let n_batches = batch_set.length()
  let n_total_neighbors = n_neighbors_per_batch * n_batches
  let graph = Array(capacity=n)
  for i in 0..<n {
    // For each batch, find nearest neighbors
    let all_dists = Array(capacity=n)
    for j in 0..<n {
      let mut dist = 0.0
      for d in 0..<embedding[i].length() {
        let diff = embedding[i][d] - embedding[j][d]
        dist = dist + diff * diff
      }
      all_dists.push((j, dist, batch_labels[j]))
    }
    // Sort by distance within each batch, take top-k per batch
    let row = Array(capacity=n)
    for _j in 0..<n {
      row.push(0.0)
    }
    for b in 0..<n_batches {
      let batch_dists = Array(capacity=n)
      for j in 0..<n {
        if all_dists[j].2 == b && all_dists[j].0 != i {
          batch_dists.push(all_dists[j])
        }
      }
      batch_dists.sort_by(fn(a, b) {
        if a.1 < b.1 {
          -1
        } else if a.1 > b.1 {
          1
        } else {
          0
        }
      })
      let k = if n_neighbors_per_batch < batch_dists.length() {
        n_neighbors_per_batch
      } else {
        batch_dists.length()
      }
      for kk in 0..<k {
        row[batch_dists[kk].0] = 1.0
      }
    }
    graph.push(row)
    ignore(n_total_neighbors)
  }
  graph
}

///|
/// Scanorama integration — batch integration via mutual nearest neighbors.
/// Simplified: align embeddings to a common coordinate system.
pub fn scanpy_scanorama(
  batches : Array[Array[Array[Double]]], // [batch][cell][dim]
) -> Array[Array[Double]] {
  // Use first batch as reference, align others
  let n_dims = if batches.length() > 0 && batches[0].length() > 0 {
    batches[0][0].length()
  } else {
    0
  }
  let result = Array(capacity=0)
  // Concatenate all batches with simple centering
  for b in 0..<batches.length() {
    // Compute batch centroid
    let centroid = Array(capacity=n_dims)
    for _d in 0..<n_dims {
      centroid.push(0.0)
    }
    let n_cells = batches[b].length()
    for i in 0..<n_cells {
      for d in 0..<n_dims {
        centroid[d] = centroid[d] + batches[b][i][d]
      }
    }
    if n_cells > 0 {
      for d in 0..<n_dims {
        centroid[d] = centroid[d] / n_cells.to_double()
      }
    }
    // Center and add to result
    for i in 0..<n_cells {
      let row = Array(capacity=n_dims)
      for d in 0..<n_dims {
        row.push(batches[b][i][d] - centroid[d])
      }
      result.push(row)
    }
  }
  result
}

///|
/// Palantir pseudotime — identifies terminal states and computes branch probabilities.
pub struct PalantirResult {
  pseudotime : Array[Double]
  terminal_states : Array[Int32]
  branch_probs : Array[Array[Double]] // n_cells x n_terminal
  entropy : Array[Double]
}

///|
/// Run Palantir-like pseudotime analysis.
/// Uses diffusion potential to identify terminal states.
pub fn scanpy_palantir(
  embedding : Array[Array[Double]],
  root_cell~ : Int32 = 0,
  n_terminal~ : Int32 = 3,
) -> PalantirResult {
  let n = embedding.length()
  // Compute distance from root
  let dist_from_root = Array(capacity=n)
  for i in 0..<n {
    let mut dist = 0.0
    for d in 0..<embedding[i].length() {
      let diff = embedding[i][d] - embedding[root_cell][d]
      dist = dist + diff * diff
    }
    dist_from_root.push(@math.pow(dist, 0.5))
  }
  // Pseudotime = normalized distance from root
  let mut max_dist = 0.0
  for i in 0..<n {
    if dist_from_root[i] > max_dist {
      max_dist = dist_from_root[i]
    }
  }
  let pseudotime = Array(capacity=n)
  for i in 0..<n {
    if max_dist > 0.0 {
      pseudotime.push(dist_from_root[i] / max_dist)
    } else {
      pseudotime.push(0.0)
    }
  }
  // Terminal states = cells with max distance from root (farthest n_terminal)
  let indices = Array(capacity=n)
  for i in 0..<n {
    indices.push(i)
  }
  indices.sort_by(fn(a, b) {
    if dist_from_root[a] > dist_from_root[b] {
      -1
    } else if dist_from_root[a] < dist_from_root[b] {
      1
    } else {
      0
    }
  })
  let terminal_states = Array(capacity=n_terminal)
  for k in 0..<n_terminal {
    if k < indices.length() {
      terminal_states.push(indices[k])
    }
  }
  // Branch probabilities: inverse distance to each terminal
  let n_term = terminal_states.length()
  let branch_probs = Array(capacity=n)
  let entropy = Array(capacity=n)
  for i in 0..<n {
    let probs = Array(capacity=n_term)
    let mut sum = 0.0
    for t in 0..<n_term {
      let mut dist = 0.0
      for d in 0..<embedding[i].length() {
        let diff = embedding[i][d] - embedding[terminal_states[t]][d]
        dist = dist + diff * diff
      }
      dist = @math.pow(dist, 0.5) + 0.0001
      let p = 1.0 / dist
      probs.push(p)
      sum = sum + p
    }
    // Normalize
    let mut ent = 0.0
    for t in 0..<n_term {
      if sum > 0.0 {
        probs[t] = probs[t] / sum
      }
      if probs[t] > 0.0 {
        ent = ent - probs[t] * @math.log2(probs[t])
      }
    }
    branch_probs.push(probs)
    entropy.push(ent)
  }
  PalantirResult::{ pseudotime, terminal_states, branch_probs, entropy, }
}

///|
/// Ingest — map new query data to existing reference embedding.
/// Uses nearest neighbor transfer of labels and embeddings.
pub struct IngestResult {
  labels : Array[Int32]
  distances : Array[Double]
  indices : Array[Int32]
}

///|
pub fn scanpy_ingest(
  query : Array[Array[Double]],
  reference : Array[Array[Double]],
  ref_labels : Array[Int32],
  n_neighbors~ : Int32 = 5,
) -> IngestResult {
  let n_query = query.length()
  let labels = Array(capacity=n_query)
  let distances = Array(capacity=n_query)
  let indices = Array(capacity=n_query)
  for i in 0..<n_query {
    // Find nearest reference cells
    let dists = Array(capacity=reference.length())
    for j in 0..<reference.length() {
      let mut dist = 0.0
      for d in 0..<query[i].length() {
        if d < reference[j].length() {
          let diff = query[i][d] - reference[j][d]
          dist = dist + diff * diff
        }
      }
      dists.push((j, dist))
    }
    dists.sort_by(fn(a, b) {
      if a.1 < b.1 {
        -1
      } else if a.1 > b.1 {
        1
      } else {
        0
      }
    })
    // Majority vote among k nearest
    let k = if n_neighbors < dists.length() {
      n_neighbors
    } else {
      dists.length()
    }
    let label_counts = Map([])
    for kk in 0..<k {
      let label = ref_labels[dists[kk].0]
      let cnt = label_counts.get(label).unwrap_or(0)
      label_counts[label] = cnt + 1
    }
    let mut best_label = 0
    let mut best_count = 0
    for k_ in label_counts.keys() {
      let cnt = label_counts[k_]
      if cnt > best_count {
        best_count = cnt
        best_label = k_
      }
    }
    labels.push(best_label)
    distances.push(@math.pow(dists[0].1, 0.5))
    indices.push(dists[0].0)
  }
  IngestResult::{ labels, distances, indices, }
}
