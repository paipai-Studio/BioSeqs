/// Scanpy-inspired: pp.scrublet — doublet detection
/// Simulates doublets and classifies cells as singlets/doublets based on
/// nearest-neighbor overlap in gene expression space.
///
/// Reference: Wolock et al. (2019) Genome Biology

///|
/// Result of Scrublet doublet detection.
pub struct ScrubletResult {
  doublet_scores : Array[Double] // per-cell doublet probability
  predicted_doublets : Array[Bool] // true if predicted doublet
  threshold : Double // score threshold for calling doublets
  n_doublets : Int32 // number of predicted doublets
  mean_scores : Double // mean doublet score
}

///|
/// Run Scrublet-like doublet detection.
///
/// Parameters:
/// - `x`: cells × genes count matrix
/// - `sim_doublet_ratio`: ratio of simulated doublets to real cells (default 2.0)
/// - `n_neighbors`: number of nearest neighbors (default 10)
/// - `threshold`: doublet score threshold (default 0.25; if 0, auto-compute)
pub fn scanpy_scrublet(
  x : Array[Array[Double]],
  sim_doublet_ratio~ : Double = 2.0,
  n_neighbors~ : Int32 = 10,
  threshold~ : Double = 0.0,
) -> ScrubletResult {
  let n_cells = x.length()
  if n_cells == 0 {
    return ScrubletResult::{
      doublet_scores: [],
      predicted_doublets: [],
      threshold: 0.0,
      n_doublets: 0,
      mean_scores: 0.0,
    }
  }
  let n_genes = x[0].length()
  // Normalize each cell to unit total counts
  let normalized = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let total = scanpy_scrublet_row_sum(x[i])
    let row = Array(capacity=n_genes)
    for j in 0..<n_genes {
      if total > 0.0 {
        row.push(x[i][j] / total)
      } else {
        row.push(0.0)
      }
    }
    normalized.push(row)
  }
  // Simulate doublets by averaging random pairs
  let n_sim = (n_cells.to_double() * sim_doublet_ratio).floor().to_int()
  let simulated = Array(capacity=n_sim)
  let mut seed : Int32 = 42
  for _ in 0..<n_sim {
    // Simple LCG for reproducible random pairs
    seed = (seed * 1103515245 + 12345) % 2147483647
    let i1 = seed.reinterpret_as_uint() % n_cells.reinterpret_as_uint()
    seed = (seed * 1103515245 + 12345) % 2147483647
    let i2 = seed.reinterpret_as_uint() % n_cells.reinterpret_as_uint()
    let i1s = i1
    let i2s = i2
    let dbl = Array(capacity=n_genes)
    for j in 0..<n_genes {
      dbl.push(
        (
          normalized[i1s.reinterpret_as_int()][j] +
          normalized[i2s.reinterpret_as_int()][j]
        ) /
        2.0,
      )
    }
    simulated.push(dbl)
  }
  // Compute doublet score: fraction of k-nearest neighbors that are simulated
  let k = if n_neighbors < n_cells { n_neighbors } else { n_cells }
  let scores = Array(capacity=n_cells)
  for i in 0..<n_cells {
    // Compute distances to all simulated doublets
    let dists = Array(capacity=n_sim)
    for s in 0..<n_sim {
      let d = scanpy_scrublet_euclidean(normalized[i], simulated[s])
      dists.push((d, true)) // true = simulated
    }
    // Also include distances to real cells
    for j in 0..<n_cells {
      if j != i {
        let d = scanpy_scrublet_euclidean(normalized[i], normalized[j])
        dists.push((d, false)) // false = real
      }
    }
    // Sort by distance and count simulated in top-k
    dists.sort_by(fn(a, b) {
      if a.0 < b.0 {
        -1
      } else if a.0 > b.0 {
        1
      } else {
        0
      }
    })
    let mut sim_count = 0
    for nn in 0..<k {
      if nn < dists.length() && dists[nn].1 {
        sim_count = sim_count + 1
      }
    }
    scores.push(sim_count.to_double() / k.to_double())
  }
  // Determine threshold
  let actual_threshold = if threshold > 0.0 {
    threshold
  } else {
    // Auto: use mean + 2*std
    let mean = scanpy_scrublet_mean(scores)
    let variance = scanpy_scrublet_variance(scores, mean)
    mean + 2.0 * variance.sqrt()
  }
  // Predict doublets
  let predicted = Array(capacity=n_cells)
  let mut n_dbl = 0
  for i in 0..<n_cells {
    let is_dbl = scores[i] >= actual_threshold
    predicted.push(is_dbl)
    if is_dbl {
      n_dbl = n_dbl + 1
    }
  }
  ScrubletResult::{
    doublet_scores: scores,
    predicted_doublets: predicted,
    threshold: actual_threshold,
    n_doublets: n_dbl,
    mean_scores: scanpy_scrublet_mean(scores),
  }
}

///|
fn scanpy_scrublet_row_sum(row : Array[Double]) -> Double {
  let mut s = 0.0
  for v in row {
    s += v
  }
  s
}

///|
fn scanpy_scrublet_euclidean(a : Array[Double], b : Array[Double]) -> Double {
  let mut s = 0.0
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  for i in 0..<n {
    let diff = a[i] - b[i]
    s += diff * diff
  }
  s.sqrt()
}

///|
fn scanpy_scrublet_mean(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let mut s = 0.0
  for v in arr {
    s += v
  }
  s / arr.length().to_double()
}

///|
fn scanpy_scrublet_variance(arr : Array[Double], mean : Double) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let mut s = 0.0
  for v in arr {
    let diff = v - mean
    s += diff * diff
  }
  s / arr.length().to_double()
}
