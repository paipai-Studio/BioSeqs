///|
/// Harmony - Fast iterative batch correction for single-cell data integration.
/// Inspired by the Harmony R package (Korsunsky et al., 2019).
/// Implements iterative soft k-means clustering with linear batch correction.

///|
/// Input data for Harmony batch correction.
pub struct HarmonyData {
  pub embeddings : Array[Array[Double]]
  pub batch_labels : Array[String]
  pub cell_ids : Array[String]
}

///|
/// Create HarmonyData from embeddings, batch labels, and cell IDs.
pub fn HarmonyData::new(
  embeddings : Array[Array[Double]],
  batch_labels : Array[String],
  cell_ids : Array[String],
) -> HarmonyData {
  HarmonyData::{ embeddings, batch_labels, cell_ids }
}

///|
/// Get embeddings matrix.
pub fn HarmonyData::get_embeddings(self : HarmonyData) -> Array[Array[Double]] {
  self.embeddings
}

///|
/// Get batch labels.
pub fn HarmonyData::get_batch_labels(self : HarmonyData) -> Array[String] {
  self.batch_labels
}

///|
/// Get cell IDs.
pub fn HarmonyData::get_cell_ids(self : HarmonyData) -> Array[String] {
  self.cell_ids
}

///|
/// Harmony hyperparameters.
pub struct HarmonyParams {
  n_clusters : Int
  max_iterations : Int
  tolerance : Double
  lambda : Double
  sigma : Double
  ridge_penalty : Double
  seed : Int
}

///|
/// Create default Harmony parameters.
pub fn HarmonyParams::new() -> HarmonyParams {
  HarmonyParams::{
    n_clusters: 20,
    max_iterations: 20,
    tolerance: 0.01,
    lambda: 1.0,
    sigma: 0.1,
    ridge_penalty: 1.0,
    seed: 42,
  }
}

///|
/// Create Harmony parameters with custom values.
pub fn HarmonyParams::create(
  n_clusters? : Int = 20,
  max_iterations? : Int = 20,
  tolerance? : Double = 0.01,
  lambda? : Double = 1.0,
  sigma? : Double = 0.1,
  ridge_penalty? : Double = 1.0,
  seed? : Int = 42,
) -> HarmonyParams {
  HarmonyParams::{
    n_clusters: n_clusters,
    max_iterations: max_iterations,
    tolerance: tolerance,
    lambda: lambda,
    sigma: sigma,
    ridge_penalty: ridge_penalty,
    seed: seed,
  }
}

///|
/// Result of Harmony batch correction.
pub struct HarmonyResult {
  pub corrected : Array[Array[Double]]
  pub membership : Array[Array[Double]]
  pub centroids : Array[Array[Double]]
  pub n_iterations : Int
  pub converged : Bool
  var_explained : Array[Double]
}

///|
/// Get corrected embeddings.
pub fn HarmonyResult::get_corrected(self : HarmonyResult) -> Array[Array[Double]] {
  self.corrected
}

///|
/// Get soft cluster membership matrix.
pub fn HarmonyResult::get_membership(self : HarmonyResult) -> Array[Array[Double]] {
  self.membership
}

///|
/// Get cluster centroids.
pub fn HarmonyResult::get_centroids(self : HarmonyResult) -> Array[Array[Double]] {
  self.centroids
}

///|
/// Get number of iterations.
pub fn HarmonyResult::get_n_iterations(self : HarmonyResult) -> Int {
  self.n_iterations
}

///|
/// Check if algorithm converged.
pub fn HarmonyResult::is_converged(self : HarmonyResult) -> Bool {
  self.converged
}

///|
/// Run Harmony batch correction on embedded data.
pub fn harmony_run(data : HarmonyData, params : HarmonyParams) -> HarmonyResult {
  let n_cells = data.embeddings.length()
  let n_dims = if n_cells > 0 { data.embeddings[0].length() } else { 0 }
  let n_clusters = if params.n_clusters > 0 { params.n_clusters } else { 1 }
  
  let centroids = harmony_init_centroids(data, n_clusters, params.seed)
  let membership = harmony_compute_membership(data.embeddings, centroids, params.lambda)
  let mut corrected = data.embeddings
  
  let mut iter = 0
  let mut converged = false
  let mut prev_obj = harmony_objective(data.embeddings, centroids, membership)
  
  while iter < params.max_iterations {
    let new_centroids = harmony_update_centroids(corrected, membership, n_clusters)
    let corrections = harmony_batch_corrections(
      HarmonyData::{ embeddings: corrected, batch_labels: data.batch_labels, cell_ids: data.cell_ids },
      membership,
      new_centroids,
      n_clusters,
      params.ridge_penalty,
    )
    corrected = harmony_apply_correction(corrected, membership, corrections, n_clusters)
    let new_membership = harmony_compute_membership(corrected, new_centroids, params.lambda)
    let obj = harmony_objective(corrected, new_centroids, new_membership)
    let delta = if prev_obj.abs() > 1.0e-10 { (prev_obj - obj).abs() / prev_obj.abs() } else { 0.0 }
    
    if delta < params.tolerance {
      converged = true
      break
    }
    prev_obj = obj
    iter = iter + 1
  }
  
  HarmonyResult::{
    corrected,
    membership,
    centroids,
    n_iterations: iter,
    converged,
    var_explained: Array::make(n_dims, 0.0),
  }
}

///|
/// Initialize cluster centroids using k-means++ style seeding.
fn harmony_init_centroids(data : HarmonyData, n_clusters : Int, seed : Int) -> Array[Array[Double]] {
  let n_cells = data.embeddings.length()
  if n_cells == 0 {
    return Array::make(n_clusters, Array::make(0, 0.0))
  }
  let n_dims = data.embeddings[0].length()
  
  let centroids : Array[Array[Double]] = Array::new()
  let mut rng_state = seed
  
  let first_idx = harmony_next_rand(rng_state) % n_cells
  rng_state = harmony_next_rand(rng_state)
  centroids.push(data.embeddings[first_idx])
  
  let mut c = 1
  while c < n_clusters {
    let dists : Array[Double] = Array::make(n_cells, 0.0)
    let mut i = 0
    while i < n_cells {
      let mut min_dist = 1.0e30
      let mut k = 0
      while k < c {
        let d = harmony_sq_dist(data.embeddings[i], centroids[k])
        if d < min_dist {
          min_dist = d
        }
        k = k + 1
      }
      dists[i] = min_dist
      i = i + 1
    }
    
    let total_dist = harmony_sum(dists)
    if total_dist > 0.0 {
      let r = harmony_next_rand(rng_state).to_double() / 2147483647.0 * total_dist
      rng_state = harmony_next_rand(rng_state)
      let mut best_idx = 0
      let mut cum = 0.0
      let mut j = 0
      while j < n_cells {
        cum = cum + dists[j]
        if cum >= r {
          best_idx = j
          break
        }
        j = j + 1
      }
      centroids.push(data.embeddings[best_idx])
    } else {
      let idx = harmony_next_rand(rng_state) % n_cells
      rng_state = harmony_next_rand(rng_state)
      centroids.push(data.embeddings[idx])
    }
    c = c + 1
  }
  centroids
}

///|
/// Compute soft membership (responsibilities) for each cell to each cluster.
fn harmony_compute_membership(
  embeddings : Array[Array[Double]],
  centroids : Array[Array[Double]],
  lambda : Double,
) -> Array[Array[Double]] {
  let n_cells = embeddings.length()
  let n_clusters = centroids.length()
  let membership : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n_cells {
    let responsibilities : Array[Double] = Array::make(n_clusters, 0.0)
    let mut j = 0
    let mut max_logit = -1.0e30
    while j < n_clusters {
      let d = harmony_sq_dist(embeddings[i], centroids[j])
      let logit = -lambda * d
      if logit > max_logit {
        max_logit = logit
      }
      responsibilities[j] = logit
      j = j + 1
    }
    
    let mut sum_exp = 0.0
    let mut k = 0
    while k < n_clusters {
      responsibilities[k] = @math.exp(responsibilities[k] - max_logit)
      sum_exp = sum_exp + responsibilities[k]
      k = k + 1
    }
    
    if sum_exp > 0.0 {
      let mut m = 0
      while m < n_clusters {
        responsibilities[m] = responsibilities[m] / sum_exp
        m = m + 1
      }
    }
    
    membership.push(responsibilities)
    i = i + 1
  }
  membership
}

///|
/// Update centroids as weighted means.
fn harmony_update_centroids(
  embeddings : Array[Array[Double]],
  membership : Array[Array[Double]],
  n_clusters : Int,
) -> Array[Array[Double]] {
  let n_cells = embeddings.length()
  let n_dims = if n_cells > 0 { embeddings[0].length() } else { 0 }
  let centroids : Array[Array[Double]] = Array::new()
  
  let mut k = 0
  while k < n_clusters {
    let mut weight_sum = 0.0
    let new_centroid : Array[Double] = Array::make(n_dims, 0.0)
    let mut i = 0
    while i < n_cells {
      let w = membership[i][k]
      weight_sum = weight_sum + w
      let mut d = 0
      while d < n_dims {
        new_centroid[d] = new_centroid[d] + w * embeddings[i][d]
        d = d + 1
      }
      i = i + 1
    }
    if weight_sum > 0.0 {
      let mut d = 0
      while d < n_dims {
        new_centroid[d] = new_centroid[d] / weight_sum
        d = d + 1
      }
    }
    centroids.push(new_centroid)
    k = k + 1
  }
  centroids
}

///|
/// Compute batch correction vectors per cluster and batch.
fn harmony_batch_corrections(
  data : HarmonyData,
  membership : Array[Array[Double]],
  centroids : Array[Array[Double]],
  n_clusters : Int,
  ridge_penalty : Double,
) -> Array[Array[Array[Double]]] {
  let n_cells = data.embeddings.length()
  let n_dims = if n_cells > 0 { data.embeddings[0].length() } else { 0 }
  
  // Get unique batches
  let unique_batches : Array[String] = Array::new()
  {
    let seen : Map[String, Bool] = Map([], capacity=16)
    let mut i = 0
    while i < n_cells {
      let batch = data.batch_labels[i]
      if !seen.get_or_default(batch, false) {
        seen.set(batch, true)
        unique_batches.push(batch)
      }
      i = i + 1
    }
  }
  
  let n_batches = unique_batches.length()
  let corrections : Array[Array[Array[Double]]] = Array::make(n_clusters, Array::make(n_batches, Array::make(n_dims, 0.0)))
  
  let mut k = 0
  while k < n_clusters {
    let mut b = 0
    while b < n_batches {
      let batch_id = unique_batches[b]
      // Compute batch mean for this cluster
      let mut weight_sum = 0.0
      let batch_mean : Array[Double] = Array::make(n_dims, 0.0)
      let global_mean : Array[Double] = Array::make(n_dims, 0.0)
      let mut gw = 0.0
      
      let mut i = 0
      while i < n_cells {
        let w = membership[i][k]
        if data.batch_labels[i] == batch_id {
          weight_sum = weight_sum + w
          let mut d = 0
          while d < n_dims {
            batch_mean[d] = batch_mean[d] + w * data.embeddings[i][d]
            d = d + 1
          }
        }
        gw = gw + w
        let mut d = 0
        while d < n_dims {
          global_mean[d] = global_mean[d] + w * data.embeddings[i][d]
          d = d + 1
        }
        i = i + 1
      }
      
      if weight_sum > 0.0 {
        let mut d = 0
        while d < n_dims {
          batch_mean[d] = batch_mean[d] / weight_sum
          d = d + 1
        }
      }
      if gw > 0.0 {
        let mut d = 0
        while d < n_dims {
          global_mean[d] = global_mean[d] / gw
          d = d + 1
        }
      }
      
      // Correction = batch_mean - global_mean (with ridge penalty)
      let correction : Array[Double] = Array::make(n_dims, 0.0)
      let effective_weight = weight_sum / (weight_sum + ridge_penalty)
      let mut d = 0
      while d < n_dims {
        correction[d] = effective_weight * (batch_mean[d] - global_mean[d])
        d = d + 1
      }
      corrections[k][b] = correction
      b = b + 1
    }
    k = k + 1
  }
  corrections
}

///|
/// Apply batch correction vectors to each cell.
fn harmony_apply_correction(
  embeddings : Array[Array[Double]],
  membership : Array[Array[Double]],
  corrections : Array[Array[Array[Double]]],
  n_clusters : Int,
) -> Array[Array[Double]] {
  let n_cells = embeddings.length()
  let n_dims = if n_cells > 0 { embeddings[0].length() } else { 0 }
  let n_batches = if n_clusters > 0 { corrections[0].length() } else { 0 }
  
  // Get unique batches from correction dimensions
  let corrected : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_cells {
    let new_cell : Array[Double] = Array::make(n_dims, 0.0)
    let mut d = 0
    while d < n_dims {
      new_cell[d] = embeddings[i][d]
      d = d + 1
    }
    
    let mut k = 0
    while k < n_clusters {
      let w = membership[i][k]
      let mut b = 0
      while b < n_batches {
        let mut d = 0
        while d < n_dims {
          new_cell[d] = new_cell[d] - w * corrections[k][b][d]
          d = d + 1
        }
        b = b + 1
      }
      k = k + 1
    }
    corrected.push(new_cell)
    i = i + 1
  }
  corrected
}

///|
/// Compute the Harmony objective function.
fn harmony_objective(
  embeddings : Array[Array[Double]],
  centroids : Array[Array[Double]],
  membership : Array[Array[Double]],
) -> Double {
  let n_cells = embeddings.length()
  let n_clusters = centroids.length()
  let mut obj = 0.0
  let mut i = 0
  while i < n_cells {
    let mut k = 0
    while k < n_clusters {
      let w = membership[i][k]
      let d = harmony_sq_dist(embeddings[i], centroids[k])
      obj = obj + w * d
      k = k + 1
    }
    i = i + 1
  }
  obj
}

///|
/// Squared Euclidean distance.
fn harmony_sq_dist(a : Array[Double], b : Array[Double]) -> Double {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    let diff = a[i] - b[i]
    sum = sum + diff * diff
    i = i + 1
  }
  sum
}

///|
/// Sum of array elements.
fn harmony_sum(arr : Array[Double]) -> Double {
  let mut s = 0.0
  for v in arr {
    s = s + v
  }
  s
}

///|
/// Simple linear congruential generator.
fn harmony_next_rand(seed : Int) -> Int {
  let state = seed * 1103515245 + 12345
  state & 0x7FFFFFFF
}

///|
/// Create example data for Harmony testing.
pub fn harmony_create_example(n_cells : Int, n_batches : Int, n_dims : Int) -> HarmonyData {
  let embeddings : Array[Array[Double]] = Array::new()
  let batch_labels : Array[String] = Array::new()
  let cell_ids : Array[String] = Array::new()
  let mut rng = 42
  
  let mut i = 0
  while i < n_cells {
    let batch_idx = i * n_batches / n_cells
    let batch_id = "batch_" + batch_idx.to_string()
    
    // Simulate batch effect
    let batch_effect = batch_idx.to_double() * 2.0
    
    let cell : Array[Double] = Array::make(n_dims, 0.0)
    let mut d = 0
    while d < n_dims {
      rng = harmony_next_rand(rng)
      let noise = (rng.to_double() / 2147483647.0 - 0.5) * 2.0
      cell[d] = batch_effect + noise
      d = d + 1
    }
    
    embeddings.push(cell)
    batch_labels.push(batch_id)
    cell_ids.push("cell_" + i.to_string())
    i = i + 1
  }
  
  HarmonyData::{ embeddings, batch_labels, cell_ids }
}
