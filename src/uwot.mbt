///|
/// uwot - Uniform Manifold Approximation and Projection.
/// Inspired by uwot Bioconductor/R package.
/// Implements the UMAP algorithm for dimensionality reduction.

///|
/// UMAP configuration parameters.
pub struct UmapConfig {
  n_neighbors : Int
  n_components : Int
  n_epochs : Int
  min_dist : Double
  spread : Double
  learning_rate : Double
  repulsion_strength : Double
  negative_sample_rate : Double
  random_seed : Int
  metric : String
}

///|
pub fn UmapConfig::new() -> UmapConfig {
  UmapConfig::{
    n_neighbors: 15,
    n_components: 2,
    n_epochs: 200,
    min_dist: 0.1,
    spread: 1.0,
    learning_rate: 1.0,
    repulsion_strength: 1.0,
    negative_sample_rate: 5.0,
    random_seed: 42,
    metric: "euclidean",
  }
}

///|
pub fn UmapConfig::new_custom(
  n_neighbors : Int,
  n_components : Int,
  n_epochs : Int,
  random_seed : Int,
) -> UmapConfig {
  UmapConfig::{
    n_neighbors: n_neighbors,
    n_components: n_components,
    n_epochs: n_epochs,
    min_dist: 0.1,
    spread: 1.0,
    learning_rate: 1.0,
    repulsion_strength: 1.0,
    negative_sample_rate: 5.0,
    random_seed: random_seed,
    metric: "euclidean",
  }
}

///|
/// UMAP result structure.
pub struct UmapResult {
  embedding : Array[Array[Double]]
  n_epochs : Int
}

///|
pub fn UmapResult::new() -> UmapResult {
  UmapResult::{
    embedding: Array::new(),
    n_epochs: 0,
  }
}

// ============================================================
// Distance and knn functions
// ============================================================

fn uwot_euclidean_distance(x : Array[Double], y : Array[Double]) -> Double {
  let mut sum_sq = 0.0
  let n = if x.length() < y.length() { x.length() } else { y.length() }
  let mut i = 0
  while i < n {
    let diff = x[i] - y[i]
    sum_sq = sum_sq + diff * diff
    i = i + 1
  }
  sum_sq.sqrt()
}

///|
/// Compute pairwise distance matrix.
pub fn uwot_distance_matrix(data : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = data.length()
  let dist_matrix : Array[Array[Double]] = Array::new()
  
  if n == 0 {
    return dist_matrix
  }
  
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      row.push(0.0)
      j = j + 1
    }
    dist_matrix.push(row)
    i = i + 1
  }
  
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let dist = uwot_euclidean_distance(data[i], data[j])
      dist_matrix[i][j] = dist
      dist_matrix[j][i] = dist
      j = j + 1
    }
    i = i + 1
  }
  
  dist_matrix
}

///|
/// Find k-nearest neighbors for each point.
/// Returns: for each point, array of (neighbor_index, distance) pairs
fn find_knn(dist_matrix : Array[Array[Double]], k : Int) -> Array[Array[(Int, Double)]] {
  let n = dist_matrix.length()
  let result : Array[Array[(Int, Double)]] = Array::new()
  
  if n == 0 || k <= 0 {
    return result
  }
  
  let actual_k = if k < n { k } else { n - 1 }
  
  let mut i = 0
  while i < n {
    // Create list of (index, distance) for all other points
    let neighbors : Array[(Int, Double)] = Array::new()
    let mut j = 0
    while j < n {
      if i != j {
        neighbors.push((j, dist_matrix[i][j]))
      }
      j = j + 1
    }
    
    // Simple selection sort to find k nearest
    let knn : Array[(Int, Double)] = Array::new()
    let mut selected = 0
    while selected < actual_k && selected < neighbors.length() {
      let mut min_dist = 1.0e300
      let mut min_idx = -1
      let mut j = 0
      while j < neighbors.length() {
        let mut already_selected = false
        let mut s = 0
        while s < knn.length() {
          if knn[s].0 == neighbors[j].0 {
            already_selected = true
            break
          }
          s = s + 1
        }
        
        if !already_selected && neighbors[j].1 < min_dist {
          min_dist = neighbors[j].1
          min_idx = j
        }
        j = j + 1
      }
      
      if min_idx >= 0 {
        knn.push(neighbors[min_idx])
      }
      selected = selected + 1
    }
    
    result.push(knn)
    i = i + 1
  }
  
  result
}

// ============================================================
// Fuzzy simplicial set construction
// ============================================================

///|
/// Compute the membership strength for a distance using the smooth k-nn distance.
fn compute_membership(distance : Double, rho : Double, sigma : Double) -> Double {
  if distance <= rho {
    return 1.0
  }
  @math.exp(-(distance - rho) / sigma)
}

///|
/// Find sigma (normalization factor) for each point using binary search.
fn find_sigma(
  distances : Array[Double],
  rho : Double,
  target_log2 : Double,
  n_neighbors : Int,
) -> Double {
  let mut lo = 0.0
  let mut hi = 1000.0
  
  let mut mid = 1.0
  let mut iter = 0
  while iter < 64 {
    mid = (lo + hi) / 2.0
    let mut sum_membership = 0.0
    let mut i = 0
    while i < n_neighbors && i < distances.length() {
      sum_membership = sum_membership + compute_membership(distances[i], rho, mid)
      i = i + 1
    }
    
    let log_sum = @math.log2(sum_membership)
    
    if log_sum > target_log2 {
      hi = mid
    } else {
      lo = mid
    }
    
    iter = iter + 1
  }
  
  mid
}

///|
/// Build the fuzzy simplicial set (high-dimensional graph).
fn build_graph(
  data : Array[Array[Double]],
  dist_matrix : Array[Array[Double]],
  knn_indices : Array[Array[Int]],
  knn_distances : Array[Array[Double]],
  config : UmapConfig,
) -> Array[Array[Double]] {
  let n = data.length()
  let n_neighbors = config.n_neighbors
  
  // Initialize graph
  let graph : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      row.push(0.0)
      j = j + 1
    }
    graph.push(row)
    i = i + 1
  }
  
  // Target entropy = log2(n_neighbors)
  let target_log2 = @math.log2(n_neighbors.to_double())
  
  // For each point, compute rho (distance to nearest neighbor) and sigma
  let sigmas : Array[Double] = Array::new()
  let rhos : Array[Double] = Array::new()
  
  let mut i = 0
  while i < n {
    // rho = distance to nearest neighbor
    let rho = if knn_distances[i].length() > 0 { knn_distances[i][0] } else { 0.0 }
    rhos.push(rho)
    
    // Find sigma using binary search
    let distances : Array[Double] = Array::new()
    let mut j = 0
    while j < knn_distances[i].length() {
      distances.push(knn_distances[i][j])
      j = j + 1
    }
    
    let sigma = find_sigma(distances, rho, target_log2, n_neighbors)
    sigmas.push(sigma)
    
    i = i + 1
  }
  
  // Build the graph with asymmetric memberships
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < knn_indices[i].length() {
      let neighbor_idx = knn_indices[i][j]
      let dist = knn_distances[i][j]
      
      let membership_i_to_j = compute_membership(dist, rhos[i], sigmas[i])
      graph[i][neighbor_idx] = membership_i_to_j
      
      j = j + 1
    }
    i = i + 1
  }
  
  // Symmetrize the graph: w_ij = w_ij + w_ji - w_ij * w_ji
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let w_ij = graph[i][j]
      let w_ji = graph[j][i]
      let combined = w_ij + w_ji - w_ij * w_ji
      graph[i][j] = combined
      graph[j][i] = combined
      j = j + 1
    }
    i = i + 1
  }
  
  graph
}

// ============================================================
// Low-dimensional embedding
// ============================================================

///|
/// Compute the low-dimensional probability (using curve from min_dist and spread).
fn low_dim_probability(q_dist : Double, a : Double, b : Double) -> Double {
  1.0 / (1.0 + a * @math.pow(q_dist, 2.0 * b))
}

///|
/// Find a and b parameters based on min_dist and spread.
fn find_ab(min_dist : Double, spread : Double) -> (Double, Double) {
  // Curve: 1 / (1 + a * x^(2b))
  // At x = min_dist: y = 0.5 (approximately)
  // At x = spread: y is small
  
  // Simplified parameter estimation
  let a = 1.576942937
  let b = 0.895061728
  
  // Scale based on min_dist and spread
  let scaled_a = a / @math.pow(spread, 2.0 * b)
  (scaled_a, b)
}

///|
/// Initialize the embedding with small random values.
fn initialize_embedding(n : Int, n_components : Int, seed : Int) -> Array[Array[Double]] {
  let embedding : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n_components {
      let seed_val = (i * n_components + j + seed).to_double()
      let rand_val = (umap_pseudo_random(seed_val) - 0.5) * 0.01
      row.push(rand_val)
      j = j + 1
    }
    embedding.push(row)
    i = i + 1
  }
  
  embedding
}

///|
/// Main UMAP algorithm.
pub fn uwot_umap(data : Array[Array[Double]], config : UmapConfig) -> UmapResult {
  let n = data.length()
  
  if n == 0 {
    return UmapResult::new()
  }
  
  if n < 2 {
    let result = UmapResult::new()
    let embedding : Array[Array[Double]] = Array::new()
    let row : Array[Double] = Array::new()
    let mut i = 0
    while i < config.n_components {
      row.push(0.0)
      i = i + 1
    }
    embedding.push(row)
    return result
  }
  
  // Step 1: Compute distance matrix
  let dist_matrix = uwot_distance_matrix(data)
  
  // Step 2: Find k-nearest neighbors
  let k = if config.n_neighbors < n - 1 { config.n_neighbors } else { n - 1 }
  let knn = find_knn(dist_matrix, k)
  
  // Extract knn indices and distances
  let knn_indices : Array[Array[Int]] = Array::new()
  let knn_distances : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n {
    let indices : Array[Int] = Array::new()
    let distances : Array[Double] = Array::new()
    let mut j = 0
    while j < knn[i].length() {
      indices.push(knn[i][j].0)
      distances.push(knn[i][j].1)
      j = j + 1
    }
    knn_indices.push(indices)
    knn_distances.push(distances)
    i = i + 1
  }
  
  // Step 3: Build fuzzy simplicial set
  let graph = build_graph(data, dist_matrix, knn_indices, knn_distances, config)
  
  // Step 4: Find a and b parameters
  let (a, b) = find_ab(config.min_dist, config.spread)
  
  // Step 5: Initialize embedding
  let embedding = initialize_embedding(n, config.n_components, config.random_seed)
  
  // Step 6: Optimize embedding using stochastic gradient descent
  let epochs_per_sample : Array[Array[Int]] = Array::new()
  let epoch_of_next_sample : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n {
    let eps_row : Array[Int] = Array::new()
    let eons_row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      eps_row.push(0)
      eons_row.push(0.0)
      j = j + 1
    }
    epochs_per_sample.push(eps_row)
    epoch_of_next_sample.push(eons_row)
    i = i + 1
  }
  
  // Compute epochs per sample based on graph weights
  let max_epochs = config.n_epochs.to_double()
  let mut max_weight = 0.0
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      if graph[i][j] > max_weight {
        max_weight = graph[i][j]
      }
      j = j + 1
    }
    i = i + 1
  }
  
  if max_weight == 0.0 {
    max_weight = 1.0
  }
  
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      if i != j && graph[i][j] > 0.0 {
        let epochs = max_epochs * graph[i][j] / max_weight
        epochs_per_sample[i][j] = epochs.to_int()
        if epochs_per_sample[i][j] < 1 {
          epochs_per_sample[i][j] = 1
        }
        epoch_of_next_sample[i][j] = epochs
      } else {
        epochs_per_sample[i][j] = -1
      }
      j = j + 1
    }
    i = i + 1
  }
  
  // SGD optimization
  let mut epoch = 0
  while epoch < config.n_epochs {
    let mut i = 0
    while i < n {
      let mut j = 0
      while j < n {
        if i != j && epochs_per_sample[i][j] > 0 {
          if epoch.to_double() >= epoch_of_next_sample[i][j] {
            // Positive sample: move points together
            let mut dist_sq = 0.0
            let mut d = 0
            while d < config.n_components {
              let diff = embedding[i][d] - embedding[j][d]
              dist_sq = dist_sq + diff * diff
              d = d + 1
            }
            
            let grad_coeff = (-2.0 * a * b * @math.pow(dist_sq, b - 1.0)) / (1.0 + a * @math.pow(dist_sq, b))
            
            let mut d = 0
            while d < config.n_components {
              let grad = config.learning_rate * grad_coeff * (embedding[i][d] - embedding[j][d])
              embedding[i][d] = embedding[i][d] + grad
              embedding[j][d] = embedding[j][d] - grad
              d = d + 1
            }
            
            // Negative samples: move points apart
            let n_neg = config.negative_sample_rate.to_int()
            let mut neg = 0
            while neg < n_neg {
              // Pick random negative sample
              let rand_val = umap_pseudo_random((epoch * n + i * n_neg + neg + config.random_seed).to_double())
              let k_random = (rand_val * n.to_double()).to_int()
              
              if k_random != i && k_random != j {
                let mut dist_sq_neg = 0.0
                let mut d = 0
                while d < config.n_components {
                  let diff = embedding[i][d] - embedding[k_random][d]
                  dist_sq_neg = dist_sq_neg + diff * diff
                  d = d + 1
                }
                
                let grad_coeff_neg = (2.0 * b * config.repulsion_strength) / ((0.001 + dist_sq_neg) * (1.0 + a * @math.pow(dist_sq_neg, b)))
                
                let mut d = 0
                while d < config.n_components {
                  let grad = config.learning_rate * grad_coeff_neg * (embedding[i][d] - embedding[k_random][d])
                  embedding[i][d] = embedding[i][d] + grad
                  d = d + 1
                }
              }
              
              neg = neg + 1
            }
            
            epoch_of_next_sample[i][j] = epoch_of_next_sample[i][j] + epochs_per_sample[i][j].to_double()
          }
        }
        j = j + 1
      }
      i = i + 1
    }
    
    epoch = epoch + 1
  }
  
  UmapResult::{
    embedding: embedding,
    n_epochs: config.n_epochs,
  }
}

// ============================================================
// Helper functions
// ============================================================

fn umap_pseudo_random(seed : Double) -> Double {
  let x = seed * 12.9898
  let y_val = seed * 78.233
  let z = x * y_val
  let sin_val = @math.sin(z)
  let result = sin_val * 43758.5453
  result - result.floor()
}

///|
/// Create sample data for testing/demo.
pub fn create_umap_test_data(n_samples : Int, n_features : Int) -> Array[Array[Double]] {
  let data : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n_samples {
    let row : Array[Double] = Array::new()
    let group = i % 4  // 4 clusters
    
    let mut j = 0
    while j < n_features {
      let seed_val = (i * n_features + j).to_double()
      let noise = (umap_pseudo_random(seed_val) - 0.5) * 2.0
      let cluster_center = group.to_double() * 4.0
      row.push(cluster_center + noise)
      j = j + 1
    }
    
    data.push(row)
    i = i + 1
  }
  
  data
}
