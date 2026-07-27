///|
/// Rtsne - t-Distributed Stochastic Neighbor Embedding.
/// Inspired by Rtsne Bioconductor package.
/// Implements the t-SNE algorithm for dimensionality reduction.

///|
/// t-SNE configuration parameters.
pub struct TsneConfig {
  perplexity : Double
  theta : Double
  max_iter : Int
  dims : Int
  eta : Double
  exaggeration_factor : Double
  momentum : Double
  final_momentum : Double
  mom_switch_iter : Int
  stop_lying_iter : Int
  random_seed : Int
}

///|
pub fn TsneConfig::new() -> TsneConfig {
  TsneConfig::{
    perplexity: 30.0,
    theta: 0.5,
    max_iter: 1000,
    dims: 2,
    eta: 200.0,
    exaggeration_factor: 12.0,
    momentum: 0.5,
    final_momentum: 0.8,
    mom_switch_iter: 250,
    stop_lying_iter: 250,
    random_seed: 42,
  }
}

///|
pub fn TsneConfig::new_custom(
  perplexity : Double,
  max_iter : Int,
  dims : Int,
  random_seed : Int,
) -> TsneConfig {
  TsneConfig::{
    perplexity: perplexity,
    theta: 0.5,
    max_iter: max_iter,
    dims: dims,
    eta: 200.0,
    exaggeration_factor: 12.0,
    momentum: 0.5,
    final_momentum: 0.8,
    mom_switch_iter: 250,
    stop_lying_iter: 250,
    random_seed: random_seed,
  }
}

///|
/// t-SNE result structure.
pub struct TsneResult {
  embedding : Array[Array[Double]]
  costs : Array[Double]
  n_iter : Int
}

///|
pub fn TsneResult::new() -> TsneResult {
  TsneResult::{
    embedding: Array::new(),
    costs: Array::new(),
    n_iter: 0,
  }
}

// ============================================================
// Distance calculation
// ============================================================

///|
/// Calculate pairwise Euclidean distances between data points.
/// Input: data matrix (rows = samples, cols = features)
/// Output: distance matrix (n x n)
pub fn calc_distance_matrix(data : Array[Array[Double]]) -> Array[Array[Double]] {
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
      let dist = rtsne_euclidean_distance(data[i], data[j])
      dist_matrix[i][j] = dist
      dist_matrix[j][i] = dist
      j = j + 1
    }
    i = i + 1
  }
  
  dist_matrix
}

fn rtsne_euclidean_distance(x : Array[Double], y : Array[Double]) -> Double {
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

// ============================================================
// Perplexity and conditional probabilities
// ============================================================

///|
/// Compute the probability distribution p_{j|i} using binary search for perplexity.
fn compute_p_ji(distances : Array[Double], perplexity : Double, n : Int) -> Array[Double] {
  let p : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    p.push(0.0)
    i = i + 1
  }
  
  let log_perp = @math.ln(perplexity)
  
  // Binary search for sigma (precision)
  let mut beta = 1.0  // beta = 1 / (2 * sigma^2)
  let mut min_beta = -1.0e300
  let mut max_beta = 1.0e300
  
  let mut iter_count = 0
  while iter_count < 50 {
    // Compute P values
    let mut sum_p = 0.0
    let mut i = 0
    while i < n {
      let p_val = @math.exp(-distances[i] * distances[i] * beta)
      p[i] = p_val
      sum_p = sum_p + p_val
      i = i + 1
    }
    
    if sum_p == 0.0 {
      sum_p = 1.0e-10
    }
    
    // Normalize
    let mut h = 0.0
    i = 0
    while i < n {
      p[i] = p[i] / sum_p
      if p[i] > 0.0 {
        h = h - p[i] * @math.ln(p[i])
      }
      i = i + 1
    }
    
    let h_diff = h - log_perp
    
    if h_diff.abs() < 1.0e-5 {
      break
    }
    
    if h_diff > 0.0 {
      // Entropy is too high, need to increase beta (decrease sigma)
      min_beta = beta
      if max_beta == 1.0e300 {
        beta = beta * 2.0
      } else {
        beta = (beta + max_beta) / 2.0
      }
    } else {
      // Entropy is too low, need to decrease beta (increase sigma)
      max_beta = beta
      if min_beta == -1.0e300 {
        beta = beta / 2.0
      } else {
        beta = (beta + min_beta) / 2.0
      }
    }
    
    iter_count = iter_count + 1
  }
  
  p
}

///|
/// Compute joint probability matrix P from distance matrix.
fn compute_p_matrix(dist_matrix : Array[Array[Double]], perplexity : Double) -> Array[Array[Double]] {
  let n = dist_matrix.length()
  
  let p_matrix : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      row.push(0.0)
      j = j + 1
    }
    p_matrix.push(row)
    i = i + 1
  }
  
  // Compute P_{j|i} for each i
  let mut i = 0
  while i < n {
    let p_ji = compute_p_ji(dist_matrix[i], perplexity, n)
    let mut j = 0
    while j < n {
      p_matrix[i][j] = p_ji[j]
      j = j + 1
    }
    i = i + 1
  }
  
  // Symmetrize: P_{ij} = (P_{j|i} + P_{i|j}) / (2n)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      p_matrix[i][j] = (p_matrix[i][j] + p_matrix[j][i]) / (2.0 * n.to_double())
      j = j + 1
    }
    i = i + 1
  }
  
  p_matrix
}

// ============================================================
// t-SNE main algorithm
// ============================================================

///|
/// Run t-SNE on input data.
/// Input: data matrix (rows = samples, cols = features)
/// Output: t-SNE embedding (n x dims)
pub fn tsne(data : Array[Array[Double]], config : TsneConfig) -> TsneResult {
  let n = data.length()
  
  if n == 0 {
    return TsneResult::new()
  }
  
  if n < 2 {
    let result = TsneResult::new()
    let embedding : Array[Array[Double]] = Array::new()
    let row : Array[Double] = Array::new()
    let mut i = 0
    while i < config.dims {
      row.push(0.0)
      i = i + 1
    }
    embedding.push(row)
    return result
  }
  
  // Step 1: Compute pairwise distances
  let dist_matrix = calc_distance_matrix(data)
  
  // Step 2: Compute joint probability matrix P
  let p_matrix = compute_p_matrix(dist_matrix, config.perplexity)
  
  // Step 3: Initialize Y (the embedding) with small random values
  let y : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < config.dims {
      // Simple pseudo-random initialization (deterministic for reproducibility)
      let seed_val = (i * config.dims + j + config.random_seed).to_double()
      let rand_val = pseudo_random(seed_val) * 0.0001
      row.push(rand_val)
      j = j + 1
    }
    y.push(row)
    i = i + 1
  }
  
  // Step 4: Gradient descent
  let d_y : Array[Array[Double]] = Array::new()  // gradient
  let u_y : Array[Array[Double]] = Array::new()  // update (momentum)
  let gains : Array[Array[Double]] = Array::new()  // learning rate gains
  
  let mut i = 0
  while i < n {
    let d_row : Array[Double] = Array::new()
    let u_row : Array[Double] = Array::new()
    let g_row : Array[Double] = Array::new()
    let mut j = 0
    while j < config.dims {
      d_row.push(0.0)
      u_row.push(0.0)
      g_row.push(1.0)
      j = j + 1
    }
    d_y.push(d_row)
    u_y.push(u_row)
    gains.push(g_row)
    i = i + 1
  }
  
  let costs : Array[Double] = Array::new()
  let mut momentum = config.momentum
  
  let mut iter = 0
  while iter < config.max_iter {
    // Compute Q matrix (low-dimensional similarities)
    let q_matrix : Array[Array[Double]] = Array::new()
    let mut sum_q = 0.0
    
    let mut i = 0
    while i < n {
      let row : Array[Double] = Array::new()
      let mut j = 0
      while j < n {
        row.push(0.0)
        j = j + 1
      }
      q_matrix.push(row)
      i = i + 1
    }
    
    let mut i = 0
    while i < n {
      let mut j = i + 1
      while j < n {
        // Student-t distribution (1 degree of freedom)
        let mut sq_dist = 0.0
        let mut d = 0
        while d < config.dims {
          let diff = y[i][d] - y[j][d]
          sq_dist = sq_dist + diff * diff
          d = d + 1
        }
        
        let q_val = 1.0 / (1.0 + sq_dist)
        q_matrix[i][j] = q_val
        q_matrix[j][i] = q_val
        sum_q = sum_q + 2.0 * q_val
        
        j = j + 1
      }
      i = i + 1
    }
    
    if sum_q == 0.0 {
      sum_q = 1.0e-10
    }
    
    // Normalize Q
    let mut i = 0
    while i < n {
      let mut j = 0
      while j < n {
        q_matrix[i][j] = q_matrix[i][j] / sum_q
        j = j + 1
      }
      i = i + 1
    }
    
    // Compute cost (KL divergence)
    let mut cost = 0.0
    let mut i = 0
    while i < n {
      let mut j = 0
      while j < n {
        if i != j && p_matrix[i][j] > 0.0 {
          cost = cost + p_matrix[i][j] * @math.ln(p_matrix[i][j] / q_matrix[i][j])
        }
        j = j + 1
      }
      i = i + 1
    }
    costs.push(cost)
    
    // Compute gradient
    let mut i = 0
    while i < n {
      let mut d = 0
      while d < config.dims {
        d_y[i][d] = 0.0
        d = d + 1
      }
      
      let mut j = 0
      while j < n {
        if i != j {
          let mut sq_dist = 0.0
          let mut d = 0
          while d < config.dims {
            let diff = y[i][d] - y[j][d]
            sq_dist = sq_dist + diff * diff
            d = d + 1
          }
          
          let mult = (p_matrix[i][j] - q_matrix[i][j]) * (1.0 / (1.0 + sq_dist))
          
          let mut d = 0
          while d < config.dims {
            d_y[i][d] = d_y[i][d] + 4.0 * mult * (y[i][d] - y[j][d])
            d = d + 1
          }
        }
        j = j + 1
      }
      i = i + 1
    }
    
    // Update momentum
    if iter == config.mom_switch_iter {
      momentum = config.final_momentum
    }
    
    // Stop exaggeration
    if iter == config.stop_lying_iter {
      let mut i = 0
      while i < n {
        let mut j = 0
        while j < n {
          p_matrix[i][j] = p_matrix[i][j] / config.exaggeration_factor
          j = j + 1
        }
        i = i + 1
      }
    }
    
    // Apply early exaggeration
    if iter < config.stop_lying_iter {
      let mut i = 0
      while i < n {
        let mut j = 0
        while j < n {
          p_matrix[i][j] = p_matrix[i][j] * config.exaggeration_factor
          j = j + 1
        }
        i = i + 1
      }
    }
    
    // Update gains and apply momentum
    let mut i = 0
    while i < n {
      let mut d = 0
      while d < config.dims {
        // Update gain
        if (d_y[i][d] > 0.0 && u_y[i][d] > 0.0) || (d_y[i][d] < 0.0 && u_y[i][d] < 0.0) {
          gains[i][d] = gains[i][d] * 0.8
        } else {
          gains[i][d] = gains[i][d] + 0.2
        }
        
        if gains[i][d] < 0.01 {
          gains[i][d] = 0.01
        }
        
        // Update with momentum
        u_y[i][d] = momentum * u_y[i][d] - config.eta * gains[i][d] * d_y[i][d]
        y[i][d] = y[i][d] + u_y[i][d]
        
        d = d + 1
      }
      i = i + 1
    }
    
    iter = iter + 1
  }
  
  TsneResult::{
    embedding: y,
    costs: costs,
    n_iter: config.max_iter,
  }
}

// ============================================================
// Helper functions
// ============================================================

fn pseudo_random(seed : Double) -> Double {
  // Simple deterministic pseudo-random number generator
  let x = seed * 12.9898
  let y_val = seed * 78.233
  let z = x * y_val
  let sin_val = @math.sin(z)
  let result = sin_val * 43758.5453
  result - result.floor()
}

///|
/// Create sample data for testing/demo.
pub fn create_tsne_test_data(n_samples : Int, n_features : Int) -> Array[Array[Double]] {
  let data : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n_samples {
    let row : Array[Double] = Array::new()
    let group = i % 3  // 3 clusters
    
    let mut j = 0
    while j < n_features {
      let seed_val = (i * n_features + j).to_double()
      let noise = (pseudo_random(seed_val) - 0.5) * 2.0
      let cluster_center = group.to_double() * 3.0
      row.push(cluster_center + noise)
      j = j + 1
    }
    
    data.push(row)
    i = i + 1
  }
  
  data
}
