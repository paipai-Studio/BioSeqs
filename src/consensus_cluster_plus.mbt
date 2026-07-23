///|
/// ConsensusClusterPlus module - Consensus clustering for cancer subtyping.
/// Inspired by Bioconductor ConsensusClusterPlus.

///|
/// Represents a consensus clustering result.
pub struct CCPResult {
  consensus_matrix : Array[Array[Double]]
  cluster_labels : Array[Int]
  k : Int
  consensus_score : Double
  stability_scores : Array[Double]
}

///|
/// Create a new CCPResult.
pub fn CCPResult::new(
  consensus_matrix : Array[Array[Double]],
  cluster_labels : Array[Int],
  k : Int,
  consensus_score : Double,
  stability_scores : Array[Double]
) -> CCPResult {
  CCPResult::{ consensus_matrix, cluster_labels, k, consensus_score, stability_scores }
}

///|
/// Calculate consensus matrix from clustering results.
pub fn ccp_calculate_consensus_matrix(
  data : Array[Array[Double]],
  k : Int,
  n_iterations : Int,
  sample_fraction : Double
) -> Array[Array[Double]] {
  let n = data.length()
  
  let consensus : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      row.push(0.0)
      j = j + 1
    }
    consensus.push(row)
    i = i + 1
  }
  
  let counts : Array[Array[Int]] = Array::new()
  let mut m = 0
  while m < n {
    let row : Array[Int] = Array::new()
    let mut p = 0
    while p < n {
      row.push(0)
      p = p + 1
    }
    counts.push(row)
    m = m + 1
  }
  
  let mut iter = 0
  while iter < n_iterations {
    let sample_size = (n.to_double() * sample_fraction).to_int()
    let sample_indices = ccp_deterministic_sample(n, sample_size, iter)
    
    let sample_data : Array[Array[Double]] = Array::new()
    let mut idx = 0
    while idx < sample_indices.length() {
      sample_data.push(data[sample_indices[idx]])
      idx = idx + 1
    }
    
    let labels = ccp_kmeans_cluster(sample_data, k)
    
    let mut a = 0
    while a < sample_indices.length() {
      let mut b = 0
      while b < sample_indices.length() {
        let si = sample_indices[a]
        let sj = sample_indices[b]
        counts[si].set(sj, counts[si][sj] + 1)
        if labels[a] == labels[b] {
          consensus[si].set(sj, consensus[si][sj] + 1.0)
        }
        b = b + 1
      }
      a = a + 1
    }
    
    iter = iter + 1
  }
  
  let mut x = 0
  while x < n {
    let mut y = 0
    while y < n {
      if counts[x][y] > 0 {
        consensus[x].set(y, consensus[x][y] / counts[x][y].to_double())
      } else {
        consensus[x].set(y, 0.0)
      }
      y = y + 1
    }
    x = x + 1
  }
  
  consensus
}

///|
/// Generate deterministic sample indices using seed.
pub fn ccp_deterministic_sample(n : Int, size : Int, seed : Int) -> Array[Int] {
  let indices : Array[Int] = Array::new()
  let mut i = 0
  while i < n {
    indices.push(i)
    i = i + 1
  }
  
  let mut state = seed
  let mut j = n - 1
  while j > 0 {
    let k_val = ((state * 1103515245 + 12345) >> 16) & 0x7fff
    let idx = (k_val.to_double() / 32767.0 * (j + 1).to_double()).to_int()
    state = k_val
    let temp = indices[j]
    indices.set(j, indices[idx])
    indices.set(idx, temp)
    j = j - 1
  }
  
  let result : Array[Int] = Array::new()
  let mut m = 0
  while m < size && m < n {
    result.push(indices[m])
    m = m + 1
  }
  
  result
}

///|
/// Perform k-means clustering with deterministic initialization.
pub fn ccp_kmeans_cluster(data : Array[Array[Double]], k : Int) -> Array[Int] {
  let n = data.length()
  let dim = if n > 0 { data[0].length() } else { 0 }
  
  let centroids : Array[Array[Double]] = Array::new()
  let mut c = 0
  while c < k && c < n {
    centroids.push(data[c])
    c = c + 1
  }
  
  let labels : Array[Int] = Array::new()
  let mut j = 0
  while j < n {
    labels.push(0)
    j = j + 1
  }
  
  let mut changed = true
  let mut iterations = 0
  while changed && iterations < 100 {
    changed = false
    
    let mut p = 0
    while p < n {
      let mut min_dist = 10000000000.0
      let mut best_cluster = 0
      
      let mut q = 0
      while q < centroids.length() {
        let dist = ccp_euclidean_distance(data[p], centroids[q])
        if dist < min_dist {
          min_dist = dist
          best_cluster = q
        }
        q = q + 1
      }
      
      if labels[p] != best_cluster {
        labels.set(p, best_cluster)
        changed = true
      }
      p = p + 1
    }
    
    let mut r = 0
    while r < centroids.length() {
      let sum : Array[Double] = Array::new()
      let mut s = 0
      while s < dim {
        sum.push(0.0)
        s = s + 1
      }
      
      let mut count = 0
      let mut t = 0
      while t < n {
        if labels[t] == r {
          let mut u = 0
          while u < dim {
            sum.set(u, sum[u] + data[t][u])
            u = u + 1
          }
          count = count + 1
        }
        t = t + 1
      }
      
      if count > 0 {
        let mut v = 0
        while v < dim {
          centroids[r].set(v, sum[v] / count.to_double())
          v = v + 1
        }
      }
      r = r + 1
    }
    
    iterations = iterations + 1
  }
  
  labels
}

///|
/// Calculate Euclidean distance between two vectors.
pub fn ccp_euclidean_distance(v1 : Array[Double], v2 : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < v1.length() && i < v2.length() {
    let diff = v1[i] - v2[i]
    sum = sum + diff * diff
    i = i + 1
  }
  sum.sqrt()
}

///|
/// Calculate consensus score from consensus matrix.
pub fn ccp_calculate_consensus_score(consensus_matrix : Array[Array[Double]]) -> Double {
  let n = consensus_matrix.length()
  if n == 0 {
    return 0.0
  }
  
  let mut sum = 0.0
  let mut count = 0
  
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      sum = sum + consensus_matrix[i][j]
      count = count + 1
      j = j + 1
    }
    i = i + 1
  }
  
  if count > 0 { sum / count.to_double() } else { 0.0 }
}

///|
/// Calculate stability scores for each cluster.
pub fn ccp_calculate_stability_scores(consensus_matrix : Array[Array[Double]], cluster_labels : Array[Int], k : Int) -> Array[Double] {
  let n = consensus_matrix.length()
  let scores : Array[Double] = Array::new()
  
  let mut c = 0
  while c < k {
    let members : Array[Int] = Array::new()
    let mut i = 0
    while i < n {
      if cluster_labels[i] == c {
        members.push(i)
      }
      i = i + 1
    }
    
    if members.length() < 2 {
      scores.push(0.0)
      c = c + 1
      continue
    }
    
    let mut sum = 0.0
    let mut count = 0
    let mut a = 0
    while a < members.length() {
      let mut b = a + 1
      while b < members.length() {
        sum = sum + consensus_matrix[members[a]][members[b]]
        count = count + 1
        b = b + 1
      }
      a = a + 1
    }
    
    scores.push(if count > 0 { sum / count.to_double() } else { 0.0 })
    c = c + 1
  }
  
  scores
}

///|
/// Perform consensus clustering.
pub fn ccp_consensus_cluster(
  data : Array[Array[Double]],
  k : Int,
  n_iterations? : Int = 100,
  sample_fraction? : Double = 0.8
) -> CCPResult {
  let consensus_matrix = ccp_calculate_consensus_matrix(data, k, n_iterations, sample_fraction)
  let cluster_labels = ccp_kmeans_cluster(data, k)
  let consensus_score = ccp_calculate_consensus_score(consensus_matrix)
  let stability_scores = ccp_calculate_stability_scores(consensus_matrix, cluster_labels, k)
  
  CCPResult::new(consensus_matrix, cluster_labels, k, consensus_score, stability_scores)
}

///|
/// Find optimal k by evaluating consensus scores across a range.
pub fn ccp_find_optimal_k(
  data : Array[Array[Double]],
  min_k : Int,
  max_k : Int,
  n_iterations? : Int = 100,
  sample_fraction? : Double = 0.8
) -> (Int, Array[Double]) {
  let scores : Array[Double] = Array::new()
  let mut best_k = min_k
  let mut best_score = 0.0
  
  let mut k = min_k
  while k <= max_k {
    let result = ccp_consensus_cluster(data, k, n_iterations = n_iterations, sample_fraction = sample_fraction)
    scores.push(result.consensus_score)
    
    if result.consensus_score > best_score {
      best_score = result.consensus_score
      best_k = k
    }
    
    k = k + 1
  }
  
  (best_k, scores)
}

///|
/// Bio API functions.

pub fn bio_consensus_cluster(data : Array[Array[Double]], k : Int) -> CCPResult {
  ccp_consensus_cluster(data, k)
}

pub fn bio_consensus_cluster_find_optimal_k(data : Array[Array[Double]], min_k : Int, max_k : Int) -> (Int, Array[Double]) {
  ccp_find_optimal_k(data, min_k, max_k)
}

pub fn bio_consensus_cluster_consensus_score(consensus_matrix : Array[Array[Double]]) -> Double {
  ccp_calculate_consensus_score(consensus_matrix)
}