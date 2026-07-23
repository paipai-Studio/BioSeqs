///|
/// SC3 module - Single-Cell Consensus Clustering.
/// Inspired by Bioconductor SC3.

///|
/// Represents an SC3 clustering result.
pub struct SC3Result {
  cluster_labels : Array[Int]
  consensus_matrix : Array[Array[Double]]
  k : Int
  silhouette_scores : Array[Double]
  gap_statistics : Array[Double]
  optimal_k : Int
}

///|
/// Create a new SC3Result.
pub fn SC3Result::new(
  cluster_labels : Array[Int],
  consensus_matrix : Array[Array[Double]],
  k : Int,
  silhouette_scores : Array[Double],
  gap_statistics : Array[Double],
  optimal_k : Int
) -> SC3Result {
  SC3Result::{ cluster_labels, consensus_matrix, k, silhouette_scores, gap_statistics, optimal_k }
}

///|
/// Preprocess single-cell data.
pub fn sc3_preprocess(data : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = data.length()
  let dim = if n > 0 { data[0].length() } else { 0 }
  
  let processed : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut sum = 0.0
    let mut j = 0
    while j < dim {
      sum = sum + data[i][j]
      j = j + 1
    }
    
    let mut k = 0
    while k < dim {
      row.push(if sum > 0.0 { data[i][k] / sum * 1000000.0 } else { 0.0 })
      k = k + 1
    }
    
    processed.push(row)
    i = i + 1
  }
  
  processed
}

///|
/// Compute PCA on preprocessed data.
pub fn sc3_pca(data : Array[Array[Double]], n_components : Int) -> Array[Array[Double]] {
  let n = data.length()
  let dim = if n > 0 { data[0].length() } else { 0 }
  
  let means : Array[Double] = Array::new()
  let mut i = 0
  while i < dim {
    let mut sum = 0.0
    let mut j = 0
    while j < n {
      sum = sum + data[j][i]
      j = j + 1
    }
    means.push(sum / n.to_double())
    i = i + 1
  }
  
  let centered : Array[Array[Double]] = Array::new()
  let mut k = 0
  while k < n {
    let row : Array[Double] = Array::new()
    let mut l = 0
    while l < dim {
      row.push(data[k][l] - means[l])
      l = l + 1
    }
    centered.push(row)
    k = k + 1
  }
  
  let covariance : Array[Array[Double]] = Array::new()
  let mut m = 0
  while m < dim {
    let row : Array[Double] = Array::new()
    let mut p = 0
    while p < dim {
      let mut sum = 0.0
      let mut q = 0
      while q < n {
        sum = sum + centered[q][m] * centered[q][p]
        q = q + 1
      }
      row.push(sum / (n - 1).to_double())
      p = p + 1
    }
    covariance.push(row)
    m = m + 1
  }
  
  let eigenvalues = sc3_power_iteration(covariance, n_components)
  
  let components : Array[Array[Double]] = Array::new()
  let mut r = 0
  while r < n_components && r < eigenvalues.length() {
    let component : Array[Double] = Array::new()
    let mut s = 0
    while s < n {
      let mut sum = 0.0
      let mut t = 0
      while t < dim {
        sum = sum + centered[s][t] * eigenvalues[r][t]
        t = t + 1
      }
      component.push(sum)
      s = s + 1
    }
    components.push(component)
    r = r + 1
  }
  
  let transposed : Array[Array[Double]] = Array::new()
  let mut u = 0
  while u < n {
    let row : Array[Double] = Array::new()
    let mut v = 0
    while v < components.length() {
      row.push(components[v][u])
      v = v + 1
    }
    transposed.push(row)
    u = u + 1
  }
  
  transposed
}

///|
/// Power iteration for eigenvalue decomposition.
pub fn sc3_power_iteration(matrix : Array[Array[Double]], n_components : Int) -> Array[Array[Double]] {
  let n = matrix.length()
  let components : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n_components {
    let v : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      v.push(((i + 1).to_double() / (j + 1).to_double()) - 0.5)
      j = j + 1
    }
    
    let mut iter = 0
    while iter < 100 {
      let v_new : Array[Double] = Array::new()
      let mut k = 0
      while k < n {
        let mut sum = 0.0
        let mut l = 0
        while l < n {
          sum = sum + matrix[k][l] * v[l]
          l = l + 1
        }
        v_new.push(sum)
        k = k + 1
      }
      
      let mut norm = 0.0
      let mut m = 0
      while m < n {
        norm = norm + v_new[m] * v_new[m]
        m = m + 1
      }
      norm = norm.sqrt()
      
      let mut p = 0
      while p < n {
        v.set(p, v_new[p] / norm)
        p = p + 1
      }
      
      iter = iter + 1
    }
    
    components.push(v)
    
    let mut q = 0
    while q < n {
      let mut r = 0
      while r < n {
        let mut sum = 0.0
        let mut s = 0
        while s < n {
          sum = sum + v[s] * matrix[s][r]
          s = s + 1
        }
        matrix[q].set(r, matrix[q][r] - v[q] * sum)
        r = r + 1
      }
      q = q + 1
    }
    
    i = i + 1
  }
  
  components
}

///|
/// Calculate consensus matrix.
pub fn sc3_calculate_consensus_matrix(
  data : Array[Array[Double]],
  k : Int,
  n_iterations : Int
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
  
  let mut iter = 0
  while iter < n_iterations {
    let labels = sc3_kmeans(data, k)
    
    let mut a = 0
    while a < n {
      let mut b = 0
      while b < n {
        if labels[a] == labels[b] {
          consensus[a].set(b, consensus[a][b] + 1.0)
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
      consensus[x].set(y, consensus[x][y] / n_iterations.to_double())
      y = y + 1
    }
    x = x + 1
  }
  
  consensus
}

///|
/// Perform k-means clustering.
pub fn sc3_kmeans(data : Array[Array[Double]], k : Int) -> Array[Int] {
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
        let dist = sc3_euclidean_distance(data[p], centroids[q])
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
/// Calculate Euclidean distance.
pub fn sc3_euclidean_distance(v1 : Array[Double], v2 : Array[Double]) -> Double {
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
/// Calculate silhouette scores.
pub fn sc3_calculate_silhouette(
  data : Array[Array[Double]],
  labels : Array[Int],
  k : Int
) -> Array[Double] {
  let n = data.length()
  let scores : Array[Double] = Array::new()
  
  let mut i = 0
  while i < n {
    let label = labels[i]
    
    let mut a_sum = 0.0
    let mut a_count = 0
    let mut j = 0
    while j < n {
      if labels[j] == label && i != j {
        a_sum = a_sum + sc3_euclidean_distance(data[i], data[j])
        a_count = a_count + 1
      }
      j = j + 1
    }
    let a = if a_count > 0 { a_sum / a_count.to_double() } else { 0.0 }
    
    let mut b_min = 10000000000.0
    let mut c = 0
    while c < k {
      if c != label {
        let mut b_sum = 0.0
        let mut b_count = 0
        let mut d = 0
        while d < n {
          if labels[d] == c {
            b_sum = b_sum + sc3_euclidean_distance(data[i], data[d])
            b_count = b_count + 1
          }
          d = d + 1
        }
        if b_count > 0 {
          let b_val = b_sum / b_count.to_double()
          if b_val < b_min {
            b_min = b_val
          }
        }
      }
      c = c + 1
    }
    
    let max_ab = if a > b_min { a } else { b_min }
    if max_ab > 0.0 {
      scores.push((b_min - a) / max_ab)
    } else {
      scores.push(0.0)
    }
    
    i = i + 1
  }
  
  scores
}

///|
/// Calculate gap statistics.
pub fn sc3_calculate_gap_statistics(
  data : Array[Array[Double]],
  k : Int
) -> Array[Double] {
  let n = data.length()
  
  let gaps : Array[Double] = Array::new()
  
  let mut current_k = 1
  while current_k <= k {
    let w = sc3_calculate_wcss(data, current_k)
    
    let mut w_null_sum = 0.0
    let mut b = 0
    while b < 10 {
      let null_data = sc3_generate_null_data(data, b)
      let w_null = sc3_calculate_wcss(null_data, current_k)
      w_null_sum = w_null_sum + @math.ln(w_null)
      b = b + 1
    }
    
    gaps.push(w_null_sum / 10.0 - @math.ln(w))
    
    current_k = current_k + 1
  }
  
  gaps
}

///|
/// Calculate within-cluster sum of squares.
pub fn sc3_calculate_wcss(data : Array[Array[Double]], k : Int) -> Double {
  let labels = sc3_kmeans(data, k)
  let n = data.length()
  let dim = if n > 0 { data[0].length() } else { 0 }
  
  let centroids : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < k {
    let centroid : Array[Double] = Array::new()
    let mut j = 0
    while j < dim {
      centroid.push(0.0)
      j = j + 1
    }
    centroids.push(centroid)
    i = i + 1
  }
  
  let counts : Array[Int] = Array::new()
  let mut c = 0
  while c < k {
    counts.push(0)
    c = c + 1
  }
  
  let mut m = 0
  while m < n {
    let label = labels[m]
    let mut p = 0
    while p < dim {
      centroids[label].set(p, centroids[label][p] + data[m][p])
      p = p + 1
    }
    counts.set(label, counts[label] + 1)
    m = m + 1
  }
  
  let mut q = 0
  while q < k {
    if counts[q] > 0 {
      let mut r = 0
      while r < dim {
        centroids[q].set(r, centroids[q][r] / counts[q].to_double())
        r = r + 1
      }
    }
    q = q + 1
  }
  
  let mut wcss = 0.0
  let mut s = 0
  while s < n {
    let label = labels[s]
    wcss = wcss + sc3_euclidean_distance(data[s], centroids[label]).pow(2)
    s = s + 1
  }
  
  wcss
}

///|
/// Generate null data for gap statistics.
pub fn sc3_generate_null_data(data : Array[Array[Double]], seed : Int) -> Array[Array[Double]] {
  let n = data.length()
  let dim = if n > 0 { data[0].length() } else { 0 }
  
  let mins : Array[Double] = Array::new()
  let maxes : Array[Double] = Array::new()
  
  let mut j = 0
  while j < dim {
    let mut min_val = 10000000000.0
    let mut max_val = -10000000000.0
    let mut i = 0
    while i < n {
      let val = data[i][j]
      min_val = if val < min_val { val } else { min_val }
      max_val = if val > max_val { val } else { max_val }
      i = i + 1
    }
    mins.push(min_val)
    maxes.push(max_val)
    j = j + 1
  }
  
  let null_data : Array[Array[Double]] = Array::new()
  let mut k = 0
  while k < n {
    let row : Array[Double] = Array::new()
    let mut l = 0
    while l < dim {
      let state = (seed * 1103515245 + k * 12345 + l * 98765) >> 16 & 0x7fff
      let rand_val = state.to_double() / 32767.0
      row.push(mins[l] + rand_val * (maxes[l] - mins[l]))
      l = l + 1
    }
    null_data.push(row)
    k = k + 1
  }
  
  null_data
}

///|
/// Perform SC3 consensus clustering.
pub fn sc3_consensus_cluster(
  data : Array[Array[Double]],
  k : Int,
  n_iterations? : Int = 100,
  n_pcs? : Int = 50
) -> SC3Result {
  let processed = sc3_preprocess(data)
  let pcs = sc3_pca(processed, n_pcs)
  
  let consensus_matrix = sc3_calculate_consensus_matrix(pcs, k, n_iterations)
  let cluster_labels = sc3_kmeans(pcs, k)
  let silhouette_scores = sc3_calculate_silhouette(pcs, cluster_labels, k)
  let gap_statistics = sc3_calculate_gap_statistics(pcs, k)
  
  let mut optimal_k = k
  let mut max_gap_diff = 0.0
  let mut g = 1
  while g < gap_statistics.length() {
    let diff = gap_statistics[g] - gap_statistics[g-1]
    if diff > max_gap_diff {
      max_gap_diff = diff
      optimal_k = g + 1
    }
    g = g + 1
  }
  
  SC3Result::new(cluster_labels, consensus_matrix, k, silhouette_scores, gap_statistics, optimal_k)
}

///|
/// Bio API functions.

pub fn bio_sc3_cluster(data : Array[Array[Double]], k : Int) -> SC3Result {
  sc3_consensus_cluster(data, k)
}

pub fn bio_sc3_pca(data : Array[Array[Double]], n_components : Int) -> Array[Array[Double]] {
  let processed = sc3_preprocess(data)
  sc3_pca(processed, n_components)
}

pub fn bio_sc3_silhouette(data : Array[Array[Double]], labels : Array[Int], k : Int) -> Array[Double] {
  sc3_calculate_silhouette(data, labels, k)
}