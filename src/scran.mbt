///|
/// scran module, inspired by Bioconductor scran.
/// Provides functions for single-cell RNA-seq data normalization and clustering.

///|
pub struct SizeFactors {
  factors : Array[Double]
  method : String
}

///|
pub struct LeidenClusterResult {
  clusters : Array[Int]
  resolution : Double
  n_clusters : Int
}

///|
pub struct ScoreMarkersResult {
  gene : Array[String]
  cluster : Array[Int]
  avg_log2_fc : Array[Double]
  p_value : Array[Double]
  fdr : Array[Double]
}

///|
pub fn compute_sum_factors(
  counts : Array[Array[Int]],
  sizes : Array[Int],
) -> SizeFactors {
  let n_cells = counts.length()
  let n_genes = if n_cells > 0 { counts[0].length() } else { 0 }
  
  let size_factors : Array[Double] = Array::make(n_cells, 1.0)
  
  let effective_sizes = if sizes.length() == 0 {
    let default_size = if n_cells > 0 { (n_cells / 20).max(10).min(100) } else { 10 }
    [default_size]
  } else {
    sizes
  }
  
  let target_size = effective_sizes[0]
  
  let all_sfs : Array[Array[Double]] = Array::new()
  
  let mut iteration = 0
  while iteration < 3 {
    let current_sfs : Array[Double] = Array::make(n_cells, 1.0)
    
    if iteration == 0 {
      let mut i = 0
      while i < n_cells {
        let mut sum = 0.0
        let mut j = 0
        while j < n_genes {
          sum = sum + counts[i][j].to_double()
          j = j + 1
        }
        current_sfs[i] = if sum > 0.0 { sum / 10000.0 } else { 1.0 }
        i = i + 1
      }
    } else {
      let mut i = 0
      while i < n_cells {
        current_sfs[i] = size_factors[i]
        i = i + 1
      }
    }
    
    let ranks : Array[Array[(Double, Int)]] = Array::new()
    let mut i = 0
    while i < n_genes {
      let gene_ranks : Array[(Double, Int)] = Array::new()
      let mut j = 0
      while j < n_cells {
        let norm_val = if current_sfs[j] > 0.0 { counts[j][i].to_double() / current_sfs[j] } else { 0.0 }
        gene_ranks.push((norm_val, j))
        j = j + 1
      }
      ranks.push(gene_ranks)
      i = i + 1
    }
    
    let mut j = 0
    while j < ranks.length() {
      let mut k = 0
      while k < ranks[j].length() - 1 {
        let mut l = k + 1
        while l < ranks[j].length() {
          if ranks[j][k].0 < ranks[j][l].0 {
            let tmp = ranks[j][k]
            ranks[j][k] = ranks[j][l]
            ranks[j][l] = tmp
          }
          l = l + 1
        }
        k = k + 1
      }
      j = j + 1
    }
    
    let subsets : Array[Array[Int]] = Array::new()
    let mut seed = 0
    while subsets.length() < 100 {
      let subset : Array[Int] = Array::new()
      let used : Map[Int, Bool] = Map([], capacity=n_cells)
      
      let mut m = 0
      while m < target_size && m < n_cells {
        let idx = (seed * 9301 + 49297) % 233280 % n_cells
        let used_idx = used.get(idx)
        if used_idx.is_none() || !used_idx.unwrap() {
          subset.push(idx)
          used[idx] = true
        }
        seed = seed + 1
        m = m + 1
      }
      
      if subset.length() == target_size {
        subsets.push(subset)
      }
    }
    
    let sf_subset : Array[Double] = Array::make(subsets.length(), 0.0)
    let mut m = 0
    while m < subsets.length() {
      let subset = subsets[m]
      let mut sum_subset = 0.0
      let mut n_non_zero = 0
      
      let mut j = 0
      while j < n_genes {
        let mut min_val = 1000000000000000000000000000000.0
        let mut k = 0
        while k < subset.length() {
          let cell_idx = subset[k]
          let norm_val = if current_sfs[cell_idx] > 0.0 { 
            counts[cell_idx][j].to_double() / current_sfs[cell_idx] 
          } else { 
            0.0 
          }
          if norm_val < min_val && norm_val > 0.0 {
            min_val = norm_val
          }
          k = k + 1
        }
        
        if min_val < 1000000000000000000000000000000.0 {
          sum_subset = sum_subset + min_val
          n_non_zero = n_non_zero + 1
        }
        j = j + 1
      }
      
      sf_subset[m] = if n_non_zero > 0 { sum_subset / n_non_zero.to_double() } else { 1.0 }
      m = m + 1
    }
    
    let mut i = 0
    while i < n_cells {
      let mut sum_sf = 0.0
      let mut count = 0
      
      let mut m = 0
      while m < subsets.length() {
        let subset = subsets[m]
        let mut found = false
        let mut k = 0
        while k < subset.length() && !found {
          if subset[k] == i {
            found = true
          }
          k = k + 1
        }
        
        if found {
          sum_sf = sum_sf + sf_subset[m]
          count = count + 1
        }
        m = m + 1
      }
      
      if count > 0 {
        size_factors[i] = sum_sf / count.to_double()
      } else {
        size_factors[i] = 1.0
      }
      i = i + 1
    }
    
    all_sfs.push(size_factors.clone())
    iteration = iteration + 1
  }
  
  let final_sfs : Array[Double] = Array::make(n_cells, 0.0)
  let mut i = 0
  while i < n_cells {
    let mut sum = 0.0
    let mut j = 0
    while j < all_sfs.length() {
      sum = sum + all_sfs[j][i]
      j = j + 1
    }
    final_sfs[i] = sum / all_sfs.length().to_double()
    i = i + 1
  }
  
  let ref_factor = final_sfs[0]
  i = 0
  while i < n_cells {
    if ref_factor > 0.0 {
      final_sfs[i] = final_sfs[i] / ref_factor
    }
    i = i + 1
  }
  
  SizeFactors::{
    factors: final_sfs,
    method: "sum_factors".to_string(),
  }
}

///|
pub fn scran_log_normalize(
  counts : Array[Array[Int]],
  size_factors : Array[Double],
) -> Array[Array[Double]] {
  let n_cells = counts.length()
  let n_genes = if n_cells > 0 { counts[0].length() } else { 0 }
  
  let normalized : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n_cells {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n_genes {
      let sf = if i < size_factors.length() { size_factors[i] } else { 1.0 }
      let norm = if sf > 0.0 { counts[i][j].to_double() / sf } else { 0.0 }
      let log_val = if norm > 0.0 { @math.ln(norm + 1.0) } else { 0.0 }
      row.push(log_val)
      j = j + 1
    }
    normalized.push(row)
    i = i + 1
  }
  
  normalized
}

///|
pub fn build_snn_graph(
  normalized : Array[Array[Double]],
  k : Int,
) -> (Array[Array[Int]], Array[Array[Double]]) {
  let n_cells = normalized.length()
  let n_genes = if n_cells > 0 { normalized[0].length() } else { 0 }
  
  let adjacency : Array[Array[Int]] = Array::make(n_cells, [])
  let weights : Array[Array[Double]] = Array::make(n_cells, [])
  
  let mut i = 0
  while i < n_cells {
    let distances : Array[(Double, Int)] = Array::new()
    
    let mut j = 0
    while j < n_cells {
      let mut sum_sq = 0.0
      let mut k_idx = 0
      while k_idx < n_genes {
        let diff = normalized[i][k_idx] - normalized[j][k_idx]
        sum_sq = sum_sq + diff * diff
        k_idx = k_idx + 1
      }
      distances.push((sum_sq.sqrt(), j))
      j = j + 1
    }
    
    let mut j = 0
    while j < distances.length() - 1 {
      let mut k_idx = j + 1
      while k_idx < distances.length() {
        if distances[j].0 > distances[k_idx].0 {
          let tmp = distances[j]
          distances[j] = distances[k_idx]
          distances[k_idx] = tmp
        }
        k_idx = k_idx + 1
      }
      j = j + 1
    }
    
    let k_neighbors = if k < distances.length() { k } else { distances.length() }
    let neighbors : Array[Int] = Array::new()
    let neighbor_weights : Array[Double] = Array::new()
    
    let mut j = 1
    while j < k_neighbors {
      let idx = distances[j].1
      neighbors.push(idx)
      neighbor_weights.push(1.0 - distances[j].0 / distances[k_neighbors - 1].0)
      j = j + 1
    }
    
    adjacency[i] = neighbors
    weights[i] = neighbor_weights
    
    i = i + 1
  }
  
  (adjacency, weights)
}

///|
pub fn leiden_clustering(
  adjacency : Array[Array[Int]],
  weights : Array[Array[Double]],
  resolution : Double,
) -> LeidenClusterResult {
  let n_cells = adjacency.length()
  
  let clusters : Array[Int] = Array::make(n_cells, 0)
  let mut i = 0
  while i < n_cells {
    clusters[i] = i
    i = i + 1
  }
  
  let mut improved = true
  let mut iter = 0
  while improved && iter < 100 {
    improved = false
    
    let mut i = 0
    while i < n_cells {
      let current_cluster = clusters[i]
      let mut best_cluster = current_cluster
      let mut best_gain = 0.0
      
      let considered_clusters : Map[Int, Bool] = Map([], capacity=10)
      considered_clusters[current_cluster] = true
      
      let mut j = 0
      while j < adjacency[i].length() {
        let neighbor = adjacency[i][j]
        considered_clusters[clusters[neighbor]] = true
        j = j + 1
      }
      
      let keys_array : Array[Int] = Array::new()
      let keys_iter = considered_clusters.keys().collect()
      let mut k = 0
      while k < keys_iter.length() {
        keys_array.push(keys_iter[k])
        k = k + 1
      }
      let mut m = 0
      while m < keys_array.length() {
        let cluster_id = keys_array[m]
        let mut gain = 0.0
        let mut n_in_cluster = 0
        
        let mut k = 0
        while k < n_cells {
          if clusters[k] == cluster_id {
            n_in_cluster = n_in_cluster + 1
            
            let mut found = false
            let mut l = 0
            while l < adjacency[i].length() && !found {
              if adjacency[i][l] == k {
                gain = gain + weights[i][l] * resolution
                found = true
              }
              l = l + 1
            }
          }
          k = k + 1
        }
        
        if n_in_cluster > 0 {
          let cluster_factor = n_in_cluster.to_double() / n_cells.to_double()
          gain = gain - cluster_factor * gain
        }
        
        if cluster_id == current_cluster {
          let mut n_current = 0
          let mut k = 0
          while k < n_cells {
            if clusters[k] == current_cluster {
              n_current = n_current + 1
            }
            k = k + 1
          }
          
          let mut current_gain = 0.0
          let mut k = 0
          while k < adjacency[i].length() {
            let neighbor = adjacency[i][k]
            if clusters[neighbor] == current_cluster {
              current_gain = current_gain + weights[i][k] * resolution
            }
            k = k + 1
          }
          
          gain = current_gain - gain
        }
        
        if gain > best_gain {
          best_gain = gain
          best_cluster = cluster_id
        }
        m = m + 1
      }
      
      if best_cluster != current_cluster {
        clusters[i] = best_cluster
        improved = true
      }
      
      i = i + 1
    }
    
    iter = iter + 1
  }
  
  let cluster_ids : Map[Int, Int] = Map([], capacity=n_cells)
  let mut next_id = 0
  let mut i = 0
  while i < n_cells {
    let mapped = cluster_ids.get(clusters[i])
    if mapped.is_none() {
      cluster_ids[clusters[i]] = next_id
      next_id = next_id + 1
    }
    clusters[i] = cluster_ids[clusters[i]]
    i = i + 1
  }
  
  LeidenClusterResult::{
    clusters,
    resolution,
    n_clusters: next_id,
  }
}

///|
pub fn score_markers(
  counts : Array[Array[Int]],
  clusters : Array[Int],
  gene_names : Array[String],
) -> ScoreMarkersResult {
  let n_cells = counts.length()
  let n_genes = if n_cells > 0 { counts[0].length() } else { 0 }
  
  let cluster_sizes : Map[Int, Int] = Map([], capacity=10)
  let mut i = 0
  while i < n_cells {
    let current = cluster_sizes.get(clusters[i])
    if current.is_some() {
      cluster_sizes[clusters[i]] = current.unwrap() + 1
    } else {
      cluster_sizes[clusters[i]] = 1
    }
    i = i + 1
  }
  
  let gene_result : Array[String] = Array::new()
  let cluster_result : Array[Int] = Array::new()
  let avg_log2_fc_result : Array[Double] = Array::new()
  let p_value_result : Array[Double] = Array::new()
  
  i = 0
  while i < n_genes {
    let mut j = 0
    while j < cluster_sizes.length() {
      let cluster_id = j
      
      let mut cluster_sum = 0.0
      let mut other_sum = 0.0
      let mut cluster_count = 0
      let mut other_count = 0
      
      let mut k = 0
      while k < n_cells {
        if clusters[k] == cluster_id {
          cluster_sum = cluster_sum + counts[k][i].to_double()
          cluster_count = cluster_count + 1
        } else {
          other_sum = other_sum + counts[k][i].to_double()
          other_count = other_count + 1
        }
        k = k + 1
      }
      
      let cluster_avg = if cluster_count > 0 { cluster_sum / cluster_count.to_double() } else { 0.0 }
      let other_avg = if other_count > 0 { other_sum / other_count.to_double() } else { 0.0 }
      
      let log2_fc = if other_avg > 0.0 { 
        @math.ln((cluster_avg + 1.0) / (other_avg + 1.0)) / @math.ln(2.0) 
      } else if cluster_avg > 0.0 {
        10.0
      } else {
        0.0
      }
      
      let p_val = scran_hypergeometric_test(counts, i, clusters, cluster_id)
      
      if i < gene_names.length() {
        gene_result.push(gene_names[i])
      } else {
        gene_result.push("Gene_" + i.to_string())
      }
      cluster_result.push(cluster_id)
      avg_log2_fc_result.push(log2_fc)
      p_value_result.push(p_val)
      
      j = j + 1
    }
    i = i + 1
  }
  
  let fdr_result = benjamini_hochberg_markers(p_value_result)
  
  ScoreMarkersResult::{
    gene: gene_result,
    cluster: cluster_result,
    avg_log2_fc: avg_log2_fc_result,
    p_value: p_value_result,
    fdr: fdr_result,
  }
}

///|
fn scran_hypergeometric_test(
  counts : Array[Array[Int]],
  gene_idx : Int,
  clusters : Array[Int],
  cluster_id : Int,
) -> Double {
  let n_cells = counts.length()
  
  let mut cluster_total = 0
  let mut other_total = 0
  let mut cluster_gene = 0
  let mut other_gene = 0
  
  let mut i = 0
  while i < n_cells {
    let count = counts[i][gene_idx]
    if clusters[i] == cluster_id {
      cluster_total = cluster_total + 1
      if count > 0 {
        cluster_gene = cluster_gene + 1
      }
    } else {
      other_total = other_total + 1
      if count > 0 {
        other_gene = other_gene + 1
      }
    }
    i = i + 1
  }
  
  let total = cluster_total + other_total
  let gene_total = cluster_gene + other_gene
  
  if gene_total == 0 || gene_total == total {
    return 1.0
  }
  
  let p_val = if cluster_gene > gene_total * cluster_total / total {
    let mut p = 0.0
    let mut k = cluster_gene
    while k <= gene_total {
      let term = combination(gene_total, k) * combination(total - gene_total, cluster_total - k) / combination(total, cluster_total)
      p = p + term
      k = k + 1
    }
    p
  } else {
    let mut p = 1.0
    let mut k = 0
    while k < cluster_gene {
      let term = combination(gene_total, k) * combination(total - gene_total, cluster_total - k) / combination(total, cluster_total)
      p = p - term
      k = k + 1
    }
    p
  }
  
  if p_val > 1.0 { 1.0 } else if p_val < 0.0 { 0.0 } else { p_val }
}

///|
fn combination(n : Int, k : Int) -> Double {
  if k < 0 || k > n {
    return 0.0
  }
  if k == 0 || k == n {
    return 1.0
  }
  
  let k = k.min(n - k)
  
  let mut result = 1.0
  let mut i = 0
  while i < k {
    result = result * (n - i).to_double() / (i + 1).to_double()
    i = i + 1
  }
  
  result
}

///|
fn benjamini_hochberg_markers(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  
  let indexed : Array[(Double, Int)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((p_values[i], i))
    i = i + 1
  }
  
  let mut j = 0
  while j < indexed.length() - 1 {
    let mut k = j + 1
    while k < indexed.length() {
      if indexed[j].0 > indexed[k].0 {
        let tmp = indexed[j]
        indexed[j] = indexed[k]
        indexed[k] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }
  
  let fdr : Array[Double] = Array::make(n, 1.0)
  let mut min_fdr = 1.0
  
  i = n - 1
  while i >= 0 {
    let p = indexed[i].0
    let rank = (i + 1).to_double()
    let current_fdr = p * n.to_double() / rank
    
    if current_fdr < min_fdr {
      min_fdr = current_fdr
    }
    
    fdr[indexed[i].1] = min_fdr
    
    i = i - 1
  }
  
  fdr
}

///|
pub fn create_example_scran_data() -> (Array[Array[Int]], Array[String]) {
  let n_cells = 100
  let n_genes = 500
  
  let counts : Array[Array[Int]] = Array::new()
  let gene_names : Array[String] = Array::new()
  
  let mut i = 0
  while i < n_genes {
    gene_names.push("Gene_" + i.to_string())
    i = i + 1
  }
  
  let mut i = 0
  while i < n_cells {
    let row : Array[Int] = Array::new()
    
    let cluster = i / 20
    
    let mut j = 0
    while j < n_genes {
      let mut base_prob = 0.01
      if j >= cluster * 50 && j < (cluster + 1) * 50 {
        base_prob = 0.1
      }
      
      let count = if (j * 7 + i * 13).to_double() % 10.0 / 10.0 < base_prob {
        ((j * 7 + i * 13).to_double() % 10.0 / 10.0 * 20.0).to_int()
      } else {
        0
      }
      row.push(count)
      j = j + 1
    }
    
    counts.push(row)
    i = i + 1
  }
  
  (counts, gene_names)
}
