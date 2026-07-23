///|
/// GENESIS module - Population structure and kinship analysis.
/// Inspired by Bioconductor GENESIS.

///|
/// Represents a kinship matrix result.
pub struct GSKinshipResult {
  kinship_matrix : Array[Array[Double]]
  sample_ids : Array[String]
  n_samples : Int
}

///|
/// Create a new GSKinshipResult.
pub fn GSKinshipResult::new(
  kinship_matrix : Array[Array[Double]],
  sample_ids : Array[String],
  n_samples : Int
) -> GSKinshipResult {
  GSKinshipResult::{ kinship_matrix, sample_ids, n_samples }
}

///|
/// Represents a PCA result.
pub struct GSPcaResult {
  eigenvalues : Array[Double]
  eigenvectors : Array[Array[Double]]
  var_explained : Array[Double]
  sample_ids : Array[String]
}

///|
/// Create a new GSPcaResult.
pub fn GSPcaResult::new(
  eigenvalues : Array[Double],
  eigenvectors : Array[Array[Double]],
  var_explained : Array[Double],
  sample_ids : Array[String]
) -> GSPcaResult {
  GSPcaResult::{ eigenvalues, eigenvectors, var_explained, sample_ids }
}

///|
/// Represents a genetic distance matrix result.
pub struct GSDistanceResult {
  distance_matrix : Array[Array[Double]]
  sample_ids : Array[String]
  distance_type : String
}

///|
/// Create a new GSDistanceResult.
pub fn GSDistanceResult::new(
  distance_matrix : Array[Array[Double]],
  sample_ids : Array[String],
  distance_type : String
) -> GSDistanceResult {
  GSDistanceResult::{ distance_matrix, sample_ids, distance_type }
}

///|
/// Estimate kinship matrix from genotype data.
pub fn gs_estimate_kinship(genotypes : Array[Array[Double]]) -> GSKinshipResult {
  let n = genotypes.length()
  let m = if n > 0 { genotypes[0].length() } else { 0 }
  
  let sample_ids : Array[String] = Array::new()
  let mut i = 0
  while i < n {
    sample_ids.push("Sample" + i.to_string())
    i = i + 1
  }
  
  let kinship : Array[Array[Double]] = Array::new()
  let mut j = 0
  while j < n {
    let row : Array[Double] = Array::new()
    let mut k = 0
    while k < n {
      row.push(0.0)
      k = k + 1
    }
    kinship.push(row)
    j = j + 1
  }
  
  let allele_freqs : Array[Double] = Array::new()
  let mut l = 0
  while l < m {
    let mut sum = 0.0
    let mut n_valid = 0
    let mut p = 0
    while p < n {
      if genotypes[p][l] >= 0.0 {
        sum = sum + genotypes[p][l]
        n_valid = n_valid + 1
      }
      p = p + 1
    }
    allele_freqs.push(if n_valid > 0 { sum / (2.0 * n_valid.to_double()) } else { 0.5 })
    l = l + 1
  }
  
  let mut a = 0
  while a < n {
    let mut b = a
    while b < n {
      let mut sum = 0.0
      let mut count = 0
      let mut c = 0
      while c < m {
        if genotypes[a][c] >= 0.0 && genotypes[b][c] >= 0.0 {
          let freq = allele_freqs[c]
          if freq > 0.0 && freq < 1.0 {
            let term = (genotypes[a][c] - 2.0 * freq) * (genotypes[b][c] - 2.0 * freq) / (2.0 * freq * (1.0 - freq))
            sum = sum + term
            count = count + 1
          }
        }
        c = c + 1
      }
      let val = if count > 0 { sum / count.to_double() } else { 0.0 }
      kinship[a].set(b, val)
      kinship[b].set(a, val)
      b = b + 1
    }
    a = a + 1
  }
  
  GSKinshipResult::new(kinship, sample_ids, n)
}

///|
/// Perform PCA on genotype data.
pub fn gs_pca(genotypes : Array[Array[Double]], n_components : Int) -> GSPcaResult {
  let n = genotypes.length()
  let m = if n > 0 { genotypes[0].length() } else { 0 }
  
  let sample_ids : Array[String] = Array::new()
  let mut i = 0
  while i < n {
    sample_ids.push("Sample" + i.to_string())
    i = i + 1
  }
  
  let allele_freqs : Array[Double] = Array::new()
  let mut j = 0
  while j < m {
    let mut sum = 0.0
    let mut n_valid = 0
    let mut k = 0
    while k < n {
      if genotypes[k][j] >= 0.0 {
        sum = sum + genotypes[k][j]
        n_valid = n_valid + 1
      }
      k = k + 1
    }
    allele_freqs.push(if n_valid > 0 { sum / (2.0 * n_valid.to_double()) } else { 0.5 })
    j = j + 1
  }
  
  let standardized : Array[Array[Double]] = Array::new()
  let mut l = 0
  while l < n {
    let row : Array[Double] = Array::new()
    let mut p = 0
    while p < m {
      if genotypes[l][p] >= 0.0 {
        let freq = allele_freqs[p]
        if freq > 0.0 && freq < 1.0 {
          row.push((genotypes[l][p] - 2.0 * freq) / (2.0 * freq * (1.0 - freq)).sqrt())
        } else {
          row.push(0.0)
        }
      } else {
        row.push(0.0)
      }
      p = p + 1
    }
    standardized.push(row)
    l = l + 1
  }
  
  let covariance : Array[Array[Double]] = Array::new()
  let mut q = 0
  while q < n {
    let row : Array[Double] = Array::new()
    let mut r = 0
    while r < n {
      let mut sum = 0.0
      let mut s = 0
      while s < m {
        sum = sum + standardized[q][s] * standardized[r][s]
        s = s + 1
      }
      row.push(sum / m.to_double())
      r = r + 1
    }
    covariance.push(row)
    q = q + 1
  }
  
  let components = gs_power_iteration(covariance, n_components)
  
  let eigenvalues : Array[Double] = Array::new()
  let mut t = 0
  while t < components.length() {
    let mut sum = 0.0
    let mut u = 0
    while u < n {
      sum = sum + components[t][u] * components[t][u]
      u = u + 1
    }
    eigenvalues.push(sum)
    t = t + 1
  }
  
  let mut total_var = 0.0
  let mut idx = 0
  while idx < eigenvalues.length() {
    total_var = total_var + eigenvalues[idx]
    idx = idx + 1
  }
  
  let var_explained : Array[Double] = Array::new()
  let mut v = 0
  while v < eigenvalues.length() {
    var_explained.push(if total_var > 0.0 { eigenvalues[v] / total_var } else { 0.0 })
    v = v + 1
  }
  
  GSPcaResult::new(eigenvalues, components, var_explained, sample_ids)
}

///|
/// Power iteration for eigenvalue decomposition.
pub fn gs_power_iteration(matrix : Array[Array[Double]], n_components : Int) -> Array[Array[Double]] {
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
/// Calculate genetic distance matrix.
pub fn gs_calculate_genetic_distance(
  genotypes : Array[Array[Double]],
  distance_type? : String = "euclidean"
) -> GSDistanceResult {
  let n = genotypes.length()
  
  let dt = distance_type
  
  let sample_ids : Array[String] = Array::new()
  let mut i = 0
  while i < n {
    sample_ids.push("Sample" + i.to_string())
    i = i + 1
  }
  
  let distance : Array[Array[Double]] = Array::new()
  let mut j = 0
  while j < n {
    let row : Array[Double] = Array::new()
    let mut k = 0
    while k < n {
      row.push(0.0)
      k = k + 1
    }
    distance.push(row)
    j = j + 1
  }
  
  let mut a = 0
  while a < n {
    let mut b = a + 1
    while b < n {
      let dist = if dt == "euclidean" {
        gs_euclidean_distance(genotypes[a], genotypes[b])
      } else if dt == "manhattan" {
        gs_manhattan_distance(genotypes[a], genotypes[b])
      } else if dt == "identity_by_state" {
        gs_ibs_distance(genotypes[a], genotypes[b])
      } else {
        gs_euclidean_distance(genotypes[a], genotypes[b])
      }
      
      distance[a].set(b, dist)
      distance[b].set(a, dist)
      b = b + 1
    }
    a = a + 1
  }
  
  GSDistanceResult::new(distance, sample_ids, dt)
}

///|
/// Calculate Euclidean distance.
pub fn gs_euclidean_distance(v1 : Array[Double], v2 : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < v1.length() && i < v2.length() {
    if v1[i] >= 0.0 && v2[i] >= 0.0 {
      let diff = v1[i] - v2[i]
      sum = sum + diff * diff
    }
    i = i + 1
  }
  sum.sqrt()
}

///|
/// Calculate Manhattan distance.
pub fn gs_manhattan_distance(v1 : Array[Double], v2 : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < v1.length() && i < v2.length() {
    if v1[i] >= 0.0 && v2[i] >= 0.0 {
      sum = sum + (v1[i] - v2[i]).abs()
    }
    i = i + 1
  }
  sum
}

///|
/// Calculate Identity-by-State distance.
pub fn gs_ibs_distance(v1 : Array[Double], v2 : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut count = 0
  let mut i = 0
  while i < v1.length() && i < v2.length() {
    if v1[i] >= 0.0 && v2[i] >= 0.0 {
      let diff = (v1[i] - v2[i]).abs()
      sum = sum + diff
      count = count + 1
    }
    i = i + 1
  }
  if count > 0 { sum / count.to_double() } else { 0.0 }
}

///|
/// Bio API functions.

pub fn bio_estimate_kinship(genotypes : Array[Array[Double]]) -> GSKinshipResult {
  gs_estimate_kinship(genotypes)
}

pub fn bio_pca(genotypes : Array[Array[Double]], n_components : Int) -> GSPcaResult {
  gs_pca(genotypes, n_components)
}

pub fn bio_genetic_distance(genotypes : Array[Array[Double]], distance_type? : String) -> GSDistanceResult {
  gs_calculate_genetic_distance(genotypes)
}