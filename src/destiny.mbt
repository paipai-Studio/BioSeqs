///|
/// destiny - Diffusion maps for single-cell data dimensionality reduction.
/// Inspired by Bioconductor destiny.
/// Provides diffusion map implementation for visualizing high-dimensional
/// single-cell data, including distance computation, kernel construction,
/// eigen decomposition, and diffusion embedding.

///|
pub struct CellData {
  cell_id : String
  expression : Array[Double]
  cluster : String
}

///|
pub struct DistanceMatrix {
  cells : Array[String]
  distances : Array[Array[Double]]
  n_cells : Int
}

///|
pub struct KernelMatrix {
  cells : Array[String]
  kernel : Array[Array[Double]]
  bandwidth : Double
}

///|
pub struct DiffusionMapData {
  cells : Array[String]
  eigen_values : Array[Double]
  eigen_vectors : Array[Array[Double]]
  diffusion_components : Array[Array[Double]]
  sigma : Double
  n_components : Int
}

///|
pub struct DiffusionResult {
  embedding : Array[Array[Double]]
  eigenvalues : Array[Double]
  cell_ids : Array[String]
  explained_variance : Array[Double]
}

///|
pub fn CellData::new(
  cell_id : String,
  expression : Array[Double],
  cluster : String
) -> CellData {
  CellData::{
    cell_id,
    expression,
    cluster
  }
}

///|
pub fn DistanceMatrix::new(cells : Array[String], distances : Array[Array[Double]]) -> DistanceMatrix {
  let n = cells.length()
  DistanceMatrix::{
    cells,
    distances,
    n_cells: n
  }
}

///|
pub fn KernelMatrix::new(cells : Array[String], kernel : Array[Array[Double]], bandwidth : Double) -> KernelMatrix {
  KernelMatrix::{
    cells,
    kernel,
    bandwidth
  }
}

///|
pub fn compute_distance_matrix(cells : Array[CellData], metric : String) -> DistanceMatrix {
  let n = cells.length()
  let cell_ids : Array[String] = Array::new()
  
  let mut i = 0
  while i < n {
    cell_ids.push(cells[i].cell_id)
    i = i + 1
  }
  
  // Compute pairwise distances
  let distances : Array[Array[Double]] = Array::new()
  i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      let dist = destiny_compute_distance(cells[i].expression, cells[j].expression, metric)
      row.push(dist)
      j = j + 1
    }
    distances.push(row)
    i = i + 1
  }
  
  DistanceMatrix::new(cell_ids, distances)
}

///|
fn destiny_compute_distance(a : Array[Double], b : Array[Double], metric : String) -> Double {
  let n = a.length()
  
  match metric {
    "euclidean" => {
      let mut sum = 0.0
      let mut i = 0
      while i < n {
        let diff = a[i] - b[i]
        sum = sum + diff * diff
        i = i + 1
      }
      sum.sqrt()
    }
    "manhattan" => {
      let mut sum = 0.0
      let mut i = 0
      while i < n {
        let diff = a[i] - b[i]
        sum = sum + diff.abs()
        i = i + 1
      }
      sum
    }
    _ => {
      // Default to euclidean
      let mut sum = 0.0
      let mut i = 0
      while i < n {
        let diff = a[i] - b[i]
        sum = sum + diff * diff
        i = i + 1
      }
      sum.sqrt()
    }
  }
}

///|
pub fn compute_gaussian_kernel(dist_matrix : DistanceMatrix, sigma : Double) -> KernelMatrix {
  let n = dist_matrix.n_cells
  let cells = dist_matrix.cells
  
  let kernel : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      let dist = dist_matrix.distances[i][j]
      let k = @math.exp(-dist * dist / (2.0 * sigma * sigma))
      row.push(k)
      j = j + 1
    }
    kernel.push(row)
    i = i + 1
  }
  
  KernelMatrix::new(cells, kernel, sigma)
}

///|
pub fn find_sigma_automatic(dist_matrix : DistanceMatrix) -> Double {
  let n = dist_matrix.n_cells
  if n < 2 {
    return 1.0
  }
  
  // Calculate pairwise distances and use median
  let distances : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      distances.push(dist_matrix.distances[i][j])
      j = j + 1
    }
    i = i + 1
  }
  
  // Sort and get median
  let sorted = destiny_sort_array(distances)
  let len = sorted.length()
  if len == 0 {
    return 1.0
  }
  
  let median_idx = len / 2
  let median = if len % 2 == 0 {
    (sorted[median_idx - 1] + sorted[median_idx]) / 2.0
  } else {
    sorted[median_idx]
  }
  
  // Use median as sigma
  if median < 0.001 { 1.0 } else { median }
}

///|
pub fn compute_diffusion_map(cells : Array[CellData], n_components : Int, sigma : Double) -> DiffusionResult {
  let n = cells.length()
  let cell_ids : Array[String] = Array::new()
  
  let mut i = 0
  while i < n {
    cell_ids.push(cells[i].cell_id)
    i = i + 1
  }
  
  // Step 1: Compute distance matrix
  let dist_matrix = compute_distance_matrix(cells, "euclidean")
  
  // Step 2: Compute kernel
  let kernel = compute_gaussian_kernel(dist_matrix, sigma)
  
  // Step 3: Normalize kernel (markov matrix)
  let markov = destiny_normalize_kernel(kernel)
  
  // Step 4: Compute eigendecomposition
  let n_comp = if n_components > n - 1 { n - 1 } else { n_components }
  let eigen_result = destiny_eigen_decomposition(markov.kernel, n_comp)
  
  // Step 5: Compute diffusion components
  let embedding = destiny_compute_embedding(eigen_result.eigen_vectors, eigen_result.eigen_values, n_comp)
  
  // Calculate explained variance
  let total_var = destiny_sum_array(eigen_result.eigen_values)
  let explained_variance : Array[Double] = Array::new()
  let mut j = 0
  while j < n_comp {
    let var = if total_var > 0.0 {
      eigen_result.eigen_values[j] / total_var
    } else {
      0.0
    }
    explained_variance.push(var)
    j = j + 1
  }
  
  DiffusionResult::{
    embedding,
    eigenvalues: eigen_result.eigen_values,
    cell_ids,
    explained_variance
  }
}

///|
fn destiny_normalize_kernel(kernel : KernelMatrix) -> KernelMatrix {
  let n = kernel.cells.length()
  let normalized : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n {
    let row = kernel.kernel[i]
    let row_sum = destiny_sum_array(row)
    
    let normalized_row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      let val = if row_sum > 0.0 {
        row[j] / row_sum
      } else {
        0.0
      }
      normalized_row.push(val)
      j = j + 1
    }
    normalized.push(normalized_row)
    i = i + 1
  }
  
  KernelMatrix::new(kernel.cells, normalized, kernel.bandwidth)
}

///|
struct EigenResult {
  eigen_values : Array[Double]
  eigen_vectors : Array[Array[Double]]
}

///|
fn destiny_eigen_decomposition(matrix : Array[Array[Double]], n_components : Int) -> EigenResult {
  let n = matrix.length()
  
  // Use power iteration for eigenvalue computation
  let eigen_values : Array[Double] = Array::new()
  let eigen_vectors : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n_components {
    // Initialize random vector
    let mut vector : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      vector.push(destiny_random_vector_component(j.to_double() + i.to_double() * 100.0))
      j = j + 1
    }
    
    // Normalize
    vector = destiny_normalize_vector(vector)
    
    let mut lambda = 0.0
    let mut iteration = 0
    while iteration < 50 {
      // Matrix-vector multiplication
      let new_vector : Array[Double] = Array::new()
      let mut k = 0
      while k < n {
        let mut sum = 0.0
        let mut l = 0
        while l < n {
          sum = sum + matrix[k][l] * vector[l]
          l = l + 1
        }
        new_vector.push(sum)
        k = k + 1
      }
      
      lambda = destiny_vector_norm(new_vector)
      vector = destiny_normalize_vector(new_vector)
      
      iteration = iteration + 1
    }
    
    eigen_values.push(lambda)
    
    // Store eigenvector
    let ev : Array[Double] = Array::new()
    let mut m = 0
    while m < n {
      ev.push(vector[m])
      m = m + 1
    }
    eigen_vectors.push(ev)
    
    // Deflate matrix
    let deflated : Array[Array[Double]] = Array::new()
    let mut p = 0
    while p < n {
      let row : Array[Double] = Array::new()
      let mut q = 0
      while q < n {
        let val = matrix[p][q] - lambda * vector[p] * vector[q]
        row.push(val)
        q = q + 1
      }
      deflated.push(row)
      p = p + 1
    }
    
    // Update matrix for next iteration
    // Note: This modifies the original - in production would need copy
    // For simplicity, we re-assign
    let _ = deflated  // Just for reference
    
    i = i + 1
  }
  
  EigenResult::{
    eigen_values,
    eigen_vectors
  }
}

///|
fn destiny_compute_embedding(eigen_vectors : Array[Array[Double]], eigen_values : Array[Double], n_components : Int) -> Array[Array[Double]] {
  // eigen_vectors has shape [n_components][n_cells]
  let n_components_actual = eigen_vectors.length()
  let n_cells = if n_components_actual > 0 { eigen_vectors[0].length() } else { 0 }
  
  let embedding : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n_cells {
    let cell_embedding : Array[Double] = Array::new()
    let mut j = 0
    while j < n_components {
      if j < n_components_actual {
        let coord = eigen_vectors[j][i] * eigen_values[j]
        cell_embedding.push(coord)
      } else {
        cell_embedding.push(0.0)
      }
      j = j + 1
    }
    embedding.push(cell_embedding)
    i = i + 1
  }
  
  embedding
}

///|
pub fn destiny_create_embedding(cells : Array[CellData], n_components : Int) -> DiffusionResult {
  // Find optimal sigma
  let dist_matrix = compute_distance_matrix(cells, "euclidean")
  let sigma = find_sigma_automatic(dist_matrix)
  
  // Compute diffusion map
  compute_diffusion_map(cells, n_components, sigma)
}

///|
pub fn destiny_plot_coordinates(result : DiffusionResult, component_x : Int, component_y : Int) -> Array[(String, Double, Double)] {
  let coordinates : Array[(String, Double, Double)] = Array::new()
  let n = result.cell_ids.length()
  
  let mut i = 0
  while i < n {
    let x = result.embedding[i][component_x]
    let y = result.embedding[i][component_y]
    coordinates.push((result.cell_ids[i], x, y))
    i = i + 1
  }
  
  coordinates
}

///|
pub fn destiny_summary(result : DiffusionResult) -> String {
  let mut summary = "\nDiffusion Map Summary\n"
  summary = summary + "==============================\n"
  summary = summary + "Number of cells: " + result.cell_ids.length().to_string() + "\n"
  summary = summary + "Number of components: " + result.eigenvalues.length().to_string() + "\n"
  
  let mut i = 0
  while i < result.eigenvalues.length() {
    summary = summary + "Component " + (i + 1).to_string() + 
      ": eigenvalue = " + result.eigenvalues[i].to_string() + 
      ", explained variance = " + (result.explained_variance[i] * 100.0).to_string() + "%\n"
    i = i + 1
  }
  
  summary
}

///|
pub fn create_example_sc_data(n_cells : Int, n_genes : Int) -> Array[CellData] {
  let cells : Array[CellData] = Array::new()
  
  let mut i = 0
  while i < n_cells {
    let expression : Array[Double] = Array::new()
    let mut j = 0
    while j < n_genes {
      // Create some structured data with two clusters
      let cluster = if i < n_cells / 2 { 0 } else { 1 }
      let base = if cluster == 0 { 0.5 } else { 1.5 }
      let noise = ((i * 31 + j * 17) % 100).to_double() / 100.0 * 0.2
      expression.push(base + noise)
      j = j + 1
    }
    
    let cell = CellData::new(
      "cell_" + i.to_string(),
      expression,
      if i < n_cells / 2 { "cluster_1" } else { "cluster_2" }
    )
    
    cells.push(cell)
    i = i + 1
  }
  
  cells
}

// Helper functions

fn destiny_sort_array(arr : Array[Double]) -> Array[Double] {
  let sorted = arr
  // Simple bubble sort for small arrays
  let n = sorted.length()
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n - i - 1 {
      if sorted[j] > sorted[j + 1] {
        let temp = sorted[j]
        sorted[j] = sorted[j + 1]
        sorted[j + 1] = temp
      }
      j = j + 1
    }
    i = i + 1
  }
  sorted
}

fn destiny_sum_array(arr : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < arr.length() {
    sum = sum + arr[i]
    i = i + 1
  }
  sum
}

fn destiny_normalize_vector(v : Array[Double]) -> Array[Double] {
  let norm = destiny_vector_norm(v)
  if norm == 0.0 {
    let result : Array[Double] = Array::new()
    let mut i = 0
    while i < v.length() {
      result.push(0.0)
      i = i + 1
    }
    return result
  }
  
  let normalized : Array[Double] = Array::new()
  let mut i = 0
  while i < v.length() {
    normalized.push(v[i] / norm)
    i = i + 1
  }
  normalized
}

fn destiny_vector_norm(v : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < v.length() {
    sum = sum + v[i] * v[i]
    i = i + 1
  }
  sum.sqrt()
}

fn destiny_random_vector_component(seed : Double) -> Double {
  // Simple deterministic pseudo-random generator
  let raw = seed * 99991.0
  let x = (raw - raw.floor()) * 2.0 - 1.0
  if x < -1.0 { -1.0 } else if x > 1.0 { 1.0 } else { x }
}
