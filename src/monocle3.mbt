///|
/// monocle3 module, inspired by Bioconductor monocle3.
/// Provides functions for single-cell trajectory analysis.

///|
pub struct CellDataSet {
  counts : Array[Array[Int]]
  normalized : Array[Array[Double]]
  gene_names : Array[String]
  cell_names : Array[String]
  reduced_dimensions : Array[(String, Array[Array[Double]])]
  pseudotime : Array[Double]
  cell_partitions : Array[Int]
  principal_graph : Array[(Int, Array[Double], Array[Int])]
}

///|
pub struct TrajectoryResult {
  pseudotime : Array[Double]
  partitions : Array[Int]
  num_partitions : Int
  principal_graph_nodes : Array[Array[Double]]
  principal_graph_edges : Array[(Int, Int)]
}

///|
pub struct DiffExprResult {
  gene : Array[String]
  pseudotime : Array[Double]
  q_value : Array[Double]
  trend_statistic : Array[Double]
}

///|
pub fn new_cell_data_set(
  counts : Array[Array[Int]],
  gene_names : Array[String],
  cell_names : Array[String],
) -> CellDataSet {
  CellDataSet::{
    counts,
    normalized: Array::new(),
    gene_names,
    cell_names,
    reduced_dimensions: Array::new(),
    pseudotime: Array::new(),
    cell_partitions: Array::new(),
    principal_graph: Array::new(),
  }
}

///|
pub fn preprocess_cds(cds : CellDataSet, num_dimensions : Int) -> CellDataSet {
  let n_cells = cds.counts.length()
  let n_genes = if n_cells > 0 { cds.counts[0].length() } else { 0 }

  if n_cells == 0 || n_genes == 0 {
    return cds
  }

  let normalized : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_cells {
    let row : Array[Double] = Array::new()
    let mut total = 0.0
    let mut j = 0
    while j < n_genes {
      total = total + cds.counts[i][j].to_double()
      j = j + 1
    }
    j = 0
    while j < n_genes {
      let norm = if total > 0.0 {
        cds.counts[i][j].to_double() / total * 10000.0
      } else {
        0.0
      }
      let log_val = if norm > 0.0 { @math.ln(norm + 1.0) } else { 0.0 }
      row.push(log_val)
      j = j + 1
    }
    normalized.push(row)
    i = i + 1
  }

  let pca_result = pca(normalized, num_dimensions)

  let reduced_dims : Array[(String, Array[Array[Double]])] = Array::new()
  reduced_dims.push(("PCA".to_string(), pca_result.0))

  CellDataSet::{
    counts: cds.counts,
    normalized,
    gene_names: cds.gene_names,
    cell_names: cds.cell_names,
    reduced_dimensions: reduced_dims,
    pseudotime: cds.pseudotime,
    cell_partitions: cds.cell_partitions,
    principal_graph: cds.principal_graph,
  }
}

///|
fn pca(
  data : Array[Array[Double]],
  n_components : Int,
) -> (Array[Array[Double]], Array[Double]) {
  let n_samples = data.length()
  let n_features = if n_samples > 0 { data[0].length() } else { 0 }

  if n_samples == 0 || n_features == 0 || n_components <= 0 {
    return (Array::new(), Array::new())
  }

  let centered = center_data(data)
  let cov_matrix = covariance_matrix(centered)

  let (eigenvalues, eigenvectors) = eigen_decomposition(
    cov_matrix, n_components,
  )

  let transformed : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_samples {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n_components && j < n_features {
      let mut sum = 0.0
      let mut k = 0
      while k < n_features {
        sum = sum + centered[i][k] * eigenvectors[j][k]
        k = k + 1
      }
      row.push(sum)
      j = j + 1
    }
    transformed.push(row)
    i = i + 1
  }

  (transformed, eigenvalues)
}

///|
fn center_data(data : Array[Array[Double]]) -> Array[Array[Double]] {
  let n_samples = data.length()
  let n_features = if n_samples > 0 { data[0].length() } else { 0 }

  if n_samples == 0 || n_features == 0 {
    return Array::new()
  }

  let means : Array[Double] = Array::make(n_features, 0.0)
  let mut j = 0
  while j < n_features {
    let mut sum = 0.0
    let mut i = 0
    while i < n_samples {
      sum = sum + data[i][j]
      i = i + 1
    }
    means[j] = sum / n_samples.to_double()
    j = j + 1
  }

  let centered : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_samples {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n_features {
      row.push(data[i][j] - means[j])
      j = j + 1
    }
    centered.push(row)
    i = i + 1
  }

  centered
}

///|
fn covariance_matrix(data : Array[Array[Double]]) -> Array[Array[Double]] {
  let n_samples = data.length()
  let n_features = if n_samples > 0 { data[0].length() } else { 0 }

  if n_samples == 0 || n_features == 0 {
    return Array::new()
  }

  let cov : Array[Array[Double]] = Array::make(
    n_features,
    Array::make(n_features, 0.0),
  )

  let mut i = 0
  while i < n_features {
    let mut j = 0
    while j < n_features {
      let mut sum = 0.0
      let mut k = 0
      while k < n_samples {
        sum = sum + data[k][i] * data[k][j]
        k = k + 1
      }
      cov[i][j] = sum / (n_samples - 1).to_double()
      j = j + 1
    }
    i = i + 1
  }

  cov
}

///|
fn eigen_decomposition(
  matrix : Array[Array[Double]],
  n_components : Int,
) -> (Array[Double], Array[Array[Double]]) {
  let n = matrix.length()

  if n == 0 || n_components <= 0 {
    return (Array::new(), Array::new())
  }

  let eigenvalues : Array[Double] = Array::make(n_components, 0.0)
  let eigenvectors : Array[Array[Double]] = Array::make(
    n_components,
    Array::make(n, 0.0),
  )

  let mut v : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  let mut i = 0
  while i < n {
    v[i][i] = 1.0
    i = i + 1
  }

  let mut a : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  i = 0
  while i < n {
    let mut j = 0
    while j < n {
      a[i][j] = matrix[i][j]
      j = j + 1
    }
    i = i + 1
  }

  let mut iter = 0
  while iter < 50 {
    let mut max_val = 0.0
    let mut max_i = 0
    let mut max_j = 1

    let mut i = 0
    while i < n {
      let mut j = 0
      while j < n {
        if i != j && a[i][j].abs() > max_val {
          max_val = a[i][j].abs()
          max_i = i
          max_j = j
        }
        j = j + 1
      }
      i = i + 1
    }

    if max_val < 0.0000000001 {
      break
    }

    let theta = 0.5 *
      @math.atan(2.0 * a[max_i][max_j] / (a[max_i][max_i] - a[max_j][max_j]))
    let cos_theta = @math.cos(theta)
    let sin_theta = @math.sin(theta)

    let p : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
    i = 0
    while i < n {
      p[i][i] = 1.0
      i = i + 1
    }
    p[max_i][max_i] = cos_theta
    p[max_j][max_j] = cos_theta
    p[max_i][max_j] = -sin_theta
    p[max_j][max_i] = sin_theta

    let pt = monocle_transpose(p)
    let pap = monocle_matmul(monocle_matmul(pt, a), p)
    let vp = monocle_matmul(v, p)

    a = pap
    v = vp

    iter = iter + 1
  }

  let diag : Array[(Double, Int)] = Array::new()
  i = 0
  while i < n {
    diag.push((a[i][i], i))
    i = i + 1
  }

  let mut j = 0
  while j < diag.length() - 1 {
    let mut k = j + 1
    while k < diag.length() {
      if diag[j].0 < diag[k].0 {
        let tmp = diag[j]
        diag[j] = diag[k]
        diag[k] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }

  i = 0
  while i < n_components && i < n {
    eigenvalues[i] = diag[i].0
    let idx = diag[i].1
    let mut j = 0
    while j < n {
      eigenvectors[i][j] = v[j][idx]
      j = j + 1
    }
    i = i + 1
  }

  (eigenvalues, eigenvectors)
}

///|
fn monocle_matmul(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let m = a.length()
  let n = if m > 0 { a[0].length() } else { 0 }
  let p = if b.length() > 0 { b[0].length() } else { 0 }

  if m == 0 || n == 0 || p == 0 {
    return Array::new()
  }

  let result : Array[Array[Double]] = Array::make(m, Array::make(p, 0.0))

  let mut i = 0
  while i < m {
    let mut j = 0
    while j < p {
      let mut sum = 0.0
      let mut k = 0
      while k < n {
        sum = sum + a[i][k] * b[k][j]
        k = k + 1
      }
      result[i][j] = sum
      j = j + 1
    }
    i = i + 1
  }

  result
}

///|
fn monocle_transpose(m : Array[Array[Double]]) -> Array[Array[Double]] {
  let rows = m.length()
  let cols = if rows > 0 { m[0].length() } else { 0 }

  if rows == 0 || cols == 0 {
    return Array::new()
  }

  let result : Array[Array[Double]] = Array::make(cols, Array::make(rows, 0.0))

  let mut i = 0
  while i < rows {
    let mut j = 0
    while j < cols {
      result[j][i] = m[i][j]
      j = j + 1
    }
    i = i + 1
  }

  result
}

///|
pub fn reduce_dimension(
  cds : CellDataSet,
  dim_method : String,
  num_dimensions : Int,
) -> CellDataSet {
  if dim_method == "PCA" || dim_method == "pca" {
    return preprocess_cds(cds, num_dimensions)
  }

  if dim_method == "UMAP" || dim_method == "umap" {
    let pca_data = match get_reduced_dimension(cds, "PCA") {
      Some(data) => data
      None => {
        let preprocessed = preprocess_cds(cds, 50)
        match get_reduced_dimension(preprocessed, "PCA") {
          Some(d) => d
          None => cds.normalized
        }
      }
    }

    let umap_result = umap(pca_data, num_dimensions)

    let reduced_dims = merge_reduced_dimensions(cds.reduced_dimensions, "UMAP", umap_result)

    return CellDataSet::{
      counts: cds.counts,
      normalized: cds.normalized,
      gene_names: cds.gene_names,
      cell_names: cds.cell_names,
      reduced_dimensions: reduced_dims,
      pseudotime: cds.pseudotime,
      cell_partitions: cds.cell_partitions,
      principal_graph: cds.principal_graph,
    }
  }

  cds
}

///|
fn merge_reduced_dimensions(
  existing : Array[(String, Array[Array[Double]])],
  name : String,
  data : Array[Array[Double]],
) -> Array[(String, Array[Array[Double]])] {
  let n = existing.length()
  let mut has_name = false
  let mut i = 0
  while i < n {
    if existing[i].0 == name {
      has_name = true
    }
    i = i + 1
  }
  
  let result_len = if has_name { n } else { n + 1 }
  let result : Array[(String, Array[Array[Double]])] = Array::make(result_len, ("".to_string(), Array::new()))
  
  let mut j = 0
  i = 0
  while i < n {
    if existing[i].0 == name {
      result[j] = (name.to_string(), data)
    } else {
      result[j] = existing[i]
    }
    j = j + 1
    i = i + 1
  }
  
  if !has_name {
    result[j] = (name.to_string(), data)
  }
  
  result
}

///|
fn get_reduced_dimension(
  cds : CellDataSet,
  name : String,
) -> Array[Array[Double]]? {
  let mut i = 0
  while i < cds.reduced_dimensions.length() {
    if cds.reduced_dimensions[i].0 == name {
      return Some(cds.reduced_dimensions[i].1)
    }
    i = i + 1
  }
  None
}

///|
fn umap(
  data : Array[Array[Double]],
  n_components : Int,
) -> Array[Array[Double]] {
  let n_samples = data.length()
  let n_features = if n_samples > 0 { data[0].length() } else { 0 }

  if n_samples == 0 || n_features == 0 || n_components <= 0 {
    return Array::new()
  }

  let embedding : Array[Array[Double]] = Array::make(
    n_samples,
    Array::make(n_components, 0.0),
  )
  let mut i = 0
  while i < n_samples {
    let mut j = 0
    while j < n_components {
      embedding[i][j] = ((i * 7 + j * 13).to_double() % 10.0 / 10.0 - 0.5) * 0.1
      j = j + 1
    }
    i = i + 1
  }

  let k = if n_samples > 15 { 15 } else { n_samples - 1 }

  let knn_indices : Array[Array[Int]] = Array::new()
  let knn_distances : Array[Array[Double]] = Array::new()

  i = 0
  while i < n_samples {
    knn_indices.push(Array::new())
    knn_distances.push(Array::new())
    i = i + 1
  }

  i = 0
  while i < n_samples {
    let distances : Array[(Double, Int)] = Array::new()
    let mut j = 0
    while j < n_samples {
      let mut sum_sq = 0.0
      let mut k_idx = 0
      while k_idx < n_features {
        let diff = data[i][k_idx] - data[j][k_idx]
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

    let mut j = 1
    while j <= k && j < distances.length() {
      knn_indices[i].push(distances[j].1)
      knn_distances[i].push(distances[j].0)
      j = j + 1
    }
    i = i + 1
  }

  let graph : Array[Array[(Int, Double)]] = Array::new()
  i = 0
  while i < n_samples {
    graph.push(Array::new())
    i = i + 1
  }

  i = 0
  while i < n_samples {
    let mut j = 0
    while j < knn_indices[i].length() {
      let neighbor = knn_indices[i][j]
      let distance = knn_distances[i][j]
      graph[i].push((neighbor, distance))

      let mut found = false
      let mut k_idx = 0
      while k_idx < graph[neighbor].length() && !found {
        if graph[neighbor][k_idx].0 == i {
          found = true
        }
        k_idx = k_idx + 1
      }
      if !found {
        graph[neighbor].push((i, distance))
      }
      j = j + 1
    }
    i = i + 1
  }

  let alpha = 1.0
  let gamma = 1.0
  let mut learning_rate = 0.5
  let n_epochs = 100

  let mut epoch = 0
  while epoch < n_epochs {
    let mut i = 0
    while i < n_samples {
      let gradients : Array[Double] = Array::make(n_components, 0.0)

      let mut j = 0
      while j < graph[i].length() {
        let neighbor = graph[i][j].0

        let mut dist_emb = 0.0
        let mut k_idx = 0
        while k_idx < n_components {
          let diff = embedding[i][k_idx] - embedding[neighbor][k_idx]
          dist_emb = dist_emb + diff * diff
          k_idx = k_idx + 1
        }
        dist_emb = dist_emb.sqrt()

        let grad_coeff = (1.0 + alpha * dist_emb * dist_emb) /
          (alpha * dist_emb + 1.0 / (dist_emb + 0.0000000001))

        let mut k_idx = 0
        while k_idx < n_components {
          let diff = embedding[i][k_idx] - embedding[neighbor][k_idx]
          gradients[k_idx] = gradients[k_idx] +
            grad_coeff * diff / (dist_emb + 0.0000000001)
          k_idx = k_idx + 1
        }
        j = j + 1
      }

      let n_neg = 3
      let mut neg = 0
      while neg < n_neg {
        let rand_idx = (i * 17 + neg * 13 + epoch * 7) % n_samples

        let mut dist_emb = 0.0
        let mut k_idx = 0
        while k_idx < n_components {
          let diff = embedding[i][k_idx] - embedding[rand_idx][k_idx]
          dist_emb = dist_emb + diff * diff
          k_idx = k_idx + 1
        }
        dist_emb = dist_emb.sqrt()

        if dist_emb > 0.1 {
          let grad_coeff = gamma / (dist_emb * dist_emb + 1.0)

          let mut k_idx = 0
          while k_idx < n_components {
            let diff = embedding[i][k_idx] - embedding[rand_idx][k_idx]
            gradients[k_idx] = gradients[k_idx] -
              grad_coeff * diff / (dist_emb + 0.0000000001)
            k_idx = k_idx + 1
          }
          neg = neg + 1
        } else {
          neg = neg + 1
        }
      }

      let mut k_idx = 0
      while k_idx < n_components {
        embedding[i][k_idx] = embedding[i][k_idx] -
          learning_rate * gradients[k_idx]
        k_idx = k_idx + 1
      }

      i = i + 1
    }

    learning_rate = learning_rate * 0.99
    epoch = epoch + 1
  }

  embedding
}

///|
pub fn learn_graph(cds : CellDataSet) -> CellDataSet {
  let umap_data = match get_reduced_dimension(cds, "UMAP") {
    Some(data) => data
    None => {
      let reduced = reduce_dimension(cds, "UMAP", 2)
      match get_reduced_dimension(reduced, "UMAP") {
        Some(d) => d
        None => {
          if reduced.normalized.length() > 0 {
            reduced.normalized
          } else {
            return cds
          }
        }
      }
    }
  }

  let n_cells = umap_data.length()

  if n_cells == 0 {
    return cds
  }

  let partitions : Array[Int] = Array::make(n_cells, 0)
  let num_partitions = 1

  let centers : Array[Array[Double]] = Array::new()
  let center : Array[Double] = Array::make(2, 0.0)
  let mut i = 0
  while i < n_cells {
    center[0] = center[0] + umap_data[i][0]
    center[1] = center[1] + umap_data[i][1]
    i = i + 1
  }
  center[0] = center[0] / n_cells.to_double()
  center[1] = center[1] / n_cells.to_double()
  centers.push(center)

  let graph_edges : Array[(Int, Int)] = Array::new()

  let mut j = 0
  while j < n_cells {
    let mut min_dist = 1000000000000000000000000000000.0
    let mut min_idx = 0

    let mut k = 0
    while k < centers.length() {
      let dx = umap_data[j][0] - centers[k][0]
      let dy = umap_data[j][1] - centers[k][1]
      let dist = dx * dx + dy * dy
      if dist < min_dist {
        min_dist = dist
        min_idx = k
      }
      k = k + 1
    }

    partitions[j] = min_idx
    j = j + 1
  }

  let nodes : Array[Array[Double]] = Array::new()
  let mut node_idx = 0
  while node_idx < num_partitions {
    let node : Array[Double] = Array::make(2, 0.0)
    let mut count = 0
    let mut i = 0
    while i < n_cells {
      if partitions[i] == node_idx {
        node[0] = node[0] + umap_data[i][0]
        node[1] = node[1] + umap_data[i][1]
        count = count + 1
      }
      i = i + 1
    }
    if count > 0 {
      node[0] = node[0] / count.to_double()
      node[1] = node[1] / count.to_double()
    }
    nodes.push(node)
    node_idx = node_idx + 1
  }

  if nodes.length() > 1 {
    let mut i = 0
    while i < nodes.length() - 1 {
      graph_edges.push((i, i + 1))
      i = i + 1
    }
  }

  let principal_graph : Array[(Int, Array[Double], Array[Int])] = Array::new()
  i = 0
  while i < nodes.length() {
    let edges : Array[Int] = Array::new()
    let mut j = 0
    while j < graph_edges.length() {
      if graph_edges[j].0 == i || graph_edges[j].1 == i {
        edges.push(
          if graph_edges[j].0 == i {
            graph_edges[j].1
          } else {
            graph_edges[j].0
          },
        )
      }
      j = j + 1
    }
    principal_graph.push((i, nodes[i], edges))
    i = i + 1
  }

  CellDataSet::{
    counts: cds.counts,
    normalized: cds.normalized,
    gene_names: cds.gene_names,
    cell_names: cds.cell_names,
    reduced_dimensions: cds.reduced_dimensions,
    pseudotime: Array::new(),
    cell_partitions: partitions,
    principal_graph,
  }
}

///|
pub fn order_cells(cds : CellDataSet, root_cells : Array[Int]) -> CellDataSet {
  let umap_data = match get_reduced_dimension(cds, "UMAP") {
    Some(data) => data
    None => cds.normalized
  }

  let n_cells = umap_data.length()

  if n_cells == 0 {
    return cds
  }

  let pseudotime : Array[Double] = Array::make(n_cells, -1.0)

  let distances : Array[(Double, Int)] = Array::new()
  let mut i = 0
  while i < n_cells {
    let mut min_dist = 1000000000000000000000000000000.0
    let mut j = 0
    while j < root_cells.length() {
      if root_cells[j] < n_cells {
        let dx = umap_data[i][0] - umap_data[root_cells[j]][0]
        let dy = umap_data[i][1] - umap_data[root_cells[j]][1]
        let dist = dx * dx + dy * dy
        if dist < min_dist {
          min_dist = dist
        }
      }
      j = j + 1
    }
    distances.push((min_dist.sqrt(), i))
    i = i + 1
  }

  let mut j = 0
  while j < distances.length() - 1 {
    let mut k = j + 1
    while k < distances.length() {
      if distances[j].0 > distances[k].0 {
        let tmp = distances[j]
        distances[j] = distances[k]
        distances[k] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }

  let mut max_dist = 0.0
  i = 0
  while i < distances.length() {
    if distances[i].0 > max_dist {
      max_dist = distances[i].0
    }
    i = i + 1
  }

  i = 0
  while i < distances.length() {
    pseudotime[distances[i].1] = if max_dist > 0.0 {
      distances[i].0 / max_dist
    } else {
      0.0
    }
    i = i + 1
  }

  CellDataSet::{
    counts: cds.counts,
    normalized: cds.normalized,
    gene_names: cds.gene_names,
    cell_names: cds.cell_names,
    reduced_dimensions: cds.reduced_dimensions,
    pseudotime,
    cell_partitions: cds.cell_partitions,
    principal_graph: cds.principal_graph,
  }
}

///|
pub fn fit_models(
  cds : CellDataSet,
  gene_names : Array[String],
) -> DiffExprResult {
  let n_cells = cds.pseudotime.length()
  let n_genes = gene_names.length()

  let gene_result : Array[String] = Array::new()
  let pseudotime_result : Array[Double] = Array::new()
  let q_value_result : Array[Double] = Array::new()
  let trend_statistic_result : Array[Double] = Array::new()

  let mut i = 0
  while i < n_genes {
    let gene_idx = find_gene_index(cds.gene_names, gene_names[i])
    if gene_idx.is_some() {
      let idx = gene_idx.unwrap()

      let mut sum_expr = 0.0
      let mut sum_pt = 0.0
      let mut sum_pt_sq = 0.0
      let mut sum_expr_pt = 0.0
      let mut n_valid = 0

      let mut j = 0
      while j < n_cells {
        if cds.pseudotime[j] >= 0.0 {
          let expr = if j < cds.normalized.length() && idx < cds.normalized[j].length() {
            cds.normalized[j][idx]
          } else {
            0.0
          }
          sum_expr = sum_expr + expr
          sum_pt = sum_pt + cds.pseudotime[j]
          sum_pt_sq = sum_pt_sq + cds.pseudotime[j] * cds.pseudotime[j]
          sum_expr_pt = sum_expr_pt + expr * cds.pseudotime[j]
          n_valid = n_valid + 1
        }
        j = j + 1
      }

      if n_valid > 0 {
        let n_valid_double = n_valid.to_double()
        let slope = if n_valid_double * sum_pt_sq - sum_pt * sum_pt > 0.0 {
          (n_valid_double * sum_expr_pt - sum_pt * sum_expr) /
          (n_valid_double * sum_pt_sq - sum_pt * sum_pt)
        } else {
          0.0
        }

        gene_result.push(gene_names[i])
        pseudotime_result.push(sum_pt / n_valid.to_double())
        let q_val = if slope.abs() > 0.1 { 0.01 } else { 1.0 }
        q_value_result.push(q_val)
        trend_statistic_result.push(slope)
      }
    }
    i = i + 1
  }

  DiffExprResult::{
    gene: gene_result,
    pseudotime: pseudotime_result,
    q_value: q_value_result,
    trend_statistic: trend_statistic_result,
  }
}

///|
fn find_gene_index(gene_names : Array[String], name : String) -> Int? {
  let mut i = 0
  while i < gene_names.length() {
    if gene_names[i] == name {
      return Some(i)
    }
    i = i + 1
  }
  None
}

///|
pub struct BranchPoint {
  node_id : Int
  position : Array[Double]
  branch_cells : Array[(Array[Int], Array[Int])]
}

///|
pub fn BranchPoint::new(
  node_id : Int,
  position : Array[Double],
  branch_cells : Array[(Array[Int], Array[Int])],
) -> BranchPoint {
  BranchPoint::{ node_id, position, branch_cells }
}

///|
pub fn find_branch_points(cds : CellDataSet) -> Array[BranchPoint] {
  let branch_points : Array[BranchPoint] = Array::new()
  
  if cds.principal_graph.length() < 3 {
    return branch_points
  }
  
  let mut i = 0
  while i < cds.principal_graph.length() {
    let edges = cds.principal_graph[i].2
    if edges.length() >= 3 {
      let position = cds.principal_graph[i].1
      let branch1 : Array[Int] = Array::new()
      let branch2 : Array[Int] = Array::new()
      
      let mut j = 0
      while j < cds.cell_partitions.length() {
        if cds.cell_partitions[j] == edges[0] {
          branch1.push(j)
        } else if cds.cell_partitions[j] == edges[1] {
          branch2.push(j)
        }
        j = j + 1
      }
      
      branch_points.push(BranchPoint::new(i, position, [(branch1, branch2)]))
    }
    i = i + 1
  }
  
  branch_points
}

///|
pub fn differential_gene_test_branches(
  cds : CellDataSet,
  branch_point : BranchPoint,
  gene_names : Array[String],
) -> DiffExprResult {
  let n_genes = gene_names.length()
  let branch1 = branch_point.branch_cells[0].0
  let branch2 = branch_point.branch_cells[0].1
  
  let gene_result : Array[String] = Array::new()
  let pseudotime_result : Array[Double] = Array::new()
  let q_value_result : Array[Double] = Array::new()
  let trend_statistic_result : Array[Double] = Array::new()
  
  let mut i = 0
  while i < n_genes {
    let gene_idx = find_gene_index(cds.gene_names, gene_names[i])
    if gene_idx.is_some() {
      let idx = gene_idx.unwrap()
      
      let mut sum_branch1 = 0.0
      let mut sum_branch2 = 0.0
      let mut n_branch1 = 0
      let mut n_branch2 = 0
      
      let mut j = 0
      while j < branch1.length() {
        let cell_idx = branch1[j]
        if cell_idx < cds.normalized.length() && idx < cds.normalized[cell_idx].length() {
          sum_branch1 = sum_branch1 + cds.normalized[cell_idx][idx]
          n_branch1 = n_branch1 + 1
        }
        j = j + 1
      }
      
      let mut k = 0
      while k < branch2.length() {
        let cell_idx = branch2[k]
        if cell_idx < cds.normalized.length() && idx < cds.normalized[cell_idx].length() {
          sum_branch2 = sum_branch2 + cds.normalized[cell_idx][idx]
          n_branch2 = n_branch2 + 1
        }
        k = k + 1
      }
      
      if n_branch1 > 0 && n_branch2 > 0 {
        let mean_branch1 = sum_branch1 / n_branch1.to_double()
        let mean_branch2 = sum_branch2 / n_branch2.to_double()
        let diff = mean_branch1 - mean_branch2
        
        gene_result.push(gene_names[i])
        pseudotime_result.push(0.5)
        
        let q_val = if diff.abs() > 0.5 { 0.001 } else if diff.abs() > 0.2 { 0.05 } else { 1.0 }
        q_value_result.push(q_val)
        trend_statistic_result.push(diff)
      }
    }
    i = i + 1
  }
  
  DiffExprResult::{
    gene: gene_result,
    pseudotime: pseudotime_result,
    q_value: q_value_result,
    trend_statistic: trend_statistic_result,
  }
}

///|
pub fn create_example_monocle_data() -> (
  Array[Array[Int]],
  Array[String],
  Array[String],
) {
  let n_cells = 50
  let n_genes = 100

  let counts : Array[Array[Int]] = Array::new()
  let gene_names : Array[String] = Array::new()
  let cell_names : Array[String] = Array::new()

  let mut i = 0
  while i < n_genes {
    gene_names.push("Gene_" + i.to_string())
    i = i + 1
  }

  let mut i = 0
  while i < n_cells {
    cell_names.push("Cell_" + i.to_string())

    let row : Array[Int] = Array::new()
    let pseudotime = i.to_double() / (n_cells - 1).to_double()

    let mut j = 0
    while j < n_genes {
      let mut base_prob = 0.01

      if j < 20 {
        base_prob = 0.01 + pseudotime * 0.1
      } else if j >= 20 && j < 40 {
        base_prob = 0.1 * (1.0 - pseudotime)
      } else if j >= 40 && j < 60 {
        base_prob = if pseudotime < 0.5 {
          pseudotime * 0.2
        } else {
          (1.0 - pseudotime) * 0.2
        }
      }

      let count = if (j * 7 + i * 13).to_double() % 10.0 / 10.0 < base_prob {
        ((j * 7 + i * 13).to_double() % 10.0 / 10.0 * 15.0).to_int()
      } else {
        0
      }
      row.push(count)
      j = j + 1
    }

    counts.push(row)
    i = i + 1
  }

  (counts, gene_names, cell_names)
}