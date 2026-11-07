// scvi_cellassign.mbt - CellAssign: probabilistic cell type assignment
// Inspired by scvi-tools CellAssign model (Zhang et al. 2019)
// Uses marker gene x cell type binary matrix for annotation

///|
pub struct CellAssignConfig {
  learning_rate : Double
  max_epochs : Int32
  delta : Double
  s : Double // Scale parameter
}

///|
pub struct CellAssignModel {
  config : CellAssignConfig
  // Marker matrix: gene x cell_type (binary)
  marker_matrix : Array[Array[Double]]
  // Cell type prior log-probabilities
  log_prior : Array[Double]
  // Per-gene baseline expression
  log_baseline : Array[Double]
  // Per-cell size factors
  size_factors : Array[Double]
  n_genes : Int32
  n_cell_types : Int32
  mut trained : Bool
}

///|
pub fn cellassign_init(
  marker_matrix : Array[Array[Double]],
  n_cells : Int32,
  config : CellAssignConfig,
) -> CellAssignModel {
  let n_genes = marker_matrix.length()
  let n_cell_types = if n_genes > 0 { marker_matrix[0].length() } else { 0 }
  let log_prior = Array(capacity=n_cell_types)
    |> a => {
      for _i in 0..<n_cell_types {
        a.push(-@math.ln(n_cell_types.to_double())) // Uniform prior
      }
      a
    }
  let log_baseline = Array(capacity=n_genes)
    |> a => {
      for _i in 0..<n_genes {
        a.push(0.0)
      }
      a
    }
  let size_factors = Array(capacity=n_cells)
    |> a => {
      for _i in 0..<n_cells {
        a.push(1.0)
      }
      a
    }
  {
    config,
    marker_matrix,
    log_prior,
    log_baseline,
    size_factors,
    n_genes,
    n_cell_types,
    trained: false,
  }
}

///|
/// Compute log-likelihood for each cell type assignment
pub fn cellassign_log_likelihood(
  model : CellAssignModel,
  x : Array[Double],
  cell_idx : Int32,
) -> Array[Double] {
  let n_ct = model.n_cell_types
  let log_liks = Array(capacity=n_ct)
  for ct in 0..<n_ct {
    let mut ll = model.log_prior[ct]
    for g in 0..<model.n_genes {
      let marker = model.marker_matrix[g][ct]
      // Expected expression: baseline * s * (1 + delta * marker)
      let log_mu = model.log_baseline[g] +
        @math.ln(model.size_factors[cell_idx]) +
        @math.ln(1.0 + model.config.delta * marker)
      let mu = @math.exp(log_mu)
      // Poisson log-likelihood
      if x[g] > 0.0 {
        ll += x[g] * log_mu - mu
      } else {
        ll += -mu
      }
    }
    log_liks.push(ll)
  }
  log_liks
}

///|
/// Predict cell types
pub fn cellassign_predict(
  model : CellAssignModel,
  x : Array[Array[Double]],
) -> Array[Int32] {
  let result = Array(capacity=x.length())
  for i in 0..<x.length() {
    let log_liks = cellassign_log_likelihood(model, x[i], i)
    let mut best_idx = 0
    let mut best_ll = log_liks[0]
    for ct in 1..<log_liks.length() {
      if log_liks[ct] > best_ll {
        best_ll = log_liks[ct]
        best_idx = ct
      }
    }
    result.push(best_idx)
  }
  result
}

///|
/// Get assignment probabilities
pub fn cellassign_predict_proba(
  model : CellAssignModel,
  x : Array[Double],
  cell_idx : Int32,
) -> Array[Double] {
  let log_liks = cellassign_log_likelihood(model, x, cell_idx)
  // Softmax
  let max_ll = scvi_array_max(log_liks)
  let mut sum_exp = 0.0
  for ct in 0..<log_liks.length() {
    sum_exp += @math.exp(log_liks[ct] - max_ll)
  }
  let result = Array(capacity=log_liks.length())
  for ct in 0..<log_liks.length() {
    result.push(@math.exp(log_liks[ct] - max_ll) / sum_exp)
  }
  result
}

///|
/// Set size factors from library sizes
pub fn cellassign_set_size_factors(
  model : CellAssignModel,
  x : Array[Array[Double]],
) -> CellAssignModel {
  let m = model
  for i in 0..<x.length() {
    m.size_factors[i] = scvi_array_sum(x[i])
  }
  m
}

///|
/// Estimate baseline expression from data
pub fn cellassign_estimate_baseline(
  model : CellAssignModel,
  x : Array[Array[Double]],
) -> CellAssignModel {
  let m = model
  let n_cells = x.length()
  for g in 0..<m.n_genes {
    let mut sum = 0.0
    for i in 0..<n_cells {
      sum += x[i][g]
    }
    let mean = sum / n_cells.to_double()
    m.log_baseline[g] = @math.ln(mean + 0.0000000001)
  }
  m
}

///|
/// Train CellAssign model
pub fn cellassign_train(
  model : CellAssignModel,
  x : Array[Array[Double]],
) -> CellAssignModel {
  let mut m = model
  // Estimate baseline
  m = cellassign_estimate_baseline(m, x)
  // Set size factors
  m = cellassign_set_size_factors(m, x)
  // Update priors based on data
  let n_cells = x.length()
  let predictions = cellassign_predict(m, x)
  let ct_counts = Array(capacity=m.n_cell_types)
    |> a => {
      for _i in 0..<m.n_cell_types {
        a.push(0.0)
      }
      a
    }
  for i in 0..<n_cells {
    ct_counts[predictions[i]] += 1.0
  }
  for ct in 0..<m.n_cell_types {
    m.log_prior[ct] = @math.ln(
      ct_counts[ct] / n_cells.to_double() + 0.0000000001,
    )
  }
  m.trained = true
  m
}

///|
/// Create marker matrix from gene/cell-type associations
pub fn cellassign_create_marker_matrix(
  gene_names : Array[String],
  marker_dict : Map[String, Array[String]],
  cell_type_names : Array[String],
) -> Array[Array[Double]] {
  let n_genes = gene_names.length()
  let n_ct = cell_type_names.length()
  let ct_index = Map([], capacity=n_ct)
    |> m => {
      for i in 0..<n_ct {
        m[cell_type_names[i]] = i
      }
      m
    }
  let result = Array(capacity=n_genes)
    |> a => {
      for g in 0..<n_genes {
        let row = Array(capacity=n_ct)
          |> r => {
            for _i in 0..<n_ct {
              r.push(0.0)
            }
            r
          }
        match marker_dict.get(gene_names[g]) {
          Some(cts) =>
            for j in 0..<cts.length() {
              match ct_index.get(cts[j]) {
                Some(idx) => row[idx] = 1.0
                None => ()
              }
            }
          None => ()
        }
        a.push(row)
      }
      a
    }
  result
}

///|
pub fn CellAssignConfig::default() -> CellAssignConfig {
  { learning_rate: 0.0, max_epochs: 0, delta: 0.0, s: 0.0, }
}
