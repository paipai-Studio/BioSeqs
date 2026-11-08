// scvi_stereoscope.mbt - Stereoscope-style spatial deconvolution
// Mirrors scvi.external.RNAStereoscope and SpatialStereoscope: learn cell-type
// expression profiles from a single-cell reference, then estimate per-spot
// cell-type proportions with multiplicative (NMF-style) updates.

///|
pub struct ScviStereoscopeConfig {
  n_epochs : Int32
  learning_rate : Double
  tolerance : Double
}

///|
pub fn ScviStereoscopeConfig::new(
  n_epochs : Int32,
  learning_rate : Double,
  tolerance : Double,
) -> ScviStereoscopeConfig {
  { n_epochs, learning_rate, tolerance, }
}

///|
pub fn ScviStereoscopeConfig::default() -> ScviStereoscopeConfig {
  { n_epochs: 200, learning_rate: 1.0, tolerance: 1.0e-6, }
}

///|
pub struct ScviStereoscopeModel {
  profiles : Array[Array[Double]] // [cell_type][gene], rows sum to one
  n_cell_types : Int32
  n_genes : Int32
}

///|
pub fn ScviStereoscopeModel::n_cell_types(self : ScviStereoscopeModel) -> Int32 {
  self.n_cell_types
}

///|
pub fn ScviStereoscopeModel::n_genes(self : ScviStereoscopeModel) -> Int32 {
  self.n_genes
}

///|
fn scvi_stereoscope_normalise(row : Array[Double]) -> Array[Double] {
  let mut s = 0.0
  for i in 0..<row.length() {
    let v = if row[i] < 0.0 { 0.0 } else { row[i] }
    s += v
  }
  let denom = if s <= 0.0 { 1.0 } else { s }
  let out = Array(capacity=row.length())
  for i in 0..<row.length() {
    let v = if row[i] < 0.0 { 0.0 } else { row[i] }
    out.push(v / denom)
  }
  out
}

///|
/// Estimate reference profiles as the mean expression of each cell type.
/// A small pseudo-count keeps profiles strictly positive.
pub fn scvi_stereoscope_fit(
  reference : Array[Array[Double]],
  cell_type : Array[Int32],
  n_cell_types : Int32,
  config : ScviStereoscopeConfig,
) -> ScviStereoscopeModel {
  let n_genes = if reference.length() == 0 { 0 } else { reference[0].length() }
  let sums = Array(capacity=n_cell_types)
  let counts = Array::make(n_cell_types, 0)
  for _t in 0..<n_cell_types {
    sums.push(Array::make(n_genes, 0.0))
  }
  for i in 0..<reference.length() {
    if i < cell_type.length() {
      let t = cell_type[i]
      if t >= 0 && t < n_cell_types {
        for g in 0..<n_genes {
          if g < reference[i].length() {
            sums[t][g] = sums[t][g] + reference[i][g]
          }
        }
        counts[t] = counts[t] + 1
      }
    }
  }
  let profiles = Array(capacity=n_cell_types)
  let pseudo = 1.0e-8
  for t in 0..<n_cell_types {
    let denom = if counts[t] == 0 { 1.0 } else { counts[t].to_double() }
    let row = Array(capacity=n_genes)
    for g in 0..<n_genes {
      row.push(sums[t][g] / denom + pseudo)
    }
    profiles.push(scvi_stereoscope_normalise(row))
  }
  ignore(config)
  { profiles, n_cell_types, n_genes, }
}

///|
/// Deconvolve one bulk / spot vector into cell-type proportions.
pub fn scvi_stereoscope_deconvolve(
  model : ScviStereoscopeModel,
  bulk : Array[Double],
  config : ScviStereoscopeConfig,
) -> Array[Double] {
  let n_types = model.n_cell_types
  if n_types <= 0 {
    return Array(capacity=0)
  }
  let w = Array::make(n_types, 1.0 / n_types.to_double())
  let n_genes = model.n_genes
  let eps = 1.0e-9
  for _iter in 0..<config.n_epochs {
    let pred = Array::make(n_genes, eps)
    for t in 0..<n_types {
      let p = model.profiles[t]
      for g in 0..<n_genes {
        if g < p.length() && g < bulk.length() {
          pred[g] = pred[g] + w[t] * p[g]
        }
      }
    }
    for t in 0..<n_types {
      let p = model.profiles[t]
      let mut num = 0.0
      for g in 0..<n_genes {
        if g < p.length() && g < bulk.length() {
          let y = if bulk[g] < 0.0 { 0.0 } else { bulk[g] }
          num += p[g] * y / pred[g]
        }
      }
      w[t] = w[t] * num
    }
    let mut s = 0.0
    for t in 0..<n_types {
      s += w[t]
    }
    let denom = if s <= 0.0 { 1.0 } else { s }
    let mut max_change = 0.0
    for t in 0..<n_types {
      let old = w[t]
      w[t] = w[t] / denom
      let change = (w[t] - old / denom).abs()
      if change > max_change {
        max_change = change
      }
    }
    if max_change < config.tolerance {
      break
    }
  }
  w
}

///|
/// Deconvolve a matrix of spatial spots (`[spot][gene]`).
pub fn scvi_stereoscope_spatial_deconvolve(
  model : ScviStereoscopeModel,
  spots : Array[Array[Double]],
  config : ScviStereoscopeConfig,
) -> Array[Array[Double]] {
  let out = Array(capacity=spots.length())
  for i in 0..<spots.length() {
    out.push(scvi_stereoscope_deconvolve(model, spots[i], config))
  }
  out
}

///|
/// Predicted bulk expression from cell-type proportions.
pub fn scvi_stereoscope_predict_expression(
  model : ScviStereoscopeModel,
  proportions : Array[Double],
) -> Array[Double] {
  let out = Array::make(model.n_genes, 0.0)
  for t in 0..<model.n_cell_types {
    let w = if t < proportions.length() { proportions[t] } else { 0.0 }
    let p = model.profiles[t]
    for g in 0..<model.n_genes {
      if g < p.length() {
        out[g] = out[g] + w * p[g]
      }
    }
  }
  out
}

///|
/// Cell type with the largest estimated proportion for each spot.
pub fn scvi_stereoscope_dominant_type(
  proportions : Array[Array[Double]],
) -> Array[Int32] {
  let out = Array(capacity=proportions.length())
  for i in 0..<proportions.length() {
    let row = proportions[i]
    let mut best = 0
    let mut best_v = if row.length() == 0 { 0.0 } else { row[0] }
    for t in 1..<row.length() {
      if row[t] > best_v {
        best_v = row[t]
        best = t
      }
    }
    out.push(best)
  }
  out
}

// ===========================================================================
// scvi.external aliases
// ===========================================================================

///|
pub fn scvi_rnastereoscope_fit(
  reference : Array[Array[Double]],
  cell_type : Array[Int32],
  n_cell_types : Int32,
  config : ScviStereoscopeConfig,
) -> ScviStereoscopeModel {
  scvi_stereoscope_fit(reference, cell_type, n_cell_types, config)
}

///|
pub fn scvi_spatialstereoscope_fit(
  reference : Array[Array[Double]],
  cell_type : Array[Int32],
  n_cell_types : Int32,
  config : ScviStereoscopeConfig,
) -> ScviStereoscopeModel {
  scvi_stereoscope_fit(reference, cell_type, n_cell_types, config)
}

///|
pub fn scvi_rnastereoscope_deconvolve(
  model : ScviStereoscopeModel,
  bulk : Array[Double],
  config : ScviStereoscopeConfig,
) -> Array[Double] {
  scvi_stereoscope_deconvolve(model, bulk, config)
}

///|
pub fn scvi_spatialstereoscope_deconvolve(
  model : ScviStereoscopeModel,
  spots : Array[Array[Double]],
  config : ScviStereoscopeConfig,
) -> Array[Array[Double]] {
  scvi_stereoscope_spatial_deconvolve(model, spots, config)
}
