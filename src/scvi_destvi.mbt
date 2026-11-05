// scvi_destvi.mbt - Spatial deconvolution (DestVI + Stereoscope)
// Inspired by scvi-tools DestVI and Stereoscope models
// Deconvolve spatial spots into cell type proportions

///|
pub struct DestVIConfig {
  n_latent : Int32
  n_cell_types : Int32
  learning_rate : Double
  max_epochs : Int32
}

///|
pub struct DestVIModel {
  config : DestVIConfig
  // Reference cell-type signatures (gene x cell_type)
  cell_type_profiles : Array[Array[Double]]
  // Spatial spot proportions (spot x cell_type)
  mut spot_proportions : Array[Array[Double]]
  // Cell-type scale factors
  scale_factors : Array[Double]
  n_genes : Int32
  n_spots : Int32
  mut trained : Bool
}

///|
/// Build reference cell-type profiles from single-cell data
pub fn destvi_build_reference(
  sc_x : Array[Array[Double]],
  cell_type_labels : Array[Int32],
  n_cell_types : Int32,
) -> Array[Array[Double]] {
  let n_genes = sc_x[0].length()
  let profiles = Array(capacity=n_genes)
    |> a => {
      for _g in 0..<n_genes {
        a.push(
          Array(capacity=n_cell_types)
          |> ct => {
            for _c in 0..<n_cell_types {
              ct.push(0.0)
            }
            ct
          },
        )
      }
      a
    }
  let counts = Array(capacity=n_cell_types)
    |> a => {
      for _i in 0..<n_cell_types {
        a.push(0)
      }
      a
    }
  for i in 0..<sc_x.length() {
    let ct = cell_type_labels[i]
    counts[ct] += 1
    for g in 0..<n_genes {
      profiles[g][ct] += sc_x[i][g]
    }
  }
  // Normalize to get mean expression per cell type
  for ct in 0..<n_cell_types {
    if counts[ct] > 0 {
      for g in 0..<n_genes {
        profiles[g][ct] /= counts[ct].to_double()
      }
    }
  }
  profiles
}

///|
/// Initialize DestVI model
pub fn destvi_init(
  cell_type_profiles : Array[Array[Double]],
  n_spots : Int32,
  config : DestVIConfig,
) -> DestVIModel {
  let n_genes = cell_type_profiles.length()
  let n_ct = config.n_cell_types
  // Initialize proportions uniformly
  let proportions = Array(capacity=n_spots)
    |> a => {
      for _s in 0..<n_spots {
        let row = Array(capacity=n_ct)
        for _c in 0..<n_ct {
          row.push(1.0 / n_ct.to_double())
        }
        a.push(row)
      }
      a
    }
  let scales = Array(capacity=n_ct)
    |> a => {
      for _i in 0..<n_ct {
        a.push(1.0)
      }
      a
    }
  {
    config,
    cell_type_profiles,
    spot_proportions: proportions,
    scale_factors: scales,
    n_genes,
    n_spots,
    trained: false,
  }
}

///|
/// Deconvolve spatial spots
pub fn destvi_deconvolve(
  model : DestVIModel,
  spatial_x : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n_spots = spatial_x.length()
  let n_ct = model.config.n_cell_types
  let result = Array(capacity=n_spots)
  for s in 0..<n_spots {
    let spot_counts = spatial_x[s]
    let lib_size = scvi_array_sum(spot_counts)
    // Estimate proportions by non-negative least squares (simplified)
    let proportions = Array(capacity=n_ct)
      |> a => {
        for _c in 0..<n_ct {
          a.push(0.0)
        }
        a
      }
    // For each cell type, compute correlation with spot profile
    let scores = Array(capacity=n_ct)
    for ct in 0..<n_ct {
      let mut dot = 0.0
      let mut norm_sq = 0.0
      for g in 0..<model.n_genes {
        let expected = model.cell_type_profiles[g][ct] * model.scale_factors[ct]
        dot += spot_counts[g] * expected
        norm_sq += expected * expected
      }
      if norm_sq > 0.0 {
        scores.push(dot / norm_sq)
      } else {
        scores.push(0.0)
      }
    }
    // Normalize scores to get proportions
    let mut total_score = 0.0
    for ct in 0..<n_ct {
      if scores[ct] < 0.0 {
        scores[ct] = 0.0
      }
      total_score += scores[ct]
    }
    if total_score > 0.0 {
      for ct in 0..<n_ct {
        proportions[ct] = scores[ct] / total_score
      }
    } else {
      for ct in 0..<n_ct {
        proportions[ct] = 1.0 / n_ct.to_double()
      }
    }
    ignore(lib_size)
    result.push(proportions)
  }
  result
}

///|
/// Get spatial factor (expected expression per spot)
pub fn destvi_get_spatial_factor(
  model : DestVIModel,
  proportions : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n_spots = proportions.length()
  let result = Array(capacity=n_spots)
  for s in 0..<n_spots {
    let expected = Array(capacity=model.n_genes)
      |> a => {
        for _g in 0..<model.n_genes {
          a.push(0.0)
        }
        a
      }
    for g in 0..<model.n_genes {
      for ct in 0..<model.config.n_cell_types {
        expected[g] += proportions[s][ct] *
          model.cell_type_profiles[g][ct] *
          model.scale_factors[ct]
      }
    }
    result.push(expected)
  }
  result
}

///|
/// Train DestVI model
pub fn destvi_train(
  model : DestVIModel,
  spatial_x : Array[Array[Double]],
  n_epochs : Int32,
) -> DestVIModel {
  let m = model
  for epoch in 0..<n_epochs {
    // Update proportions
    m.spot_proportions = destvi_deconvolve(m, spatial_x)
    // Update scale factors
    for ct in 0..<m.config.n_cell_types {
      let mut total_expected = 0.0
      let mut total_observed = 0.0
      for s in 0..<spatial_x.length() {
        for g in 0..<m.n_genes {
          total_expected += m.spot_proportions[s][ct] *
            m.cell_type_profiles[g][ct]
          total_observed += spatial_x[s][g]
        }
      }
      if total_expected > 0.0 {
        m.scale_factors[ct] = total_observed /
          (total_expected * m.config.n_cell_types.to_double())
      }
    }
    ignore(epoch)
  }
  m.trained = true
  m
}

///|
/// Stereoscope-style deconvolution (simpler NB-based)
pub fn stereoscope_deconvolve(
  reference_profiles : Array[Array[Double]],
  spatial_counts : Array[Double],
  n_cell_types : Int32,
) -> Array[Double] {
  let n_genes = reference_profiles.length()
  let lib_size = scvi_array_sum(spatial_counts)
  // NB-based proportion estimation
  let scores = Array(capacity=n_cell_types)
  for ct in 0..<n_cell_types {
    let mut log_lik = 0.0
    for g in 0..<n_genes {
      let mu = reference_profiles[g][ct] * lib_size / n_genes.to_double()
      if spatial_counts[g] > 0.0 && mu > 0.0 {
        log_lik += spatial_counts[g] * @math.ln(mu + 0.0000000001) - mu
      }
    }
    scores.push(log_lik)
  }
  // Softmax to get proportions
  let max_score = scvi_array_max(scores)
  let mut sum_exp = 0.0
  for ct in 0..<n_cell_types {
    scores[ct] = @math.exp(scores[ct] - max_score)
    sum_exp += scores[ct]
  }
  for ct in 0..<n_cell_types {
    scores[ct] /= sum_exp
  }
  scores
}

///|
pub fn DestVIConfig::default() -> DestVIConfig {
  { n_latent: 0, n_cell_types: 0, learning_rate: 0.0, max_epochs: 0, }
}
