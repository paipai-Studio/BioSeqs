// scvi_tangram.mbt - Tangram: spatial mapping of single cells
// Inspired by scvi-tools Tangram model (Biancalani et al. 2021)
// Maps single-cell profiles onto spatial transcriptomics via optimal transport

///|
pub struct TangramConfig {
  learning_rate : Double
  max_epochs : Int32
  lambda_g1 : Double // Gene expression similarity weight
  lambda_g2 : Double // Marker gene weight
  lambda_s : Double // Sparsity weight
}

///|
pub struct TangramModel {
  config : TangramConfig
  // Mapping matrix: cell x spot (probability)
  mapping : Array[Array[Double]]
  // Spatial gene expression (spot x gene)
  spatial_genes : Array[Array[Double]]
  // Single-cell gene expression (cell x gene)
  sc_genes : Array[Array[Double]]
  n_cells : Int32
  n_spots : Int32
  n_genes : Int32
  mut trained : Bool
}

///|
pub fn tangram_init(
  sc_genes : Array[Array[Double]],
  spatial_genes : Array[Array[Double]],
  config : TangramConfig,
) -> TangramModel {
  let n_cells = sc_genes.length()
  let n_spots = spatial_genes.length()
  let n_genes = if n_cells > 0 { sc_genes[0].length() } else { 0 }
  // Initialize mapping uniformly
  let mapping = Array(capacity=n_cells)
    |> a => {
      for _c in 0..<n_cells {
        let row = Array(capacity=n_spots)
        for _s in 0..<n_spots {
          row.push(1.0 / n_spots.to_double())
        }
        a.push(row)
      }
      a
    }
  {
    config,
    mapping,
    spatial_genes,
    sc_genes,
    n_cells,
    n_spots,
    n_genes,
    trained: false,
  }
}

///|
/// Compute mapping loss (gene expression similarity)
pub fn tangram_compute_loss(model : TangramModel) -> Double {
  // Mapped expression = mapping^T @ sc_genes
  let mapped = Array(capacity=model.n_spots)
    |> a => {
      for _s in 0..<model.n_spots {
        a.push(
          Array(capacity=model.n_genes)
          |> g => {
            for _i in 0..<model.n_genes {
              g.push(0.0)
            }
            g
          },
        )
      }
      a
    }
  for s in 0..<model.n_spots {
    for c in 0..<model.n_cells {
      for g in 0..<model.n_genes {
        mapped[s][g] += model.mapping[c][s] * model.sc_genes[c][g]
      }
    }
  }
  // Cosine similarity loss
  let mut loss = 0.0
  for g in 0..<model.n_genes {
    let mut dot = 0.0
    let mut norm_mapped = 0.0
    let mut norm_spatial = 0.0
    for s in 0..<model.n_spots {
      dot += mapped[s][g] * model.spatial_genes[s][g]
      norm_mapped += mapped[s][g] * mapped[s][g]
      norm_spatial += model.spatial_genes[s][g] * model.spatial_genes[s][g]
    }
    if norm_mapped > 0.0 && norm_spatial > 0.0 {
      loss += 1.0 - dot / (norm_mapped.sqrt() * norm_spatial.sqrt())
    }
  }
  loss
}

///|
/// Update mapping (gradient-free optimization)
pub fn tangram_update_mapping(model : TangramModel) -> TangramModel {
  let m = model
  for c in 0..<m.n_cells {
    for s in 0..<m.n_spots {
      // Compute similarity between cell c and spot s
      let mut dot = 0.0
      let mut norm_c = 0.0
      let mut norm_s = 0.0
      for g in 0..<m.n_genes {
        dot += m.sc_genes[c][g] * m.spatial_genes[s][g]
        norm_c += m.sc_genes[c][g] * m.sc_genes[c][g]
        norm_s += m.spatial_genes[s][g] * m.spatial_genes[s][g]
      }
      let sim = if norm_c > 0.0 && norm_s > 0.0 {
        dot / (norm_c.sqrt() * norm_s.sqrt())
      } else {
        0.0
      }
      // Update mapping based on similarity
      m.mapping[c][s] = @math.exp(sim * 10.0) // Soft assignment
    }
    // Normalize rows
    let mut row_sum = 0.0
    for s in 0..<m.n_spots {
      row_sum += m.mapping[c][s]
    }
    if row_sum > 0.0 {
      for s in 0..<m.n_spots {
        m.mapping[c][s] /= row_sum
      }
    }
  }
  m
}

///|
/// Train Tangram model
pub fn tangram_train(model : TangramModel, n_epochs : Int32) -> TangramModel {
  let mut m = model
  for epoch in 0..<n_epochs {
    m = tangram_update_mapping(m)
    let loss = tangram_compute_loss(m)
    ignore(epoch)
    ignore(loss)
  }
  m.trained = true
  m
}

///|
/// Get mapped single cells onto spatial coordinates
pub fn tangram_get_mapped_cells(model : TangramModel) -> Array[Array[Double]] {
  model.mapping
}

///|
/// Get cell density per spot
pub fn tangram_cell_density(model : TangramModel) -> Array[Double] {
  let result = Array(capacity=model.n_spots)
  for s in 0..<model.n_spots {
    let mut density = 0.0
    for c in 0..<model.n_cells {
      density += model.mapping[c][s]
    }
    result.push(density)
  }
  result
}

///|
/// Get spatial gene expression from mapped cells
pub fn tangram_get_spatial_expression(
  model : TangramModel,
) -> Array[Array[Double]] {
  let result = Array(capacity=model.n_spots)
    |> a => {
      for _s in 0..<model.n_spots {
        a.push(
          Array(capacity=model.n_genes)
          |> g => {
            for _i in 0..<model.n_genes {
              g.push(0.0)
            }
            g
          },
        )
      }
      a
    }
  for s in 0..<model.n_spots {
    for c in 0..<model.n_cells {
      for g in 0..<model.n_genes {
        result[s][g] += model.mapping[c][s] * model.sc_genes[c][g]
      }
    }
  }
  result
}

///|
/// Compare mapped vs observed spatial expression
pub fn tangram_compare_expression(model : TangramModel) -> Array[Double] {
  let mapped = tangram_get_spatial_expression(model)
  let result = Array(capacity=model.n_genes)
  for g in 0..<model.n_genes {
    let mut mapped_sum = 0.0
    let mut obs_sum = 0.0
    for s in 0..<model.n_spots {
      mapped_sum += mapped[s][g]
      obs_sum += model.spatial_genes[s][g]
    }
    // Correlation-like score
    if obs_sum > 0.0 {
      result.push(mapped_sum / obs_sum)
    } else {
      result.push(0.0)
    }
  }
  result
}

///|
pub fn TangramConfig::default() -> TangramConfig {
  {
    learning_rate: 0.0,
    max_epochs: 0,
    lambda_g1: 0.0,
    lambda_g2: 0.0,
    lambda_s: 0.0,
  }
}
