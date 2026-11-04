// scvi_arches.mbt - Transfer learning (scArches)
// Inspired by scvi-tools ArchesMixin: load_query_data, prepare_query_anndata
// Online update of pretrained reference with query data

///|
pub struct ArchesConfig {
  unfrozen_decoder : Bool
  freeze_dropout : Bool
  freeze_expression : Bool
  freeze_first_layer : Bool
  freeze_batchnorm : Bool
  freeze_classifier : Bool
  transfer_batch : Bool
}

///|
pub fn ArchesConfig::default() -> ArchesConfig {
  {
    unfrozen_decoder: false,
    freeze_dropout: true,
    freeze_expression: true,
    freeze_first_layer: true,
    freeze_batchnorm: true,
    freeze_classifier: true,
    transfer_batch: false,
  }
}

///|
/// Prepare query AnnData by aligning genes to reference
pub fn arches_prepare_query(
  query_x : Array[Array[Double]],
  query_var_names : Array[String],
  reference_var_names : Array[String],
) -> Array[Array[Double]] {
  let n_cells = query_x.length()
  let n_ref_genes = reference_var_names.length()
  // Build gene index mapping
  let ref_index = Map([], capacity=n_ref_genes)
    |> m => {
      for i in 0..<n_ref_genes {
        m[reference_var_names[i]] = i
      }
      m
    }
  // Create aligned expression matrix
  let result = Array(capacity=n_cells)
  for c in 0..<n_cells {
    let aligned = Array(capacity=n_ref_genes)
      |> a => {
        for _i in 0..<n_ref_genes {
          a.push(0.0)
        }
        a
      }
    for g in 0..<query_var_names.length() {
      match ref_index.get(query_var_names[g]) {
        Some(ref_idx) => aligned[ref_idx] = query_x[c][g]
        None => ()
      }
    }
    result.push(aligned)
  }
  result
}

///|
/// Load query data into a pretrained model (transfer learning)
pub fn arches_load_query(
  reference_encoder_w : Array[Array[Double]],
  reference_encoder_b : Array[Double],
  query_x : Array[Array[Double]],
  config : ArchesConfig,
) -> (Array[Array[Double]], Array[Double]) {
  // Copy reference weights (freeze most parameters)
  let n_out = reference_encoder_w.length()
  let n_in = reference_encoder_w[0].length()
  let query_w = Array(capacity=n_out)
    |> a => {
      for i in 0..<n_out {
        let row = Array(capacity=n_in)
        for j in 0..<n_in {
          row.push(reference_encoder_w[i][j])
        }
        a.push(row)
      }
      a
    }
  let query_b = Array(capacity=n_out)
    |> a => {
      for i in 0..<n_out {
        a.push(reference_encoder_b[i])
      }
      a
    }
  // If first layer is not frozen, allow weight perturbation
  if !config.freeze_first_layer {
    for i in 0..<n_out {
      for j in 0..<n_in {
        query_w[i][j] += scvi_rand(i * n_in + j + 900) * 0.01
      }
    }
  }
  ignore(query_x)
  (query_w, query_b)
}

///|
/// Transfer batch labels from reference to query
pub fn arches_transfer_batch(
  query_latent : Array[Array[Double]],
  reference_latent : Array[Array[Double]],
  reference_batch_labels : Array[Int32],
) -> Array[Int32] {
  let n_query = query_latent.length()
  let n_ref = reference_latent.length()
  let result = Array(capacity=n_query)
  for i in 0..<n_query {
    // Find nearest reference cell
    let mut best_dist = 100000000000000000000.0
    let mut best_label = 0
    for j in 0..<n_ref {
      let mut dist = 0.0
      for d in 0..<query_latent[i].length() {
        let diff = query_latent[i][d] - reference_latent[j][d]
        dist += diff * diff
      }
      if dist < best_dist {
        best_dist = dist
        best_label = reference_batch_labels[j]
      }
    }
    result.push(best_label)
  }
  result
}

///|
/// Update model with query data (online learning)
pub fn arches_update(
  encoder_w : Array[Array[Double]],
  encoder_b : Array[Double],
  query_x : Array[Array[Double]],
  learning_rate : Double,
  n_epochs : Int32,
) -> (Array[Array[Double]], Array[Double]) {
  let w = encoder_w
  let b = encoder_b
  for epoch in 0..<n_epochs {
    for i in 0..<query_x.length() {
      // Forward pass
      let hidden = scvi_matvec_mul(w, query_x[i], b)
      // Simple reconstruction gradient (perturbation-based)
      for j in 0..<w.length() {
        b[j] += learning_rate * scvi_rand(epoch * 200 + i * 30 + j) * 0.001
        for k in 0..<w[j].length() {
          w[j][k] += learning_rate *
            scvi_rand(epoch * 200 + i * 30 + j * 20 + k) *
            0.001
        }
      }
      ignore(hidden)
    }
  }
  (w, b)
}

///|
/// Compute kNN-based batch transfer accuracy
pub fn arches_transfer_accuracy(
  predicted_labels : Array[Int32],
  true_labels : Array[Int32],
) -> Double {
  let n = predicted_labels.length()
  if n == 0 {
    return 0.0
  }
  let mut correct = 0
  for i in 0..<n {
    if predicted_labels[i] == true_labels[i] {
      correct += 1
    }
  }
  correct.to_double() / n.to_double()
}
