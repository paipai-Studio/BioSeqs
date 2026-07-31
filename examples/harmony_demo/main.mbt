///|
/// Harmony batch correction demo.
/// Demonstrates iterative batch correction on synthetic single-cell data.

fn main {
  let n_cells = 100
  let n_batches = 3
  let n_dims = 10

  println("=== Harmony Batch Correction Demo ===")
  println("")

  // Create example data with batch effects
  let data = @bio.harmony_create_example(n_cells, n_batches, n_dims)

  println("Input data: \{n_cells.to_string()} cells, \{n_batches.to_string()} batches, \{n_dims.to_string()} dimensions")

  // Run Harmony
  let params = @bio.HarmonyParams::create(
    n_clusters = 10,
    max_iterations = 30,
    tolerance = 0.001,
  )

  println("")
  println("Running Harmony...")
  let result = @bio.harmony_run(data, params)

  println("")
  println("Harmony results:")
  println("  Iterations: \{result.n_iterations.to_string()}")
  println("  Converged: \{result.converged.to_string()}")

  // Print batch centroids before/after (first 3 dims)
  println("")
  println("Batch centroids before correction (dim 0-2):")
  {
    let batch_sums : Array[Array[Double]] = Array::make(n_batches, Array::make(3, 0.0))
    let batch_counts : Array[Int] = Array::make(n_batches, 0)
    let mut i = 0
    while i < n_cells {
      let batch = data.batch_labels[i]
      let batch_idx = if batch == "batch_0" { 0 } else if batch == "batch_1" { 1 } else { 2 }
      batch_counts[batch_idx] = batch_counts[batch_idx] + 1
      batch_sums[batch_idx][0] = batch_sums[batch_idx][0] + data.embeddings[i][0]
      batch_sums[batch_idx][1] = batch_sums[batch_idx][1] + data.embeddings[i][1]
      batch_sums[batch_idx][2] = batch_sums[batch_idx][2] + data.embeddings[i][2]
      i = i + 1
    }
    let mut b = 0
    while b < n_batches {
      let name = "batch_" + b.to_string()
      let c0 = if batch_counts[b] > 0 { batch_sums[b][0] / batch_counts[b].to_double() } else { 0.0 }
      let c1 = if batch_counts[b] > 0 { batch_sums[b][1] / batch_counts[b].to_double() } else { 0.0 }
      let c2 = if batch_counts[b] > 0 { batch_sums[b][2] / batch_counts[b].to_double() } else { 0.0 }
      println("  \{name}: [\{c0.to_string()}, \{c1.to_string()}, \{c2.to_string()}]")
      b = b + 1
    }
  }

  // Compute batch centroids after correction
  println("")
  println("Batch centroids after correction (dim 0-2):")
  {
    let batch_sums : Array[Array[Double]] = Array::make(n_batches, Array::make(3, 0.0))
    let batch_counts : Array[Int] = Array::make(n_batches, 0)
    let mut i = 0
    while i < n_cells {
      let batch = data.batch_labels[i]
      let batch_idx = if batch == "batch_0" { 0 } else if batch == "batch_1" { 1 } else { 2 }
      batch_counts[batch_idx] = batch_counts[batch_idx] + 1
      batch_sums[batch_idx][0] = batch_sums[batch_idx][0] + result.corrected[i][0]
      batch_sums[batch_idx][1] = batch_sums[batch_idx][1] + result.corrected[i][1]
      batch_sums[batch_idx][2] = batch_sums[batch_idx][2] + result.corrected[i][2]
      i = i + 1
    }
    let mut b = 0
    while b < n_batches {
      let name = "batch_" + b.to_string()
      let c0 = if batch_counts[b] > 0 { batch_sums[b][0] / batch_counts[b].to_double() } else { 0.0 }
      let c1 = if batch_counts[b] > 0 { batch_sums[b][1] / batch_counts[b].to_double() } else { 0.0 }
      let c2 = if batch_counts[b] > 0 { batch_sums[b][2] / batch_counts[b].to_double() } else { 0.0 }
      println("  \{name}: [\{c0.to_string()}, \{c1.to_string()}, \{c2.to_string()}]")
      b = b + 1
    }
  }

  println("")
  println("Done!")
}
