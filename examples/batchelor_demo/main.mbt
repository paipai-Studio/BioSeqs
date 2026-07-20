///|
/// Batchelor Demo - Single-cell batch correction methods.
/// Inspired by Bioconductor batchelor package.

///|
fn main {
  println("=== Batchelor Single-cell Batch Correction Demo ===")
  println("Inspired by Bioconductor batchelor package")
  println("")
  
  println("1. Generating simulated single-cell data with batch effects...")
  let (batches, batch_names) = @bio.batchelor_create_example_data(2, 100, 500, 2.0)
  
  let mut i = 0
  while i < batches.length() {
    println("   \{batch_names[i]}: \{batches[i].length().to_string()} cells, \{batches[i][0].length().to_string()} genes")
    i = i + 1
  }
  println("")
  
  println("2. Evaluating batch mixing before correction...")
  let raw_data : Array[Array[Double]] = Array::new()
  for cell in batches[0] {
    raw_data.push(cell)
  }
  for cell in batches[1] {
    raw_data.push(cell)
  }
  
  let raw_batch_indices : Array[Int] = Array::new()
  let mut idx = 0
  while idx < 100 {
    raw_batch_indices.push(0)
    idx = idx + 1
  }
  idx = 0
  while idx < 100 {
    raw_batch_indices.push(1)
    idx = idx + 1
  }
  
  let raw_score = @bio.compute_batch_mixing_score(raw_data, raw_batch_indices)
  println("   Batch mixing score (higher = better mixing): \{raw_score.to_string()}")
  println("")
  
  println("3. Method 1: rescaleBatches (simple scaling correction)...")
  let rescale_result = @bio.rescale_batches(batches, batch_names)
  let rescale_score = @bio.compute_batch_mixing_score(rescale_result.corrected, rescale_result.batch_indices)
  println("   Corrected dimensions: \{rescale_result.corrected.length().to_string()} cells, \{rescale_result.corrected[0].length().to_string()} genes")
  println("   Batch mixing score after rescaleBatches: \{rescale_score.to_string()}")
  println("   Improvement: +\{((rescale_score - raw_score) / raw_score * 100.0).to_string()}%")
  println("")
  
  println("4. Method 2: fastMNN (Mutual Nearest Neighbors correction)...")
  println("   Performing PCA and MNN correction...")
  let fastmnn_result = @bio.fast_mnn(batches, batch_names, 5, 50)
  println("   Corrected dimensions: \{fastmnn_result.corrected.length().to_string()} cells, \{fastmnn_result.corrected[0].length().to_string()} dimensions")
  println("   Number of MNN pairs found: \{fastmnn_result.mnn_pairs.length().to_string()}")
  println("   Variance explained by top 5 PCs:")
  let mut j = 0
  while j < 5 && j < fastmnn_result.var_explained.length() {
    println("     PC\{(j + 1).to_string()}: \{(fastmnn_result.var_explained[j] * 100.0).to_string()}%")
    j = j + 1
  }
  
  let fastmnn_score = @bio.compute_batch_mixing_score(fastmnn_result.corrected, fastmnn_result.batch_indices)
  println("   Batch mixing score after fastMNN: \{fastmnn_score.to_string()}")
  println("   Improvement: +\{((fastmnn_score - raw_score) / raw_score * 100.0).to_string()}%")
  println("")
  
  println("5. Comparing batch correction methods:")
  println("   Method               | Batch Mixing Score")
  println("   -------------------------------------------")
  println("   Raw data             | \{raw_score.to_string()}")
  println("   rescaleBatches       | \{rescale_score.to_string()}")
  println("   fastMNN              | \{fastmnn_score.to_string()}")
  println("")
  
  println("6. Finding mutual nearest neighbors between batches...")
  let mnn_pairs = @bio.find_mutual_nearest_neighbors(batches[0], batches[1], 5)
  println("   Found \{mnn_pairs.length().to_string()} MNN pairs between Batch_0 and Batch_1")
  if mnn_pairs.length() > 0 {
    println("   First 5 MNN pairs:")
    let mut k = 0
    while k < 5 && k < mnn_pairs.length() {
      let pair = mnn_pairs[k]
      println("     (Batch_0 cell \{pair.0.to_string()}, Batch_1 cell \{pair.1.to_string()})")
      k = k + 1
    }
  }
  println("")
  
  println("=== Demo completed successfully! ===")
  println("")
  println("Key takeaways:")
  println("- rescaleBatches: Simple and fast, good for global scaling differences")
  println("- fastMNN: More sophisticated, handles complex batch effects")
  println("- Batch mixing score quantifies how well batches are mixed after correction")
}
