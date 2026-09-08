///|
/// scMerge Demo - Single-cell batch correction using RUV4.
/// Demonstrates removal of unwanted variation (batch effects)
/// from synthetic single-cell expression data.

fn main {
  println("=== scMerge Demo ===")
  println("")

  // 1. Create synthetic single-cell data with batch effects
  println("--- Step 1: Create Synthetic Data ---")

  // Batch 1: 5 cells, 4 genes, base expression level ~1.0
  let batch1 = [
    [1.0, 2.0, 3.0, 4.0],
    [1.1, 2.1, 3.1, 4.1],
    [0.9, 1.9, 2.9, 3.9],
    [1.2, 2.2, 3.2, 4.2],
    [0.8, 1.8, 2.8, 3.8],
  ]

  // Batch 2: 5 cells, 4 genes, shifted by ~4.0 (batch effect)
  let batch2 = [
    [5.0, 6.0, 7.0, 8.0],
    [5.1, 6.1, 7.1, 8.1],
    [4.9, 5.9, 6.9, 7.9],
    [5.2, 6.2, 7.2, 8.2],
    [4.8, 5.8, 6.8, 7.8],
  ]

  // Batch 3: 5 cells, 4 genes, shifted by ~-2.0 (batch effect)
  let batch3 = [
    [-1.0, 0.0, 1.0, 2.0],
    [-0.9, 0.1, 1.1, 2.1],
    [-1.1, -0.1, 0.9, 1.9],
    [-0.8, 0.2, 1.2, 2.2],
    [-1.2, -0.2, 0.8, 1.8],
  ]

  let batches = [batch1, batch2, batch3]
  let n_batches = batches.length()
  let n_genes = batches[0][0].length()

  println("Batches: \{n_batches.to_string()}")
  println("Genes per cell: \{n_genes.to_string()}")

  // Show batch means before correction
  println("")
  println("Batch means before correction (gene 0):")
  let mut b = 0
  while b < n_batches {
    let batch = batches[b]
    let mut sum = 0.0
    let mut c = 0
    while c < batch.length() {
      sum = sum + batch[c][0]
      c = c + 1
    }
    let mean = sum / batch.length().to_double()
    println("  Batch \{b.to_string()}: mean = \{mean.to_string()}")
    b = b + 1
  }

  // 2. Find stably expressed genes (SEGs)
  println("")
  println("--- Step 2: Find Stably Expressed Genes (SEGs) ---")

  let seg = @src.scm_find_seg(batches, n_seg=2)
  println("SEG indices: \{seg.to_string()}")
  // Interpret: which genes are most stable across batches
  let mut i = 0
  while i < seg.length() {
    println("  Gene \{seg[i].to_string()} is stably expressed")
    i = i + 1
  }

  // 3. Run scMerge batch correction
  println("")
  println("--- Step 3: Run scMerge Batch Correction ---")

  let result = @src.scmerge(batches, n_seg=2, n_factors=2)

  println("Corrected data dimensions: \{result.corrected_data.length().to_string()} cells x \{result.corrected_data[0].length().to_string()} genes")

  // 4. Show corrected batch means
  println("")
  println("--- Step 4: Batch Means After Correction ---")

  let mut offset = 0
  b = 0
  while b < n_batches {
    let batch_size = batches[b].length()
    let mut gene0_sum = 0.0
    let mut c = 0
    while c < batch_size {
      gene0_sum = gene0_sum + result.corrected_data[offset + c][0]
      c = c + 1
    }
    let mean = gene0_sum / batch_size.to_double()
    println("  Batch \{b.to_string()}: corrected mean (gene 0) = \{mean.to_string()}")
    offset = offset + batch_size
    b = b + 1
  }

  // 5. Show cell type labels
  println("")
  println("--- Step 5: Cell Type Labels ---")

  println("Cell type labels: \{result.cell_type_labels.to_string()}")

  // 6. Show estimated batch factors
  println("")
  println("--- Step 6: Estimated Batch Factors ---")

  println("First 3 cells' batch factors:")
  let mut j = 0
  while j < 3 {
    println(
      "  Cell \{j.to_string()}: factors = \{result.estimated_batch_factors[j].to_string()}",
    )
    j = j + 1
  }

  // 7. Print sample values from corrected data
  println("")
  println("--- Step 7: Sample Corrected Values ---")

  println("First 3 cells, first 2 genes:")
  let mut k = 0
  while k < 3 {
    println(
      "  Cell \{k.to_string()}: gene 0 = \{result.corrected_data[k][0].to_string()}, gene 1 = \{result.corrected_data[k][1].to_string()}",
    )
    k = k + 1
  }

  println("")
  println("Summary:")
  println("  Total cells corrected: \{result.corrected_data.length().to_string()}")
  println("  SEGs used: \{result.seg_genes.length().to_string()}")
  println("  Batch correction completed successfully")

  println("")
  println("Done!")
}