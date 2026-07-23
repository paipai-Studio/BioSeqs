///|
/// EDASeq RNA-seq Exploratory Data Analysis Demo.
fn main {
  println("=== EDASeq RNA-seq Exploratory Data Analysis Demo ===\n")

  // Create simulated RNA-seq data
  println("1. Creating simulated RNA-seq data...")
  let gene_ids = ["Gene_A", "Gene_B", "Gene_C", "Gene_D", "Gene_E", "Gene_F"]
  let sample_ids = ["sample_1", "sample_2"]
  // Counts: genes (rows) x samples (columns)
  let counts : Array[Array[Double]] = [
    [100.0, 150.0],   // Gene_A
    [200.0, 180.0],   // Gene_B
    [500.0, 600.0],   // Gene_C
    [50.0, 80.0],     // Gene_D
    [300.0, 250.0],   // Gene_E
    [150.0, 200.0]    // Gene_F
  ]

  // Create gene annotations with GC content and length
  let annotations : Array[@src.EDASeqGeneAnno] = [
    @src.EDASeqGeneAnno::new("Gene_A", 0.45, 2000.0),
    @src.EDASeqGeneAnno::new("Gene_B", 0.65, 3500.0),
    @src.EDASeqGeneAnno::new("Gene_C", 0.55, 1200.0),
    @src.EDASeqGeneAnno::new("Gene_D", 0.72, 4500.0),
    @src.EDASeqGeneAnno::new("Gene_E", 0.38, 800.0),
    @src.EDASeqGeneAnno::new("Gene_F", 0.50, 2800.0)
  ]
  println("   Genes: \{gene_ids.length()}, Samples: \{sample_ids.length()}")

  // Create dataset
  println("\n2. Creating EDASeqDataSet...")
  let data = @src.EDASeqDataSet::new(gene_ids, sample_ids, counts, annotations)
  println("   Dataset created: \{data.eda_n_genes()} genes, \{data.eda_n_samples()} samples")

  // Configure parameters
  let params = @src.EDASeqParams::new()

  // Within-lane normalization
  println("\n3. Within-lane normalization (GC content)...")
  let within_result = @src.eda_normalize_within_lane(data, params)
  println("   Within-lane normalized counts computed")
  println("   GC fits converged: \{within_result.gc_fits.length()} samples processed")

  // Between-lane normalization
  println("\n4. Between-lane normalization (median ratio)...")
  let between_result = @src.eda_normalize_between_lane(data, params)
  println("   Size factors per sample:")
  let mut sf_idx = 0
  while sf_idx < between_result.size_factors.length() {
    let sf = between_result.size_factors[sf_idx].to_string()
    println("     \{sample_ids[sf_idx]}: \{sf}")
    sf_idx = sf_idx + 1
  }

  // Full normalization
  println("\n5. Running full normalization pipeline...")
  let full_result = @src.eda_run(data, params)
  println("   Final normalized counts computed")

  // Compare original vs normalized counts
  println("\n6. Original vs Normalized counts (Gene_A):")
  let mut s = 0
  while s < sample_ids.length() {
    let orig = counts[0][s].to_string()
    let norm = full_result.norm_counts[0][s].to_string()
    println("     \{sample_ids[s]}: original=\{orig}, normalized=\{norm}")
    s = s + 1
  }

  // RPKM calculation
  println("\n7. Computing RPKM values...")
  let rpkm = @src.eda_rpkm(counts, annotations, data.lib_sizes)
  println("   RPKM values for Gene_C:")
  let mut rpkm_s = 0
  while rpkm_s < sample_ids.length() {
    let val = rpkm[2][rpkm_s].to_string()
    println("     \{sample_ids[rpkm_s]}: \{val}")
    rpkm_s = rpkm_s + 1
  }

  println("\n=== EDASeq Demo Complete ===")
}
