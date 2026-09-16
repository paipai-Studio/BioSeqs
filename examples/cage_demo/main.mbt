///|
/// Demo for cage module: CAGE data TSS analysis.

fn main {
  println("=== CAGE Demo: TSS Identification & Clustering ===")
  println("")

  // Example 1: Create CTSS data
  println("--- Step 1: CTSS Data ---")
  let ctss = [
    @src.CageCtss::new("chr1", 1000, "+", 50),
    @src.CageCtss::new("chr1", 1001, "+", 30),
    @src.CageCtss::new("chr1", 1002, "+", 20),
    @src.CageCtss::new("chr1", 1005, "+", 15),
    @src.CageCtss::new("chr1", 1010, "+", 10),
    @src.CageCtss::new("chr1", 2000, "-", 40),
    @src.CageCtss::new("chr1", 2001, "-", 35),
  ]
  let sample = @src.CageSample::new("sample1", ctss)
  println("Sample: " + sample.sample_id + ", total tags=" + sample.total_tags.to_string())
  println("")

  // Example 2: Create CAGEset
  println("--- Step 2: CAGEset ---")
  let cs0 = @src.cage_new("my_cage")
  let cs = @src.cage_add_sample(cs0, sample)
  println("CAGEset: " + cs.name + ", samples=" + @src.cage_n_samples(cs).to_string())
  println("")

  // Example 3: TPM normalization
  println("--- Step 3: TPM Normalization ---")
  let norm = @src.cage_normalize_tpm(cs)
  println("Normalization: " + norm.normalization)
  let s1 = norm.samples[0]
  println("  CTSS 1000+: normalized=" + s1.ctss[0].normalized.to_string())
  println("")

  // Example 4: CTSS clustering
  println("--- Step 4: CTSS Clustering ---")
  let clusters = @src.cage_cluster(norm, max_dist=10, min_tags=2)
  println("Clusters found: " + clusters.length().to_string())
  for i in 0..<clusters.length() {
    let cl = clusters[i]
    let shape = @src.cage_classify_tss_shape(cl)
    println("  Cluster " + i.to_string() + ": " + cl.chr + ":" + cl.start.to_string() + "-" + cl.end.to_string() +
            " tags=" + cl.total_tags.to_string() + " shape=" + shape +
            " score=" + cl.shape_score.to_string())
  }
  println("")

  // Example 5: Cluster summary
  println("--- Step 5: Cluster Summary ---")
  let (n, total_tags, avg_tags, avg_ctss, avg_shape) = @src.cage_cluster_summary(clusters)
  println("Total clusters: " + n.to_string())
  println("Total tags: " + total_tags.to_string())
  println("Avg tags/cluster: " + avg_tags.to_string())
  println("Avg shape score: " + avg_shape.to_string())
  println("")

  // Example 6: Power-law normalization
  println("--- Step 6: Power-law Normalization ---")
  let pl = @src.cage_normalize_powerlaw(cs, alpha=1.2)
  println("Normalization: " + pl.normalization)
  println("")

  // Example 7: Differential TSS usage
  println("--- Step 7: Differential TSS Usage ---")
  let ctss2 = [
    @src.CageCtss::new("chr1", 1000, "+", 100),
    @src.CageCtss::new("chr1", 1001, "+", 60),
    @src.CageCtss::new("chr1", 2000, "-", 20),
  ]
  let sample2 = @src.CageSample::new("sample2", ctss2)
  let diff = @src.cage_diff_tss(sample, sample2)
  println("Differential TSS: " + diff.length().to_string() + " positions")
  if diff.length() > 0 {
    let (chr, pos, strand, log2fc, total) = diff[0]
    println("  " + chr + ":" + pos.to_string() + " " + strand + " log2FC=" + log2fc.to_string())
  }
  println("")

  // Example 8: TSS annotation
  println("--- Step 8: TSS Annotation ---")
  let genes = [
    ("GENE1", 1000, "+", "chr1"),
    ("GENE2", 2000, "-", "chr1"),
  ]
  let ann = @src.cage_annotate_tss(clusters, genes)
  println("Annotations: " + ann.length().to_string())
  for i in 0..<ann.length() {
    let a = ann[i]
    println("  TSS " + a.tss_pos.to_string() + " -> " + a.gene_symbol + " (" + a.feature + ")")
  }
  println("")

  println("=== CAGE Demo Complete ===")
}
