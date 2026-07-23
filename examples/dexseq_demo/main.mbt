///|
/// DEXSeq Demo - Demonstrates differential exon usage analysis.

fn main {
  println("=== DEXSeq Demo ===")
  
  let ds = @bio.create_example_dexseq_dataset()
  println("Number of genes: \{ds.gene_ids.length()}")
  println("Number of exons: \{ds.exon_counts.length()}")
  
  println("")
  println("Normalizing counts...")
  let normalized = @bio.dexseq_normalize_counts(ds)
  println("Normalized exons: \{normalized.exon_counts.length()}")
  
  println("")
  println("Testing for differential exon usage...")
  let results = @bio.dexseq_test_for_exon_usage(ds)
  println("Results: \{results.length()} exons tested")
  
  println("")
  println("Filtering results (log2FC >= 0.5, padj <= 0.1)...")
  let filtered = @bio.dexseq_filter_results(results, 0.5, 0.1)
  println("Significant exons: \{filtered.length()}")
}
