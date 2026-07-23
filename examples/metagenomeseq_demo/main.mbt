///|
/// metagenomeSeq Demo - Demonstrates zero-inflated model for microbiome analysis.

fn main {
  println("=== metagenomeSeq Demo ===")
  
  let obj = @bio.create_example_mrexperiment()
  println("Number of taxa: \{obj.taxa_names.length()}")
  println("Number of samples: \{obj.sample_names.length()}")
  
  println("")
  println("Normalizing counts...")
  let normalized = @bio.mg_normalize_counts(obj)
  println("Normalized counts: \{normalized.counts.length()} taxa")
  
  println("")
  println("Calculating zero inflation...")
  let zero_probs = @bio.mg_calculate_zero_inflation(obj)
  println("Zero inflation probabilities:")
  let mut i = 0
  while i < obj.taxa_names.length() {
    println("  \{obj.taxa_names[i]}: \{zero_probs[i]}")
    i = i + 1
  }
  
  println("")
  println("Testing differential abundance...")
  let group = ["control", "control", "control", "treatment", "treatment", "treatment"]
  let results = @bio.mg_test_zero_inflated(obj, group)
  println("Results: \{results.length()} taxa tested")
}
