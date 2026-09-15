///|
/// primersearch demo - EMBOSS primer matching and PCR product prediction.
/// Inspired by Biopython Bio.Emboss.PrimerSearch.

fn main {
  println("=== primersearch Demo: Primer Matching & PCR ===")
  println("")

  // Step 1: Define primer pairs
  println("--- Step 1: Primer Pairs ---")
  let pairs = [
    @src.PrimerPair::new("GeneA", "ATGAAACGTT", "TTACAGCTAA"),
    @src.PrimerPair::new("GeneB", "GGCCTTAACC", "GGAATTCCGG"),
  ]
  let mut i = 0
  while i < pairs.length() {
    println(pairs[i].name + ": " + pairs[i].forward + " / " + pairs[i].reverse)
    i = i + 1
  }
  println("")

  // Step 2: Target sequences
  println("--- Step 2: Target Sequences ---")
  let seqs : Array[(String, String)] = [
    ("chr1", "GGGATGAAACGTTCCCTTTACAGCTAAAGGG"),
    ("chr2", "TTTGGCCGGCCTTAACCCGGGGAATTCCGGTTT"),
  ]
  let mut j = 0
  while j < seqs.length() {
    println(seqs[j].0 + ": " + seqs[j].1)
    j = j + 1
  }
  println("")

  // Step 3: Run primersearch
  println("--- Step 3: Primersearch Results ---")
  let result = @src.primersearch_run(pairs, seqs)
  println(@src.primersearch_to_string(result))
  println("")

  // Step 4: Filter by size
  println("--- Step 4: Filter by Amplicon Size ---")
  if result.products.length() > 0 {
    let size = result.products[0].product_size
    let filtered = @src.primersearch_filter_by_size(result, size - 5, size + 5)
    println("Products within " + (size - 5).to_string() + "-" + (size + 5).to_string() + " bp: " + filtered.length().to_string())
  }
  println("")

  // Step 5: Mismatch-tolerant search
  println("--- Step 5: Mismatch-Tolerant Search ---")
  let seq = "GGGATGAAACGTTCCCTTTACAGCTAAAGGG"
  let hits0 = @src.primersearch_find_hits(seq, "ATGAAACGTA", max_mismatches=0)
  let hits1 = @src.primersearch_find_hits(seq, "ATGAAACGTA", max_mismatches=1)
  println("Primer 'ATGAAACGTA' (1 mismatch):")
  println("  0 mismatches: " + hits0.length().to_string() + " hits")
  println("  1 mismatch:   " + hits1.length().to_string() + " hits")
  println("")

  println("=== primersearch Demo Complete ===")
}