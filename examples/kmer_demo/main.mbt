///|
/// Kmer Demo - Demonstrates k-mer counting and analysis.

fn main {
  println("=== Kmer Demo ===")
  
  let counter = @bio.create_example_kmer_counter()
  println("K-mer size: \{counter.k}")
  println("Total kmers: \{counter.total_kmers()}")
  println("Unique kmers: \{counter.unique_kmers()}")
  
  let top = counter.get_top_kmers(3)
  println("")
  println("Top 3 kmers:")
  let mut i = 0
  while i < top.length() {
    println("  \{top[i].0}: \{top[i].1}")
    i = i + 1
  }
  
  let spectrum = @bio.kmer_spectrum("ATCGATCG", 2, 4)
  println("")
  println("K-mer spectrum (k=2 to 4):")
  let mut k = 2
  let mut j = 0
  while j < spectrum.length() {
    println("  k=\{k}: \{spectrum[j]} unique kmers")
    k = k + 1
    j = j + 1
  }
}