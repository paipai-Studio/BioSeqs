// Demo for Bio.FreqAnalysis module

fn main {
  println("=== Bio.FreqAnalysis Module Demo ===\n")
  
  // Nucleotide frequency
  let seq = "ATGCGTACGTTAGCTA"
  println("Sequence: \{seq}")
  
  let nuc_freq = @bio.fa_nucleotide_frequency(seq)
  println("\nNucleotide frequencies:")
  println("  A: \{nuc_freq.fa_get_frequency("A")}")
  println("  T: \{nuc_freq.fa_get_frequency("T")}")
  println("  G: \{nuc_freq.fa_get_frequency("G")}")
  println("  C: \{nuc_freq.fa_get_frequency("C")}")
  println("  Total count: \{nuc_freq.fa_get_total_count()}")
  
  // GC content
  let gc = @bio.fa_gc_content(seq)
  println("\nGC content: \{gc}")
  
  // AT content
  let at = @bio.fa_at_content(seq)
  println("AT content: \{at}")
  
  // K-mer counting
  println("\n=== K-mer Analysis ===")
  let di_freq = @bio.fa_dinucleotide_frequency(seq)
  let patterns = di_freq.fa_get_patterns()
  println("Dinucleotide frequencies:")
  let mut i = 0
  while i < patterns.length() {
    let p = patterns[i]
    println("  \{p}: \{di_freq.fa_get_frequency(p)}")
    i = i + 1
  }
  
  // Codon usage
  println("\n=== Codon Usage ===")
  let coding_seq = "ATGAAACCCGGGTTTAAA"
  println("Coding sequence: \{coding_seq}")
  let codon_freq = @bio.fa_codon_usage(coding_seq)
  let codon_patterns = codon_freq.fa_get_patterns()
  let mut j = 0
  while j < codon_patterns.length() {
    let c = codon_patterns[j]
    let count = codon_freq.fa_get_count(c)
    let freq = codon_freq.fa_get_frequency(c)
    println("  \{c}: count=\{count}, freq=\{freq}")
    j = j + 1
  }
  
  // Sequence complexity
  println("\n=== Sequence Complexity ===")
  let entropy = @bio.fa_sequence_complexity(seq)
  println("Shannon entropy: \{entropy}")
  
  // Find motifs
  println("\n=== Motif Finding ===")
  let positions = @bio.fa_find_motif(seq, "TAC")
  if positions.length() > 0 {
    println("Found 'TAC' at positions:")
    let mut k = 0
    while k < positions.length() {
      println("  \{positions[k]}")
      k = k + 1
    }
  } else {
    println("'TAC' not found in sequence")
  }
  
  println("\nDemo completed successfully!")
}
