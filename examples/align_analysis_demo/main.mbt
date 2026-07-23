// Demo for Bio.Align.analysis module

fn main {
  println("=== Bio.Align.analysis Module Demo ===\n")
  
  // Create test sequences
  let seq1 = "ATGAAACCC"
  let seq2 = "ATGAAGCCC"
  let seq3 = "ATGTACCCC"
  
  println("Sequence 1: \{seq1}")
  println("Sequence 2: \{seq2}")
  println("Sequence 3: \{seq3}")
  
  // Amino acid translation
  println("\n=== Amino Acid Translation ===")
  println("ATG -> \{@bio.aln_aa_from_codon("ATG")}")
  println("AAA -> \{@bio.aln_aa_from_codon("AAA")}")
  println("AAG -> \{@bio.aln_aa_from_codon("AAG")}")
  println("TAC -> \{@bio.aln_aa_from_codon("TAC")}")
  println("CCC -> \{@bio.aln_aa_from_codon("CCC")}")
  
  // dn/ds analysis
  println("\n=== dn/ds Analysis (Nei-Gojobori Method) ===")
  
  // Compare seq1 and seq2 (synonymous change)
  let result_12 = @bio.aln_analyze_dn_ds(seq1, seq2)
  println("Comparison seq1 vs seq2 (synonymous):")
  println("  dn: \{result_12.aln_get_dn()}")
  println("  ds: \{result_12.aln_get_ds()}")
  println("  dn/ds ratio: \{result_12.aln_get_ratio()}")
  println("  Positive selection: \{result_12.aln_has_positive_selection()}")
  println("  Purifying selection: \{result_12.aln_has_purifying_selection()}")
  
  // Compare seq1 and seq3 (non-synonymous change)
  let result_13 = @bio.aln_analyze_dn_ds(seq1, seq3)
  println("\nComparison seq1 vs seq3 (non-synonymous):")
  println("  dn: \{result_13.aln_get_dn()}")
  println("  ds: \{result_13.aln_get_ds()}")
  println("  dn/ds ratio: \{result_13.aln_get_ratio()}")
  println("  Positive selection: \{result_13.aln_has_positive_selection()}")
  
  // Evolutionary distances
  println("\n=== Evolutionary Distances ===")
  
  let sequences = [seq1, seq2, seq3]
  
  // P-distance matrix
  println("\nP-distance matrix:")
  let n = sequences.length()
  let mut i = 0
  while i < n {
    let mut j = 0
    let row: Array[Double] = []
    while j < n {
      let dist = @bio.aln_p_distance(sequences[i], sequences[j])
      row.push(dist)
      j = j + 1
    }
    println("  Row \{i}: \{row}")
    i = i + 1
  }
  
  // Jukes-Cantor distance matrix
  println("\nJukes-Cantor distance matrix:")
  let jc_matrix = @bio.aln_jukes_cantor_matrix(sequences)
  let mut x = 0
  while x < n {
    println("  Row \{x}: \{jc_matrix[x]}")
    x = x + 1
  }
  
  // Kimura 2-parameter distance matrix
  println("\nKimura 2-parameter distance matrix:")
  let k2p_matrix = @bio.aln_kimura_2p_matrix(sequences)
  let mut y = 0
  while y < n {
    println("  Row \{y}: \{k2p_matrix[y]}")
    y = y + 1
  }
  
  // Count sites
  println("\n=== Site Counting ===")
  let sites = @bio.aln_count_sites(seq1)
  println("Synonymous sites: \{sites[0]}")
  println("Non-synonymous sites: \{sites[1]}")
  
  println("\nDemo completed successfully!")
}
