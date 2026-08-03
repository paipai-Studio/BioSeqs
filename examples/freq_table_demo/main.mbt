///|
/// FreqTable module demo - Frequency table for sequence analysis.
///
/// Run: moon run examples/freq_table_demo
fn main {
  println("=== FreqTable Module Demo ===")
  println("Frequency table for sequence alphabet analysis")
  println("Reference: Biopython Bio.SubsMat.FreqTable module")
  println("")

  // ------------------------------------------------------------------
  // 1. Read a count file and display frequencies
  // ------------------------------------------------------------------
  println("=== 1. Read Count File ===")
  let count_text = @bio.freq_table_sample_count_text()
  println("Input count file:")
  println(count_text)
  let ft = @bio.freq_table_read_count(count_text)
  println("Parsed frequency table:")
  println(ft.to_string())
  println("Total count: \{ft.total_count()}")
  println("Has counts: \{ft.has_counts()}")
  println("")

  // ------------------------------------------------------------------
  // 2. Read a frequency file
  // ------------------------------------------------------------------
  println("=== 2. Read Frequency File ===")
  let freq_text = @bio.freq_table_sample_freq_text()
  println("Input frequency file:")
  println(freq_text)
  let ft2 = @bio.freq_table_read_freq(freq_text)
  println("Parsed frequency table:")
  println(ft2.to_string())
  println("Has counts: \{ft2.has_counts()}")
  println("")

  // ------------------------------------------------------------------
  // 3. Build from a sequence
  // ------------------------------------------------------------------
  println("=== 3. Build From Sequence ===")
  let sequence = "ATGCGTACGTAGATGCGTACGTACGTACGTACGATCGATCG"
  let ft3 = @bio.freq_table_from_sequence(sequence, None)
  println("Sequence: \{sequence}")
  println("Length: \{sequence.length()}")
  println("Letter counts:")
  let letters = ft3.letters()
  for i in 0..<letters.length() {
    let l = letters[i]
    println("  \{l}: count=\{ft3.count_of(l)}, freq=\{ft3.frequency(l) * 100.0}%")
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Shannon entropy
  // ------------------------------------------------------------------
  println("=== 4. Shannon Entropy ===")
  let h = @bio.freq_table_entropy(ft)
  println("DNA sample table entropy: \{h} bits")
  let max_entropy = 2.0
  println("Maximum entropy (uniform 4-letter): \{max_entropy} bits")
  println("Relative entropy: \{h / max_entropy * 100.0}%")
  println("")

  // ------------------------------------------------------------------
  // 5. Compare two distributions with Jensen-Shannon distance
  // ------------------------------------------------------------------
  println("=== 5. Compare Distributions (JS Distance) ===")
  let ft_dna = @bio.freq_table_sample_dna()
  let ft_protein = @bio.freq_table_sample_protein()
  println("DNA table size: \{ft_dna.size()} letters")
  println("Protein table size: \{ft_protein.size()} letters")
  // JS distance only works on shared alphabet, so compare two DNA tables.
  let dict : Map[String, Double] = Map([], capacity=4)
  dict["A"] = 0.4
  dict["C"] = 0.1
  dict["G"] = 0.1
  dict["T"] = 0.4
  let ft_skewed = @bio.freq_table_from_freqs(dict)
  let js = @bio.freq_table_js_distance(ft_dna, ft_skewed)
  println("JS distance (sample vs skewed A/T-rich): \{js}")
  println("")

  // ------------------------------------------------------------------
  // 6. Normalize a frequency table
  // ------------------------------------------------------------------
  println("=== 6. Normalize ===")
  let dict2 : Map[String, Double] = Map([], capacity=3)
  dict2["A"] = 2.0
  dict2["C"] = 3.0
  dict2["G"] = 5.0
  let ft_unnorm = @bio.freq_table_from_freqs(dict2)
  let ft_norm = @bio.freq_table_normalize(ft_unnorm)
  println("Before normalization:")
  let ls1 = ft_unnorm.letters()
  for i in 0..<ls1.length() {
    println("  \{ls1[i]}: \{ft_unnorm.frequency(ls1[i])}")
  }
  println("After normalization:")
  let ls2 = ft_norm.letters()
  for i in 0..<ls2.length() {
    println("  \{ls2[i]}: \{ft_norm.frequency(ls2[i])}")
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Protein frequency table
  // ------------------------------------------------------------------
  println("=== 7. Protein Frequency Table ===")
  let ft_p = @bio.freq_table_sample_protein()
  println("Protein table: \{ft_p.size()} amino acids")
  println("Total count: \{ft_p.total_count()}")
  let h_p = @bio.freq_table_entropy(ft_p)
  println("Shannon entropy: \{h_p} bits")
  let max_p = @math.ln(20.0) / @math.ln(2.0)
  println("Maximum entropy (uniform 20-letter): \{max_p} bits")
  println("Relative entropy: \{h_p / max_p * 100.0}%")
  println("")

  // ------------------------------------------------------------------
  // 8. Relative entropy (KL divergence)
  // ------------------------------------------------------------------
  println("=== 8. KL Divergence ===")
  let kl = @bio.freq_table_relative_entropy(ft_dna, ft_skewed)
  println("KL(sample DNA || skewed A/T-rich): \{kl} bits")
  println("")

  println("=== Demo complete ===")
}
