///|
/// Advanced Codon Analysis Demo - Comprehensive codon usage analysis example.

///|
fn main {
  println("=== Advanced Codon Analysis Demo ===\n")

  let dna_seq = "ATGGCCATGAAATGGTGAATGGCCGCCGCCAAAAAA"

  println("DNA sequence: " + dna_seq)
  println("Sequence length: " + dna_seq.length().to_string() + "\n")

  println("=== Codon Usage Table ===")
  let usage = @bio.codon_usage_from_sequence(dna_seq)
  println("Total codons: " + @bio.CodonUsageTable::total(usage).to_string())
  println("")

  let counts = @bio.CodonUsageTable::counts(usage)
  let keys = counts.keys().collect()
  let mut ci = 0
  while ci < keys.length() {
    let codon = keys[ci]
    let count = counts.get(codon).unwrap_or(0)
    let freq = @bio.CodonUsageTable::get_frequency(usage, codon)
    println("  " + codon + ": count=" + count.to_string() +
      ", freq=" + freq.to_string())
    ci = ci + 1
  }

  println("\n=== Amino Acid Translation ===")
  println("  ATG -> " + @bio.get_amino_acid("ATG") + " (start)")
  println("  GCC -> " + @bio.get_amino_acid("GCC"))
  println("  AAA -> " + @bio.get_amino_acid("AAA"))
  println("  TGG -> " + @bio.get_amino_acid("TGG"))
  println("  TGA -> " + @bio.get_amino_acid("TGA") + " (stop)")

  println("\n=== RSCU Analysis ===")
  let rscu = @bio.calculate_rscu(usage)
  let rscu_values = @bio.RSCUResult::values(rscu)
  let rscu_keys = rscu_values.keys().collect()
  let mut ri = 0
  while ri < rscu_keys.length() {
    let codon = rscu_keys[ri]
    let val = rscu_values.get(codon).unwrap_or(0.0)
    if val > 1.5 || val < 0.5 {
      println("  " + codon + ": " + val.to_string() + " (bias)")
    }
    ri = ri + 1
  }

  println("\n=== Codon Adaptation Index ===")
  let ref_seq = "ATGATGATGGCCGCCGCCAAAAAA"
  let ref_usage = @bio.codon_usage_from_sequence(ref_seq)
  let cai_val = @bio.calculate_cai(dna_seq, ref_usage)
  println("CAI value: " + cai_val.to_string())
  if cai_val > 0.8 {
    println("Interpretation: High codon adaptation")
  } else if cai_val > 0.5 {
    println("Interpretation: Moderate codon adaptation")
  } else {
    println("Interpretation: Low codon adaptation")
  }

  let cai_detailed = @bio.calculate_cai_detailed(dna_seq, ref_usage)
  println("Detailed CAI score: " + @bio.CAIResult::score(cai_detailed).to_string())

  println("\n=== Effective Number of Codons ===")
  let enc_val = @bio.calculate_enc(usage)
  println("ENC value: " + enc_val.to_string())
  if enc_val < 35.0 {
    println("Interpretation: High codon usage bias")
  } else if enc_val < 50.0 {
    println("Interpretation: Moderate codon usage bias")
  } else {
    println("Interpretation: Low codon usage bias")
  }

  println("\n=== GC Content Analysis ===")
  let gc3 = @bio.calculate_gc3(dna_seq)
  println("GC3 (third position): " + gc3.to_string())
  let gc1 = @bio.gc_first_fraction(dna_seq)
  println("GC first position: " + gc1.to_string())
  let gc2 = @bio.gc_second_fraction(dna_seq)
  println("GC second position: " + gc2.to_string())
  let at3 = @bio.at_third_fraction(dna_seq)
  println("AT third position: " + at3.to_string())

  println("\n=== Optimal Codons ===")
  let optimal = @bio.optimal_codons(usage)
  let mut oi = 0
  while oi < optimal.length() {
    println("  " + optimal[oi].0 + " -> " + optimal[oi].1)
    oi = oi + 1
  }

  println("\n=== Rare Codons ===")
  let rare = @bio.rare_codons(usage)
  if rare.length() > 0 {
    let mut rri = 0
    while rri < rare.length() {
      println("  " + rare[rri])
      rri = rri + 1
    }
  } else {
    println("  No rare codons detected")
  }

  println("\n=== Comprehensive Bias Analysis ===")
  let bias = @bio.codon_bias_analysis(dna_seq)
  println("  ENC: " + bias.enc.to_string())
  println("  CAI: " + bias.cai.to_string())
  println("  GC3: " + bias.gc3.to_string())
  println("  Bias Score: " + bias.bias_score.to_string())
  println("  Optimal Fraction: " + bias.optimal_fraction.to_string())
  println("  Rare Fraction: " + bias.rare_fraction.to_string())
  println("  Total Codons Analyzed: " + bias.total_codons_analyzed.to_string())

  println("\n=== Organism-Specific Reference Tables ===")
  let cai_high = @bio.calculate_cai_with_org(dna_seq, "high_gc")
  println("  CAI (high GC): " + cai_high.to_string())
  let cai_low = @bio.calculate_cai_with_org(dna_seq, "low_gc")
  println("  CAI (low GC): " + cai_low.to_string())
  let cai_gen = @bio.calculate_cai_with_org(dna_seq, "general")
  println("  CAI (general): " + cai_gen.to_string())

  println("\n=== Amino Acid Composition ===")
  let comp = @bio.get_amino_acid_composition(usage)
  let comp_keys = comp.keys().collect()
  let mut cpi = 0
  while cpi < comp_keys.length() {
    let aa = comp_keys[cpi]
    let fraction = comp.get(aa).unwrap_or(0.0)
    println("  " + aa + ": " + fraction.to_string())
    cpi = cpi + 1
  }

  println("\n=== Multiple Sequence Usage Table ===")
  let seqs = [
    "ATGGCCATGAAATGGTGA",
    "GCCGCCGCCAAAAAAATG",
    "ATGATGGCCAAAGCCGCC",
  ]
  let multi_usage = @bio.codon_usage_from_sequences(seqs)
  println("Total codons (multi): " +
    @bio.CodonUsageTable::total(multi_usage).to_string())

  println("\n=== Codon Usage Table String ===")
  let table_str = @bio.codon_usage_to_string(usage)
  println(table_str)

  println("=== Advanced Codon Analysis Done ===")
}