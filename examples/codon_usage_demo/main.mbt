///|
/// Codon Usage Analysis Demo - DNA sequence codon usage analysis example.

///|
fn main {
  println("=== Codon Usage Analysis Demo ===\n")

  let dna_seq = "ATGGCCATGAAATGGTGAATGGCCGCCGCCAAATGGG"

  println("DNA sequence: " + dna_seq)
  println("Sequence length: " + dna_seq.length().to_string() + "\n")

  println("=== Count Codons ===")
  let counts = @bio.count_codons(dna_seq)
  println("Total unique codons: " + counts.length().to_string())

  println("\n=== RSCU Analysis ===")
  let rscu_map = @bio.codon_rscu(counts.copy())
  let top_codons = ["ATG", "GCC", "AAA", "TGG"]
  let mut i = 0
  while i < top_codons.length() {
    let codon = top_codons[i]
    match rscu_map.get(codon) {
      Some(val) => println("  " + codon + ": " + val.to_string())
      None => println("  " + codon + ": not found")
    }
    i = i + 1
  }

  println("\n=== Optimal Codons ===")
  let optimal = @bio.codon_optimal(rscu_map.copy())
  println("Number of optimal codons: " + optimal.length().to_string())

  println("\n=== CAI (Codon Adaptation Index) ===")
  let cai_val = @bio.codon_cai(dna_seq, optimal.copy())
  println("CAI value: " + cai_val.to_string())
  if cai_val > 0.8 {
    println("Interpretation: High codon adaptation")
  } else if cai_val > 0.5 {
    println("Interpretation: Moderate codon adaptation")
  } else {
    println("Interpretation: Low codon adaptation")
  }

  println("\n=== ENC (Effective Number of Codons) ===")
  let enc_val = @bio.codon_enc(counts.copy())
  println("ENC value: " + enc_val.to_string())
  if enc_val < 35.0 {
    println("Interpretation: High codon usage bias")
  } else if enc_val < 50.0 {
    println("Interpretation: Moderate codon usage bias")
  } else {
    println("Interpretation: Low codon usage bias")
  }

  println("\n=== GC3 (GC content at third position) ===")
  let gc3_val = @bio.codon_gc3(dna_seq)
  println("GC3 value: " + gc3_val.to_string())

  println("\n=== CBI (Codon Bias Index) ===")
  let opt_array = ["ATG", "GCC", "AAA", "TGG"]
  let cbi_val = @bio.codon_cbi(dna_seq, opt_array)
  println("CBI value: " + cbi_val.to_string())

  println("\n=== Fop (Frequency of Optimal Codons) ===")
  let fop_val = @bio.codon_fop(dna_seq, opt_array)
  println("Fop value: " + fop_val.to_string())

  println("\n=== Codon Translation ===")
  println("  ATG -> " + @bio.codon_translate("ATG"))
  println("  GCC -> " + @bio.codon_translate("GCC"))
  println("  TAA -> " + @bio.codon_translate("TAA"))

  println("\n=== Reverse Complement ===")
  println("  ATG -> " + @bio.codon_reverse_complement("ATG"))
  println("  GCA -> " + @bio.codon_reverse_complement("GCA"))

  println("\n=== Complete Analysis ===")
  let result = @bio.codon_usage_analysis(dna_seq)
  println("CAI: " + result.cai.to_string())
  println("ENC: " + result.enc.to_string())
  println("GC3: " + result.gc3.to_string())
  println("Optimal codons count: " + result.optimal_codons.length().to_string())

  println("\n=== Analysis with Reference Sequence ===")
  let ref_seq = "ATGATGATGGCCGCCGCCAAA"
  let ref_result = @bio.codon_usage_with_ref(dna_seq, ref_seq)
  println("CAI (with reference): " + ref_result.cai.to_string())

  println("\n=== Codon Usage Analysis Done ===")
}
