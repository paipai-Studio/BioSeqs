///|
/// Codon alignment and dN/dS analysis demo inspired by Biopython Bio.codonalign.
///
/// Run: moon run examples/codon_align_demo/main.mbt
fn main {
  println("=== Bio.codonalign Demo ===")

  let seq1 = "ATGTCTGAAGTGGAA"
  let seq2 = "ATGTCTGAAGTGGAA"
  let seq3 = "ATGTCCGAAGTGGAA"

  println("\n1. Input sequences:")
  println("   Seq1: \{seq1}")
  println("   Seq2: \{seq2}")
  println("   Seq3: \{seq3}")

  println("\n2. Translate codons:")
  println("   ATG -> \{@bio.codon_to_aa("ATG")}")
  println("   TCT -> \{@bio.codon_to_aa("TCT")}")
  println("   GAA -> \{@bio.codon_to_aa("GAA")}")

  println("\n3. Classify codon substitution:")
  println("   TCT vs TCC: \{@bio.classify_codon_substitution("TCT", "TCC")}")
  println("   TCT vs TTT: \{@bio.classify_codon_substitution("TCT", "TTT")}")

  println("\n4. Count codon substitutions:")
  let counts = @bio.count_codon_substitutions(seq1, seq3)
  println("   Synonymous: \{counts.synonymous}")
  println("   Non-synonymous: \{counts.nonsynonymous}")
  println("   Total: \{counts.total}")

  println("\n5. dN/dS ratio (Nei-Gojobori):")
  let dnds = @bio.calc_dnds_nei_gojobori(seq1, seq3)
  println("   dN: \{dnds.dn}")
  println("   dS: \{dnds.ds}")
  println("   dN/dS: \{dnds.dnds}")

  println("\n6. Codon usage bias (RSCU):")
  let bias = @bio.codon_usage_bias(seq1)
  for kv in bias {
    println("   \{kv.0}: \{kv.1}")
  }

  println("\n7. Total codon sites:")
  let sites = @bio.total_codon_sites(seq1)
  println("   Synonymous sites: \{sites.0}")
  println("   Non-synonymous sites: \{sites.1}")

  println("\n8. Effective number of codons (ENC):")
  let enc = @bio.effective_number_of_codons(seq1)
  println("   ENC: \{enc}")

  println("\n9. Jukes-Cantor correction:")
  println("   p=0.1 -> d=\{@bio.jukes_cantor_correction(0.1)}")
  println("   p=0.3 -> d=\{@bio.jukes_cantor_correction(0.3)}")
}
