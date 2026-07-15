///|
/// Sequence utilities demo.
fn main {
  println("=== Sequence Utilities Demo ===\n")

  let dna_seq = "ATCGATCGATCG"
  let protein_seq = "MALWMRLLPLLALLALWGPDPAAAFVN"

  println("1. Molecular weight:")
  let dna_weight = @src.sequtils_molecular_weight(dna_seq, "dna")
  let protein_weight = @src.sequtils_molecular_weight(protein_seq, "protein")
  println("   DNA weight: \{dna_weight}")
  println("   Protein weight: \{protein_weight}")

  println("\n2. GC content:")
  let gc = @src.gc_content(dna_seq)
  println("   GC content: \{gc}")

  println("\n3. Melting temperature:")
  let tm = @src.melting_temperature(dna_seq)
  println("   Tm: \{tm}")

  println("\n4. Sequence type detection:")
  println("   Is DNA: \{@src.is_dna(dna_seq)}")
  println("   Is RNA: \{@src.is_rna(dna_seq)}")
  println("   Is protein: \{@src.is_protein(protein_seq)}")

  println("\n5. Amino acid conversion:")
  let seq3 = @src.seq3("ARND")
  println("   seq3(ARND): \{seq3}")

  println("\n6. K-mer counting:")
  let kmers = @src.count_kmers(dna_seq, 2)
  println("   2-mer count: \{kmers.size()}")

  println("\n7. Shannon entropy:")
  let entropy = @src.sequtils_shannon_entropy(dna_seq)
  println("   Entropy: \{entropy}")

  println("\n8. Reverse complement:")
  let rc = @src.reverse_complement(dna_seq)
  println("   Original: \{dna_seq}")
  println("   Reverse complement: \{rc}")

  println("\n=== Demo Complete ===")
}
