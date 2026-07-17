///|
/// Sequence utilities demo.
///
/// Run: moon run examples/seq_utils_demo/main.mbt
fn main {
  println("=== Sequence Utilities Demo ===\n")

  let dna_seq = "ATCGATCGATCG"
  let protein_seq = "MALWMRLLPLLALLALWGPDPAAAFVN"

  println("1. Molecular weight:")
  let protein = @bio.ProteinAnalysis::new(protein_seq)
  let protein_weight = protein.molecular_weight()
  println("   Protein weight: \{protein_weight}")

  println("\n2. GC content:")
  let gc = @bio.gc(dna_seq)
  println("   GC content: \{gc}%")

  println("\n3. Melting temperature (Wallace):")
  let tm = @bio.tm_wallace(dna_seq)
  println("   Tm: \{tm}")

  println("\n4. Protein aromaticity:")
  let aromaticity = protein.aromaticity()
  println("   Aromaticity: \{aromaticity}")

  println("\n5. Protein GRAVY score:")
  let gravy = protein.gravy()
  println("   GRAVY: \{gravy}")

  println("\n6. Amino acid composition:")
  let aa_counts = protein.count_amino_acids()
  println("   Amino acids: \{aa_counts.size()}")

  println("\n7. Amino acid percentage:")
  let aa_percent = protein.get_amino_acids_percent()
  println("   Amino acid percent: \{aa_percent.size()}")

  println("\n8. Reverse complement:")
  let dna = @bio.Seq::new(dna_seq)
  let rc = dna.reverse_complement()
  println("   Original: \{dna.data}")
  println("   Reverse complement: \{rc.data}")

  println("\n9. Complement:")
  let comp = dna.complement()
  println("   Complement: \{comp.data}")

  println("\n10. Transcribe:")
  let rna = dna.transcribe()
  println("   RNA: \{rna.data}")

  println("\n=== Demo Complete ===")
}
