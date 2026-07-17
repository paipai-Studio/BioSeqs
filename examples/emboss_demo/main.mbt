///|
/// EMBOSS tools example.
///
/// Run: moon run examples/emboss_demo/main.mbt
fn main {
  println("=== EMBOSS Tools Demo ===")
  println("")

  let dna_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
  let prot_seq = "MKKLLLISVLLFLSSAYSR"

  println("=== Sequence Analysis ===")
  println("DNA: \{dna_seq}")
  println("Protein: \{prot_seq}")
  println("")

  println("=== GC/AT Skew ===")
  let gc_skew_val = @bio.emboss_gc_skew(dna_seq)
  println("GC Skew: \{gc_skew_val}")
  let at_skew_val = @bio.emboss_at_skew(dna_seq)
  println("AT Skew: \{at_skew_val}")
  println("")

  println("=== Molecular Weight ===")
  let mw = @bio.dna_molecular_weight(dna_seq)
  println("DNA Molecular Weight: \{mw}")
  println("")

  println("=== Melting Temperature ===")
  let tm_s = @bio.tm_simple("ATGC")
  println("Tm (simple): \{tm_s}")
  println("")

  println("=== Reverse Complement ===")
  let rc = @bio.emboss_reverse_complement("ATGC")
  println("ATGC -> \{rc}")
  println("")

  println("=== Translation ===")
  let protein = @bio.translate("ATG")
  println("ATG -> \{protein}")
  println("")

  println("=== ORF Finding ===")
  let orfs = @bio.find_orfs("ATGAAATAA", min_length=3)
  println("ORFs found: \{orfs.length()}")
  println("")

  println("=== Distance Calculation ===")
  let hd = @bio.hamming_distance("ATGC", "ATGG")
  println("Hamming distance (ATGC vs ATGG): \{hd}")
  let ld = @bio.levenshtein_distance("ATGC", "AGC")
  println("Levenshtein distance (ATGC vs AGC): \{ld}")
  println("")

  println("=== k-mer Counting ===")
  let kmers = @bio.emboss_count_kmers("ATGCAT", 2)
  println("2-mers found: \{kmers.length()}")
  println("")

  println("=== Protein Analysis ===")
  let pi = @bio.isoelectric_point("K")
  println("Isoelectric point (K): \{pi}")
  let ai = @bio.aliphatic_index(prot_seq)
  println("Aliphatic index: \{ai}")
  let gravy = @bio.gravy_score(prot_seq)
  println("GRAVY score: \{gravy}")
  let ii = @bio.instability_index(prot_seq)
  println("Instability index: \{ii}")
  println("")

  println("=== Amino Acid Composition ===")
  let comp = @bio.emboss_amino_acid_composition(prot_seq)
  println("Composition entries: \{comp.length()}")
  println("")

  println("=== Done ===")
}
