///|
/// Machine learning feature extraction example.
///
/// Demonstrates extracting features from DNA and protein sequences
/// for machine learning tasks.
///
/// Run: moon run examples/ml_features/main.mbt
fn main {
  let dna_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

  println("=== DNA Sequence Features ===")
  println("Sequence: \{dna_seq}")
  println("Length: \{dna_seq.length()}")

  let comp = @bio.nucleotide_composition(dna_seq)
  println("\nNucleotide Composition:")
  println("  A: \{comp.a}")
  println("  T: \{comp.t}")
  println("  G: \{comp.g}")
  println("  C: \{comp.c}")
  println("  GC content: \{comp.g + comp.c}")

  let kmer2 = @bio.kmer_frequency(dna_seq, 2, true)
  let kmer3 = @bio.kmer_frequency(dna_seq, 3, true)
  let kmer4 = @bio.kmer_frequency(dna_seq, 4, true)

  println("\nk-mer Frequency:")
  println("  Dinucleotides: \{kmer2.keys().count()} unique")
  println("  Trinucleotides: \{kmer3.keys().count()} unique")
  println("  Tetranucleotides: \{kmer4.keys().count()} unique")

  let dinuc_keys = kmer2.keys().collect()
  let dinuc_limit = min_int(dinuc_keys.length(), 5)
  for i in 0..<dinuc_limit {
    println("    \{dinuc_keys[i]}: \{kmer2[dinuc_keys[i]]}")
  }

  let canonical = @bio.dna_canonical_kmer_frequency(dna_seq, 3, true)
  println("\nCanonical k-mer Frequency (k=3):")
  println("  Unique canonical: \{canonical.keys().count()}")

  let dna_features = @bio.dna_feature_vector(dna_seq)
  let dna_names = @bio.dna_feature_names()
  println("\nDNA Feature Vector:")
  println("  Dimensions: \{dna_features.length()}")
  let dna_limit = min_int(dna_features.length(), 10)
  for i in 0..<dna_limit {
    println("  \{dna_names[i]}: \{dna_features[i]}")
  }

  let prot_seq = "MKKLLLISVLLFLSSAYSR"

  println("\n=== Protein Sequence Features ===")
  println("Sequence: \{prot_seq}")
  println("Length: \{prot_seq.length()}")

  let aa_comp = @bio.amino_acid_composition(prot_seq)
  println("\nAmino Acid Composition:")
  let aa_order : Array[UInt16] = [
    'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
    'T', 'W', 'Y', 'V',
  ]
  for i in 0..<aa_order.length() {
    let val = aa_comp.get(aa_order[i]).unwrap_or(0.0)
    if val > 0.0 {
      println("  \{aa_order[i].to_char().unwrap_or('?')}: \{val}")
    }
  }

  let dipep = @bio.dipeptide_composition(prot_seq)
  let tripep = @bio.tripeptide_composition(prot_seq)
  println("\nDipeptide Composition: \{dipep.keys().count()} unique")
  println("Tripeptide Composition: \{tripep.keys().count()} unique")

  let avg_h = @bio.avg_hydrophobicity(prot_seq)
  let avg_c = @bio.avg_charge(prot_seq)
  let avg_p = @bio.avg_polarity(prot_seq)
  let mw = @bio.molecular_weight(prot_seq)

  println("\nPhysical-Chemical Properties:")
  println("  Average Hydrophobicity: \{avg_h}")
  println("  Average Charge: \{avg_c}")
  println("  Average Polarity: \{avg_p}")
  println("  Molecular Weight: \{mw}")

  let (h, s, c) = @bio.secondary_structure_propensity(prot_seq)
  println("\nSecondary Structure Propensity (Chou-Fasman):")
  println("  Helix: \{h}")
  println("  Sheet: \{s}")
  println("  Coil: \{c}")

  let pos_feat = @bio.position_specific_features(prot_seq, 5)
  println("\nPosition-Specific Features (5 positions):")
  println("  Features: \{pos_feat.length()}")

  let prot_features = @bio.protein_feature_vector(prot_seq)
  let prot_names = @bio.protein_feature_names()
  println("\nProtein Feature Vector:")
  println("  Dimensions: \{prot_features.length()}")

  println("\nFirst 20 features (amino acid composition):")
  let prot_limit = min_int(20, prot_features.length())
  for i in 0..<prot_limit {
    println("  \{prot_names[i]}: \{prot_features[i]}")
  }

  println("\nPhysical-chemical features:")
  println("  \{prot_names[40]}: \{prot_features[40]}  (avg_hydrophobicity)")
  println("  \{prot_names[41]}: \{prot_features[41]}  (avg_charge)")
  println("  \{prot_names[42]}: \{prot_features[42]}  (avg_polarity)")
  println("  \{prot_names[43]}: \{prot_features[43]}  (molecular_weight)")
}

///|
fn min_int(a : Int, b : Int) -> Int {
  if a < b {
    a
  } else {
    b
  }
}
