///|
test "iupac_ambiguous_dna_map" {
  let m = @src.iupac_ambiguous_dna_map()
  assert_eq(m.size(), 11)
  let r = m.get("R").unwrap()
  assert_eq(r, ["G", "A"])
  let y = m.get("Y").unwrap()
  assert_eq(y, ["C", "T"])
  let n = m.get("N").unwrap()
  assert_eq(n.length(), 4)
}

///|
test "iupac_ambiguous_rna_map" {
  let m = @src.iupac_ambiguous_rna_map()
  assert_eq(m.size(), 11)
  let u = m.get("Y").unwrap()
  assert_eq(u, ["C", "U"])
}

///|
test "is_ambiguous_dna" {
  assert_true(@src.is_ambiguous_dna("R"))
  assert_true(@src.is_ambiguous_dna("N"))
  assert_true(!@src.is_ambiguous_dna("A"))
  assert_true(!@src.is_ambiguous_dna("T"))
  assert_true(@src.is_ambiguous_dna("r")) // case insensitive
}

///|
test "is_ambiguous_rna" {
  assert_true(@src.is_ambiguous_rna("Y"))
  assert_true(!@src.is_ambiguous_rna("A"))
}

///|
test "expand_ambiguous_dna" {
  let r = @src.expand_ambiguous_dna("R")
  assert_eq(r, ["G", "A"])
  let a = @src.expand_ambiguous_dna("A")
  assert_eq(a, ["A"])
  let n = @src.expand_ambiguous_dna("N")
  assert_eq(n.length(), 4)
}

///|
test "expand_ambiguous_rna" {
  let y = @src.expand_ambiguous_rna("Y")
  assert_eq(y, ["C", "U"])
  let a = @src.expand_ambiguous_rna("A")
  assert_eq(a, ["A"])
}

///|
test "amino_3_to_1" {
  let m = @src.amino_3_to_1()
  assert_eq(m.get("Ala").unwrap(), "A")
  assert_eq(m.get("Ter").unwrap(), "*")
  assert_eq(m.size(), 21) // 20 standard + Ter
}

///|
test "amino_1_to_3" {
  let m = @src.amino_1_to_3()
  assert_eq(m.get("A").unwrap(), "Ala")
  assert_eq(m.get("*").unwrap(), "Ter")
  assert_eq(m.size(), 21)
}

///|
test "amino3_to_1" {
  assert_eq(@src.amino3_to_1("Ala"), "A")
  assert_eq(@src.amino3_to_1("Ter"), "*")
  assert_eq(@src.amino3_to_1("XXX"), "?")
}

///|
test "amino1_to_3" {
  assert_eq(@src.amino1_to_3("A"), "Ala")
  assert_eq(@src.amino1_to_3("*"), "Ter")
  assert_eq(@src.amino1_to_3("X"), "???")
}

///|
test "amino_acid_properties" {
  let ala = @src.amino_acid_properties("A")
  assert_eq(ala.amino_acid, "A")
  assert_eq(ala.property, "hydrophobic")
  assert_eq(ala.polarity, "nonpolar")
  assert_true((ala.hydropathy - 1.8).abs() < 0.01)

  let arg = @src.amino_acid_properties("R")
  assert_eq(arg.charge, 1.0)
  assert_eq(arg.property, "positive")

  let unknown = @src.amino_acid_properties("X")
  assert_eq(unknown.property, "unknown")
}

///|
test "standard_codon_table" {
  let t = @src.standard_codon_table()
  assert_eq(t.size(), 64)
  assert_eq(t.get("ATG").unwrap(), "M") // Start codon
  assert_eq(t.get("TAA").unwrap(), "*") // Stop
  assert_eq(t.get("TAG").unwrap(), "*") // Stop
  assert_eq(t.get("TGA").unwrap(), "*") // Stop
}

///|
test "codon_lookup" {
  assert_eq(@src.codon_lookup("ATG"), "M")
  assert_eq(@src.codon_lookup("TAA"), "*")
  assert_eq(@src.codon_lookup("atg"), "M") // case insensitive
  assert_eq(@src.codon_lookup("XXX"), "?")
}

///|
test "reverse_codon_table" {
  let rev = @src.reverse_codon_table()
  assert_eq(rev.get("M").unwrap(), ["ATG"]) // Met has only 1 codon
  let leu = rev.get("L").unwrap()
  assert_eq(leu.length(), 6) // Leu has 6 codons
}

///|
test "num_synonymous_codons" {
  assert_eq(@src.num_synonymous_codons("M"), 1)
  assert_eq(@src.num_synonymous_codons("L"), 6)
  assert_eq(@src.num_synonymous_codons("*"), 3)
  assert_eq(@src.num_synonymous_codons("X"), 0)
}

///|
test "reverse_complement_iupac" {
  assert_eq(@src.reverse_complement_iupac("ATCG"), "CGAT")
  assert_eq(@src.reverse_complement_iupac("AATTCC"), "GGAATT")
  assert_eq(@src.reverse_complement_iupac("R"), "Y")
  assert_eq(@src.reverse_complement_iupac("Y"), "R")
  assert_eq(@src.reverse_complement_iupac("N"), "N")
  assert_eq(@src.reverse_complement_iupac("K"), "M")
  assert_eq(@src.reverse_complement_iupac("M"), "K")
  assert_eq(@src.reverse_complement_iupac("B"), "V")
  assert_eq(@src.reverse_complement_iupac("D"), "H")
  assert_eq(@src.reverse_complement_iupac("H"), "D")
  assert_eq(@src.reverse_complement_iupac("S"), "S")
  assert_eq(@src.reverse_complement_iupac("W"), "W")
}

///|
test "reverse_complement_iupac_with_uracil" {
  assert_eq(@src.reverse_complement_iupac("AUCG"), "CGAT") // U -> A
}
