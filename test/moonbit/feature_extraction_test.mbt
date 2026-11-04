///|
test "kmer_frequency" {
  let seq = "ATGC"
  let freq = @src.kmer_frequency(seq, 2, true)

  assert_eq(freq.get("AT").unwrap_or(0.0), 0.3333333333333333)
  assert_eq(freq.get("TG").unwrap_or(0.0), 0.3333333333333333)
  assert_eq(freq.get("GC").unwrap_or(0.0), 0.3333333333333333)
  assert_eq(freq.keys().count(), 3)
}

///|
test "kmer_frequency_non_normalized" {
  let seq = "AAAAA"
  let freq = @src.kmer_frequency(seq, 2, false)

  assert_eq(freq.get("AA").unwrap_or(0.0), 4.0)
}

///|
test "kmer_frequency_empty" {
  let freq = @src.kmer_frequency("", 2, true)
  assert_eq(freq.keys().count(), 0)
}

///|
test "dna_canonical_kmer_frequency" {
  let seq = "ATGC"
  let freq = @src.dna_canonical_kmer_frequency(seq, 2, true)

  assert_eq(freq.get("AT").unwrap_or(0.0), 0.3333333333333333)
}

///|
test "nucleotide_composition" {
  let seq = "ATGC"
  let comp = @src.nucleotide_composition(seq)

  assert_eq(comp.a, 0.25)
  assert_eq(comp.t, 0.25)
  assert_eq(comp.g, 0.25)
  assert_eq(comp.c, 0.25)
}

///|
test "nucleotide_composition_empty" {
  let comp = @src.nucleotide_composition("")

  assert_eq(comp.a, 0.0)
  assert_eq(comp.t, 0.0)
  assert_eq(comp.g, 0.0)
  assert_eq(comp.c, 0.0)
}

///|
test "amino_acid_composition" {
  let seq = "MRK"
  let comp = @src.amino_acid_composition(seq)

  assert_eq(comp.get('M').unwrap_or(0.0), 0.3333333333333333)
  assert_eq(comp.get('R').unwrap_or(0.0), 0.3333333333333333)
  assert_eq(comp.get('K').unwrap_or(0.0), 0.3333333333333333)
}

///|
test "dipeptide_composition" {
  let seq = "MRK"
  let comp = @src.dipeptide_composition(seq)

  assert_eq(comp.get("MR").unwrap_or(0.0), 0.5)
  assert_eq(comp.get("RK").unwrap_or(0.0), 0.5)
}

///|
test "tripeptide_composition" {
  let seq = "MRKK"
  let comp = @src.tripeptide_composition(seq)

  assert_eq(comp.get("MRK").unwrap_or(0.0), 0.5)
  assert_eq(comp.get("RKK").unwrap_or(0.0), 0.5)
}

///|
test "avg_hydrophobicity" {
  let seq = "ILV"
  let h = @src.avg_hydrophobicity(seq)

  assert_true(h > 0.0)
}

///|
test "avg_charge" {
  let seq = "RK"
  let c = @src.avg_charge(seq)

  assert_eq(c, 1.0)
}

///|
test "avg_polarity" {
  let seq = "RK"
  let p = @src.avg_polarity(seq)

  assert_eq(p, 1.0)
}

///|
test "molecular_weight" {
  let seq = "A"
  let mw = @src.molecular_weight(seq)

  assert_true(mw > 80.0)
  assert_true(mw < 100.0)
}

///|
test "secondary_structure_propensity" {
  let seq = "AE"
  let ss = @src.secondary_structure_propensity(seq)

  assert_true(ss.0 > 0.0)
  assert_true(ss.1 > 0.0)
}

///|
test "position_specific_features" {
  let seq = "ILV"
  let features = @src.position_specific_features(seq, 2)

  assert_eq(features.length(), 6)
}

///|
test "protein_feature_vector" {
  let seq = "ILV"
  let features = @src.protein_feature_vector(seq)

  assert_true(features.length() > 0)
}

///|
test "dna_feature_vector" {
  let seq = "ATGC"
  let features = @src.dna_feature_vector(seq)

  assert_true(features.length() > 0)
}

///|
test "protein_feature_names" {
  let names = @src.protein_feature_names()

  assert_true(names.length() > 0)
}

///|
test "dna_feature_names" {
  let names = @src.dna_feature_names()

  assert_true(names.length() > 0)
}
