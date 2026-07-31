///|
/// Tests for mol_wt module (molecular weight calculation).

///|
test "nucleotide_weight_a" {
  let w = @src.nucleotide_weight("A")
  assert_true(w > 300.0)
  assert_true(w < 400.0)
}

///|
test "nucleotide_weight_c" {
  let w = @src.nucleotide_weight("C")
  assert_true(w > 280.0)
  assert_true(w < 380.0)
}

///|
test "nucleotide_weight_g" {
  let w = @src.nucleotide_weight("G")
  assert_true(w > 320.0)
  assert_true(w < 400.0)
}

///|
test "nucleotide_weight_t" {
  let w = @src.nucleotide_weight("T")
  assert_true(w > 300.0)
  assert_true(w < 380.0)
}

///|
test "nucleotide_weight_u" {
  let w = @src.nucleotide_weight("U")
  assert_true(w > 300.0)
  assert_true(w < 380.0)
}

///|
test "nucleotide_weight_unknown" {
  let w = @src.nucleotide_weight("Z")
  assert_true(w > 0.0) // Should fall back to default
}

///|
test "amino_acid_weight_a" {
  let w = @src.amino_acid_weight("A")
  assert_true(w > 70.0)
  assert_true(w < 100.0)
}

///|
test "amino_acid_weight_g" {
  let w = @src.amino_acid_weight("G")
  assert_true(w > 60.0)
  assert_true(w < 90.0)
}

///|
test "amino_acid_weight_w" {
  let w = @src.amino_acid_weight("W")
  assert_true(w > 180.0)
  assert_true(w < 230.0)
}

///|
test "amino_acid_weight_unknown" {
  let w = @src.amino_acid_weight("Z")
  assert_true(w > 0.0) // Should fall back to average
}

///|
test "mol_weight_dna_empty" {
  let mw = @src.mol_weight_dna("")
  assert_eq(mw, 0.0)
}

///|
test "mol_weight_dna_single" {
  let mw = @src.mol_weight_dna("A")
  assert_true(mw > 300.0) // Should be approximately 335 - 0*18 = 335
  assert_true(mw < 400.0)
}

///|
test "mol_weight_dna_known" {
  // ACGT: 335 + 311 + 351 + 326 - 3*18.015
  let mw = @src.mol_weight_dna("ACGT")
  assert_true(mw > 1200.0)
  assert_true(mw < 1400.0)
}

///|
test "mol_weight_rna" {
  let mw = @src.mol_weight_rna("ACGT")
  // RNA has U instead of T, but weights are similar
  assert_true(mw > 1200.0)
}

///|
test "mol_weight_protein_empty" {
  let mw = @src.mol_weight_protein("")
  assert_eq(mw, 0.0)
}

///|
test "mol_weight_protein_single" {
  let mw = @src.mol_weight_protein("G")
  assert_true(mw > 50.0) // Glycine - 0*water = 75
  assert_true(mw < 100.0)
}

///|
test "mol_weight_protein_known" {
  // "MV" = Met(149) + Val(117) - 1*18 = 248
  let mw = @src.mol_weight_protein("MV")
  assert_true(mw > 200.0)
  assert_true(mw < 280.0)
}

///|
test "extinction_coefficient_280_no_trp" {
  let ext = @src.extinction_coefficient_280("A")
  assert_eq(ext, 0.0)
}

///|
test "extinction_coefficient_280_trp" {
  let ext = @src.extinction_coefficient_280("W")
  assert_eq(ext, 5500.0)
}

///|
test "extinction_coefficient_280_tyr" {
  let ext = @src.extinction_coefficient_280("Y")
  assert_eq(ext, 1490.0)
}

///|
test "extinction_coefficient_280_cys" {
  let ext = @src.extinction_coefficient_280("C")
  assert_eq(ext, 0.0) // Single Cys = 0 cystines
}

///|
test "extinction_coefficient_280_two_cys" {
  let ext = @src.extinction_coefficient_280("CC")
  assert_eq(ext, 125.0) // 1 cystine
}

///|
test "extinction_coefficient_280_combined" {
  let ext = @src.extinction_coefficient_280("WYCC")
  assert_eq(ext, 5500.0 + 1490.0 + 125.0)
}

///|
test "absorbance_280_zero_mw" {
  let abs = @src.absorbance_280("")
  assert_eq(abs, 0.0)
}

///|
test "isoelectric_point_empty" {
  let pI = @src.mol_wt_isoelectric_point("")
  assert_eq(pI, 7.0)
}

///|
test "isoelectric_point_basic" {
  let pI = @src.mol_wt_isoelectric_point("KKK")
  assert_true(pI > 7.0) // Basic protein
}

///|
test "isoelectric_point_acidic" {
  let pI = @src.mol_wt_isoelectric_point("DEE")
  assert_true(pI < 7.0) // Acidic protein
}

///|
test "mol_weight_summary_dna" {
  let s = @src.mol_weight_summary("ACGT")
  assert_true(s.length() > 0)
  assert_true(s.contains("DNA"))
}

///|
test "mol_weight_summary_protein" {
  let s = @src.mol_weight_summary("MVL", is_protein=true)
  assert_true(s.length() > 0)
  assert_true(s.contains("Protein"))
}

///|
test "mol_wt_sample_dna" {
  let dna = @src.mol_wt_sample_dna()
  assert_true(dna.length() > 0)
}

///|
test "mol_wt_sample_protein" {
  let prot = @src.mol_wt_sample_protein()
  assert_true(prot.length() > 0)
}

///|
test "atomic_weights_default" {
  let aw = @src.default_atomic_weights()
  assert_true(aw.h > 1.0)
  assert_true(aw.c > 10.0)
  assert_true(aw.n > 10.0)
  assert_true(aw.o > 10.0)
  assert_true(aw.p > 20.0)
  assert_true(aw.s > 20.0)
}
