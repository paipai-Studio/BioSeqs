///|
/// Tests for seq_descriptors.mbt — batch protein / nucleotide descriptors.

///|
test "seq_protein_mw_single_glycine" {
  // residue mass of G (57.0519) + water (18.01528)
  let d = @src.seq_calc_protein_descriptors(@src.Seq::new("G"))
  assert_true(d["mol_weight"] > 75.06 && d["mol_weight"] < 75.08)
  assert_eq(d["length"], 1.0)
}

///|
test "seq_protein_gravy_and_aromaticity" {
  let d = @src.seq_calc_protein_descriptors(@src.Seq::new("AAAA"))
  assert_eq(d["gravy"], 1.8)
  let d2 = @src.seq_calc_protein_descriptors(@src.Seq::new("IIII"))
  assert_eq(d2["gravy"], 4.5)
  let d3 = @src.seq_calc_protein_descriptors(@src.Seq::new("FWYA"))
  assert_eq(d3["aromaticity"], 0.75)
}

///|
test "seq_protein_composition_fractions" {
  let d = @src.seq_calc_protein_descriptors(@src.Seq::new("AAAG"))
  assert_eq(d["A_A"], 0.75)
  assert_eq(d["A_G"], 0.25)
  assert_eq(d["A_W"], 0.0)
}

///|
test "seq_protein_pi_neutral_polypeptide" {
  // AAAA has no ionizable side chains: pI between Cterm and Nterm pKa
  let d = @src.seq_calc_protein_descriptors(@src.Seq::new("AAAA"))
  assert_true(d["isoelectric_point"] > 3.55)
  assert_true(d["isoelectric_point"] < 7.5)
}

///|
test "seq_protein_charge_at_ph7" {
  let basic = @src.seq_protein_charge_at_ph(@src.Seq::new("KKKK"), 7.0)
  let acidic = @src.seq_protein_charge_at_ph(@src.Seq::new("DDDD"), 7.0)
  assert_true(basic > 3.0)
  assert_true(acidic < -3.0)
  // acidic protein has a low pI, basic protein a high one
  let pi_acid = @src.seq_protein_isoelectric_point(@src.Seq::new("DDDD"))
  let pi_basic = @src.seq_protein_isoelectric_point(@src.Seq::new("KKKK"))
  assert_true(pi_acid < 5.0)
  assert_true(pi_basic > 9.0)
}

///|
test "seq_nucleotide_gc_and_cpg" {
  let d = @src.seq_calc_nucleotide_descriptors(@src.Seq::new("GCGCGC"))
  assert_eq(d["GC"], 1.0)
  assert_eq(d["AT"], 0.0)
  assert_eq(d["cpg_count"], 2.0)
  assert_eq(d["purine"], 0.5)
}

///|
test "seq_nucleotide_rna_counts" {
  let d = @src.seq_calc_nucleotide_descriptors(@src.Seq::new("AUUU"))
  assert_eq(d["U"], 0.75)
  assert_eq(d["T"], 0.0)
  assert_eq(d["A"], 0.25)
}

///|
test "seq_nucleotide_shannon_homopolymer" {
  let d = @src.seq_calc_nucleotide_descriptors(@src.Seq::new("AAAAAAAA"))
  assert_eq(d["shannon_bits"], 0.0)
  // equal mix of 4 bases -> 2 bits
  let d2 = @src.seq_calc_nucleotide_descriptors(@src.Seq::new("ACGT"))
  assert_eq(d2["shannon_bits"], 2.0)
}
