///|
/// Tests for ProtParam protein analysis module.

///|
test "protparam_molecular_weight" {
  let protein = @src.ProteinAnalysis::new("ACDEFGHIKLMNPQRSTVWY")
  let mw = protein.molecular_weight()
  assert_true(mw > 2000.0 && mw < 2500.0)
}

///|
test "protparam_count_amino_acids" {
  let protein = @src.ProteinAnalysis::new("AABBCC")
  let counts = protein.count_amino_acids()
  assert_eq(counts.get('A'), Some(2))
  assert_eq(counts.get('B'), Some(2))
  assert_eq(counts.get('C'), Some(2))
}

///|
test "protparam_instability_index" {
  let stable = @src.ProteinAnalysis::new("AVGLEALVLGLLAVGLLPGL")
  let index = stable.instability_index()
  assert_true(index < 40.0)
}

///|
test "protparam_gravy" {
  let hydrophobic = @src.ProteinAnalysis::new("LLLLLL")
  let gravy = hydrophobic.gravy()
  assert_true(gravy > 3.0)
  
  let hydrophilic = @src.ProteinAnalysis::new("RRRRRR")
  let gravy2 = hydrophilic.gravy()
  assert_true(gravy2 < -3.0)
}

///|
test "protparam_aromaticity" {
  let protein = @src.ProteinAnalysis::new("FFWWYY")
  let aroma = protein.aromaticity()
  assert_eq(aroma, 1.0)
}

///|
test "protparam_isoelectric_point" {
  let acidic = @src.ProteinAnalysis::new("EEEEEE")
  let pi = acidic.isoelectric_point()
  assert_true(pi < 5.0)
  
  let basic = @src.ProteinAnalysis::new("KKKKKK")
  let pi2 = basic.isoelectric_point()
  assert_true(pi2 > 9.0)
}

///|
test "protparam_charge_at_ph" {
  let neutral = @src.ProteinAnalysis::new("G")
  let charge = neutral.charge_at_ph(7.0)
  assert_true(charge > -0.5 && charge < 0.5)
}

///|
test "protparam_secondary_structure" {
  let protein = @src.ProteinAnalysis::new("AEAEAE")
  let (h, s, c) = protein.secondary_structure()
  assert_true(h > 1.0)
}

///|
test "protparam_amino_acid_composition" {
  let protein = @src.ProteinAnalysis::new("AABB")
  let comp = protein.amino_acid_composition()
  assert_eq(comp.get('A'), Some(0.5))
  assert_eq(comp.get('B'), Some(0.5))
}

///|
test "protparam_signal_peptide" {
  let sp_seq = "MRVLLFLLLIAAPVLAAPVTLSRGA"
  let protein = @src.ProteinAnalysis::new(sp_seq)
  let result = protein.predict_signal_peptide()
  assert_true(result.has_signal_peptide)
}

///|
test "protparam_length" {
  let protein = @src.ProteinAnalysis::new("ABCDE")
  assert_eq(protein.length(), 5)
}

///|
test "protparam_get_amino_acids_percent" {
  let protein = @src.ProteinAnalysis::new("AABBCC")
  let percent = protein.get_amino_acids_percent()
  assert_eq(percent.get('A'), Some(33.33333333333333))
}

///|
test "protparam_empty_sequence" {
  let protein = @src.ProteinAnalysis::new("")
  assert_eq(protein.length(), 0)
  assert_eq(protein.molecular_weight(), 0.0)
  assert_eq(protein.instability_index(), 0.0)
}

///|
test "protparam_small_sequence" {
  let protein = @src.ProteinAnalysis::new("A")
  assert_eq(protein.length(), 1)
  assert_true(protein.molecular_weight() > 70.0)
}

///|
test "protparam_instability_index_small" {
  let protein = @src.ProteinAnalysis::new("AB")
  let index = protein.instability_index()
  assert_true(index >= 0.0)
}