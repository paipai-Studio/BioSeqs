///|
test "polypeptide_calculate_composition" {
  let sequence = "MKLILVLLVSLSL"
  let comp = @src.calculate_composition(sequence)
  assert_eq(comp.total_residues, sequence.length())
}

///|
test "polypeptide_calculate_percent_composition" {
  let sequence = "AAAA"
  let comp = @src.calculate_composition(sequence)
  let percent = @src.calculate_percent_composition(comp)
  assert_eq(percent.get("A").unwrap_or(0.0), 100.0)
}

///|
test "polypeptide_calculate_hydrophobicity" {
  let sequence = "MKLILVLLVSLSL"
  let profile = @src.calculate_hydrophobicity(sequence, "kyte-doolittle")
  assert_eq(profile.values.length(), sequence.length())
}

///|
test "polypeptide_calculate_hydrophobicity_window" {
  let sequence = "MKLILVLLVSLSL"
  let profile = @src.calculate_hydrophobicity_window(sequence, 5, "kyte-doolittle")
  assert_eq(profile.values.length(), sequence.length())
}

///|
test "polypeptide_predict_transmembrane_regions" {
  let sequence = "MVLILVLLVSLSLILVLLV"
  let regions = @src.predict_transmembrane_regions(sequence, 7, 2.0)
  assert_eq(regions.length() >= 1, true)
}

///|
test "polypeptide_calculate_isoelectric_point" {
  let sequence = "MKLILVLLVSLSL"
  let pi = @src.calculate_isoelectric_point(sequence)
  assert_eq(pi >= 0.0, true)
  assert_eq(pi <= 14.0, true)
}

///|
test "polypeptide_calculate_aromaticity" {
  let sequence = "FFYYWW"
  let aromaticity = @src.calculate_aromaticity(sequence)
  assert_eq(aromaticity, 1.0)
}

///|
test "polypeptide_create_example_data" {
  let (seq, comp, prof) = @src.create_example_polypeptide_data()
  assert_eq(seq.length() > 0, true)
  assert_eq(comp.total_residues > 0, true)
  assert_eq(prof.values.length() > 0, true)
}