///|
test "dssp_predict_secondary_structure" {
  let sequence = "MKLILVLLVSLSL"
  let ss = @src.predict_secondary_structure(sequence)
  assert_eq(ss.length(), sequence.length())
}

///|
test "dssp_calculate_solvent_accessibility" {
  let sequence = "MKLILVLLVSLSL"
  let structure = "HHHHHHHHHHHH"
  let acc = @src.calculate_solvent_accessibility(sequence, structure)
  assert_eq(acc.length(), sequence.length())
}

///|
test "dssp_calculate_ramachandran_plot" {
  let data = @src.create_example_dssp_data()
  let plot = @src.calculate_ramachandran_plot(data.records)
  assert_eq(plot.length(), data.records.length())
}

///|
test "dssp_analyze_structure_composition" {
  let data = @src.create_example_dssp_data()
  let comp = @src.analyze_structure_composition(data.records)
  assert_eq(comp.get("H").unwrap_or(0), 3)
}

///|
test "dssp_create_example_data" {
  let data = @src.create_example_dssp_data()
  assert_eq(data.records.length(), 8)
}