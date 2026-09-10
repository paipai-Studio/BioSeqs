// Whitebox tests for the cif-seqres and cif-atom formats (spec 5.2.3).

///|
let cif_seqres_sample : String = "data_TEST1\n" +
  "#\n" +
  "loop_\n" +
  "_entity_poly_seq.entity_id\n" +
  "_entity_poly_seq.num\n" +
  "_entity_poly_seq.mon_id\n" +
  "1 1 MET\n" +
  "1 2 ALA\n" +
  "1 3 LYS\n" +
  "1 4 GLY\n" +
  "2 1 GLU\n" +
  "2 2 VAL\n" +
  "#\n"

///|
let cif_atom_sample : String = "data_TEST2\n" +
  "#\n" +
  "loop_\n" +
  "_atom_site.group_PDB\n" +
  "_atom_site.id\n" +
  "_atom_site.type_symbol\n" +
  "_atom_site.label_atom_id\n" +
  "_atom_site.label_comp_id\n" +
  "_atom_site.label_asym_id\n" +
  "_atom_site.label_seq_id\n" +
  "_atom_site.Cartn_x\n" +
  "_atom_site.Cartn_y\n" +
  "_atom_site.Cartn_z\n" +
  "_atom_site.occupancy\n" +
  "_atom_site.B_iso_or_equiv\n" +
  "ATOM 1 N N MET A 1 10.0 20.0 30.0 1.0 20.0\n" +
  "ATOM 2 C CA MET A 1 11.0 21.0 31.0 1.0 21.0\n" +
  "ATOM 3 N N ALA A 2 12.0 22.0 32.0 1.0 22.0\n" +
  "ATOM 4 N N GLY B 1 13.0 23.0 33.0 1.0 23.0\n" +
  "ATOM 5 N N PRO B 2 14.0 24.0 34.0 1.0 24.0\n" +
  "#\n"

///|
test "cif_seqres_parses_two_entities" {
  let records = @src.parse_cif_seqres(cif_seqres_sample)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "1")
  assert_eq(records[0].seq.to_string(), "MAKG")
  assert_eq(records[1].id, "2")
  assert_eq(records[1].seq.to_string(), "EV")
}

///|
test "cif_seqres_orders_by_residue_number" {
  // Rows deliberately out of residue-number order
  let content = "data_T\n" +
    "#\n" +
    "loop_\n" +
    "_entity_poly_seq.entity_id\n" +
    "_entity_poly_seq.num\n" +
    "_entity_poly_seq.mon_id\n" +
    "1 3 LYS\n" +
    "1 1 MET\n" +
    "1 2 ALA\n" +
    "#\n"
  let records = @src.parse_cif_seqres(content)
  assert_eq(records.length(), 1)
  assert_eq(records[0].seq.to_string(), "MAK")
}

///|
test "cif_seqres_unknown_residue_maps_to_x" {
  let content = "data_T\n" +
    "#\n" +
    "loop_\n" +
    "_entity_poly_seq.entity_id\n" +
    "_entity_poly_seq.num\n" +
    "_entity_poly_seq.mon_id\n" +
    "1 1 MET\n" +
    "1 2 XYZ\n" +
    "#\n"
  let records = @src.parse_cif_seqres(content)
  assert_eq(records[0].seq.to_string(), "MX")
}

///|
test "cif_seqres_missing_category_returns_empty" {
  let records = @src.parse_cif_seqres("data_EMPTY\n#\n_entry.id EMPTY\n")
  assert_eq(records.length(), 0)
}

///|
test "cif_seqres_malformed_returns_empty" {
  // parse_mmcif is lenient: unrecognised input yields no entity_poly_seq
  // category, so the cif-seqres reader returns an empty collection.
  let records = @src.parse_cif_seqres("this is not mmCIF at all\nbroken {{{")
  assert_eq(records.length(), 0)
}

///|
test "cif_atom_parses_and_dedups_two_chains" {
  let records = @src.parse_cif_atom(cif_atom_sample)
  assert_eq(records.length(), 2)
  // Chain A: MET (deduped by seq_id 1 despite two atom rows) + ALA
  assert_eq(records[0].id, "A")
  assert_eq(records[0].seq.to_string(), "MA")
  assert_eq(records[1].id, "B")
  assert_eq(records[1].seq.to_string(), "GP")
}

///|
test "cif_atom_no_atoms_returns_empty" {
  let records = @src.parse_cif_atom("data_EMPTY\n#\n_entry.id EMPTY\n")
  assert_eq(records.length(), 0)
}

///|
test "cif_atom_skips_water" {
  let content = cif_atom_sample +
    "HETATM 6 O O HOH A . 15.0 25.0 35.0 1.0 30.0\n"
  let records = @src.parse_cif_atom(content)
  assert_eq(records.length(), 2)
  assert_eq(records[0].seq.to_string(), "MA")
}