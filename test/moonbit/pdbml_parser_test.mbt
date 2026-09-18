///|
/// Tests for the Bio.PDB.PDBMLParser port. Expected values come from running
/// Biopython 1.88's PDBMLParser on the same fixture. Note: resolution is
/// expected to be Some(1.5) here because this port fixes Biopython's
/// `_parse_resolution_from` presence-check bug (Biopython always returns None).


///|
fn pm_xml() -> String {
  "<PDBx:datablock xmlns:PDBx='http://pdbml.pdb.org/schema/pdbx-v50.xsd'>" +
  "<PDBx:entryCategory><PDBx:entry id='1ABC'/></PDBx:entryCategory>" +
  "<PDBx:structCategory><PDBx:struct><PDBx:title>TEST STRUCTURE</PDBx:title></PDBx:struct></PDBx:structCategory>" +
  "<PDBx:struct_keywordsCategory><PDBx:struct_keywords><PDBx:text>test keyword</PDBx:text></PDBx:struct_keywords></PDBx:struct_keywordsCategory>" +
  "<PDBx:pdbx_database_statusCategory><PDBx:pdbx_database_status><PDBx:recvd_initial_deposition_date>2020-01-02</PDBx:recvd_initial_deposition_date></PDBx:pdbx_database_status></PDBx:pdbx_database_statusCategory>" +
  "<PDBx:exptlCategory><PDBx:exptl method='X-RAY DIFFRACTION'/></PDBx:exptlCategory>" +
  "<PDBx:refineCategory><PDBx:refine><PDBx:ls_d_res_high>1.50</PDBx:ls_d_res_high></PDBx:refine></PDBx:refineCategory>" +
  "<PDBx:atom_siteCategory>" +
  "<PDBx:atom_site id='1'><PDBx:group_PDB>ATOM</PDBx:group_PDB><PDBx:label_comp_id>ALA</PDBx:label_comp_id><PDBx:auth_seq_id>1</PDBx:auth_seq_id><PDBx:pdbx_PDB_ins_code>?</PDBx:pdbx_PDB_ins_code><PDBx:auth_asym_id>A</PDBx:auth_asym_id><PDBx:pdbx_PDB_model_num>1</PDBx:pdbx_PDB_model_num><PDBx:label_atom_id>CA</PDBx:label_atom_id><PDBx:type_symbol>C</PDBx:type_symbol><PDBx:label_alt_id>.</PDBx:label_alt_id><PDBx:Cartn_x>1.100</PDBx:Cartn_x><PDBx:Cartn_y>2.200</PDBx:Cartn_y><PDBx:Cartn_z>3.300</PDBx:Cartn_z><PDBx:occupancy>1.00</PDBx:occupancy><PDBx:B_iso_or_equiv>20.50</PDBx:B_iso_or_equiv></PDBx:atom_site>" +
  "<PDBx:atom_site id='2'><PDBx:group_PDB>ATOM</PDBx:group_PDB><PDBx:label_comp_id>ALA</PDBx:label_comp_id><PDBx:auth_seq_id>1</PDBx:auth_seq_id><PDBx:pdbx_PDB_ins_code>?</PDBx:pdbx_PDB_ins_code><PDBx:auth_asym_id>A</PDBx:auth_asym_id><PDBx:pdbx_PDB_model_num>1</PDBx:pdbx_PDB_model_num><PDBx:label_atom_id>CB</PDBx:label_atom_id><PDBx:type_symbol>C</PDBx:type_symbol><PDBx:label_alt_id>.</PDBx:label_alt_id><PDBx:Cartn_x>1.500</PDBx:Cartn_x><PDBx:Cartn_y>2.600</PDBx:Cartn_y><PDBx:Cartn_z>3.700</PDBx:Cartn_z><PDBx:occupancy>1.00</PDBx:occupancy><PDBx:B_iso_or_equiv>22.10</PDBx:B_iso_or_equiv></PDBx:atom_site>" +
  "<PDBx:atom_site id='3'><PDBx:group_PDB>HETATM</PDBx:group_PDB><PDBx:label_comp_id>HOH</PDBx:label_comp_id><PDBx:auth_seq_id>2</PDBx:auth_seq_id><PDBx:pdbx_PDB_ins_code>?</PDBx:pdbx_PDB_ins_code><PDBx:auth_asym_id>B</PDBx:auth_asym_id><PDBx:pdbx_PDB_model_num>1</PDBx:pdbx_PDB_model_num><PDBx:label_atom_id>O</PDBx:label_atom_id><PDBx:type_symbol>O</PDBx:type_symbol><PDBx:label_alt_id>.</PDBx:label_alt_id><PDBx:Cartn_x>4.000</PDBx:Cartn_x><PDBx:Cartn_y>5.000</PDBx:Cartn_y><PDBx:Cartn_z>6.000</PDBx:Cartn_z><PDBx:occupancy>0.50</PDBx:occupancy><PDBx:B_iso_or_equiv>30.00</PDBx:B_iso_or_equiv></PDBx:atom_site>" +
  "</PDBx:atom_siteCategory></PDBx:datablock>"
}

///|
test "pdbml_header" {
  let s = @src.pdbml_get_structure(pm_xml())
  assert_eq(s.header.name, "TEST STRUCTURE")
  assert_eq(s.header.head, "test keyword")
  assert_eq(s.header.idcode, "1ABC")
  assert_eq(s.header.deposition_date, "2020-01-02")
  assert_eq(s.header.structure_method, "X-RAY DIFFRACTION")
  match s.header.resolution {
    Some(v) => assert_true(assert_close(v, 1.5, 1.0e-9))
    None => abort("resolution should be parsed")
  }
  assert_eq(s.idcode, "1ABC")
}

///|
test "pdbml_structure_hierarchy" {
  let s = @src.pdbml_get_structure(pm_xml())
  assert_eq(s.models.length(), 1)
  assert_eq(s.models[0].serial_number, 1)
  assert_eq(s.models[0].chains.length(), 2)

  let c0 = s.models[0].chains[0]
  assert_eq(c0.chain_id, "A")
  assert_eq(c0.residues.length(), 1)
  let r0 = c0.residues[0]
  assert_eq(r0.component_id, "ALA")
  assert_eq(r0.hetero_field, " ")
  assert_eq(r0.sequence_id, 1)
  assert_eq(r0.insertion_code, "?")
  assert_eq(r0.atoms.length(), 2)

  let a0 = r0.atoms[0]
  assert_eq(a0.name, "CA")
  assert_eq(a0.fullname, "CA")
  assert_eq(a0.serial_number, 1)
  assert_eq(a0.element, "C")
  assert_eq(a0.altloc, ".")
  assert_true(assert_close(a0.coord[0], 1.1, 1.0e-9))
  assert_true(assert_close(a0.coord[1], 2.2, 1.0e-9))
  assert_true(assert_close(a0.coord[2], 3.3, 1.0e-9))
  assert_true(assert_close(a0.b_factor, 20.5, 1.0e-9))
  assert_true(assert_close(a0.occupancy, 1.0, 1.0e-9))

  let a1 = r0.atoms[1]
  assert_eq(a1.name, "CB")
  assert_eq(a1.serial_number, 2)
  assert_true(assert_close(a1.coord[0], 1.5, 1.0e-9))
  assert_true(assert_close(a1.b_factor, 22.1, 1.0e-9))

  let c1 = s.models[0].chains[1]
  assert_eq(c1.chain_id, "B")
  let r1 = c1.residues[0]
  assert_eq(r1.component_id, "HOH")
  assert_eq(r1.hetero_field, "W")
  assert_eq(r1.sequence_id, 2)
  assert_eq(r1.insertion_code, "?")
  assert_eq(r1.atoms.length(), 1)
  let a2 = r1.atoms[0]
  assert_eq(a2.name, "O")
  assert_eq(a2.serial_number, 3)
  assert_eq(a2.element, "O")
  assert_true(assert_close(a2.coord[2], 6.0, 1.0e-9))
  assert_true(assert_close(a2.b_factor, 30.0, 1.0e-9))
  assert_true(assert_close(a2.occupancy, 0.5, 1.0e-9))
}

///|
test "pdbml_missing_resolution" {
  let xml = "<PDBx:datablock xmlns:PDBx='x'><PDBx:entryCategory><PDBx:entry id='X'/></PDBx:entryCategory></PDBx:datablock>"
  let s = @src.pdbml_get_structure(xml)
  match s.header.resolution {
    None => ()
    Some(_) => abort("resolution should be none when absent")
  }
}