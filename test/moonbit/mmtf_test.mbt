///|
/// Tests for Bio.PDB.mmtf module.

// ============================================================================
// MmtfSecStruct enum
// ============================================================================

///|
test "mmtf_sec_struct_from_int_pii_helix" {
  let s = @src.MmtfSecStruct::from_int(0)
  assert_true(s is @src.MmtfSecStruct::MmtfSecPpiHelix)
  assert_eq(s.to_string(), "pi_helix")
}

///|
test "mmtf_sec_struct_from_int_bend" {
  let s = @src.MmtfSecStruct::from_int(1)
  assert_true(s is @src.MmtfSecStruct::MmtfSecBend)
  assert_eq(s.to_string(), "bend")
}

///|
test "mmtf_sec_struct_from_int_alpha_helix" {
  let s = @src.MmtfSecStruct::from_int(2)
  assert_true(s is @src.MmtfSecStruct::MmtfSecAlphaHelix)
  assert_eq(s.to_string(), "alpha_helix")
}

///|
test "mmtf_sec_struct_from_int_extended" {
  let s = @src.MmtfSecStruct::from_int(3)
  assert_true(s is @src.MmtfSecStruct::MmtfSecExtended)
  assert_eq(s.to_string(), "extended")
}

///|
test "mmtf_sec_struct_from_int_three_ten_helix" {
  let s = @src.MmtfSecStruct::from_int(4)
  assert_true(s is @src.MmtfSecStruct::MmtfSecThreeTenHelix)
  assert_eq(s.to_string(), "310_helix")
}

///|
test "mmtf_sec_struct_from_int_bridge" {
  let s = @src.MmtfSecStruct::from_int(5)
  assert_true(s is @src.MmtfSecStruct::MmtfSecBridge)
  assert_eq(s.to_string(), "bridge")
}

///|
test "mmtf_sec_struct_from_int_turn" {
  let s = @src.MmtfSecStruct::from_int(6)
  assert_true(s is @src.MmtfSecStruct::MmtfSecTurn)
  assert_eq(s.to_string(), "turn")
}

///|
test "mmtf_sec_struct_from_int_coil" {
  let s = @src.MmtfSecStruct::from_int(7)
  assert_true(s is @src.MmtfSecStruct::MmtfSecCoil)
  assert_eq(s.to_string(), "coil")
}

///|
test "mmtf_sec_struct_from_int_undefined" {
  let s = @src.MmtfSecStruct::from_int(99)
  assert_true(s is @src.MmtfSecStruct::MmtfSecUndefined)
  assert_eq(s.to_string(), "undefined")
}

///|
test "mmtf_sec_struct_to_string_round_trip" {
  let codes = [0, 1, 2, 3, 4, 5, 6, 7, 99]
  let names = [
    "pi_helix",
    "bend",
    "alpha_helix",
    "extended",
    "310_helix",
    "bridge",
    "turn",
    "coil",
    "undefined",
  ]
  for i in 0..<codes.length() {
    let s = @src.MmtfSecStruct::from_int(codes[i])
    assert_eq(s.to_string(), names[i])
  }
}

// ============================================================================
// MmtfEntityType
// ============================================================================

///|
test "mmtf_entity_type_to_string" {
  // Construct via from_string since variants cannot be created cross-package
  assert_eq(@src.MmtfEntityType::from_string("polymer").to_string(), "polymer")
  assert_eq(
    @src.MmtfEntityType::from_string("non-polymer").to_string(),
    "non-polymer",
  )
  assert_eq(@src.MmtfEntityType::from_string("water").to_string(), "water")
  assert_eq(
    @src.MmtfEntityType::from_string("branched").to_string(),
    "branched",
  )
  assert_eq(
    @src.MmtfEntityType::from_string("undefined").to_string(),
    "undefined",
  )
}

///|
test "mmtf_entity_type_from_string_round_trip" {
  let names = [
    "polymer",
    "non-polymer",
    "water",
    "branched",
    "undefined",
  ]
  for n in names {
    let t = @src.MmtfEntityType::from_string(n)
    assert_eq(t.to_string(), n)
  }
}

///|
test "mmtf_entity_type_from_string_unknown" {
  let t = @src.MmtfEntityType::from_string("nonsense")
  assert_true(t is @src.MmtfEntityType::MmtfEntityUndefined)
  assert_eq(t.to_string(), "undefined")
}

// ============================================================================
// MmtfAtom
// ============================================================================

///|
test "mmtf_atom_new_defaults" {
  let a = @src.MmtfAtom::new(name="CA", serial=1, x=1.0, y=2.0, z=3.0)
  assert_eq(a.name, "CA")
  assert_eq(a.serial, 1)
  assert_eq(a.alt_loc_id, "")
  assert_eq(a.x, 1.0)
  assert_eq(a.y, 2.0)
  assert_eq(a.z, 3.0)
  assert_true((a.occupancy - 1.0).abs() < 0.001)
  assert_true((a.b_factor - 0.0).abs() < 0.001)
  assert_eq(a.element, "")
  assert_eq(a.charge, 0)
}

///|
test "mmtf_atom_new_full" {
  let a = @src.MmtfAtom::new(
    name="FE",
    serial=42,
    alt_loc_id="A",
    x=10.5,
    y=20.5,
    z=30.5,
    occupancy=0.5,
    b_factor=40.0,
    element="FE",
    charge=2,
  )
  assert_eq(a.name, "FE")
  assert_eq(a.serial, 42)
  assert_eq(a.alt_loc_id, "A")
  assert_eq(a.x, 10.5)
  assert_eq(a.y, 20.5)
  assert_eq(a.z, 30.5)
  assert_eq(a.occupancy, 0.5)
  assert_eq(a.b_factor, 40.0)
  assert_eq(a.element, "FE")
  assert_eq(a.charge, 2)
}

// ============================================================================
// MmtfGroup
// ============================================================================

///|
test "mmtf_group_new_defaults" {
  let g = @src.MmtfGroup::new(group_name="ALA", group_number=1)
  assert_eq(g.group_name, "ALA")
  assert_eq(g.group_number, 1)
  assert_eq(g.insertion_code, "")
  assert_eq(g.group_type, "")
  assert_eq(g.single_letter_code, "")
  assert_true(g.sec_struct is @src.MmtfSecStruct::MmtfSecUndefined)
  assert_eq(g.sequence_index, -1)
  assert_eq(g.get_num_atoms(), 0)
}

///|
test "mmtf_group_add_atom_and_get_atoms" {
  let g = @src.MmtfGroup::new(
    group_name="ALA",
    group_number=1,
    single_letter_code="A",
    sec_struct=@src.MmtfSecStruct::from_int(2),
    sequence_index=0,
  )
  g.add_atom(@src.MmtfAtom::new(name="N", serial=1, x=0.0, y=0.0, z=0.0))
  g.add_atom(@src.MmtfAtom::new(name="CA", serial=2, x=1.0, y=0.0, z=0.0))
  assert_eq(g.get_num_atoms(), 2)
  let atoms = g.get_atoms()
  assert_eq(atoms.length(), 2)
  assert_eq(atoms[0].name, "N")
  assert_eq(atoms[1].name, "CA")
}

// ============================================================================
// MmtfChain
// ============================================================================

///|
test "mmtf_chain_new_defaults" {
  let c = @src.MmtfChain::new()
  assert_eq(c.chain_id, "")
  assert_eq(c.chain_name, "")
  assert_eq(c.get_num_groups(), 0)
}

///|
test "mmtf_chain_add_group_and_get_groups" {
  let c = @src.MmtfChain::new(chain_id="A", chain_name="A")
  let g1 = @src.MmtfGroup::new(group_name="ALA", group_number=1)
  let g2 = @src.MmtfGroup::new(group_name="GLY", group_number=2)
  c.add_group(g1)
  c.add_group(g2)
  assert_eq(c.get_num_groups(), 2)
  let groups = c.get_groups()
  assert_eq(groups.length(), 2)
  assert_eq(groups[0].group_name, "ALA")
  assert_eq(groups[1].group_name, "GLY")
}

// ============================================================================
// MmtfModel
// ============================================================================

///|
test "mmtf_model_new_default_id" {
  let m = @src.MmtfModel::new()
  assert_eq(m.model_id, 0)
  assert_eq(m.get_num_chains(), 0)
}

///|
test "mmtf_model_add_chain_and_get_chains" {
  let m = @src.MmtfModel::new(model_id=3)
  let c1 = @src.MmtfChain::new(chain_id="A", chain_name="A")
  let c2 = @src.MmtfChain::new(chain_id="B", chain_name="B")
  m.add_chain(c1)
  m.add_chain(c2)
  assert_eq(m.model_id, 3)
  assert_eq(m.get_num_chains(), 2)
  let chains = m.get_chains()
  assert_eq(chains[0].chain_id, "A")
  assert_eq(chains[1].chain_id, "B")
}

// ============================================================================
// MmtfStructure - construction and accessors
// ============================================================================

///|
test "mmtf_structure_new_empty" {
  let s = @src.MmtfStructure::new()
  assert_eq(s.structure_id, "")
  assert_eq(s.get_num_models(), 0)
  assert_eq(s.count_chains(), 0)
  assert_eq(s.count_groups(), 0)
  assert_eq(s.count_atoms(), 0)
  assert_eq(s.get_chains().length(), 0)
  assert_eq(s.get_groups().length(), 0)
  assert_eq(s.get_atoms().length(), 0)
}

///|
test "mmtf_structure_add_model_and_iterate" {
  let s = @src.MmtfStructure::new()
  let model = @src.MmtfModel::new(model_id=0)
  let chain = @src.MmtfChain::new(chain_id="A", chain_name="A")
  let group = @src.MmtfGroup::new(group_name="ALA", group_number=1)
  group.add_atom(@src.MmtfAtom::new(name="CA", serial=1, x=1.0, y=2.0, z=3.0))
  chain.add_group(group)
  model.add_chain(chain)
  s.add_model(model)
  assert_eq(s.get_num_models(), 1)
  assert_eq(s.count_chains(), 1)
  assert_eq(s.count_groups(), 1)
  assert_eq(s.count_atoms(), 1)
  assert_eq(s.get_models().length(), 1)
  assert_eq(s.get_chains().length(), 1)
  assert_eq(s.get_groups().length(), 1)
  let atoms = s.get_atoms()
  assert_eq(atoms.length(), 1)
  assert_eq(atoms[0].name, "CA")
}

///|
test "mmtf_structure_header_get_set" {
  // pub mut fields cannot be assigned cross-package, so we transfer a header
  // from the sample structure to a fresh structure via get/set.
  let sample = @src.mmtf_sample()
  let sample_header = sample.get_header()
  let s = @src.MmtfStructure::new()
  s.set_header(sample_header)
  let got = s.get_header()
  assert_eq(got.title, "Sample protein structure")
  assert_eq(got.deposition_date, "2024-01-01")
  assert_eq(got.release_date, "2024-06-01")
  assert_true((got.resolution - 1.5).abs() < 0.001)
  assert_true((got.r_work - 0.18).abs() < 0.001)
  assert_true((got.r_free - 0.22).abs() < 0.001)
  assert_eq(got.experimental_methods.length(), 1)
  assert_eq(got.experimental_methods[0], "X-RAY DIFFRACTION")
}

///|
test "mmtf_structure_xtal_info_get_set" {
  // pub mut fields cannot be assigned cross-package, so we transfer xtal_info
  // from the sample structure to a fresh structure via get/set.
  let sample = @src.mmtf_sample()
  let sample_xtal = sample.get_xtal_info()
  let s = @src.MmtfStructure::new()
  s.set_xtal_info(sample_xtal)
  let got = s.get_xtal_info()
  assert_eq(got.space_group, "P 21 21 21")
  assert_eq(got.unit_cell.length(), 6)
  assert_eq(got.unit_cell[0], 30.0)
  assert_eq(got.unit_cell[1], 40.0)
  assert_eq(got.unit_cell[2], 50.0)
  assert_eq(got.unit_cell[3], 90.0)
  assert_eq(got.unit_cell[4], 90.0)
  assert_eq(got.unit_cell[5], 90.0)
}

///|
test "mmtf_header_new_defaults" {
  let h = @src.MmtfHeader::new()
  assert_eq(h.title, "")
  assert_eq(h.deposition_date, "")
  assert_eq(h.release_date, "")
  assert_true((h.resolution - 0.0).abs() < 0.001)
  assert_true((h.r_work - 0.0).abs() < 0.001)
  assert_true((h.r_free - 0.0).abs() < 0.001)
  assert_eq(h.experimental_methods.length(), 0)
}

///|
test "mmtf_xtal_info_new_defaults" {
  let x = @src.MmtfXtalInfo::new()
  assert_eq(x.space_group, "")
  assert_eq(x.unit_cell.length(), 6)
  for i in 0..<6 {
    assert_true((x.unit_cell[i] - 0.0).abs() < 0.001)
  }
}

///|
test "mmtf_bond_new" {
  let b = @src.MmtfBond::new(atom_index1=3, atom_index2=7, bond_order=2)
  assert_eq(b.atom_index1, 3)
  assert_eq(b.atom_index2, 7)
  assert_eq(b.bond_order, 2)
}

///|
test "mmtf_bond_new_default_order" {
  let b = @src.MmtfBond::new(atom_index1=0, atom_index2=1)
  assert_eq(b.bond_order, 1)
}

///|
test "mmtf_entity_new" {
  let e = @src.MmtfEntity::new(
    chain_indices=[0, 1],
    sequence="MKL",
    description="test entity",
    entity_type=@src.MmtfEntityType::from_string("polymer"),
  )
  assert_eq(e.chain_indices.length(), 2)
  assert_eq(e.chain_indices[0], 0)
  assert_eq(e.chain_indices[1], 1)
  assert_eq(e.sequence, "MKL")
  assert_eq(e.description, "test entity")
  assert_eq(e.entity_type.to_string(), "polymer")
}

///|
test "mmtf_group_type_new" {
  let g = @src.MmtfGroupType::new(
    group_name="ALA",
    chem_comp_type="L-PEPTIDE LINKING",
    single_letter_code="A",
  )
  assert_eq(g.group_name, "ALA")
  assert_eq(g.chem_comp_type, "L-PEPTIDE LINKING")
  assert_eq(g.single_letter_code, "A")
  assert_eq(g.atom_name_list.length(), 0)
  assert_eq(g.element_list.length(), 0)
  assert_eq(g.bond_atom_list.length(), 0)
  assert_eq(g.bond_order_list.length(), 0)
  assert_eq(g.formal_charge_list.length(), 0)
}

// ============================================================================
// mmtf_sample
// ============================================================================

///|
test "mmtf_sample_structure_id" {
  let s = @src.mmtf_sample()
  assert_eq(s.structure_id, "1XXX")
}

///|
test "mmtf_sample_header_fields" {
  let s = @src.mmtf_sample()
  let h = s.get_header()
  assert_eq(h.title, "Sample protein structure")
  assert_eq(h.deposition_date, "2024-01-01")
  assert_eq(h.release_date, "2024-06-01")
  assert_true((h.resolution - 1.5).abs() < 0.001)
  assert_true((h.r_work - 0.18).abs() < 0.001)
  assert_true((h.r_free - 0.22).abs() < 0.001)
  assert_eq(h.experimental_methods.length(), 1)
  assert_eq(h.experimental_methods[0], "X-RAY DIFFRACTION")
}

///|
test "mmtf_sample_xtal_info_fields" {
  let s = @src.mmtf_sample()
  let x = s.get_xtal_info()
  assert_eq(x.space_group, "P 21 21 21")
  assert_eq(x.unit_cell.length(), 6)
  assert_true((x.unit_cell[0] - 30.0).abs() < 0.001)
  assert_true((x.unit_cell[1] - 40.0).abs() < 0.001)
  assert_true((x.unit_cell[2] - 50.0).abs() < 0.001)
  assert_true((x.unit_cell[3] - 90.0).abs() < 0.001)
  assert_true((x.unit_cell[4] - 90.0).abs() < 0.001)
  assert_true((x.unit_cell[5] - 90.0).abs() < 0.001)
}

///|
test "mmtf_sample_counts" {
  let s = @src.mmtf_sample()
  assert_eq(s.get_num_models(), 1)
  assert_eq(s.count_chains(), 2)
  assert_eq(s.count_groups(), 3)
  assert_eq(s.count_atoms(), 6)
}

///|
test "mmtf_sample_total_num_metadata" {
  let s = @src.mmtf_sample()
  assert_eq(s.total_num_models, 1)
  assert_eq(s.total_num_chains, 2)
  assert_eq(s.total_num_groups, 4)
  assert_eq(s.total_num_atoms, 10)
  assert_eq(s.total_num_bonds, 2)
}

///|
test "mmtf_sample_entities" {
  let s = @src.mmtf_sample()
  assert_eq(s.entities.length(), 1)
  let e = s.entities[0]
  assert_eq(e.sequence, "ACDE")
  assert_eq(e.description, "sample protein")
  assert_eq(e.entity_type.to_string(), "polymer")
  assert_eq(e.chain_indices.length(), 1)
  assert_eq(e.chain_indices[0], 0)
}

///|
test "mmtf_sample_bonds" {
  let s = @src.mmtf_sample()
  assert_eq(s.bonds.length(), 2)
  assert_eq(s.bonds[0].atom_index1, 0)
  assert_eq(s.bonds[0].atom_index2, 1)
  assert_eq(s.bonds[0].bond_order, 1)
  assert_eq(s.bonds[1].atom_index1, 1)
  assert_eq(s.bonds[1].atom_index2, 2)
}

///|
test "mmtf_sample_chain_contents" {
  let s = @src.mmtf_sample()
  let models = s.get_models()
  let model = models[0]
  let chains = model.get_chains()
  let chain_a = chains[0]
  assert_eq(chain_a.chain_id, "A")
  assert_eq(chain_a.chain_name, "A")
  assert_eq(chain_a.get_num_groups(), 2)
  let groups = chain_a.get_groups()
  assert_eq(groups[0].group_name, "ALA")
  assert_eq(groups[0].group_number, 1)
  assert_eq(groups[0].single_letter_code, "A")
  assert_eq(groups[0].sec_struct.to_string(), "alpha_helix")
  assert_eq(groups[0].get_num_atoms(), 3)
  assert_eq(groups[1].group_name, "CYS")
  assert_eq(groups[1].group_number, 2)
  assert_eq(groups[1].single_letter_code, "C")
  assert_eq(groups[1].sec_struct.to_string(), "extended")
  assert_eq(groups[1].get_num_atoms(), 2)
  let chain_b = chains[1]
  assert_eq(chain_b.chain_id, "B")
  assert_eq(chain_b.get_num_groups(), 1)
  let b_groups = chain_b.get_groups()
  assert_eq(b_groups[0].group_name, "HOH")
  assert_eq(b_groups[0].group_number, 100)
  assert_eq(b_groups[0].get_num_atoms(), 1)
}

// ============================================================================
// Run-length encoding
// ============================================================================

///|
test "mmtf_encode_run_length_basic" {
  let data = [1, 1, 1, 2, 2, 3]
  let encoded = @src.mmtf_encode_run_length(data)
  // (1,3), (2,2), (3,1)
  assert_eq(encoded, [1, 3, 2, 2, 3, 1])
}

///|
test "mmtf_decode_run_length_basic" {
  let encoded = [1, 3, 2, 2, 3, 1]
  let decoded = @src.mmtf_decode_run_length(encoded)
  assert_eq(decoded, [1, 1, 1, 2, 2, 3])
}

///|
test "mmtf_run_length_round_trip" {
  let data = [5, 5, 5, 5, 7, 7, 9, 9, 9, 4]
  let encoded = @src.mmtf_encode_run_length(data)
  let decoded = @src.mmtf_decode_run_length(encoded)
  assert_eq(decoded, data)
}

///|
test "mmtf_encode_run_length_empty" {
  let encoded = @src.mmtf_encode_run_length([])
  assert_eq(encoded.length(), 0)
}

///|
test "mmtf_decode_run_length_empty" {
  let decoded = @src.mmtf_decode_run_length([])
  assert_eq(decoded.length(), 0)
}

///|
test "mmtf_encode_run_length_single" {
  let encoded = @src.mmtf_encode_run_length([42])
  assert_eq(encoded, [42, 1])
}

// ============================================================================
// Delta encoding
// ============================================================================

///|
test "mmtf_encode_delta_basic" {
  let data = [10, 11, 13, 16, 20]
  let encoded = @src.mmtf_encode_delta(data)
  // 10, 11-10=1, 13-11=2, 16-13=3, 20-16=4
  assert_eq(encoded, [10, 1, 2, 3, 4])
}

///|
test "mmtf_decode_delta_basic" {
  let encoded = [10, 1, 2, 3, 4]
  let decoded = @src.mmtf_decode_delta(encoded)
  assert_eq(decoded, [10, 11, 13, 16, 20])
}

///|
test "mmtf_delta_round_trip" {
  let data = [100, 95, 90, 85, 80, 75]
  let encoded = @src.mmtf_encode_delta(data)
  let decoded = @src.mmtf_decode_delta(encoded)
  assert_eq(decoded, data)
}

///|
test "mmtf_encode_delta_empty" {
  let encoded = @src.mmtf_encode_delta([])
  assert_eq(encoded.length(), 0)
}

///|
test "mmtf_decode_delta_empty" {
  let decoded = @src.mmtf_decode_delta([])
  assert_eq(decoded.length(), 0)
}

///|
test "mmtf_decode_delta_single" {
  let decoded = @src.mmtf_decode_delta([42])
  assert_eq(decoded, [42])
}

// ============================================================================
// Recursive index decode (currently a passthrough)
// ============================================================================

///|
test "mmtf_recursive_index_decode_passthrough" {
  let result = @src.mmtf_recursive_index_decode([1, 2, 3, 4, 5], 100)
  assert_eq(result, [1, 2, 3, 4, 5])
}

// ============================================================================
// mmtf_chain_sequence
// ============================================================================

///|
test "mmtf_chain_sequence_polymer" {
  let c = @src.MmtfChain::new(chain_id="A", chain_name="A")
  c.add_group(
    @src.MmtfGroup::new(
      group_name="ALA",
      group_number=1,
      single_letter_code="A",
    ),
  )
  c.add_group(
    @src.MmtfGroup::new(
      group_name="CYS",
      group_number=2,
      single_letter_code="C",
    ),
  )
  c.add_group(
    @src.MmtfGroup::new(
      group_name="ASP",
      group_number=3,
      single_letter_code="D",
    ),
  )
  assert_eq(@src.mmtf_chain_sequence(c), "ACD")
}

///|
test "mmtf_chain_sequence_skips_empty_and_unknown" {
  let c = @src.MmtfChain::new(chain_id="A", chain_name="A")
  c.add_group(
    @src.MmtfGroup::new(
      group_name="ALA",
      group_number=1,
      single_letter_code="A",
    ),
  )
  c.add_group(
    @src.MmtfGroup::new(
      group_name="HOH",
      group_number=2,
      single_letter_code="",
    ),
  )
  c.add_group(
    @src.MmtfGroup::new(
      group_name="UNK",
      group_number=3,
      single_letter_code="?",
    ),
  )
  c.add_group(
    @src.MmtfGroup::new(
      group_name="GLY",
      group_number=4,
      single_letter_code="G",
    ),
  )
  assert_eq(@src.mmtf_chain_sequence(c), "AG")
}

///|
test "mmtf_chain_sequence_empty" {
  let c = @src.MmtfChain::new()
  assert_eq(@src.mmtf_chain_sequence(c), "")
}

///|
test "mmtf_chain_sequence_sample" {
  let s = @src.mmtf_sample()
  let chains = s.get_chains()
  // Chain A: ALA (A), CYS (C) -> "AC"
  assert_eq(@src.mmtf_chain_sequence(chains[0]), "AC")
  // Chain B: HOH (no single letter code) -> ""
  assert_eq(@src.mmtf_chain_sequence(chains[1]), "")
}

// ============================================================================
// mmtf_chain_num_atoms / mmtf_model_num_atoms
// ============================================================================

///|
test "mmtf_chain_num_atoms" {
  let c = @src.MmtfChain::new(chain_id="A", chain_name="A")
  let g1 = @src.MmtfGroup::new(group_name="ALA", group_number=1)
  g1.add_atom(@src.MmtfAtom::new(name="N", serial=1, x=0.0, y=0.0, z=0.0))
  g1.add_atom(@src.MmtfAtom::new(name="CA", serial=2, x=1.0, y=0.0, z=0.0))
  let g2 = @src.MmtfGroup::new(group_name="GLY", group_number=2)
  g2.add_atom(@src.MmtfAtom::new(name="N", serial=3, x=2.0, y=0.0, z=0.0))
  c.add_group(g1)
  c.add_group(g2)
  assert_eq(@src.mmtf_chain_num_atoms(c), 3)
}

///|
test "mmtf_model_num_atoms" {
  let m = @src.MmtfModel::new(model_id=0)
  let c1 = @src.MmtfChain::new(chain_id="A", chain_name="A")
  let g1 = @src.MmtfGroup::new(group_name="ALA", group_number=1)
  g1.add_atom(@src.MmtfAtom::new(name="N", serial=1, x=0.0, y=0.0, z=0.0))
  c1.add_group(g1)
  let c2 = @src.MmtfChain::new(chain_id="B", chain_name="B")
  let g2 = @src.MmtfGroup::new(group_name="HOH", group_number=1)
  g2.add_atom(@src.MmtfAtom::new(name="O", serial=2, x=0.0, y=0.0, z=0.0))
  c2.add_group(g2)
  m.add_chain(c1)
  m.add_chain(c2)
  assert_eq(@src.mmtf_model_num_atoms(m), 2)
}

// ============================================================================
// mmtf_summary
// ============================================================================

///|
test "mmtf_summary_contains_key_fields" {
  let s = @src.mmtf_sample()
  let summary = @src.mmtf_summary(s)
  assert_true(summary.contains("MMTF Structure Summary:"))
  assert_true(summary.contains("1XXX"))
  assert_true(summary.contains("Sample protein structure"))
  assert_true(summary.contains("P 21 21 21"))
  assert_true(summary.contains("Models: 1"))
  assert_true(summary.contains("Chains: 2"))
  assert_true(summary.contains("Groups: 3"))
  assert_true(summary.contains("Atoms: 6"))
  assert_true(summary.contains("Bonds: 2"))
  assert_true(summary.contains("Entities: 1"))
}

///|
test "mmtf_summary_empty_structure" {
  let s = @src.MmtfStructure::new()
  let summary = @src.mmtf_summary(s)
  assert_true(summary.contains("Models: 0"))
  assert_true(summary.contains("Chains: 0"))
  assert_true(summary.contains("Groups: 0"))
  assert_true(summary.contains("Atoms: 0"))
}

// ============================================================================
// mmtf_write_text / mmtf_read_text round-trip
// ============================================================================

///|
test "mmtf_write_text_starts_with_mmtf" {
  let s = @src.mmtf_sample()
  let text = @src.mmtf_write_text(s)
  assert_true(text.has_prefix("MMTF\n"))
  assert_true(text.contains("structureId=1XXX"))
  assert_true(text.contains("title=Sample protein structure"))
  assert_true(text.contains("spaceGroup=P 21 21 21"))
  assert_true(text.contains("MODEL 0"))
  assert_true(text.contains("CHAIN 0 A A"))
  assert_true(text.contains("GROUP ALA 1 A alpha_helix"))
  assert_true(text.contains("ATOM N 1"))
}

///|
test "mmtf_read_text_round_trip_basic" {
  let s = @src.mmtf_sample()
  let text = @src.mmtf_write_text(s)
  let parsed = @src.mmtf_read_text(text)
  assert_eq(parsed.structure_id, "1XXX")
  assert_eq(parsed.total_num_models, 1)
  assert_eq(parsed.total_num_chains, 2)
  assert_eq(parsed.total_num_groups, 4)
  assert_eq(parsed.total_num_atoms, 10)
  assert_eq(parsed.total_num_bonds, 2)
}

///|
test "mmtf_read_text_round_trip_header" {
  let s = @src.mmtf_sample()
  let text = @src.mmtf_write_text(s)
  let parsed = @src.mmtf_read_text(text)
  let h = parsed.get_header()
  assert_eq(h.title, "Sample protein structure")
  assert_eq(h.deposition_date, "2024-01-01")
  assert_eq(h.release_date, "2024-06-01")
  assert_true((h.resolution - 1.5).abs() < 0.001)
  assert_true((h.r_work - 0.18).abs() < 0.001)
  assert_true((h.r_free - 0.22).abs() < 0.001)
}

///|
test "mmtf_read_text_round_trip_xtal" {
  let s = @src.mmtf_sample()
  let text = @src.mmtf_write_text(s)
  let parsed = @src.mmtf_read_text(text)
  let x = parsed.get_xtal_info()
  assert_eq(x.space_group, "P 21 21 21")
}

///|
test "mmtf_read_text_round_trip_hierarchy" {
  let s = @src.mmtf_sample()
  let text = @src.mmtf_write_text(s)
  let parsed = @src.mmtf_read_text(text)
  assert_eq(parsed.get_num_models(), 1)
  assert_eq(parsed.count_chains(), 2)
  assert_eq(parsed.count_groups(), 3)
  assert_eq(parsed.count_atoms(), 6)
  let models = parsed.get_models()
  let model = models[0]
  assert_eq(model.model_id, 0)
  let chains = model.get_chains()
  assert_eq(chains[0].chain_id, "A")
  assert_eq(chains[0].chain_name, "A")
  assert_eq(chains[0].get_num_groups(), 2)
  let groups = chains[0].get_groups()
  assert_eq(groups[0].group_name, "ALA")
  assert_eq(groups[0].group_number, 1)
  assert_eq(groups[0].single_letter_code, "A")
  assert_eq(groups[0].sec_struct.to_string(), "alpha_helix")
  assert_eq(groups[0].get_num_atoms(), 3)
  let atoms = groups[0].get_atoms()
  assert_eq(atoms[0].name, "N")
  assert_eq(atoms[0].serial, 1)
  assert_eq(atoms[1].name, "CA")
  assert_eq(atoms[2].name, "C")
  assert_eq(chains[1].chain_id, "B")
  assert_eq(chains[1].get_num_groups(), 1)
}

///|
test "mmtf_read_text_round_trip_bonds" {
  let s = @src.mmtf_sample()
  let text = @src.mmtf_write_text(s)
  let parsed = @src.mmtf_read_text(text)
  assert_eq(parsed.bonds.length(), 2)
  assert_eq(parsed.bonds[0].atom_index1, 0)
  assert_eq(parsed.bonds[0].atom_index2, 1)
  assert_eq(parsed.bonds[0].bond_order, 1)
  assert_eq(parsed.bonds[1].atom_index1, 1)
  assert_eq(parsed.bonds[1].atom_index2, 2)
}

///|
test "mmtf_read_text_empty" {
  let parsed = @src.mmtf_read_text("")
  assert_eq(parsed.structure_id, "")
  assert_eq(parsed.get_num_models(), 0)
}

///|
test "mmtf_write_text_empty_structure" {
  let s = @src.MmtfStructure::new()
  let text = @src.mmtf_write_text(s)
  assert_true(text.has_prefix("MMTF\n"))
  assert_true(text.contains("structureId="))
  assert_true(text.contains("totalNumModels=0"))
}

// ============================================================================
// mmtf_to_pdb_structure
// ============================================================================

///|
test "mmtf_to_pdb_structure_basic" {
  let s = @src.mmtf_sample()
  let pdb = @src.mmtf_to_pdb_structure(s)
  assert_eq(pdb.get_id(), "1XXX")
  assert_eq(pdb.get_num_models(), 1)
}

///|
test "mmtf_to_pdb_structure_counts" {
  let s = @src.mmtf_sample()
  let pdb = @src.mmtf_to_pdb_structure(s)
  // Sample has 2 chains, 3 residues (groups), 6 atoms
  assert_eq(pdb.get_num_chains(), 2)
  assert_eq(pdb.get_num_residues(), 3)
  assert_eq(pdb.get_num_atoms(), 6)
}

///|
test "mmtf_to_pdb_structure_chain_ids" {
  let s = @src.mmtf_sample()
  let pdb = @src.mmtf_to_pdb_structure(s)
  let chains = pdb.get_chains()
  assert_eq(chains[0].get_id(), 'A')
  assert_eq(chains[1].get_id(), 'B')
}

///|
test "mmtf_to_pdb_structure_residues" {
  let s = @src.mmtf_sample()
  let pdb = @src.mmtf_to_pdb_structure(s)
  let residues = pdb.get_residues()
  assert_eq(residues[0].resname, "ALA")
  assert_eq(residues[0].resseq, 1)
  assert_eq(residues[0].chainid, 'A')
  assert_eq(residues[0].get_num_atoms(), 3)
  assert_eq(residues[1].resname, "CYS")
  assert_eq(residues[1].resseq, 2)
  assert_eq(residues[2].resname, "HOH")
  assert_eq(residues[2].resseq, 100)
}

///|
test "mmtf_to_pdb_structure_atoms" {
  let s = @src.mmtf_sample()
  let pdb = @src.mmtf_to_pdb_structure(s)
  let atoms = pdb.get_atoms()
  assert_eq(atoms.length(), 6)
  assert_eq(atoms[0].name, "N")
  let coord = atoms[0].get_coord()
  assert_true((coord.x - 10.0).abs() < 0.001)
  assert_true((coord.y - 20.0).abs() < 0.001)
  assert_true((coord.z - 30.0).abs() < 0.001)
}

///|
test "mmtf_to_pdb_structure_het_residues" {
  let s = @src.mmtf_sample()
  let pdb = @src.mmtf_to_pdb_structure(s)
  let residues = pdb.get_residues()
  // ALA and CYS are not HET
  assert_false(residues[0].is_het())
  assert_false(residues[1].is_het())
  // HOH is a HET residue
  assert_true(residues[2].is_het())
}
