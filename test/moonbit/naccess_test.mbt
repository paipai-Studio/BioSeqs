///|
/// Tests for Bio.PDB.NACCESS module.
///
/// NACCESS (Atomic Solvent Accessible Area Calculations) parses the .rsa
/// (residue-level) and .asa (atom-level) output files produced by the
/// NACCESS program.

///|
/// Tolerance for floating-point comparisons.
let naccess_eps : Double = 1.0e-6

// ============================================================================
// NaccessAtom::new - construction and field access
// ============================================================================

///|
test "naccess_atom_new_basic" {
  let a = @src.NaccessAtom::new(
    atom_name="CA",
    atom_num=42,
    res_name="ALA",
    res_num=7,
    chain_id="A",
    abs_asa=15.5,
    rel_asa=30.2,
    element="C",
  )
  assert_eq(a.atom_name, "CA")
  assert_eq(a.atom_num, 42)
  assert_eq(a.res_name, "ALA")
  assert_eq(a.res_num, 7)
  assert_eq(a.chain_id, "A")
  assert_true((a.abs_asa - 15.5).abs() < naccess_eps)
  assert_true((a.rel_asa - 30.2).abs() < naccess_eps)
  assert_eq(a.element, "C")
}

///|
test "naccess_atom_new_defaults" {
  let a = @src.NaccessAtom::new(
    atom_name="N",
    atom_num=1,
    res_name="MET",
    res_num=1,
    chain_id="A",
    abs_asa=10.0,
  )
  // rel_asa defaults to 0.0, element defaults to "".
  assert_true(a.rel_asa == 0.0)
  assert_eq(a.element, "")
}

// ============================================================================
// NaccessResidue::new - construction and field access
// ============================================================================

///|
test "naccess_residue_new_basic" {
  let r = @src.NaccessResidue::new(
    res_name="MET",
    res_num=1,
    chain_id="A",
    abs_all=45.3,
    rel_all=22.1,
    abs_side_chain=30.1,
    rel_side_chain=18.5,
    abs_main_chain=15.2,
    rel_main_chain=35.2,
    non_polar_abs=25.1,
    non_polar_rel=20.1,
    polar_abs=20.2,
    polar_rel=25.1,
  )
  assert_eq(r.res_name, "MET")
  assert_eq(r.res_num, 1)
  assert_eq(r.chain_id, "A")
  assert_true((r.abs_all - 45.3).abs() < naccess_eps)
  assert_true((r.rel_all - 22.1).abs() < naccess_eps)
  assert_true((r.abs_side_chain - 30.1).abs() < naccess_eps)
  assert_true((r.rel_side_chain - 18.5).abs() < naccess_eps)
  assert_true((r.abs_main_chain - 15.2).abs() < naccess_eps)
  assert_true((r.rel_main_chain - 35.2).abs() < naccess_eps)
  assert_true((r.non_polar_abs - 25.1).abs() < naccess_eps)
  assert_true((r.non_polar_rel - 20.1).abs() < naccess_eps)
  assert_true((r.polar_abs - 20.2).abs() < naccess_eps)
  assert_true((r.polar_rel - 25.1).abs() < naccess_eps)
}

///|
test "naccess_residue_new_defaults" {
  let r = @src.NaccessResidue::new(res_name="GLY", res_num=5, chain_id="B")
  assert_eq(r.res_name, "GLY")
  assert_eq(r.res_num, 5)
  assert_eq(r.chain_id, "B")
  assert_true(r.abs_all == 0.0)
  assert_true(r.rel_all == 0.0)
  assert_true(r.abs_side_chain == 0.0)
  assert_true(r.rel_side_chain == 0.0)
  assert_true(r.abs_main_chain == 0.0)
  assert_true(r.rel_main_chain == 0.0)
  assert_true(r.non_polar_abs == 0.0)
  assert_true(r.non_polar_rel == 0.0)
  assert_true(r.polar_abs == 0.0)
  assert_true(r.polar_rel == 0.0)
}

// ============================================================================
// NaccessResult::new - empty construction
// ============================================================================

///|
test "naccess_result_new_empty" {
  let r = @src.NaccessResult::new()
  assert_eq(r.get_num_residues(), 0)
  assert_eq(r.get_num_atoms(), 0)
  assert_eq(r.structure_id, "")
  assert_true(r.total_abs_asa == 0.0)
  assert_true(r.total_rel_asa == 0.0)
  assert_eq(r.chain_totals.length(), 0)
}

// ============================================================================
// NaccessResult - add_residue, add_atom, get_residues, get_atoms
// ============================================================================

///|
test "naccess_result_add_residue" {
  let r = @src.NaccessResult::new()
  assert_eq(r.get_num_residues(), 0)
  r.add_residue(@src.NaccessResidue::new(res_name="ALA", res_num=1, chain_id="A"))
  assert_eq(r.get_num_residues(), 1)
  r.add_residue(@src.NaccessResidue::new(res_name="GLY", res_num=2, chain_id="A"))
  assert_eq(r.get_num_residues(), 2)
}

///|
test "naccess_result_add_atom" {
  let r = @src.NaccessResult::new()
  assert_eq(r.get_num_atoms(), 0)
  r.add_atom(
    @src.NaccessAtom::new(
      atom_name="CA",
      atom_num=1,
      res_name="ALA",
      res_num=1,
      chain_id="A",
      abs_asa=10.0,
    ),
  )
  assert_eq(r.get_num_atoms(), 1)
  r.add_atom(
    @src.NaccessAtom::new(
      atom_name="N",
      atom_num=2,
      res_name="ALA",
      res_num=1,
      chain_id="A",
      abs_asa=5.0,
    ),
  )
  assert_eq(r.get_num_atoms(), 2)
}

///|
test "naccess_result_get_residues" {
  let r = @src.NaccessResult::new()
  r.add_residue(@src.NaccessResidue::new(res_name="MET", res_num=1, chain_id="A"))
  r.add_residue(@src.NaccessResidue::new(res_name="ALA", res_num=2, chain_id="A"))
  let residues = r.get_residues()
  assert_eq(residues.length(), 2)
  assert_eq(residues[0].res_name, "MET")
  assert_eq(residues[1].res_name, "ALA")
}

///|
test "naccess_result_get_atoms" {
  let r = @src.NaccessResult::new()
  r.add_atom(
    @src.NaccessAtom::new(
      atom_name="N",
      atom_num=1,
      res_name="MET",
      res_num=1,
      chain_id="A",
      abs_asa=15.0,
    ),
  )
  r.add_atom(
    @src.NaccessAtom::new(
      atom_name="CA",
      atom_num=2,
      res_name="MET",
      res_num=1,
      chain_id="A",
      abs_asa=10.0,
    ),
  )
  let atoms = r.get_atoms()
  assert_eq(atoms.length(), 2)
  assert_eq(atoms[0].atom_name, "N")
  assert_eq(atoms[1].atom_name, "CA")
}

///|
test "naccess_result_get_num_residues" {
  let r = @src.NaccessResult::new()
  assert_eq(r.get_num_residues(), 0)
  for i in 0..<5 {
    r.add_residue(@src.NaccessResidue::new(res_name="ALA", res_num=i, chain_id="A"))
  }
  assert_eq(r.get_num_residues(), 5)
}

///|
test "naccess_result_get_num_atoms" {
  let r = @src.NaccessResult::new()
  assert_eq(r.get_num_atoms(), 0)
  for i in 0..<3 {
    r.add_atom(
      @src.NaccessAtom::new(
        atom_name="CA",
        atom_num=i,
        res_name="ALA",
        res_num=1,
        chain_id="A",
        abs_asa=1.0,
      ),
    )
  }
  assert_eq(r.get_num_atoms(), 3)
}

// ============================================================================
// naccess_parse_rsa - parse RSA format text
// ============================================================================

///|
test "naccess_parse_rsa_count" {
  let residues = @src.naccess_parse_rsa(@src.naccess_sample_rsa())
  // MET, ALA, CYS, HOH - header should be skipped.
  assert_eq(residues.length(), 4)
}

///|
test "naccess_parse_rsa_met_fields" {
  let residues = @src.naccess_parse_rsa(@src.naccess_sample_rsa())
  let r = residues[0]
  assert_eq(r.res_name, "MET")
  assert_eq(r.res_num, 1)
  assert_eq(r.chain_id, "A")
  assert_true((r.abs_all - 45.3).abs() < naccess_eps)
  assert_true((r.rel_all - 22.1).abs() < naccess_eps)
  assert_true((r.abs_side_chain - 30.1).abs() < naccess_eps)
  assert_true((r.rel_side_chain - 18.5).abs() < naccess_eps)
  assert_true((r.abs_main_chain - 15.2).abs() < naccess_eps)
  assert_true((r.rel_main_chain - 35.2).abs() < naccess_eps)
  assert_true((r.non_polar_abs - 25.1).abs() < naccess_eps)
  assert_true((r.non_polar_rel - 20.1).abs() < naccess_eps)
  assert_true((r.polar_abs - 20.2).abs() < naccess_eps)
  assert_true((r.polar_rel - 25.1).abs() < naccess_eps)
}

///|
test "naccess_parse_rsa_ala_fields" {
  let residues = @src.naccess_parse_rsa(@src.naccess_sample_rsa())
  let r = residues[1]
  assert_eq(r.res_name, "ALA")
  assert_eq(r.res_num, 2)
  assert_eq(r.chain_id, "A")
  assert_true((r.abs_all - 20.1).abs() < naccess_eps)
  assert_true((r.rel_all - 17.8).abs() < naccess_eps)
  assert_true((r.abs_side_chain - 5.1).abs() < naccess_eps)
  assert_true((r.rel_side_chain - 8.2).abs() < naccess_eps)
}

///|
test "naccess_parse_rsa_cys_fields" {
  let residues = @src.naccess_parse_rsa(@src.naccess_sample_rsa())
  let r = residues[2]
  assert_eq(r.res_name, "CYS")
  assert_eq(r.res_num, 3)
  assert_eq(r.chain_id, "A")
  assert_true((r.abs_all - 5.2).abs() < naccess_eps)
  assert_true((r.rel_all - 3.7).abs() < naccess_eps)
}

///|
test "naccess_parse_rsa_hoh_fields" {
  let residues = @src.naccess_parse_rsa(@src.naccess_sample_rsa())
  let r = residues[3]
  assert_eq(r.res_name, "HOH")
  assert_eq(r.res_num, 1)
  assert_eq(r.chain_id, "B")
  assert_true((r.abs_all - 25.0).abs() < naccess_eps)
  assert_true(r.rel_all == 0.0)
}

///|
test "naccess_parse_rsa_skips_header" {
  // Header line "RES NUM CHN ..." must not be parsed as a residue.
  let residues = @src.naccess_parse_rsa(@src.naccess_sample_rsa())
  for r in residues {
    assert_false(r.res_name == "NUM")
  }
}

///|
test "naccess_parse_rsa_empty" {
  let residues = @src.naccess_parse_rsa("")
  assert_eq(residues.length(), 0)
}

///|
test "naccess_parse_rsa_header_only" {
  let residues = @src.naccess_parse_rsa("RES NUM CHN ABS_ALL REL_ALL\n")
  assert_eq(residues.length(), 0)
}

// ============================================================================
// naccess_parse_asa - parse ASA format text
// ============================================================================

///|
test "naccess_parse_asa_count" {
  let atoms = @src.naccess_parse_asa(@src.naccess_sample_asa())
  assert_eq(atoms.length(), 6)
}

///|
test "naccess_parse_asa_first_atom" {
  let atoms = @src.naccess_parse_asa(@src.naccess_sample_asa())
  let a = atoms[0]
  assert_eq(a.atom_name, "N")
  assert_eq(a.atom_num, 1)
  assert_eq(a.res_name, "MET")
  assert_eq(a.chain_id, "A")
  assert_true((a.abs_asa - 15.2).abs() < naccess_eps)
  assert_true((a.rel_asa - 35.2).abs() < naccess_eps)
}

///|
test "naccess_parse_asa_atom_names" {
  let atoms = @src.naccess_parse_asa(@src.naccess_sample_asa())
  assert_eq(atoms[0].atom_name, "N")
  assert_eq(atoms[1].atom_name, "CA")
  assert_eq(atoms[2].atom_name, "N")
  assert_eq(atoms[3].atom_name, "CA")
  assert_eq(atoms[4].atom_name, "N")
  assert_eq(atoms[5].atom_name, "O")
}

///|
test "naccess_parse_asa_residues" {
  let atoms = @src.naccess_parse_asa(@src.naccess_sample_asa())
  assert_eq(atoms[0].res_name, "MET")
  assert_eq(atoms[1].res_name, "MET")
  assert_eq(atoms[2].res_name, "ALA")
  assert_eq(atoms[3].res_name, "ALA")
  assert_eq(atoms[4].res_name, "CYS")
  assert_eq(atoms[5].res_name, "HOH")
}

///|
test "naccess_parse_asa_chain_ids" {
  let atoms = @src.naccess_parse_asa(@src.naccess_sample_asa())
  assert_eq(atoms[0].chain_id, "A")
  assert_eq(atoms[1].chain_id, "A")
  assert_eq(atoms[2].chain_id, "A")
  assert_eq(atoms[3].chain_id, "A")
  assert_eq(atoms[4].chain_id, "A")
  assert_eq(atoms[5].chain_id, "B")
}

///|
test "naccess_parse_asa_last_atom" {
  let atoms = @src.naccess_parse_asa(@src.naccess_sample_asa())
  let a = atoms[5]
  assert_eq(a.atom_name, "O")
  assert_eq(a.atom_num, 6)
  assert_eq(a.res_name, "HOH")
  assert_eq(a.chain_id, "B")
  assert_true((a.abs_asa - 25.0).abs() < naccess_eps)
  assert_true(a.rel_asa == 0.0)
}

///|
test "naccess_parse_asa_skips_header" {
  let atoms = @src.naccess_parse_asa(@src.naccess_sample_asa())
  for a in atoms {
    assert_false(a.atom_name == "ATOM_NAME")
  }
}

///|
test "naccess_parse_asa_empty" {
  let atoms = @src.naccess_parse_asa("")
  assert_eq(atoms.length(), 0)
}

///|
test "naccess_parse_asa_header_only" {
  let atoms = @src.naccess_parse_asa("ATM NUM RES CHN NAME ABS REL\n")
  assert_eq(atoms.length(), 0)
}

// ============================================================================
// naccess_parse - combined parse, verify totals
// ============================================================================

///|
test "naccess_parse_combined_counts" {
  let result = @src.naccess_parse(@src.naccess_sample_rsa(), @src.naccess_sample_asa())
  assert_eq(result.get_num_residues(), 4)
  assert_eq(result.get_num_atoms(), 6)
}

///|
test "naccess_parse_total_abs_asa" {
  let result = @src.naccess_parse(@src.naccess_sample_rsa(), @src.naccess_sample_asa())
  // 45.3 + 20.1 + 5.2 + 25.0 = 95.6
  assert_true((result.total_abs_asa - 95.6).abs() < naccess_eps)
}

///|
test "naccess_parse_chain_totals" {
  let result = @src.naccess_parse(@src.naccess_sample_rsa(), @src.naccess_sample_asa())
  assert_eq(result.chain_totals.length(), 2)
  // Chain A: 45.3 + 20.1 + 5.2 = 70.6
  assert_true((result.chain_totals.get("A").unwrap() - 70.6).abs() < naccess_eps)
  // Chain B: 25.0
  assert_true((result.chain_totals.get("B").unwrap() - 25.0).abs() < naccess_eps)
}

///|
test "naccess_parse_empty_strings" {
  let result = @src.naccess_parse("", "")
  assert_eq(result.get_num_residues(), 0)
  assert_eq(result.get_num_atoms(), 0)
  assert_true(result.total_abs_asa == 0.0)
  assert_eq(result.chain_totals.length(), 0)
}

// ============================================================================
// naccess_sample_rsa / naccess_sample_asa / naccess_sample
// ============================================================================

///|
test "naccess_sample_rsa_has_res_prefix" {
  let s = @src.naccess_sample_rsa()
  // Data lines start with "RES ".
  assert_true(s.contains("RES MET"))
  assert_true(s.contains("RES ALA"))
  assert_true(s.contains("RES CYS"))
  assert_true(s.contains("RES HOH"))
}

///|
test "naccess_sample_asa_has_atm_prefix" {
  let s = @src.naccess_sample_asa()
  assert_true(s.contains("ATM"))
  assert_true(s.contains("MET"))
  assert_true(s.contains("ALA"))
  assert_true(s.contains("CYS"))
  assert_true(s.contains("HOH"))
}

///|
test "naccess_sample_result" {
  let result = @src.naccess_sample()
  assert_eq(result.get_num_residues(), 4)
  assert_eq(result.get_num_atoms(), 6)
  assert_true((result.total_abs_asa - 95.6).abs() < naccess_eps)
}

// ============================================================================
// naccess_exposed_residues / naccess_buried_residues - threshold filtering
// ============================================================================

///|
test "naccess_exposed_residues" {
  let result = @src.naccess_sample()
  // rel_all > 20: only MET (22.1)
  let exposed = @src.naccess_exposed_residues(result, 20.0)
  assert_eq(exposed.length(), 1)
  assert_eq(exposed[0].res_name, "MET")
}

///|
test "naccess_buried_residues" {
  let result = @src.naccess_sample()
  // rel_all < 20: ALA (17.8), CYS (3.7), HOH (0.0)
  let buried = @src.naccess_buried_residues(result, 20.0)
  assert_eq(buried.length(), 3)
  assert_eq(buried[0].res_name, "ALA")
  assert_eq(buried[1].res_name, "CYS")
  assert_eq(buried[2].res_name, "HOH")
}

///|
test "naccess_exposed_residues_threshold_zero" {
  let result = @src.naccess_sample()
  // rel_all > 0: MET, ALA, CYS (HOH has rel_all == 0, not > 0)
  let exposed = @src.naccess_exposed_residues(result, 0.0)
  assert_eq(exposed.length(), 3)
}

///|
test "naccess_exposed_residues_empty_result" {
  let result = @src.NaccessResult::new()
  let exposed = @src.naccess_exposed_residues(result, 20.0)
  assert_eq(exposed.length(), 0)
}

// ============================================================================
// naccess_mean_rel_asa - average relative ASA
// ============================================================================

///|
test "naccess_mean_rel_asa" {
  let result = @src.naccess_sample()
  // (22.1 + 17.8 + 3.7 + 0.0) / 4 = 10.9
  assert_true((@src.naccess_mean_rel_asa(result) - 10.9).abs() < naccess_eps)
}

///|
test "naccess_mean_rel_asa_empty" {
  let result = @src.NaccessResult::new()
  assert_true(@src.naccess_mean_rel_asa(result) == 0.0)
}

// ============================================================================
// naccess_get_chain_residues - filter by chain
// ============================================================================

///|
test "naccess_get_chain_residues_a" {
  let result = @src.naccess_sample()
  let chain_a = @src.naccess_get_chain_residues(result, "A")
  assert_eq(chain_a.length(), 3)
  assert_eq(chain_a[0].res_name, "MET")
  assert_eq(chain_a[1].res_name, "ALA")
  assert_eq(chain_a[2].res_name, "CYS")
}

///|
test "naccess_get_chain_residues_b" {
  let result = @src.naccess_sample()
  let chain_b = @src.naccess_get_chain_residues(result, "B")
  assert_eq(chain_b.length(), 1)
  assert_eq(chain_b[0].res_name, "HOH")
}

///|
test "naccess_get_chain_residues_nonexistent" {
  let result = @src.naccess_sample()
  let chain_c = @src.naccess_get_chain_residues(result, "C")
  assert_eq(chain_c.length(), 0)
}

// ============================================================================
// naccess_count_exposed - count exposed
// ============================================================================

///|
test "naccess_count_exposed" {
  let result = @src.naccess_sample()
  assert_eq(@src.naccess_count_exposed(result, 20.0), 1)
  assert_eq(@src.naccess_count_exposed(result, 10.0), 2)
  assert_eq(@src.naccess_count_exposed(result, 0.0), 3)
  assert_eq(@src.naccess_count_exposed(result, 100.0), 0)
}

// ============================================================================
// naccess_chain_total_asa - chain total
// ============================================================================

///|
test "naccess_chain_total_asa_a" {
  let result = @src.naccess_sample()
  // 45.3 + 20.1 + 5.2 = 70.6
  assert_true((@src.naccess_chain_total_asa(result, "A") - 70.6).abs() < naccess_eps)
}

///|
test "naccess_chain_total_asa_b" {
  let result = @src.naccess_sample()
  assert_true((@src.naccess_chain_total_asa(result, "B") - 25.0).abs() < naccess_eps)
}

///|
test "naccess_chain_total_asa_nonexistent" {
  let result = @src.naccess_sample()
  assert_true(@src.naccess_chain_total_asa(result, "Z") == 0.0)
}

// ============================================================================
// naccess_is_exposed - classification
// ============================================================================

///|
test "naccess_is_exposed_true" {
  let r = @src.NaccessResidue::new(
    res_name="MET",
    res_num=1,
    chain_id="A",
    rel_all=22.1,
  )
  assert_true(@src.naccess_is_exposed(r, 20.0))
}

///|
test "naccess_is_exposed_false" {
  let r = @src.NaccessResidue::new(
    res_name="ALA",
    res_num=2,
    chain_id="A",
    rel_all=17.8,
  )
  assert_false(@src.naccess_is_exposed(r, 20.0))
}

///|
test "naccess_is_exposed_boundary" {
  // rel_all == threshold: not exposed (> is strict).
  let r = @src.NaccessResidue::new(
    res_name="X",
    res_num=1,
    chain_id="A",
    rel_all=20.0,
  )
  assert_false(@src.naccess_is_exposed(r, 20.0))
}

// ============================================================================
// naccess_standard_asa - reference values for all 20 amino acids
// ============================================================================

///|
test "naccess_standard_asa_ala" {
  assert_true(@src.naccess_standard_asa("ALA") == 113.0)
}

///|
test "naccess_standard_asa_all_twenty" {
  assert_true(@src.naccess_standard_asa("ALA") == 113.0)
  assert_true(@src.naccess_standard_asa("ARG") == 241.0)
  assert_true(@src.naccess_standard_asa("ASN") == 158.0)
  assert_true(@src.naccess_standard_asa("ASP") == 151.0)
  assert_true(@src.naccess_standard_asa("CYS") == 140.0)
  assert_true(@src.naccess_standard_asa("GLN") == 189.0)
  assert_true(@src.naccess_standard_asa("GLU") == 183.0)
  assert_true(@src.naccess_standard_asa("GLY") == 85.0)
  assert_true(@src.naccess_standard_asa("HIS") == 194.0)
  assert_true(@src.naccess_standard_asa("ILE") == 182.0)
  assert_true(@src.naccess_standard_asa("LEU") == 180.0)
  assert_true(@src.naccess_standard_asa("LYS") == 211.0)
  assert_true(@src.naccess_standard_asa("MET") == 204.0)
  assert_true(@src.naccess_standard_asa("PHE") == 218.0)
  assert_true(@src.naccess_standard_asa("PRO") == 143.0)
  assert_true(@src.naccess_standard_asa("SER") == 122.0)
  assert_true(@src.naccess_standard_asa("THR") == 146.0)
  assert_true(@src.naccess_standard_asa("TRP") == 259.0)
  assert_true(@src.naccess_standard_asa("TYR") == 229.0)
  assert_true(@src.naccess_standard_asa("VAL") == 160.0)
}

///|
test "naccess_standard_asa_unknown" {
  assert_true(@src.naccess_standard_asa("HOH") == 0.0)
  assert_true(@src.naccess_standard_asa("XXX") == 0.0)
  assert_true(@src.naccess_standard_asa("") == 0.0)
}

// ============================================================================
// naccess_residue_to_rsa / naccess_atom_to_asa - serialization
// ============================================================================

///|
test "naccess_residue_to_rsa_format" {
  let r = @src.NaccessResidue::new(
    res_name="MET",
    res_num=1,
    chain_id="A",
    abs_all=45.3,
    rel_all=22.1,
  )
  let s = @src.naccess_residue_to_rsa(r)
  assert_true(s.has_prefix("RES  MET"))
  assert_true(s.contains("  A  "))
  assert_true(s.contains("1"))
}

///|
test "naccess_residue_to_rsa_roundtrip" {
  let r = @src.NaccessResidue::new(
    res_name="MET",
    res_num=1,
    chain_id="A",
    abs_all=45.3,
    rel_all=22.1,
    abs_side_chain=30.1,
    rel_side_chain=18.5,
    abs_main_chain=15.2,
    rel_main_chain=35.2,
    non_polar_abs=25.1,
    non_polar_rel=20.1,
    polar_abs=20.2,
    polar_rel=25.1,
  )
  let s = @src.naccess_residue_to_rsa(r)
  let parsed = @src.naccess_parse_rsa(s)
  assert_eq(parsed.length(), 1)
  let p = parsed[0]
  assert_eq(p.res_name, "MET")
  assert_eq(p.res_num, 1)
  assert_eq(p.chain_id, "A")
  assert_true((p.abs_all - 45.3).abs() < naccess_eps)
  assert_true((p.rel_all - 22.1).abs() < naccess_eps)
  assert_true((p.abs_side_chain - 30.1).abs() < naccess_eps)
  assert_true((p.rel_side_chain - 18.5).abs() < naccess_eps)
  assert_true((p.abs_main_chain - 15.2).abs() < naccess_eps)
  assert_true((p.rel_main_chain - 35.2).abs() < naccess_eps)
}

///|
test "naccess_atom_to_asa_format" {
  let a = @src.NaccessAtom::new(
    atom_name="N",
    atom_num=1,
    res_name="MET",
    res_num=1,
    chain_id="A",
    abs_asa=15.2,
    rel_asa=35.2,
  )
  let s = @src.naccess_atom_to_asa(a)
  assert_true(s.has_prefix("ATM  1"))
  assert_true(s.contains("MET"))
  assert_true(s.contains("A"))
  assert_true(s.contains("N"))
}

///|
test "naccess_atom_to_asa_roundtrip" {
  let a = @src.NaccessAtom::new(
    atom_name="CA",
    atom_num=2,
    res_name="ALA",
    res_num=5,
    chain_id="A",
    abs_asa=10.5,
    rel_asa=22.3,
  )
  let s = @src.naccess_atom_to_asa(a)
  let parsed = @src.naccess_parse_asa(s)
  assert_eq(parsed.length(), 1)
  let p = parsed[0]
  assert_eq(p.atom_name, "CA")
  assert_eq(p.atom_num, 2)
  assert_eq(p.res_name, "ALA")
  assert_eq(p.chain_id, "A")
  assert_true((p.abs_asa - 10.5).abs() < naccess_eps)
  assert_true((p.rel_asa - 22.3).abs() < naccess_eps)
}

// ============================================================================
// naccess_write_rsa / naccess_write_asa - file serialization
// ============================================================================

///|
test "naccess_write_rsa_format" {
  let result = @src.naccess_sample()
  let s = @src.naccess_write_rsa(result)
  assert_true(s.has_prefix("NACCESS"))
  assert_true(s.contains("RES"))
  assert_true(s.contains("MET"))
  assert_true(s.contains("ALA"))
  assert_true(s.contains("CYS"))
  assert_true(s.contains("HOH"))
}

///|
test "naccess_write_rsa_roundtrip" {
  let result = @src.naccess_sample()
  let s = @src.naccess_write_rsa(result)
  let reparsed = @src.naccess_parse_rsa(s)
  assert_eq(reparsed.length(), result.get_num_residues())
  // Verify first residue round-trips.
  assert_eq(reparsed[0].res_name, "MET")
  assert_eq(reparsed[0].res_num, 1)
  assert_true((reparsed[0].abs_all - 45.3).abs() < naccess_eps)
}

///|
test "naccess_write_asa_format" {
  let result = @src.naccess_sample()
  let s = @src.naccess_write_asa(result)
  assert_true(s.has_prefix("NACCESS"))
  assert_true(s.contains("ATM"))
  assert_true(s.contains("MET"))
  assert_true(s.contains("ALA"))
  assert_true(s.contains("HOH"))
}

///|
test "naccess_write_asa_roundtrip" {
  let result = @src.naccess_sample()
  let s = @src.naccess_write_asa(result)
  let reparsed = @src.naccess_parse_asa(s)
  assert_eq(reparsed.length(), result.get_num_atoms())
  assert_eq(reparsed[0].atom_name, "N")
  assert_eq(reparsed[0].res_name, "MET")
  assert_true((reparsed[0].abs_asa - 15.2).abs() < naccess_eps)
}

// ============================================================================
// naccess_summary - summary string
// ============================================================================

///|
test "naccess_summary_content" {
  let result = @src.naccess_sample()
  let s = @src.naccess_summary(result)
  assert_true(s.has_prefix("NACCESS Result Summary:"))
  assert_true(s.contains("Residues: 4"))
  assert_true(s.contains("Atoms: 6"))
  assert_true(s.contains("Total ASA:"))
  assert_true(s.contains("Mean Rel ASA:"))
  assert_true(s.contains("Exposed (>20%): 1"))
  assert_true(s.contains("Buried"))
  assert_true(s.contains("Chains: 2"))
  assert_true(s.contains("Chain A:"))
  assert_true(s.contains("Chain B:"))
}

///|
test "naccess_summary_empty" {
  let result = @src.NaccessResult::new()
  let s = @src.naccess_summary(result)
  assert_true(s.has_prefix("NACCESS Result Summary:"))
  assert_true(s.contains("Residues: 0"))
  assert_true(s.contains("Atoms: 0"))
  assert_true(s.contains("Chains: 0"))
}
