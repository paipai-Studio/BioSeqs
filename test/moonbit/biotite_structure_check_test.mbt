///|
/// Tests for biotite_structure_check.mbt.
fn bsc_chk_atom(
  name : String,
  res_name : String,
  res_id : Int32,
  x : Double,
  y : Double,
  z : Double,
  element : String,
) -> @src.Atom {
  @src.Atom::new(
    name~,
    coord=@src.Vector3::new(x, y, z),
    resname=res_name,
    chainid='A',
    resseq=res_id,
    element~,
  )
}

///|
/// Two alanine residues with a sane backbone (C-N = 1.2 A, N-CA-C
/// angles around 105 degrees).
fn bsc_sane_backbone() -> @src.AtomArray {
  @src.AtomArray::from_atoms([
    bsc_chk_atom("N", "ALA", 1, 0.0, 0.0, 0.0, "N"),
    bsc_chk_atom("CA", "ALA", 1, 0.0, 1.5, 0.4, "C"),
    bsc_chk_atom("C", "ALA", 1, 0.0, 1.5, 1.9, "C"),
    bsc_chk_atom("N", "ALA", 2, 0.0, 1.5, 3.1, "N"),
    bsc_chk_atom("CA", "ALA", 2, 0.0, 2.9, 3.3, "C"),
    bsc_chk_atom("C", "ALA", 2, 0.0, 3.1, 4.8, "C"),
  ])
}

///|
test "backbone check accepts sane backbone" {
  let arr = bsc_sane_backbone()
  let violations = @src.biotite_check_backbone_structure(arr)
  assert_eq(violations.length(), 0)
}

///|
test "backbone check reports missing CA" {
  let arr = bsc_sane_backbone()
  // remove the CA of residue 2 (index 4)
  let mask = [true, true, true, true, false, true]
  let truncated = arr.filter(mask)
  let violations = @src.biotite_check_backbone_structure(truncated)
  let mut found_missing = false
  for v in violations {
    if v.kind == "missing_backbone_atom" {
      found_missing = true
    }
  }
  assert_true(found_missing)
}

///|
test "backbone check reports broken peptide bond" {
  let arr = @src.AtomArray::from_atoms([
    bsc_chk_atom("N", "ALA", 1, 0.0, 0.0, 0.0, "N"),
    bsc_chk_atom("CA", "ALA", 1, 0.0, 1.5, 0.4, "C"),
    bsc_chk_atom("C", "ALA", 1, 0.0, 1.5, 1.9, "C"),
    // N of residue 2 is 5 A away -> broken bond
    bsc_chk_atom("N", "ALA", 2, 0.0, 1.5, 6.9, "N"),
    bsc_chk_atom("CA", "ALA", 2, 0.0, 2.9, 7.1, "C"),
    bsc_chk_atom("C", "ALA", 2, 0.0, 3.1, 8.6, "C"),
  ])
  let violations = @src.biotite_check_backbone_structure(arr)
  let mut found_bond = false
  for v in violations {
    if v.kind == "invalid_peptide_bond" {
      found_bond = true
    }
  }
  assert_true(found_bond)
}

///|
test "backbone check reports bad N-CA-C angle" {
  // linear N-CA-C arrangement gives a 180 degree angle
  let arr = @src.AtomArray::from_atoms([
    bsc_chk_atom("N", "ALA", 1, 0.0, 0.0, 0.0, "N"),
    bsc_chk_atom("CA", "ALA", 1, 0.0, 1.5, 0.0, "C"),
    bsc_chk_atom("C", "ALA", 1, 0.0, 3.0, 0.0, "C"),
    bsc_chk_atom("N", "ALA", 2, 0.0, 4.2, 0.0, "N"),
    bsc_chk_atom("CA", "ALA", 2, 0.0, 5.7, 0.0, "C"),
    bsc_chk_atom("C", "ALA", 2, 0.0, 7.2, 0.0, "C"),
  ])
  let violations = @src.biotite_check_backbone_structure(arr)
  let mut found_angle = false
  for v in violations {
    if v.kind == "invalid_bond_angle" {
      found_angle = true
    }
  }
  assert_true(found_angle)
}

///|
test "steric clash detection" {
  let arr = @src.AtomArray::from_atoms([
    bsc_chk_atom("O1", "HOH", 1, 0.0, 0.0, 0.0, "O"),
    bsc_chk_atom("O2", "HOH", 2, 0.0, 0.0, 1.0, "O"),
    bsc_chk_atom("O3", "HOH", 3, 0.0, 0.0, 5.0, "O"),
  ])
  let check = @src.biotite_steric_clash_check(arr)
  let clashes = check.get_clashes()
  assert_eq(clashes.length(), 1)
  assert_eq(clashes[0].atom_i, 0)
  assert_eq(clashes[0].atom_j, 1)
  assert_true(@moonbit.assert_num_close(clashes[0].distance, 1.0, 1.0e-12))
  // mean distance = (1.52 + 1.52) / 2
  assert_true(
    @moonbit.assert_num_close(clashes[0].mean_distance, 1.52, 1.0e-12),
  )
  let indices = check.get_clash_indices()
  assert_eq(indices.length(), 2)
  assert_eq(indices[0], 0)
  assert_eq(indices[1], 1)
}

///|
test "steric clash respects stereochemistry" {
  let arr = @src.AtomArray::from_atoms([
    bsc_chk_atom("C1", "LIG", 1, 0.0, 0.0, 0.0, "C"),
    bsc_chk_atom("C2", "LIG", 1, 0.0, 0.0, 1.0, "C"),
  ])
  let bonds = @src.BondList::new(n_atoms=2).add(0, 1)
  // bonded pair: no clash reported when respecting stereochemistry
  let check_respect = @src.biotite_steric_clash_check(
    arr,
    bonds=Some(bonds),
    respect_stereochemistry=true,
  )
  assert_eq(check_respect.get_clashes().length(), 0)
  // without respect: 1.0 A < (1.7+1.7)*0.7 -> clash
  let check_plain = @src.biotite_steric_clash_check(arr)
  assert_eq(check_plain.get_clashes().length(), 1)
}

///|
test "vdw radius table" {
  assert_true(
    @moonbit.assert_num_close(@src.biotite_vdw_radius("H"), 1.10, 1.0e-12),
  )
  assert_true(
    @moonbit.assert_num_close(@src.biotite_vdw_radius("C"), 1.70, 1.0e-12),
  )
  assert_true(
    @moonbit.assert_num_close(@src.biotite_vdw_radius("S"), 1.80, 1.0e-12),
  )
  assert_true(
    @moonbit.assert_num_close(@src.biotite_vdw_radius("XX"), 1.50, 1.0e-12),
  )
}
