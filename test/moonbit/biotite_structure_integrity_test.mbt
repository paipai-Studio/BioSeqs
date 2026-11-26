///|
/// Tests for biotite_structure_integrity.mbt (official
/// `biotite.structure.integrity` semantics: continuity checks,
/// duplicate detection and the in-box check, ported from the official
/// integrity.py and test_integrity.py).
fn big_atom(
  chain : Char,
  res_id : Int32,
  res_name : String,
  name : String,
  x : Double,
) -> @src.Atom {
  @src.Atom::new(
    name~,
    coord=@src.Vector3::new(x, 0.0, 0.0),
    resname=res_name,
    chainid=chain,
    resseq=res_id,
    element="C",
  )
}

///|
test "integrity atom id continuity official semantics" {
  // Official _check_continuity: the returned index points to the atom
  // after a gap; decrements count as discontinuities too.
  assert_eq(@src.biotite_check_atom_id_continuity([1, 2, 3, 5, 6]), [3])
  assert_eq(@src.biotite_check_atom_id_continuity([1, 2, 3]), [])
  assert_eq(@src.biotite_check_atom_id_continuity([3, 2, 1]), [1, 2])
  // Repeated IDs (diff == 0) are continuous, e.g. altloc duplicates.
  assert_eq(@src.biotite_check_atom_id_continuity([1, 2, 2, 3]), [])
}

///|
test "integrity res id continuity" {
  let arr = @src.AtomArray::from_atoms([
    big_atom('A', 1, "ALA", "N", 0.0),
    big_atom('A', 2, "ALA", "CA", 1.0),
    big_atom('A', 3, "ALA", "C", 2.0),
    big_atom('A', 5, "ALA", "N", 3.0),
  ])
  let discon = @src.biotite_check_res_id_continuity(arr)
  assert_eq(discon.length(), 1)
  assert_eq(discon[0], 3)
}

///|
test "integrity linear continuity official bounds" {
  // Atom chain with distances 1.5, 1.5, 5.0, 1.5 along x: the atom
  // after the too-long bond (index 3) is flagged with the default
  // bounds [1.2, 1.8]. Note that bonds *shorter* than min_len count
  // as discontinuities too.
  let arr = @src.AtomArray::from_atoms([
    big_atom('A', 1, "ALA", "C1", 0.0),
    big_atom('A', 1, "ALA", "C2", 1.5),
    big_atom('A', 1, "ALA", "C4", 3.0),
    big_atom('A', 1, "ALA", "C3", 8.0),
    big_atom('A', 1, "ALA", "C5", 9.5),
  ])
  let discon = @src.biotite_check_linear_continuity(arr)
  assert_eq(discon.length(), 1)
  assert_eq(discon[0], 3)
  // The 1.0-spaced chain is entirely outside [1.2, 1.8]: every atom
  // except the first is flagged.
  let tight = @src.AtomArray::from_atoms([
    big_atom('A', 1, "ALA", "C1", 0.0),
    big_atom('A', 1, "ALA", "C2", 1.0),
    big_atom('A', 1, "ALA", "C4", 2.0),
    big_atom('A', 1, "ALA", "C3", 3.0),
  ])
  let short_bonds = @src.biotite_check_linear_continuity(tight)
  assert_eq(short_bonds.length(), 3)
  assert_eq(short_bonds[0], 1)
  // A custom window flags everything outside it.
  let wide = @src.biotite_check_linear_continuity(arr, min_len=0.5, max_len=6.0)
  assert_eq(wide.length(), 0)
}

///|
test "integrity backbone continuity peptide and phosphate" {
  // Peptide backbone atoms ("N", "CA", "C" of amino acids) form a
  // subsequence; the jump from 3.0 A to 8.0 A between two backbone
  // atoms is flagged, while non-backbone atoms (e.g. "O" is not in
  // the official peptide set) are skipped.
  let arr = @src.AtomArray::from_atoms([
    big_atom('A', 1, "ALA", "N", 0.0),
    big_atom('A', 1, "ALA", "CA", 1.5),
    big_atom('A', 1, "ALA", "C", 3.0),
    big_atom('A', 1, "ALA", "O", 3.5),
    big_atom('A', 2, "ALA", "N", 8.0),
  ])
  let discon = @src.biotite_check_backbone_continuity(arr)
  assert_eq(discon.length(), 1)
  assert_eq(discon[0], 4)
  // Phosphate backbone: "P" and "O5'" of a nucleotide participate.
  let na = @src.AtomArray::from_atoms([
    @src.Atom::new(
      name="P",
      coord=@src.Vector3::new(0.0, 0.0, 0.0),
      resname="DA",
      chainid='A',
      resseq=1,
      element="P",
    ),
    @src.Atom::new(
      name="O5'",
      coord=@src.Vector3::new(1.6, 0.0, 0.0),
      resname="DA",
      chainid='A',
      resseq=1,
      element="O",
    ),
    @src.Atom::new(
      name="O5'",
      coord=@src.Vector3::new(9.0, 0.0, 0.0),
      resname="DT",
      chainid='A',
      resseq=2,
      element="O",
    ),
  ])
  let na_discon = @src.biotite_check_backbone_continuity(na)
  assert_eq(na_discon.length(), 1)
  assert_eq(na_discon[0], 2)
}

///|
test "integrity duplicate atoms annotations only" {
  // Duplicate = all annotations equal (coordinates may differ); the
  // first occurrence is not reported (official test_duplicate_atoms).
  let dup = @src.Atom::new(
    name="CA",
    coord=@src.Vector3::new(1.0, 2.0, 3.0),
    resname="ALA",
    chainid='A',
    resseq=1,
    element="C",
  )
  let arr = @src.AtomArray::from_atoms([
    @src.Atom::new(
      name="CA",
      coord=@src.Vector3::new(0.0, 0.0, 0.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="C",
    ),
    @src.Atom::new(
      name="CA",
      coord=@src.Vector3::new(0.0, 0.0, 0.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="C",
    ),
    dup,
  ])
  let discon = @src.biotite_check_duplicate_atoms(arr)
  assert_eq(discon.length(), 2)
  assert_eq(discon[0], 1)
  assert_eq(discon[1], 2)
  // A different atom name breaks the annotation equality.
  let arr2 = @src.AtomArray::from_atoms([
    @src.Atom::new(
      name="CA",
      coord=@src.Vector3::new(0.0, 0.0, 0.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="C",
    ),
    @src.Atom::new(
      name="N",
      coord=@src.Vector3::new(0.0, 0.0, 0.0),
      resname="ALA",
      chainid='A',
      resseq=1,
      element="N",
    ),
  ])
  assert_eq(@src.biotite_check_duplicate_atoms(arr2).length(), 0)
}

///|
test "integrity in box returns inside indices" {
  // Official implementation semantics: indices of atoms with
  // fractional coordinates within [0, 1) on all axes.
  let arr = @src.AtomArray::from_atoms([
    big_atom('A', 1, "ALA", "N", 5.0),
    big_atom('A', 2, "ALA", "CA", 15.0),
  ])
  let box = @src.BiotiteBox::from_vectors([10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [
    0.0, 0.0, 10.0,
  ])
  let inside = @src.biotite_check_in_box(arr, box=Some(box))
  assert_eq(inside.length(), 1)
  assert_eq(inside[0], 0)
  // No box: official raises TypeError, here the result is empty.
  assert_eq(@src.biotite_check_in_box(arr).length(), 0)
}

///|
test "info vdw radius single mantina table" {
  // Official vdw_radius_single: Mantina et al. 2009 main group values
  // and the RDKit transition metal subset.
  match @src.biotite_vdw_radius_single("C") {
    Some(r) => assert_true(@moonbit.assert_num_close(r, 1.70, 1.0e-12))
    None => fail("C must have a vdW radius")
  }
  match @src.biotite_vdw_radius_single("BR") {
    Some(r) => assert_true(@moonbit.assert_num_close(r, 1.83, 1.0e-12))
    None => fail("BR must have a vdW radius")
  }
  // Case-insensitive lookup (official upper-cases the input).
  match @src.biotite_vdw_radius_single("fe") {
    Some(r) => assert_true(@moonbit.assert_num_close(r, 2.05, 1.0e-12))
    None => fail("FE must have a vdW radius")
  }
  match @src.biotite_vdw_radius_single("Xx") {
    Some(_) => fail("unknown elements must have no radius")
    None => ()
  }
}
