///|
/// Tests for biotite_structure_repair.mbt (official
/// `biotite.structure.repair` semantics, verified against the official
/// doctest values of repair.py).
fn brp_atom(
  chain : Char,
  res_id : Int32,
  name : String,
  x : Double,
) -> @src.Atom {
  @src.Atom::new(
    name~,
    coord=@src.Vector3::new(x, 0.0, 0.0),
    resname="ALA",
    chainid=chain,
    resseq=res_id,
    element="C",
  )
}

///|
test "repair create continuous res ids single chain" {
  // Official doctest: removing residue 5 from [1..20] and re-running
  // yields contiguous IDs; here with three atoms per residue and a
  // missing residue 3.
  let arr = @src.AtomArray::from_atoms([
    brp_atom('A', 1, "N", 0.0),
    brp_atom('A', 1, "CA", 1.0),
    brp_atom('A', 1, "C", 2.0),
    brp_atom('A', 2, "N", 3.0),
    brp_atom('A', 2, "CA", 4.0),
    brp_atom('A', 2, "C", 5.0),
    brp_atom('A', 4, "N", 6.0),
    brp_atom('A', 4, "CA", 7.0),
    brp_atom('A', 4, "C", 8.0),
  ])
  let ids = @src.biotite_create_continuous_res_ids(arr)
  assert_eq(ids.length(), 9)
  // Residues 1, 2, 4 become 1, 2, 3 (three atoms each).
  for i = 0; i < 3; i = i + 1 {
    assert_eq(ids[i], 1)
    assert_eq(ids[3 + i], 2)
    assert_eq(ids[6 + i], 3)
  }
  // restart_each_chain=false keeps the global numbering (same for a
  // single chain).
  let ids_global = @src.biotite_create_continuous_res_ids(
    arr,
    restart_each_chain=false,
  )
  for i = 0; i < 9; i = i + 1 {
    assert_eq(ids_global[i], ids[i])
  }
}

///|
test "repair create continuous res ids restarts per chain" {
  // Chain A has residues 1, 2; chain B has residues 5, 7 (gaps). With
  // restart_each_chain=true chain B restarts at 1 and becomes 1, 2.
  let arr = @src.AtomArray::from_atoms([
    brp_atom('A', 1, "N", 0.0),
    brp_atom('A', 2, "N", 1.0),
    brp_atom('B', 5, "N", 2.0),
    brp_atom('B', 7, "N", 3.0),
  ])
  let ids = @src.biotite_create_continuous_res_ids(arr)
  assert_eq(ids.length(), 4)
  assert_eq(ids[0], 1)
  assert_eq(ids[1], 2)
  assert_eq(ids[2], 1)
  assert_eq(ids[3], 2)
  // Without restart the IDs run 1..4 across both chains.
  let ids_global = @src.biotite_create_continuous_res_ids(
    arr,
    restart_each_chain=false,
  )
  assert_eq(ids_global[0], 1)
  assert_eq(ids_global[1], 2)
  assert_eq(ids_global[2], 3)
  assert_eq(ids_global[3], 4)
}

///|
test "repair infer elements official doctest" {
  // Official doctest: infer_elements(["CA", "C", "C1", "OD1", "HD21",
  // "1H", "FE"]) == ['C', 'C', 'C', 'O', 'H', 'H', 'FE'].
  let names = ["CA", "C", "C1", "OD1", "HD21", "1H", "FE"]
  let elems = @src.biotite_infer_elements(names)
  assert_eq(elems.length(), 7)
  assert_eq(elems[0], "C")
  assert_eq(elems[1], "C")
  assert_eq(elems[2], "C")
  assert_eq(elems[3], "O")
  assert_eq(elems[4], "H")
  assert_eq(elems[5], "H")
  assert_eq(elems[6], "FE")
  // Two-character elements behind the biomolecule-first-char rule:
  // "SE" starts with S -> "S" (official rule order), while "ZN" has
  // no biomolecule first character and matches the table exactly.
  assert_eq(@src.biotite_infer_elements(["SE"])[0], "S")
  assert_eq(@src.biotite_infer_elements(["ZN"])[0], "ZN")
  // Unresolvable names yield the empty string (official warns first).
  assert_eq(@src.biotite_infer_elements(["XQ"])[0], "")
}

///|
test "repair create atom names official doctest" {
  // Official doctest (uracil): elements
  // ['N','C','O','N','C','O','C','C','H','H','H','H'] map to
  // ['N1','C1','O1','N2','C2','O2','C3','C4','H1','H2','H3','H4'].
  let elements = ["N", "C", "O", "N", "C", "O", "C", "C", "H", "H", "H", "H"]
  let names = @src.biotite_create_atom_names(elements)
  let expected = [
    "N1", "C1", "O1", "N2", "C2", "O2", "C3", "C4", "H1", "H2", "H3", "H4",
  ]
  assert_eq(names.length(), expected.length())
  for i = 0; i < expected.length(); i = i + 1 {
    assert_eq(names[i], expected[i])
  }
  // For a single atom the element itself is the name (official
  // residue("CL") doctest).
  let single = @src.biotite_create_atom_names(["CL"])
  assert_eq(single.length(), 1)
  assert_eq(single[0], "CL")
}
