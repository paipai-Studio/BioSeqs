///|
/// Tests for biotite_structure_charges.mbt (PEOE / Gasteiger-Marsili,
/// verified against the official Biotite Rust implementation: the
/// fluoromethane doctest, the 17 molecules of table 3 of the
/// Gasteiger-Marsili publication, and the NaN semantics for
/// unparametrized atoms).
fn bch_gst_atom(
  element : String,
  x : Double,
  y : Double,
  z : Double,
) -> @src.Atom {
  @src.Atom::new(
    name="A",
    coord=@src.Vector3::new(x, y, z),
    resname="LIG",
    chainid='A',
    resseq=1,
    element~,
  )
}

///|
/// Fluoromethane topology (official doctest molecule CF0):
/// C bonded to F and three H.
fn bch_fluoromethane() -> (@src.AtomArray, @src.BondList) {
  let arr = @src.AtomArray::from_atoms([
    bch_gst_atom("C", 0.0, 0.0, 0.0),
    bch_gst_atom("F", 1.3, 0.0, 0.0),
    bch_gst_atom("H", -0.5, 0.9, 0.0),
    bch_gst_atom("H", -0.5, -0.45, 0.78),
    bch_gst_atom("H", -0.5, -0.45, -0.78),
  ])
  let bonds = @src.BondList::new(n_atoms=5)
    .add(0, 1)
    .add(0, 2)
    .add(0, 3)
    .add(0, 4)
  (arr, bonds)
}

///|
/// Bond type codes: 0 = single, 1 = double, 2 = triple, 3 = any.
fn bch_bond(code : Int32) -> @src.BondType {
  match code {
    1 => @src.biotite_bond_double()
    2 => @src.biotite_bond_triple()
    3 => @src.biotite_bond_any()
    _ => @src.biotite_bond_single()
  }
}

///|
/// Build a molecule from an element list and (atom1, atom2, code)
/// bond specs; coordinates are arbitrary (only the BondList matters).
fn bch_mol(
  elements : Array[String],
  bond_specs : Array[(Int32, Int32, Int32)],
) -> (@src.AtomArray, @src.BondList) {
  let atoms : Array[@src.Atom] = []
  for i = 0; i < elements.length(); i = i + 1 {
    atoms.push(bch_gst_atom(elements[i], i.to_double(), 0.0, 0.0))
  }
  let arr = @src.AtomArray::from_atoms(atoms)
  let mut bonds = @src.BondList::new(n_atoms=elements.length())
  for (a, b, code) in bond_specs {
    bonds = bonds.add(a, b, bond_type=bch_bond(code))
  }
  (arr, bonds)
}

///|
test "peoe official doctest single iteration" {
  // Official docstring: partial_charges(fluoromethane, iteration_step_num=1)
  // == [0.115, -0.175, 0.020, 0.020, 0.020]
  let (arr, bonds) = bch_fluoromethane()
  let q = @src.biotite_partial_charges_peoe(arr, bonds, iteration_step_num=1)
  assert_eq(q.length(), 5)
  assert_true(@moonbit.assert_num_close(q[0], 0.115, 2.0e-3))
  assert_true(@moonbit.assert_num_close(q[1], -0.175, 2.0e-3))
  assert_true(@moonbit.assert_num_close(q[2], 0.020, 2.0e-3))
  assert_true(@moonbit.assert_num_close(q[3], 0.020, 2.0e-3))
  assert_true(@moonbit.assert_num_close(q[4], 0.020, 2.0e-3))
}

///|
test "peoe official doctest six iterations" {
  // Official docstring: partial_charges(fluoromethane, iteration_step_num=6)
  // == [0.079, -0.253, 0.058, 0.058, 0.058]
  let (arr, bonds) = bch_fluoromethane()
  let q = @src.biotite_partial_charges_peoe(arr, bonds, iteration_step_num=6)
  assert_true(@moonbit.assert_num_close(q[0], 0.079, 5.0e-3))
  assert_true(@moonbit.assert_num_close(q[1], -0.253, 5.0e-3))
  assert_true(@moonbit.assert_num_close(q[2], 0.058, 5.0e-3))
  assert_true(@moonbit.assert_num_close(q[3], 0.058, 5.0e-3))
  assert_true(@moonbit.assert_num_close(q[4], 0.058, 5.0e-3))
}

///|
test "peoe default iteration count is six" {
  let (arr, bonds) = bch_fluoromethane()
  let q_default = @src.biotite_partial_charges_peoe(arr, bonds)
  let q_six = @src.biotite_partial_charges_peoe(
    arr,
    bonds,
    iteration_step_num=6,
  )
  assert_eq(q_default.length(), q_six.length())
  for i = 0; i < q_default.length(); i = i + 1 {
    assert_true(@moonbit.assert_num_close(q_default[i], q_six[i], 1.0e-12))
  }
}

///|
test "peoe conserves total charge" {
  let (arr, bonds) = bch_fluoromethane()
  let q = @src.biotite_partial_charges_peoe(arr, bonds)
  let mut total = 0.0
  for v in q {
    total = total + v
  }
  assert_true(@moonbit.assert_num_close(total, 0.0, 1.0e-9))
}

///|
test "peoe formal charges seed and conserve" {
  // Ammonium-like: N with formal +1 bonded to 3 H.
  let arr = @src.AtomArray::from_atoms([
    bch_gst_atom("N", 0.0, 0.0, 0.0),
    bch_gst_atom("H", 1.0, 0.0, 0.0),
    bch_gst_atom("H", 0.0, 1.0, 0.0),
    bch_gst_atom("H", 0.0, 0.0, 1.0),
  ])
  let bonds = @src.BondList::new(n_atoms=4).add(0, 1).add(0, 2).add(0, 3)
  let q = @src.biotite_partial_charges_peoe(
    arr,
    bonds,
    formal_charges=Some([1, 0, 0, 0]),
  )
  let mut total = 0.0
  for v in q {
    total = total + v
  }
  assert_true(@moonbit.assert_num_close(total, 1.0, 1.0e-9))
  // PEOE semantics: the formal positive charge raises the nitrogen
  // electronegativity (chi(q=1) = A + B + C), so it pulls electron
  // density back from its hydrogens and ends up negative while the
  // hydrogens turn positive.
  assert_true(q[0] < 0.0)
  assert_true(q[1] > 0.0)
}

///|
test "peoe electronegative oxygen" {
  // Water topology: O bonded to two H.
  let arr = @src.AtomArray::from_atoms([
    bch_gst_atom("O", 0.0, 0.0, 0.0),
    bch_gst_atom("H", 0.96, 0.0, 0.0),
    bch_gst_atom("H", 0.0, 0.96, 0.0),
  ])
  let bonds = @src.BondList::new(n_atoms=3).add(0, 1).add(0, 2)
  let q = @src.biotite_partial_charges_peoe(arr, bonds)
  assert_true(q[0] < 0.0)
  assert_true(q[1] > 0.0)
  assert_true(q[2] > 0.0)
}

///|
test "peoe without bonds returns formal charges" {
  let arr = @src.AtomArray::from_atoms([
    bch_gst_atom("Na", 0.0, 0.0, 0.0),
    bch_gst_atom("Cl", 2.0, 0.0, 0.0),
  ])
  let bonds = @src.BondList::new(n_atoms=2)
  let q = @src.biotite_partial_charges_peoe(
    arr,
    bonds,
    formal_charges=Some([1, -1]),
  )
  assert_true(@moonbit.assert_num_close(q[0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(q[1], -1.0, 1.0e-12))
}

///|
test "peoe unparametrized element yields NaN charges" {
  // Fe has no Gasteiger parameters: like the official Rust
  // implementation its charge becomes NaN (and stays NaN), while the
  // parametrized O-H part still polarizes normally.
  let arr = @src.AtomArray::from_atoms([
    bch_gst_atom("Fe", 0.0, 0.0, 0.0),
    bch_gst_atom("O", 1.5, 0.0, 0.0),
    bch_gst_atom("H", 2.4, 0.0, 0.0),
  ])
  let bonds = @src.BondList::new(n_atoms=3).add(0, 1).add(1, 2)
  let q = @src.biotite_partial_charges_peoe(arr, bonds)
  assert_true(q[0].is_nan())
  assert_true(q[1] < 0.0)
  assert_true(q[2] > 0.0)
}

///|
test "peoe unparametrized sulfur double bond yields NaN" {
  // Thioformaldehyde (official test_valence_state_not_parametrized):
  // the C=S double bond implies sp2 sulfur, which has no parameters.
  // The sulfur charge is NaN while the carbon stays parametrized and
  // pulls density from its hydrogens (C more negative than H).
  let (arr, bonds) = bch_mol(["C", "S", "H", "H"], [
    (0, 1, 1),
    (0, 2, 0),
    (0, 3, 0),
  ])
  let q = @src.biotite_partial_charges_peoe(arr, bonds)
  assert_true(q[1].is_nan())
  assert_true(q[0] < q[2])
  assert_true(q[0] < q[3])
}

///|
test "peoe sulfur sp3 single bonds" {
  // S with only single bonds uses the sp3 parameters: oxygen pulls
  // density away (chi_O(0) = 14.18 > chi_S(0) = 10.14) while the
  // hydrogen donates, so sulfur ends positive, oxygens negative.
  let (arr, bonds) = bch_mol(["S", "O", "O", "H"], [
    (0, 1, 0),
    (0, 2, 0),
    (0, 3, 0),
  ])
  let q = @src.biotite_partial_charges_peoe(arr, bonds)
  assert_true(q[0] > 0.0)
  assert_true(q[1] < 0.0)
  assert_true(q[2] < 0.0)
  assert_true(q[3] > 0.0)
}

///|
test "peoe sp2 carbonyl oxygen more negative than sp3 ether" {
  let arr = @src.AtomArray::from_atoms([
    bch_gst_atom("C", 0.0, 0.0, 0.0),
    bch_gst_atom("O", 1.3, 0.0, 0.0),
  ])
  let bonds_double = @src.BondList::new(n_atoms=2).add(
    0,
    1,
    bond_type=@src.biotite_bond_double(),
  )
  let q_double = @src.biotite_partial_charges_peoe(arr, bonds_double)
  let bonds_single = @src.BondList::new(n_atoms=2).add(0, 1)
  let q_single = @src.biotite_partial_charges_peoe(arr, bonds_single)
  // Both polarize C positive / O negative; the double bond transfers
  // more charge per the PEOE parameterization.
  assert_true(q_double[1] < 0.0)
  assert_true(q_single[1] < 0.0)
  assert_true(q_double[1] < q_single[1])
}

///|
test "peoe parameter lookup table" {
  match @src.biotite_poeo_params("C", "sp3") {
    Some((a, b, c)) => {
      assert_true(@moonbit.assert_num_close(a, 7.98, 1.0e-12))
      assert_true(@moonbit.assert_num_close(b, 9.18, 1.0e-12))
      assert_true(@moonbit.assert_num_close(c, 1.88, 1.0e-12))
    }
    None => fail("C sp3 must have parameters")
  }
  // Official Rust: carbon sp2 shares B/C with sp3.
  match @src.biotite_poeo_params("C", "sp2") {
    Some((a, b, c)) => {
      assert_true(@moonbit.assert_num_close(a, 8.79, 1.0e-12))
      assert_true(@moonbit.assert_num_close(b, 9.18, 1.0e-12))
      assert_true(@moonbit.assert_num_close(c, 1.88, 1.0e-12))
    }
    None => fail("C sp2 must have parameters")
  }
  match @src.biotite_poeo_params("Fe", "sp3") {
    Some(_) => fail("Fe should have no parameters")
    None => ()
  }
  // P/Si/B/Be/Mg/Al are not in the official parameter table.
  match @src.biotite_poeo_params("P", "sp3") {
    Some(_) => fail("P should have no parameters")
    None => ()
  }
  match @src.biotite_poeo_params("H", "*") {
    Some((a, _b, _c)) =>
      assert_true(@moonbit.assert_num_close(a, 7.17, 1.0e-12))
    None => fail("H must have parameters")
  }
}

///|
test "peoe any bond type falls back to partner count" {
  // With every bond typed `Any`, the hybridization is inferred from
  // the binding-partner count (C: 4 -> sp3, H: 1 -> sp3), which for
  // methane matches the typed topology exactly.
  let (arr, bonds) = bch_mol(["C", "H", "H", "H", "H"], [
    (0, 1, 3),
    (0, 2, 3),
    (0, 3, 3),
    (0, 4, 3),
  ])
  let q_any = @src.biotite_partial_charges_peoe(arr, bonds)
  let (arr2, bonds2) = bch_mol(["C", "H", "H", "H", "H"], [
    (0, 1, 0),
    (0, 2, 0),
    (0, 3, 0),
    (0, 4, 0),
  ])
  let q_typed = @src.biotite_partial_charges_peoe(arr2, bonds2)
  for i = 0; i < 5; i = i + 1 {
    assert_true(@moonbit.assert_num_close(q_any[i], q_typed[i], 1.0e-12))
  }
}

///|
test "peoe gasteiger-marsili table 3 molecules" {
  // The 17 molecules of table 3 of the Gasteiger-Marsili publication,
  // comparing the computed carbon charges against the official test
  // expectations (absolute tolerance 1e-2, like upstream).
  let cases : Array[
    (Array[String], Array[(Int32, Int32, Int32)], Array[Int32], Array[Double]),
  ] = [
    (
      // methane
      ["C", "H", "H", "H", "H"],
      [(0, 1, 0), (0, 2, 0), (0, 3, 0), (0, 4, 0)],
      [0],
      [-0.078],
    ),
    (
      // ethane
      ["C", "C", "H", "H", "H", "H", "H", "H"],
      [
        (0, 1, 0),
        (0, 2, 0),
        (0, 3, 0),
        (0, 4, 0),
        (1, 5, 0),
        (1, 6, 0),
        (1, 7, 0),
      ],
      [0, 1],
      [-0.068, -0.068],
    ),
    (
      // ethylene
      ["C", "C", "H", "H", "H", "H"],
      [(0, 1, 1), (0, 2, 0), (0, 3, 0), (1, 4, 0), (1, 5, 0)],
      [0, 1],
      [-0.106, -0.106],
    ),
    (
      // acetylene
      ["C", "C", "H", "H"],
      [(0, 1, 2), (0, 2, 0), (1, 3, 0)],
      [0, 1],
      [-0.122, -0.122],
    ),
    (
      // difluoromethane
      ["C", "F", "F", "H", "H"],
      [(0, 1, 0), (0, 2, 0), (0, 3, 0), (0, 4, 0)],
      [0],
      [0.23],
    ),
    (
      // trifluoromethane
      ["C", "F", "F", "F", "H"],
      [(0, 1, 0), (0, 2, 0), (0, 3, 0), (0, 4, 0)],
      [0],
      [0.38],
    ),
    (
      // tetrafluoromethane
      ["C", "F", "F", "F", "F"],
      [(0, 1, 0), (0, 2, 0), (0, 3, 0), (0, 4, 0)],
      [0],
      [0.561],
    ),
    (
      // fluoroethane
      ["C", "C", "F", "H", "H", "H", "H", "H"],
      [
        (0, 1, 0),
        (0, 2, 0),
        (0, 3, 0),
        (0, 4, 0),
        (1, 5, 0),
        (1, 6, 0),
        (1, 7, 0),
      ],
      [0, 1],
      [0.087, -0.037],
    ),
    (
      // trifluoroethane
      ["C", "C", "F", "F", "F", "H", "H", "H"],
      [
        (0, 1, 0),
        (0, 2, 0),
        (0, 3, 0),
        (0, 4, 0),
        (1, 5, 0),
        (1, 6, 0),
        (1, 7, 0),
      ],
      [0, 1],
      [0.387, 0.039],
    ),
    (
      // methanol
      ["C", "O", "H", "H", "H", "H"],
      [(0, 1, 0), (0, 2, 0), (0, 3, 0), (0, 4, 0), (1, 5, 0)],
      [0],
      [0.033],
    ),
    (
      // dimethyl ether
      ["C", "C", "O", "H", "H", "H", "H", "H", "H"],
      [
        (0, 2, 0),
        (1, 2, 0),
        (0, 3, 0),
        (0, 4, 0),
        (0, 5, 0),
        (1, 6, 0),
        (1, 7, 0),
        (1, 8, 0),
      ],
      [0, 1],
      [0.036, 0.036],
    ),
    (
      // formaldehyde
      ["C", "O", "H", "H"],
      [(0, 1, 1), (0, 2, 0), (0, 3, 0)],
      [0],
      [0.115],
    ),
    (
      // acetaldehyde
      ["C", "C", "O", "H", "H", "H", "H"],
      [(0, 1, 0), (1, 2, 1), (0, 3, 0), (0, 4, 0), (0, 5, 0), (1, 6, 0)],
      [0, 1],
      [-0.009, 0.123],
    ),
    (
      // acetone
      ["C", "C", "C", "O", "H", "H", "H", "H", "H", "H"],
      [
        (0, 1, 0),
        (1, 2, 0),
        (1, 3, 1),
        (0, 4, 0),
        (0, 5, 0),
        (0, 6, 0),
        (2, 7, 0),
        (2, 8, 0),
        (2, 9, 0),
      ],
      [0, 1, 2],
      [-0.006, 0.131, -0.006],
    ),
    (
      // hydrogen cyanide
      ["C", "N", "H"],
      [(0, 1, 2), (0, 2, 0)],
      [0],
      [0.051],
    ),
    (
      // acetonitrile
      ["C", "C", "N", "H", "H", "H"],
      [(0, 1, 0), (1, 2, 2), (0, 3, 0), (0, 4, 0), (0, 5, 0)],
      [0, 1],
      [0.023, 0.06],
    ),
  ]
  for case in cases {
    let (elements, bond_specs, carbon_indices, expected) = case
    let (arr, bonds) = bch_mol(elements, bond_specs)
    let q = @src.biotite_partial_charges_peoe(arr, bonds)
    for k = 0; k < expected.length(); k = k + 1 {
      let idx = carbon_indices[k]
      assert_true(@moonbit.assert_num_close(q[idx], expected[k], 1.0e-2))
    }
    // Total charge is conserved (official test_total_charge_zero).
    let mut total = 0.0
    for v in q {
      total = total + v
    }
    assert_true(@moonbit.assert_num_close(total, 0.0, 1.0e-6))
  }
}

///|
test "bond type constructors" {
  let bt = @src.biotite_bond_double()
  match bt {
    _ => ()
  }
  let bt3 = @src.biotite_bond_triple()
  match bt3 {
    _ => ()
  }
  let bt1 = @src.biotite_bond_single()
  match bt1 {
    _ => ()
  }
  let bta = @src.biotite_bond_any()
  match bta {
    _ => ()
  }
}
