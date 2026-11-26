///|
/// Tests for biotite_structure_charges.mbt (PEOE / Gasteiger-Marsili,
/// verified against the official Biotite/RDKit doctest values for
/// fluoromethane).
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
test "peoe unknown element falls back to X parameters" {
  // Fe has no Gasteiger parameters: falls back to (0, 0, 0), so it
  // neither receives nor donates charge through its bonds.
  let arr = @src.AtomArray::from_atoms([
    bch_gst_atom("Fe", 0.0, 0.0, 0.0),
    bch_gst_atom("O", 1.5, 0.0, 0.0),
    bch_gst_atom("H", 2.4, 0.0, 0.0),
  ])
  let bonds = @src.BondList::new(n_atoms=3).add(0, 1).add(1, 2)
  let q = @src.biotite_partial_charges_peoe(arr, bonds)
  // Fe stays at 0 (X parameters), O pulls charge from H.
  assert_true(@moonbit.assert_num_close(q[0], 0.0, 1.0e-12))
  assert_true(q[1] < 0.0)
  assert_true(q[2] > 0.0)
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
  match @src.biotite_poeo_params("Fe", "sp3") {
    Some(_) => fail("Fe should have no parameters")
    None => ()
  }
  match @src.biotite_poeo_params("H", "*") {
    Some((a, _b, _c)) =>
      assert_true(@moonbit.assert_num_close(a, 7.17, 1.0e-12))
    None => fail("H must have parameters")
  }
}

///|
test "peoe sulfur so2 special case" {
  // S with two oxygen neighbors uses the so2 parameters (A=12.00).
  let arr = @src.AtomArray::from_atoms([
    bch_gst_atom("S", 0.0, 0.0, 0.0),
    bch_gst_atom("O", 1.5, 0.0, 0.0),
    bch_gst_atom("O", -1.5, 0.0, 0.0),
    bch_gst_atom("H", 0.0, 1.2, 0.0),
  ])
  let bonds = @src.BondList::new(n_atoms=4).add(0, 1).add(0, 2).add(0, 3)
  let q = @src.biotite_partial_charges_peoe(arr, bonds)
  // The hydrogens donate density to sulfur (chi_S(0) = 12.0 > chi_H),
  // while both oxygens pull density away; sulfur ends slightly
  // negative, oxygens negative, hydrogen positive.
  assert_true(q[0] < 0.0)
  assert_true(q[1] < 0.0)
  assert_true(q[2] < 0.0)
  assert_true(q[3] > 0.0)
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
}
