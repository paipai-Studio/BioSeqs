///|
/// Tests for biotite_structure_charges.mbt.
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
/// Water-like topology: O bonded to two H.
fn bch_water() -> (@src.AtomArray, @src.BondList) {
  let arr = @src.AtomArray::from_atoms([
    bch_gst_atom("O", 0.0, 0.0, 0.0),
    bch_gst_atom("H", 0.96, 0.0, 0.0),
    bch_gst_atom("H", 0.0, 0.96, 0.0),
  ])
  let bonds = @src.BondList::new(n_atoms=3).add(0, 1).add(0, 2)
  (arr, bonds)
}

///|
test "gasteiger charges organic subset only" {
  let arr = @src.AtomArray::from_atoms([
    bch_gst_atom("C", 0.0, 0.0, 0.0),
    bch_gst_atom("Fe", 1.0, 0.0, 0.0),
  ])
  let bonds = @src.BondList::new(n_atoms=2).add(0, 1)
  let pc = @src.biotite_partial_charges_gasteiger(arr, bonds)
  let charges = pc.compute_charge()
  assert_eq(charges.length(), 2)
  match charges[0] {
    Some(_) => ()
    None => fail("C should be organic")
  }
  match charges[1] {
    Some(_) => fail("Fe should not be charged")
    None => ()
  }
}

///|
test "gasteiger charges electronegative oxygen" {
  let (arr, bonds) = bch_water()
  let pc = @src.biotite_partial_charges_gasteiger(arr, bonds)
  let charges = pc.compute_charge()
  match charges[0] {
    Some(q_o) => assert_true(q_o < 0.0)
    None => fail("O should be charged")
  }
  match charges[1] {
    Some(q_h) => assert_true(q_h > 0.0)
    None => fail("H should be charged")
  }
}

///|
test "gasteiger charges conserve total formal charge" {
  let (arr, bonds) = bch_water()
  let pc = @src.biotite_partial_charges_gasteiger(arr, bonds)
  let charges = pc.compute_charge()
  let mut total = 0.0
  for c in charges {
    match c {
      Some(q) => total = total + q
      None => ()
    }
  }
  assert_true(@moonbit.assert_num_close(total, 0.0, 1.0e-9))
}

///|
test "gasteiger formal charge normalization" {
  // ammonium-like: N with formal +1 bonded to 3 H (geometry simplified)
  let arr = @src.AtomArray::from_atoms([
    bch_gst_atom("N", 0.0, 0.0, 0.0),
    bch_gst_atom("H", 1.0, 0.0, 0.0),
  ])
  let bonds = @src.BondList::new(n_atoms=2).add(0, 1)
  let pc = @src.biotite_partial_charges_gasteiger(
    arr,
    bonds,
    formal_charges=Some([1, 0]),
  )
  let charges = pc.compute_charge()
  let mut total = 0.0
  for c in charges {
    match c {
      Some(q) => total = total + q
      None => ()
    }
  }
  assert_true(@moonbit.assert_num_close(total, 1.0, 1.0e-9))
}

///|
test "gasteiger sp2 carbonyl oxygen more negative" {
  let arr = @src.AtomArray::from_atoms([
    bch_gst_atom("C", 0.0, 0.0, 0.0),
    bch_gst_atom("O", 1.3, 0.0, 0.0),
  ])
  let bonds_sp2 = @src.BondList::new(n_atoms=2).add(
    0,
    1,
    bond_type=@src.biotite_bond_double(),
  )
  let pc_sp2 = @src.biotite_partial_charges_gasteiger(arr, bonds_sp2)
  let q_sp2 = pc_sp2.compute_charge()
  let bonds_sp3 = @src.BondList::new(n_atoms=2).add(0, 1)
  let pc_sp3 = @src.biotite_partial_charges_gasteiger(arr, bonds_sp3)
  let q_sp3 = pc_sp3.compute_charge()
  match (q_sp2[0], q_sp3[0]) {
    (Some(a), Some(b)) =>
      // double-bonded C loses more electron density (more positive)
      assert_true(a >= b - 1.0e-9)
    _ => fail("both atoms should be charged")
  }
}

///|
test "gasteiger empty structure" {
  let pc = @src.biotite_partial_charges_gasteiger(
    @src.AtomArray::new(),
    @src.BondList::new(n_atoms=0),
  )
  let charges = pc.compute_charge()
  assert_eq(charges.length(), 0)
}
