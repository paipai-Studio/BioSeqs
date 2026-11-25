///|
/// Tests for biotite_formal_charges.mbt — formal and partial charges.
test "biotite_formal_charge_carbon_and_oxygen" {
  assert_eq(@src.biotite_formal_charge("C", 4.0), 0)
  assert_eq(@src.biotite_formal_charge("O", 2.0, n_lone_pairs=2), 0)
  assert_eq(@src.biotite_formal_charge("N", 4.0, n_lone_pairs=0), 1)
}

///|
test "biotite_formal_charges_batch_and_verify" {
  let atoms : Array[@src.BiotiteFormalChargeAtom] = [
    { element: "C", bond_order_sum: 4.0, n_lone_pairs: -1, },
    { element: "O", bond_order_sum: 2.0, n_lone_pairs: 2, },
  ]
  let charges = @src.biotite_formal_charges(atoms)
  assert_eq(charges.length(), 2)
  assert_eq(charges[0], 0)
  assert_eq(charges[1], 0)
  assert_true(@src.biotite_verify_total_charge(atoms, 0))
  assert_true(!@src.biotite_verify_total_charge(atoms, 2))
}

///|
test "biotite_partial_charges_length" {
  let charges = @src.biotite_partial_charges(["C", "O"], [(0, 1, 2.0)])
  assert_eq(charges.length(), 2)
}

///|
test "biotite_formal_charge_unknown_element" {
  assert_eq(@src.biotite_formal_charge("XX", 1.0), 0)
}
