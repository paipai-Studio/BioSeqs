///|
/// Chemical utilities module tests (Bio.PDB.chem_utils equivalent).

///|
test "chem_vdw_radius_carbon" {
  let r = @src.chem_vdw_radius("C")
  assert_true(r > 1.5)
  assert_true(r < 2.0)
}

///|
test "chem_vdw_radius_hydrogen" {
  let r = @src.chem_vdw_radius("H")
  assert_true(r > 1.0)
  assert_true(r < 1.5)
}

///|
test "chem_covalent_radius_carbon" {
  let r = @src.chem_covalent_radius("C")
  assert_true(r > 0.7)
  assert_true(r < 0.8)
}

///|
test "chem_covalent_radius_unknown" {
  let r = @src.chem_covalent_radius("X")
  assert_true(r > 0.5) // should get default (carbon-like)
  assert_true(r < 1.0)
}

///|
test "chem_bond_length_same_point" {
  let coords = [0.0, 0.0, 0.0]
  let dist = @src.chem_bond_length(coords, coords)
  assert_true(dist < 0.001)
}

///|
test "chem_bond_length_known" {
  let a = [0.0, 0.0, 0.0]
  let b = [1.0, 0.0, 0.0]
  let dist = @src.chem_bond_length(a, b)
  assert_true(dist > 0.999)
  assert_true(dist < 1.001)
}

///|
test "chem_bond_length_3d" {
  let a = [0.0, 0.0, 0.0]
  let b = [3.0, 4.0, 0.0]
  let dist = @src.chem_bond_length(a, b)
  assert_true(dist > 4.99)
  assert_true(dist < 5.01)
}

///|
test "chem_bond_angle_collinear" {
  let a = [0.0, 0.0, 0.0]
  let b = [1.0, 0.0, 0.0]
  let c = [2.0, 0.0, 0.0]
  let angle = @src.chem_bond_angle(a, b, c)
  assert_true(angle > 179.0)
  assert_true(angle < 181.0)
}

///|
test "chem_bond_angle_right" {
  let a = [1.0, 0.0, 0.0]
  let b = [0.0, 0.0, 0.0]
  let c = [0.0, 1.0, 0.0]
  let angle = @src.chem_bond_angle(a, b, c)
  assert_true(angle > 89.0)
  assert_true(angle < 91.0)
}

///|
test "chem_dihedral_angle_zero" {
  let a = [1.0, 0.0, 0.0]
  let b = [0.0, 0.0, 0.0]
  let c = [0.0, 1.0, 0.0]
  let d = [1.0, 1.0, 0.0]
  let angle = @src.chem_dihedral_angle(a, b, c, d)
  assert_true(angle > -1.0)
  assert_true(angle < 1.0)
}

///|
test "chem_dihedral_angle_180" {
  let a = [1.0, 0.0, 0.0]
  let b = [0.0, 0.0, 0.0]
  let c = [0.0, 1.0, 0.0]
  let d = [-1.0, 1.0, 0.0]
  let angle = @src.chem_dihedral_angle(a, b, c, d)
  assert_true(angle > 179.0)
}

///|
test "chem_empirical_formula_water" {
  let counts : Map[String, Int] = Map([], capacity=8)
  counts["H"] = 2
  counts["O"] = 1
  let formula = @src.chem_empirical_formula(counts)
  assert_true(formula.length() > 0)
}

///|
test "chem_formula_weight_water" {
  let weight = @src.chem_formula_weight("H2O")
  assert_true(weight > 18.0)
  assert_true(weight < 19.0)
}

///|
test "chem_formula_weight_glucose" {
  let weight = @src.chem_formula_weight("C6H12O6")
  assert_true(weight > 179.0)
  assert_true(weight < 181.0)
}

///|
test "chem_isoelectric_point_lysine" {
  let pI = @src.chem_isoelectric_point("K")
  assert_true(pI > 7.0)
}

///|
test "chem_isoelectric_point_unknown" {
  let pI = @src.chem_isoelectric_point("X")
  assert_true(pI > 0.0)
  assert_true(pI < 14.0)
}

///|
test "chem_ionic_radius_default" {
  let r = @src.chem_ionic_radius("Na")
  assert_true(r > 0.9)
  assert_true(r < 1.1)
}

///|
test "chem_hydrogen_bond_length" {
  let d = @src.chem_hydrogen_bond_length("OH", "OH")
  assert_true(d > 2.5)
  assert_true(d < 3.5)
}
