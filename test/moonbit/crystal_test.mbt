///|
/// Tests for the Bio.Crystal module (crystal.mbt).
///
/// Covers unit-cell geometry, space groups, atom sites, CIF parsing, and
/// crystal-structure analysis.

// ============================================================================
// Unit cell
// ============================================================================

///|
test "unit_cell_construction" {
  let cell = @src.UnitCell::new(
    a=5.0,
    b=6.0,
    c=7.0,
    alpha=90.0,
    beta=90.0,
    gamma=90.0,
  )
  assert_true((cell.a - 5.0).abs() < 0.001)
  assert_true((cell.b - 6.0).abs() < 0.001)
  assert_true((cell.c - 7.0).abs() < 0.001)
  assert_true((cell.alpha - 90.0).abs() < 0.001)
}

///|
test "unit_cell_volume_orthorhombic" {
  // Orthorhombic: alpha=beta=gamma=90 -> V = a*b*c
  let cell = @src.UnitCell::new(
    a=5.0,
    b=6.0,
    c=7.0,
    alpha=90.0,
    beta=90.0,
    gamma=90.0,
  )
  let v = @src.unit_cell_volume(cell)
  assert_true((v - 210.0).abs() < 0.001)
}

///|
test "unit_cell_volume_hexagonal" {
  // Hexagonal: alpha=beta=90, gamma=120 -> V = a*b*c*sqrt(3)/2
  let cell = @src.UnitCell::new(
    a=5.0,
    b=5.0,
    c=7.0,
    alpha=90.0,
    beta=90.0,
    gamma=120.0,
  )
  let v = @src.unit_cell_volume(cell)
  let expected = 5.0 * 5.0 * 7.0 * (3.0).sqrt() / 2.0
  assert_true((v - expected).abs() < 0.01)
}

///|
test "unit_cell_volume_triclinic" {
  // Triclinic with non-trivial angles.
  let cell = @src.UnitCell::new(
    a=5.0,
    b=6.0,
    c=7.0,
    alpha=80.0,
    beta=85.0,
    gamma=75.0,
  )
  let v = @src.unit_cell_volume(cell)
  assert_true(v > 0.0)
  // V should be less than the orthogonal volume (210).
  assert_true(v < 210.0)
}

///|
test "crystal_density_basic" {
  // Orthorhombic cell, Z=2, MW=78.11 (benzene).
  let cell = @src.UnitCell::new(
    a=5.0,
    b=6.0,
    c=7.0,
    alpha=90.0,
    beta=90.0,
    gamma=90.0,
  )
  let rho = @src.crystal_density(cell, 2, 78.11)
  assert_true(rho > 0.0)
  // Typical organic crystal densities are 1.0-1.5 g/cm^3.
  assert_true(rho > 0.5)
  assert_true(rho < 2.0)
}

///|
test "crystal_density_zero_volume" {
  // Degenerate cell with zero volume.
  let cell = @src.UnitCell::new(
    a=0.0,
    b=0.0,
    c=0.0,
    alpha=90.0,
    beta=90.0,
    gamma=90.0,
  )
  let rho = @src.crystal_density(cell, 2, 78.11)
  assert_true((rho - 0.0).abs() < 0.001)
}

// ============================================================================
// Space groups
// ============================================================================

///|
test "space_group_construction" {
  let sg = @src.SpaceGroup::new(number=19, symbol="P 21 21 21")
  assert_eq(sg.number, 19)
  assert_eq(sg.symbol, "P 21 21 21")
}

///|
test "lookup_space_group_common" {
  let sg1 = @src.lookup_space_group(1)
  assert_eq(sg1.number, 1)
  assert_eq(sg1.symbol, "P 1")
  let sg19 = @src.lookup_space_group(19)
  assert_eq(sg19.number, 19)
  assert_eq(sg19.symbol, "P 21 21 21")
  let sg225 = @src.lookup_space_group(225)
  assert_eq(sg225.number, 225)
  assert_eq(sg225.symbol, "F m -3 m")
}

///|
test "lookup_space_group_unknown" {
  let sg = @src.lookup_space_group(999)
  assert_eq(sg.number, 999)
  assert_eq(sg.symbol, "")
}

// ============================================================================
// Crystal atoms and structures
// ============================================================================

///|
test "crystal_atom_construction" {
  let a = @src.CrystalAtom::new(
    label="C1",
    element="C",
    x=0.5,
    y=0.5,
    z=0.5,
  )
  assert_eq(a.label, "C1")
  assert_eq(a.element, "C")
  assert_true((a.x - 0.5).abs() < 0.001)
  assert_true((a.occupancy - 1.0).abs() < 0.001)
  assert_eq(a.is_fractional, true)
}

///|
test "crystal_atom_cartesian" {
  let a = @src.CrystalAtom::new(
    label="C1",
    element="C",
    x=1.5,
    y=2.5,
    z=3.5,
    occupancy=0.8,
    is_fractional=false,
  )
  assert_eq(a.is_fractional, false)
  assert_true((a.occupancy - 0.8).abs() < 0.001)
}

///|
test "crystal_bond_construction" {
  let b = @src.CrystalBond::new(
    atom1="C1",
    atom2="C2",
    bond_type="DOUBLE",
    distance=1.34,
  )
  assert_eq(b.atom1, "C1")
  assert_eq(b.atom2, "C2")
  assert_eq(b.bond_type, "DOUBLE")
  assert_true((b.distance - 1.34).abs() < 0.001)
}

///|
test "crystal_bond_default_type" {
  let b = @src.CrystalBond::new(atom1="A", atom2="B")
  assert_eq(b.bond_type, "SINGLE")
  assert_true((b.distance - 0.0).abs() < 0.001)
}

///|
test "crystal_structure_construction" {
  let s = @src.crystal_sample_structure()
  assert_eq(s.name, "benzene")
  assert_eq(s.atoms.length(), 6)
  assert_eq(s.bonds.length(), 6)
  assert_eq(s.z, 2)
  assert_true((s.molecular_weight - 78.11).abs() < 0.001)
}

// ============================================================================
// Coordinate transformations
// ============================================================================

///|
test "fractional_to_cartesian_orthorhombic" {
  let cell = @src.UnitCell::new(
    a=10.0,
    b=20.0,
    c=30.0,
    alpha=90.0,
    beta=90.0,
    gamma=90.0,
  )
  let(x, y, z) = @src.fractional_to_cartesian(cell, 0.5, 0.25, 0.1)
  assert_true((x - 5.0).abs() < 0.001)
  assert_true((y - 5.0).abs() < 0.001)
  assert_true((z - 3.0).abs() < 0.001)
}

///|
test "cartesian_to_fractional_orthorhombic" {
  let cell = @src.UnitCell::new(
    a=10.0,
    b=20.0,
    c=30.0,
    alpha=90.0,
    beta=90.0,
    gamma=90.0,
  )
  let(xf, yf, zf) = @src.cartesian_to_fractional(cell, 5.0, 5.0, 3.0)
  assert_true((xf - 0.5).abs() < 0.001)
  assert_true((yf - 0.25).abs() < 0.001)
  assert_true((zf - 0.1).abs() < 0.001)
}

///|
test "coord_conversion_roundtrip_orthorhombic" {
  let cell = @src.UnitCell::new(
    a=10.0,
    b=20.0,
    c=30.0,
    alpha=90.0,
    beta=90.0,
    gamma=90.0,
  )
  let xf0 = 0.3
  let yf0 = 0.6
  let zf0 = 0.9
  let(x, y, z) = @src.fractional_to_cartesian(cell, xf0, yf0, zf0)
  let(xf, yf, zf) = @src.cartesian_to_fractional(cell, x, y, z)
  assert_true((xf - xf0).abs() < 0.01)
  assert_true((yf - yf0).abs() < 0.01)
  assert_true((zf - zf0).abs() < 0.01)
}

///|
test "atoms_to_fractional_conversion" {
  let s = @src.crystal_sample_structure()
  let frac_atoms = @src.atoms_to_fractional(s)
  // Sample atoms are already fractional; should be unchanged.
  assert_eq(frac_atoms.length(), 6)
  assert_true((frac_atoms[0].x - s.atoms[0].x).abs() < 0.001)
}

///|
test "atoms_to_cartesian_conversion" {
  let s = @src.crystal_sample_structure()
  let cart_atoms = @src.atoms_to_cartesian(s)
  assert_eq(cart_atoms.length(), 6)
  for i in 0..<cart_atoms.length() {
    assert_eq(cart_atoms[i].is_fractional, false)
  }
}

// ============================================================================
// Distance and center of mass
// ============================================================================

///|
test "crystal_atom_distance_basic" {
  let s = @src.crystal_sample_structure()
  let d = @src.crystal_atom_distance(s, "C1", "C2")
  assert_true(d.is_some())
  let dist = d.unwrap()
  assert_true(dist > 0.0)
}

///|
test "crystal_atom_distance_missing_atom" {
  let s = @src.crystal_sample_structure()
  let d = @src.crystal_atom_distance(s, "C1", "NONEXISTENT")
  assert_eq(d.is_empty(), true)
}

///|
test "center_of_mass_basic" {
  let s = @src.crystal_sample_structure()
  let(cx, cy, cz) = @src.center_of_mass(s)
  // Center of mass should be a 3-tuple of finite doubles.
  assert_true(cx >= 0.0 || cx < 0.0) // just check it's a number
  assert_true(cy >= 0.0 || cy < 0.0)
  assert_true(cz >= 0.0 || cz < 0.0)
}

///|
test "center_of_mass_empty_structure" {
  let cell = @src.UnitCell::new(
    a=5.0,
    b=5.0,
    c=5.0,
    alpha=90.0,
    beta=90.0,
    gamma=90.0,
  )
  let s = @src.CrystalStructure::new(
    name="empty",
    cell=cell,
    space_group=@src.SpaceGroup::new(number=1, symbol="P 1"),
    atoms=[],
    bonds=[],
    z=1,
    molecular_weight=0.0,
  )
  let(cx, cy, cz) = @src.center_of_mass(s)
  assert_true((cx - 0.0).abs() < 0.001)
  assert_true((cy - 0.0).abs() < 0.001)
  assert_true((cz - 0.0).abs() < 0.001)
}

// ============================================================================
// CIF parsing
// ============================================================================

///|
test "cif_block_construction" {
  let block = @src.CifBlock::new(name="test")
  assert_eq(block.name, "test")
}

///|
test "cif_block_scalar_set_get" {
  let block = @src.CifBlock::new(name="test")
  block.set_scalar("_cell_length_a", "5.418")
  let v = block.get_scalar("_cell_length_a")
  assert_eq(v, "5.418")
}

///|
test "cif_block_scalar_missing_returns_empty" {
  let block = @src.CifBlock::new(name="test")
  let v = block.get_scalar("_nonexistent")
  assert_eq(v, "")
}

///|
test "cif_block_loop_set_get" {
  let block = @src.CifBlock::new(name="test")
  let headers = ["_atom_site_label", "_atom_site_x"]
  let rows = [["C1", "0.5"], ["C2", "0.6"]]
  block.set_loop("_atom_site_label", headers, rows)
  let(h, r) = block.get_loop("_atom_site_label")
  assert_eq(h.length(), 2)
  assert_eq(r.length(), 2)
  assert_eq(r[0][0], "C1")
}

///|
test "cif_block_loop_missing_returns_empty" {
  let block = @src.CifBlock::new(name="test")
  let(h, r) = block.get_loop("_nonexistent")
  assert_eq(h.length(), 0)
  assert_eq(r.length(), 0)
}

///|
test "parse_cif_basic" {
  let text =
    #|data_test
    #|_cell_length_a 5.418
    #|_cell_length_b 5.418
    #|_cell_length_c 7.543
  let blocks = @src.parse_cif(text)
  assert_eq(blocks.length(), 1)
  assert_eq(blocks[0].name, "test")
  assert_eq(blocks[0].get_scalar("_cell_length_a"), "5.418")
  assert_eq(blocks[0].get_scalar("_cell_length_b"), "5.418")
}

///|
test "parse_cif_with_loop" {
  let text =
    #|data_test
    #|loop_
    #|_atom_site_label
    #|_atom_site_x
    #|C1 0.5
    #|C2 0.6
    #|C3 0.7
  let blocks = @src.parse_cif(text)
  assert_eq(blocks.length(), 1)
  let(h, r) = blocks[0].get_loop("_atom_site_label")
  assert_eq(h.length(), 2)
  assert_eq(r.length(), 3)
  assert_eq(r[0][0], "C1")
  assert_eq(r[2][0], "C3")
}

///|
test "parse_cif_multiple_blocks" {
  let text =
    #|data_block1
    #|_cell_length_a 5.0
    #|data_block2
    #|_cell_length_a 6.0
  let blocks = @src.parse_cif(text)
  assert_eq(blocks.length(), 2)
  assert_eq(blocks[0].name, "block1")
  assert_eq(blocks[0].get_scalar("_cell_length_a"), "5.0")
  assert_eq(blocks[1].name, "block2")
  assert_eq(blocks[1].get_scalar("_cell_length_a"), "6.0")
}

///|
test "parse_cif_quoted_values" {
  let text =
    #|data_test
    #|_space_group_name_H-M 'P 21/c'
    #|_chemical_name 'Benzene, C6H6'
  let blocks = @src.parse_cif(text)
  assert_eq(blocks[0].get_scalar("_space_group_name_H-M"), "P 21/c")
  assert_eq(blocks[0].get_scalar("_chemical_name"), "Benzene, C6H6")
}

///|
test "parse_cif_comments" {
  let text =
    #|# top-level comment
    #|data_test
    #|# comment in block
    #|_cell_length_a 5.0
  let blocks = @src.parse_cif(text)
  assert_eq(blocks.length(), 1)
  assert_eq(blocks[0].get_scalar("_cell_length_a"), "5.0")
}

///|
test "cif_block_to_crystal_basic" {
  let text = @src.crystal_sample_cif()
  let blocks = @src.parse_cif(text)
  assert_eq(blocks.length(), 1)
  let s = @src.cif_block_to_crystal(blocks[0])
  assert_eq(s.name, "benzene")
  assert_true((s.cell.a - 5.418).abs() < 0.001)
  assert_true((s.cell.b - 5.418).abs() < 0.001)
  assert_true((s.cell.c - 7.543).abs() < 0.001)
  assert_eq(s.space_group.number, 14)
  assert_eq(s.z, 2)
  assert_true((s.molecular_weight - 78.11).abs() < 0.001)
  // 12 atoms: 6 C + 6 H
  assert_eq(s.atoms.length(), 12)
}

///|
test "cif_block_to_crystal_atom_parsing" {
  let text = @src.crystal_sample_cif()
  let blocks = @src.parse_cif(text)
  let s = @src.cif_block_to_crystal(blocks[0])
  assert_eq(s.atoms[0].label, "C1")
  assert_eq(s.atoms[0].element, "C")
  assert_true((s.atoms[0].x - 0.1234).abs() < 0.001)
  assert_eq(s.atoms[0].is_fractional, true)
  // Check the 7th atom (H1)
  assert_eq(s.atoms[6].label, "H1")
  assert_eq(s.atoms[6].element, "H")
}

// ============================================================================
// Format and summary
// ============================================================================

///|
test "format_crystal_structure_output" {
  let s = @src.crystal_sample_structure()
  let out = @src.format_crystal_structure(s)
  assert_true(out.contains("benzene"))
  assert_true(out.contains("Cell:"))
  assert_true(out.contains("Volume:"))
  assert_true(out.contains("Space group:"))
  assert_true(out.contains("P 21/c"))
  assert_true(out.contains("Atoms (6):"))
  assert_true(out.contains("C1"))
}

///|
test "crystal_summary_basic" {
  let s = @src.crystal_sample_structure()
  let summary = @src.crystal_summary(s)
  assert_eq(summary.name, "benzene")
  assert_eq(summary.n_atoms, 6)
  assert_eq(summary.n_bonds, 6)
  assert_eq(summary.z, 2)
  assert_eq(summary.space_group, "P 21/c")
  assert_true(summary.volume > 0.0)
  assert_true(summary.density > 0.0)
}

///|
test "crystal_summary_density_zero_when_mw_zero" {
  let cell = @src.UnitCell::new(
    a=5.0,
    b=5.0,
    c=5.0,
    alpha=90.0,
    beta=90.0,
    gamma=90.0,
  )
  let s = @src.CrystalStructure::new(
    name="test",
    cell=cell,
    space_group=@src.SpaceGroup::new(number=1, symbol="P 1"),
    atoms=[],
    bonds=[],
    z=1,
    molecular_weight=0.0,
  )
  let summary = @src.crystal_summary(s)
  assert_true((summary.density - 0.0).abs() < 0.001)
}

// ============================================================================
// Sample data validation
// ============================================================================

///|
test "crystal_sample_cif_content" {
  let s = @src.crystal_sample_cif()
  assert_true(s.contains("data_benzene"))
  assert_true(s.contains("_cell_length_a"))
  assert_true(s.contains("_atom_site_label"))
  assert_true(s.contains("C1"))
}

///|
test "crystal_sample_structure_bonds" {
  let s = @src.crystal_sample_structure()
  // 6 bonds in benzene ring
  assert_eq(s.bonds.length(), 6)
  assert_eq(s.bonds[0].atom1, "C1")
  assert_eq(s.bonds[0].atom2, "C2")
  assert_true((s.bonds[0].distance - 1.39).abs() < 0.001)
}

///|
test "orthogonalization_matrix_shape" {
  let cell = @src.UnitCell::new(
    a=5.0,
    b=6.0,
    c=7.0,
    alpha=90.0,
    beta=90.0,
    gamma=90.0,
  )
  let m = @src.orthogonalization_matrix(cell)
  assert_eq(m.length(), 3)
  assert_eq(m[0].length(), 3)
  // For orthogonal cell, m is diagonal.
  assert_true((m[0][0] - 5.0).abs() < 0.001)
  assert_true((m[1][1] - 6.0).abs() < 0.001)
  assert_true((m[2][2] - 7.0).abs() < 0.001)
}

///|
test "cif_keys_returns_all_items" {
  let block = @src.CifBlock::new(name="test")
  block.set_scalar("_a", "1")
  block.set_scalar("_b", "2")
  let keys = block.keys()
  assert_eq(keys.length(), 2)
}
