test "gro_parse_basic" {
  // Proper GRO fixed-column: 5+5+5+5+8+8+8 = 44 chars per atom line
  let content = "TITLE\n2\n    1ALA  CA        1   0.000   0.000   0.000\n    2ALA  CB        2   0.150   0.000   0.000\n1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000 1.00000\n"
  let frame = @src.gro_parse(content)
  assert_eq(frame.atoms.length(), 2)
  assert_eq(frame.atoms[0].res_name, "ALA")
  assert_eq(frame.atoms[0].atom_name, "CA")
  assert_eq(frame.atoms[1].res_id, 2)
}

test "gro_parse_coords" {
  let content = "Test\n1\n    1GLY  N        1   1.234   2.345   3.456\n1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000 1.00000\n"
  let frame = @src.gro_parse(content)
  assert_eq(frame.atoms[0].x, 1.234)
  assert_eq(frame.atoms[0].y, 2.345)
  assert_eq(frame.atoms[0].z, 3.456)
}

test "gro_parse_box" {
  let content = "Box\n1\n    1ALA  CA        1   0.000   0.000   0.000\n2.00000 0.00000 0.00000 0.00000 2.00000 0.00000 0.00000 0.00000 2.00000\n"
  let frame = @src.gro_parse(content)
  assert_true(frame.box.length() == 3)
  assert_eq(frame.box[0][0], 2.0)
  assert_eq(frame.box[1][1], 2.0)
  assert_eq(frame.box[2][2], 2.0)
}

test "gro_nm_to_angstrom" {
  let nm = [1.0, 2.0, 3.0]
  let ang = @src.gro_nm_to_angstrom(nm)
  assert_eq(ang[0], 10.0)
  assert_eq(ang[1], 20.0)
  assert_eq(ang[2], 30.0)
}

test "gro_angstrom_to_nm" {
  let ang = [10.0, 20.0, 30.0]
  let nm = @src.gro_angstrom_to_nm(ang)
  assert_eq(nm[0], 1.0)
  assert_eq(nm[1], 2.0)
  assert_eq(nm[2], 3.0)
}
