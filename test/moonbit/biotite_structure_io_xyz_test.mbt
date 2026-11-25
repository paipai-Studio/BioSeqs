///|
/// Tests for biotite_structure_io_xyz.mbt.
test "biotite_xyz_parse" {
  let text = "2\ncomment\nC 0.0 0.0 0.0\nO 1.5 2.0 3.0\n"
  let arr = @src.biotite_xyz_parse(text)
  assert_eq(arr.len(), 2)
  assert_eq(arr.get_name(0), "C")
  assert_true(@moonbit.assert_num_close(arr.get_coord(1).x, 1.5, 1.0e-12))
  assert_true(@moonbit.assert_num_close(arr.get_coord(1).y, 2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(arr.get_coord(1).z, 3.0, 1.0e-12))
}

///|
test "biotite_xyz_write_roundtrip" {
  let arr = @src.AtomArray::from_atoms([
    @src.Atom::new(
      name="C1",
      coord=@src.Vector3::new(0.0, 0.0, 0.0),
      element="C",
    ),
    @src.Atom::new(
      name="O1",
      coord=@src.Vector3::new(1.0, 2.0, 3.0),
      element="O",
    ),
  ])
  let text = @src.biotite_xyz_write(arr)
  assert_true(text.has_prefix("2\n"))
  let parsed = @src.biotite_xyz_parse(text)
  assert_eq(parsed.len(), 2)
  assert_eq(parsed.get_name(0), "C")
}

///|
test "biotite_xyz_short_input" {
  assert_eq(@src.biotite_xyz_parse("").len(), 0)
  assert_eq(@src.biotite_xyz_parse("3\nonly-comment\n").len(), 0)
}
