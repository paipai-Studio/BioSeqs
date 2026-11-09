///|
/// Tests for PQR (topology + charge + radius) reading and writing.

///|
test "pqr_parse_whitespace" {
  let text = "REMARKS test pqr\nATOM      1  N   ALA A   1      -0.415   7.123  12.345  -0.415700  14.0100\nATOM      2  CA  ALA A   1       0.021   8.234  11.111   0.021300  12.0110\nTER\nEND\n"
  let s = @src.pqr_parse(text)
  assert_eq(s.natoms(), 2)
  let a = s.atoms[0]
  assert_eq(a.serial, 1)
  assert_eq(a.name, "N")
  assert_eq(a.resname, "ALA")
  assert_eq(a.chain, "A")
  assert_eq(a.resid, 1)
  assert_true(@moonbit.assert_num_close(a.x, -0.415, 1.0e-9))
  assert_true(@moonbit.assert_num_close(a.y, 7.123, 1.0e-9))
  assert_true(@moonbit.assert_num_close(a.z, 12.345, 1.0e-9))
  assert_true(@moonbit.assert_num_close(a.charge, -0.4157, 1.0e-9))
  assert_true(@moonbit.assert_num_close(a.radius, 14.01, 1.0e-9))
  let b = s.atoms[1]
  assert_eq(b.name, "CA")
  assert_eq(b.resid, 1)
  assert_true(@moonbit.assert_num_close(b.x, 0.021, 1.0e-9))
  assert_true(@moonbit.assert_num_close(b.charge, 0.0213, 1.0e-9))
}

///|
test "pqr_total_charge" {
  let text = "ATOM      1  N   ALA A   1      -0.415   7.123  12.345  -0.415700  14.0100\nATOM      2  CA  ALA A   1       0.021   8.234  11.111   0.021300  12.0110\n"
  let s = @src.pqr_parse(text)
  // -0.4157 + 0.0213 = -0.3944
  assert_true(@moonbit.assert_num_close(s.total_charge(), -0.3944, 1.0e-9))
}

///|
test "pqr_write_roundtrip" {
  let atoms : Array[@src.PqrAtom] = [
    @src.PqrAtom::{
      serial: 1,
      name: "N",
      resname: "ALA",
      chain: "A",
      resid: 1,
      x: -0.415,
      y: 7.123,
      z: 12.345,
      charge: -0.4157,
      radius: 1.55,
    },
    @src.PqrAtom::{
      serial: 2,
      name: "CA",
      resname: "ALA",
      chain: "A",
      resid: 1,
      x: 0.021,
      y: 8.234,
      z: 11.111,
      charge: 0.0213,
      radius: 1.7,
    },
  ]
  let s = @src.PqrStructure::new(atoms)
  let text = @src.pqr_write(s)
  let back = @src.pqr_parse(text)
  assert_eq(back.natoms(), 2)
  let a = back.atoms[0]
  assert_eq(a.serial, 1)
  assert_eq(a.name, "N")
  assert_eq(a.resname, "ALA")
  assert_eq(a.chain, "A")
  assert_eq(a.resid, 1)
  assert_true(@moonbit.assert_num_close(a.x, -0.415, 1.0e-9))
  assert_true(@moonbit.assert_num_close(a.y, 7.123, 1.0e-9))
  assert_true(@moonbit.assert_num_close(a.z, 12.345, 1.0e-9))
  assert_true(@moonbit.assert_num_close(a.charge, -0.4157, 1.0e-9))
  assert_true(@moonbit.assert_num_close(a.radius, 1.55, 1.0e-9))
  assert_true(@moonbit.assert_num_close(back.atoms[1].z, 11.111, 1.0e-9))
  assert_true(@moonbit.assert_num_close(back.atoms[1].radius, 1.7, 1.0e-9))
}

///|
test "pqr_bounds_and_coords" {
  let atoms : Array[@src.PqrAtom] = [
    @src.PqrAtom::{
      serial: 1,
      name: "N",
      resname: "ALA",
      chain: "A",
      resid: 1,
      x: -1.0,
      y: 2.0,
      z: 0.5,
      charge: 0.0,
      radius: 1.0,
    },
    @src.PqrAtom::{
      serial: 2,
      name: "O",
      resname: "HOH",
      chain: "B",
      resid: 2,
      x: 3.0,
      y: -4.0,
      z: 2.5,
      charge: 0.0,
      radius: 1.0,
    },
  ]
  let s = @src.PqrStructure::new(atoms)
  let (mins, maxs) = s.bounds()
  assert_eq(mins.length(), 3)
  assert_eq(maxs.length(), 3)
  assert_true(@moonbit.assert_num_close(mins[0], -1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(mins[1], -4.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(mins[2], 0.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(maxs[0], 3.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(maxs[1], 2.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(maxs[2], 2.5, 1.0e-9))
  let coords = @src.pqr_coords(s)
  assert_eq(coords.length(), 2)
  assert_true(@moonbit.assert_num_close(coords[1][2], 2.5, 1.0e-9))
}

///|
test "pqr_parse_no_chain" {
  // 10 whitespace fields: no chain column between resname and resid.
  let text = "ATOM      1  N   ALA   1      -0.415   7.123  12.345  -0.415700  14.0100\n"
  let s = @src.pqr_parse(text)
  assert_eq(s.natoms(), 1)
  let a = s.atoms[0]
  assert_eq(a.chain, "")
  assert_eq(a.resid, 1)
  assert_eq(a.resname, "ALA")
  assert_true(@moonbit.assert_num_close(a.x, -0.415, 1.0e-9))
  assert_true(@moonbit.assert_num_close(a.charge, -0.4157, 1.0e-9))
  assert_true(@moonbit.assert_num_close(a.radius, 14.01, 1.0e-9))
}

///|
test "pqr_parse_missing_columns_raises" {
  // Neither 10 whitespace fields nor 54 fixed columns available.
  let text = "ATOM      1  N   ALA\n"
  let raised = try {
    let _ = @src.pqr_parse(text)
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}
