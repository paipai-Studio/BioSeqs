///|
/// Tests for the LAMMPS dump text reader (structure_lammps.mbt).
test "lammps_dump_read_orthogonal" {
  let text = "ITEM: TIMESTEP\n1000\nITEM: NUMBER OF ATOMS\n4\nITEM: BOX BOUNDS\n0.0000000000000000E+00 1.0000000000000000E+01\n0.0000000000000000E+00 2.0000000000000000E+01\n0.0000000000000000E+00 3.0000000000000000E+01\nITEM: ATOMS id mol type q x y z ix iy iz\n1 1 1 0.5 0.5 0.5 0.5 0 0 0\n2 1 2 0.5 1.5 2.5 3.5 0 0 1\n3 2 1 -0.5 4.5 5.5 6.5 1 0 0\n4 2 2 0.5 7.5 8.5 9.5 1 1 1\n"
  let frames = @src.lammps_dump_read(text)
  assert_eq(frames.length(), 1)
  let f = frames[0]
  assert_eq(f.timestep, 1000)
  assert_eq(f.natoms, 4)
  assert_eq(f.box_lo.length(), 3)
  assert_eq(f.box_hi.length(), 3)
  assert_true(@moonbit.assert_num_close(f.box_hi[0], 10.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(f.box_hi[2], 30.0, 1.0e-9))
  assert_eq(f.tilt.length(), 0)
  assert_eq(f.atoms.length(), 4)
  assert_eq(f.atoms[0].id, 1)
  assert_eq(f.atoms[0].mol, 1)
  assert_eq(f.atoms[1].atom_type, 2)
  assert_true(@moonbit.assert_num_close(f.atoms[1].y, 2.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(f.atoms[3].z, 9.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(f.atoms[2].x, 4.5, 1.0e-9))
  // q + ix + iy + iz end up in extra
  assert_eq(f.atoms[0].extra.length(), 4)
  assert_true(@moonbit.assert_num_close(f.atoms[0].extra[0], 0.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(f.atoms[2].extra[0], -0.5, 1.0e-9))
}

///|
test "lammps_dump_read_triclinic" {
  let text = "ITEM: TIMESTEP\n5\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS xy xz yz\n0.0000000000000000E+00 1.0000000000000000E+01 5.0000000000000000E-01\n0.0000000000000000E+00 2.0000000000000000E+01 6.0000000000000000E-01\n0.0000000000000000E+00 3.0000000000000000E+01 7.0000000000000000E-01\nITEM: ATOMS id type x y z\n1 1 0.5 0.5 0.5\n2 2 1.5 2.5 3.5\n"
  let frames = @src.lammps_dump_read(text)
  assert_eq(frames.length(), 1)
  let f = frames[0]
  assert_eq(f.timestep, 5)
  assert_eq(f.natoms, 2)
  assert_eq(f.tilt.length(), 3)
  assert_true(@moonbit.assert_num_close(f.tilt[0], 0.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(f.tilt[1], 0.6, 1.0e-9))
  assert_true(@moonbit.assert_num_close(f.tilt[2], 0.7, 1.0e-9))
  assert_true(@moonbit.assert_num_close(f.box_hi[1], 20.0, 1.0e-9))
  // mol is absent, defaults to 0
  assert_eq(f.atoms[0].mol, 0)
  assert_eq(f.atoms[1].atom_type, 2)
  assert_true(@moonbit.assert_num_close(f.atoms[1].x, 1.5, 1.0e-9))
  assert_eq(f.atoms[0].extra.length(), 0)
}

///|
test "lammps_dump_read_multiple_frames" {
  let text = "ITEM: TIMESTEP\n100\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS\n0.0000000000000000E+00 1.0000000000000000E+00\n0.0000000000000000E+00 1.0000000000000000E+00\n0.0000000000000000E+00 1.0000000000000000E+00\nITEM: ATOMS id x y z\n1 0.1 0.2 0.3\n\nITEM: TIMESTEP\n200\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS\n0.0000000000000000E+00 2.0000000000000000E+00\n0.0000000000000000E+00 2.0000000000000000E+00\n0.0000000000000000E+00 2.0000000000000000E+00\nITEM: ATOMS id x y z\n1 1.1 1.2 1.3\n2 1.4 1.5 1.6\n"
  let frames = @src.lammps_dump_read(text)
  assert_eq(frames.length(), 2)
  assert_eq(frames[0].timestep, 100)
  assert_eq(frames[0].natoms, 1)
  assert_eq(frames[1].timestep, 200)
  assert_eq(frames[1].natoms, 2)
  assert_true(@moonbit.assert_num_close(frames[0].atoms[0].z, 0.3, 1.0e-9))
  assert_true(@moonbit.assert_num_close(frames[1].atoms[1].z, 1.6, 1.0e-9))
  assert_eq(frames[1].atoms[1].id, 2)
}

///|
test "lammps_dump_read_scaled_and_minimal_columns" {
  let text = "ITEM: TIMESTEP\n7\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS\n0.0000000000000000E+00 1.0000000000000000E+01\n0.0000000000000000E+00 1.0000000000000000E+01\n0.0000000000000000E+00 1.0000000000000000E+01\nITEM: ATOMS id type xs ys zs\n7 3 0.05 0.06 0.07\n8 1 0.15 0.16 0.17\n"
  let frames = @src.lammps_dump_read(text)
  assert_eq(frames.length(), 1)
  let f = frames[0]
  // xs/ys/zs are treated as cartesian coordinates
  assert_true(@moonbit.assert_num_close(f.atoms[0].x, 0.05, 1.0e-9))
  assert_true(@moonbit.assert_num_close(f.atoms[1].z, 0.17, 1.0e-9))
  assert_eq(f.atoms[0].id, 7)
  assert_eq(f.atoms[0].atom_type, 3)
  assert_eq(f.atoms[0].mol, 0)
  assert_eq(f.atoms[0].extra.length(), 0)
}

///|
test "lammps_to_trajectory" {
  let text = "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS\n0.0000000000000000E+00 1.0000000000000000E+01\n0.0000000000000000E+00 1.0000000000000000E+01\n0.0000000000000000E+00 1.0000000000000000E+01\nITEM: ATOMS id mol type x y z\n1 1 1 0.5 0.5 0.5\n2 1 2 1.5 2.5 3.5\nITEM: TIMESTEP\n2\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS\n0.0000000000000000E+00 1.0000000000000000E+01\n0.0000000000000000E+00 1.0000000000000000E+01\n0.0000000000000000E+00 1.0000000000000000E+01\nITEM: ATOMS id mol type x y z\n1 1 1 1.5 1.5 1.5\n2 1 2 2.5 3.5 4.5\n"
  let frames = @src.lammps_dump_read(text)
  let traj = @src.lammps_to_trajectory(frames)
  assert_eq(traj.num_frames(), 2)
  assert_eq(traj.num_atoms(), 2)
  let frame0 = traj.get_frame(0)
  assert_true(@moonbit.assert_num_close(frame0[0][0], 0.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(frame0[1][2], 3.5, 1.0e-9))
  let frame1 = traj.get_frame(1)
  assert_true(@moonbit.assert_num_close(frame1[1][0], 2.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(frame1[1][2], 4.5, 1.0e-9))
}

///|
test "lammps_dump_read_rejects_missing_coords" {
  let text = "ITEM: TIMESTEP\n1\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS\n0.0000000000000000E+00 1.0000000000000000E+00\n0.0000000000000000E+00 1.0000000000000000E+00\n0.0000000000000000E+00 1.0000000000000000E+00\nITEM: ATOMS id type q\n1 1 0.5\n"
  let mut raised = false
  try {
    let r = @src.lammps_dump_read(text)
    ignore(r)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}
