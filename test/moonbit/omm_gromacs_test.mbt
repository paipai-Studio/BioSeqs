///|
test "omm gromacs top parse atoms bonds angles dihedrals" {
  let content = "; topology\n[ atoms ]\n1 C 1 SOL OW 1 0.0 16.0\n2 C 1 SOL HW 1 0.4 1.0\n[ bonds ]\n1 2 1 300.0 0.1\n[ angles ]\n1 2 3 1 100.0 50.0\n[ dihedrals ]\n1 2 3 4 1 0.0 5.0 3\n"
  let top = @src.omm_parse_gromacs_top(content)
  assert_eq(top.n_atoms, 2)
  assert_eq(top.atom_names[0], "OW")
  assert_eq(top.atom_names[1], "HW")
  assert_eq(top.atom_types[0], "C")
  assert_eq(top.residues[0], "SOL")
  assert_eq(top.charges[0], 0.0)
  assert_eq(top.charges[1], 0.4)
  assert_eq(top.masses[0], 16.0)
  assert_eq(top.masses[1], 1.0)
  assert_eq(top.bonds.length(), 1)
  let (i, j, func, k, r0) = top.bonds[0]
  assert_eq(i, 0)
  assert_eq(j, 1)
  assert_eq(func, 1)
  assert_eq(k, 300.0)
  assert_eq(r0, 0.1)
  assert_eq(top.angles.length(), 1)
  let (ai, aj, ak, afunc, theta0, k_theta) = top.angles[0]
  assert_eq(ai, 0)
  assert_eq(aj, 1)
  assert_eq(ak, 2)
  assert_eq(afunc, 1)
  assert_eq(theta0, 100.0)
  assert_eq(k_theta, 50.0)
  assert_eq(top.dihedrals.length(), 1)
  let (d1, d2, d3, d4, dfunc, phi0, k_phi, n) = top.dihedrals[0]
  assert_eq(d1, 0)
  assert_eq(d2, 1)
  assert_eq(d3, 2)
  assert_eq(d4, 3)
  assert_eq(dfunc, 1)
  assert_eq(phi0, 0.0)
  assert_eq(k_phi, 5.0)
  assert_eq(n, 3)
}

///|
test "omm gromacs top parse skips comments and unknown sections" {
  let content = "; only a comment\n[ system ]\nwater\n[ atoms ]\n1 C 1 SOL OW 1 0.0 16.0\n"
  let top = @src.omm_parse_gromacs_top(content)
  assert_eq(top.n_atoms, 1)
  assert_eq(top.bonds.length(), 0)
}

///|
test "omm gromacs gro parse coordinates and box" {
  // standard gro fixed columns: resnr(5) resname(5) atomname(5) atomnr(5)
  // then x(8) y(8) z(8)
  let content = "test\n    1\n    1SOL     OW    1   1.000   2.000   3.000\n   10.0 10.0 10.0\n"
  let gro = @src.omm_parse_gro(content)
  assert_eq(gro.n_atoms, 1)
  assert_eq(gro.coords.length(), 1)
  let (x, y, z) = gro.coords[0]
  assert_true(@moonbit.assert_num_close(x, 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(y, 2.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(z, 3.0, 1.0e-9))
  let (vx, vy, vz) = gro.velocities[0]
  assert_eq(vx, 0.0)
  assert_eq(vy, 0.0)
  assert_eq(vz, 0.0)
  let (bx, by, bz) = gro.box
  assert_eq(bx, 10.0)
  assert_eq(by, 10.0)
  assert_eq(bz, 10.0)
}

///|
test "omm gromacs gro parse with velocities" {
  let content = "test\n    1\n    1SOL     OW    1   1.000   2.000   3.000  0.5000  1.0000  2.0000\n   10.0 10.0 10.0\n"
  let gro = @src.omm_parse_gro(content)
  let (vx, vy, vz) = gro.velocities[0]
  assert_true(@moonbit.assert_num_close(vx, 0.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(vy, 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(vz, 2.0, 1.0e-9))
}

///|
test "omm gromacs gro parse multiple atoms" {
  let content = "water\n    2\n    1SOL     OW    1   0.000   0.000   0.000\n    1SOL     HW    2   1.000   0.000   0.000\n   20.0 20.0 20.0\n"
  let gro = @src.omm_parse_gro(content)
  assert_eq(gro.n_atoms, 2)
  assert_eq(gro.coords.length(), 2)
  let (x1, _y1, _z1) = gro.coords[0]
  let (x2, _y2, _z2) = gro.coords[1]
  assert_true(@moonbit.assert_num_close(x1, 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(x2, 1.0, 1.0e-9))
}

///|
test "omm gromacs bond energy at stretched pair" {
  // implementation compares squared distance: E = 0.5*k*(r^2 - r0)^2
  // r=2 -> r^2=4 = r0 -> zero energy
  let top = @src.omm_parse_gromacs_top(
    "[ atoms ]\n1 C 1 SOL OW 1 0.0 16.0\n2 C 1 SOL HW 1 0.4 1.0\n[ bonds ]\n1 2 1 300.0 4.0\n",
  )
  let e = @src.omm_gromacs_bond_energy(top, [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
  assert_true(@moonbit.assert_num_close(e, 0.0, 1.0e-12))
}

///|
test "omm gromacs bond energy nonzero displacement" {
  // r=2, r0=1.5 -> r^2=4, dr=2.5, E = 0.5*300*6.25 = 937.5
  let top = @src.omm_parse_gromacs_top(
    "[ atoms ]\n1 C 1 SOL OW 1 0.0 16.0\n2 C 1 SOL HW 1 0.4 1.0\n[ bonds ]\n1 2 1 300.0 1.5\n",
  )
  let e = @src.omm_gromacs_bond_energy(top, [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
  assert_true(@moonbit.assert_num_close(e, 937.5, 1.0e-9))
}

///|
test "omm gromacs bond energy ignores out of range indices" {
  let top = @src.omm_parse_gromacs_top(
    "[ atoms ]\n1 C 1 SOL OW 1 0.0 16.0\n2 C 1 SOL HW 1 0.4 1.0\n[ bonds ]\n1 9 1 300.0 1.5\n",
  )
  let e = @src.omm_gromacs_bond_energy(top, [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
  assert_eq(e, 0.0)
}
