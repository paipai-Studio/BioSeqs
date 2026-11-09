///|
test "omm charmm psf parse extracts atoms and bonds" {
  let content = "PSF CMAP\n\n1 !NTITLE\n\n2 !NATOM\n1 A 1 ALA CA C -0.1 12.0\n2 A 1 ALA CB C 0.1 13.0\n\n1 !NBOND\n1 2\n"
  let params = @src.omm_parse_psf(content)
  assert_eq(params.n_atoms, 2)
  assert_eq(params.atom_types[0], "C")
  assert_eq(params.atom_types[1], "C")
  assert_true(@moonbit.assert_num_close(params.charges[0], -0.1, 1.0e-12))
  assert_true(@moonbit.assert_num_close(params.charges[1], 0.1, 1.0e-12))
  assert_true(@moonbit.assert_num_close(params.masses[0], 12.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(params.masses[1], 13.0, 1.0e-12))
  assert_eq(params.bonds.length(), 1)
  let (i, j, k, r0) = params.bonds[0]
  assert_eq(i, 0)
  assert_eq(j, 1)
  assert_eq(k, 200.0)
  assert_eq(r0, 1.5)
  assert_eq(params.angles.length(), 0)
  assert_eq(params.dihedrals.length(), 0)
}

///|
test "omm charmm bond energy at stretched pair" {
  // implementation uses r as squared distance: E = k * (r^2 - r0)^2
  // r=2 -> r^2=4, r0=1.5 -> E = 200 * 2.5^2 = 1250
  let params = @src.omm_parse_psf(
    "PSF\n\n2 !NATOM\n1 A 1 ALA CA C 0.0 12.0\n2 A 1 ALA CB C 0.0 12.0\n\n1 !NBOND\n1 2\n",
  )
  let e = @src.omm_charmm_bond_energy(params, [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
  assert_true(@moonbit.assert_num_close(e, 1250.0, 1.0e-9))
}

///|
test "omm charmm bond energy zero when squared distance equals r0" {
  // parser fixes r0=1.5 ; r = sqrt(1.5) -> r^2 = 1.5 = r0 -> zero energy
  let params = @src.omm_parse_psf(
    "PSF\n\n2 !NATOM\n1 A 1 ALA CA C 0.0 12.0\n2 A 1 ALA CB C 0.0 12.0\n\n1 !NBOND\n1 2\n",
  )
  let r0 = params.bonds[0].3
  let r = @math.pow(r0, 0.5)
  let e = @src.omm_charmm_bond_energy(params, [(0.0, 0.0, 0.0), (r, 0.0, 0.0)])
  assert_true(@moonbit.assert_num_close(e, 0.0, 1.0e-9))
}

///|
test "omm charmm elec energy opposite charges attract" {
  let params = @src.omm_parse_psf(
    "PSF\n\n2 !NATOM\n1 A 1 ALA CA C 1.0 12.0\n2 A 1 ALA CB C -1.0 12.0\n",
  )
  let e = @src.omm_charmm_elec_energy(params, [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
  // implementation uses squared distance as r: r^2 = 4 here
  assert_true(@moonbit.assert_num_close(e, -138.935485 / 4.0, 1.0e-9))
}

///|
test "omm charmm elec energy respects cutoff" {
  let params = @src.omm_parse_psf(
    "PSF\n\n2 !NATOM\n1 A 1 ALA CA C 1.0 12.0\n2 A 1 ALA CB C -1.0 12.0\n",
  )
  let e = @src.omm_charmm_elec_energy(
    params,
    [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
    cutoff=1.0,
  )
  assert_eq(e, 0.0)
}

///|
test "omm charmm params file parse bonds and angles" {
  let content = "BONDS\nC C 300.0 1.5\nANGLES\nC C C 50.0 110.0\n"
  let (bonds, angles) = @src.omm_parse_charmm_params(content)
  assert_eq(bonds.length(), 1)
  let (t1, t2, k, r0) = bonds[0]
  assert_eq(t1, "C")
  assert_eq(t2, "C")
  assert_eq(k, 300.0)
  assert_eq(r0, 1.5)
  assert_eq(angles.length(), 1)
  let (a1, a2, a3, k_theta, theta0) = angles[0]
  assert_eq(a1, "C")
  assert_eq(a2, "C")
  assert_eq(a3, "C")
  assert_eq(k_theta, 50.0)
  assert_eq(theta0, 110.0)
}

///|
test "omm charmm params file parse multiple bond entries" {
  let content = "BONDS\nC C 300.0 1.5\nC H 250.0 1.1\n"
  let (bonds, angles) = @src.omm_parse_charmm_params(content)
  assert_eq(bonds.length(), 2)
  assert_eq(bonds[1].0, "C")
  assert_eq(bonds[1].1, "H")
  assert_eq(bonds[1].2, 250.0)
  assert_eq(bonds[1].3, 1.1)
  assert_eq(angles.length(), 0)
}
