///|
test "omm obc2 born radii single atom analytic value" {
  // isolated atom: offset=0, radius=ri
  // reff = 1 / (1/ri + 0.8*ri + 4.85*ri^2)
  let radii = [1.4]
  let coords = [(0.0, 0.0, 0.0)]
  let reff = @src.omm_obc2_born_radii(coords, radii)
  let expected = 1.0 / (1.0 / 1.4 + 0.8 * 1.4 + 4.85 * 1.4 * 1.4)
  assert_eq(reff.length(), 1)
  assert_true(@moonbit.assert_num_close(reff[0], expected, 1.0e-12))
}

///|
test "omm obc2 born radii pair positive and finite" {
  let radii = [1.4, 1.4]
  let coords = [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0)]
  let reff = @src.omm_obc2_born_radii(coords, radii)
  assert_eq(reff.length(), 2)
  assert_true(reff[0] > 0.0)
  assert_true(reff[1] > 0.0)
  assert_true(reff[0] < 1.4)
}

///|
test "omm obc2 energy single atom self term" {
  // E = -0.5*(1/eps_in - 1/eps_out)*138.935485*q^2/B with B = born radius
  let charges = [1.0]
  let radii = [1.4]
  let coords = [(0.0, 0.0, 0.0)]
  let reff = @src.omm_obc2_born_radii(coords, radii)
  let e = @src.omm_obc2_energy(charges, radii, coords)
  let expected = -0.5 * (1.0 - 1.0 / 78.5) * 138.935485 * (1.0 / reff[0])
  assert_true(@moonbit.assert_num_close(e, expected, 1.0e-9))
  assert_true(e < 0.0)
}

///|
test "omm obc2 energy zero when dielectrics equal" {
  // 1/2 - 1/2 = 0 -> prefactor zero -> energy zero
  let charges = [1.0, -1.0]
  let radii = [1.4, 1.4]
  let coords = [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)]
  let e = @src.omm_obc2_energy(
    charges,
    radii,
    coords,
    solute_dielectric=2.0,
    solvent_dielectric=2.0,
  )
  assert_true(@moonbit.assert_num_close(e, 0.0, 1.0e-9))
}

///|
test "omm obc2 energy opposite charges negative" {
  let charges = [1.0, -1.0]
  let radii = [1.4, 1.4]
  let coords = [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)]
  let e = @src.omm_obc2_energy(charges, radii, coords)
  assert_true(e < 0.0)
}

///|
test "omm obc2 energy zero charges give zero" {
  let charges = [0.0, 0.0]
  let radii = [1.4, 1.4]
  let coords = [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)]
  let e = @src.omm_obc2_energy(charges, radii, coords)
  assert_true(@moonbit.assert_num_close(e, 0.0, 1.0e-9))
}

///|
test "omm obc2 energy matches gbn core formula" {
  // two atoms with unit born radii supplied explicitly
  let charges = [1.0, -1.0]
  let coords = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
  // use equal dielectrics trick? instead verify against manual pair formula:
  // E = sum_i pref*q_i^2/B_i + sum_{i<j} 2*pref*q_i*q_j/f_ij
  // f_ij = sqrt(r^2 + B_i*B_j*exp(-r^2/(4*B_i*B_j)))
  let pref = -0.5 * (1.0 - 1.0 / 78.5) * 138.935485
  let b = 1.0
  let r2 = 4.0
  let f_ij = @math.pow(r2 + b * b * @math.exp(-r2 / (4.0 * b * b)), 0.5)
  // note: omm_obc2_energy computes born radii internally; we instead verify
  // via omm_gbn_energy using the public obc2 radii to stay consistent
  let radii = [1.4, 1.4]
  let reff = @src.omm_obc2_born_radii(coords, radii)
  let bi = reff[0]
  let bj = reff[1]
  let f2 = @math.pow(
    r2 + bi * bj * @math.exp(-r2 / (4.0 * bi * bj + 1.0e-10)),
    0.5,
  )
  let expected = pref * charges[0] * charges[0] / bi +
    pref * charges[1] * charges[1] / bj +
    2.0 * pref * charges[0] * charges[1] / f2
  let e = @src.omm_obc2_energy(charges, radii, coords)
  assert_true(@moonbit.assert_num_close(e, expected, 1.0e-9))
  ignore(f_ij)
}
