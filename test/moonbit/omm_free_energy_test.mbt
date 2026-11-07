///|
test "omm thermodynamic integration linear energy function" {
  // U(lambda) = 10*lambda -> dU/dlambda = 10 -> dG = 10
  let energy_fn = fn(_p : Array[Array[Double]], l : Double) -> Double {
    10.0 * l
  }
  let dg = @src.omm_thermodynamic_integration([[0.0, 0.0, 0.0]], energy_fn)
  assert_true(@moonbit.assert_num_close(dg, 10.0, 1.0e-4))
}

///|
test "omm thermodynamic integration quadratic energy function" {
  // U(lambda) = 4*lambda^2 -> dU/dlambda = 8*lambda -> dG = 4
  // central finite difference is exact for quadratics
  let energy_fn = fn(_p : Array[Array[Double]], l : Double) -> Double {
    4.0 * l * l
  }
  let dg = @src.omm_thermodynamic_integration([[0.0, 0.0, 0.0]], energy_fn)
  assert_true(@moonbit.assert_num_close(dg, 4.0, 1.0e-6))
}

///|
test "omm thermodynamic integration degenerate window count" {
  let energy_fn = fn(_p : Array[Array[Double]], l : Double) -> Double {
    10.0 * l
  }
  let dg = @src.omm_thermodynamic_integration(
    [[0.0, 0.0, 0.0]],
    energy_fn,
    n_lambda=1,
  )
  assert_eq(dg, 0.0)
}

///|
test "omm fep identical states give zero free energy" {
  let e_a = [1.0, 2.0, 3.0]
  let e_b = [1.0, 2.0, 3.0]
  let dg = @src.omm_fep(e_a, e_b)
  assert_true(@moonbit.assert_num_close(dg, 0.0, 1.0e-12))
}

///|
test "omm fep constant perturbation recovers offset" {
  // dU = c for every sample -> dG = -kT*ln(exp(-beta*c)) = c
  let e_a = [1.0, 2.0, 3.0, 4.0]
  let e_b = [3.0, 4.0, 5.0, 6.0]
  let dg = @src.omm_fep(e_a, e_b)
  assert_true(@moonbit.assert_num_close(dg, 2.0, 1.0e-9))
}

///|
test "omm fep empty sample list returns zero" {
  let dg = @src.omm_fep([], [])
  assert_eq(dg, 0.0)
}

///|
test "omm bar symmetric perturbations give zero" {
  // dU_A = +1 from A samples, dU_B = -1 from B samples, equal counts -> dG = 0
  let dG = @src.omm_bar([0.0], [1.0], [1.0], [0.0])
  assert_true(@moonbit.assert_num_close(dG, 0.0, 1.0e-3))
}

///|
test "omm bar constant shift is recovered" {
  // every sample has dU = 2 in both directions -> BAR converges to 2
  let dG = @src.omm_bar([0.0], [2.0], [0.0], [2.0])
  assert_true(@moonbit.assert_num_close(dG, 2.0, 1.0e-3))
}

///|
test "omm bar equal offsets in both ensembles give unit shift" {
  // dU_A = dU_B = 1 with equal counts -> BAR converges to 1
  let dG = @src.omm_bar([0.0], [1.0], [-1.0], [0.0])
  assert_true(@moonbit.assert_num_close(dG, 1.0, 1.0e-3))
}

///|
test "omm bar empty samples return zero" {
  let dG = @src.omm_bar([], [], [], [])
  assert_eq(dG, 0.0)
}

///|
test "omm alchemical energy endpoints reduce to pure states" {
  // lambda=0: U = u_A(r) with sigma=1, eps=1 at r=1 -> 0
  let u0 = @src.omm_alchemical_energy(1.0, 1.0, 1.0, 3.0, 2.0, 0.0)
  assert_true(@moonbit.assert_num_close(u0, 0.0, 1.0e-12))
  // lambda=1: U = u_B ; r at LJ minimum of B (sigma=1, eps=1) -> -1
  let rmin = 1.122462048309373
  let u1 = @src.omm_alchemical_energy(rmin, 1.0, 1.0, 1.0, 1.0, 1.0)
  assert_true(@moonbit.assert_num_close(u1, -1.0, 1.0e-9))
}

///|
test "omm alchemical energy symmetric in lambda for identical states" {
  let u03 = @src.omm_alchemical_energy(1.3, 1.0, 1.0, 1.0, 1.0, 0.3)
  let u07 = @src.omm_alchemical_energy(1.3, 1.0, 1.0, 1.0, 1.0, 0.7)
  assert_true(@moonbit.assert_num_close(u03, u07, 1.0e-12))
}

///|
test "omm alchemical energy softened at intermediate lambda" {
  // endpoints reduce to pure LJ at r=1 (zero) ; soft-core inflates r_eff at
  // lambda=0.5 so the pair enters the attractive well -> negative energy
  let u0 = @src.omm_alchemical_energy(1.0, 1.0, 1.0, 1.0, 1.0, 0.0)
  let u1 = @src.omm_alchemical_energy(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
  let u_mid = @src.omm_alchemical_energy(1.0, 1.0, 1.0, 1.0, 1.0, 0.5)
  assert_true(@moonbit.assert_num_close(u0, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(u1, 0.0, 1.0e-12))
  assert_true(u_mid < 0.0)
}
