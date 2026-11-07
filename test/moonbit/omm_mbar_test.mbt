///|
test "softcore_lj_limits" {
  // lambda = 0 -> no interaction.
  assert_true(@src.omm_softcore_lj(1.0, 0.2, 1.0, 0.0) == 0.0)
  // lambda = 1 with tiny alpha -> standard LJ: 4*eps*((s/r)^12 - (s/r)^6).
  let e = @src.omm_softcore_lj(1.0, 0.2, 1.0, 1.0, alpha_sc=1.0e-9)
  let standard = 4.0 * 0.2 * (1.0 - 1.0)
  ignore(standard)
  // At r = sigma the LJ energy is exactly 0 regardless of soft-core.
  assert_true(@moonbit.assert_num_close(e, 0.0, 1.0e-9))
  // No singularity: energy stays finite at r -> 0.
  let e0 = @src.omm_softcore_lj(1.0, 0.2, 1.0e-6, 0.5)
  assert_true(e0 == e0)
  assert_true(e0.abs() < 1.0e6)
  // Monotone switching: |E| grows toward the standard value as lambda rises.
  let e25 = @src.omm_softcore_lj(1.0, 0.2, 2.0, 0.25)
  let e75 = @src.omm_softcore_lj(1.0, 0.2, 2.0, 0.75)
  // At r=2sigma LJ is attractive (negative); deeper at larger lambda.
  assert_true(e75 < e25)
  assert_true(e25 < 0.0)
}

///|
test "softcore_coulomb_limits" {
  assert_true(@src.omm_softcore_coulomb(1.0, 2.0, 0.0) == 0.0)
  // lambda = 1, tiny alpha -> bare Coulomb 332.06371 * q1q2 / r.
  let e = @src.omm_softcore_coulomb(1.0, 2.0, 1.0, alpha_sc=1.0e-12)
  assert_true(@moonbit.assert_num_close(e, 332.06371 / 2.0, 1.0e-6))
  // Repulsive pair (+/+): larger lambda -> larger energy.
  let e25 = @src.omm_softcore_coulomb(1.0, 2.0, 0.25)
  let e75 = @src.omm_softcore_coulomb(1.0, 2.0, 0.75)
  assert_true(e75 > e25 && e25 > 0.0)
  // No singularity at r = 0.
  let e0 = @src.omm_softcore_coulomb(1.0, 0.0, 0.5)
  assert_true(e0 == e0)
}

///|
test "mbar_two_harmonic_states" {
  // Two 1D harmonic states with different force constants:
  // U_A = 0.5 kA x^2, U_B = 0.5 kB x^2. Exact dF = 0.5 ln(kB/kA).
  // Samples are drawn from each state's Boltzmann distribution.
  let k_a = 1.0
  let k_b = 4.0
  let n_samples = 600
  @src.md_seed_rng(12345)
  let u_nk : Array[Array[Array[Double]]] = []
  for k = 0; k < 2; k = k + 1 {
    let rows : Array[Array[Double]] = []
    for _n = 0; _n < n_samples; _n = _n + 1 {
      let kk = if k == 0 { k_a } else { k_b }
      let x = @src.md_gaussian() / @math.pow(kk, 0.5)
      rows.push([0.5 * k_a * x * x, 0.5 * k_b * x * x])
    }
    u_nk.push(rows)
  }
  let log_f = @src.omm_mbar_solve(u_nk, max_iter=500, tol=1.0e-10)
  let exact = 0.5 * @math.ln(k_b / k_a)
  assert_true(@moonbit.assert_num_close(log_f[0], 0.0, 1.0e-12))
  assert_true((log_f[1] - exact).abs() < 0.12)
}

///|
test "mbar_three_states_consistency" {
  // Three states: A (k=1), B (k=2), C (k=4). Exact f_B = 0.5 ln 2,
  // f_C = 0.5 ln 4 and f_C = f_A + f_B (path independence holds since
  // all are relative to A).
  let ks = [1.0, 2.0, 4.0]
  let n_samples = 500
  @src.md_seed_rng(777)
  let u_nk : Array[Array[Array[Double]]] = []
  for k = 0; k < 3; k = k + 1 {
    let rows : Array[Array[Double]] = []
    for _n = 0; _n < n_samples; _n = _n + 1 {
      let x = @src.md_gaussian() / @math.pow(ks[k], 0.5)
      let row : Array[Double] = []
      for j = 0; j < 3; j = j + 1 {
        row.push(0.5 * ks[j] * x * x)
      }
      rows.push(row)
    }
    u_nk.push(rows)
  }
  let log_f = @src.omm_mbar_solve(u_nk, max_iter=500, tol=1.0e-10)
  assert_true((log_f[1] - 0.5 * @math.ln(2.0)).abs() < 0.15)
  assert_true((log_f[2] - 0.5 * @math.ln(4.0)).abs() < 0.15)
}

///|
test "lambda_schedules" {
  // Uniform schedule.
  let lam = @src.omm_lambda_linear(n=5)
  assert_eq(lam.length(), 5)
  assert_true(@moonbit.assert_num_close(lam[0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(lam[4], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(lam[2], 0.5, 1.0e-12))
  // Rebalance with a constant gradient keeps the uniform grid when the
  // target equals the actual per-interval dF.
  let lams = [0.0, 0.25, 0.5, 0.75, 1.0]
  let grads = [1.0, 1.0, 1.0, 1.0, 1.0]
  let reb = @src.omm_lambda_rebalance(lams, grads, target=0.25)
  assert_true(reb.length() == 5)
  for i = 0; i < 5; i = i + 1 {
    assert_true(@moonbit.assert_num_close(reb[i], lams[i], 1.0e-9))
  }
  // Rebalance with a growing gradient concentrates points at large lambda.
  let grads2 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let reb2 = @src.omm_lambda_rebalance(lams, grads2, target=0.5)
  assert_true(reb2.length() >= 3)
  assert_true(@moonbit.assert_num_close(reb2[0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(reb2[reb2.length() - 1], 1.0, 1.0e-12))
  // Interior points shift toward high-lambda region: first interval is
  // wider than the last.
  let n2 = reb2.length()
  assert_true(reb2[1] - reb2[0] > reb2[n2 - 1] - reb2[n2 - 2])
}
