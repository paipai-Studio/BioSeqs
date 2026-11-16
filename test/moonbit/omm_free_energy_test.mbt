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

// --- A10 高层便利入口 ---

///|
test "omm_free_energy_result_fields" {
  let r = @src.OmmFreeEnergyResult::new(-1.23, "MBAR", 1000, 5, true)
  assert_true(@moonbit.assert_num_close(r.dG, -1.23, 1.0e-12))
  assert_eq(r.alg, "MBAR")
  assert_eq(r.n_samples, 1000)
  assert_eq(r.n_states, 5)
  assert_true(r.converged)
}

///|
test "omm_fep_batch_two_windows" {
  // 两窗口 FEP batch 返回 1 条 dG
  // convention: omm_fep(A_samples_in_A, A_samples_in_B) → dG_{A→B}
  let e0_in_0 = [1.0, 2.0, 3.0, 4.0] // 状态 A 样本在 A 的能量
  let e0_in_1 = [3.0, 4.0, 5.0, 6.0] // 状态 A 样本在 B 的能量 (= +2 恒偏移)
  let results = @src.omm_fep_batch([e0_in_0, e0_in_1])
  assert_eq(results.length(), 1)
  assert_eq(results[0].alg, "FEP")
  // dU = e0_in_1 - e0_in_0 = +2，FEP 算 dG = +2
  assert_true(@moonbit.assert_num_close(results[0].dG, 2.0, 1.0e-9))
}

///|
test "omm_fep_batch_single_window_returns_empty" {
  let results = @src.omm_fep_batch([[1.0, 2.0]])
  assert_eq(results.length(), 0)
}

///|
test "omm_mbar_from_energy_trajectories_two_windows_constant_shift" {
  // 两个 lambda 窗口，每个 4 个样本，窗口 1 比窗口 0 高 2.0
  let e0 = [1.0, 2.0, 3.0, 4.0]
  let e1 = [3.0, 4.0, 5.0, 6.0]
  let result = @src.omm_mbar_from_energy_trajectories(
    [e0, e1],
    lambdas=[0.0, 1.0],
    temperature=300.0,
  )
  assert_eq(result.alg, "MBAR")
  assert_eq(result.n_states, 2)
  assert_true(result.converged)
}

///|
test "omm_bar_batch_two_pairs" {
  // BAR batch 输入 k=3 个窗口 → 2 对相邻对 → 输出 2 条结果
  let eA0_A0 = [0.0, 0.0, 0.0]
  let eB0_A0 = [1.0, 1.0, 1.0]
  let eA0_B0 = [0.0, 0.0, 0.0]
  let eB0_B0 = [1.0, 1.0, 1.0]
  let eA1_A1 = [0.0, 0.0, 0.0]
  let eB1_A1 = [1.0, 1.0, 1.0]
  let eA1_B1 = [0.0, 0.0, 0.0]
  let eB1_B1 = [1.0, 1.0, 1.0]
  let eA2_A2 = [0.0, 0.0, 0.0]
  let eB2_A2 = [1.0, 1.0, 1.0]
  let eA2_B2 = [0.0, 0.0, 0.0]
  let eB2_B2 = [1.0, 1.0, 1.0]
  let results = @src.omm_bar_batch(
    [eA0_A0, eA1_A1, eA2_A2],
    [eB0_A0, eB1_A1, eB2_A2],
    [eA0_B0, eA1_B1, eA2_B2],
    [eB0_B0, eB1_B1, eB2_B2],
  )
  assert_eq(results.length(), 2)
  assert_eq(results[0].alg, "BAR")
  assert_eq(results[1].alg, "BAR")
  // BAR 恒偏移 1.0 → dG = 1.0
  assert_true(@moonbit.assert_num_close(results[0].dG, 1.0, 1.0e-3))
  assert_true(@moonbit.assert_num_close(results[1].dG, 1.0, 1.0e-3))
}

///|
test "omm_wham_from_cvs_three_windows" {
  // 三个 umbrella 窗口的简单测试
  let cvs0 = [-0.1, 0.0, 0.1, -0.05]
  let cvs1 = [0.9, 1.0, 1.1, 0.95]
  let cvs2 = [1.9, 2.0, 2.1, 1.95]
  let (bins, pmf) = @src.omm_wham_from_cvs(
    [cvs0, cvs1, cvs2],
    [0.0, 1.0, 2.0],
    [50.0, 50.0, 50.0],
    n_bins=50,
    range_min=-0.5,
    range_max=2.5,
  )
  assert_eq(bins.length(), 50)
  assert_eq(pmf.length(), 50)
}
