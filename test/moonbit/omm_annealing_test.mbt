///|
test "omm anneal linear schedule endpoints and midpoint" {
  let s = @src.omm_linear_schedule(300.0, 100.0, 10)
  assert_true(
    @moonbit.assert_num_close(@src.omm_anneal_temperature(s, 0), 300.0, 1.0e-12),
  )
  assert_true(
    @moonbit.assert_num_close(@src.omm_anneal_temperature(s, 5), 200.0, 1.0e-12),
  )
  assert_true(
    @moonbit.assert_num_close(@src.omm_anneal_temperature(s, 9), 120.0, 1.0e-12),
  )
}

///|
test "omm anneal linear schedule clamps beyond last step" {
  let s = @src.omm_linear_schedule(300.0, 100.0, 10)
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_anneal_temperature(s, 10),
      100.0,
      1.0e-12,
    ),
  )
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_anneal_temperature(s, 100),
      100.0,
      1.0e-12,
    ),
  )
}

///|
test "omm anneal exponential schedule geometric interpolation" {
  // T(2) with T0=100, T1=400, n=4: 100 * 4^(2/4) = 200
  let s = @src.omm_exponential_schedule(100.0, 400.0, 4)
  assert_true(
    @moonbit.assert_num_close(@src.omm_anneal_temperature(s, 0), 100.0, 1.0e-12),
  )
  assert_true(
    @moonbit.assert_num_close(@src.omm_anneal_temperature(s, 2), 200.0, 1.0e-9),
  )
  assert_true(
    @moonbit.assert_num_close(@src.omm_anneal_temperature(s, 4), 400.0, 1.0e-12),
  )
  assert_true(
    @moonbit.assert_num_close(@src.omm_anneal_temperature(s, 9), 400.0, 1.0e-12),
  )
}

///|
test "omm anneal exponential schedule nonpositive start falls back to end" {
  let s = @src.omm_exponential_schedule(0.0, 250.0, 5)
  assert_true(
    @moonbit.assert_num_close(@src.omm_anneal_temperature(s, 0), 250.0, 1.0e-12),
  )
}

///|
test "omm simulated annealing runs all steps and is deterministic" {
  let force_fn = fn(p : Array[Array[Double]]) -> Array[Array[Double]] {
    [[-p[0][0], -p[0][1], -p[0][2]]]
  }
  let schedule = @src.omm_linear_schedule(50.0, 10.0, 20)
  @src.md_seed_rng(7L)
  let (pos1, e1) = @src.omm_simulated_annealing(
    [[1.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    force_fn,
    [1.0],
    0.001,
    schedule,
  )
  @src.md_seed_rng(7L)
  let (pos2, e2) = @src.omm_simulated_annealing(
    [[1.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    force_fn,
    [1.0],
    0.001,
    schedule,
  )
  assert_eq(pos1.length(), 1)
  assert_eq(pos1[0].length(), 3)
  assert_eq(e1.length(), 20)
  assert_true(@moonbit.assert_num_close(pos1[0][0], pos2[0][0], 1.0e-12))
  assert_true(@moonbit.assert_num_close(e1[0], e2[0], 1.0e-12))
  assert_true(@moonbit.assert_num_close(e1[19], e2[19], 1.0e-12))
}

///|
test "omm simulated annealing exponential schedule counts total steps" {
  let force_fn = fn(p : Array[Array[Double]]) -> Array[Array[Double]] {
    [[-p[0][0], -p[0][1], -p[0][2]]]
  }
  let schedule = @src.omm_exponential_schedule(100.0, 10.0, 7)
  @src.md_seed_rng(11L)
  let (_pos, energies) = @src.omm_simulated_annealing(
    [[0.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    force_fn,
    [1.0],
    0.001,
    schedule,
  )
  assert_eq(energies.length(), 7)
}
