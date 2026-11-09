///|
test "omm_create_replicas_geometric_temperatures" {
  let pos = [[0.0, 0.0, 0.0]]
  let reps = @src.omm_create_replicas(4, 300.0, 3000.0, pos)
  assert_eq(reps.length(), 4)
  assert_true(@moonbit.assert_num_close(reps[0].temperature, 300.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(reps[3].temperature, 3000.0, 1.0e-9))
  let t1 = 300.0 * @math.pow(10.0, 1.0 / 3.0)
  assert_true(@moonbit.assert_num_close(reps[1].temperature, t1, 1.0e-9))
  let t2 = 300.0 * @math.pow(10.0, 2.0 / 3.0)
  assert_true(@moonbit.assert_num_close(reps[2].temperature, t2, 1.0e-9))
  assert_eq(reps[2].positions.length(), 1)
  assert_true(@moonbit.assert_num_close(reps[2].positions[0][0], 0.0, 1.0e-12))
  // a single replica sits at t_min
  let single = @src.omm_create_replicas(1, 310.0, 500.0, pos)
  assert_eq(single.length(), 1)
  assert_true(@moonbit.assert_num_close(single[0].temperature, 310.0, 1.0e-12))
  // non-positive count yields no replicas
  assert_eq(@src.omm_create_replicas(0, 300.0, 3000.0, pos).length(), 0)
}

///|
test "omm_attempt_replica_swap_with_equal_energies" {
  // OmmReplica carries energy = 0 here, so the Metropolis delta
  // (beta_i - beta_j) * (E_j - E_i) is exactly 0 and exp(0) = 1 > rand in [0,1):
  // every attempt must be accepted, in both directions and at equal temperature.
  let pos = [[0.0, 0.0, 0.0]]
  let reps = @src.omm_create_replicas(2, 300.0, 600.0, pos)
  assert_true(@src.omm_attempt_replica_swap(reps[0], reps[1]))
  assert_true(@src.omm_attempt_replica_swap(reps[1], reps[0]))
  let same = @src.omm_create_replicas(2, 400.0, 400.0, pos)
  assert_true(@src.omm_attempt_replica_swap(same[0], same[1]))
}

///|
test "omm_replica_exchange_round_swaps_temperatures" {
  let pos = [[0.0, 0.0, 0.0]]
  let reps = @src.omm_create_replicas(2, 300.0, 600.0, pos)
  let swaps = @src.omm_replica_exchange_round(reps)
  assert_eq(swaps.length(), 1)
  assert_true(swaps[0])
  // configurations stay, temperatures swap
  assert_true(@moonbit.assert_num_close(reps[0].temperature, 600.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(reps[1].temperature, 300.0, 1.0e-12))
  // odd parity with only two replicas attempts no pair
  let reps2 = @src.omm_create_replicas(2, 300.0, 600.0, pos)
  let swaps2 = @src.omm_replica_exchange_round(reps2, parity=1)
  assert_eq(swaps2.length(), 1)
  assert_false(swaps2[0])
  assert_true(@moonbit.assert_num_close(reps2[0].temperature, 300.0, 1.0e-12))
}

///|
test "omm_run_remd_identity_step_counts_attempts" {
  let pos = [[0.0, 0.0, 0.0]]
  let reps = @src.omm_create_replicas(2, 300.0, 600.0, pos)
  let identity = fn(r : @src.OmmReplica, _dt : Double) -> @src.OmmReplica { r }
  // attempts at step 100 (odd parity: no pair) and step 200 (even: swap)
  let (final_reps, n_swaps, n_attempts) = @src.omm_run_remd(
    reps,
    201,
    swap_interval=100,
    identity,
  )
  assert_eq(n_attempts, 2)
  assert_eq(n_swaps, 1)
  assert_true(
    @moonbit.assert_num_close(final_reps[0].temperature, 600.0, 1.0e-12),
  )
  assert_true(
    @moonbit.assert_num_close(final_reps[1].temperature, 300.0, 1.0e-12),
  )
}

///|
test "omm_optimal_temperatures_and_acceptance_rate" {
  let temps = @src.omm_optimal_temperatures(3, 300.0, 100)
  assert_eq(temps.length(), 3)
  assert_true(@moonbit.assert_num_close(temps[0], 300.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(temps[1], 330.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(temps[2], 363.0, 1.0e-9))
  assert_eq(@src.omm_optimal_temperatures(0, 300.0, 100).length(), 0)
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_swap_acceptance_rate(5, 10),
      0.5,
      1.0e-12,
    ),
  )
  assert_true(
    @moonbit.assert_num_close(@src.omm_swap_acceptance_rate(0, 0), 0.0, 1.0e-12),
  )
}
