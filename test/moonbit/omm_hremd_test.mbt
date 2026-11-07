///|
test "hremd_log_accept_symmetric" {
  // Identical Hamiltonians (pure temperature REMD with equal temps):
  // swapping identical energies gives Delta = 0 -> log_prob = 0 -> accept.
  let lp = @src.omm_hremd_log_accept(2.0, 2.0, 2.0, 2.0, 1.0, 1.0)
  assert_true(@moonbit.assert_num_close(lp, 0.0, 1.0e-12))
  assert_true(@src.omm_hremd_attempt(2.0, 2.0, 2.0, 2.0, 1.0, 1.0, uniform=0.5))
}

///|
test "hremd_temperature_remd_matches_known_formula" {
  // Same Hamiltonian U, different betas:
  // Delta = (b_i - b_j)(U_j - U_i); verify the sign and magnitude.
  let lp = @src.omm_hremd_log_accept(1.0, 3.0, 1.0, 3.0, 2.0, 1.0)
  // Delta = 2*(3-1) + 1*(1-3) = 2 -> log_p = -2.
  assert_true(@moonbit.assert_num_close(lp, -2.0, 1.0e-12))
  // Uniform above exp(-2) rejects, below accepts.
  assert_false(
    @src.omm_hremd_attempt(1.0, 3.0, 1.0, 3.0, 2.0, 1.0, uniform=0.9),
  )
  assert_true(
    @src.omm_hremd_attempt(1.0, 3.0, 1.0, 3.0, 2.0, 1.0, uniform=0.01),
  )
}

///|
test "hremd_swap_cycle_occupancy" {
  // Two replicas, two configs: replica 0 prefers config 1 and replica 1
  // prefers config 0 (energy ordering drives an accepted swap).
  let energies : Array[Array[Double]] = [[5.0, 1.0], [1.0, 5.0]]
  let betas = [1.0, 1.0]
  let swaps = @src.omm_hremd_swap_cycle(energies, betas, uniforms=[0.5])
  assert_eq(swaps.length(), 1)
  assert_true(swaps[0].accepted)
  let occ = @src.omm_hremd_final_occupancy(2, swaps)
  // Config 0 moved to replica 1, config 1 to replica 0.
  assert_eq(occ[0], 1)
  assert_eq(occ[1], 0)
  // With an unfavorable energy landscape (diagonal preferred) the swap
  // is rejected at Delta = +8.
  let energies_r : Array[Array[Double]] = [[1.0, 5.0], [5.0, 1.0]]
  let swaps2 = @src.omm_hremd_swap_cycle(energies_r, betas, uniforms=[0.5])
  let occ2 = @src.omm_hremd_final_occupancy(2, swaps2)
  assert_eq(occ2[0], 0)
  assert_eq(occ2[1], 1)
}

///|
test "hremd_three_replicas_even_pairs" {
  // Three replicas: only pairs (0,1) attempted in an even cycle.
  let energies : Array[Array[Double]] = [
    [1.0, 2.0, 3.0],
    [1.0, 2.0, 3.0],
    [1.0, 2.0, 3.0],
  ]
  let betas = [1.0, 1.0, 1.0]
  let swaps = @src.omm_hremd_swap_cycle(energies, betas, uniforms=[0.5, 0.5])
  assert_eq(swaps.length(), 1)
  assert_eq(swaps[0].i, 0)
  assert_eq(swaps[0].j, 1)
  // Identical rows -> swap always accepted at Delta = 0.
  assert_true(swaps[0].accepted)
}
