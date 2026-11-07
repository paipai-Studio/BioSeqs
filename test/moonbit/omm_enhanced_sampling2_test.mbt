///|
test "omm_steered_md_energy_and_force" {
  // center(t) = center0 + v * t ; F = -k * (cv - center)
  // at t = 2 with v = 0.5 the center sits at 1.0
  let f0 = @src.omm_steered_md_force(1.0, 0.0, 0.5, 2.0, 10.0)
  assert_true(@moonbit.assert_num_close(f0, 0.0, 1.0e-12))
  let f1 = @src.omm_steered_md_force(1.5, 0.0, 0.5, 2.0, 10.0)
  assert_true(@moonbit.assert_num_close(f1, -5.0, 1.0e-12))
  let f2 = @src.omm_steered_md_force(0.5, 0.0, 0.5, 2.0, 10.0)
  assert_true(@moonbit.assert_num_close(f2, 5.0, 1.0e-12))
  // E = 0.5 * k * (cv - center)^2 = 0.5 * 10 * 0.25
  let e1 = @src.omm_steered_md_energy(1.5, 0.0, 0.5, 2.0, 10.0)
  assert_true(@moonbit.assert_num_close(e1, 1.25, 1.0e-12))
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_steered_md_energy(1.0, 0.0, 0.5, 2.0, 10.0),
      0.0,
      1.0e-12,
    ),
  )
}

///|
test "omm_steered_md_work_accumulation" {
  // work = sum of f_prev * dcv; cv moves 1.0 -> 1.5 while center moves 0 -> 1
  // f_prev = -2 * (1.0 - 0.0) = -2, dcv = 0.5 -> work = -1
  let work = @src.omm_steered_md_work([1.0, 1.5], [0.0, 1.0], 0.0, 1.0, 2.0)
  assert_true(@moonbit.assert_num_close(work, -1.0, 1.0e-12))
  // two segments: f at t=0 is -2 (dcv 0.5), f at t=1 is -2*(1.5-1.0)=-1 (dcv 0.5)
  let work2 = @src.omm_steered_md_work(
    [1.0, 1.5, 2.0],
    [0.0, 1.0, 2.0],
    0.0,
    1.0,
    2.0,
  )
  assert_true(@moonbit.assert_num_close(work2, -1.5, 1.0e-12))
  // fewer than two points yields no work
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_steered_md_work([1.0], [0.0], 0.0, 1.0, 2.0),
      0.0,
      1.0e-12,
    ),
  )
}

///|
test "omm_rest2_scale_factors_and_energy" {
  // ratio = 300/600 = 0.5: solute-solute scales as ratio^2, solute-solvent as ratio
  let (ss, sv) = @src.omm_rest2_scale_factors(300.0, 600.0)
  assert_true(@moonbit.assert_num_close(ss, 0.25, 1.0e-12))
  assert_true(@moonbit.assert_num_close(sv, 0.5, 1.0e-12))
  let (ss1, sv1) = @src.omm_rest2_scale_factors(300.0, 300.0)
  assert_true(@moonbit.assert_num_close(ss1, 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(sv1, 1.0, 1.0e-12))
  // E' = ss*E_ss + sv*E_sv + E_vv
  let e = @src.omm_rest2_scaled_energy(4.0, 2.0, 1.0, 300.0, 600.0)
  assert_true(
    @moonbit.assert_num_close(e, 0.25 * 4.0 + 0.5 * 2.0 + 1.0, 1.0e-12),
  )
}

///|
test "omm_rest2_swap_accept_direction" {
  // replica 1 (cold, T=300) carries solute-solute energy 0 and replica 2
  // (hot, T=600) carries 100: swapping lowers the scaled total
  // (delta = 100 - 25 = 75 > 0 with beta_1 > beta_2) -> accepted
  let accepted = @src.omm_rest2_swap_accept(
    0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 300.0, 600.0,
  )
  assert_true(accepted)
  // mirrored energies reverse the sign of delta while betas keep their order
  let rejected = @src.omm_rest2_swap_accept(
    100.0, 0.0, 0.0, 0.0, 0.0, 0.0, 300.0, 600.0,
  )
  assert_false(rejected)
  // identical configurations give delta = 0 -> arg not positive -> rejected
  let neutral = @src.omm_rest2_swap_accept(
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 300.0, 600.0,
  )
  assert_false(neutral)
}
