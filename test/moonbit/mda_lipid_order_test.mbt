///|
/// Tests for mda_lipid_order.mbt (deuterium order parameter SCD).

///|
test "mda_lipid_scd_perpendicular_bond" {
  // C-H perpendicular to the z normal: cos^2 = 0 -> SCD = -0.5*(0-1) = 0.5.
  let scd = @src.mda_lipid_scd([[1.0, 0.0, 0.0]])
  assert_true(@moonbit.assert_num_close(scd, 0.5, 1.0e-12))
}

///|
test "mda_lipid_scd_parallel_bond" {
  // C-H along the z normal: cos^2 = 1 -> SCD = -0.5*(3-1) = -1.
  let scd = @src.mda_lipid_scd([[0.0, 0.0, 2.0]])
  assert_true(@moonbit.assert_num_close(scd, -1.0, 1.0e-12))
}

///|
test "mda_lipid_scd_magic_angle" {
  // Vector (1,1,1): cos^2 = 1/3 -> 3*(1/3)-1 = 0 -> SCD = 0.
  let scd = @src.mda_lipid_scd([[1.0, 1.0, 1.0]])
  assert_true(@moonbit.assert_num_close(scd, 0.0, 1.0e-12))
}

///|
test "mda_lipid_scd_averages_over_frames" {
  // Frame 0 perpendicular (contribution 3*0-1 = -1), frame 1 parallel
  // (contribution 2): mean = 0.5 -> SCD = -0.25.
  let scd = @src.mda_lipid_scd([[1.0, 0.0, 0.0], [0.0, 0.0, 3.0]])
  assert_true(@moonbit.assert_num_close(scd, -0.25, 1.0e-12))
}

///|
test "mda_lipid_scd_zero_vector_counts_in_denominator" {
  // Zero-length vector is skipped in the sum but still counted in n.
  // sum = -1 (from the perpendicular frame), n = 2 -> SCD = 0.25.
  let scd = @src.mda_lipid_scd([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
  assert_true(@moonbit.assert_num_close(scd, 0.25, 1.0e-12))
}

///|
test "mda_lipid_scd_empty_and_degenerate_normal" {
  assert_true(@moonbit.assert_num_close(@src.mda_lipid_scd([]), 0.0, 1.0e-12))
  // Zero-length normal -> 0.
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_lipid_scd([[1.0, 0.0, 0.0]], normal=[0.0, 0.0, 0.0]),
      0.0,
      1.0e-12,
    ),
  )
}

///|
test "mda_lipid_scd_custom_normal" {
  // Normal along x, bond along x: fully ordered -> SCD = -1.
  let scd = @src.mda_lipid_scd([[3.0, 0.0, 0.0]], normal=[2.0, 0.0, 0.0])
  assert_true(@moonbit.assert_num_close(scd, -1.0, 1.0e-12))
}

///|
test "mda_lipid_order_profile_and_average" {
  // Carbon 0: perpendicular bond (SCD = 0.5); carbon 1: parallel (SCD = -1).
  let profile = @src.mda_lipid_order_profile([
    [[1.0, 0.0, 0.0]],
    [[0.0, 0.0, 1.0]],
  ])
  assert_eq(profile.length(), 2)
  assert_true(@moonbit.assert_num_close(profile[0], 0.5, 1.0e-12))
  assert_true(@moonbit.assert_num_close(profile[1], -1.0, 1.0e-12))
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_lipid_average_order(profile),
      -0.25,
      1.0e-12,
    ),
  )
}

///|
test "mda_lipid_average_order_empty" {
  assert_true(
    @moonbit.assert_num_close(@src.mda_lipid_average_order([]), 0.0, 1.0e-12),
  )
}

///|
test "mda_lipid_scd_from_positions" {
  // Carbon at origin, hydrogen at (1,0,0): perpendicular -> SCD = 0.5.
  let scd = @src.mda_lipid_scd_from_positions([[0.0, 0.0, 0.0]], [
    [1.0, 0.0, 0.0],
  ])
  assert_true(@moonbit.assert_num_close(scd, 0.5, 1.0e-12))
}

///|
test "mda_lipid_order_state_classification" {
  assert_eq(@src.mda_lipid_order_state(0.5), "Ordered")
  assert_eq(@src.mda_lipid_order_state(-0.7), "Ordered")
  assert_eq(@src.mda_lipid_order_state(0.4), "Intermediate")
  assert_eq(@src.mda_lipid_order_state(-0.3), "Intermediate")
  assert_eq(@src.mda_lipid_order_state(0.2), "Intermediate")
  assert_eq(@src.mda_lipid_order_state(0.1), "Disordered")
  assert_eq(@src.mda_lipid_order_state(-0.19), "Disordered")
}
