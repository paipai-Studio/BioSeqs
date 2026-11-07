///|
/// Tests for hbond_frequency.mbt (HBond frequency across models).

///|
fn hb_freq_sample_obs(
  d_res : Int32,
  a_res : Int32,
  d_atom : String,
  a_atom : String,
  dist : Double,
  ang : Double,
) -> @src.HBondObservation {
  @src.HBondObservation::new(
    donor_residue=d_res,
    acceptor_residue=a_res,
    donor_atom=d_atom,
    acceptor_atom=a_atom,
    distance=dist,
    angle=ang,
  )
}

///|
test "hbond_frequency_counts_across_models" {
  // HBond A (5N -> 10O) appears in all 3 models; HBond B (7NE2 -> 12OD1)
  // only in the second model.
  let obs_a0 = hb_freq_sample_obs(5, 10, "N", "O", 3.0, 160.0)
  let obs_a1 = hb_freq_sample_obs(5, 10, "N", "O", 3.2, 150.0)
  let obs_b1 = hb_freq_sample_obs(7, 12, "NE2", "OD1", 2.8, 170.0)
  let obs_a2 = hb_freq_sample_obs(5, 10, "N", "O", 3.1, 155.0)
  let result = @src.hbond_frequency(
    [[obs_a0], [obs_a1, obs_b1], [obs_a2]],
    total_models=3,
  )
  assert_eq(result.length(), 2)
  // First-seen bond first: A with count 3, frequency 1.0.
  assert_eq(@src.HBondFrequency::get_donor_residue(result[0]), 5)
  assert_eq(@src.HBondFrequency::get_acceptor_residue(result[0]), 10)
  assert_eq(result[0].donor_atom, "N")
  assert_eq(result[0].acceptor_atom, "O")
  assert_eq(result[0].count, 3)
  assert_eq(result[0].total_models, 3)
  assert_true(
    @moonbit.assert_num_close(
      @src.HBondFrequency::get_frequency(result[0]),
      1.0,
      1.0e-12,
    ),
  )
  assert_true(@moonbit.assert_num_close(result[0].avg_distance, 3.1, 1.0e-9))
  assert_true(@moonbit.assert_num_close(result[0].avg_angle, 155.0, 1.0e-9))
  // Second bond: count 1, frequency 1/3.
  assert_eq(result[1].count, 1)
  assert_true(
    @moonbit.assert_num_close(
      @src.HBondFrequency::get_frequency(result[1]),
      1.0 / 3.0,
      1.0e-12,
    ),
  )
  assert_true(@moonbit.assert_num_close(result[1].avg_distance, 2.8, 1.0e-9))
  assert_true(@moonbit.assert_num_close(result[1].avg_angle, 170.0, 1.0e-9))
}

///|
test "hbond_frequency_empty" {
  assert_eq(@src.hbond_frequency([], total_models=5).length(), 0)
  // Models present but no observations -> no bonds.
  assert_eq(@src.hbond_frequency([[], []], total_models=2).length(), 0)
}

///|
test "filter_hbonds_by_frequency_threshold" {
  let obs_a0 = hb_freq_sample_obs(5, 10, "N", "O", 3.0, 160.0)
  let obs_a1 = hb_freq_sample_obs(5, 10, "N", "O", 3.2, 150.0)
  let obs_b1 = hb_freq_sample_obs(7, 12, "NE2", "OD1", 2.8, 170.0)
  let result = @src.hbond_frequency(
    [[obs_a0], [obs_a1, obs_b1]],
    total_models=2,
  )
  // Frequencies: A = 1.0, B = 0.5.
  let strong = @src.filter_hbonds_by_frequency(result, 0.75)
  assert_eq(strong.length(), 1)
  assert_eq(strong[0].count, 2)
  // Threshold 0 keeps both.
  assert_eq(@src.filter_hbonds_by_frequency(result, 0.0).length(), 2)
  // Threshold above 1 keeps none.
  assert_eq(@src.filter_hbonds_by_frequency(result, 1.5).length(), 0)
}

///|
test "sort_hbonds_by_frequency_descending" {
  let obs_a0 = hb_freq_sample_obs(5, 10, "N", "O", 3.0, 160.0)
  let obs_a1 = hb_freq_sample_obs(5, 10, "N", "O", 3.2, 150.0)
  let obs_b1 = hb_freq_sample_obs(7, 12, "NE2", "OD1", 2.8, 170.0)
  let result = @src.hbond_frequency(
    [[obs_a0], [obs_a1, obs_b1]],
    total_models=2,
  )
  let sorted = @src.sort_hbonds_by_frequency(result)
  assert_eq(sorted.length(), 2)
  // Descending: A (1.0) before B (0.5).
  assert_true(
    @moonbit.assert_num_close(
      @src.HBondFrequency::get_frequency(sorted[0]),
      1.0,
      1.0e-12,
    ),
  )
  assert_true(
    @moonbit.assert_num_close(
      @src.HBondFrequency::get_frequency(sorted[1]),
      0.5,
      1.0e-12,
    ),
  )
  // Sorting an empty list stays empty.
  assert_eq(@src.sort_hbonds_by_frequency([]).length(), 0)
}
