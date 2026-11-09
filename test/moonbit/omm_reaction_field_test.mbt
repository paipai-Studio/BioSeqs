///|
test "omm reaction field energy analytic value" {
  // U(r) = qi*qj*(1/r + k_rf*r^2 - c_rf), cutoff=10, eps=78.5, r=1, q=+1/+1
  let cutoff = 10.0
  let eps = 78.5
  let k_rf = 1.0 / (cutoff * cutoff * cutoff) * (eps - 1.0) / (2.0 * eps + 1.0)
  let c_rf = 1.0 / cutoff * (3.0 * eps) / (2.0 * eps + 1.0)
  let u = @src.omm_reaction_field_energy(1.0, 1.0, 1.0, cutoff)
  assert_true(@moonbit.assert_num_close(u, 1.0 + k_rf - c_rf, 1.0e-12))
}

///|
test "omm reaction field energy zero beyond cutoff and at origin" {
  assert_eq(@src.omm_reaction_field_energy(10.0, 1.0, 1.0, 10.0), 0.0)
  assert_eq(@src.omm_reaction_field_energy(11.0, 1.0, 1.0, 10.0), 0.0)
  assert_eq(@src.omm_reaction_field_energy(0.0, 1.0, 1.0, 10.0), 0.0)
}

///|
test "omm reaction field energy opposite charges attractive" {
  let u = @src.omm_reaction_field_energy(1.0, 1.0, -1.0, 10.0)
  assert_true(u < 0.0)
}

///|
test "omm reaction field force analytic value" {
  // F(r) = qi*qj*(1/r^2 - 2*k_rf*r)
  let cutoff = 10.0
  let eps = 78.5
  let k_rf = 1.0 / (cutoff * cutoff * cutoff) * (eps - 1.0) / (2.0 * eps + 1.0)
  let f = @src.omm_reaction_field_force(1.0, 1.0, 1.0, cutoff)
  assert_true(@moonbit.assert_num_close(f, 1.0 - 2.0 * k_rf, 1.0e-12))
  assert_true(f > 0.0)
}

///|
test "omm reaction field force zero beyond cutoff" {
  assert_eq(@src.omm_reaction_field_force(10.0, 1.0, 1.0, 10.0), 0.0)
}

///|
test "omm reaction field total open boundary pair" {
  // two atoms at r=2, q=+1/-1, cutoff=10, box=0
  let cutoff = 10.0
  let eps = 78.5
  let k_rf = 1.0 / (cutoff * cutoff * cutoff) * (eps - 1.0) / (2.0 * eps + 1.0)
  let c_rf = 1.0 / cutoff * (3.0 * eps) / (2.0 * eps + 1.0)
  let e = @src.omm_reaction_field_total(
    [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
    [1.0, -1.0],
    cutoff,
  )
  let expected = -1.0 * (0.5 + k_rf * 4.0 - c_rf)
  assert_true(@moonbit.assert_num_close(e, expected, 1.0e-12))
}

///|
test "omm reaction field total applies minimum image" {
  // separation 9 in box 10 wraps to 1: same charges pair at r=1
  let e_wrapped = @src.omm_reaction_field_total(
    [[0.0, 0.0, 0.0], [9.0, 0.0, 0.0]],
    [1.0, 1.0],
    10.0,
    box=[10.0, 10.0, 10.0],
  )
  let e_direct = @src.omm_reaction_field_total(
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    [1.0, 1.0],
    10.0,
  )
  assert_true(@moonbit.assert_num_close(e_wrapped, e_direct, 1.0e-12))
}

///|
test "omm reaction field total zero charges" {
  let e = @src.omm_reaction_field_total(
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    [0.0, 0.0],
    10.0,
  )
  assert_eq(e, 0.0)
}

///|
test "omm reaction field shifted energy vanishes at cutoff" {
  // just below the cutoff the shifted potential tends to zero
  let u = @src.omm_reaction_field_shifted_energy(9.999999999, 1.0, 1.0, 10.0)
  assert_true(u.abs() < 1.0e-6)
}

///|
test "omm reaction field shifted energy equals unshifted minus cutoff value" {
  let cutoff = 10.0
  let eps = 78.5
  let k_rf = 1.0 / (cutoff * cutoff * cutoff) * (eps - 1.0) / (2.0 * eps + 1.0)
  let c_rf = 1.0 / cutoff * (3.0 * eps) / (2.0 * eps + 1.0)
  let u_plain = @src.omm_reaction_field_energy(1.0, 1.0, 1.0, cutoff)
  let u_shifted = @src.omm_reaction_field_shifted_energy(1.0, 1.0, 1.0, cutoff)
  let u_cut = 1.0 / cutoff + k_rf * cutoff * cutoff - c_rf
  assert_true(@moonbit.assert_num_close(u_shifted, u_plain - u_cut, 1.0e-12))
}

///|
test "omm reaction field vs pme error nonnegative" {
  let err1 = @src.omm_rf_vs_pme_error(1.0, 1.0, 1.0, 10.0)
  let err2 = @src.omm_rf_vs_pme_error(2.0, 1.0, -1.0, 10.0)
  assert_true(err1 >= 0.0)
  assert_true(err2 >= 0.0)
  assert_true(err1.is_finite())
  assert_true(err2.is_finite())
}
