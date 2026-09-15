///|
/// Tests for ashr module (Adaptive Shrinkage).

test "ashr_erf" {
  let e0 = @src.ashr_erf(0.0)
  assert_true(e0.abs() < 1.0e-9)
  let e1 = @src.ashr_erf(1.0)
  assert_true((e1 - 0.8427).abs() < 0.001)
  let eneg = @src.ashr_erf(-1.0)
  assert_true((eneg + e1).abs() < 1.0e-9)
}

test "ashr_pnorm" {
  let p0 = @src.ashr_pnorm(0.0)
  assert_true((p0 - 0.5).abs() < 1.0e-9)
  let p_pos = @src.ashr_pnorm(1.96)
  assert_true(p_pos > 0.97 && p_pos < 0.98)
  let p_neg = @src.ashr_pnorm(-1.96)
  assert_true(p_neg < 0.03 && p_neg > 0.02)
  assert_true((p_pos + p_neg - 1.0).abs() < 1.0e-9)
}

test "ashr_dnorm" {
  let d0 = @src.ashr_dnorm(0.0, 0.0, 1.0)
  assert_true((d0 - 0.3989423).abs() < 1.0e-6)
  let d_far = @src.ashr_dnorm(10.0, 0.0, 1.0)
  assert_true(d_far < 1.0e-20)
  let d_shift = @src.ashr_dnorm(2.0, 2.0, 1.0)
  assert_true((d_shift - d0).abs() < 1.0e-9)
}

test "ashr_default_grid" {
  let grid = @src.ashr_default_grid()
  assert_true(grid.length() == 7)
  assert_eq(grid[0], 0.0)
  assert_true(grid[1] > 0.0)
}

test "ashr_shrink_basic" {
  let estimates = [2.0, -1.8, 0.3, 0.05]
  let ses = [0.2, 0.2, 0.2, 0.18]
  let fit = @src.ashr_shrink(estimates, ses)
  assert_eq(fit.n, 4)
  assert_eq(fit.posterior_mean.length(), 4)
  assert_eq(fit.posterior_sd.length(), 4)
  assert_eq(fit.lfsr.length(), 4)
  assert_true(fit.iterations > 0)
  assert_true(fit.log_likelihood > -1.0e300)
  let mut wsum = 0.0
  for wv in fit.prior_weights {
    wsum = wsum + wv
  }
  assert_true((wsum - 1.0).abs() < 1.0e-6)
}

test "ashr_shrink_z" {
  let z = [2.5, -2.0, 0.5, 0.1]
  let fit = @src.ashr_shrink_z(z)
  assert_eq(fit.n, 4)
  assert_true(fit.log_likelihood.abs() < 1.0e12)
}

test "ashr_count_significant" {
  let estimates = [3.0, -2.5, 0.2, 0.1]
  let ses = [0.3, 0.3, 0.3, 0.3]
  let fit = @src.ashr_shrink(estimates, ses)
  let n_sig = @src.ashr_count_significant(fit, 0.05)
  assert_true(n_sig >= 0 && n_sig <= 4)
  let n_all = @src.ashr_count_significant(fit, 1.01)
  assert_eq(n_all, 4)
}

test "ashr_shrinkage_factor" {
  let estimates = [2.0, -1.5, 0.3]
  let ses = [0.2, 0.2, 0.2]
  let fit = @src.ashr_shrink(estimates, ses)
  let factors = @src.ashr_shrinkage_factor(fit, estimates)
  assert_eq(factors.length(), 3)
  for f in factors {
    assert_true(f >= 0.0 && f <= 1.0)
  }
}

test "ashr_mode_component" {
  let estimates = [2.0, -1.5, 0.3]
  let ses = [0.2, 0.2, 0.2]
  let fit = @src.ashr_shrink(estimates, ses)
  let mode = @src.ashr_mode_component(fit)
  assert_true(mode >= 0 && mode < fit.prior_grid.length())
}

test "ashr_summary" {
  let estimates = [2.0, -1.5]
  let ses = [0.2, 0.2]
  let fit = @src.ashr_shrink(estimates, ses)
  let s = @src.ashr_summary(fit)
  assert_true(s.contains("AshrFit"))
  assert_true(s.contains("lfsr"))
}

test "ashr_example_data" {
  let (est, ses) = @src.ashr_example_data()
  assert_eq(est.length(), 8)
  assert_eq(ses.length(), 8)
  let fit = @src.ashr_shrink(est, ses)
  assert_eq(fit.n, 8)
}

test "ashr_posterior_shrinks_toward_zero" {
  let estimates = [0.5, 0.3, 0.1]
  let ses = [0.5, 0.5, 0.5]
  let fit = @src.ashr_shrink(estimates, ses)
  let mut i = 0
  while i < fit.n {
    assert_true(fit.posterior_mean[i].abs() <= estimates[i].abs() + 1.0e-9)
    i = i + 1
  }
}

test "ashr_lfsr_extremes" {
  let estimates = [5.0, 0.0]
  let ses = [0.1, 0.5]
  let fit = @src.ashr_shrink(estimates, ses)
  assert_true(fit.lfsr[0] < fit.lfsr[1])
}

test "ashr_posterior_sd_positive" {
  let estimates = [2.0, -1.0, 0.5]
  let ses = [0.2, 0.3, 0.4]
  let fit = @src.ashr_shrink(estimates, ses)
  let mut i = 0
  while i < fit.n {
    assert_true(fit.posterior_sd[i] >= 0.0)
    assert_true(fit.lfsr[i] >= 0.0 && fit.lfsr[i] <= 1.0)
    i = i + 1
  }
}
