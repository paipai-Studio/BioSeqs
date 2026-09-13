///|
/// Tests for affyPLM — probe-level model and QC.

///|
test "affy_plm_sample_probeset construction" {
  let ps = @src.affy_plm_sample_probeset()
  assert_eq(ps.name, "ps1")
  assert_eq(ps.probes.length(), 4)
  assert_eq(ps.probes[0].length(), 3)
}

///|
test "affy_plm_fit_plm basic" {
  let ps = @src.affy_plm_sample_probeset()
  let fit = @src.fit_plm(ps)
  assert_eq(fit.probeset_name, "ps1")
  assert_eq(fit.chip_effects.length(), 3)
  assert_eq(fit.probe_effects.length(), 4)
  assert_true(fit.n_iterations > 0)
  assert_true(fit.sigma > 0.0)
}

///|
test "affy_plm_fit_plm residuals dimensions" {
  let ps = @src.affy_plm_sample_probeset()
  let fit = @src.fit_plm(ps)
  assert_eq(fit.residuals.length(), 4)
  for row in fit.residuals {
    assert_eq(row.length(), 3)
  }
}

///|
test "affy_plm_fit_plm weights dimensions" {
  let ps = @src.affy_plm_sample_probeset()
  let fit = @src.fit_plm(ps)
  assert_eq(fit.weights.length(), 4)
  for row in fit.weights {
    assert_eq(row.length(), 3)
    for w in row {
      assert_true(w >= 0.0)
      assert_true(w <= 1.0)
    }
  }
}

///|
test "affy_plm_fit_plm probe effects sum to zero" {
  let ps = @src.affy_plm_sample_probeset()
  let fit = @src.fit_plm(ps)
  let mut sum = 0.0
  for pe in fit.probe_effects {
    sum += pe
  }
  assert_true(sum.abs() < 1.0e-6)
}

///|
test "affy_plm_nuse basic" {
  let ps = @src.affy_plm_sample_probeset()
  let fit = @src.fit_plm(ps)
  let nuse_vals = @src.nuse([fit])
  assert_eq(nuse_vals.length(), 3)
  for v in nuse_vals {
    assert_true(v >= 0.0)
  }
}

///|
test "affy_plm_rle basic" {
  let ps = @src.affy_plm_sample_probeset()
  let fit = @src.fit_plm(ps)
  let rle_vals = @src.rle([fit])
  assert_eq(rle_vals.length(), 3)
}

///|
test "affy_plm_qc report" {
  let ps = @src.affy_plm_sample_probeset()
  let fit = @src.fit_plm(ps)
  let report = @src.affy_plm_qc([fit])
  assert_eq(report.n_probesets, 1)
  assert_eq(report.nuse.length(), 3)
  assert_eq(report.rle.length(), 3)
}

///|
test "affy_plm_fit_plm huber estimator" {
  let ps = @src.affy_plm_sample_probeset()
  let fit = @src.fit_plm(ps, estimator=@src.PlmEstimator::Huber(1.345))
  assert_eq(fit.chip_effects.length(), 3)
}

///|
test "affy_plm_fit_plm tukey estimator" {
  let ps = @src.affy_plm_sample_probeset()
  let fit = @src.fit_plm(ps, estimator=@src.PlmEstimator::Tukey(4.685))
  assert_eq(fit.chip_effects.length(), 3)
}

///|
test "affy_plm_fit_plm empty probeset" {
  let ps = @src.PlmProbeset::new(name="empty", probes=[])
  let fit = @src.fit_plm(ps)
  assert_eq(fit.chip_effects.length(), 0)
  assert_eq(fit.probe_effects.length(), 0)
}

///|
test "affy_plm_fit_plm converges with more iterations" {
  let ps = @src.affy_plm_sample_probeset()
  let fit = @src.fit_plm(ps, max_iter=50, tol=1.0e-8)
  assert_true(fit.n_iterations <= 50)
}