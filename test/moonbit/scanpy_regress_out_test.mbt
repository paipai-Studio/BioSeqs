///|
/// Tests for scanpy `pp.regress_out` (linear-regression confounder removal).

///|
test "scanpy_regress_out" {
  let x = [[1.0, 5.0], [2.0, 7.0], [3.0, 6.0]]
  let cov = [[1.0], [2.0], [3.0]]
  let resid = @src.scanpy_regress_out(x, cov)
  // var0 == covariate perfectly -> residual 0
  for i = 0; i < 3; i = i + 1 {
    assert_true(assert_close(resid[i][0], 0.0, 1.0e-6))
  }
  // var1 residual: zero mean and uncorrelated with the covariate
  let mut mean = 0.0
  let mut cov_sum = 0.0
  for i = 0; i < 3; i = i + 1 {
    mean = mean + resid[i][1]
    cov_sum = cov_sum + resid[i][1] * cov[i][0]
  }
  assert_true(assert_close(mean / 3.0, 0.0, 1.0e-6))
  assert_true(assert_close(cov_sum, 0.0, 1.0e-6))
}