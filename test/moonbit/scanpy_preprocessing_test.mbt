///|
/// Tests for scanpy preprocessing core (normalize_total / log1p / scale /
/// highly_variable_genes).

///|
fn sp_close(a : Double, b : Double) -> Bool {
  let d = a - b
  let ad = if d < 0.0 { 0.0 - d } else { d }
  ad < 1.0e-9
}

///|
test "scanpy_normalize_total" {
  let x = [[1.0, 2.0, 3.0], [4.0, 0.0, 0.0]]
  let out = @src.scanpy_normalize_total(x, 10.0)
  // cell0 sum=6 -> *10/6 ; cell1 sum=4 -> *10/4
  assert_eq(out[0][2], 5.0)
  assert_eq(out[1][0], 10.0)
  // each row sums to target_sum
  let mut s0 = 0.0
  for v in out[0] {
    s0 = s0 + v
  }
  let mut s1 = 0.0
  for v in out[1] {
    s1 = s1 + v
  }
  assert_eq(sp_close(s0, 10.0), true)
  assert_eq(sp_close(s1, 10.0), true)
}

///|
test "scanpy_log1p" {
  let out = @src.scanpy_log1p([[1.0, 3.0]])
  assert_eq(sp_close(out[0][0], 0.6931471805599453), true)
  assert_eq(sp_close(out[0][1], 1.3862943611198906), true)
}

///|
test "scanpy_scale" {
  let x = [[0.0, 1.0], [2.0, 1.0], [4.0, 1.0]]
  let out = @src.scanpy_scale(x, None)
  // gene0 [0,2,4] mean=2 std=2 -> [-1,0,1]; gene1 constant -> [0,0,0]
  assert_eq(out[0][0], -1.0)
  assert_eq(out[1][0], 0.0)
  assert_eq(out[2][0], 1.0)
  assert_eq(out[0][1], 0.0)
  assert_eq(out[2][1], 0.0)
}

///|
test "scanpy_scale_clip" {
  let x = [[0.0], [4.0]]
  let out = @src.scanpy_scale(x, Some(0.5))
  // unscaled would be [-1, 1]; clipped to [-0.5, 0.5]
  assert_eq(sp_close(out[0][0], -0.5), true)
  assert_eq(sp_close(out[1][0], 0.5), true)
}

///|
test "scanpy_highly_variable_genes" {
  // 4 cells × 4 genes; gene0 constant (zero dispersion), gene1/3 variable,
  // gene2 high mean (log1p mean > max_mean); gene1/3 binned together, gene3 wins
  // the lower-dispersion tie so only gene1 is flagged
  let x = [
    [5.0, 0.0, 10.0, 1.0],
    [5.0, 2.0, 20.0, 2.0],
    [5.0, 4.0, 30.0, 3.0],
    [5.0, 6.0, 40.0, 4.0],
  ]
  let hv = @src.scanpy_highly_variable_genes(x, 4, 0.0125, 3.0, 0.5)
  assert_eq(hv.highly_variable.length(), 4)
  assert_eq(hv.means.length(), 4)
  assert_eq(hv.highly_variable[0], false)
  assert_eq(hv.highly_variable[1], true)
  assert_eq(hv.highly_variable[2], false)
  assert_eq(hv.highly_variable[3], false)
}