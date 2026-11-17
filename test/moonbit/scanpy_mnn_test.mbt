///|
/// Tests for MNN batch correction.

///|
test "scanpy_mnn_correct" {
  // reference with two clusters; batch shifted +5 on gene 0 near both clusters
  let ref_mat = [
    [10.0, 0.0],
    [10.0, 1.0],
    [10.0, 2.0],
    [0.0, 10.0],
    [0.0, 11.0],
  ]
  let bat_mat = [[15.0, 0.0], [15.0, 1.0], [5.0, 10.0], [5.0, 11.0]]
  let res = @src.scanpy_mnn_correct([ref_mat, bat_mat], k=1, sigma=1.0)
  assert_eq(res.corrected.length(), 2)
  // four mutual pairs: (r0,b0),(r1,b1),(r3,b2),(r4,b3)
  assert_eq(res.n_mnn[1], 4)
  // reference unchanged
  assert_eq(res.corrected[0][0][0], 10.0)
  // b0 pulled back to (10, 0)
  assert_true(@moonbit.assert_num_close(res.corrected[1][0][0], 10.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(res.corrected[1][0][1], 0.0, 1.0e-9))
  // b2 pulled back to (0, 10)
  assert_true(@moonbit.assert_num_close(res.corrected[1][2][0], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(res.corrected[1][2][1], 10.0, 1.0e-9))
}

///|
test "scanpy_mnn_correct_adatas" {
  let a1 = @src.AnnData::new([[10.0, 0.0], [10.0, 1.0]], ["r0", "r1"], [
    "G0", "G1",
  ])
  let a2 = @src.AnnData::new([[15.0, 0.0], [15.0, 1.0]], ["b0", "b1"], [
    "G0", "G1",
  ])
  let _ = a1.set_obs("batch", ["ref", "ref"])
  let _ = a2.set_obs("batch", ["b", "b"])
  let res = @src.scanpy_mnn_correct_adatas([a1, a2], k=1, sigma=1.0)
  assert_eq(res.corrected.length(), 2)
  // batch cells corrected in place toward the reference
  assert_true(@moonbit.assert_num_close(a2.x[0][0], 10.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(a2.x[1][0], 10.0, 1.0e-9))
  // reference untouched
  assert_eq(a1.x[0][0], 10.0)
  // obs annotation preserved
  match a2.obs_column("batch") {
    Some(v) => assert_eq(v[0], "b")
    None => assert_true(false)
  }
}
