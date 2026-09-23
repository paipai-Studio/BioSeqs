///|
/// Tests for scikit-allel site frequency spectra and Patterson Fst.

///|
test "skallel_sfs" {
  let dac = [0, 1, 1, 2, 2, 2]
  assert_eq(@src.skallel_sfs(dac, Some(3)), [1, 2, 3, 0])
  // scaled: s[k] * k
  assert_eq(@src.skallel_sfs_scaled(dac, Some(3)), [0.0, 2.0, 6.0, 0.0])
}

///|
test "skallel_sfs_folded" {
  let ac = [[0, 4], [1, 3], [2, 2], [4, 0], [0, 4]]
  // minor allele counts: 0,1,2,0,0  -> counts [3,1,1]
  assert_eq(@src.skallel_sfs_folded(ac, Some(4)), [3, 1, 1])
  // folded scaled: s[k]*k*(n-k)/n, n=4
  let sc = @src.skallel_sfs_folded_scaled(ac, Some(4))
  assert_true(assert_close(sc[0], 0.0, 1.0e-9))
  assert_true(assert_close(sc[1], 0.75, 1.0e-9))
  assert_true(assert_close(sc[2], 1.0, 1.0e-9))
}

///|
test "skallel_joint_sfs" {
  let d1 = [0, 1, 2]
  let d2 = [1, 1, 2]
  let j = @src.skallel_joint_sfs(d1, d2, Some(2), Some(2))
  assert_eq(j.length(), 3)
  assert_eq(j[0], [0, 1, 0])
  assert_eq(j[1], [0, 1, 0])
  assert_eq(j[2], [0, 0, 1])
}

///|
test "skallel_joint_sfs_folded" {
  let ac1 = [[0, 4], [2, 2], [4, 0]] // mac 0,2,0
  let ac2 = [[2, 2], [0, 4], [1, 3]] // mac 2,0,1
  let j = @src.skallel_joint_sfs_folded(ac1, ac2, Some(4), Some(4))
  // dims (4//2+1) x (4//2+1) = 3x3; entries at (mac1,mac2): (0,2),(2,0),(0,1)
  assert_eq(j.length(), 3)
  assert_eq(j[0], [0, 1, 1])  // (0,1)=1 and (0,2)=1
  assert_eq(j[1], [0, 0, 0])
  assert_eq(j[2], [1, 0, 0])  // (2,0)=1
}

///|
test "skallel_fold_sfs" {
  assert_eq(@src.skallel_fold_sfs([1, 2, 3, 0], 3), [1, 5])
}

///|
test "skallel_patterson_fst" {
  let aca = [[0, 6], [6, 0]]
  let acb = [[6, 0], [0, 6]]
  let (num, den) = @src.skallel_patterson_fst(aca, acb)
  assert_eq(num.length(), 2)
  assert_true(assert_close(num[0], 1.0, 1.0e-9))
  assert_true(assert_close(num[1], 1.0, 1.0e-9))
  assert_true(assert_close(den[0], 1.0, 1.0e-9))
  assert_true(assert_close(den[1], 1.0, 1.0e-9))
}