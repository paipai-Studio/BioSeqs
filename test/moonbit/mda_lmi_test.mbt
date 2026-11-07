///|
/// Tests for linear mutual information (mda_lmi.mbt).
test "mda_lmi_independent" {
  // 二元独立变量：联合分布 = 边缘之积 → I(X;Y) = 0
  let x = [0.0, 0.0, 1.0, 1.0]
  let y = [0.0, 1.0, 0.0, 1.0]
  let mi = @src.mda_lmi(x, y)
  assert_true(@moonbit.assert_num_close(mi, 0.0, 1.0e-9))
}

///|
test "mda_lmi_perfect_correlation" {
  // y = x 完全相关 → I = H(X) = ln 2
  let x = [0.0, 0.0, 1.0, 1.0]
  let mi = @src.mda_lmi(x, x)
  assert_true(@moonbit.assert_num_close(mi, 0.6931471805599453, 1.0e-9))
}

///|
test "mda_lmi_constant" {
  // 常数序列 H(X)=0 → 返回 0.0
  let mi = @src.mda_lmi([1.0, 1.0, 1.0], [0.0, 1.0, 0.0])
  assert_true(@moonbit.assert_num_close(mi, 0.0, 1.0e-9))
}

///|
test "mda_lmi_entropy_hist" {
  // H([0,0,1,1]) = ln 2
  let h = @src.mda_entropy_hist([0.0, 0.0, 1.0, 1.0])
  assert_true(@moonbit.assert_num_close(h, 0.6931471805599453, 1.0e-9))
  // 常数序列熵为 0
  let h0 = @src.mda_entropy_hist([2.0, 2.0, 2.0])
  assert_true(@moonbit.assert_num_close(h0, 0.0, 1.0e-9))
}

///|
test "mda_lmi_matrix_and_normalized" {
  let s1 = [0.0, 0.0, 1.0, 1.0]
  let s2 = [0.0, 1.0, 0.0, 1.0]
  let m = @src.mda_lmi_matrix([s1, s2])
  assert_eq(m.length(), 2)
  // 对角 = H(X) = ln 2
  assert_true(@moonbit.assert_num_close(m[0][0], 0.6931471805599453, 1.0e-9))
  assert_true(@moonbit.assert_num_close(m[1][1], 0.6931471805599453, 1.0e-9))
  // 对称
  assert_true(@moonbit.assert_num_close(m[0][1], m[1][0], 1.0e-12))
  // 独立信号对的互信息 ≈ 0
  assert_true(@moonbit.assert_num_close(m[0][1], 0.0, 1.0e-9))
  // 归一化：完全相关 → 1；零熵分母 → 0
  let n1 = @src.mda_lmi_normalized(
    0.6931471805599453, 0.6931471805599453, 0.6931471805599453,
  )
  assert_true(@moonbit.assert_num_close(n1, 1.0, 1.0e-9))
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_lmi_normalized(0.5, 0.0, 1.0),
      0.0,
      1.0e-12,
    ),
  )
}
