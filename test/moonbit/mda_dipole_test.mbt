///|
/// Tests for trajectory dipole moments (mda_dipole.mbt).
test "mda_dipole_single_charge" {
  // 单电荷 q=1 位于 (1,0,0) → μ=[1,0,0]
  let mu = @src.mda_dipole([[1.0, 0.0, 0.0]], [1.0])
  assert_true(@moonbit.assert_num_close(mu[0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(mu[1], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(mu[2], 0.0, 1.0e-9))
  // Debye 转换系数 4.80320425
  let d = @src.mda_dipole_debye(mu)
  assert_true(@moonbit.assert_num_close(d[0], 4.80320425, 1.0e-9))
  assert_true(@moonbit.assert_num_close(d[1], 0.0, 1.0e-12))
  // 模长 = sqrt(mu . mu)
  assert_true(
    @moonbit.assert_num_close(@src.mda_dipole_magnitude(mu), 1.0, 1.0e-9),
  )
}

///|
test "mda_dipole_neutral_pair" {
  // +1 在 (1,0,0)、-1 在原点 → μ=[1,0,0]
  let mu = @src.mda_dipole([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]], [1.0, -1.0])
  assert_true(@moonbit.assert_num_close(mu[0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(mu[1], 0.0, 1.0e-9))
  assert_true(
    @moonbit.assert_num_close(@src.mda_dipole_magnitude(mu), 1.0, 1.0e-9),
  )
}

///|
test "mda_dipole_traj" {
  // 两帧相同构型 → 每帧偶极均为 [1,0,0]
  let traj = [
    [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
  ]
  let mus = @src.mda_dipole_traj(traj, [1.0, -1.0])
  assert_eq(mus.length(), 2)
  assert_true(@moonbit.assert_num_close(mus[0][0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(mus[1][0], 1.0, 1.0e-9))
}

///|
test "mda_dipole_correlation_constant" {
  // 平移不变序列：C(t) = 1 对所有 t（含 C(0)=1）
  let dipoles = [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let c = @src.mda_dipole_correlation(dipoles)
  assert_eq(c.length(), 3)
  for v in c {
    assert_true(@moonbit.assert_num_close(v, 1.0, 1.0e-9))
  }
}

///|
test "mda_dipole_correlation_zero" {
  // 全零偶极 → 分母为 0 → 全 0
  let c = @src.mda_dipole_correlation([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
  assert_eq(c.length(), 2)
  for v in c {
    assert_true(@moonbit.assert_num_close(v, 0.0, 1.0e-12))
  }
}
