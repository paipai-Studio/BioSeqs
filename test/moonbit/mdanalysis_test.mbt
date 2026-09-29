///|
/// Tests for MDAnalysis-style trajectory analysis (mdanalysis.mbt).

test "mdanalysis_distance_array" {
  let c1 = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let c2 = [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let d = @src.mdanalysis_distance_array(c1, c2)
  assert_true(assert_close(d[0][0], 0.0, 1.0e-9))
  assert_true(assert_close(d[0][1], 1.0, 1.0e-9))
  assert_true(assert_close(d[1][0], 1.0, 1.0e-9))
  assert_true(assert_close(d[1][1], 2.0.sqrt(), 1.0e-9))
}

test "mdanalysis_rmsd_translation_invariant" {
  // 非共面 4 点，b 平移 +5 → 叠加后 RMSD = 0
  let a = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  let b = [[5.0, 5.0, 5.0], [6.0, 5.0, 5.0], [5.0, 6.0, 5.0], [5.0, 5.0, 6.0]]
  assert_true(assert_close(@src.mdanalysis_rmsd(a, b), 0.0, 1.0e-6))
}

test "mdanalysis_rmsf" {
  let traj = [
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
  ]
  let rmsf = @src.mdanalysis_rmsf(traj)
  // 每原子偏离均值 ±0.5 → rmsf = 0.5
  assert_true(assert_close(rmsf[0], 0.5, 1.0e-9))
  assert_true(assert_close(rmsf[1], 0.5, 1.0e-9))
}

test "mdanalysis_rdf_peak" {
  // 原点一个原子，距离 1 处两个原子；range_max=2, nbins=2 → 峰在 bin 1
  let c1 = [[0.0, 0.0, 0.0]]
  let c2 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let (centers, g) = @src.mdanalysis_rdf(c1, c2, 2, 2.0)
  assert_eq(centers.length(), 2)
  assert_true(assert_close(g[0], 0.0, 1.0e-9))
  assert_true(g[1] > 0.0)
}

test "mdanalysis_msd_uniform" {
  // 匀速运动 x(t)=t → MSD = t²
  let traj = [
    [[0.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0]],
    [[2.0, 0.0, 0.0]],
    [[3.0, 0.0, 0.0]],
  ]
  let msd = @src.mdanalysis_msd(traj)
  assert_true(assert_close(msd[0], 0.0, 1.0e-9))
  assert_true(assert_close(msd[1], 1.0, 1.0e-9))
  assert_true(assert_close(msd[2], 4.0, 1.0e-9))
  assert_true(assert_close(msd[3], 9.0, 1.0e-9))
}
