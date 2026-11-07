///|
/// Tests for site-specific RDF and 2D RDF (mda_rdf_site.mbt).
test "mda_rdf_sitewise_peak" {
  // 单帧 1 对原子相距 1.5：bin 宽 1.0 → 计数落在 bin 1
  let traj = [[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]]
  let (centers, counts) = @src.mda_rdf_sitewise(traj, [0], [1], 3, 3.0)
  assert_eq(centers.length(), 3)
  assert_true(@moonbit.assert_num_close(centers[1], 1.5, 1.0e-9))
  assert_eq(counts.length(), 1)
  assert_true(@moonbit.assert_num_close(counts[0][1], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(counts[0][0], 0.0, 1.0e-9))
  // 归一化 g(r)：无计数 bin 为 0，峰 bin > 0
  let g = @src.mda_rdf_sitewise_normalize(counts, 1, 1, 1, 3, 3.0)
  assert_true(@moonbit.assert_num_close(g[0][0], 0.0, 1.0e-9))
  assert_true(g[0][1] > 0.0)
}

///|
test "mda_rdf_sitewise_pair_order" {
  // sel1=[0,1] 外层、sel2=[2,3] 内层：d(0,2)=2 → counts[0][4]=1；
  // d(0,3)=3 >= 2.5 超出；d(1,*)=10 超出 → 其余全 0
  let frame = [
    [0.0, 0.0, 0.0],
    [10.0, 0.0, 0.0],
    [0.0, 0.0, 2.0],
    [0.0, 0.0, 3.0],
  ]
  let (_, counts) = @src.mda_rdf_sitewise([frame], [0, 1], [2, 3], 5, 2.5)
  assert_eq(counts.length(), 4)
  assert_true(@moonbit.assert_num_close(counts[0][4], 1.0, 1.0e-9))
  let mut total_rest = 0.0
  for p = 1; p < 4; p = p + 1 {
    for b = 0; b < 5; b = b + 1 {
      total_rest = total_rest + counts[p][b]
    }
  }
  assert_true(@moonbit.assert_num_close(total_rest, 0.0, 1.0e-12))
}

///|
test "mda_rdf_2d_fold" {
  // 原子 (0,0,0) 与 (1,1,0)：|dx|=|dy|=1，bin 宽 1.0 → 象限折叠后落 bin (1,1)
  let traj = [[[0.0, 0.0, 0.0], [1.0, 1.0, 0.0]]]
  let (xc, yc, hist) = @src.mda_rdf_2d(traj, [0], [1], 2.0, 2.0, 2, 2)
  assert_true(@moonbit.assert_num_close(xc[0], 0.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(yc[1], 1.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(hist[1][1], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(hist[0][0], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(hist[0][1], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(hist[1][0], 0.0, 1.0e-9))
}

///|
test "mda_rdf_2d_pbc" {
  // box=[10,10,10]：dx=-9.8 最小镜像折叠为 0.2 → bin 0
  let traj = [[[0.0, 0.0, 0.0], [9.8, 0.0, 0.0]]]
  let (_, _, hist) = @src.mda_rdf_2d(
    traj,
    [0],
    [1],
    2.0,
    2.0,
    2,
    2,
    box=Some([10.0, 10.0, 10.0]),
  )
  assert_true(@moonbit.assert_num_close(hist[0][0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(hist[1][0], 0.0, 1.0e-9))
}
