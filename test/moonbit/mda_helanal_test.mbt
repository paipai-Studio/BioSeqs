///|
/// Tests for HELANAL helix axis ensemble analysis (mda_helanal.mbt).
test "mda_helix_axis_window_straight" {
  // 4 点沿 x 轴排列 → 拟合单位轴 [1, 0, 0]（符号规范化后）
  let pts = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
  let axis = @src.mda_helix_axis_window(pts)
  assert_true(@moonbit.assert_num_close(axis[0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(axis[1], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(axis[2], 0.0, 1.0e-9))
}

///|
test "mda_helanal_ideal_helix" {
  // 理想陡螺旋：r=2.3，每残基转 100°（dt=1.7453 rad），z 每弧度上升 3.0
  let dt = 1.7453292519943295
  let n_res = 10
  let frame : Array[Array[Double]] = Array(capacity=n_res)
  for i = 0; i < n_res; i = i + 1 {
    let t = i.to_double() * dt
    frame.push([2.3 * @math.cos(t), 2.3 * @math.sin(t), 3.0 * t])
  }
  let frames = @src.mda_helanal([frame])
  assert_eq(frames.length(), 1)
  let f = frames[0]
  // 全局轴接近 z 方向
  assert_true(f.global_axis[2] > 0.9)
  // 窗口数 = n_res - 3 = 7；twist 数 = 6
  assert_eq(f.window_axes.length(), 7)
  assert_eq(f.bend_angles.length(), 7)
  assert_eq(f.twist_angles.length(), 6)
  // 理想螺旋的 bend 角应很小
  for b in f.bend_angles {
    assert_true(b < 15.0)
  }
  // twist 角度合理（> 0）
  for tw in f.twist_angles {
    assert_true(tw > 0.0)
  }
}

///|
test "mda_helanal_short_and_straight" {
  // 直线 5 点：窗口 2 个、twist 1 个、bend 全 0；3 点片段（< 4 残基）→ 空结果
  let straight : Array[Array[Double]] = Array(capacity=5)
  for i = 0; i < 5; i = i + 1 {
    straight.push([0.0, 0.0, i.to_double()])
  }
  let short = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
  let frames = @src.mda_helanal([straight, short])
  assert_eq(frames.length(), 2)
  let f0 = frames[0]
  assert_eq(f0.window_axes.length(), 2)
  assert_eq(f0.twist_angles.length(), 1)
  for b in f0.bend_angles {
    assert_true(@moonbit.assert_num_close(b, 0.0, 1.0e-9))
  }
  let f1 = frames[1]
  assert_eq(f1.window_axes.length(), 0)
  assert_eq(f1.bend_angles.length(), 0)
  assert_eq(f1.twist_angles.length(), 0)
  assert_true(@moonbit.assert_num_close(f1.global_axis[0], 0.0, 1.0e-9))
}
