///|
/// Tests for Trajectory + DCD write/read (dcd.mbt).

test "trajectory_accessors" {
  let frames = [
    [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]],
    [[10.0, 11.0, 12.0], [13.0, 14.0, 15.0]],
  ]
  let t = @src.Trajectory::new(frames)
  assert_eq(t.num_frames(), 2)
  assert_eq(t.num_atoms(), 2)
  assert_eq(t.get_frame(1)[0][0], 10.0)
}

test "dcd_roundtrip_integer_coords" {
  // 整数坐标在 float32 里精确表示 → 写读应逐位无损
  let frames = [
    [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]],
    [[10.0, -11.0, 12.0], [13.0, 14.0, -15.0]],
  ]
  let t = @src.Trajectory::new(frames)
  let bytes = @src.write_dcd(t)
  let back = @src.parse_dcd(bytes)
  assert_eq(back.num_frames(), 2)
  assert_eq(back.num_atoms(), 2)
  let mut f = 0
  while f < 2 {
    let mut a = 0
    while a < 2 {
      let mut d = 0
      while d < 3 {
        assert_eq(back.get_frame(f)[a][d], frames[f][a][d])
        d = d + 1
      }
      a = a + 1
    }
    f = f + 1
  }
}

test "dcd_roundtrip_fractional" {
  // 0.5 等二进分数也能精确表示
  let frames = [
    [[0.5, -0.25, 1.75], [2.5, -3.5, 0.0625]],
  ]
  let t = @src.Trajectory::new(frames)
  let back = @src.parse_dcd(@src.write_dcd(t))
  assert_eq(back.num_atoms(), 2)
  assert_true(assert_close(back.get_frame(0)[0][1], -0.25, 1.0e-9))
  assert_true(assert_close(back.get_frame(0)[1][2], 0.0625, 1.0e-9))
}
