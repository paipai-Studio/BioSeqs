test "xtc_write_read_roundtrip" {
  // Create a simple frame with 3 atoms
  let frame = @src.XtcFrame::{
    step: 100,
    time: 10.0,
    natoms: 3,
    box: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0],
    coords: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
  }
  let bytes = @src.xtc_write_frame(frame)
  // Read it back
  let frames = @src.xtc_read(bytes)
  assert_eq(frames.length(), 1)
  let f = frames[0]
  assert_eq(f.step, 100)
  assert_eq(f.natoms, 3)
  assert_true(@moonbit.assert_close(f.time, 10.0, 0.001))
  // Check coords (precision = 1000, so tolerance ~0.001)
  assert_true(@moonbit.assert_close(f.coords[0], 0.1, 0.002))
  assert_true(@moonbit.assert_close(f.coords[1], 0.2, 0.002))
  assert_true(@moonbit.assert_close(f.coords[8], 0.9, 0.002))
}

test "xtc_read_multiple_frames" {
  let frame1 = @src.XtcFrame::{
    step: 0,
    time: 0.0,
    natoms: 2,
    box: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
    coords: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
  }
  let frame2 = @src.XtcFrame::{
    step: 1000,
    time: 100.0,
    natoms: 2,
    box: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
    coords: [0.5, 0.5, 0.5, 1.5, 1.5, 1.5],
  }
  let bytes1 = @src.xtc_write_frame(frame1)
  let bytes2 = @src.xtc_write_frame(frame2)
  // Concatenate
  let combined : Array[Int] = Array(capacity=bytes1.length() + bytes2.length())
  for i = 0; i < bytes1.length(); i = i + 1 {
    combined.push(bytes1[i])
  }
  for i = 0; i < bytes2.length(); i = i + 1 {
    combined.push(bytes2[i])
  }
  let frames = @src.xtc_read(combined)
  assert_eq(frames.length(), 2)
  assert_eq(frames[0].step, 0)
  assert_eq(frames[1].step, 1000)
}

test "xtc_frame_box" {
  let frame = @src.XtcFrame::{
    step: 0,
    time: 0.0,
    natoms: 1,
    box: [3.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 5.0],
    coords: [0.0, 0.0, 0.0],
  }
  let bytes = @src.xtc_write_frame(frame)
  let frames = @src.xtc_read(bytes)
  assert_true(@moonbit.assert_close(frames[0].box[0], 3.0, 0.001))
  assert_true(@moonbit.assert_close(frames[0].box[4], 4.0, 0.001))
  assert_true(@moonbit.assert_close(frames[0].box[8], 5.0, 0.001))
}
