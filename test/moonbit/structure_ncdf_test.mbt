///|
/// Tests for the AMBER NetCDF (CDF-1) reader/writer (structure_ncdf.mbt).
test "ncdf_write_read_roundtrip_with_box" {
  let frames = [
    @src.NcdfFrame::{
      coords: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
      box_lengths: [10.0, 11.0, 12.0],
      box_angles: [90.0, 90.0, 90.0],
      time: 0.0,
    },
    @src.NcdfFrame::{
      coords: [1.1, 1.2, 1.3, 1.4, 1.5, 1.6],
      box_lengths: [10.5, 11.5, 12.5],
      box_angles: [89.0, 90.0, 91.0],
      time: 1.0,
    },
  ]
  let data = @src.ncdf_write_amber(frames, 2)
  let back = @src.ncdf_read_amber(data)
  assert_eq(back.length(), 2)
  // float32 payload: 1e-5 tolerance
  assert_true(@moonbit.assert_num_close(back[0].coords[0], 0.1, 1.0e-5))
  assert_true(@moonbit.assert_num_close(back[0].coords[5], 0.6, 1.0e-5))
  assert_true(@moonbit.assert_num_close(back[1].coords[0], 1.1, 1.0e-5))
  assert_true(@moonbit.assert_num_close(back[1].coords[5], 1.6, 1.0e-5))
  // double payload: exact
  assert_true(@moonbit.assert_num_close(back[0].time, 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(back[1].time, 1.0, 1.0e-9))
  assert_eq(back[0].box_lengths.length(), 3)
  assert_eq(back[0].box_angles.length(), 3)
  assert_true(@moonbit.assert_num_close(back[1].box_lengths[0], 10.5, 1.0e-5))
  assert_true(@moonbit.assert_num_close(back[1].box_lengths[2], 12.5, 1.0e-5))
  assert_true(@moonbit.assert_num_close(back[1].box_angles[2], 91.0, 1.0e-5))
}

///|
test "ncdf_write_read_no_box" {
  let frames = [
    @src.NcdfFrame::{
      coords: [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
      box_lengths: [],
      box_angles: [],
      time: 0.0,
    },
    @src.NcdfFrame::{
      coords: [2.0, 2.0, 2.0, 3.0, 3.0, 3.0],
      box_lengths: [],
      box_angles: [],
      time: 2.0,
    },
  ]
  let data = @src.ncdf_write_amber(frames, 2)
  let back = @src.ncdf_read_amber(data)
  assert_eq(back.length(), 2)
  assert_eq(back[0].box_lengths.length(), 0)
  assert_eq(back[0].box_angles.length(), 0)
  assert_true(@moonbit.assert_num_close(back[1].coords[0], 2.0, 1.0e-5))
  assert_true(@moonbit.assert_num_close(back[1].coords[5], 3.0, 1.0e-5))
  assert_true(@moonbit.assert_num_close(back[1].time, 2.0, 1.0e-9))
}

///|
test "ncdf_to_trajectory" {
  let frames = [
    @src.NcdfFrame::{
      coords: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
      box_lengths: [],
      box_angles: [],
      time: 0.0,
    },
    @src.NcdfFrame::{ coords: [], box_lengths: [], box_angles: [], time: 1.0, },
  ]
  let traj = @src.ncdf_to_trajectory(frames)
  // the frame without coordinates is skipped
  assert_eq(traj.num_frames(), 1)
  assert_eq(traj.num_atoms(), 2)
  let frame0 = traj.get_frame(0)
  assert_true(@moonbit.assert_num_close(frame0[1][2], 0.6, 1.0e-9))
}

///|
test "ncdf_read_rejects_bad_magic" {
  let frames = [
    @src.NcdfFrame::{
      coords: [0.0, 0.0, 0.0],
      box_lengths: [],
      box_angles: [],
      time: 0.0,
    },
  ]
  let data = @src.ncdf_write_amber(frames, 1)
  data[0] = 88
  let mut raised = false
  try {
    let r = @src.ncdf_read_amber(data)
    ignore(r)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "ncdf_read_rejects_missing_coordinates" {
  // A truncated stream: valid magic but no variable list.
  let data : Array[Int32] = [67, 68, 70, 1, 0, 0, 0, 0]
  let mut raised = false
  try {
    let r = @src.ncdf_read_amber(data)
    ignore(r)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}
