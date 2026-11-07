///|
/// Tests for the GROMACS TRR reader/writer (structure_trr.mbt).
test "trr_write_read_roundtrip_with_velocities" {
  let f1 = @src.TrrFrame::{
    step: 10,
    time: 0.5,
    lambda: 0.25,
    box: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0],
    coords: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
    velocities: [0.01, -0.02, 0.03, -0.04, 0.05, -0.06],
    forces: [],
  }
  let f2 = @src.TrrFrame::{
    step: 20,
    time: 1.5,
    lambda: 0.75,
    box: [],
    coords: [1.1, 1.2, 1.3, 1.4, 1.5, 1.6],
    velocities: [],
    forces: [2.0, -2.5, 3.0, -3.5, 4.0, -4.5],
  }
  let b1 = @src.trr_write_frame(f1)
  let b2 = @src.trr_write_frame(f2)
  let combined : Array[Int32] = Array(capacity=b1.length() + b2.length())
  for b in b1 {
    combined.push(b)
  }
  for b in b2 {
    combined.push(b)
  }
  let frames = @src.trr_read(combined)
  assert_eq(frames.length(), 2)
  assert_eq(frames[0].step, 10)
  assert_eq(frames[1].step, 20)
  assert_true(@moonbit.assert_num_close(frames[0].time, 0.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(frames[0].lambda, 0.25, 1.0e-9))
  assert_true(@moonbit.assert_num_close(frames[1].time, 1.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(frames[1].lambda, 0.75, 1.0e-9))
  // float32 payload: 1e-5 tolerance
  assert_true(@moonbit.assert_num_close(frames[0].coords[0], 0.1, 1.0e-5))
  assert_true(@moonbit.assert_num_close(frames[0].coords[5], 0.6, 1.0e-5))
  assert_true(@moonbit.assert_num_close(frames[1].coords[3], 1.4, 1.0e-5))
  assert_true(@moonbit.assert_num_close(frames[0].velocities[1], -0.02, 1.0e-5))
  assert_true(@moonbit.assert_num_close(frames[0].box[8], 2.0, 1.0e-5))
  assert_true(@moonbit.assert_num_close(frames[1].forces[3], -3.5, 1.0e-5))
  assert_eq(frames[0].forces.length(), 0)
  assert_eq(frames[1].box.length(), 0)
  assert_eq(frames[1].velocities.length(), 0)
}

///|
test "trr_read_rejects_bad_magic" {
  let f = @src.TrrFrame::{
    step: 0,
    time: 0.0,
    lambda: 0.0,
    box: [],
    coords: [0.0, 0.0, 0.0],
    velocities: [],
    forces: [],
  }
  let buf = @src.trr_write_frame(f)
  buf[0] = 199
  let mut raised = false
  try {
    let r = @src.trr_read(buf)
    ignore(r)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "trr_empty_frame_roundtrip" {
  let f = @src.TrrFrame::{
    step: 1,
    time: 2.0,
    lambda: 0.0,
    box: [],
    coords: [],
    velocities: [],
    forces: [],
  }
  let buf = @src.trr_write_frame(f)
  let frames = @src.trr_read(buf)
  assert_eq(frames.length(), 1)
  assert_eq(frames[0].step, 1)
  assert_true(@moonbit.assert_num_close(frames[0].time, 2.0, 1.0e-9))
  assert_eq(frames[0].coords.length(), 0)
  assert_eq(frames[0].box.length(), 0)
}

///|
test "trr_to_trajectory" {
  let f1 = @src.TrrFrame::{
    step: 0,
    time: 0.0,
    lambda: 0.0,
    box: [],
    coords: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
    velocities: [],
    forces: [],
  }
  let f2 = @src.TrrFrame::{
    step: 1,
    time: 1.0,
    lambda: 0.0,
    box: [],
    coords: [],
    velocities: [],
    forces: [],
  }
  let f3 = @src.TrrFrame::{
    step: 2,
    time: 2.0,
    lambda: 0.0,
    box: [],
    coords: [1.1, 1.2, 1.3, 1.4, 1.5, 1.6],
    velocities: [],
    forces: [],
  }
  let traj = @src.trr_to_trajectory([f1, f2, f3])
  // the frame without coordinates is skipped
  assert_eq(traj.num_frames(), 2)
  assert_eq(traj.num_atoms(), 2)
  let frame0 = traj.get_frame(0)
  assert_true(@moonbit.assert_num_close(frame0[0][0], 0.1, 1.0e-9))
  assert_true(@moonbit.assert_num_close(frame0[1][2], 0.6, 1.0e-9))
  let frame1 = traj.get_frame(1)
  assert_true(@moonbit.assert_num_close(frame1[1][0], 1.4, 1.0e-9))
}

///|
test "trajectory_to_trr_and_back" {
  let frames_data = [
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
  ]
  let traj = @src.Trajectory::new(frames_data)
  let trr_frames = @src.trajectory_to_trr(traj, step=100, time_step=0.5)
  assert_eq(trr_frames.length(), 2)
  assert_eq(trr_frames[0].step, 100)
  assert_eq(trr_frames[1].step, 101)
  assert_true(@moonbit.assert_num_close(trr_frames[1].time, 0.5, 1.0e-9))
  assert_eq(trr_frames[0].coords.length(), 6)
  let buf : Array[Int32] = Array(capacity=0)
  for f in trr_frames {
    let b = @src.trr_write_frame(f)
    for x in b {
      buf.push(x)
    }
  }
  let back = @src.trr_read(buf)
  assert_eq(back.length(), 2)
  assert_true(@moonbit.assert_num_close(back[0].coords[3], 1.0, 1.0e-5))
  assert_true(@moonbit.assert_num_close(back[1].coords[0], 1.0, 1.0e-5))
}
