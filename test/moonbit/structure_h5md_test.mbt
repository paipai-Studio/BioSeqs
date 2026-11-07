///|
/// Tests for the H5MD object layer (structure_h5md.mbt).

///|
test "h5md_validate_ok" {
  let f = @src.h5md_file_new(author="Tester")
  let pos = @src.H5mdDataset::{
    name: "position",
    step: [0, 1],
    time: [0.0, 1.0],
    value: [
      [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
      [[1.0, 0.0, 0.0], [2.0, 1.0, 1.0]],
    ],
  }
  let vel = @src.H5mdDataset::{
    name: "velocity",
    step: [0, 1],
    time: [0.0, 1.0],
    value: [
      [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0]],
      [[0.2, 0.0, 0.0], [0.0, 0.2, 0.0]],
    ],
  }
  let box_edges = [
    [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
    [[12.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 12.0]],
  ]
  let idx = @src.h5md_add_particles(
    f,
    "prot",
    pos,
    velocity=Some(vel),
    box_edges=Some(box_edges),
  )
  assert_eq(idx, 0)
  let found = @src.h5md_find_particles(f, "prot")
  match found {
    Some(pg) => assert_eq(pg.name, "prot")
    None => fail("particle group 'prot' not found")
  }
  assert_eq(@src.h5md_find_particles(f, "missing") is None, true)
  let errors = @src.h5md_validate(f)
  assert_eq(errors.length(), 0)
}

///|
test "h5md_validate_frame_mismatch" {
  let f = @src.h5md_file_new()
  // step/time carry 3 entries but value has only 2 frames.
  let pos = @src.H5mdDataset::{
    name: "position",
    step: [0, 1, 2],
    time: [0.0, 1.0, 2.0],
    value: [[[0.0, 0.0, 0.0]], [[1.0, 0.0, 0.0]]],
  }
  let _ = @src.h5md_add_particles(f, "prot", pos)
  let errors = @src.h5md_validate(f)
  assert_eq(errors.length() > 0, true)
  assert_true(errors[0].contains("position"))
}

///|
test "h5md_validate_time_length_mismatch" {
  let f = @src.h5md_file_new()
  let pos = @src.H5mdDataset::{
    name: "position",
    step: [0, 1],
    time: [0.0, 1.0],
    value: [[[0.0, 0.0, 0.0]], [[1.0, 0.0, 0.0]]],
  }
  // velocity: two steps but only one time entry.
  let vel = @src.H5mdDataset::{
    name: "velocity",
    step: [0, 1],
    time: [0.0],
    value: [[[0.1, 0.0, 0.0]], [[0.2, 0.0, 0.0]]],
  }
  let _ = @src.h5md_add_particles(f, "prot", pos, velocity=Some(vel))
  let errors = @src.h5md_validate(f)
  assert_eq(errors.length() > 0, true)
  let mut mentions_velocity = false
  for e in errors {
    if e.contains("velocity") {
      mentions_velocity = true
    }
  }
  assert_eq(mentions_velocity, true)
}

///|
test "h5md_trajectory_roundtrip" {
  let frames = [[[1.0, 2.0, 3.0]], [[4.0, 5.0, 6.0]], [[7.0, 8.0, 9.0]]]
  let traj = @src.Trajectory::new(frames)
  let ds = @src.trajectory_to_h5md_dataset("position", traj, time_step=0.5)
  assert_eq(ds.step.length(), 3)
  assert_eq(ds.step[2], 2)
  assert_true(@moonbit.assert_num_close(ds.time[2], 1.0, 1.0e-9))
  let pg = @src.H5mdParticleGroup::{
    name: "p",
    position: ds,
    velocity: None,
    force: None,
    box_edges: None,
  }
  let traj2 = @src.h5md_to_trajectory(pg)
  assert_eq(traj2.num_frames(), 3)
  assert_eq(traj2.num_atoms(), 1)
  let frame = traj2.get_frame(1)
  assert_true(@moonbit.assert_num_close(frame[0][1], 5.0, 1.0e-9))
  let frame2 = traj2.get_frame(2)
  assert_true(@moonbit.assert_num_close(frame2[0][0], 7.0, 1.0e-9))
}

///|
test "h5md_to_json_fields" {
  let f = @src.h5md_file_new(author="Ada", creator_version="1.2.3")
  let pos = @src.H5mdDataset::{
    name: "position",
    step: [0, 1],
    time: [0.0, 0.5],
    value: [
      [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
      [[0.5, 0.0, 0.0], [1.5, 0.0, 0.0], [2.5, 0.0, 0.0], [3.5, 0.0, 0.0]],
    ],
  }
  let vel = @src.H5mdDataset::{
    name: "velocity",
    step: [0, 1],
    time: [0.0, 0.5],
    value: [
      [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
      [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
    ],
  }
  let _ = @src.h5md_add_particles(f, "prot", pos, velocity=Some(vel))
  let json = @src.h5md_to_json(f)
  assert_true(json.contains("\"h5md_version\":[1,1]"))
  assert_true(json.contains("\"author\":\"Ada\""))
  assert_true(json.contains("\"creator_version\":\"1.2.3\""))
  assert_true(json.contains("\"name\":\"prot\""))
  assert_true(json.contains("\"n_frames\":2"))
  assert_true(json.contains("\"n_atoms\":4"))
  assert_true(json.contains("\"has_velocity\":true"))
  assert_true(json.contains("\"has_force\":false"))
}

///|
test "h5md_add_duplicate_raises" {
  let f = @src.h5md_file_new()
  let pos = @src.H5mdDataset::{
    name: "position",
    step: [0],
    time: [0.0],
    value: [[[0.0, 0.0, 0.0]]],
  }
  let idx = @src.h5md_add_particles(f, "prot", pos)
  assert_eq(idx, 0)
  let raised = try {
    let _ = @src.h5md_add_particles(f, "prot", pos)
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}
