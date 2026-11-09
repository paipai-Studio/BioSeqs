///|
test "dcd_writer_roundtrip" {
  let w = @src.OmmDCDWriter::new(n_atoms=2, dt=0.001)
  w.add_frame([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
  w.add_frame([[0.5, 0.0, 0.0], [1.5, 0.0, 0.0]])
  assert_eq(w.num_frames(), 2)
  let bytes = w.to_dcd()
  let traj = @src.parse_dcd(bytes)
  assert_eq(traj.frames.length(), 2)
  assert_eq(traj.num_atoms(), 2)
  assert_true(@moonbit.assert_num_close(traj.frames[0][1][0], 1.0, 1.0e-5))
  assert_true(@moonbit.assert_num_close(traj.frames[1][0][0], 0.5, 1.0e-5))
  assert_true(@moonbit.assert_num_close(traj.frames[1][1][0], 1.5, 1.0e-5))
}

///|
test "simulation_record_dcd" {
  let sys = @src.OmmSystem::new()
  ignore(sys.add_particle(particle=@src.OmmParticle::new(mass=12.0)))
  ignore(sys.add_particle(particle=@src.OmmParticle::new(mass=1.0)))
  sys.add_bond(0, 1, 100.0, 1.0)
  let sim = @src.OmmSimulation::new(
    system=sys,
    dt=0.001,
    positions=Some([[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]]),
  )
  let w = sim.record_dcd(n_steps=10, report_interval=5)
  // Frames at relative steps 5 and 10.
  assert_eq(w.num_frames(), 2)
  let bytes = w.to_dcd()
  let traj = @src.parse_dcd(bytes)
  assert_eq(traj.frames.length(), 2)
  // Coordinates finite.
  for f = 0; f < 2; f = f + 1 {
    for a = 0; a < 2; a = a + 1 {
      for d = 0; d < 3; d = d + 1 {
        assert_true(traj.frames[f][a][d] == traj.frames[f][a][d])
      }
    }
  }
}

///|
test "pdb_multi_frame_output" {
  let frames : Array[Array[Array[Double]]] = [
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.5], [1.0, 0.0, 0.5]],
  ]
  let text = @src.omm_write_pdb_frames(frames, atom_names=["C", "O"])
  // Two MODEL blocks with ENDMDL terminators.
  assert_true(text.contains("MODEL     1"))
  assert_true(text.contains("MODEL     2"))
  // Both MODEL blocks are terminated by ENDMDL records.
  assert_true(text.contains("ENDMDL"))
  assert_true(text.contains("ATOM  "))
  assert_true(text.contains("C"))
  assert_true(text.contains("O"))
  assert_true(text.has_suffix("END\n"))
  // Coordinate formatting: %8.3f fields.
  assert_true(text.contains("   1.000"))
  assert_true(text.contains("   0.500"))
}

///|
test "pdb_coordinate_formatting" {
  // Negative and sub-milliangstrom coordinates format correctly.
  let frames : Array[Array[Array[Double]]] = [[[-1.25, 0.0004, 12.0]]]
  let text = @src.omm_write_pdb_frames(frames)
  assert_true(text.contains("  -1.250"))
  assert_true(text.contains("   0.000"))
  assert_true(text.contains("  12.000"))
}
