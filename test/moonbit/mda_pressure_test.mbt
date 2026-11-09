///|
/// Tests for mda_pressure.mbt (pressure tensor, virial, derived quantities).

///|
test "mda_pressure_tensor_kinetic_only" {
  // Zero positions/forces -> virial 0. Velocities (1,0,0) m=2 and
  // (0,2,0) m=1 give kinetic tensor [[2,0,0],[0,4,0],[0,0,0]]; V = 2.
  let nul = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let p = @src.mda_pressure_tensor(
    nul,
    [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
    nul,
    [2.0, 1.0],
    2.0,
  )
  assert_true(@moonbit.assert_num_close(p[0][0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[1][1], 2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[2][2], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[0][1], 0.0, 1.0e-12))
  assert_true(
    @moonbit.assert_num_close(@src.mda_scalar_pressure(p), 1.0, 1.0e-12),
  )
  assert_true(
    @moonbit.assert_num_close(@src.mda_pressure_anisotropy(p), 1.5, 1.0e-12),
  )
}

///|
test "mda_pressure_tensor_virial_only" {
  // Zero velocities -> kinetic 0. One atom r=(1,0,0), F=(3,2,1):
  // virial = [[3,2,1],[0,0,0],[0,0,0]]; V = 1.
  let nul = [[0.0, 0.0, 0.0]]
  let p = @src.mda_pressure_tensor(
    [[1.0, 0.0, 0.0]],
    nul,
    [[3.0, 2.0, 1.0]],
    [1.0],
    1.0,
  )
  assert_true(@moonbit.assert_num_close(p[0][0], 3.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[0][1], 2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[0][2], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[1][1], 0.0, 1.0e-12))
}

///|
test "mda_pressure_tensor_volume_scaling" {
  // Same system with V = 4: every component is scaled by 1/4.
  let nul = [[0.0, 0.0, 0.0]]
  let p = @src.mda_pressure_tensor(
    [[1.0, 0.0, 0.0]],
    nul,
    [[3.0, 2.0, 1.0]],
    [1.0],
    4.0,
  )
  assert_true(@moonbit.assert_num_close(p[0][0], 0.75, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[0][1], 0.5, 1.0e-12))
}

///|
test "mda_pressure_tensor_zero_volume_guard" {
  let p = @src.mda_pressure_tensor(
    [[1.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0]],
    [1.0],
    0.0,
  )
  assert_true(@moonbit.assert_num_close(p[0][0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[2][2], 0.0, 1.0e-12))
}

///|
test "mda_pressure_tensor_symmetric_offdiagonal" {
  // Velocity (1,1,0) m=1, V=1: kinetic [[1,1,0],[1,1,0],[0,0,0]].
  let nul = [[0.0, 0.0, 0.0]]
  let p = @src.mda_pressure_tensor(nul, [[1.0, 1.0, 0.0]], nul, [1.0], 1.0)
  assert_true(@moonbit.assert_num_close(p[0][1], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[1][0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[0][0], 1.0, 1.0e-12))
}

///|
test "mda_surface_tension_membrane" {
  // diag(2,2,1) tensor, box_z = 4: gamma = 4 * (2 - 1) / 2 = 2.
  let tensor = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.0]]
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_surface_tension(tensor, 4.0),
      2.0,
      1.0e-12,
    ),
  )
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_pressure_anisotropy(tensor),
      1.0,
      1.0e-12,
    ),
  )
}

///|
test "mda_von_mises_stress_hydrostatic" {
  // Isotropic pressure has zero deviatoric part -> stress 0.
  let tensor = [[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]]
  assert_true(
    @moonbit.assert_num_close(@src.mda_von_mises_stress(tensor), 0.0, 1.0e-12),
  )
}

///|
test "mda_von_mises_stress_uniaxial" {
  // diag(3,0,0): p = 1; sxx = 2, syy = szz = -1.
  // vms = sqrt(0.5*(9 + 0 + 9)) = 3.
  let tensor = [[3.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  assert_true(
    @moonbit.assert_num_close(@src.mda_von_mises_stress(tensor), 3.0, 1.0e-9),
  )
}

///|
test "mda_von_mises_stress_pure_shear" {
  // Only off-diagonal sxy = 1: vms = sqrt(3).
  let tensor = [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_von_mises_stress(tensor),
      3.0.sqrt(),
      1.0e-9,
    ),
  )
}
