///|
test "omm periodic torsion trans conformation energy" {
  // phi = pi ; E = k*(1 + cos(n*phi - phi0))
  // n=1, phase=0, k=2 -> 2*(1+cos(pi)) = 0
  let t = @src.omm_amber_periodic_torsion(0, 1, 2, 3, 2.0, 0.0, 1)
  let coords = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (2.0, 1.0, 0.0),
  ]
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_periodic_torsion_energy(t, coords),
      0.0,
      1.0e-10,
    ),
  )
}

///|
test "omm periodic torsion trans periodicity two" {
  // n=2, phase=0, k=1.5 -> 1.5*(1+cos(2*pi)) = 3
  let t = @src.omm_amber_periodic_torsion(0, 1, 2, 3, 1.5, 0.0, 2)
  let coords = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (2.0, 1.0, 0.0),
  ]
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_periodic_torsion_energy(t, coords),
      3.0,
      1.0e-10,
    ),
  )
}

///|
test "omm periodic torsion trans with 180 degree phase" {
  // n=1, phase=180, k=1 -> 1*(1+cos(pi - pi)) = 2
  let t = @src.omm_amber_periodic_torsion(0, 1, 2, 3, 1.0, 180.0, 1)
  let coords = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (2.0, 1.0, 0.0),
  ]
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_periodic_torsion_energy(t, coords),
      2.0,
      1.0e-10,
    ),
  )
}

///|
test "omm periodic torsion cis conformation energy" {
  // phi = 0 ; n=1, phase=0, k=1 -> 2 ; n=1, phase=180, k=2 -> 0
  let coords = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
  ]
  let t1 = @src.omm_amber_periodic_torsion(0, 1, 2, 3, 1.0, 0.0, 1)
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_periodic_torsion_energy(t1, coords),
      2.0,
      1.0e-10,
    ),
  )
  let t2 = @src.omm_amber_periodic_torsion(0, 1, 2, 3, 2.0, 180.0, 1)
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_periodic_torsion_energy(t2, coords),
      0.0,
      1.0e-10,
    ),
  )
}

///|
test "omm periodic torsion orthogonal conformation energy" {
  // phi = pi/2 ; n=1, phase=0, k=1 -> 1 + cos(pi/2) = 1
  let coords = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (1.0, 1.0, 1.0),
  ]
  let t1 = @src.omm_amber_periodic_torsion(0, 1, 2, 3, 1.0, 0.0, 1)
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_periodic_torsion_energy(t1, coords),
      1.0,
      1.0e-10,
    ),
  )
  // n=4, phase=0, k=1 -> 1 + cos(4*pi/2) = 2
  let t2 = @src.omm_amber_periodic_torsion(0, 1, 2, 3, 1.0, 0.0, 4)
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_periodic_torsion_energy(t2, coords),
      2.0,
      1.0e-10,
    ),
  )
}

///|
test "omm periodic torsion total sums terms" {
  let coords = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
  ]
  let t1 = @src.omm_amber_periodic_torsion(0, 1, 2, 3, 1.0, 0.0, 1)
  let t2 = @src.omm_amber_periodic_torsion(0, 1, 2, 3, 2.0, 180.0, 1)
  let total = @src.omm_periodic_torsion_total([t1, t2], coords)
  // t1 -> 2, t2 -> 0
  assert_true(@moonbit.assert_num_close(total, 2.0, 1.0e-10))
}

///|
test "omm periodic torsion forces vanish at stationary phases" {
  // phi = pi, n=1, phase=0 -> sin(n*phi - phi0) = 0 -> zero forces
  let t = @src.omm_amber_periodic_torsion(0, 1, 2, 3, 1.0, 0.0, 1)
  let coords = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (2.0, 1.0, 0.0),
  ]
  let f = @src.omm_periodic_torsion_forces(t, coords)
  assert_eq(f.length(), 12)
  for i = 0; i < 12; i = i + 1 {
    assert_true(@moonbit.assert_num_close(f[i], 0.0, 1.0e-12))
  }
}

///|
test "omm periodic torsion forces orthogonal conformation" {
  // phi = pi/2, n=1, phase=0, k=1: dE/dphi = -sin(pi/2) = -1
  // f = dE/dphi * 0.25 = -0.25 ; forces i,j = +0.25 ; k,l = -0.25 (x component)
  let t = @src.omm_amber_periodic_torsion(0, 1, 2, 3, 1.0, 0.0, 1)
  let coords = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (1.0, 1.0, 1.0),
  ]
  let f = @src.omm_periodic_torsion_forces(t, coords)
  assert_true(@moonbit.assert_num_close(f[0], 0.25, 1.0e-12))
  assert_true(@moonbit.assert_num_close(f[3], 0.25, 1.0e-12))
  assert_true(@moonbit.assert_num_close(f[6], -0.25, 1.0e-12))
  assert_true(@moonbit.assert_num_close(f[9], -0.25, 1.0e-12))
  // all y/z components zero and total force zero
  let mut sum = 0.0
  for i = 0; i < 12; i = i + 1 {
    sum = sum + f[i]
  }
  assert_true(@moonbit.assert_num_close(sum, 0.0, 1.0e-12))
}

///|
test "omm amber periodic torsion converts degrees to radians" {
  // phase 90 deg -> phi0 = pi/2 ; at phi=0: E = k*(1+cos(-pi/2)) = k
  let t = @src.omm_amber_periodic_torsion(0, 1, 2, 3, 1.0, 90.0, 1)
  let coords = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
  ]
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_periodic_torsion_energy(t, coords),
      1.0,
      1.0e-10,
    ),
  )
  assert_true(
    @moonbit.assert_num_close(t.phi0, 3.141592653589793 / 2.0, 1.0e-12),
  )
}
