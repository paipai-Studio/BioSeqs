///|
test "omm pme params alpha from cutoff" {
  let params = @src.OmmPMEParams::new(8.0, box=[10.0, 10.0, 10.0])
  assert_true(
    @moonbit.assert_num_close(
      params.alpha,
      3.141592653589793.sqrt() / 8.0,
      1.0e-12,
    ),
  )
  assert_eq(params.cutoff, 8.0)
  assert_eq(params.grid_x, 32)
  assert_eq(params.grid_y, 32)
  assert_eq(params.grid_z, 32)
  assert_eq(params.box[0], 10.0)
}

///|
test "omm pme self energy analytic value" {
  // E_self = -alpha/sqrt(pi) * sum(q^2) ; alpha = sqrt(pi)/8, charges +1/-1
  // -> -alpha/sqrt(pi)*2 = -2/8 = -0.25
  let params = @src.OmmPMEParams::new(8.0)
  let e = @src.omm_pme_self_energy([1.0, -1.0], params)
  assert_true(
    @moonbit.assert_num_close(
      e,
      -(3.141592653589793.sqrt() / 8.0) / 3.141592653589793.sqrt() * 2.0,
      1.0e-12,
    ),
  )
  assert_true(@moonbit.assert_num_close(e, -0.25, 1.0e-12))
}

///|
test "omm pme self energy zero charges" {
  let params = @src.OmmPMEParams::new(8.0)
  let e = @src.omm_pme_self_energy([0.0, 0.0], params)
  assert_eq(e, 0.0)
}

///|
test "omm pme real space like charges repulsive and decaying" {
  let params = @src.OmmPMEParams::new(8.0, box=[0.0, 0.0, 0.0])
  let e_near = @src.omm_pme_real_space(
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    [1.0, 1.0],
    params,
  )
  let e_far = @src.omm_pme_real_space(
    [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
    [1.0, 1.0],
    params,
  )
  assert_true(e_near > 0.0)
  assert_true(e_far > 0.0)
  assert_true(e_near > e_far)
  // erfc weighted Coulomb must stay below bare Coulomb
  assert_true(e_near < 1.0)
}

///|
test "omm pme real space opposite charges attractive" {
  let params = @src.OmmPMEParams::new(8.0, box=[0.0, 0.0, 0.0])
  let e = @src.omm_pme_real_space(
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    [1.0, -1.0],
    params,
  )
  assert_true(e < 0.0)
}

///|
test "omm pme real space respects cutoff" {
  let params = @src.OmmPMEParams::new(8.0, box=[0.0, 0.0, 0.0])
  let e = @src.omm_pme_real_space(
    [[0.0, 0.0, 0.0], [9.0, 0.0, 0.0]],
    [1.0, 1.0],
    params,
  )
  assert_eq(e, 0.0)
}

///|
test "omm pme real space minimum image convention" {
  // separation 9 in a 10 A box wraps to 1 -> same as direct distance 1
  let params_boxed = @src.OmmPMEParams::new(8.0, box=[10.0, 10.0, 10.0])
  let e_wrapped = @src.omm_pme_real_space(
    [[0.0, 0.0, 0.0], [9.0, 0.0, 0.0]],
    [1.0, 1.0],
    params_boxed,
  )
  let params_open = @src.OmmPMEParams::new(8.0, box=[0.0, 0.0, 0.0])
  let e_direct = @src.omm_pme_real_space(
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    [1.0, 1.0],
    params_open,
  )
  assert_true(@moonbit.assert_num_close(e_wrapped, e_direct, 1.0e-12))
}

///|
test "omm pme reciprocal space zero box returns zero" {
  let params = @src.OmmPMEParams::new(8.0, box=[0.0, 0.0, 0.0])
  let e = @src.omm_pme_reciprocal_space(
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    [1.0, 1.0],
    params,
  )
  assert_eq(e, 0.0)
}

///|
test "omm pme reciprocal space positive for like charges" {
  let params = @src.OmmPMEParams::new(8.0, box=[10.0, 10.0, 10.0])
  let e = @src.omm_pme_reciprocal_space(
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
    [1.0, 1.0],
    params,
  )
  assert_true(e > 0.0)
}

///|
test "omm pme full energy equals sum of components" {
  let params = @src.OmmPMEParams::new(8.0, box=[10.0, 10.0, 10.0])
  let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 1.0, 0.5]]
  let charges = [0.5, -0.3, -0.2]
  let e_total = @src.omm_pme_energy(positions, charges, params)
  let e_real = @src.omm_pme_real_space(positions, charges, params)
  let e_recip = @src.omm_pme_reciprocal_space(positions, charges, params)
  let e_self = @src.omm_pme_self_energy(charges, params)
  assert_true(
    @moonbit.assert_num_close(e_total, e_real + e_recip + e_self, 1.0e-12),
  )
}
