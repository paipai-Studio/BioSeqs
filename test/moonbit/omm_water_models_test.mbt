///|
test "omm_water_spc_spce_parameters" {
  let spc = @src.omm_water_spc()
  assert_eq(spc.name, "SPC")
  assert_eq(spc.n_sites, 3)
  assert_true(@moonbit.assert_num_close(spc.dOH, 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(spc.dHH, 1.63299, 1.0e-12))
  assert_true(@moonbit.assert_num_close(spc.angle_HOH, 109.47, 1.0e-12))
  assert_true(@moonbit.assert_num_close(spc.qO, -0.82, 1.0e-12))
  assert_true(@moonbit.assert_num_close(spc.qH, 0.41, 1.0e-12))
  assert_true(@moonbit.assert_num_close(spc.sigma_O, 3.166, 1.0e-12))
  assert_true(@moonbit.assert_num_close(spc.epsilon_O, 0.1553, 1.0e-12))
  assert_true(@moonbit.assert_num_close(spc.sigma_H, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(spc.mass_O, 15.9994, 1.0e-12))
  assert_true(@moonbit.assert_num_close(spc.mass_H, 1.008, 1.0e-12))
  assert_false(spc.has_virtual_site)
  let spce = @src.omm_water_spce()
  assert_eq(spce.name, "SPC/E")
  assert_eq(spce.n_sites, 3)
  assert_true(@moonbit.assert_num_close(spce.qO, -0.8476, 1.0e-12))
  assert_true(@moonbit.assert_num_close(spce.qH, 0.4238, 1.0e-12))
  assert_true(@moonbit.assert_num_close(spce.sigma_O, 3.166, 1.0e-12))
}

///|
test "omm_water_4site_5site_parameters" {
  let ew = @src.omm_water_tip4p_ew()
  assert_eq(ew.name, "TIP4P-Ew")
  assert_eq(ew.n_sites, 4)
  assert_true(ew.has_virtual_site)
  assert_true(@moonbit.assert_num_close(ew.qO, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(ew.qH, 0.52422, 1.0e-12))
  assert_true(@moonbit.assert_num_close(ew.qV, -1.04844, 1.0e-12))
  assert_true(@moonbit.assert_num_close(ew.dOV, 0.125, 1.0e-12))
  assert_true(@moonbit.assert_num_close(ew.dOH, 0.9572, 1.0e-12))
  assert_true(@moonbit.assert_num_close(ew.sigma_O, 3.16435, 1.0e-12))
  assert_true(@moonbit.assert_num_close(ew.epsilon_O, 0.16275, 1.0e-12))
  let t5 = @src.omm_water_tip5p()
  assert_eq(t5.name, "TIP5P")
  assert_eq(t5.n_sites, 5)
  assert_true(t5.has_virtual_site)
  assert_true(@moonbit.assert_num_close(t5.qH, 0.241, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t5.qV, -0.241, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t5.dOV, 0.70, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t5.sigma_O, 3.097, 1.0e-12))
  assert_true(@moonbit.assert_num_close(t5.epsilon_O, 0.0989, 1.0e-12))
  let fb = @src.omm_water_tip3p_fb()
  assert_eq(fb.name, "TIP3P-FB")
  assert_eq(fb.n_sites, 3)
  assert_false(fb.has_virtual_site)
  assert_true(@moonbit.assert_num_close(fb.qO, -0.410, 1.0e-12))
  assert_true(@moonbit.assert_num_close(fb.qH, 0.205, 1.0e-12))
  assert_true(@moonbit.assert_num_close(fb.sigma_O, 3.18195, 1.0e-12))
  assert_true(@moonbit.assert_num_close(fb.epsilon_O, 0.0901, 1.0e-12))
}

///|
test "omm_water_model_by_name_lookup" {
  match @src.omm_water_model_by_name("spc") {
    Some(m) => assert_eq(m.name, "SPC")
    None => fail("spc should resolve")
  }
  match @src.omm_water_model_by_name("TIP4P-EW") {
    Some(m) => assert_eq(m.name, "TIP4P-Ew")
    None => fail("tip4p-ew should resolve")
  }
  match @src.omm_water_model_by_name("spce") {
    Some(m) => assert_eq(m.name, "SPC/E")
    None => fail("spce should resolve")
  }
  match @src.omm_water_model_by_name("SPC/E") {
    Some(m) => assert_eq(m.name, "SPC/E")
    None => fail("spc/e should resolve")
  }
  match @src.omm_water_model_by_name("tip3p-fb") {
    Some(m) => assert_eq(m.name, "TIP3P-FB")
    None => fail("tip3p-fb should resolve")
  }
  match @src.omm_water_model_by_name("mystery") {
    Some(_) => fail("unknown model must not resolve")
    None => ()
  }
  let names = @src.omm_water_model_names()
  assert_eq(names.length(), 5)
  assert_eq(names[0], "SPC")
  assert_eq(names[4], "TIP5P")
}

///|
test "omm_water_coords_charges_masses" {
  let spc = @src.omm_water_spc()
  let coords = @src.omm_water_coords(spc)
  assert_eq(coords.length(), 3)
  let (ox, oy, oz) = coords[0]
  assert_true(@moonbit.assert_num_close(ox, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(oy, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(oz, 0.0, 1.0e-12))
  // O-H = 1.0 with the 109.47 degree angle in the xz plane
  let (h1x, h1y, h1z) = coords[1]
  assert_true(@moonbit.assert_num_close(h1x, 0.816497, 1.0e-5))
  assert_true(@moonbit.assert_num_close(h1y, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(h1z, 0.577350, 1.0e-5))
  // H2 mirrors H1
  let (h2x, _, h2z) = coords[2]
  assert_true(@moonbit.assert_num_close(h2x, -h1x, 1.0e-12))
  assert_true(@moonbit.assert_num_close(h2z, h1z, 1.0e-12))
  // H-H distance matches dHH
  let dhh = @math.pow(
    (h2x - h1x) * (h2x - h1x) + (h2z - h1z) * (h2z - h1z),
    0.5,
  )
  assert_true(@moonbit.assert_num_close(dhh, 1.63299, 1.0e-5))
  let q3 = @src.omm_water_charges(spc)
  assert_eq(q3.length(), 3)
  assert_true(@moonbit.assert_num_close(q3[0], -0.82, 1.0e-12))
  assert_true(@moonbit.assert_num_close(q3[1], 0.41, 1.0e-12))
  let m3 = @src.omm_water_masses(spc)
  assert_eq(m3.length(), 3)
  assert_true(@moonbit.assert_num_close(m3[0], 15.9994, 1.0e-12))
  assert_true(@moonbit.assert_num_close(m3[1], 1.008, 1.0e-12))
  // TIP4P-Ew adds a massless, negatively charged M site at (0,0,-0.125)
  let ew = @src.omm_water_tip4p_ew()
  let c4 = @src.omm_water_coords(ew)
  assert_eq(c4.length(), 4)
  let (_, _, mz) = c4[3]
  assert_true(@moonbit.assert_num_close(mz, -0.125, 1.0e-12))
  let q4 = @src.omm_water_charges(ew)
  assert_eq(q4.length(), 4)
  assert_true(@moonbit.assert_num_close(q4[3], -1.04844, 1.0e-12))
  let m4 = @src.omm_water_masses(ew)
  assert_eq(m4.length(), 4)
  assert_true(@moonbit.assert_num_close(m4[3], 0.0, 1.0e-12))
  let mut total_q = 0.0
  for qi in q4 {
    total_q = total_q + qi
  }
  assert_true(@moonbit.assert_num_close(total_q, 0.0, 1.0e-12))
  // TIP5P adds two lone pairs at (0,0,-0.70)
  let t5 = @src.omm_water_tip5p()
  let c5 = @src.omm_water_coords(t5)
  assert_eq(c5.length(), 5)
  let (_, _, lz) = c5[4]
  assert_true(@moonbit.assert_num_close(lz, -0.70, 1.0e-12))
  let q5 = @src.omm_water_charges(t5)
  assert_eq(q5.length(), 5)
  assert_true(@moonbit.assert_num_close(q5[3], -0.241, 1.0e-12))
  assert_true(@moonbit.assert_num_close(q5[4], -0.241, 1.0e-12))
}
