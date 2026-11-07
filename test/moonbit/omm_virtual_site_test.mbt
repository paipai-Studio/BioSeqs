///|
test "omm_tip4p_virtual_site_position" {
  // general check of the M-site weights (-0.8, 0.9, 0.9):
  // M = 0.9*H1 + 0.9*H2 - 0.8*O
  let pos = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
  let m = @src.omm_tip4p_virtual_site(0, 1, 2, 3)
  assert_eq(m.index(), 3)
  let p = @src.omm_virtual_site_position(m, pos)
  assert_true(@moonbit.assert_num_close(p[0], 0.9, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[1], 0.9, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[2], 0.0, 1.0e-12))
  // symmetric water: M stays on the H-O-H bisector
  let half = 52.26 * 3.141592653589793 / 180.0
  let hx = 0.9572 * @math.sin(half)
  let hz = 0.9572 * @math.cos(half)
  let water = [[0.0, 0.0, 0.0], [hx, 0.0, hz], [-hx, 0.0, hz]]
  let m2 = @src.omm_tip4p_virtual_site(0, 1, 2, 3)
  let p2 = @src.omm_virtual_site_position(m2, water)
  assert_true(@moonbit.assert_num_close(p2[0], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(p2[1], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(p2[2], 1.8 * hz, 1.0e-9))
}

///|
test "omm_drude_virtual_site_tracks_parent" {
  let pos = [[1.5, -2.0, 3.0], [0.0, 0.0, 0.0]]
  let d = @src.omm_drude_virtual_site(0, 1)
  assert_eq(d.index(), 1)
  let p = @src.omm_virtual_site_position(d, pos)
  // TwoParticleAverage(parent, parent, 1.0, 0.0) collapses onto the parent
  assert_true(@moonbit.assert_num_close(p[0], 1.5, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[1], -2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[2], 3.0, 1.0e-12))
}

///|
test "omm_update_virtual_sites_batch" {
  // 4 particles: O, H1, H2 and a placeholder for the M site
  let half = 52.26 * 3.141592653589793 / 180.0
  let hx = 0.9572 * @math.sin(half)
  let hz = 0.9572 * @math.cos(half)
  let pos = [[0.0, 0.0, 0.0], [hx, 0.0, hz], [-hx, 0.0, hz], [9.0, 9.0, 9.0]]
  let m = @src.omm_tip4p_virtual_site(0, 1, 2, 3)
  let out = @src.omm_update_virtual_sites(pos, [m])
  assert_eq(out.length(), 4)
  // parent atoms untouched
  assert_true(@moonbit.assert_num_close(out[0][0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(out[1][0], hx, 1.0e-12))
  assert_true(@moonbit.assert_num_close(out[2][0], -hx, 1.0e-12))
  // placeholder replaced by the M position on the bisector
  assert_true(@moonbit.assert_num_close(out[3][0], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(out[3][2], 1.8 * hz, 1.0e-9))
  // input array not mutated
  assert_true(@moonbit.assert_num_close(pos[3][0], 9.0, 1.0e-12))
}
