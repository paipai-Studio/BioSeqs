///|
test "omm_parse_tinker_params_sections" {
  let content = "atom\n1 C 12.011 -0.5 1\n2 H 1.008 0.25 1\nbond\n1 2 5.0 1.09\nangle\n1 2 2 3.0 109.5\ntorsion\n1 2 3 4 1.5 30.0 3\nvdw\n1 3.5 0.1\n"
  let p = @src.omm_parse_tinker_params(content)
  assert_eq(p.atom_types.length(), 2)
  let (tid, el, mass, chg) = p.atom_types[0]
  assert_eq(tid, 1)
  assert_eq(el, "C")
  assert_true(@moonbit.assert_num_close(mass, 12.011, 1.0e-12))
  assert_true(@moonbit.assert_num_close(chg, -0.5, 1.0e-12))
  let (_, el1, mass1, chg1) = p.atom_types[1]
  assert_eq(el1, "H")
  assert_true(@moonbit.assert_num_close(mass1, 1.008, 1.0e-12))
  assert_true(@moonbit.assert_num_close(chg1, 0.25, 1.0e-12))
  assert_eq(p.bond_params.length(), 1)
  let (ti, tj, k, r0) = p.bond_params[0]
  assert_eq(ti, 1)
  assert_eq(tj, 2)
  assert_true(@moonbit.assert_num_close(k, 5.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(r0, 1.09, 1.0e-12))
  assert_eq(p.angle_params.length(), 1)
  let (_, _, _, ak, ath) = p.angle_params[0]
  assert_true(@moonbit.assert_num_close(ak, 3.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(ath, 109.5, 1.0e-12))
  assert_eq(p.torsion_params.length(), 1)
  let (_, _, _, _, tk, tph, tn) = p.torsion_params[0]
  assert_true(@moonbit.assert_num_close(tk, 1.5, 1.0e-12))
  assert_true(@moonbit.assert_num_close(tph, 30.0, 1.0e-12))
  assert_eq(tn, 3)
  assert_eq(p.vdw_params.length(), 1)
  let (vt, vr, ve) = p.vdw_params[0]
  assert_eq(vt, 1)
  assert_true(@moonbit.assert_num_close(vr, 3.5, 1.0e-12))
  assert_true(@moonbit.assert_num_close(ve, 0.1, 1.0e-12))
}

///|
test "omm_parse_tinker_xyz_water" {
  let content = "3\nwater\n1 O 0.0 0.0 0.0 1 2 3\n2 H 0.9572 0.0 0.0 2 1\n3 H -0.9572 0.0 0.0 2 1\n"
  let xyz = @src.omm_parse_tinker_xyz(content)
  assert_eq(xyz.n_atoms, 3)
  assert_eq(xyz.title, "water")
  assert_eq(xyz.atoms.length(), 3)
  let (el, x, y, z, tid) = xyz.atoms[0]
  assert_eq(el, "O")
  assert_true(@moonbit.assert_num_close(x, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(y, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(z, 0.0, 1.0e-12))
  assert_eq(tid, 1)
  let (_, hx, _, _, _) = xyz.atoms[1]
  assert_true(@moonbit.assert_num_close(hx, 0.9572, 1.0e-12))
  let (_, h2x, _, _, _) = xyz.atoms[2]
  assert_true(@moonbit.assert_num_close(h2x, -0.9572, 1.0e-12))
  // connectivity list converted to 0-based indices
  // oxygen line carries bonded partners "2 3" -> indices 1 and 2
  assert_eq(xyz.bonds[0].length(), 2)
  assert_eq(xyz.bonds[0][0], 1)
  assert_eq(xyz.bonds[0][1], 2)
  assert_eq(xyz.bonds[1].length(), 1)
  assert_eq(xyz.bonds[1][0], 0)
  assert_eq(xyz.bonds[2].length(), 1)
  assert_eq(xyz.bonds[2][0], 0)
}

///|
test "omm_tinker_bond_energy" {
  let params = @src.omm_parse_tinker_params("bond\n1 2 5.0 1.0\n")
  let xyz = @src.omm_parse_tinker_xyz(
    "2\ntest\n1 O 0.0 0.0 0.0 1 2\n2 H 2.0 0.0 0.0 2 1\n",
  )
  // implementation squares the raw distance component sum: r = 4.0, dr = 3.0
  let e = @src.omm_tinker_bond_energy(params, xyz)
  assert_true(@moonbit.assert_num_close(e, 5.0 * 9.0, 1.0e-12))
  // equilibrium geometry gives zero energy regardless of convention
  let xyz0 = @src.omm_parse_tinker_xyz(
    "2\ntest\n1 O 0.0 0.0 0.0 1 2\n2 H 1.0 0.0 0.0 2 1\n",
  )
  let e0 = @src.omm_tinker_bond_energy(params, xyz0)
  assert_true(@moonbit.assert_num_close(e0, 0.0, 1.0e-12))
  // no matching bond parameters -> zero energy
  let empty = @src.omm_parse_tinker_params("bond\n3 4 5.0 1.0\n")
  let e2 = @src.omm_tinker_bond_energy(empty, xyz)
  assert_true(@moonbit.assert_num_close(e2, 0.0, 1.0e-12))
}
