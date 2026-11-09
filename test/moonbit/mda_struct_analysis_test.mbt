///|
/// Tests for mda_struct_analysis.mbt (ANM, DSSP-like assignment, Voronoi,
/// bond geometry, protein-membrane analysis).

///|
test "mda_anm_two_atoms_along_x" {
  // Two residues 1 A apart: g = -1/d^2 = -1, so h[0][0] = h[3][3] = 1
  // and h[0][3] = h[3][0] = -1 (xx block); all other blocks vanish.
  let coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
  let r = @src.mda_anm(coords)
  assert_eq(r.n_modes, 6)
  assert_eq(r.hessian.length(), 6)
  assert_true(@moonbit.assert_num_close(r.hessian[0][0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r.hessian[3][3], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r.hessian[0][3], -1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r.hessian[3][0], -1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r.hessian[0][4], 0.0, 1.0e-12))
  // Eigenvalues are the Hessian diagonal entries.
  assert_true(@moonbit.assert_num_close(r.eigenvalues[0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r.eigenvalues[3], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r.eigenvalues[1], 0.0, 1.0e-12))
  // Fluctuations are 1/ev for nonzero entries.
  assert_true(@moonbit.assert_num_close(r.mode_fluctuations[0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r.mode_fluctuations[1], 0.0, 1.0e-12))
}

///|
test "mda_anm_cutoff_excludes_distant_pairs" {
  // Pair separated by 15 A exceeds the 12 A cutoff: hessian stays zero.
  let coords = [(0.0, 0.0, 0.0), (15.0, 0.0, 0.0)]
  let r = @src.mda_anm(coords)
  assert_eq(r.n_modes, 6)
  assert_true(@moonbit.assert_num_close(r.hessian[0][0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(r.hessian[0][3], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(r.eigenvalues[0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(r.mode_fluctuations[0], 0.0, 1.0e-12))
}

///|
test "mda_dssp_assign_alpha_helix" {
  // Ideal alpha helix: rise 1.5 A, 100 deg per residue, radius 2.3.
  // CA(i)-CA(i+3) distance ~5.05 A lies in (4.5, 6) -> all H.
  let omega = 1.7453292519943295
  let coords : Array[(Double, Double, Double)] = Array(capacity=8)
  for k = 0; k < 8; k = k + 1 {
    let t = omega * k.to_double()
    coords.push((2.3 * @math.cos(t), 2.3 * @math.sin(t), 1.5 * k.to_double()))
  }
  let res_names = ["ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA", "ALA"]
  let ss = @src.mda_dssp_assign(coords, res_names)
  assert_eq(ss.length(), 8)
  for s in ss {
    assert_eq(s, "H")
  }
}

///|
test "mda_dssp_assign_extended_strand" {
  // Four residues on a 3.8 x 6.8 grid: CA(i)-CA(i+2) = 6.8 in (6.5, 8)
  // and CA(i)-CA(i+3) = sqrt(60.68) ~ 7.79 outside the helix window
  // -> all residues marked E.
  let coords = [
    (0.0, 0.0, 0.0),
    (3.8, 0.0, 0.0),
    (0.0, 6.8, 0.0),
    (3.8, 6.8, 0.0),
  ]
  let ss = @src.mda_dssp_assign(coords, ["GLY", "GLY", "GLY", "GLY"])
  assert_eq(ss.length(), 4)
  for s in ss {
    assert_eq(s, "E")
  }
}

///|
test "mda_dssp_assign_turn_and_coil" {
  // 2 x 4 grid: CA(i)-CA(i+2) = 4 < 5.5 -> turn at i and i+1 only.
  let coords = [
    (0.0, 0.0, 0.0),
    (2.0, 0.0, 0.0),
    (0.0, 4.0, 0.0),
    (2.0, 4.0, 0.0),
  ]
  let ss = @src.mda_dssp_assign(coords, ["GLY", "GLY", "GLY", "GLY"])
  assert_eq(ss[0], "T")
  assert_eq(ss[1], "T")
  assert_eq(ss[2], "C")
  assert_eq(ss[3], "C")
}

///|
test "mda_dssp_assign_short_chain" {
  // Fewer than 4 residues: no helix/strand/turn detection.
  let ss = @src.mda_dssp_assign([(0.0, 0.0, 0.0), (3.8, 0.0, 0.0)], [
    "GLY", "GLY",
  ])
  assert_eq(ss.length(), 2)
  assert_eq(ss[0], "C")
  assert_eq(ss[1], "C")
}

///|
test "mda_voronoi_volumes_linear_chain" {
  // Atoms at x = 0, 1, 3 (cutoff 5): each has 2 neighbours, so the
  // volume is min_dist^3 / 6 -> 1/6, 1/6, 8/6.
  let vols = @src.mda_voronoi_volumes([
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (3.0, 0.0, 0.0),
  ])
  assert_eq(vols.length(), 3)
  assert_true(@moonbit.assert_num_close(vols[0], 1.0 / 6.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(vols[1], 1.0 / 6.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(vols[2], 8.0 / 6.0, 1.0e-9))
}

///|
test "mda_voronoi_volumes_isolated_atom" {
  // A single atom has no neighbours -> volume 0.
  let vols = @src.mda_voronoi_volumes([(1.0, 2.0, 3.0)])
  assert_eq(vols.length(), 1)
  assert_true(@moonbit.assert_num_close(vols[0], 0.0, 1.0e-12))
}

///|
test "mda_bond_length_stats_two_frames" {
  // Bond (0,1) is 1 A in frame 0 and 2 A in frame 1.
  let stats = @src.mda_bond_length_stats(
    [[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)], [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)]],
    [(0, 1)],
  )
  assert_eq(stats.n_samples, 2)
  assert_true(@moonbit.assert_num_close(stats.mean, 1.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(stats.std, 0.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(stats.min_, 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(stats.max_, 2.0, 1.0e-9))
}

///|
test "mda_bond_length_stats_empty" {
  let stats = @src.mda_bond_length_stats([], [])
  assert_eq(stats.n_samples, 0)
  assert_true(@moonbit.assert_num_close(stats.mean, 0.0, 1.0e-12))
}

///|
test "mda_bond_angle_stats_right_angle_to_zero" {
  // Angle (0,1,2): 90 deg in frame 0, 0 deg in frame 1.
  let frames = [
    [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
    [(2.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
  ]
  let stats = @src.mda_bond_angle_stats(frames, [(0, 1, 2)])
  assert_eq(stats.n_samples, 2)
  assert_true(@moonbit.assert_num_close(stats.mean, 45.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(stats.min_, 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(stats.max_, 90.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(stats.std, 45.0, 1.0e-9))
}

///|
test "mda_protein_membrane_transmembrane_span" {
  // CA z = -20..20 at x=0 with a 30 A membrane centered at 0:
  // |z| <= 15 selects residues at z = -10, 0, 10.
  let coords = [
    (0.0, 0.0, -20.0),
    (0.0, 0.0, -10.0),
    (0.0, 0.0, 0.0),
    (0.0, 0.0, 10.0),
    (0.0, 0.0, 20.0),
  ]
  let res_ids = [1, 2, 3, 4, 5]
  let r = @src.mda_protein_membrane(coords, res_ids)
  assert_eq(r.n_transmembrane, 3)
  assert_eq(r.transmembrane_residues[0], 2)
  assert_eq(r.transmembrane_residues[1], 3)
  assert_eq(r.transmembrane_residues[2], 4)
  assert_true(@moonbit.assert_num_close(r.insertion_depth, 0.0, 1.0e-9))
  // x range is 0 -> tilt angle atan(0/40) = 0.
  assert_true(@moonbit.assert_num_close(r.tilt_angle, 0.0, 1.0e-9))
}

///|
test "mda_protein_membrane_empty" {
  let r = @src.mda_protein_membrane([], [])
  assert_eq(r.n_transmembrane, 0)
  assert_true(@moonbit.assert_num_close(r.insertion_depth, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(r.tilt_angle, 0.0, 1.0e-12))
}
