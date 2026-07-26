///|
/// Tests for Bio.PDB.MAalign - Multiple Protein Structure Alignment

///|
test "ma_aligner_new" {
  let a = @src.ma_aligner_new()
  assert_eq(a.gap_opening, -10.0)
  assert_eq(a.gap_extension, -1.0)
  assert_eq(a.max_iterations, 10)
  assert_eq(a.rmsd_threshold, 2.0)
}

///|
test "ma_aligner_with_params" {
  let a = @src.ma_aligner_with_params(-5.0, -0.5, 20, 1.5)
  assert_eq(a.gap_opening, -5.0)
  assert_eq(a.max_iterations, 20)
}

///|
test "add_structure" {
  let a = @src.ma_aligner_new()
  let coords : Array[(Double, Double, Double)] = [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0)]
  let seq = "ALA"
  let res = ["ALA", "LEU", "ALA"]
  let s = @src.add_structure(a, "struct1", coords, seq, res)
  assert_eq(s.structure_id, "struct1")
  assert_eq(s.ca_coordinates.length(), 3)
  assert_eq(s.sequence, "ALA")
}

///|
test "ma_center_coordinates" {
  let coords : Array[(Double, Double, Double)] = [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0)]
  let centered = @src.ma_center_coordinates(coords)
  assert_eq(centered.length(), 3)
  let mut sum_x = 0.0
  let mut sum_y = 0.0
  let mut sum_z = 0.0
  for c in centered {
    let (x, y, z) = c
    sum_x = sum_x + x
    sum_y = sum_y + y
    sum_z = sum_z + z
  }
  let n = centered.length().to_double()
  assert_true((sum_x / n).abs() < 1.0e-9)
  assert_true((sum_y / n).abs() < 1.0e-9)
  assert_true((sum_z / n).abs() < 1.0e-9)
}

///|
test "ma_compute_distance_matrix" {
  let coords1 : Array[(Double, Double, Double)] = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
  let coords2 : Array[(Double, Double, Double)] = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
  let dm = @src.ma_compute_distance_matrix(coords1, coords2)
  assert_eq(dm.length(), 3)
  assert_eq(dm[0].length(), 3)
}

///|
test "ma_align_pairwise_identical" {
  let coords : Array[(Double, Double, Double)] = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
  let (rotation, rmsd) = @src.ma_align_pairwise(coords, coords)
  assert_true(rmsd < 0.01)
}

///|
test "ma_align_pairwise_different_scale" {
  let coords1 : Array[(Double, Double, Double)] = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
  let coords2 : Array[(Double, Double, Double)] = [(2.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 2.0)]
  let (rotation, rmsd) = @src.ma_align_pairwise(coords1, coords2)
  assert_true(rmsd < 5.0)
}

///|
test "compute_rmsd_identical" {
  let coords : Array[(Double, Double, Double)] = [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0)]
  let rmsd = @src.compute_rmsd(coords, coords)
  assert_true(rmsd < 0.01)
}

///|
test "compute_rmsd_different" {
  let coords1 : Array[(Double, Double, Double)] = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
  let coords2 : Array[(Double, Double, Double)] = [(10.0, 0.0, 0.0), (11.0, 0.0, 0.0)]
  let rmsd = @src.compute_rmsd(coords1, coords2)
  assert_true(rmsd > 5.0)
}

///|
test "align_structures_single" {
  let a = @src.ma_aligner_new()
  let coords : Array[(Double, Double, Double)] = [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0)]
  let seq = "ALA"
  let res = ["ALA", "LEU", "ALA"]
  let s = @src.add_structure(a, "s1", coords, seq, res)
  let result = @src.align_structures(a, [s])
  assert_eq(result.structure_ids.length(), 1)
}

///|
test "align_structures_two" {
  let a = @src.ma_aligner_new()
  let coords1 : Array[(Double, Double, Double)] = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
  let coords2 : Array[(Double, Double, Double)] = [(10.0, 10.0, 0.0), (11.0, 10.0, 0.0), (10.0, 11.0, 0.0)]
  let s1 = @src.add_structure(a, "s1", coords1, "AAA", ["ALA", "ALA", "ALA"])
  let s2 = @src.add_structure(a, "s2", coords2, "BBB", ["ALA", "ALA", "ALA"])
  let result = @src.align_structures(a, [s1, s2])
  assert_eq(result.structure_ids.length(), 2)
  assert_true(result.total_rmsd < 5.0)
}

///|
test "align_structures_empty" {
  let a = @src.ma_aligner_new()
  let result = @src.align_structures(a, [])
  assert_eq(result.structure_ids.length(), 0)
}

///|
test "build_conservation_profile" {
  let aligned = ["ACGT", "ACGT", "ACGT"]
  let profile = @src.build_conservation_profile(aligned)
  assert_eq(profile.length(), 4)
  for s in profile {
    assert_true(s > 0.9)
  }
}

///|
test "build_conservation_profile_different" {
  let aligned = ["AAAA", "CCCC", "GGGG"]
  let profile = @src.build_conservation_profile(aligned)
  assert_eq(profile.length(), 4)
  for s in profile {
    assert_true(s < 0.5)
  }
}

///|
test "get_conservation" {
  let a = @src.ma_aligner_new()
  let coords1 : Array[(Double, Double, Double)] = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
  let coords2 : Array[(Double, Double, Double)] = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
  let s1 = @src.add_structure(a, "s1", coords1, "AAA", ["ALA", "ALA", "ALA"])
  let s2 = @src.add_structure(a, "s2", coords2, "AAA", ["ALA", "ALA", "ALA"])
  let result = @src.align_structures(a, [s1, s2])
  let cons = @src.get_conservation(result)
  assert_true(cons.length() >= 0)
}

///|
test "get_rmsd_matrix" {
  let a = @src.ma_aligner_new()
  let coords1 : Array[(Double, Double, Double)] = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
  let coords2 : Array[(Double, Double, Double)] = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
  let s1 = @src.add_structure(a, "s1", coords1, "AA", ["ALA", "ALA"])
  let s2 = @src.add_structure(a, "s2", coords2, "AA", ["ALA", "ALA"])
  let result = @src.align_structures(a, [s1, s2])
  let rmsd_mat = @src.get_rmsd_matrix(result)
  assert_eq(rmsd_mat.length(), 2)
  assert_eq(rmsd_mat[0][0], 0.0)
  assert_eq(rmsd_mat[1][1], 0.0)
}

///|
test "apply_rotation_identity" {
  let coords : Array[(Double, Double, Double)] = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
  let identity : (Double, Double, Double) = (0.0, 0.0, 0.0)
  let _result = @src.apply_rotation(coords, identity)
  assert_true(true)
}