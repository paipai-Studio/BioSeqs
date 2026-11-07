///|
/// Tests for rdkit_3d_align.mbt — PMI/NPR, PBF, USR, GetBestRMS and
/// grid shape alignment.

///|
test "pmi_mass_weighted_linear_chain" {
  // 5 carbons along the x axis: smallest moment about x
  let atoms = ["C", "C", "C", "C", "C"]
  let coords = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (2.0, 0.0, 0.0),
    (3.0, 0.0, 0.0),
    (4.0, 0.0, 0.0),
  ]
  let conf = @src.RdkitConformer::from_coords(atoms, coords)
  let (i1, i2, i3) = @src.rdkit_pmi_mass_weighted(conf)
  // rod along x: I_y = I_z >> I_x = 0
  assert_eq(i1, 0.0)
  assert_eq(i2, i3)
  let (npr1, _npr2, _npr3) = @src.rdkit_npr_ratios(conf)
  assert_eq(npr1, 0.0)
}

///|
test "pmi_mass_weighted_symmetric_square" {
  // square of carbons in the xy plane: I_z > I_x = I_y
  let atoms = ["C", "C", "C", "C"]
  let coords = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
  ]
  let conf = @src.RdkitConformer::from_coords(atoms, coords)
  let (i1, i2, i3) = @src.rdkit_pmi_mass_weighted(conf)
  assert_eq(i1, i2)
  assert_true(i3 > i1)
}

///|
test "pbf_planar_vs_out_of_plane" {
  // planar set -> ~0; same set displaced along z -> larger
  let atoms = ["C", "C", "C", "C"]
  let planar = @src.RdkitConformer::from_coords(atoms, [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
  ])
  let buckled = @src.RdkitConformer::from_coords(atoms, [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.6),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, -0.5),
  ])
  assert_true(@src.rdkit_pbf(planar) < 0.01)
  assert_true(@src.rdkit_pbf(buckled) > @src.rdkit_pbf(planar))
}

///|
test "usr_descriptors_self_distance_zero" {
  let atoms = ["C", "C", "O"]
  let coords = [(0.0, 0.0, 0.0), (1.5, 0.3, 0.1), (2.6, -0.4, 0.9)]
  let conf = @src.RdkitConformer::from_coords(atoms, coords)
  let d1 = @src.rdkit_usr_descriptors(conf)
  assert_eq(d1.length(), 12)
  assert_eq(@src.rdkit_usr_distance(d1, d1), 0.0)
  // a different geometry has a positive distance
  let coords2 = [(0.0, 0.0, 0.0), (1.5, 0.3, 0.1), (2.6, -0.4, 3.9)]
  let conf2 = @src.RdkitConformer::from_coords(atoms, coords2)
  assert_true(
    @src.rdkit_usr_distance(d1, @src.rdkit_usr_descriptors(conf2)) > 0.1,
  )
}

///|
test "get_best_rms_handles_symmetry_permutation" {
  // square of carbons: probe is a relabeling of the same geometry
  let atoms = ["C", "C", "C", "C"]
  let reference = @src.RdkitConformer::from_coords(atoms, [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
  ])
  // probe lists the same 4 points in a rotated labeling
  let probe = @src.RdkitConformer::from_coords(atoms, [
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 0.0),
  ])
  match @src.rdkit_get_best_rms(probe, reference) {
    Some((rms, _moved)) => assert_true(rms < 0.000001)
    None => fail("same-molecule alignment should succeed")
  }
}

///|
test "get_best_rms_invariant_under_rigid_motion" {
  let atoms = ["C", "C", "O"]
  let reference = @src.RdkitConformer::from_coords(atoms, [
    (0.0, 0.0, 0.0),
    (1.5, 0.0, 0.0),
    (2.4, 0.9, 0.2),
  ])
  // probe = reference translated by (10, -3, 7)
  let probe = @src.RdkitConformer::from_coords(atoms, [
    (10.0, -3.0, 7.0),
    (11.5, -3.0, 7.0),
    (12.4, -2.1, 7.2),
  ])
  match @src.rdkit_get_best_rms(probe, reference) {
    Some((rms, _moved)) => assert_true(rms < 0.000001)
    None => fail("rigid translation should align perfectly")
  }
}

///|
test "get_best_rms_rejects_different_composition" {
  let reference = @src.RdkitConformer::from_coords(["C", "C", "O"], [
    (0.0, 0.0, 0.0),
    (1.5, 0.0, 0.0),
    (2.4, 0.9, 0.2),
  ])
  let probe = @src.RdkitConformer::from_coords(["C", "C", "N"], [
    (0.0, 0.0, 0.0),
    (1.5, 0.0, 0.0),
    (2.4, 0.9, 0.2),
  ])
  assert_true(@src.rdkit_get_best_rms(probe, reference) is None)
}

///|
test "shape_align_self_overlap_high" {
  let conf = @src.rdkit_generate_conformer("CCO", seed=13)
  let result = @src.rdkit_shape_align(conf, conf)
  assert_true(result.shape_tanimoto > 0.8)
  // self volume ratio is 1
  assert_true(result.protrude > 0.9 && result.protrude < 1.1)
}
