///|
/// Tests for Biopython Bio.PDB.cealign-compatible structure alignment.

///|
fn ce_test_aligner(window_size : Int, max_gap : Int) -> @src.CEAligner {
  @src.CEAligner::new(window_size~, max_gap~) catch {
    _ => abort("valid CEAligner parameters should be accepted")
  }
}

///|
fn ce_test_curve(count : Int) -> Array[@src.Vector3] {
  let coordinates : Array[@src.Vector3] = []
  let mut index = 0
  while index < count {
    let value = index.to_double()
    coordinates.push(
      @src.Vector3::new(
        value * 1.3,
        value * value * 0.17 + (index % 3).to_double() * 0.4,
        value * 0.9 + (index % 4).to_double() * 0.6,
      ),
    )
    index = index + 1
  }
  coordinates
}

///|
fn ce_test_transform_coordinates(
  coordinates : Array[@src.Vector3],
  cosine : Double,
  sine : Double,
  translation : @src.Vector3,
) -> Array[@src.Vector3] {
  let result : Array[@src.Vector3] = []
  for coordinate in coordinates {
    result.push(
      @src.Vector3::new(
        cosine * coordinate.x - sine * coordinate.y + translation.x,
        sine * coordinate.x + cosine * coordinate.y + translation.y,
        coordinate.z + translation.z,
      ),
    )
  }
  result
}

///|
fn ce_test_structure(
  id : String,
  coordinates : Array[@src.Vector3],
  guide_name : String,
) -> @src.Structure {
  let residues : Array[@src.Residue] = []
  let mut index = 0
  while index < coordinates.length() {
    let guide = @src.Atom::new(
      name=guide_name,
      coord=coordinates[index],
      resname="GLY",
      chainid='A',
      resseq=index + 1,
      element=if guide_name == "CA" { "C" } else { "C" },
    )
    let companion = @src.Atom::new(
      name="N",
      coord=@src.Vector3::new(
        coordinates[index].x + 0.3,
        coordinates[index].y - 0.2,
        coordinates[index].z + 0.1,
      ),
      resname="GLY",
      chainid='A',
      resseq=index + 1,
      element="N",
    )
    residues.push(
      @src.Residue::new(resname="GLY", chainid='A', resseq=index + 1, atoms=[
        guide, companion,
      ]),
    )
    index = index + 1
  }
  @src.Structure::new(id~, models=[
    @src.Model::new(id=0, chains=[@src.Chain::new(id='A', residues~)]),
  ])
}

///|
fn ce_test_assert_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort("values differ beyond tolerance")
  }
}

///|
test "cealign: creates default and custom aligners" {
  let default_aligner = @src.CEAligner::new() catch {
    _ => abort("default CEAligner should be valid")
  }
  let custom = ce_test_aligner(4, 7)
  assert_eq(default_aligner.window_size, 8)
  assert_eq(default_aligner.max_gap, 30)
  assert_eq(default_aligner.reference_length(), 0)
  assert_eq(custom.window_size, 4)
  assert_eq(custom.max_gap, 7)
}

///|
test "cealign: rejects invalid parameters" {
  let invalid_window = try {
    ignore(@src.CEAligner::new(window_size=0))
    false
  } catch {
    CeAlignError(_) => true
  }
  let invalid_gap = try {
    ignore(@src.CEAligner::new(max_gap=-1))
    false
  } catch {
    CeAlignError(_) => true
  }
  assert_true(invalid_window)
  assert_true(invalid_gap)
}

///|
test "cealign: extracts protein CA guide atoms" {
  let structure = ce_test_structure("protein", ce_test_curve(4), "CA")
  let atoms = @src.cealign_get_guide_atoms(structure) catch {
    _ => abort("protein structure should provide guide atoms")
  }
  assert_eq(atoms.length(), 4)
  assert_eq(atoms[0].atom_name, "CA")
  assert_eq(atoms[0].residue_number, 1)
}

///|
test "cealign: extracts nucleic acid C4 prime guide atoms" {
  let structure = ce_test_structure("rna", ce_test_curve(4), "C4'")
  let atoms = @src.cealign_get_guide_atoms(structure) catch {
    _ => abort("nucleic acid structure should provide guide atoms")
  }
  assert_eq(atoms.length(), 4)
  assert_eq(atoms[0].atom_name, "C4'")
}

///|
test "cealign: prefers CA when a residue contains both guide atom types" {
  let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(1.0, 2.0, 3.0))
  let c4 = @src.Atom::new(name="C4'", coord=@src.Vector3::new(9.0, 9.0, 9.0))
  let residue = @src.Residue::new(resname="MIX", chainid='A', resseq=1, atoms=[
    c4, ca,
  ])
  let structure = @src.Structure::new(id="mixed", models=[
    @src.Model::new(id=0, chains=[@src.Chain::new(id='A', residues=[residue])]),
  ])
  let atoms = @src.cealign_get_guide_atoms(structure) catch {
    _ => abort("mixed residue should provide one guide atom")
  }
  assert_eq(atoms.length(), 1)
  assert_eq(atoms[0].atom_name, "CA")
  assert_eq(atoms[0].coordinate.x, 1.0)
}

///|
test "cealign: sorts guide atoms by model chain and residue" {
  let atom_a2 = @src.Atom::new(
    name="CA",
    coord=@src.Vector3::new(2.0, 0.0, 0.0),
  )
  let atom_a1 = @src.Atom::new(
    name="CA",
    coord=@src.Vector3::new(1.0, 0.0, 0.0),
  )
  let atom_b1 = @src.Atom::new(
    name="CA",
    coord=@src.Vector3::new(3.0, 0.0, 0.0),
  )
  let chain_b = @src.Chain::new(id='B', residues=[
    @src.Residue::new(resname="GLY", chainid='B', resseq=1, atoms=[atom_b1]),
  ])
  let chain_a = @src.Chain::new(id='A', residues=[
    @src.Residue::new(resname="GLY", chainid='A', resseq=2, atoms=[atom_a2]),
    @src.Residue::new(resname="GLY", chainid='A', resseq=1, atoms=[atom_a1]),
  ])
  let structure = @src.Structure::new(id="unordered", models=[
    @src.Model::new(id=0, chains=[chain_b, chain_a]),
  ])
  let atoms = @src.cealign_get_guide_atoms(structure) catch {
    _ => abort("guide atoms should be sortable")
  }
  assert_eq(atoms[0].chain_id, 'A')
  assert_eq(atoms[0].residue_number, 1)
  assert_eq(atoms[1].chain_id, 'A')
  assert_eq(atoms[1].residue_number, 2)
  assert_eq(atoms[2].chain_id, 'B')
}

///|
test "cealign: rejects structures without guide atoms" {
  let atom = @src.Atom::new(name="N", coord=@src.Vector3::new(0.0, 0.0, 0.0))
  let residue = @src.Residue::new(resname="GLY", chainid='A', resseq=1, atoms=[
    atom,
  ])
  let structure = @src.Structure::new(id="no-guide", models=[
    @src.Model::new(id=0, chains=[@src.Chain::new(id='A', residues=[residue])]),
  ])
  let raised = try {
    ignore(@src.cealign_get_guide_atoms(structure))
    false
  } catch {
    CeAlignError(_) => true
  }
  assert_true(raised)
}

///|
test "cealign: builds symmetric distance matrices" {
  let matrix = @src.cealign_distance_matrix([
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(3.0, 4.0, 0.0),
    @src.Vector3::new(3.0, 4.0, 12.0),
  ])
  assert_eq(matrix.length(), 3)
  assert_eq(matrix[0][0], 0.0)
  assert_eq(matrix[0][1], 5.0)
  assert_eq(matrix[1][0], 5.0)
  assert_eq(matrix[1][2], 12.0)
}

///|
test "cealign: identical fragments have zero dissimilarity" {
  let coordinates = ce_test_curve(5)
  let distances = @src.cealign_distance_matrix(coordinates)
  let similarity = @src.cealign_fragment_similarity(
    distances, distances, 0, 0, 4,
  )
  assert_eq(similarity, 0.0)
}

///|
test "cealign: distorted fragments have negative similarity" {
  let first = ce_test_curve(5)
  let second = first.copy()
  second[3] = @src.Vector3::new(20.0, -10.0, 5.0)
  let similarity = @src.cealign_fragment_similarity(
    @src.cealign_distance_matrix(first),
    @src.cealign_distance_matrix(second),
    0,
    0,
    4,
  )
  assert_true(similarity < 0.0)
}

///|
test "cealign: similarity matrix has all fragment starts" {
  let similarity = @src.cealign_similarity_matrix(
    ce_test_curve(7),
    ce_test_curve(6),
    3,
  ) catch {
    _ => abort("valid coordinate sets should produce a matrix")
  }
  assert_eq(similarity.length(), 5)
  assert_eq(similarity[0].length(), 4)
  assert_eq(similarity[0][0], 0.0)
}

///|
test "cealign: similarity matrix validates full windows" {
  let raised = try {
    ignore(
      @src.cealign_similarity_matrix(ce_test_curve(2), ce_test_curve(4), 3),
    )
    false
  } catch {
    CeAlignError(_) => true
  }
  assert_true(raised)
}

///|
test "cealign: identical coordinates produce a full CE path" {
  let coordinates = ce_test_curve(8)
  let paths = @src.cealign_find_paths(coordinates, coordinates, 2, 0) catch {
    _ => abort("identical coordinates should align")
  }
  assert_true(paths.length() > 0)
  assert_eq(paths[0].fragment_count, 4)
  assert_eq(paths[0].aligned_length(), 8)
  assert_eq(paths[0].reference_path(), [0, 1, 2, 3, 4, 5, 6, 7])
  assert_eq(paths[0].mobile_path(), [0, 1, 2, 3, 4, 5, 6, 7])
}

///|
test "cealign: path search returns at most twenty ranked candidates" {
  let coordinates = ce_test_curve(12)
  let paths = @src.cealign_find_paths(coordinates, coordinates, 2, 2) catch {
    _ => abort("valid coordinates should produce CE paths")
  }
  assert_true(paths.length() <= 20)
  let mut index = 1
  while index < paths.length() {
    assert_true(paths[index - 1].fragment_count >= paths[index].fragment_count)
    index = index + 1
  }
}

///|
test "cealign: CE paths are strictly monotonic" {
  let coordinates = ce_test_curve(10)
  let path = (@src.cealign_find_paths(coordinates, coordinates, 2, 1) catch {
      _ => abort("valid coordinates should produce CE paths")
    })[0]
  let mut index = 1
  while index < path.aligned_length() {
    assert_true(
      path.reference_indices[index] > path.reference_indices[index - 1],
    )
    assert_true(path.mobile_indices[index] > path.mobile_indices[index - 1])
    index = index + 1
  }
}

///|
test "cealign: path accessors return defensive copies" {
  let coordinates = ce_test_curve(8)
  let path = (@src.cealign_find_paths(coordinates, coordinates, 2, 0) catch {
      _ => abort("valid coordinates should produce CE paths")
    })[0]
  let reference_copy = path.reference_path()
  reference_copy[0] = 99
  let mobile_copy = path.mobile_path()
  mobile_copy[0] = 99
  assert_eq(path.reference_path()[0], 0)
  assert_eq(path.mobile_path()[0], 0)
}

///|
test "cealign: paths expose aligned index pairs" {
  let coordinates = ce_test_curve(6)
  let path = (@src.cealign_find_paths(coordinates, coordinates, 2, 0) catch {
      _ => abort("valid coordinates should produce CE paths")
    })[0]
  let pairs = path.aligned_pairs()
  assert_eq(pairs.length(), 6)
  assert_eq(pairs[0], (0, 0))
  assert_eq(pairs[5], (5, 5))
}

///|
test "cealign: path search validates parameters" {
  let invalid_gap = try {
    ignore(@src.cealign_find_paths(ce_test_curve(4), ce_test_curve(4), 2, -1))
    false
  } catch {
    CeAlignError(_) => true
  }
  let invalid_window = try {
    ignore(@src.cealign_find_paths(ce_test_curve(4), ce_test_curve(4), 0, 0))
    false
  } catch {
    CeAlignError(_) => true
  }
  assert_true(invalid_gap)
  assert_true(invalid_window)
}

///|
test "cealign: align requires a reference structure" {
  let aligner = ce_test_aligner(2, 2)
  let mobile = ce_test_structure("mobile", ce_test_curve(8), "CA")
  let raised = try {
    ignore(aligner.align(mobile, transform=false, final_optimization=false))
    false
  } catch {
    CeAlignError(_) => true
  }
  assert_true(raised)
}

///|
test "cealign: reference must contain two complete windows" {
  let aligner = ce_test_aligner(3, 2)
  let reference = ce_test_structure("short", ce_test_curve(5), "CA")
  let raised = try {
    ignore(aligner.set_reference(reference))
    false
  } catch {
    CeAlignError(_) => true
  }
  assert_true(raised)
}

///|
test "cealign: mobile structure must contain two complete windows" {
  let reference = ce_test_structure("reference", ce_test_curve(8), "CA")
  let mobile = ce_test_structure("short-mobile", ce_test_curve(3), "CA")
  let aligner = ce_test_aligner(2, 2).set_reference(reference) catch {
    _ => abort("reference should contain enough guide atoms")
  }
  let raised = try {
    ignore(aligner.align(mobile, transform=false, final_optimization=false))
    false
  } catch {
    CeAlignError(_) => true
  }
  assert_true(raised)
}

///|
test "cealign: set_reference retains ordered coordinates" {
  let coordinates = ce_test_curve(8)
  let reference = ce_test_structure("reference", coordinates, "CA")
  let aligner = ce_test_aligner(2, 2).set_reference(reference) catch {
    _ => abort("reference should be accepted")
  }
  let stored = aligner.reference_coordinates()
  assert_eq(aligner.reference_length(), 8)
  ce_test_assert_close(stored[0].x, coordinates[0].x, 1.0e-12)
  ce_test_assert_close(stored[7].z, coordinates[7].z, 1.0e-12)
}

///|
test "cealign: translated structures superimpose with low RMSD" {
  let reference_coordinates = ce_test_curve(12)
  let mobile_coordinates = ce_test_transform_coordinates(
    reference_coordinates,
    1.0,
    0.0,
    @src.Vector3::new(7.0, -4.0, 2.0),
  )
  let reference = ce_test_structure("reference", reference_coordinates, "CA")
  let mobile = ce_test_structure("mobile", mobile_coordinates, "CA")
  let aligner = ce_test_aligner(3, 2).set_reference(reference) catch {
    _ => abort("reference should be accepted")
  }
  let result = aligner.align(mobile, transform=false, final_optimization=false) catch {
    _ => abort("translated structure should align")
  }
  assert_true(result.rmsd < 1.0e-5)
  assert_eq(result.path.aligned_length(), 12)
}

///|
test "cealign: rotated and translated sample structures align" {
  let (reference, mobile) = @src.cealign_sample_structures()
  let aligner = @src.CEAligner::new() catch {
    _ => abort("default aligner should be valid")
  }
  let configured = aligner.set_reference(reference) catch {
    _ => abort("sample reference should be accepted")
  }
  let result = configured.align(
    mobile,
    transform=true,
    final_optimization=false,
  ) catch {
    _ => abort("sample structures should align")
  }
  assert_true(result.rmsd < 1.0e-4)
  assert_eq(result.path.aligned_length(), 32)
  assert_true(result.transformed)
}

///|
test "cealign: transform false preserves returned coordinates" {
  let coordinates = ce_test_curve(8)
  let mobile_coordinates = ce_test_transform_coordinates(
    coordinates,
    0.8,
    0.6,
    @src.Vector3::new(2.0, 3.0, -1.0),
  )
  let reference = ce_test_structure("reference", coordinates, "CA")
  let mobile = ce_test_structure("mobile", mobile_coordinates, "CA")
  let original = mobile.get_atoms()[0].coord
  let aligner = ce_test_aligner(2, 2).set_reference(reference) catch {
    _ => abort("reference should be accepted")
  }
  let result = aligner.align(mobile, transform=false, final_optimization=false) catch {
    _ => abort("mobile structure should align")
  }
  let returned = result.structure.get_atoms()[0].coord
  assert_false(result.transformed)
  ce_test_assert_close(returned.x, original.x, 1.0e-12)
  ce_test_assert_close(returned.y, original.y, 1.0e-12)
  ce_test_assert_close(returned.z, original.z, 1.0e-12)
}

///|
test "cealign: transform true moves guide and companion atoms" {
  let coordinates = ce_test_curve(10)
  let mobile_coordinates = ce_test_transform_coordinates(
    coordinates,
    0.8,
    0.6,
    @src.Vector3::new(4.0, -2.0, 3.0),
  )
  let reference = ce_test_structure("reference", coordinates, "CA")
  let mobile = ce_test_structure("mobile", mobile_coordinates, "CA")
  let aligner = ce_test_aligner(2, 2).set_reference(reference) catch {
    _ => abort("reference should be accepted")
  }
  let result = aligner.align(mobile, transform=true, final_optimization=false) catch {
    _ => abort("mobile structure should align")
  }
  let transformed_atoms = result.structure.get_atoms()
  assert_eq(transformed_atoms.length(), 20)
  assert_true(transformed_atoms[0].coord.distance(coordinates[0]) < 1.0e-4)
  assert_true(
    transformed_atoms[1].coord.distance(mobile.get_atoms()[1].coord) > 0.1,
  )
}

///|
test "cealign: alignment never mutates the input structure" {
  let (reference, mobile) = @src.cealign_sample_structures()
  let original = mobile.get_atoms()[0].coord
  let aligner = @src.CEAligner::new() catch {
    _ => abort("default aligner should be valid")
  }
  let configured = aligner.set_reference(reference) catch {
    _ => abort("sample reference should be accepted")
  }
  ignore(
    configured.align(mobile, transform=true, final_optimization=false) catch {
      _ => abort("sample structures should align")
    },
  )
  let after = mobile.get_atoms()[0].coord
  assert_eq(after.x, original.x)
  assert_eq(after.y, original.y)
  assert_eq(after.z, original.z)
}

///|
test "cealign: aligns nucleic acid structures through C4 prime atoms" {
  let coordinates = ce_test_curve(8)
  let mobile_coordinates = ce_test_transform_coordinates(
    coordinates,
    0.8,
    0.6,
    @src.Vector3::new(-2.0, 5.0, 1.0),
  )
  let reference = ce_test_structure("rna-ref", coordinates, "C4'")
  let mobile = ce_test_structure("rna-mobile", mobile_coordinates, "C4'")
  let aligner = ce_test_aligner(2, 2).set_reference(reference) catch {
    _ => abort("nucleic reference should be accepted")
  }
  let result = aligner.align(mobile, transform=false, final_optimization=false) catch {
    _ => abort("nucleic structures should align")
  }
  assert_true(result.rmsd < 1.0e-4)
  assert_eq(result.path.aligned_length(), 8)
}

///|
test "cealign: reports reference and mobile coverage" {
  let coordinates = ce_test_curve(10)
  let reference = ce_test_structure("reference", coordinates, "CA")
  let mobile = ce_test_structure("mobile", coordinates, "CA")
  let aligner = ce_test_aligner(2, 2).set_reference(reference) catch {
    _ => abort("reference should be accepted")
  }
  let result = aligner.align(mobile, transform=false, final_optimization=false) catch {
    _ => abort("identical structures should align")
  }
  assert_eq(result.reference_coverage(), 1.0)
  assert_eq(result.mobile_coverage(), 1.0)
}

///|
test "cealign: local path can bridge an insertion within max_gap" {
  let reference_coordinates = ce_test_curve(10)
  let mobile_coordinates : Array[@src.Vector3] = []
  let mut index = 0
  while index < reference_coordinates.length() {
    if index == 4 {
      mobile_coordinates.push(@src.Vector3::new(30.0, -20.0, 15.0))
    }
    mobile_coordinates.push(reference_coordinates[index])
    index = index + 1
  }
  let paths = @src.cealign_find_paths(
    reference_coordinates, mobile_coordinates, 2, 2,
  ) catch {
    _ => abort("CE should search paths across a short insertion")
  }
  assert_true(paths.length() > 0)
  assert_true(paths[0].aligned_length() >= 8)
  let pairs = paths[0].aligned_pairs()
  let mut bridges_insertion = false
  for pair in pairs {
    if pair.0 >= 4 && pair.1 == pair.0 + 1 {
      bridges_insertion = true
    }
  }
  assert_true(bridges_insertion)
}

///|
test "cealign: final optimization never increases RMSD" {
  let (reference, mobile) = @src.cealign_sample_structures()
  let configured = (@src.CEAligner::new() catch {
    _ => abort("default aligner should be valid")
  }).set_reference(reference) catch {
    _ => abort("sample reference should be accepted")
  }
  let baseline = configured.align(
    mobile,
    transform=false,
    final_optimization=false,
  ) catch {
    _ => abort("baseline alignment should succeed")
  }
  let optimized = configured.align(
    mobile,
    transform=false,
    final_optimization=true,
  ) catch {
    _ => abort("optimized alignment should succeed")
  }
  assert_true(optimized.rmsd <= baseline.rmsd + 1.0e-12)
  assert_true(optimized.path.z_score >= 3.5)
  assert_true(optimized.optimization_applied)
}

///|
test "cealign: result transform accessors are defensive" {
  let coordinates = ce_test_curve(8)
  let reference = ce_test_structure("reference", coordinates, "CA")
  let mobile = ce_test_structure("mobile", coordinates, "CA")
  let configured = ce_test_aligner(2, 2).set_reference(reference) catch {
    _ => abort("reference should be accepted")
  }
  let result = configured.align(
    mobile,
    transform=false,
    final_optimization=false,
  ) catch {
    _ => abort("identical structures should align")
  }
  let rotation = result.rotation_matrix()
  let translation = result.translation_vector()
  rotation[0][0] = 99.0
  translation[0] = 99.0
  assert_true(result.rotation_matrix()[0][0] < 2.0)
  assert_true(result.translation_vector()[0].abs() < 1.0e-6)
}

///|
test "cealign: validates explicit rigid transforms" {
  let structure = ce_test_structure("protein", ce_test_curve(4), "CA")
  let bad_rotation = try {
    ignore(
      @src.cealign_transform_structure(structure, [[1.0, 0.0], [0.0, 1.0]], [
        0.0, 0.0, 0.0,
      ]),
    )
    false
  } catch {
    CeAlignError(_) => true
  }
  let bad_translation = try {
    ignore(
      @src.cealign_transform_structure(
        structure,
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [0.0, 0.0],
      ),
    )
    false
  } catch {
    CeAlignError(_) => true
  }
  assert_true(bad_rotation)
  assert_true(bad_translation)
}

///|
test "cealign: identity transform preserves atom coordinates" {
  let structure = ce_test_structure("protein", ce_test_curve(4), "CA")
  let transformed = @src.cealign_transform_structure(
    structure,
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    [0.0, 0.0, 0.0],
  ) catch {
    _ => abort("identity transform should be valid")
  }
  let original_coordinate = structure.get_atoms()[3].coord
  let transformed_coordinate = transformed.get_atoms()[3].coord
  assert_eq(transformed_coordinate.x, original_coordinate.x)
  assert_eq(transformed_coordinate.y, original_coordinate.y)
  assert_eq(transformed_coordinate.z, original_coordinate.z)
}

///|
test "cealign: summary includes core alignment metrics" {
  let coordinates = ce_test_curve(8)
  let reference = ce_test_structure("reference", coordinates, "CA")
  let mobile = ce_test_structure("mobile", coordinates, "CA")
  let configured = ce_test_aligner(2, 2).set_reference(reference) catch {
    _ => abort("reference should be accepted")
  }
  let result = configured.align(
    mobile,
    transform=false,
    final_optimization=false,
  ) catch {
    _ => abort("identical structures should align")
  }
  let summary = result.summary()
  assert_true(summary.has_prefix("CEAlignment("))
  assert_true(summary.contains("aligned=8"))
  assert_true(summary.contains("rmsd="))
  assert_true(summary.contains("reference_coverage=1"))
}
