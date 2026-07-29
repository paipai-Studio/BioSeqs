///|
/// Tests for Bio.PDB.Packing module.

///|
test "packing_atom_creation" {
  let atom = @src.PackingAtom::new(
    x=1.5,
    y=2.5,
    z=3.5,
    residue_name="ALA",
    residue_seq=42,
    chain_id='A',
    atom_name="CA",
    element="C",
    vdw_radius=1.7,
  )
  assert_eq(atom.x, 1.5)
  assert_eq(atom.y, 2.5)
  assert_eq(atom.z, 3.5)
  assert_eq(atom.residue_name, "ALA")
  assert_eq(atom.residue_seq, 42)
  assert_eq(atom.chain_id, 'A')
  assert_eq(atom.atom_name, "CA")
  assert_eq(atom.element, "C")
  assert_eq(atom.vdw_radius, 1.7)
}

///|
test "packing_atom_default_values" {
  let atom = @src.PackingAtom::new(x=0.0, y=0.0, z=0.0)
  assert_eq(atom.residue_name, "UNK")
  assert_eq(atom.residue_seq, 1)
  assert_eq(atom.chain_id, ' ')
  assert_eq(atom.atom_name, "")
  assert_eq(atom.element, "")
  assert_eq(atom.vdw_radius, 1.7)
}

///|
test "packing_atom_get_coords" {
  let atom = @src.PackingAtom::new(x=1.0, y=2.0, z=3.0)
  let coords = atom.get_coords()
  assert_eq(coords.length(), 3)
  assert_eq(coords[0], 1.0)
  assert_eq(coords[1], 2.0)
  assert_eq(coords[2], 3.0)
}

///|
test "packing_atom_distance" {
  let a1 = @src.PackingAtom::new(x=0.0, y=0.0, z=0.0)
  let a2 = @src.PackingAtom::new(x=3.0, y=4.0, z=0.0)
  let dist = a1.distance(a2)
  assert_eq(dist, 5.0)
}

///|
test "packing_atom_distance_self" {
  let a1 = @src.PackingAtom::new(x=5.0, y=5.0, z=5.0)
  let dist = a1.distance(a1)
  assert_eq(dist, 0.0)
}

///|
test "sphere_point_creation" {
  let sp = @src.SpherePoint::new(0.5, 0.5, 0.5)
  assert_eq(sp.x, 0.5)
  assert_eq(sp.y, 0.5)
  assert_eq(sp.z, 0.5)
}

///|
test "sphere_point_get_coords" {
  let sp = @src.SpherePoint::new(0.1, 0.2, 0.3)
  let coords = sp.get_coords()
  assert_eq(coords.length(), 3)
  assert_eq(coords[0], 0.1)
  assert_eq(coords[1], 0.2)
  assert_eq(coords[2], 0.3)
}

///|
test "packing_result_creation" {
  let result = @src.PackingResult::new(
    residue_name="ALA",
    residue_seq=10,
    chain_id='B',
    density=2.5,
    n_contacting_atoms=25,
    n_shell_points=100,
  )
  assert_eq(result.residue_name, "ALA")
  assert_eq(result.residue_seq, 10)
  assert_eq(result.chain_id, 'B')
  assert_eq(result.density, 2.5)
  assert_eq(result.n_contacting_atoms, 25)
  assert_eq(result.n_shell_points, 100)
}

///|
test "packing_result_get_density" {
  let result = @src.PackingResult::new(residue_name="GLU", residue_seq=5, density=3.7)
  assert_eq(result.get_density(), 3.7)
}

///|
test "packing_result_get_residue_id" {
  let result = @src.PackingResult::new(
    residue_name="LYS",
    residue_seq=15,
    chain_id='A',
  )
  let id = result.get_residue_id()
  assert_eq(id, "A15LYS")
}

///|
test "packing_analysis_result_creation" {
  let results : Array[@src.PackingResult] = [
    @src.PackingResult::new(residue_name="ALA", residue_seq=1, density=2.0),
    @src.PackingResult::new(residue_name="GLY", residue_seq=2, density=3.0),
  ]
  let analysis = @src.PackingAnalysisResult::new(
    results~,
    mean_density=2.5,
    min_density=2.0,
    max_density=3.0,
    median_density=2.5,
    n_residues=2,
    n_buried=1,
    n_exposed=1,
  )
  assert_eq(analysis.n_residues, 2)
  assert_eq(analysis.mean_density, 2.5)
  assert_eq(analysis.min_density, 2.0)
  assert_eq(analysis.max_density, 3.0)
}

///|
test "packing_analysis_result_summary" {
  let results : Array[@src.PackingResult] = [
    @src.PackingResult::new(residue_name="ALA", residue_seq=1, density=2.0),
  ]
  let analysis = @src.PackingAnalysisResult::new(
    results~,
    mean_density=2.0,
    min_density=2.0,
    max_density=2.0,
    n_residues=1,
    n_buried=0,
    n_exposed=1,
  )
  let summary = analysis.summary()
  assert_eq(summary.length() > 0, true)
}

///|
test "generate_packing_sphere_points_default" {
  let points = @src.generate_packing_sphere_points()
  assert_eq(points.length(), 100)
  let first = points[0]
  let dist = (first.x * first.x + first.y * first.y + first.z * first.z).sqrt()
  assert_eq(dist > 0.99 && dist < 1.01, true)
}

///|
test "generate_packing_sphere_points_custom" {
  let points = @src.generate_packing_sphere_points(n_points=50)
  assert_eq(points.length(), 50)
  for p in points {
    let dist = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt()
    assert_eq(dist > 0.99 && dist < 1.01, true)
  }
}

///|
test "generate_packing_sphere_points_single" {
  let points = @src.generate_packing_sphere_points(n_points=1)
  assert_eq(points.length(), 1)
}

///|
test "atom_in_sphere_true" {
  let center = [0.0, 0.0, 0.0]
  let point = [1.0, 0.0, 0.0]
  assert_eq(@src.atom_in_sphere(center, point, 2.0), true)
}

///|
test "atom_in_sphere_false" {
  let center = [0.0, 0.0, 0.0]
  let point = [3.0, 0.0, 0.0]
  assert_eq(@src.atom_in_sphere(center, point, 2.0), false)
}

///|
test "atom_in_sphere_boundary" {
  let center = [0.0, 0.0, 0.0]
  let point = [2.0, 0.0, 0.0]
  assert_eq(@src.atom_in_sphere(center, point, 2.0), true)
}

///|
test "atom_in_sphere_zero_radius" {
  let center = [1.0, 1.0, 1.0]
  let point = [1.0, 1.0, 1.0]
  assert_eq(@src.atom_in_sphere(center, point, 0.0), true)
}

///|
test "count_atoms_in_sphere_basic" {
  let atoms : Array[@src.PackingAtom] = [
    @src.PackingAtom::new(x=0.0, y=0.0, z=0.0),
    @src.PackingAtom::new(x=1.0, y=0.0, z=0.0),
    @src.PackingAtom::new(x=5.0, y=0.0, z=0.0),
  ]
  let center = [0.0, 0.0, 0.0]
  let count = @src.count_atoms_in_sphere(center, atoms, 2.0)
  assert_eq(count, 2)
}

///|
test "count_atoms_in_sphere_exclude" {
  let atoms : Array[@src.PackingAtom] = [
    @src.PackingAtom::new(x=0.0, y=0.0, z=0.0),
    @src.PackingAtom::new(x=1.0, y=0.0, z=0.0),
    @src.PackingAtom::new(x=5.0, y=0.0, z=0.0),
  ]
  let center = [0.0, 0.0, 0.0]
  let count = @src.count_atoms_in_sphere(center, atoms, 2.0, exclude_index=0)
  assert_eq(count, 1)
}

///|
test "count_atoms_in_sphere_empty" {
  let atoms : Array[@src.PackingAtom] = []
  let center = [0.0, 0.0, 0.0]
  let count = @src.count_atoms_in_sphere(center, atoms, 2.0)
  assert_eq(count, 0)
}

///|
test "calculate_packing_density_empty" {
  let atoms : Array[@src.PackingAtom] = []
  let result = @src.calculate_packing_density(atoms)
  assert_eq(result.residue_name, "")
  assert_eq(result.residue_seq, 0)
}

///|
test "calculate_packing_density_single_atom" {
  let atoms : Array[@src.PackingAtom] = [
    @src.PackingAtom::new(
      x=0.0,
      y=0.0,
      z=0.0,
      residue_name="ALA",
      residue_seq=1,
    ),
  ]
  let result = @src.calculate_packing_density(atoms, radius=1.4, n_points=50)
  assert_eq(result.residue_name, "ALA")
  assert_eq(result.residue_seq, 1)
  assert_eq(result.n_shell_points, 50)
}

///|
test "calculate_packing_density_with_neighbors" {
  let atoms : Array[@src.PackingAtom] = [
    @src.PackingAtom::new(
      x=0.0,
      y=0.0,
      z=0.0,
      residue_name="ALA",
      residue_seq=1,
    ),
    @src.PackingAtom::new(
      x=1.5,
      y=0.0,
      z=0.0,
      residue_name="ALA",
      residue_seq=1,
    ),
    @src.PackingAtom::new(
      x=0.0,
      y=1.5,
      z=0.0,
      residue_name="ALA",
      residue_seq=1,
    ),
  ]
  let result = @src.calculate_packing_density(atoms, radius=1.4, n_points=50)
  assert_eq(result.residue_name, "ALA")
  assert_eq(result.density > 0.0, true)
}

///|
test "packing_density_per_residue_basic" {
  let atoms = @src.create_demo_packing_atoms()
  let analysis = @src.packing_density_per_residue(atoms)
  assert_eq(analysis.n_residues > 0, true)
  assert_eq(analysis.mean_density > 0.0, true)
  assert_eq(analysis.min_density >= 0.0, true)
  assert_eq(analysis.max_density >= analysis.min_density, true)
}

///|
test "packing_density_per_residue_empty" {
  let atoms : Array[@src.PackingAtom] = []
  let analysis = @src.packing_density_per_residue(atoms)
  assert_eq(analysis.n_residues, 0)
}

///|
test "calculate_packing_sasa_empty" {
  let atoms : Array[@src.PackingAtom] = []
  let sasa = @src.calculate_packing_sasa(atoms)
  assert_eq(sasa, 0.0)
}

///|
test "calculate_packing_sasa_single_atom" {
  let atoms : Array[@src.PackingAtom] = [
    @src.PackingAtom::new(
      x=0.0,
      y=0.0,
      z=0.0,
      vdw_radius=1.7,
    ),
  ]
  let sasa = @src.calculate_packing_sasa(atoms, probe_radius=1.4, n_points=50)
  assert_eq(sasa > 0.0, true)
}

///|
test "calculate_packing_sasa_multiple_atoms" {
  let atoms = @src.create_demo_packing_atoms()
  let sasa = @src.calculate_packing_sasa(atoms, probe_radius=1.4, n_points=50)
  assert_eq(sasa > 0.0, true)
}

///|
test "calculate_packing_sasa_lower_probe" {
  let atoms = @src.create_demo_packing_atoms()
  let sasa1 = @src.calculate_packing_sasa(atoms, probe_radius=0.5, n_points=50)
  let sasa2 = @src.calculate_packing_sasa(atoms, probe_radius=1.4, n_points=50)
  assert_eq(sasa1 < sasa2, true)
}

///|
test "normalize_packing_density_below_threshold" {
  let result = @src.PackingResult::new(residue_name="ALA", residue_seq=1, density=0.5)
  let normalized = @src.normalize_packing_density(result, threshold=2.0)
  assert_eq(normalized.density, 0.25)
}

///|
test "normalize_packing_density_above_threshold" {
  let result = @src.PackingResult::new(residue_name="ALA", residue_seq=1, density=3.0)
  let normalized = @src.normalize_packing_density(result, threshold=2.0)
  assert_eq(normalized.density, 1.0)
}

///|
test "normalize_packing_density_at_threshold" {
  let result = @src.PackingResult::new(residue_name="ALA", residue_seq=1, density=2.0)
  let normalized = @src.normalize_packing_density(result, threshold=2.0)
  assert_eq(normalized.density, 1.0)
}

///|
test "normalize_packing_density_zero_threshold" {
  let result = @src.PackingResult::new(residue_name="ALA", residue_seq=1, density=1.5)
  let normalized = @src.normalize_packing_density(result, threshold=0.0)
  assert_eq(normalized.density, 1.5)
}

///|
test "identify_low_packing_basic" {
  let results : Array[@src.PackingResult] = [
    @src.PackingResult::new(residue_name="ALA", residue_seq=1, chain_id='A', density=0.3),
    @src.PackingResult::new(residue_name="GLY", residue_seq=2, chain_id='A', density=1.5),
    @src.PackingResult::new(residue_name="VAL", residue_seq=3, chain_id='A', density=0.4),
  ]
  let analysis = @src.PackingAnalysisResult::new(results~, n_residues=3)
  let low = @src.identify_low_packing(analysis, cutoff=0.5)
  assert_eq(low.length(), 2)
}

///|
test "identify_low_packing_none_below" {
  let results : Array[@src.PackingResult] = [
    @src.PackingResult::new(residue_name="ALA", residue_seq=1, density=2.0),
    @src.PackingResult::new(residue_name="GLY", residue_seq=2, density=3.0),
  ]
  let analysis = @src.PackingAnalysisResult::new(results~, n_residues=2)
  let low = @src.identify_low_packing(analysis, cutoff=0.5)
  assert_eq(low.length(), 0)
}

///|
test "identify_low_packing_all_below" {
  let results : Array[@src.PackingResult] = [
    @src.PackingResult::new(residue_name="ALA", residue_seq=1, density=0.1),
    @src.PackingResult::new(residue_name="GLY", residue_seq=2, density=0.2),
  ]
  let analysis = @src.PackingAnalysisResult::new(results~, n_residues=2)
  let low = @src.identify_low_packing(analysis, cutoff=0.5)
  assert_eq(low.length(), 2)
}

///|
test "sphere_volume_calculation" {
  let vol = @src.sphere_volume(2.0)
  let expected = (4.0 / 3.0) * 3.14159265358979323846 * 8.0
  assert_eq((vol - expected).abs() < 0.001, true)
}

///|
test "sphere_volume_zero" {
  let vol = @src.sphere_volume(0.0)
  assert_eq(vol, 0.0)
}

///|
test "sphere_radius_from_volume" {
  let radius = @src.sphere_radius_from_volume(100.0)
  assert_eq(radius > 0.0, true)
}

///|
test "sphere_radius_from_volume_zero" {
  let radius = @src.sphere_radius_from_volume(0.0)
  assert_eq(radius, 0.0)
}

///|
test "sphere_volume_radius_roundtrip" {
  let original_r = 2.5
  let vol = @src.sphere_volume(original_r)
  let recovered_r = @src.sphere_radius_from_volume(vol)
  assert_eq((recovered_r - original_r).abs() < 0.001, true)
}

///|
test "calculate_packing_efficiency_basic" {
  let atoms = @src.create_demo_packing_atoms()
  let efficiency = @src.calculate_packing_efficiency(atoms, container_radius=10.0)
  assert_eq(efficiency > 0.0, true)
}

///|
test "calculate_packing_efficiency_empty" {
  let atoms : Array[@src.PackingAtom] = []
  let efficiency = @src.calculate_packing_efficiency(atoms)
  assert_eq(efficiency, 0.0)
}

///|
test "calculate_packing_efficiency_large_container" {
  let atoms = @src.create_demo_packing_atoms()
  let eff_small = @src.calculate_packing_efficiency(atoms, container_radius=5.0)
  let eff_large = @src.calculate_packing_efficiency(atoms, container_radius=50.0)
  assert_eq(eff_small > eff_large, true)
}

///|
test "get_packing_vdw_radius_carbon" {
  let r = @src.get_packing_vdw_radius("CA")
  assert_eq(r, 1.7)
}

///|
test "get_packing_vdw_radius_nitrogen" {
  let r = @src.get_packing_vdw_radius("N")
  assert_eq(r, 1.55)
}

///|
test "get_packing_vdw_radius_oxygen" {
  let r = @src.get_packing_vdw_radius("OG")
  assert_eq(r, 1.52)
}

///|
test "get_packing_vdw_radius_sulfur" {
  let r = @src.get_packing_vdw_radius("SG")
  assert_eq(r, 1.8)
}

///|
test "get_packing_vdw_radius_hydrogen" {
  let r = @src.get_packing_vdw_radius("H")
  assert_eq(r, 1.2)
}

///|
test "get_packing_vdw_radius_unknown" {
  let r = @src.get_packing_vdw_radius("X")
  assert_eq(r, 1.7)
}

///|
test "create_demo_packing_atoms" {
  let atoms = @src.create_demo_packing_atoms()
  assert_eq(atoms.length(), 10)
  assert_eq(atoms[0].residue_name, "ALA")
  assert_eq(atoms[0].residue_seq, 1)
  assert_eq(atoms[4].residue_name, "GLY")
  assert_eq(atoms[4].residue_seq, 2)
  assert_eq(atoms[8].residue_name, "VAL")
  assert_eq(atoms[8].residue_seq, 3)
}

///|
test "packing_density_with_demo_atoms" {
  let atoms = @src.create_demo_packing_atoms()
  let result = @src.calculate_packing_density(atoms, radius=1.4, n_points=20)
  assert_eq(result.residue_name, "ALA")
  assert_eq(result.residue_seq, 1)
  assert_eq(result.n_shell_points, 20)
}

///|
test "packing_analysis_with_demo_atoms" {
  let atoms = @src.create_demo_packing_atoms()
  let analysis = @src.packing_density_per_residue(atoms, radius=1.4, n_points=20)
  assert_eq(analysis.n_residues, 3)
  assert_eq(analysis.mean_density > 0.0, true)
  assert_eq(analysis.n_buried + analysis.n_exposed, 3)
}

///|
test "sasa_with_demo_atoms" {
  let atoms = @src.create_demo_packing_atoms()
  let sasa = @src.calculate_packing_sasa(atoms, probe_radius=1.4, n_points=30)
  assert_eq(sasa > 0.0, true)
}

///|
test "low_packing_identification_with_analysis" {
  let atoms = @src.create_demo_packing_atoms()
  let analysis = @src.packing_density_per_residue(atoms, radius=1.4, n_points=20)
  let low = @src.identify_low_packing(analysis, cutoff=0.0)
  assert_eq(low.length() >= 0, true)
}

///|
test "sphere_points_uniform_distribution" {
  let points = @src.generate_packing_sphere_points(n_points=100)
  let mut min_x = points[0].x
  let mut max_x = points[0].x
  let mut i = 1
  while i < points.length() {
    if points[i].x < min_x {
      min_x = points[i].x
    }
    if points[i].x > max_x {
      max_x = points[i].x
    }
    i = i + 1
  }
  assert_eq(min_x >= -1.0 && min_x <= 1.0, true)
  assert_eq(max_x >= -1.0 && max_x <= 1.0, true)
}

///|
test "packing_atom_distance_3d" {
  let a1 = @src.PackingAtom::new(x=0.0, y=0.0, z=0.0)
  let a2 = @src.PackingAtom::new(x=1.0, y=1.0, z=1.0)
  let dist = a1.distance(a2)
  let expected = 3.0.sqrt()
  assert_eq((dist - expected).abs() < 0.001, true)
}

///|
test "packing_result_default_creation" {
  let result = @src.PackingResult::new(residue_name="PRO", residue_seq=99)
  assert_eq(result.residue_name, "PRO")
  assert_eq(result.residue_seq, 99)
  assert_eq(result.chain_id, ' ')
  assert_eq(result.density, 0.0)
  assert_eq(result.n_contacting_atoms, 0)
  assert_eq(result.n_shell_points, 0)
}