///|
/// PDB Packing Density Analysis Demo.
///
/// Demonstrates protein packing density analysis including:
/// - Packing density calculation
/// - Per-residue packing analysis
/// - Solvent Accessible Surface Area (SASA)
/// - Low packing region identification
/// - Normalization and comparison
///
/// Run: moon run examples/pdb_packing_demo/main.mbt

///|
fn main {
  println("=== PDB Packing Density Analysis Demo ===\n")

  // === 1. Create demo atoms ===
  println("1. Creating Demo Packing Atoms:")
  println("================================")
  let atoms = @src.create_demo_packing_atoms()
  println("  Total atoms: " + atoms.length().to_string())

  let mut i = 0
  while i < atoms.length() {
    let atom = atoms[i]
    println(
      "  [" +
      i.to_string() +
      "] " +
      atom.atom_name +
      " (" +
      atom.residue_name +
      atom.residue_seq.to_string() +
      ") at (" +
      atom.x.to_string() +
      ", " +
      atom.y.to_string() +
      ", " +
      atom.z.to_string() +
      ") R=" +
      atom.vdw_radius.to_string(),
    )
    i = i + 1
  }

  // === 2. Generate sphere points ===
  println("\n2. Sphere Points Generation:")
  println("============================")
  let point_counts = [10, 50, 100]
  let mut pc_i = 0
  while pc_i < point_counts.length() {
    let n = point_counts[pc_i]
    let points = @src.generate_packing_sphere_points(n_points=n)
    println("  Generated " + points.length().to_string() + " points on unit sphere")
    let first = points[0]
    let dist = (first.x * first.x + first.y * first.y + first.z * first.z).sqrt()
    println("    First point distance from origin: " + dist.to_string())
    pc_i = pc_i + 1
  }

  // === 3. Test atom_in_sphere ===
  println("\n3. Atom in Sphere Tests:")
  println("========================")
  let test_center = [1.0, 2.0, 3.0]
  let test_point_inside = [1.5, 2.5, 3.5]
  let test_point_outside = [10.0, 10.0, 10.0]
  println("  Point inside sphere (r=2.0): " +
    @src.atom_in_sphere(test_center, test_point_inside, 2.0).to_string())
  println("  Point outside sphere (r=2.0): " +
    @src.atom_in_sphere(test_center, test_point_outside, 2.0).to_string())

  // === 4. Calculate packing density ===
  println("\n4. Local Packing Density Calculation:")
  println("=======================================")
  let densities = [
    @src.calculate_packing_density(atoms, radius=3.0, n_points=20),
    @src.calculate_packing_density(atoms, radius=5.0, n_points=50),
    @src.calculate_packing_density(atoms, radius=7.0, n_points=100),
  ]
  let radii = [3.0, 5.0, 7.0]
  let mut d_i = 0
  while d_i < densities.length() {
    let d = densities[d_i]
    println(
      "  Radius=" +
      radii[d_i].to_string() +
      " A: density=" +
      d.density.to_string() +
      ", shell points=" +
      d.n_shell_points.to_string(),
    )
    d_i = d_i + 1
  }

  // === 5. Per-residue packing analysis ===
  println("\n5. Per-Residue Packing Analysis:")
  println("==================================")
  let analysis = @src.packing_density_per_residue(atoms, radius=5.0, n_points=50)
  println("  Number of residues analyzed: " + analysis.n_residues.to_string())
  println("  Mean density: " + analysis.mean_density.to_string())
  println("  Min density: " + analysis.min_density.to_string())
  println("  Max density: " + analysis.max_density.to_string())
  println("  Median density: " + analysis.median_density.to_string())
  println("  Buried residues: " + analysis.n_buried.to_string())
  println("  Exposed residues: " + analysis.n_exposed.to_string())

  let results = analysis.get_results()
  let mut r_i = 0
  while r_i < results.length() {
    let r = results[r_i]
    println(
      "  Residue " +
      r.residue_seq.to_string() +
      " (" +
      r.residue_name +
      "): density=" +
      r.density.to_string() +
      ", contact atoms=" +
      r.n_contacting_atoms.to_string(),
    )
    r_i = r_i + 1
  }

  // === 6. Solvent Accessible Surface Area ===
  println("\n6. Solvent Accessible Surface Area (SASA):")
  println("==========================================")
  let probe_radii = [0.5, 1.0, 1.4, 2.0]
  let mut pr_i = 0
  while pr_i < probe_radii.length() {
    let pr = probe_radii[pr_i]
    let sasa = @src.calculate_packing_sasa(atoms, probe_radius=pr, n_points=50)
    println("  Probe radius=" + pr.to_string() + " A: SASA=" + sasa.to_string() + " A^2")
    pr_i = pr_i + 1
  }

  // === 7. Normalize packing density ===
  println("\n7. Density Normalization:")
  println("==========================")
  let original = @src.calculate_packing_density(atoms, radius=5.0, n_points=50)
  let thresholds = [0.5, 1.0, 2.0, 5.0]
  let mut th_i = 0
  while th_i < thresholds.length() {
    let t = thresholds[th_i]
    let normalized = @src.normalize_packing_density(original, threshold=t)
    println(
      "  Original=" +
      original.density.to_string() +
      ", threshold=" +
      t.to_string() +
      " -> normalized=" +
      normalized.density.to_string(),
    )
    th_i = th_i + 1
  }

  // === 8. Identify low packing regions ===
  println("\n8. Low Packing Region Identification:")
  println("=====================================")
  let cutoffs = [0.0, 0.5, 1.0, 2.0]
  let mut co_i = 0
  while co_i < cutoffs.length() {
    let c = cutoffs[co_i]
    let low = @src.identify_low_packing(analysis, cutoff=c)
    println("  Cutoff=" + c.to_string() + ": " + low.length().to_string() + " low-density residues")
    for id in low {
      println("    - " + id)
    }
    co_i = co_i + 1
  }

  // === 9. Packing efficiency ===
  println("\n9. Packing Efficiency:")
  println("=====================")
  let container_radii = [5.0, 10.0, 20.0, 50.0]
  let mut cr_i = 0
  while cr_i < container_radii.length() {
    let cr = container_radii[cr_i]
    let efficiency = @src.calculate_packing_efficiency(atoms, container_radius=cr)
    println("  Container radius=" + cr.to_string() + " A: efficiency=" + efficiency.to_string())
    cr_i = cr_i + 1
  }

  // === 10. Sphere geometry ===
  println("\n10. Sphere Geometry:")
  println("====================")
  let test_radii = [1.0, 2.0, 3.0, 5.0]
  let mut sr_i = 0
  while sr_i < test_radii.length() {
    let r = test_radii[sr_i]
    let vol = @src.sphere_volume(r)
    let recovered_r = @src.sphere_radius_from_volume(vol)
    println(
      "  R=" +
      r.to_string() +
      " A: volume=" +
      vol.to_string() +
      " A^3, recovered R=" +
      recovered_r.to_string(),
    )
    sr_i = sr_i + 1
  }

  // === 11. VDW radii lookup ===
  println("\n11. Van der Waals Radii:")
  println("=========================")
  let atom_types = ["H", "C", "CA", "CB", "N", "NZ", "O", "OG", "S", "SG", "P"]
  let mut at_i = 0
  while at_i < atom_types.length() {
    let at = atom_types[at_i]
    let r = @src.get_packing_vdw_radius(at)
    println("  " + at + ": " + r.to_string() + " A")
    at_i = at_i + 1
  }

  // === 12. Summary ===
  println("\n12. Analysis Summary:")
  println("=====================")
  println(analysis.summary())

  // === 13. Distance calculations ===
  println("\n13. Distance Calculations:")
  println("==========================")
  let atom_a = @src.PackingAtom::new(x=0.0, y=0.0, z=0.0)
  let atom_b = @src.PackingAtom::new(x=3.0, y=4.0, z=0.0)
  println("  Distance between (0,0,0) and (3,4,0): " + atom_a.distance(atom_b).to_string())
  let atom_c = @src.PackingAtom::new(x=1.0, y=1.0, z=1.0)
  println("  Distance between (0,0,0) and (1,1,1): " + atom_a.distance(atom_c).to_string())

  println("\n=== PDB Packing Density Analysis Demo Complete ===")
}