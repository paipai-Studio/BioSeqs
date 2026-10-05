///|
/// Tests for CaviFinder (cavity/channel detection).

///|
/// Build a residue with atoms at given positions.
fn cav_make_res(resname : String, resseq : Int, atoms : Array[(@src.Vector3, String)]) -> @src.Residue {
  let atom_arr = Array(capacity=atoms.length())
  for (coord, name) in atoms {
    atom_arr.push(@src.Atom::new(name=name, coord=coord, resname=resname, resseq=resseq))
  }
  @src.Residue::new(resname=resname, chainid='A', resseq=resseq, atoms=atom_arr)
}

///|
test "cavifinder: empty chain" {
  let chain = @src.Chain::new(id='A', residues=[])
  let result = @src.find_cavities(chain)
  assert_eq(result.cavities.length(), 0)
}

///|
test "cavifinder: solid protein no cavity" {
  // A compact cluster of atoms — no internal cavity
  let atoms = [
    (@src.Vector3::new(0.0, 0.0, 0.0), "CA"),
    (@src.Vector3::new(3.0, 0.0, 0.0), "CA"),
    (@src.Vector3::new(0.0, 3.0, 0.0), "CA"),
    (@src.Vector3::new(0.0, 0.0, 3.0), "CA"),
    (@src.Vector3::new(3.0, 3.0, 0.0), "CA"),
    (@src.Vector3::new(3.0, 0.0, 3.0), "CA"),
    (@src.Vector3::new(0.0, 3.0, 3.0), "CA"),
    (@src.Vector3::new(3.0, 3.0, 3.0), "CA"),
  ]
  let res = cav_make_res("ALA", 1, atoms)
  let chain = @src.Chain::new(id='A', residues=[res])
  let result = @src.find_cavities(chain, grid_spacing=1.5)
  // A compact cluster should have no significant cavities
  // (may have small gaps but not enclosed)
  let has_large_cavity = result.cavities.length() > 0 && result.cavities[0].num_points > 5
  assert_false(has_large_cavity)
}

///|
test "cavifinder: cavity detection" {
  // Create a hollow shell: atoms on the surface of a sphere
  // leaving an empty interior
  let atoms = Array(capacity=0)
  let mut theta = 0
  while theta < 360 {
    let mut phi = 0
    while phi < 180 {
      let th_rad = theta.to_double() * 3.14159265358979323846 / 180.0
      let ph_rad = phi.to_double() * 3.14159265358979323846 / 180.0
      let r = 8.0
      let x = r * @math.sin(ph_rad) * @math.cos(th_rad)
      let y = r * @math.sin(ph_rad) * @math.sin(th_rad)
      let z = r * @math.cos(ph_rad)
      atoms.push((@src.Vector3::new(x, y, z), "CA"))
      phi = phi + 20
    }
    theta = theta + 20
  }
  let res = cav_make_res("ALA", 1, atoms)
  let chain = @src.Chain::new(id='A', residues=[res])
  let result = @src.find_cavities(chain, grid_spacing=2.0, probe_radius=2.5)
  // Should detect at least one cavity inside the shell
  assert_true(result.cavities.length() > 0)
  assert_true(result.total_volume > 0.0)
}

///|
test "cavifinder: grid spacing affects resolution" {
  // Same structure, different grid spacing
  let atoms = [
    (@src.Vector3::new(0.0, 0.0, 0.0), "CA"),
    (@src.Vector3::new(5.0, 0.0, 0.0), "CA"),
    (@src.Vector3::new(0.0, 5.0, 0.0), "CA"),
    (@src.Vector3::new(0.0, 0.0, 5.0), "CA"),
    (@src.Vector3::new(5.0, 5.0, 0.0), "CA"),
    (@src.Vector3::new(5.0, 0.0, 5.0), "CA"),
    (@src.Vector3::new(0.0, 5.0, 5.0), "CA"),
    (@src.Vector3::new(5.0, 5.0, 5.0), "CA"),
  ]
  let res = cav_make_res("ALA", 1, atoms)
  let chain = @src.Chain::new(id='A', residues=[res])
  let r1 = @src.find_cavities(chain, grid_spacing=1.0)
  let r2 = @src.find_cavities(chain, grid_spacing=2.0)
  // Both should run without error
  assert_true(r1.grid_spacing == 1.0)
  assert_true(r2.grid_spacing == 2.0)
}
