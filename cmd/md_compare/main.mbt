// Companion program that evaluates the same hard-coded MD-analysis test
// system as test/python/python_mdanalysis_reference.py and emits results as
// JSON. The two outputs are diffed by test/python/mdanalysis_compare.sh to
// verify the MoonBit implementation matches MDAnalysis.
//
// Run:  moon run cmd/md_compare

///|
struct Case {
  name : String
  result : String
}

///|
/// Format a double with 6 decimal places (rounded half up).
fn fmt6(x : Double) -> String {
  let neg = x < 0.0
  let a = if neg { 0.0 - x } else { x }
  let scaled = (a * 1000000.0 + 0.5).floor()
  let ip = (scaled / 1000000.0).floor().to_int()
  let fp = (scaled - ip.to_double() * 1000000.0).to_int()
  let sb = StringBuilder()
  if neg {
    sb.write_string("-")
  }
  sb.write_string(ip.to_string())
  sb.write_string(".")
  let fs = fp.to_string()
  let mut k = fs.length()
  while k < 6 {
    sb.write_string("0")
    k = k + 1
  }
  sb.write_string(fs)
  sb.to_string()
}

///|
fn record(cases : Array[Case], name : String, values : Array[Double]) -> Unit {
  let sb = StringBuilder()
  for i = 0; i < values.length(); i = i + 1 {
    if i > 0 {
      sb.write_string(",")
    }
    sb.write_string(fmt6(values[i]))
  }
  cases.push(Case::{ name, result: sb.to_string(), })
}

///|
fn main {
  let cases : Array[Case] = Array(capacity=0)
  // Hard-coded system: 3 atoms (N, C, O), 2 frames — must match
  // test/python/python_mdanalysis_reference.py exactly.
  let x0 : Array[Array[Double]] = [
    [0.5, 0.2, 0.1],
    [1.5, 0.7, -0.3],
    [2.2, -0.4, 0.6],
  ]
  let x1 : Array[Array[Double]] = [
    [0.1, -0.3, 0.9],
    [1.8, 0.5, 0.4],
    [2.5, 0.9, -0.7],
  ]
  let masses : Array[Double] = [14.007, 12.011, 15.999]

  // 1. pairwise distance matrix (frame 0)
  let dm = @lib.mdanalysis_distance_array(x0, x0)
  let flat : Array[Double] = Array(capacity=0)
  for row in dm {
    for v in row {
      flat.push(v)
    }
  }
  record(cases, "distance_matrix", flat)

  // 2. RMSD of frame 1 onto frame 0 (centered + optimal superposition)
  record(cases, "rmsd_sup", [@lib.mdanalysis_rmsd(x0, x1)])

  // 3. per-atom RMSF over the two frames
  let traj = @lib.Trajectory::new([x0, x1])
  record(cases, "rmsf", @lib.mdanalysis_rmsf(traj.frames))

  // 4-6. radius of gyration / centers via the AtomGroup object model
  let atoms : Array[@lib.UnivAtom] = Array(capacity=0)
  let names : Array[String] = ["N", "C", "O"]
  for i = 0; i < 3; i = i + 1 {
    atoms.push(@lib.UnivAtom::{
      name: names[i],
      resname: "ALA",
      resid: 1,
      chain: "A",
      segid: "A",
      element: names[i],
      atom_type: "",
      charge: 0.0,
      mass: masses[i],
      is_protein: true,
      is_water: false,
      is_backbone: false,
    })
  }
  let u = @lib.Universe::new(atoms, traj)
  let g = u.all()
  record(cases, "radius_of_gyration", [g.radius_of_gyration()])
  record(cases, "center_of_geometry", g.center_of_geometry())
  record(cases, "center_of_mass", g.center_of_mass())

  // 7. dihedral angle of four points (radians, absolute value)
  let p1 = @lib.Vector3::from_array([0.0, 0.0, 0.0])
  let p2 = @lib.Vector3::from_array([1.0, 0.0, 0.0])
  let p3 = @lib.Vector3::from_array([1.0, 1.0, 0.0])
  let p4 = @lib.Vector3::from_array([1.0, 1.0, 1.0])
  let d = @lib.vector_dihedral_angle(p1, p2, p3, p4)
  record(cases, "dihedral_abs", [if d < 0.0 { 0.0 - d } else { d }])

  // 8. single bond distance
  record(cases, "bond_distance", [dm[0][1]])

  // emit JSON
  let sb = StringBuilder()
  sb.write_string("[\n")
  for i = 0; i < cases.length(); i = i + 1 {
    let c = cases[i]
    sb.write_string("  {\"name\": \"")
    sb.write_string(c.name)
    sb.write_string("\", \"result\": \"")
    sb.write_string(c.result)
    sb.write_string("\"}")
    if i + 1 < cases.length() {
      sb.write_string(",")
    }
    sb.write_string("\n")
  }
  sb.write_string("]\n")
  println(sb.to_string())
}
