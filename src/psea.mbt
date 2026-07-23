///|
/// Bio.PDB.PSEA module - Protein Secondary Structure Element Assignment.
/// Inspired by Biopython's Bio.PDB.PSEA wrapper and the PSEA algorithm.
/// Predicts secondary structure (H=helix, E=strand, C=coil) from
/// alpha-carbon (CA) atomic coordinates.

// ===== Core data structures =====

///|
/// A single atomic coordinate.
pub struct PseaAtom {
  /// Residue name (e.g., "ALA")
  res_name : String
  /// Residue number
  res_num : Int
  /// Atom name (e.g., "CA")
  atom_name : String
  /// X coordinate
  x : Double
  /// Y coordinate
  y : Double
  /// Z coordinate
  z : Double
}

///|
/// PSEA secondary structure assignment result.
pub struct PseaResult {
  /// Residue numbers
  res_nums : Array[Int]
  /// Residue names
  res_names : Array[String]
  /// Secondary structure per residue: "H" (helix), "E" (strand), "C" (coil)
  ss : Array[String]
  /// Confidence scores (0.0 - 1.0)
  confidence : Array[Double]
  /// Summary counts
  n_helix : Int
  /// Summary counts
  n_strand : Int
  /// Summary counts
  n_coil : Int
}

///|
/// Virtual torsion angles computed along the CA chain.
pub struct PseaGeometry {
  /// CA-CA distances (length n-1)
  ca_distances : Array[Double]
  /// Virtual bond angles (length n-2)
  bond_angles : Array[Double]
  /// Virtual dihedral angles in degrees (length n-3)
  dihedral_angles : Array[Double]
}

// ===== Geometric helpers =====

///|
/// Create a new atom.
pub fn PseaAtom::new(res_name : String, res_num : Int, atom_name : String, x : Double, y : Double, z : Double) -> PseaAtom {
  PseaAtom::{ res_name, res_num, atom_name, x, y, z }
}

///|
/// Compute Euclidean distance between two atoms.
pub fn psea_distance(a : PseaAtom, b : PseaAtom) -> Double {
  let dx = a.x - b.x
  let dy = a.y - b.y
  let dz = a.z - b.z
  (dx * dx + dy * dy + dz * dz).sqrt()
}

///|
/// Compute the virtual bond angle (in degrees) at atom b, formed by a-b-c.
pub fn psea_bond_angle(a : PseaAtom, b : PseaAtom, c : PseaAtom) -> Double {
  // Vectors ba and bc
  let ba_x = a.x - b.x
  let ba_y = a.y - b.y
  let ba_z = a.z - b.z
  let bc_x = c.x - b.x
  let bc_y = c.y - b.y
  let bc_z = c.z - b.z
  let dot = ba_x * bc_x + ba_y * bc_y + ba_z * bc_z
  let ba_len = (ba_x * ba_x + ba_y * ba_y + ba_z * ba_z).sqrt()
  let bc_len = (bc_x * bc_x + bc_y * bc_y + bc_z * bc_z).sqrt()
  if ba_len < 1.0e-9 || bc_len < 1.0e-9 {
    return 0.0
  }
  let cos_theta = dot / (ba_len * bc_len)
  // Clamp to avoid numerical issues
  let cos_clamped = if cos_theta > 1.0 { 1.0 } else if cos_theta < -1.0 { -1.0 } else { cos_theta }
  // arccos in degrees
  psea_arccos_deg(cos_clamped)
}

///|
/// Compute the virtual dihedral angle (in degrees) formed by atoms a-b-c-d.
pub fn psea_dihedral_angle(a : PseaAtom, b : PseaAtom, c : PseaAtom, d : PseaAtom) -> Double {
  // Vectors b1 = b-a, b2 = c-b, b3 = d-c
  let b1_x = b.x - a.x
  let b1_y = b.y - a.y
  let b1_z = b.z - a.z
  let b2_x = c.x - b.x
  let b2_y = c.y - b.y
  let b2_z = c.z - b.z
  let b3_x = d.x - c.x
  let b3_y = d.y - c.y
  let b3_z = d.z - c.z
  // Normal of plane 1: n1 = b1 x b2
  let n1_x = b1_y * b2_z - b1_z * b2_y
  let n1_y = b1_z * b2_x - b1_x * b2_z
  let n1_z = b1_x * b2_y - b1_y * b2_x
  // Normal of plane 2: n2 = b2 x b3
  let n2_x = b2_y * b3_z - b2_z * b3_y
  let n2_y = b2_z * b3_x - b2_x * b3_z
  let n2_z = b2_x * b3_y - b2_y * b3_x
  // Angle between normals
  let n1_len = (n1_x * n1_x + n1_y * n1_y + n1_z * n1_z).sqrt()
  let n2_len = (n2_x * n2_x + n2_y * n2_y + n2_z * n2_z).sqrt()
  if n1_len < 1.0e-9 || n2_len < 1.0e-9 {
    return 0.0
  }
  let dot = n1_x * n2_x + n1_y * n2_y + n1_z * n2_z
  let cos_val = dot / (n1_len * n2_len)
  let cos_clamped = if cos_val > 1.0 { 1.0 } else if cos_val < -1.0 { -1.0 } else { cos_val }
  let mut angle = psea_arccos_deg(cos_clamped)
  // Sign of the dihedral: (n1 x n2) . b2
  let sign_vec_x = n1_y * n2_z - n1_z * n2_y
  let sign_vec_y = n1_z * n2_x - n1_x * n2_z
  let sign_vec_z = n1_x * n2_y - n1_y * n2_x
  let sign_dot = sign_vec_x * b2_x + sign_vec_y * b2_y + sign_vec_z * b2_z
  if sign_dot < 0.0 {
    angle = -angle
  }
  angle
}

///|
/// Arccosine in degrees.
fn psea_arccos_deg(x : Double) -> Double {
  @math.acos(x) * 180.0 / 3.141592653589793
}

///|
/// Arctangent2 in degrees.
fn psea_atan2_deg(y : Double, x : Double) -> Double {
  @math.atan2(y, x) * 180.0 / 3.141592653589793
}

// ===== Geometry computation =====

///|
/// Compute geometric features (CA-CA distances, bond angles, dihedral angles)
/// for a chain of CA atoms.
pub fn psea_compute_geometry(ca_atoms : Array[PseaAtom]) -> PseaGeometry {
  let n = ca_atoms.length()
  let ca_distances : Array[Double] = Array::new()
  let bond_angles : Array[Double] = Array::new()
  let dihedral_angles : Array[Double] = Array::new()
  // CA-CA distances
  let mut i = 1
  while i < n {
    ca_distances.push(psea_distance(ca_atoms[i - 1], ca_atoms[i]))
    i = i + 1
  }
  // Bond angles (a-b-c)
  let mut j = 2
  while j < n {
    bond_angles.push(psea_bond_angle(ca_atoms[j - 2], ca_atoms[j - 1], ca_atoms[j]))
    j = j + 1
  }
  // Dihedral angles (a-b-c-d)
  let mut k = 3
  while k < n {
    dihedral_angles.push(psea_dihedral_angle(ca_atoms[k - 3], ca_atoms[k - 2], ca_atoms[k - 1], ca_atoms[k]))
    k = k + 1
  }
  PseaGeometry::{ ca_distances, bond_angles, dihedral_angles }
}

// ===== PSEA secondary structure assignment =====

///|
/// Run the PSEA secondary structure assignment algorithm.
/// Uses a window-based approach on CA-CA geometry:
/// - Helices (H): regions where virtual dihedral angles fall in alpha-helix range
///   (~ -60 degrees) for at least 4 consecutive residues.
/// - Strands (E): regions where dihedral angles fall in beta-sheet range
///   (~ -120 or +120 degrees) for at least 3 consecutive residues.
/// - Coil (C): everything else.
pub fn psea_run(ca_atoms : Array[PseaAtom]) -> PseaResult {
  let n = ca_atoms.length()
  let ss : Array[String] = Array::new()
  let confidence : Array[Double] = Array::new()
  // Initialize all to coil
  let mut idx = 0
  while idx < n {
    ss.push("C")
    confidence.push(0.5)
    idx = idx + 1
  }
  if n < 4 {
    // Not enough atoms for dihedral analysis; mark as coil
    let res_nums : Array[Int] = Array::new()
    let res_names : Array[String] = Array::new()
    let mut r = 0
    while r < n {
      res_nums.push(ca_atoms[r].res_num)
      res_names.push(ca_atoms[r].res_name)
      r = r + 1
    }
    return PseaResult::{
      res_nums,
      res_names,
      ss,
      confidence,
      n_helix: 0,
      n_strand: 0,
      n_coil: n
    }
  }
  let geom = psea_compute_geometry(ca_atoms)
  // Score each residue based on local dihedral angles
  // Dihedral angle for residues i..i+3 corresponds to central residue i+1
  // Map dihedral angle index k to residue (k+1)
  let helix_scores : Array[Double] = Array::new()
  let strand_scores : Array[Double] = Array::new()
  let mut s = 0
  while s < n {
    helix_scores.push(0.0)
    strand_scores.push(0.0)
    s = s + 1
  }
  let mut di = 0
  while di < geom.dihedral_angles.length() {
    let angle = geom.dihedral_angles[di]
    let res_idx = di + 1
    // Alpha helix: dihedral around -57 to -70 (centered at -60)
    let helix_dist = (angle - (-60.0)).abs()
    let helix_score = if helix_dist < 30.0 { 1.0 - helix_dist / 30.0 } else { 0.0 }
    // Beta strand: dihedral around -120 or +120
    let strand_dist1 = (angle - (-120.0)).abs()
    let strand_dist2 = (angle - 120.0).abs()
    let strand_dist = if strand_dist1 < strand_dist2 { strand_dist1 } else { strand_dist2 }
    let strand_score = if strand_dist < 40.0 { 1.0 - strand_dist / 40.0 } else { 0.0 }
    // Apply to residues in the window
    let mut ri = res_idx - 1
    while ri <= res_idx + 1 && ri < n {
      if ri >= 0 {
        if helix_score > helix_scores[ri] {
          helix_scores[ri] = helix_score
        }
        if strand_score > strand_scores[ri] {
          strand_scores[ri] = strand_score
        }
      }
      ri = ri + 1
    }
    di = di + 1
  }
  // Smooth and assign: a residue is helix if helix_score > strand_score and helix_score > 0.3
  let mut i = 0
  while i < n {
    let h = helix_scores[i]
    let e = strand_scores[i]
    if h > e && h > 0.3 {
      ss[i] = "H"
      confidence[i] = h
    } else if e > h && e > 0.3 {
      ss[i] = "E"
      confidence[i] = e
    } else {
      ss[i] = "C"
      let max_score = if h > e { h } else { e }
      confidence[i] = 1.0 - max_score
    }
    i = i + 1
  }
  // Enforce minimum segment lengths:
  // - Helices need at least 4 consecutive residues
  // - Strands need at least 3 consecutive residues
  psea_enforce_min_lengths(ss, 4, 3)
  // Build result
  let res_nums : Array[Int] = Array::new()
  let res_names : Array[String] = Array::new()
  let mut n_helix = 0
  let mut n_strand = 0
  let mut n_coil = 0
  let mut r = 0
  while r < n {
    res_nums.push(ca_atoms[r].res_num)
    res_names.push(ca_atoms[r].res_name)
    if ss[r] == "H" {
      n_helix = n_helix + 1
    } else if ss[r] == "E" {
      n_strand = n_strand + 1
    } else {
      n_coil = n_coil + 1
    }
    r = r + 1
  }
  PseaResult::{ res_nums, res_names, ss, confidence, n_helix, n_strand, n_coil }
}

///|
/// Enforce minimum lengths for helix and strand segments.
fn psea_enforce_min_lengths(ss : Array[String], min_helix : Int, min_strand : Int) -> Unit {
  let n = ss.length()
  let mut i = 0
  while i < n {
    let current = ss[i]
    if current == "H" || current == "E" {
      let min_len = if current == "H" { min_helix } else { min_strand }
      // Find run end
      let mut j = i
      while j < n && ss[j] == current {
        j = j + 1
      }
      let run_len = j - i
      if run_len < min_len {
        // Convert to coil
        let mut k = i
        while k < j {
          ss[k] = "C"
          k = k + 1
        }
      }
      i = j
    } else {
      i = i + 1
    }
  }
  // After pruning, there might be adjacent strands that should merge; ignore for simplicity
}

// ===== Result accessors =====

///|
/// Get the secondary structure string (concatenated, e.g., "CCCHHHHHCCCCEEECCC").
pub fn PseaResult::psea_ss_string(self : PseaResult) -> String {
  let mut s = ""
  let mut i = 0
  while i < self.ss.length() {
    s = s + self.ss[i]
    i = i + 1
  }
  s
}

///|
/// Get number of residues.
pub fn PseaResult::psea_n_residues(self : PseaResult) -> Int {
  self.ss.length()
}

///|
/// Get the secondary structure assignment at an index.
pub fn PseaResult::psea_get_ss(self : PseaResult, idx : Int) -> String {
  self.ss[idx]
}

///|
/// Get the confidence score at an index.
pub fn PseaResult::psea_get_confidence(self : PseaResult, idx : Int) -> Double {
  self.confidence[idx]
}

///|
/// Get the residue number at an index.
pub fn PseaResult::psea_get_res_num(self : PseaResult, idx : Int) -> Int {
  self.res_nums[idx]
}

///|
/// Get the residue name at an index.
pub fn PseaResult::psea_get_res_name(self : PseaResult, idx : Int) -> String {
  self.res_names[idx]
}

///|
/// Summary statistics as a string.
pub fn PseaResult::psea_summary(self : PseaResult) -> String {
  "PSEA: " + self.psea_n_residues().to_string() + " residues, " +
  "H=" + self.n_helix.to_string() + " " +
  "E=" + self.n_strand.to_string() + " " +
  "C=" + self.n_coil.to_string()
}

// ===== Three-state to eight-state conversion (DSSP-like) =====

///|
/// Convert a 3-state (H/E/C) assignment to an 8-state DSSP-like assignment.
/// H -> H (alpha helix), G (3-10 helix if short), I (pi helix if very long)
/// E -> E (extended strand), B (bridge if isolated)
/// C -> T (turn), S (bend), L (loop), C (coil)
pub fn psea_three_to_eight_state(ss : Array[String]) -> Array[String] {
  let n = ss.length()
  let result : Array[String] = Array::new()
  let mut i = 0
  while i < n {
    let s = ss[i]
    if s == "H" {
      // Count helix length
      let mut j = i
      while j < n && ss[j] == "H" {
        j = j + 1
      }
      let len = j - i
      let mut k = i
      while k < j {
        let state = if len <= 3 {
          "G"
        } else if len >= 12 {
          "I"
        } else {
          "H"
        }
        result.push(state)
        k = k + 1
      }
      i = j
    } else if s == "E" {
      let mut j = i
      while j < n && ss[j] == "E" {
        j = j + 1
      }
      let len = j - i
      let mut k = i
      while k < j {
        let state = if len < 2 { "B" } else { "E" }
        result.push(state)
        k = k + 1
      }
      i = j
    } else {
      // Coil: try to distinguish T/S/L
      let prev = if i > 0 { ss[i - 1] } else { "C" }
      let next = if i + 1 < n { ss[i + 1] } else { "C" }
      let state = if (prev == "H" || next == "H") && (prev != next) {
        "T"
      } else if prev == "E" || next == "E" {
        "T"
      } else {
        "L"
      }
      result.push(state)
      i = i + 1
    }
  }
  result
}

// ===== Example data generators =====

///|
/// Generate a small example CA chain with an alpha-helix-like geometry.
/// Produces ~12 residues in a helical arrangement.
pub fn psea_create_example_helix() -> Array[PseaAtom] {
  // Approximate alpha-helix CA coordinates: helix rise ~1.5 Å, radius ~2.3 Å
  // 100 degrees per residue (3.6 residues per turn)
  let atoms : Array[PseaAtom] = Array::new()
  let names = ["ALA", "GLY", "VAL", "LEU", "SER", "THR", "PHE", "TRP", "TYR", "CYS", "MET", "ILE"]
  let mut i = 0
  while i < 12 {
    let angle = i.to_double() * 100.0 * 3.141592653589793 / 180.0
    let radius = 2.3
    let x = radius * @math.cos(angle)
    let y = radius * @math.sin(angle)
    let z = i.to_double() * 1.5
    let rname = if i < names.length() { names[i] } else { "ALA" }
    atoms.push(PseaAtom::new(rname, i + 1, "CA", x, y, z))
    i = i + 1
  }
  atoms
}

///|
/// Generate a small example CA chain with a beta-strand-like geometry.
/// Produces ~8 residues in an extended conformation.
pub fn psea_create_example_strand() -> Array[PseaAtom] {
  let atoms : Array[PseaAtom] = Array::new()
  let names = ["ALA", "GLY", "VAL", "LEU", "SER", "THR", "PHE", "TRP"]
  // Extended strand: alternating zig-zag with ~3.8 Å CA-CA distance
  let mut i = 0
  while i < 8 {
    let x = i.to_double() * 3.4
    let y = if i % 2 == 0 { 1.0 } else { -1.0 }
    let z = 0.0
    let rname = if i < names.length() { names[i] } else { "ALA" }
    atoms.push(PseaAtom::new(rname, i + 1, "CA", x, y, z))
    i = i + 1
  }
  atoms
}

///|
/// Generate a mixed example with helix + coil + strand regions.
pub fn psea_create_example_mixed() -> Array[PseaAtom] {
  let helix_part = psea_create_example_helix()
  let strand_part = psea_create_example_strand()
  let atoms : Array[PseaAtom] = Array::new()
  // Take first 6 helix residues, then a 3-residue coil, then 6 strand residues
  let mut idx = 0
  while idx < 6 && idx < helix_part.length() {
    let a = helix_part[idx]
    atoms.push(PseaAtom::new(a.res_name, idx + 1, "CA", a.x, a.y, a.z))
    idx = idx + 1
  }
  // Coil region (random-ish positions)
  atoms.push(PseaAtom::new("PRO", 7, "CA", 12.0, 5.0, 8.0))
  atoms.push(PseaAtom::new("GLY", 8, "CA", 13.0, 6.0, 9.0))
  atoms.push(PseaAtom::new("SER", 9, "CA", 14.5, 6.5, 9.5))
  // Strand region
  let mut s = 0
  while s < 6 && s < strand_part.length() {
    let a = strand_part[s]
    atoms.push(PseaAtom::new(a.res_name, 10 + s, "CA", a.x + 18.0, a.y, a.z))
    s = s + 1
  }
  atoms
}
