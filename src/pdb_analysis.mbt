///|
/// Bio.PDB.StructureAnalysis - Advanced protein structure analysis.
/// Inspired by Bio.PDB's polypeptide, DSSP, and hydrogen bond modules.

///|
/// Dihedral angles (phi, psi, omega) for a residue.
pub struct DihedralAngles {
  phi : Double?
  psi : Double?
  omega : Double?
}

///|
/// Construct DihedralAngles.
pub fn DihedralAngles::new(
  phi? : Double? = None,
  psi? : Double? = None,
  omega? : Double? = None,
) -> DihedralAngles {
  DihedralAngles::{ phi, psi, omega }
}

///|
/// Secondary structure types (DSSP-style).
pub enum SSType {
  AlphaHelix // H
  Helix3_10 // G
  PiHelix // I
  ExtendedStrand // E
  IsolatedBridge // B
  Turn // T
  Bend // S
  Coil // -
} derive(Eq, Show)

///|
/// Construct AlphaHelix SSType.
pub fn ss_alpha_helix() -> SSType {
  AlphaHelix
}

///|
/// Construct Helix3_10 SSType.
pub fn ss_helix_3_10() -> SSType {
  Helix3_10
}

///|
/// Construct PiHelix SSType.
pub fn ss_pi_helix() -> SSType {
  PiHelix
}

///|
/// Construct ExtendedStrand SSType.
pub fn ss_extended_strand() -> SSType {
  ExtendedStrand
}

///|
/// Construct IsolatedBridge SSType.
pub fn ss_isolated_bridge() -> SSType {
  IsolatedBridge
}

///|
/// Construct Turn SSType.
pub fn ss_turn() -> SSType {
  Turn
}

///|
/// Construct Bend SSType.
pub fn ss_bend() -> SSType {
  Bend
}

///|
/// Construct Coil SSType.
pub fn ss_coil() -> SSType {
  Coil
}

///|
/// Convert SSType to DSSP single-letter code.
pub fn SSType::to_code(self : SSType) -> Char {
  match self {
    AlphaHelix => 'H'
    Helix3_10 => 'G'
    PiHelix => 'I'
    ExtendedStrand => 'E'
    IsolatedBridge => 'B'
    Turn => 'T'
    Bend => 'S'
    Coil => '-'
  }
}

///|
/// Convert DSSP single-letter code to SSType.
pub fn SSType::from_code(code : Char) -> SSType {
  match code {
    'H' => AlphaHelix
    'G' => Helix3_10
    'I' => PiHelix
    'E' => ExtendedStrand
    'B' => IsolatedBridge
    'T' => Turn
    'S' => Bend
    _ => Coil
  }
}

///|
/// Calculate the angle between three points (BAC convention: angle at A).
fn calc_angle(a : Vector3, b : Vector3, c : Vector3) -> Double {
  let ab = Vector3::new(b.x - a.x, b.y - a.y, b.z - a.z)
  let cb = Vector3::new(b.x - c.x, b.y - c.y, b.z - c.z)
  let dot = ab.x * cb.x + ab.y * cb.y + ab.z * cb.z
  let mag_ab = (ab.x * ab.x + ab.y * ab.y + ab.z * ab.z).sqrt()
  let mag_cb = (cb.x * cb.x + cb.y * cb.y + cb.z * cb.z).sqrt()
  let cos_val = dot / (mag_ab * mag_cb)
  // Clamp to [-1, 1] to avoid NaN from floating point errors
  let clamped = if cos_val > 1.0 {
    1.0
  } else if cos_val < -1.0 {
    -1.0
  } else {
    cos_val
  }
  @math.acos(clamped)
}

///|
/// Cross product of two vectors.
fn cross_product(a : Vector3, b : Vector3) -> Vector3 {
  Vector3::new(
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x,
  )
}

///|
/// Dot product of two vectors.
fn dot_product(a : Vector3, b : Vector3) -> Double {
  a.x * b.x + a.y * b.y + a.z * b.z
}

///|
/// Calculate the dihedral angle formed by four atoms (a-b-c-d).
/// Returns the angle in degrees, in range [-180, 180].
pub fn calc_dihedral(
  a : Vector3,
  b : Vector3,
  c : Vector3,
  d : Vector3,
) -> Double {
  let b1 = Vector3::new(b.x - a.x, b.y - a.y, b.z - a.z)
  let b2 = Vector3::new(c.x - b.x, c.y - b.y, c.z - b.z)
  let b3 = Vector3::new(d.x - c.x, d.y - c.y, d.z - c.z)

  let n1 = cross_product(b1, b2)
  let n2 = cross_product(b2, b3)
  let m1 = cross_product(n1, Vector3::new(b2.x, b2.y, b2.z).normalize_safe())

  let x = dot_product(n1, n2)
  let y = dot_product(m1, n2)

  let angle = @math.atan2(y, x)
  angle * 180.0 / 3.14159265358979323846
}

///|
/// Normalize a Vector3 to unit length (with safe handling of zero vector).
fn Vector3::normalize_safe(self : Vector3) -> Vector3 {
  let mag = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
  if mag == 0.0 {
    return Vector3::new(0.0, 0.0, 0.0)
  }
  Vector3::new(self.x / mag, self.y / mag, self.z / mag)
}

///|
/// A hydrogen bond detected in the structure.
pub struct HydrogenBond {
  donor_resname : String
  donor_resseq : Int
  donor_chainid : Char
  acceptor_resname : String
  acceptor_resseq : Int
  acceptor_chainid : Char
  distance : Double
  angle : Double
}

///|
/// Construct a HydrogenBond.
pub fn HydrogenBond::new(
  donor_resname~ : String,
  donor_resseq~ : Int,
  donor_chainid~ : Char,
  acceptor_resname~ : String,
  acceptor_resseq~ : Int,
  acceptor_chainid~ : Char,
  distance~ : Double,
  angle~ : Double,
) -> HydrogenBond {
  HydrogenBond::{
    donor_resname,
    donor_resseq,
    donor_chainid,
    acceptor_resname,
    acceptor_resseq,
    acceptor_chainid,
    distance,
    angle,
  }
}

///|
/// Calculate phi angle for a residue.
/// Requires C(i-1), N(i), CA(i), C(i) atoms.
pub fn calc_phi(prev_res : Residue, res : Residue) -> Double? {
  let c_prev = prev_res.find_atom("C")
  let n = res.find_atom("N")
  let ca = res.find_atom("CA")
  let c = res.find_atom("C")
  match (c_prev, n, ca, c) {
    (Some(cp), Some(n), Some(ca), Some(c)) =>
      Some(calc_dihedral(cp.coord, n.coord, ca.coord, c.coord))
    _ => None
  }
}

///|
/// Calculate psi angle for a residue.
/// Requires N(i), CA(i), C(i), N(i+1) atoms.
pub fn calc_psi(res : Residue, next_res : Residue) -> Double? {
  let n = res.find_atom("N")
  let ca = res.find_atom("CA")
  let c = res.find_atom("C")
  let n_next = next_res.find_atom("N")
  match (n, ca, c, n_next) {
    (Some(n), Some(ca), Some(c), Some(n_next)) =>
      Some(calc_dihedral(n.coord, ca.coord, c.coord, n_next.coord))
    _ => None
  }
}

///|
/// Calculate omega angle between two residues.
/// Requires CA(i), C(i), N(i+1), CA(i+1) atoms.
pub fn calc_omega(res : Residue, next_res : Residue) -> Double? {
  let ca = res.find_atom("CA")
  let c = res.find_atom("C")
  let n_next = next_res.find_atom("N")
  let ca_next = next_res.find_atom("CA")
  match (ca, c, n_next, ca_next) {
    (Some(ca), Some(c), Some(n_next), Some(ca_next)) =>
      Some(calc_dihedral(ca.coord, c.coord, n_next.coord, ca_next.coord))
    _ => None
  }
}

///|
/// Calculate all dihedral angles for a chain of residues.
/// Returns an array of (resseq, phi, psi, omega) tuples.
pub fn calc_chain_dihedrals(
  chain : Chain,
) -> Array[(Int, Double?, Double?, Double?)] {
  let result : Array[(Int, Double?, Double?, Double?)] = Array::new()
  let residues = chain.get_residues()
  let n = residues.length()
  let mut i = 0
  while i < n {
    let res = residues[i]
    let phi = if i > 0 { calc_phi(residues[i - 1], res) } else { None }
    let psi = if i < n - 1 { calc_psi(res, residues[i + 1]) } else { None }
    let omega = if i < n - 1 { calc_omega(res, residues[i + 1]) } else { None }
    result.push((res.resseq, phi, psi, omega))
    i = i + 1
  }
  result
}

///|
/// Compute the distance matrix between all CA atoms in a chain.
/// Returns a 2D array where element [i][j] is the distance
/// between CA of residue i and CA of residue j.
pub fn calc_ca_distance_matrix(chain : Chain) -> Array[Array[Double]] {
  let residues = chain.get_residues()
  let ca_atoms : Array[Vector3] = Array::new()
  for res in residues {
    match res.get_ca() {
      Some(ca) => ca_atoms.push(ca.coord)
      None => ()
    }
  }
  let n = ca_atoms.length()
  let matrix : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      row.push(ca_atoms[i].distance(ca_atoms[j]))
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  matrix
}

///|
/// Compute the contact map for CA atoms.
/// Two residues are in contact if their CA distance < threshold.
pub fn calc_contact_map(
  chain : Chain,
  threshold? : Double = 8.0,
) -> Array[Array[Bool]] {
  let matrix = calc_ca_distance_matrix(chain)
  let n = matrix.length()
  let contact : Array[Array[Bool]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Bool] = Array::new()
    let mut j = 0
    while j < n {
      row.push(matrix[i][j] < threshold && i != j)
      j = j + 1
    }
    contact.push(row)
    i = i + 1
  }
  contact
}

///|
/// Detect hydrogen bonds in a chain based on N-H...O=C geometry.
/// Uses a simple distance + angle criterion.
/// donor_max_dist defaults to 3.5 Angstroms.
/// min_angle defaults to 90 degrees.
pub fn detect_hydrogen_bonds(
  chain : Chain,
  max_dist? : Double = 3.5,
  min_angle? : Double = 90.0,
) -> Array[HydrogenBond] {
  let result : Array[HydrogenBond] = Array::new()
  let residues = chain.get_residues()
  let n = residues.length()
  let mut i = 0
  while i < n {
    let donor = residues[i]
    let n_atom = donor.find_atom("N")
    let ca_atom = donor.find_atom("CA")
    let c_atom = donor.find_atom("C")
    let mut j = 0
    while j < n {
      if i != j {
        let acceptor = residues[j]
        let o_atom = acceptor.find_atom("O")
        match (n_atom, ca_atom, c_atom, o_atom) {
          (Some(n), Some(ca), Some(c), Some(o)) => {
            let dist = n.coord.distance(o.coord)
            if dist <= max_dist {
              // Simplified angle check: N-H...O angle
              // Using N-CA-C-O as proxy
              let angle = calc_angle(ca.coord, n.coord, o.coord)
              let angle_deg = angle * 180.0 / 3.14159265358979323846
              if angle_deg >= min_angle {
                result.push(
                  HydrogenBond::new(
                    donor_resname=donor.resname,
                    donor_resseq=donor.resseq,
                    donor_chainid=donor.chainid,
                    acceptor_resname=acceptor.resname,
                    acceptor_resseq=acceptor.resseq,
                    acceptor_chainid=acceptor.chainid,
                    distance=dist,
                    angle=angle_deg,
                  ),
                )
              }
            }
          }
          _ => ()
        }
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Assign secondary structure based on phi/psi angles and hydrogen bond patterns.
/// Uses a simplified Ramachandran-based classification.
pub fn assign_secondary_structure(chain : Chain) -> Array[(Int, SSType)] {
  let result : Array[(Int, SSType)] = Array::new()
  let dihedrals = calc_chain_dihedrals(chain)
  let hbonds = detect_hydrogen_bonds(chain)

  // Build hydrogen bond map: (donor_resseq -> acceptor_resseq)
  let hbond_map : Map[Int, Array[Int]] = Map([], capacity=16)
  for hb in hbonds {
    match hbond_map.get(hb.donor_resseq) {
      Some(arr) => arr.push(hb.acceptor_resseq)
      None => {
        let arr : Array[Int] = Array::new()
        arr.push(hb.acceptor_resseq)
        hbond_map[hb.donor_resseq] = arr
      }
    }
  }

  let n = dihedrals.length()
  let mut i = 0
  while i < n {
    let (resseq, phi_opt, psi_opt, _) = dihedrals[i]
    let ss = match (phi_opt, psi_opt) {
      (Some(phi), Some(psi)) =>
        // Alpha helix region: phi ~ -60, psi ~ -45
        if phi > -90.0 && phi < -30.0 && psi > -70.0 && psi < -10.0 {
          // Check for i -> i+4 hydrogen bond (helix signature)
          let has_hbond = match hbond_map.get(resseq) {
            Some(arr) => {
              let mut found = false
              for target in arr {
                if target == resseq + 4 || target == resseq - 4 {
                  found = true
                  break
                }
              }
              found
            }
            None => false
          }
          if has_hbond {
            AlphaHelix
          } else {
            Turn
          }
          // Beta strand region: phi ~ -120, psi ~ 130
        } else if phi < -60.0 && phi > -180.0 && psi > 90.0 && psi < 180.0 {
          ExtendedStrand
          // Turn region
        } else if phi > -60.0 && phi < 60.0 && psi > 0.0 && psi < 180.0 {
          Turn
        } else {
          Coil
        }
      _ => Coil
    }
    result.push((resseq, ss))
    i = i + 1
  }
  result
}

///|
/// Count residues by secondary structure type.
pub fn count_secondary_structure(
  ss_assignments : Array[(Int, SSType)],
) -> Map[String, Int] {
  let result : Map[String, Int] = Map([], capacity=8)
  for entry in ss_assignments {
    let ss = entry.1
    let key = ss.to_code().to_string()
    match result.get(key) {
      Some(count) => result[key] = count + 1
      None => result[key] = 1
    }
  }
  result
}

///|
/// Calculate the radius of gyration of a structure (using CA atoms).
/// Measures the overall size/compactness of the protein.
pub fn calc_radius_of_gyration(chain : Chain) -> Double? {
  let residues = chain.get_residues()
  let coords : Array[Vector3] = Array::new()
  for res in residues {
    match res.get_ca() {
      Some(ca) => coords.push(ca.coord)
      None => ()
    }
  }
  let n = coords.length()
  if n == 0 {
    return None
  }
  // Calculate centroid
  let mut cx = 0.0
  let mut cy = 0.0
  let mut cz = 0.0
  for c in coords {
    cx = cx + c.x
    cy = cy + c.y
    cz = cz + c.z
  }
  cx = cx / n.to_double()
  cy = cy / n.to_double()
  cz = cz / n.to_double()
  // Calculate sum of squared distances
  let mut sum_sq = 0.0
  for c in coords {
    let dx = c.x - cx
    let dy = c.y - cy
    let dz = c.z - cz
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
  }
  Some((sum_sq / n.to_double()).sqrt())
}

///|
/// Calculate the Ramachandran plot data for a chain.
/// Returns an array of (resseq, resname, phi, psi) tuples.
pub fn ramachandran_plot(chain : Chain) -> Array[(Int, String, Double, Double)] {
  let result : Array[(Int, String, Double, Double)] = Array::new()
  let dihedrals = calc_chain_dihedrals(chain)
  let residues = chain.get_residues()
  let n = dihedrals.length()
  let mut i = 0
  while i < n {
    let (resseq, phi_opt, psi_opt, _) = dihedrals[i]
    let resname = residues[i].resname
    match (phi_opt, psi_opt) {
      (Some(phi), Some(psi)) => result.push((resseq, resname, phi, psi))
      _ => ()
    }
    i = i + 1
  }
  result
}

///|
/// Classify a residue's position in the Ramachandran plot.
pub enum RamachandranRegion {
  Favored // Most favorable region
  Allowed // Allowed region
  Generous // Generous/outlier region
} derive(Eq, Show)

///|
/// Construct Favored RamachandranRegion.
pub fn ramachandran_favored() -> RamachandranRegion {
  Favored
}

///|
/// Construct Allowed RamachandranRegion.
pub fn ramachandran_allowed() -> RamachandranRegion {
  Allowed
}

///|
/// Construct Generous RamachandranRegion.
pub fn ramachandran_generous() -> RamachandranRegion {
  Generous
}

///|
/// Classify phi/psi angle pair into Ramachandran region.
pub fn classify_ramachandran(phi : Double, psi : Double) -> RamachandranRegion {
  // Simplified classification based on standard Ramachandran plot regions
  // Alpha-helix region (favored)
  if phi > -90.0 && phi < -30.0 && psi > -70.0 && psi < -10.0 {
    Favored
    // Beta-sheet region (favored)
  } else if phi > -180.0 && phi < -60.0 && psi > 90.0 && psi < 180.0 {
    Favored
    // Left-handed helix (allowed)
  } else if phi > 30.0 && phi < 90.0 && psi > 0.0 && psi < 60.0 {
    Allowed
    // Generous regions (broader bounds)
  } else if phi > -150.0 && phi < 0.0 && psi > -120.0 && psi < 180.0 {
    Generous
  } else {
    Generous
  }
}

///|
/// Create a demo protein structure for testing.
pub fn create_demo_structure() -> Structure {
  // Create a small beta-hairpin-like structure
  let chain_a = Chain::new(id='A')
  // Residue 1: ALA
  let res1 = Residue::new(resname="ALA", chainid='A', resseq=1, atoms=[
    Atom::new(name="N", coord=Vector3::new(0.0, 0.0, 0.0)),
    Atom::new(name="CA", coord=Vector3::new(1.5, 0.0, 0.0)),
    Atom::new(name="C", coord=Vector3::new(2.0, 1.5, 0.0)),
    Atom::new(name="O", coord=Vector3::new(1.2, 2.3, 0.0)),
  ])
  // Residue 2: GLY
  let res2 = Residue::new(resname="GLY", chainid='A', resseq=2, atoms=[
    Atom::new(name="N", coord=Vector3::new(3.3, 1.8, 0.0)),
    Atom::new(name="CA", coord=Vector3::new(4.0, 3.0, 0.0)),
    Atom::new(name="C", coord=Vector3::new(5.5, 2.8, 0.0)),
    Atom::new(name="O", coord=Vector3::new(6.2, 3.7, 0.0)),
  ])
  // Residue 3: VAL
  let res3 = Residue::new(resname="VAL", chainid='A', resseq=3, atoms=[
    Atom::new(name="N", coord=Vector3::new(6.0, 1.5, 0.0)),
    Atom::new(name="CA", coord=Vector3::new(7.5, 1.3, 0.0)),
    Atom::new(name="C", coord=Vector3::new(8.0, 2.8, 0.0)),
    Atom::new(name="O", coord=Vector3::new(7.2, 3.6, 0.0)),
  ])
  // Residue 4: LYS
  let res4 = Residue::new(resname="LYS", chainid='A', resseq=4, atoms=[
    Atom::new(name="N", coord=Vector3::new(9.3, 3.1, 0.0)),
    Atom::new(name="CA", coord=Vector3::new(10.0, 4.3, 0.0)),
    Atom::new(name="C", coord=Vector3::new(11.5, 4.1, 0.0)),
    Atom::new(name="O", coord=Vector3::new(12.2, 5.0, 0.0)),
  ])
  let chain = Chain::new(id='A', residues=[res1, res2, res3, res4])
  let model = Model::new(id=0, chains=[chain])
  Structure::new(id="demo", models=[model])
}
