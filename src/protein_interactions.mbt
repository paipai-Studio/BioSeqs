///|
/// Protein interaction analysis (InSty): salt bridges, π-π stacking,
/// π-cation, hydrophobic interactions, repulsive ionic bonding,
/// and water bridges.
///
/// Ported from ProDy's InSty suite.

// ─── helpers ───────────────────────────────────────────────

fn insty_sqrt(x : Double) -> Double {
  if x < 0.0 { 0.0 } else { x.sqrt() }
}

fn insty_cross(a : Vector3, b : Vector3) -> Vector3 {
  Vector3::new(
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x,
  )
}

fn insty_dot(a : Vector3, b : Vector3) -> Double {
  a.x * b.x + a.y * b.y + a.z * b.z
}

fn insty_normalize(v : Vector3) -> Vector3 {
  let mag = insty_sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
  if mag == 0.0 { Vector3::new(0.0, 0.0, 0.0) }
  else { Vector3::new(v.x / mag, v.y / mag, v.z / mag) }
}

///|
/// Positive-charge atom names for charged residues.
fn get_positive_atoms(resname : String) -> Array[String] {
  match resname {
    "ARG" => ["NH1", "NH2", "NE"]
    "LYS" => ["NZ"]
    "HIS" => ["ND1", "NE2"]
    _ => []
  }
}

///|
/// Negative-charge atom names for charged residues.
fn get_negative_atoms(resname : String) -> Array[String] {
  match resname {
    "ASP" => ["OD1", "OD2"]
    "GLU" => ["OE1", "OE2"]
    _ => []
  }
}

///|
/// Aromatic ring atom names for aromatic residues.
fn get_aromatic_ring_atoms(resname : String) -> Array[String]? {
  match resname {
    "PHE" => Some(["CG", "CD1", "CE1", "CZ", "CE2", "CD2"])
    "TYR" => Some(["CG", "CD1", "CE1", "CZ", "CE2", "CD2"])
    "TRP" => Some(["CD2", "CE2", "CZ2", "CH2", "CZ3", "CE3"])
    "HIS" => Some(["CG", "ND1", "CE1", "NE2", "CD2"])
    _ => None
  }
}

///|
/// Hydrophobic sidechain carbon atom names.
fn get_hydrophobic_atoms(resname : String) -> Array[String] {
  match resname {
    "ALA" => ["CB"]
    "VAL" => ["CB", "CG1", "CG2"]
    "LEU" => ["CB", "CG", "CD1", "CD2"]
    "ILE" => ["CB", "CG1", "CG2", "CD1"]
    "PHE" => ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
    "MET" => ["CB", "CG", "SD", "CE"]
    "TRP" => ["CB", "CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"]
    "PRO" => ["CB", "CG", "CD"]
    _ => []
  }
}

///|
/// Check if a residue is hydrophobic.
fn is_hydrophobic(resname : String) -> Bool {
  match resname {
    "ALA" | "VAL" | "LEU" | "ILE" | "PHE" | "MET" | "TRP" | "PRO" => true
    _ => false
  }
}

///|
/// Calculate the centroid of a set of atoms in a residue.
fn calc_ring_centroid(res : Residue, atom_names : Array[String]) -> Vector3? {
  let mut sx = 0.0
  let mut sy = 0.0
  let mut sz = 0.0
  let mut count = 0
  for name in atom_names {
    match res.find_atom(name) {
      Some(a) => {
        sx = sx + a.coord.x
        sy = sy + a.coord.y
        sz = sz + a.coord.z
        count = count + 1
      }
      None => ()
    }
  }
  if count == 0 { return None }
  Some(Vector3::new(sx / count.to_double(), sy / count.to_double(), sz / count.to_double()))
}

///|
/// Calculate the normal vector of an aromatic ring plane.
fn calc_ring_normal(res : Residue, atom_names : Array[String]) -> Vector3? {
  if atom_names.length() < 3 { return None }
  let a1 = res.find_atom(atom_names[0])
  let a2 = res.find_atom(atom_names[1])
  let a3 = res.find_atom(atom_names[2])
  match (a1, a2, a3) {
    (Some(p1), Some(p2), Some(p3)) => {
      let v1 = Vector3::new(p2.coord.x - p1.coord.x, p2.coord.y - p1.coord.y, p2.coord.z - p1.coord.z)
      let v2 = Vector3::new(p3.coord.x - p1.coord.x, p3.coord.y - p1.coord.y, p3.coord.z - p1.coord.z)
      let n = insty_cross(v1, v2)
      let normal = insty_normalize(n)
      let mag = insty_sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z)
      if mag < 0.001 { return None }
      Some(normal)
    }
    _ => None
  }
}

///|
/// Angle (degrees) between two normal vectors.
fn angle_between_normals(n1 : Vector3, n2 : Vector3) -> Double {
  let dot = insty_dot(n1, n2)
  let clamped = if dot > 1.0 { 1.0 } else if dot < -1.0 { -1.0 } else { dot }
  let angle = @math.acos(clamped)
  angle * 180.0 / 3.14159265358979323846
}

// ─── Salt Bridges ──────────────────────────────────────────

///|
/// A salt bridge between oppositely charged residues.
pub struct SaltBridge {
  resname1 : String
  resseq1 : Int
  chainid1 : Char
  resname2 : String
  resseq2 : Int
  chainid2 : Char
  distance : Double
  atom_name1 : String
  atom_name2 : String
}

///|
/// Find salt bridges between oppositely charged residues.
pub fn find_salt_bridges(
  chain : Chain,
  distance_cutoff? : Double = 4.0,
) -> Array[SaltBridge] {
  let result : Array[SaltBridge] = Array(capacity=0)
  let residues = chain.get_residues()
  let n = residues.length()
  let mut i = 0
  while i < n {
    let res1 = residues[i]
    let pos_atoms = get_positive_atoms(res1.resname)
    let neg_atoms = get_negative_atoms(res1.resname)
    let mut j = i + 1
    while j < n {
      let res2 = residues[j]
      if pos_atoms.length() > 0 {
        let neg2 = get_negative_atoms(res2.resname)
        for a1 in pos_atoms {
          match res1.find_atom(a1) {
            Some(atom1) => {
              for a2 in neg2 {
                match res2.find_atom(a2) {
                  Some(atom2) => {
                    let d = atom1.coord.distance(atom2.coord)
                    if d <= distance_cutoff {
                      result.push(SaltBridge::{
                        resname1: res1.resname, resseq1: res1.resseq, chainid1: res1.chainid,
                        resname2: res2.resname, resseq2: res2.resseq, chainid2: res2.chainid,
                        distance: d, atom_name1: a1, atom_name2: a2,
                      })
                    }
                  }
                  None => ()
                }
              }
            }
            None => ()
          }
        }
      }
      if neg_atoms.length() > 0 {
        let pos2 = get_positive_atoms(res2.resname)
        for a1 in neg_atoms {
          match res1.find_atom(a1) {
            Some(atom1) => {
              for a2 in pos2 {
                match res2.find_atom(a2) {
                  Some(atom2) => {
                    let d = atom1.coord.distance(atom2.coord)
                    if d <= distance_cutoff {
                      result.push(SaltBridge::{
                        resname1: res2.resname, resseq1: res2.resseq, chainid1: res2.chainid,
                        resname2: res1.resname, resseq2: res1.resseq, chainid2: res1.chainid,
                        distance: d, atom_name1: a2, atom_name2: a1,
                      })
                    }
                  }
                  None => ()
                }
              }
            }
            None => ()
          }
        }
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

// ─── π-π Stacking ──────────────────────────────────────────

///|
/// π-π stacking interaction type.
pub enum PiStackingType {
  Parallel
  Antiparallel
  Tilted
}

///|
/// A π-π stacking interaction between aromatic rings.
pub struct PiStacking {
  resname1 : String
  resseq1 : Int
  chainid1 : Char
  resname2 : String
  resseq2 : Int
  chainid2 : Char
  centroid_distance : Double
  ring_angle : Double
  stacking_type : PiStackingType
}

///|
/// Find π-π stacking interactions between aromatic residues.
pub fn find_pi_stacking(
  chain : Chain,
  centroid_cutoff? : Double = 6.0,
  angle_tolerance? : Double = 30.0,
) -> Array[PiStacking] {
  let result : Array[PiStacking] = Array(capacity=0)
  let residues = chain.get_residues()
  let n = residues.length()
  let centroids : Array[Vector3?] = Array::make(n, None)
  let normals : Array[Vector3?] = Array::make(n, None)
  let mut i = 0
  while i < n {
    match get_aromatic_ring_atoms(residues[i].resname) {
      Some(names) => {
        centroids[i] = calc_ring_centroid(residues[i], names)
        normals[i] = calc_ring_normal(residues[i], names)
      }
      None => ()
    }
    i = i + 1
  }
  let mut j_outer = 0
  while j_outer < n {
    let mut k = j_outer + 1
    while k < n {
      match (centroids[j_outer], centroids[k], normals[j_outer], normals[k]) {
        (Some(c1), Some(c2), Some(n1), Some(n2)) => {
          let d = c1.distance(c2)
          if d <= centroid_cutoff {
            let angle = angle_between_normals(n1, n2)
            let stacking_type = if angle < angle_tolerance {
              Parallel
            } else if angle > 180.0 - angle_tolerance {
              Antiparallel
            } else {
              Tilted
            }
            result.push(PiStacking::{
              resname1: residues[j_outer].resname, resseq1: residues[j_outer].resseq, chainid1: residues[j_outer].chainid,
              resname2: residues[k].resname, resseq2: residues[k].resseq, chainid2: residues[k].chainid,
              centroid_distance: d, ring_angle: angle, stacking_type,
            })
          }
        }
        _ => ()
      }
      k = k + 1
    }
    j_outer = j_outer + 1
  }
  result
}

// ─── π-Cation ──────────────────────────────────────────────

///|
/// A π-cation interaction between a cation and an aromatic ring.
pub struct PiCation {
  cation_resname : String
  cation_resseq : Int
  cation_chainid : Char
  aromatic_resname : String
  aromatic_resseq : Int
  aromatic_chainid : Char
  distance : Double
  cation_atom : String
}

///|
/// Find π-cation interactions.
pub fn find_pi_cation(
  chain : Chain,
  distance_cutoff? : Double = 6.0,
) -> Array[PiCation] {
  let result : Array[PiCation] = Array(capacity=0)
  let residues = chain.get_residues()
  let n = residues.length()
  let centroids : Array[Vector3?] = Array::make(n, None)
  let mut i = 0
  while i < n {
    match get_aromatic_ring_atoms(residues[i].resname) {
      Some(names) => { centroids[i] = calc_ring_centroid(residues[i], names) }
      None => ()
    }
    i = i + 1
  }
  let mut j_outer = 0
  while j_outer < n {
    let cation_atoms = get_positive_atoms(residues[j_outer].resname)
    let mut k = 0
    while k < n {
      if k == j_outer { k = k + 1; continue }
      match centroids[k] {
        Some(centroid) => {
          for ca in cation_atoms {
            match residues[j_outer].find_atom(ca) {
              Some(atom) => {
                let d = atom.coord.distance(centroid)
                if d <= distance_cutoff {
                  result.push(PiCation::{
                    cation_resname: residues[j_outer].resname, cation_resseq: residues[j_outer].resseq, cation_chainid: residues[j_outer].chainid,
                    aromatic_resname: residues[k].resname, aromatic_resseq: residues[k].resseq, aromatic_chainid: residues[k].chainid,
                    distance: d, cation_atom: ca,
                  })
                }
              }
              None => ()
            }
          }
        }
        None => ()
      }
      k = k + 1
    }
    j_outer = j_outer + 1
  }
  result
}

// ─── Hydrophobic Interactions ──────────────────────────────

///|
/// A hydrophobic interaction between nonpolar sidechains.
pub struct HydrophobicInteraction {
  resname1 : String
  resseq1 : Int
  chainid1 : Char
  resname2 : String
  resseq2 : Int
  chainid2 : Char
  distance : Double
  atom_name1 : String
  atom_name2 : String
}

///|
/// Find hydrophobic interactions between nonpolar residues.
pub fn find_hydrophobic_interactions(
  chain : Chain,
  distance_cutoff? : Double = 5.0,
) -> Array[HydrophobicInteraction] {
  let result : Array[HydrophobicInteraction] = Array(capacity=0)
  let residues = chain.get_residues()
  let n = residues.length()
  let mut i = 0
  while i < n {
    let res1 = residues[i]
    if not(is_hydrophobic(res1.resname)) { i = i + 1; continue }
    let atoms1 = get_hydrophobic_atoms(res1.resname)
    let mut j = i + 1
    while j < n {
      let res2 = residues[j]
      if not(is_hydrophobic(res2.resname)) { j = j + 1; continue }
      let atoms2 = get_hydrophobic_atoms(res2.resname)
      for a1 in atoms1 {
        match res1.find_atom(a1) {
          Some(atom1) => {
            for a2 in atoms2 {
              match res2.find_atom(a2) {
                Some(atom2) => {
                  let d = atom1.coord.distance(atom2.coord)
                  if d <= distance_cutoff {
                    result.push(HydrophobicInteraction::{
                      resname1: res1.resname, resseq1: res1.resseq, chainid1: res1.chainid,
                      resname2: res2.resname, resseq2: res2.resseq, chainid2: res2.chainid,
                      distance: d, atom_name1: a1, atom_name2: a2,
                    })
                  }
                }
                None => ()
              }
            }
          }
          None => ()
        }
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

// ─── Repulsive Ionic Bonding ───────────────────────────────

///|
/// A repulsive ionic interaction (same-charge residues in proximity).
pub struct RepulsiveIonic {
  resname1 : String
  resseq1 : Int
  chainid1 : Char
  resname2 : String
  resseq2 : Int
  chainid2 : Char
  distance : Double
  atom_name1 : String
  atom_name2 : String
}

///|
/// Find repulsive ionic interactions (same-charge residues close together).
pub fn find_repulsive_ionic(
  chain : Chain,
  distance_cutoff? : Double = 4.0,
) -> Array[RepulsiveIonic] {
  let result : Array[RepulsiveIonic] = Array(capacity=0)
  let residues = chain.get_residues()
  let n = residues.length()
  let mut i = 0
  while i < n {
    let res1 = residues[i]
    let pos1 = get_positive_atoms(res1.resname)
    let neg1 = get_negative_atoms(res1.resname)
    let mut j = i + 1
    while j < n {
      let res2 = residues[j]
      // positive-positive repulsion
      let pos2 = get_positive_atoms(res2.resname)
      for a1 in pos1 {
        match res1.find_atom(a1) {
          Some(atom1) => {
            for a2 in pos2 {
              match res2.find_atom(a2) {
                Some(atom2) => {
                  let d = atom1.coord.distance(atom2.coord)
                  if d <= distance_cutoff {
                    result.push(RepulsiveIonic::{
                      resname1: res1.resname, resseq1: res1.resseq, chainid1: res1.chainid,
                      resname2: res2.resname, resseq2: res2.resseq, chainid2: res2.chainid,
                      distance: d, atom_name1: a1, atom_name2: a2,
                    })
                  }
                }
                None => ()
              }
            }
          }
          None => ()
        }
      }
      // negative-negative repulsion
      let neg2 = get_negative_atoms(res2.resname)
      for a1 in neg1 {
        match res1.find_atom(a1) {
          Some(atom1) => {
            for a2 in neg2 {
              match res2.find_atom(a2) {
                Some(atom2) => {
                  let d = atom1.coord.distance(atom2.coord)
                  if d <= distance_cutoff {
                    result.push(RepulsiveIonic::{
                      resname1: res1.resname, resseq1: res1.resseq, chainid1: res1.chainid,
                      resname2: res2.resname, resseq2: res2.resseq, chainid2: res2.chainid,
                      distance: d, atom_name1: a1, atom_name2: a2,
                    })
                  }
                }
                None => ()
              }
            }
          }
          None => ()
        }
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

// ─── Water Bridges ─────────────────────────────────────────

///|
/// A water bridge: water molecule mediating interaction between two residues.
pub struct WaterBridge {
  resname1 : String
  resseq1 : Int
  chainid1 : Char
  resname2 : String
  resseq2 : Int
  chainid2 : Char
  water_resseq : Int
  water_chainid : Char
  distance1 : Double
  distance2 : Double
}

///|
/// Find water bridges: water molecules H-bonding with two protein residues.
/// Water residues are identified by resname HOH or WAT.
pub fn find_water_bridges(
  chain : Chain,
  hbond_cutoff? : Double = 3.5,
) -> Array[WaterBridge] {
  let result : Array[WaterBridge] = Array(capacity=0)
  let residues = chain.get_residues()
  let n = residues.length()
  // Separate protein and water residues
  let protein_indices = Array(capacity=0)
  let water_indices = Array(capacity=0)
  let mut i = 0
  while i < n {
    let rn = residues[i].resname
    if rn == "HOH" || rn == "WAT" || rn == "TIP" || rn == "SOL" {
      water_indices.push(i)
    } else {
      protein_indices.push(i)
    }
    i = i + 1
  }
  // For each water, find protein residues with O/N within hbond_cutoff
  for wi in water_indices {
    let water = residues[wi]
    let water_o = water.find_atom("O")
    match water_o {
      Some(wo) => {
        // Find all protein residues near this water
        let nearby = Array(capacity=0)
        for pi in protein_indices {
          let res = residues[pi]
          // Check N and O atoms of protein residue
          let n_atom = res.find_atom("N")
          let o_atom = res.find_atom("O")
          match n_atom {
            Some(na) => {
              let d = na.coord.distance(wo.coord)
              if d <= hbond_cutoff { nearby.push((pi, d)) }
            }
            None => ()
          }
          match o_atom {
            Some(oa) => {
              let d = oa.coord.distance(wo.coord)
              if d <= hbond_cutoff { nearby.push((pi, d)) }
            }
            None => ()
          }
        }
        // If water bridges 2+ residues, create bridges for all pairs
        let nn = nearby.length()
        let mut a = 0
        while a < nn {
          let mut b = a + 1
          while b < nn {
            let (pi1, d1) = nearby[a]
            let (pi2, d2) = nearby[b]
            if pi1 != pi2 {
              result.push(WaterBridge::{
                resname1: residues[pi1].resname, resseq1: residues[pi1].resseq, chainid1: residues[pi1].chainid,
                resname2: residues[pi2].resname, resseq2: residues[pi2].resseq, chainid2: residues[pi2].chainid,
                water_resseq: water.resseq, water_chainid: water.chainid,
                distance1: d1, distance2: d2,
              })
            }
            b = b + 1
          }
          a = a + 1
        }
      }
      None => ()
    }
  }
  result
}

// ─── InSty: All-in-one Analysis ────────────────────────────

///|
/// Comprehensive protein interaction analysis result.
pub struct InStyResult {
  salt_bridges : Array[SaltBridge]
  pi_stackings : Array[PiStacking]
  pi_cations : Array[PiCation]
  hydrophobic : Array[HydrophobicInteraction]
  repulsive_ionic : Array[RepulsiveIonic]
  water_bridges : Array[WaterBridge]
}

///|
/// Run all interaction analyses on a chain.
pub fn analyze_interactions(chain : Chain) -> InStyResult {
  InStyResult::{
    salt_bridges: find_salt_bridges(chain),
    pi_stackings: find_pi_stacking(chain),
    pi_cations: find_pi_cation(chain),
    hydrophobic: find_hydrophobic_interactions(chain),
    repulsive_ionic: find_repulsive_ionic(chain),
    water_bridges: find_water_bridges(chain),
  }
}

///|
/// Format an InStyResult as a human-readable summary.
pub fn interaction_summary(result : InStyResult) -> String {
  let lines = Array(capacity=0)
  lines.push("=== Protein Interaction Analysis ===")
  lines.push("")
  lines.push("Salt Bridges: " + result.salt_bridges.length().to_string())
  for sb in result.salt_bridges {
    lines.push(
      "  " + sb.resname1 + sb.resseq1.to_string() + "-" + sb.resname2 + sb.resseq2.to_string() +
      "  d=" + sb.distance.to_string() + "  " + sb.atom_name1 + "..." + sb.atom_name2,
    )
  }
  lines.push("")
  lines.push("π-π Stacking: " + result.pi_stackings.length().to_string())
  for ps in result.pi_stackings {
    let stype = match ps.stacking_type {
      Parallel => "parallel"
      Antiparallel => "antiparallel"
      Tilted => "tilted"
    }
    lines.push(
      "  " + ps.resname1 + ps.resseq1.to_string() + "-" + ps.resname2 + ps.resseq2.to_string() +
      "  d=" + ps.centroid_distance.to_string() + "  angle=" + ps.ring_angle.to_string() + "  " + stype,
    )
  }
  lines.push("")
  lines.push("π-Cation: " + result.pi_cations.length().to_string())
  for pc in result.pi_cations {
    lines.push(
      "  " + pc.cation_resname + pc.cation_resseq.to_string() + "→" +
      pc.aromatic_resname + pc.aromatic_resseq.to_string() +
      "  d=" + pc.distance.to_string() + "  " + pc.cation_atom,
    )
  }
  lines.push("")
  lines.push("Hydrophobic: " + result.hydrophobic.length().to_string())
  for hi in result.hydrophobic {
    lines.push(
      "  " + hi.resname1 + hi.resseq1.to_string() + "-" + hi.resname2 + hi.resseq2.to_string() +
      "  d=" + hi.distance.to_string() + "  " + hi.atom_name1 + "..." + hi.atom_name2,
    )
  }
  lines.push("")
  lines.push("Repulsive Ionic: " + result.repulsive_ionic.length().to_string())
  for ri in result.repulsive_ionic {
    lines.push(
      "  " + ri.resname1 + ri.resseq1.to_string() + "-" + ri.resname2 + ri.resseq2.to_string() +
      "  d=" + ri.distance.to_string(),
    )
  }
  lines.push("")
  lines.push("Water Bridges: " + result.water_bridges.length().to_string())
  for wb in result.water_bridges {
    lines.push(
      "  " + wb.resname1 + wb.resseq1.to_string() + "-" + wb.resname2 + wb.resseq2.to_string() +
      "  via HOH" + wb.water_resseq.to_string() +
      "  d1=" + wb.distance1.to_string() + " d2=" + wb.distance2.to_string(),
    )
  }
  insty_join_strings(lines, "\n")
}

///|
/// Join an array of strings with a separator.
fn insty_join_strings(arr : Array[String], sep : String) -> String {
  let mut result = ""
  let n = arr.length()
  let mut i = 0
  while i < n {
    if i > 0 { result = result + sep }
    result = result + arr[i]
    i = i + 1
  }
  result
}
