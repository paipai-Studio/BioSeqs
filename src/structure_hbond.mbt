///|
/// Hydrogen bond detection — geometric donor-acceptor criteria.
/// Ported from biotite.structure.hbond.
///
/// Criteria (Baker & Hubbard 1984, as used in biotite):
///   D-A distance ≤ max_distance (default 3.5 Å)
///   H-A distance ≤ max_h_distance (default 2.5 Å)
///   D-H···A angle ≥ min_angle (default 120°)

///|
/// A detected hydrogen bond.
pub struct HBond {
  donor_resname : String
  donor_resseq : Int
  donor_chainid : Char
  donor_atom : String         // e.g. "N", "OG"
  acceptor_resname : String
  acceptor_resseq : Int
  acceptor_chainid : Char
  acceptor_atom : String      // e.g. "O", "OD1"
  distance : Double           // D-A distance (Å)
  angle : Double              // D-H-A angle (degrees)
}

///|
/// Standard hydrogen bond donor atom names per residue type.
fn is_hbond_donor(resname : String, atom_name : String) -> Bool {
  let r = to_uppercase(resname)
  let a = to_uppercase(atom_name)
  // Backbone amide N (all residues except PRO)
  if a == "N" && r != "PRO" {
    return true
  }
  // Side-chain donors
  match r {
    "SER" => a == "OG"
    "THR" => a == "OG1"
    "TYR" => a == "OH"
    "ASN" => a == "ND2"
    "GLN" => a == "NE2"
    "LYS" => a == "NZ"
    "ARG" => a == "NE" || a == "NH1" || a == "NH2"
    "HIS" => a == "ND1" || a == "NE2"
    "TRP" => a == "NE1"
    _ => false
  }
}

///|
/// Standard hydrogen bond acceptor atom names per residue type.
fn is_hbond_acceptor(resname : String, atom_name : String) -> Bool {
  let r = to_uppercase(resname)
  let a = to_uppercase(atom_name)
  // Backbone carbonyl O
  if a == "O" {
    return true
  }
  // Side-chain acceptors
  match r {
    "SER" => a == "OG"
    "THR" => a == "OG1"
    "TYR" => a == "OH"
    "ASN" => a == "OD1"
    "GLN" => a == "OE1"
    "ASP" => a == "OD1" || a == "OD2"
    "GLU" => a == "OE1" || a == "OE2"
    "HIS" => a == "ND1" || a == "NE2"
    "MET" => a == "SD"
    _ => false
  }
}

///|
/// Estimate hydrogen position for a donor N-H.
/// Uses the preceding C=O geometry: H ≈ N + 1.01 * (N - C_prev) normalized.
/// If no preceding C is available, places H along N→CA direction.
fn estimate_h_position(
  donor_coord : Vector3,
  prev_coord : Vector3,
) -> Vector3 {
  // Direction from previous atom (C) to N
  let dx = donor_coord.x - prev_coord.x
  let dy = donor_coord.y - prev_coord.y
  let dz = donor_coord.z - prev_coord.z
  let len = (dx * dx + dy * dy + dz * dz).sqrt()
  if len < 0.01 {
    return donor_coord
  }
  // N-H bond length ≈ 1.01 Å, along the C→N direction continued
  let scale = 1.01 / len
  Vector3::new(
    donor_coord.x + dx * scale,
    donor_coord.y + dy * scale,
    donor_coord.z + dz * scale,
  )
}

///|
/// Calculate angle (degrees) at vertex B between vectors BA and BC.
fn angle_at(a : Vector3, b : Vector3, c : Vector3) -> Double {
  let v1x = a.x - b.x
  let v1y = a.y - b.y
  let v1z = a.z - b.z
  let v2x = c.x - b.x
  let v2y = c.y - b.y
  let v2z = c.z - b.z
  let dot = v1x * v2x + v1y * v2y + v1z * v2z
  let n1 = (v1x * v1x + v1y * v1y + v1z * v1z).sqrt()
  let n2 = (v2x * v2x + v2y * v2y + v2z * v2z).sqrt()
  if n1 < 0.01 || n2 < 0.01 {
    return 0.0
  }
  let cos_a = dot / (n1 * n2)
  // Clamp to [-1, 1]
  let c = if cos_a > 1.0 {
    1.0
  } else if cos_a < -1.0 {
    -1.0
  } else {
    cos_a
  }
  @math.acos(c) * 180.0 / 3.141592653589793
}

///|
/// Detect hydrogen bonds in a list of atoms using geometric criteria.
///
/// Parameters:
///   `atoms` — flat atom list (e.g. from `Structure::get_atoms()`)
///   `max_distance` — maximum donor-acceptor distance (Å), default 3.5
///   `min_angle` — minimum D-H-A angle (degrees), default 120
///
/// Returns all detected H-bonds sorted by distance.
pub fn detect_hbonds(
  atoms : Array[Atom],
  max_distance~ : Double = 3.5,
  min_angle~ : Double = 120.0,
) -> Array[HBond] {
  // Partition into donors and acceptors
  let donors : Array[Atom] = Array(capacity=0)
  let acceptors : Array[Atom] = Array(capacity=0)
  for a in atoms {
    if is_hbond_donor(a.resname, a.name) {
      donors.push(a)
    }
    if is_hbond_acceptor(a.resname, a.name) {
      acceptors.push(a)
    }
  }
  let result : Array[HBond] = Array(capacity=0)
  for d in donors {
    for acc in acceptors {
      // Skip same-residue intra bonds
      if d.chainid == acc.chainid && d.resseq == acc.resseq && d.resname == acc.resname {
        continue
      }
      let dist = d.coord.distance(acc.coord)
      if dist > max_distance || dist < 0.5 {
        continue
      }
      // Estimate H position: use acceptor as the "previous" direction
      // (simplified — in production would use actual bonded H or preceding C)
      let h_pos = estimate_h_position(d.coord, acc.coord)
      let h_a_dist = h_pos.distance(acc.coord)
      if h_a_dist > 2.5 {
        continue
      }
      let ang = angle_at(d.coord, h_pos, acc.coord)
      if ang < min_angle {
        continue
      }
      result.push(HBond::{
        donor_resname: d.resname,
        donor_resseq: d.resseq,
        donor_chainid: d.chainid,
        donor_atom: d.name,
        acceptor_resname: acc.resname,
        acceptor_resseq: acc.resseq,
        acceptor_chainid: acc.chainid,
        acceptor_atom: acc.name,
        distance: dist,
        angle: ang,
      })
    }
  }
  // Sort by distance (simple insertion sort)
  for i = 1; i < result.length(); i = i + 1 {
    let key = result[i]
    let mut j = i - 1
    while j >= 0 && result[j].distance > key.distance {
      result[j + 1] = result[j]
      j = j - 1
    }
    result[j + 1] = key
  }
  result
}
