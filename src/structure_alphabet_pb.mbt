///|
/// Protein Blocks (PB) structural alphabet — converts protein backbone
/// conformations into a sequence of 16 structural states (a–p).
///
/// Each Protein Block is defined by a characteristic pattern of eight
/// consecutive phi/psi dihedral angles. The 16 PBs were derived from
/// non-redundant protein structure databases and capture the most
/// frequent local backbone conformations.
///
/// Reference: de Brevern et al., "Protein Blocks: a near-unique model
/// for protein backbone structural classes", Proteins 2000.

///|
/// The 16 Protein Block letters.
pub fn protein_blocks_alphabet() -> Array[Char] {
  [
    'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
  ]
}

///|
/// Result of Protein Blocks conversion for a single residue.
pub struct ProteinBlockState {
  residue_index : Int32
  residue_name : String
  block : Char
  rmsd : Double
}

///|
/// Construct a ProteinBlockState.
pub fn ProteinBlockState::new(
  residue_index~ : Int32,
  residue_name~ : String,
  block~ : Char,
  rmsd~ : Double,
) -> ProteinBlockState {
  ProteinBlockState::{ residue_index, residue_name, block, rmsd, }
}

///|
/// Get the Protein Block letter.
pub fn ProteinBlockState::get_block(self : ProteinBlockState) -> Char {
  self.block
}

///|
/// Get the RMSD to the ideal PB pattern.
pub fn ProteinBlockState::get_rmsd(self : ProteinBlockState) -> Double {
  self.rmsd
}

///|
/// Result of converting an entire chain to Protein Blocks.
pub struct ProteinBlocksResult {
  states : Array[ProteinBlockState]
  sequence : String
}

///|
/// Construct a ProteinBlocksResult.
pub fn ProteinBlocksResult::new(
  states~ : Array[ProteinBlockState],
) -> ProteinBlocksResult {
  let mut seq = ""
  for s in states {
    seq = seq + s.block.to_string()
  }
  ProteinBlocksResult::{ states, sequence: seq, }
}

///|
/// Get the PB sequence as a string.
pub fn ProteinBlocksResult::get_sequence(self : ProteinBlocksResult) -> String {
  self.sequence
}

///|
/// Get all PB states.
pub fn ProteinBlocksResult::get_states(
  self : ProteinBlocksResult,
) -> Array[ProteinBlockState] {
  self.states
}

///|
/// Get the length of the PB sequence.
pub fn ProteinBlocksResult::length(self : ProteinBlocksResult) -> Int32 {
  self.states.length()
}

///|
/// Ideal phi/psi patterns for the 16 Protein Blocks.
/// Each PB is defined by 8 consecutive (phi, psi) pairs.
fn pb_ideal_phi(index : Int32) -> Array[Double] {
  match index {
    0 => [-57.0, -57.0, -57.0, -57.0, -57.0, -57.0, -57.0, -57.0]
    1 => [-63.0, -66.0, -67.0, -64.0, -60.0, -57.0, -57.0, -57.0]
    2 => [-57.0, -57.0, -57.0, -57.0, -60.0, -64.0, -67.0, -63.0]
    3 => [-60.0, -60.0, -60.0, -60.0, -60.0, -60.0, -60.0, -60.0]
    4 => [-120.0, -120.0, -120.0, -120.0, -120.0, -120.0, -120.0, -120.0]
    5 => [-100.0, -110.0, -120.0, -120.0, -120.0, -120.0, -120.0, -120.0]
    6 => [-60.0, -90.0, -120.0, -60.0, -60.0, -60.0, -60.0, -60.0]
    7 => [-60.0, -80.0, -100.0, -80.0, -60.0, -60.0, -60.0, -60.0]
    8 => [-60.0, -100.0, -120.0, -120.0, -120.0, -60.0, -60.0, -60.0]
    9 => [-75.0, -75.0, -75.0, -75.0, -75.0, -75.0, -75.0, -75.0]
    10 => [-60.0, -60.0, -60.0, -60.0, -60.0, -60.0, -60.0, -60.0]
    11 => [-60.0, -60.0, -60.0, -60.0, -60.0, -60.0, -60.0, -60.0]
    12 => [-80.0, -70.0, -90.0, -100.0, -80.0, -70.0, -90.0, -100.0]
    13 => [-100.0, -120.0, -60.0, -120.0, -120.0, -120.0, -120.0, -120.0]
    14 => [-65.0, -62.0, -58.0, -55.0, -62.0, -65.0, -62.0, -58.0]
    _ => [-110.0, -125.0, -115.0, -130.0, -110.0, -125.0, -115.0, -130.0]
  }
}

///|
/// Ideal psi patterns for the 16 Protein Blocks.
fn pb_ideal_psi(index : Int32) -> Array[Double] {
  match index {
    0 => [-47.0, -47.0, -47.0, -47.0, -47.0, -47.0, -47.0, -47.0]
    1 => [-41.0, -38.0, -37.0, -40.0, -43.0, -46.0, -47.0, -47.0]
    2 => [-47.0, -47.0, -47.0, -47.0, -43.0, -40.0, -37.0, -41.0]
    3 => [-30.0, -30.0, -30.0, -30.0, -30.0, -30.0, -30.0, -30.0]
    4 => [130.0, 130.0, 130.0, 130.0, 130.0, 130.0, 130.0, 130.0]
    5 => [110.0, 120.0, 130.0, 130.0, 130.0, 130.0, 130.0, 130.0]
    6 => [30.0, 0.0, 130.0, -45.0, -45.0, -45.0, -45.0, -45.0]
    7 => [30.0, 10.0, 120.0, -120.0, -45.0, -45.0, -45.0, -45.0]
    8 => [30.0, 120.0, 130.0, 130.0, 130.0, -45.0, -45.0, -45.0]
    9 => [145.0, 145.0, 145.0, 145.0, 145.0, 145.0, 145.0, 145.0]
    10 => [-120.0, -120.0, -120.0, -120.0, -120.0, -120.0, -120.0, -120.0]
    11 => [120.0, 120.0, 120.0, 120.0, 120.0, 120.0, 120.0, 120.0]
    12 => [10.0, 20.0, -10.0, 30.0, 10.0, 20.0, -10.0, 30.0]
    13 => [130.0, 130.0, 30.0, 130.0, 130.0, 130.0, 130.0, 130.0]
    14 => [-42.0, -44.0, -46.0, -48.0, -44.0, -42.0, -44.0, -46.0]
    _ => [125.0, 135.0, 120.0, 140.0, 125.0, 135.0, 120.0, 140.0]
  }
}

///|
/// Compute the RMSD between observed (phi, psi) angles and an ideal PB pattern.
fn pb_rmsd(
  observed_phi : Array[Double],
  observed_psi : Array[Double],
  ideal_phi : Array[Double],
  ideal_psi : Array[Double],
) -> Double {
  let n = observed_phi.length()
  if n == 0 {
    return 999.0
  }
  let mut sum_sq = 0.0
  for i in 0..<n {
    let dp = pb_angle_diff(observed_phi[i], ideal_phi[i])
    let ds = pb_angle_diff(observed_psi[i], ideal_psi[i])
    sum_sq = sum_sq + dp * dp + ds * ds
  }
  (sum_sq / (n * 2).to_double()).sqrt()
}

///|
/// Compute the smallest difference between two angles (in degrees).
fn pb_angle_diff(a : Double, b : Double) -> Double {
  let mut d = a - b
  while d > 180.0 {
    d = d - 360.0
  }
  while d < -180.0 {
    d = d + 360.0
  }
  d
}

///|
/// Assign the best Protein Block to a window of (phi, psi) angles.
fn assign_protein_block(
  phi_window : Array[Double],
  psi_window : Array[Double],
) -> (Char, Double) {
  let letters = "abcdefghijklmnop"
  let mut best_idx = 0
  let mut best_rmsd = 999.0
  for pb_idx in 0..<16 {
    let ideal_phi = pb_ideal_phi(pb_idx)
    let ideal_psi = pb_ideal_psi(pb_idx)
    let r = pb_rmsd(phi_window, psi_window, ideal_phi, ideal_psi)
    if r < best_rmsd {
      best_rmsd = r
      best_idx = pb_idx
    }
  }
  (letters.unsafe_get(best_idx).unsafe_to_char(), best_rmsd)
}

///|
/// Convert a chain of residues to a Protein Blocks sequence.
pub fn to_protein_blocks(chain : Chain) -> ProteinBlocksResult {
  let residues = chain.get_residues()
  let n = residues.length()
  let states : Array[ProteinBlockState] = []
  if n == 0 {
    return ProteinBlocksResult::new(states=[])
  }
  let coords : Array[(Vector3?, Vector3?, Vector3?)] = []
  for r in residues {
    let n_atom = match r.find_atom("N") {
      Some(a) => Some(a.get_coord())
      None => None
    }
    let ca_atom = match r.find_atom("CA") {
      Some(a) => Some(a.get_coord())
      None => None
    }
    let c_atom = match r.find_atom("C") {
      Some(a) => Some(a.get_coord())
      None => None
    }
    coords.push((n_atom, ca_atom, c_atom))
  }
  let phis : Array[Double] = []
  let psis : Array[Double] = []
  for i in 0..<n {
    let (n_curr, ca_curr, c_curr) = coords[i]
    let n_prev = if i > 0 {
      let (np, _, _) = coords[i - 1]
      np
    } else {
      None
    }
    let c_next = if i < n - 1 {
      let (_, _, cp) = coords[i + 1]
      cp
    } else {
      None
    }
    phis.push(pb_compute_phi(n_prev, ca_curr, c_curr, n_curr))
    psis.push(pb_compute_psi(c_next, n_curr, ca_curr, c_curr))
  }
  let window_size = 8
  for i in 0..<n {
    let phi_window : Array[Double] = []
    let psi_window : Array[Double] = []
    for j in 0..<window_size {
      let idx = i + j - window_size / 2
      if idx >= 0 && idx < n {
        phi_window.push(phis[idx])
        psi_window.push(psis[idx])
      } else {
        phi_window.push(0.0)
        psi_window.push(0.0)
      }
    }
    let (block, rmsd) = assign_protein_block(phi_window, psi_window)
    states.push(
      ProteinBlockState::new(
        residue_index=i,
        residue_name=residues[i].resname,
        block~,
        rmsd~,
      ),
    )
  }
  ProteinBlocksResult::new(states~)
}

///|
/// Compute phi dihedral angle.
fn pb_compute_phi(
  n_prev : Vector3?,
  ca : Vector3?,
  c : Vector3?,
  n_curr : Vector3?,
) -> Double {
  match (n_prev, ca, c, n_curr) {
    (Some(p), Some(a), Some(cc), Some(nc)) => pb_dihedral(p, a, cc, nc)
    _ => 0.0
  }
}

///|
/// Compute psi dihedral angle.
fn pb_compute_psi(
  c_prev : Vector3?,
  n : Vector3?,
  ca : Vector3?,
  c : Vector3?,
) -> Double {
  match (c_prev, n, ca, c) {
    (Some(p), Some(nn), Some(a), Some(cc)) => pb_dihedral(p, nn, a, cc)
    _ => 0.0
  }
}

///|
/// Compute dihedral angle in degrees.
fn pb_dihedral(
  p1 : Vector3,
  p2 : Vector3,
  p3 : Vector3,
  p4 : Vector3,
) -> Double {
  let b1 = Vector3::new(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
  let b2 = Vector3::new(p3.x - p2.x, p3.y - p2.y, p3.z - p2.z)
  let b3 = Vector3::new(p4.x - p3.x, p4.y - p3.y, p4.z - p3.z)
  let n1 = pb_cross(b1, b2)
  let n2 = pb_cross(b2, b3)
  let m1 = pb_veclen(n1)
  let m2 = pb_veclen(n2)
  if m1 < 0.0000000001 || m2 < 0.0000000001 {
    return 0.0
  }
  let u1 = Vector3::new(n1.x / m1, n1.y / m1, n1.z / m1)
  let u2 = Vector3::new(n2.x / m2, n2.y / m2, n2.z / m2)
  let dot = u1.x * u2.x + u1.y * u2.y + u1.z * u2.z
  let mut d = dot
  if d > 1.0 {
    d = 1.0
  }
  if d < -1.0 {
    d = -1.0
  }
  let angle = @math.acos(d)
  let triple = b1.x * n2.x + b1.y * n2.y + b1.z * n2.z
  if triple < 0.0 {
    -angle * 180.0 / 3.141592653589793
  } else {
    angle * 180.0 / 3.141592653589793
  }
}

///|
/// Cross product.
fn pb_cross(a : Vector3, b : Vector3) -> Vector3 {
  Vector3::new(
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x,
  )
}

///|
/// Vector length.
fn pb_veclen(v : Vector3) -> Double {
  (v.x * v.x + v.y * v.y + v.z * v.z).sqrt()
}

///|
/// Convert a Protein Blocks result to a string.
pub fn protein_blocks_to_string(result : ProteinBlocksResult) -> String {
  result.get_sequence()
}
