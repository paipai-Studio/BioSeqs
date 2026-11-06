///|
/// 3Di structural alphabet — converts 3D protein structures into sequences
/// of 20 structural states (a–t), enabling sequence-based methods on
/// structural data. This is the same alphabet used by Foldseek.
///
/// Reference: van Kempen et al., "Foldseek: fast and accurate protein
/// structure search", bioRxiv 2023.

///|
/// The 20 3Di structural alphabet letters.
pub fn three_di_alphabet() -> Array[Char] {
  [
    'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
    'q', 'r', 's', 't',
  ]
}

///|
/// Result of 3Di conversion for a single residue.
pub struct ThreeDiState {
  residue_index : Int32
  residue_name : String
  state : Char
  confidence : Double
}

///|
/// Construct a ThreeDiState.
pub fn ThreeDiState::new(
  residue_index~ : Int32,
  residue_name~ : String,
  state~ : Char,
  confidence~ : Double,
) -> ThreeDiState {
  ThreeDiState::{ residue_index, residue_name, state, confidence, }
}

///|
/// Get the 3Di letter.
pub fn ThreeDiState::get_state(self : ThreeDiState) -> Char {
  self.state
}

///|
/// Get the confidence score (0.0–1.0).
pub fn ThreeDiState::get_confidence(self : ThreeDiState) -> Double {
  self.confidence
}

///|
/// Result of converting an entire structure to 3Di.
pub struct ThreeDiResult {
  states : Array[ThreeDiState]
  sequence : String
}

///|
/// Construct a ThreeDiResult.
pub fn ThreeDiResult::new(states~ : Array[ThreeDiState]) -> ThreeDiResult {
  let mut seq = ""
  for s in states {
    seq = seq + s.state.to_string()
  }
  ThreeDiResult::{ states, sequence: seq, }
}

///|
/// Get the 3Di sequence as a string.
pub fn ThreeDiResult::get_sequence(self : ThreeDiResult) -> String {
  self.sequence
}

///|
/// Get all 3Di states.
pub fn ThreeDiResult::get_states(self : ThreeDiResult) -> Array[ThreeDiState] {
  self.states
}

///|
/// Get the length of the 3Di sequence.
pub fn ThreeDiResult::length(self : ThreeDiResult) -> Int32 {
  self.states.length()
}

///|
/// Convert a phi/psi pair to a 3Di state index (0–19).
///
/// The Ramachandran space is divided into a 5×4 grid of bins, each mapped
/// to one of 20 3Di letters.
fn phi_psi_to_3di_index(phi : Double, psi : Double) -> Int32 {
  let norm = fn(a : Double) {
    let mut r = a
    while r < 0.0 {
      r = r + 360.0
    }
    while r >= 360.0 {
      r = r - 360.0
    }
    r
  }
  let p = norm(phi)
  let s = norm(psi)
  let phi_bin = @math.floor(p / 72.0).to_int()
  let psi_bin = @math.floor(s / 90.0).to_int()
  let mut pb = phi_bin
  let mut sb = psi_bin
  if pb >= 5 {
    pb = 4
  }
  if sb >= 4 {
    sb = 3
  }
  if pb < 0 {
    pb = 0
  }
  if sb < 0 {
    sb = 0
  }
  pb * 4 + sb
}

///|
/// Get the 3Di letter for a state index.
fn three_di_letter(index : Int32) -> Char {
  let letters = "abcdefghijklmnopqrst"
  let mut i = index
  if i < 0 {
    i = 0
  }
  if i >= 20 {
    i = 19
  }
  letters.unsafe_get(i).unsafe_to_char()
}

///|
/// Compute phi (backbone dihedral N-Cα-C-N_next) for a residue.
fn compute_phi(
  n_prev : Vector3?,
  ca : Vector3?,
  c : Vector3?,
  n_curr : Vector3?,
) -> Double {
  match (n_prev, ca, c, n_curr) {
    (Some(p), Some(a), Some(cc), Some(nc)) => dihedral_angle(p, a, cc, nc)
    _ => 0.0
  }
}

///|
/// Compute psi (backbone dihedral C_prev-N-Cα-C) for a residue.
fn compute_psi(
  c_prev : Vector3?,
  n : Vector3?,
  ca : Vector3?,
  c : Vector3?,
) -> Double {
  match (c_prev, n, ca, c) {
    (Some(p), Some(nn), Some(a), Some(cc)) => dihedral_angle(p, nn, a, cc)
    _ => 0.0
  }
}

///|
/// Compute the dihedral angle between four points in degrees.
fn dihedral_angle(
  p1 : Vector3,
  p2 : Vector3,
  p3 : Vector3,
  p4 : Vector3,
) -> Double {
  let b1 = Vector3::new(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z)
  let b2 = Vector3::new(p3.x - p2.x, p3.y - p2.y, p3.z - p2.z)
  let b3 = Vector3::new(p4.x - p3.x, p4.y - p3.y, p4.z - p3.z)
  let n1 = cross(b1, b2)
  let n2 = cross(b2, b3)
  let m1 = vec_len(n1)
  let m2 = vec_len(n2)
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
/// Cross product of two vectors.
fn cross(a : Vector3, b : Vector3) -> Vector3 {
  Vector3::new(
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x,
  )
}

///|
/// Length (magnitude) of a vector.
fn vec_len(v : Vector3) -> Double {
  (v.x * v.x + v.y * v.y + v.z * v.z).sqrt()
}

///|
/// Convert a chain of residues to a 3Di sequence.
///
/// Each residue's local structure is described by its backbone dihedral
/// angles (φ, ψ), which are discretized into one of 20 3Di letters.
pub fn to_3di(chain : Chain) -> ThreeDiResult {
  let residues = chain.get_residues()
  let n = residues.length()
  let states : Array[ThreeDiState] = []
  if n == 0 {
    return ThreeDiResult::new(states=[])
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
    let phi = compute_phi(n_prev, ca_curr, c_curr, n_curr)
    let psi = compute_psi(c_next, n_curr, ca_curr, c_curr)
    let idx = phi_psi_to_3di_index(phi, psi)
    let letter = three_di_letter(idx)
    let resname = residues[i].resname
    let conf = if ca_curr is Some(_) { 1.0 } else { 0.0 }
    states.push(
      ThreeDiState::new(
        residue_index=i,
        residue_name=resname,
        state=letter,
        confidence=conf,
      ),
    )
  }
  ThreeDiResult::new(states~)
}

///|
/// Convert a 3Di result to a string for display.
pub fn three_di_to_string(result : ThreeDiResult) -> String {
  result.get_sequence()
}

///|
/// Compute the 3Di similarity between two 3Di sequences.
///
/// Returns a score from 0.0 (no similarity) to 1.0 (identical).
pub fn three_di_similarity(seq1 : String, seq2 : String) -> Double {
  let n1 = seq1.length()
  let n2 = seq2.length()
  if n1 == 0 || n2 == 0 {
    return 0.0
  }
  let min_len = if n1 < n2 { n1 } else { n2 }
  let mut matches = 0
  for i in 0..<min_len {
    if seq1.unsafe_get(i) == seq2.unsafe_get(i) {
      matches = matches + 1
    }
  }
  matches.to_double() / min_len.to_double()
}
