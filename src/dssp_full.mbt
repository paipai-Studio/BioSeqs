///|
/// Full DSSP (Dictionary of Secondary Structure of Proteins) algorithm.
/// Pure algorithm implementation — no external mkdssp needed.
///
/// DSSP assigns secondary structure based on hydrogen bond patterns
/// detected using an electrostatic energy formula, plus backbone
/// geometry (phi/psi dihedrals, CA bend angles).
///
/// H-bond energy formula (Kabsch & Sander, 1983):
///   E = q1 * q2 * (1/r(ON) + 1/r(CH) - 1/r(OH) - 1/r(CN)) * f
/// where q1=0.42, q2=0.20, f=332 (kcal/mol conversion),
/// and H position is approximated from N and the preceding C.
///
/// Secondary structure priority: H > G > I > E > B > T > S > -

///|
/// DSSP 8-state secondary structure types.
pub enum DSSPType {
  DSSP_H // H
  DSSP_G // G
  DSSP_I // I
  DSSP_E // E
  DSSP_B // B
  DSSP_T // T
  DSSP_S // S
  DSSP_C // -
} derive(Eq, Debug)
pub extend DSSPType with Eq::{equal, not_equal}
pub extend DSSPType with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn DSSPType::to_char(self : DSSPType) -> Char {
  match self {
    DSSP_H => 'H'
    DSSP_G => 'G'
    DSSP_I => 'I'
    DSSP_E => 'E'
    DSSP_B => 'B'
    DSSP_T => 'T'
    DSSP_S => 'S'
    DSSP_C => '-'
  }
}

///|
pub fn DSSPType::to_str(self : DSSPType) -> String {
  match self {
    DSSP_H => "H"
    DSSP_G => "G"
    DSSP_I => "I"
    DSSP_E => "E"
    DSSP_B => "B"
    DSSP_T => "T"
    DSSP_S => "S"
    DSSP_C => "-"
  }
}

///|
pub fn dssp_type_from_char(c : Char) -> DSSPType {
  match c {
    'H' => DSSP_H
    'G' => DSSP_G
    'I' => DSSP_I
    'E' => DSSP_E
    'B' => DSSP_B
    'T' => DSSP_T
    'S' => DSSP_S
    _ => DSSP_C
  }
}

///|
/// A single DSSP entry for one residue.
pub struct DSSPEntry {
  resseq : Int
  resname : String
  chainid : Char
  ss : DSSPType
  phi : Double
  psi : Double
  asa : Double
  hbond_donors : Array[(Int, Double)]
  hbond_acceptors : Array[(Int, Double)]
}

///|
pub fn DSSPEntry::new(
  resseq~ : Int,
  resname~ : String,
  chainid~ : Char,
  ss~ : DSSPType,
  phi~ : Double = 0.0,
  psi~ : Double = 0.0,
  asa~ : Double = 0.0,
  hbond_donors? : Array[(Int, Double)] = [],
  hbond_acceptors? : Array[(Int, Double)] = [],
) -> DSSPEntry {
  DSSPEntry::{
    resseq, resname, chainid, ss, phi, psi, asa, hbond_donors, hbond_acceptors,
  }
}

// ===== Constants =====
let dssp_hbond_threshold : Double = -0.5
let dssp_q1 : Double = 0.42
let dssp_q2 : Double = 0.20
let dssp_f : Double = 332.0
let dssp_nh_length : Double = 1.01
let dssp_bend_threshold : Double = 70.0

// ===== Vector helpers =====
fn dssp_vec_sub(a : Vector3, b : Vector3) -> Vector3 {
  Vector3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

fn dssp_vec_norm(v : Vector3) -> Double {
  (v.x * v.x + v.y * v.y + v.z * v.z).sqrt()
}

fn dssp_vec_normalize(v : Vector3) -> Vector3 {
  let n = dssp_vec_norm(v)
  if n < 0.0000000001 {
    return Vector3::new(0.0, 0.0, 0.0)
  }
  Vector3::new(v.x / n, v.y / n, v.z / n)
}

fn dssp_vec_dot(a : Vector3, b : Vector3) -> Double {
  a.x * b.x + a.y * b.y + a.z * b.z
}

fn dssp_vec_cross(a : Vector3, b : Vector3) -> Vector3 {
  Vector3::new(
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x,
  )
}

fn dssp_inv_dist(a : Vector3, b : Vector3) -> Double {
  let d = a.distance(b)
  if d < 0.000001 {
    return 0.0
  }
  1.0 / d
}

// ===== H-bond Energy =====

///|
/// Approximate the H atom position on residue j.
/// H_j ≈ N_j + 1.01 * (N_j - C_{j-1}) / |N_j - C_{j-1}|
pub fn dssp_approx_h_position(
  n_j : Vector3,
  ca_j : Vector3,
  prev_c : Vector3?,
) -> Vector3 {
  match prev_c {
    Some(cp) => {
      let dir = dssp_vec_sub(n_j, cp)
      let unit = dssp_vec_normalize(dir)
      Vector3::new(
        n_j.x + dssp_nh_length * unit.x,
        n_j.y + dssp_nh_length * unit.y,
        n_j.z + dssp_nh_length * unit.z,
      )
    }
    None => {
      let dir = dssp_vec_sub(n_j, ca_j)
      let unit = dssp_vec_normalize(dir)
      Vector3::new(
        n_j.x + dssp_nh_length * unit.x,
        n_j.y + dssp_nh_length * unit.y,
        n_j.z + dssp_nh_length * unit.z,
      )
    }
  }
}

///|
/// Compute the DSSP hydrogen bond energy.
/// E = q1 * q2 * (1/r(ON) + 1/r(CH) - 1/r(OH) - 1/r(CN)) * f
pub fn dssp_hbond_energy(
  o_i : Vector3,
  c_i : Vector3,
  n_j : Vector3,
  h_j : Vector3,
) -> Double {
  let inv_on = dssp_inv_dist(o_i, n_j)
  let inv_ch = dssp_inv_dist(c_i, h_j)
  let inv_oh = dssp_inv_dist(o_i, h_j)
  let inv_cn = dssp_inv_dist(c_i, n_j)
  dssp_q1 * dssp_q2 * (inv_on + inv_ch - inv_oh - inv_cn) * dssp_f
}

///|
/// Compute H-bond energy between two residues.
/// Residue i donates C=O, residue j accepts with N-H.
pub fn dssp_residue_hbond_energy(
  o_i : Vector3,
  c_i : Vector3,
  n_j : Vector3,
  ca_j : Vector3,
  prev_c_j : Vector3?,
) -> Double {
  let h_j = dssp_approx_h_position(n_j, ca_j, prev_c_j)
  dssp_hbond_energy(o_i, c_i, n_j, h_j)
}

// ===== Backbone extraction =====

///|
pub struct BackboneAtoms {
  n : Vector3
  ca : Vector3
  c : Vector3
  o : Vector3
  prev_c : Vector3?
  has_atoms : Bool
}

///|
fn extract_backbone(residues : Array[Residue]) -> Array[BackboneAtoms] {
  let n = residues.length()
  let result : Array[BackboneAtoms] = Array(capacity=n)
  let mut prev_c : Vector3? = None
  for i = 0; i < n; i = i + 1 {
    let res = residues[i]
    let n_opt = res.find_atom("N")
    let ca_opt = res.find_atom("CA")
    let c_opt = res.find_atom("C")
    let o_opt = res.find_atom("O")
    match (n_opt, ca_opt, c_opt, o_opt) {
      (Some(na), Some(ca), Some(c), Some(o)) => {
        result.push(
          BackboneAtoms::{
            n: na.coord, ca: ca.coord, c: c.coord, o: o.coord, prev_c, has_atoms: true,
          },
        )
        prev_c = Some(c.coord)
      }
      _ => {
        let ph = Vector3::new(0.0, 0.0, 0.0)
        result.push(BackboneAtoms::{ n: ph, ca: ph, c: ph, o: ph, prev_c, has_atoms: false })
        match c_opt {
          Some(c) => prev_c = Some(c.coord)
          None => ()
        }
      }
    }
  }
  result
}

// ===== H-bond detection =====

///|
pub struct DSSPHBond {
  donor_idx : Int
  acceptor_idx : Int
  energy : Double
}

///|
/// Find all hydrogen bonds using the DSSP energy criterion (E < -0.5).
pub fn dssp_find_hbonds(backbones : Array[BackboneAtoms]) -> Array[DSSPHBond] {
  let n = backbones.length()
  let result : Array[DSSPHBond] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    if not(backbones[i].has_atoms) {
      continue
    }
    for j = 0; j < n; j = j + 1 {
      let diff = i - j
      let abs_diff = if diff < 0 { -diff } else { diff }
      if abs_diff < 2 {
        continue
      }
      if not(backbones[j].has_atoms) {
        continue
      }
      let e = dssp_residue_hbond_energy(
        backbones[i].o, backbones[i].c,
        backbones[j].n, backbones[j].ca, backbones[j].prev_c,
      )
      if e < dssp_hbond_threshold {
        result.push(DSSPHBond::{ donor_idx: i, acceptor_idx: j, energy: e })
      }
    }
  }
  result
}

fn build_donor_map(hbonds : Array[DSSPHBond]) -> Map[Int, Array[(Int, Double)]] {
  let m : Map[Int, Array[(Int, Double)]] = Map([], capacity=16)
  for hb in hbonds {
    match m.get(hb.donor_idx) {
      Some(arr) => arr.push((hb.acceptor_idx, hb.energy))
      None => {
        let arr : Array[(Int, Double)] = Array(capacity=0)
        arr.push((hb.acceptor_idx, hb.energy))
        m[hb.donor_idx] = arr
      }
    }
  }
  m
}

fn build_acceptor_map(hbonds : Array[DSSPHBond]) -> Map[Int, Array[(Int, Double)]] {
  let m : Map[Int, Array[(Int, Double)]] = Map([], capacity=16)
  for hb in hbonds {
    match m.get(hb.acceptor_idx) {
      Some(arr) => arr.push((hb.donor_idx, hb.energy))
      None => {
        let arr : Array[(Int, Double)] = Array(capacity=0)
        arr.push((hb.donor_idx, hb.energy))
        m[hb.acceptor_idx] = arr
      }
    }
  }
  m
}

fn has_hbond(donor_map : Map[Int, Array[(Int, Double)]], i : Int, j : Int) -> Bool {
  match donor_map.get(i) {
    Some(arr) => {
      for pair in arr {
        if pair.0 == j {
          return true
        }
      }
      false
    }
    None => false
  }
}

///|
/// Check if residue i has an n-turn: H-bond from C=O(i) to N-H(i+n).
pub fn dssp_has_n_turn(
  donor_map : Map[Int, Array[(Int, Double)]],
  i : Int,
  n : Int,
) -> Bool {
  has_hbond(donor_map, i, i + n)
}

// ===== β-bridge detection =====

///|
pub struct BetaBridge {
  i : Int
  j : Int
  bridge_type : Char // 'E' = antiparallel, 'P' = parallel
}

fn bridge_key(i : Int, j : Int) -> String {
  let a = if i < j { i } else { j }
  let b = if i < j { j } else { i }
  a.to_string() + "_" + b.to_string()
}

///|
/// Find all β-bridges.
/// Antiparallel (E): H-bond(i→j) AND H-bond(j→i)
/// Parallel (P): H-bond(i→j) AND H-bond((j-1)→(i+1))
pub fn dssp_find_bridges(
  donor_map : Map[Int, Array[(Int, Double)]],
  n_residues : Int,
) -> Array[BetaBridge] {
  let result : Array[BetaBridge] = Array(capacity=0)
  let seen : Map[String, Unit] = Map([], capacity=16)
  for i = 0; i < n_residues; i = i + 1 {
    for j = 0; j < n_residues; j = j + 1 {
      let diff = i - j
      let abs_diff = if diff < 0 { -diff } else { diff }
      if abs_diff < 3 {
        continue
      }
      let hb_ij = has_hbond(donor_map, i, j)
      let hb_ji = has_hbond(donor_map, j, i)
      if hb_ij && hb_ji {
        let key = bridge_key(i, j)
        if not(seen.contains(key)) {
          seen[key] = ()
          result.push(BetaBridge::{ i, j, bridge_type: 'E' })
        }
      }
      if hb_ij && j > 0 && i + 1 < n_residues {
        if has_hbond(donor_map, j - 1, i + 1) {
          let key = bridge_key(i, j)
          if not(seen.contains(key)) {
            seen[key] = ()
            result.push(BetaBridge::{ i, j, bridge_type: 'P' })
          }
        }
      }
      if hb_ji && i > 0 && j + 1 < n_residues {
        if has_hbond(donor_map, i - 1, j + 1) {
          let key = bridge_key(i, j)
          if not(seen.contains(key)) {
            seen[key] = ()
            result.push(BetaBridge::{ i, j, bridge_type: 'P' })
          }
        }
      }
    }
  }
  result
}

fn build_bridge_map(bridges : Array[BetaBridge]) -> Map[Int, Array[BetaBridge]] {
  let m : Map[Int, Array[BetaBridge]] = Map([], capacity=16)
  for b in bridges {
    match m.get(b.i) {
      Some(arr) => arr.push(b)
      None => {
        let arr : Array[BetaBridge] = Array(capacity=0)
        arr.push(b)
        m[b.i] = arr
      }
    }
    match m.get(b.j) {
      Some(arr) => arr.push(b)
      None => {
        let arr : Array[BetaBridge] = Array(capacity=0)
        arr.push(b)
        m[b.j] = arr
      }
    }
  }
  m
}

// ===== Bend detection =====

///|
/// Detect bends: residue i is a bend if the angle between
/// CA(i-2)→CA(i) and CA(i-2)→CA(i+2) exceeds 70°.
pub fn dssp_find_bends(ca_coords : Array[Vector3]) -> Array[Bool] {
  let n = ca_coords.length()
  let bends : Array[Bool] = Array::make(n, false)
  for i = 0; i < n; i = i + 1 {
    if i < 2 || i + 2 >= n {
      continue
    }
    let ca_prev2 = ca_coords[i - 2]
    let ca_i = ca_coords[i]
    let ca_next2 = ca_coords[i + 2]
    let v1 = dssp_vec_sub(ca_i, ca_prev2)
    let v2 = dssp_vec_sub(ca_next2, ca_prev2)
    let n1 = dssp_vec_norm(v1)
    let n2 = dssp_vec_norm(v2)
    if n1 < 0.000001 || n2 < 0.000001 {
      continue
    }
    let cos_angle = dssp_vec_dot(v1, v2) / (n1 * n2)
    let clamped = if cos_angle > 1.0 {
      1.0
    } else if cos_angle < -1.0 {
      -1.0
    } else {
      cos_angle
    }
    let angle_deg = @math.acos(clamped) * 180.0 / 3.14159265358979323846
    if angle_deg > dssp_bend_threshold {
      bends[i] = true
    }
  }
  bends
}

// ===== Dihedral angle computation =====

fn dssp_calc_dihedral(a : Vector3, b : Vector3, c : Vector3, d : Vector3) -> Double {
  let b1 = dssp_vec_sub(b, a)
  let b2 = dssp_vec_sub(c, b)
  let b3 = dssp_vec_sub(d, c)
  let n1 = dssp_vec_cross(b1, b2)
  let n2 = dssp_vec_cross(b2, b3)
  let m1 = dssp_vec_cross(n1, dssp_vec_normalize(b2))
  let x = dssp_vec_dot(n1, n2)
  let y = dssp_vec_dot(m1, n2)
  @math.atan2(y, x) * 180.0 / 3.14159265358979323846
}

///|
/// Compute phi/psi dihedrals for each residue from backbone atoms.
/// Returns Array[(phi, psi)] in degrees; 0.0 for terminal/missing.
pub fn dssp_calc_dihedrals(backbones : Array[BackboneAtoms]) -> Array[(Double, Double)] {
  let n = backbones.length()
  let result : Array[(Double, Double)] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let phi = if i > 0 && backbones[i - 1].has_atoms && backbones[i].has_atoms {
      dssp_calc_dihedral(
        backbones[i - 1].c, backbones[i].n, backbones[i].ca, backbones[i].c,
      )
    } else {
      0.0
    }
    let psi = if i + 1 < n && backbones[i].has_atoms && backbones[i + 1].has_atoms {
      dssp_calc_dihedral(
        backbones[i].n, backbones[i].ca, backbones[i].c, backbones[i + 1].n,
      )
    } else {
      0.0
    }
    result.push((phi, psi))
  }
  result
}

// ===== Main DSSP assignment =====

///|
/// Run the full DSSP secondary structure assignment on a chain.
///
/// Algorithm:
/// 1. Extract backbone atoms (N, CA, C, O) for each residue
/// 2. Compute all H-bonds using the DSSP energy formula
/// 3. Detect n-turns (n=3,4,5) from H-bond patterns
/// 4. Assign helices (H, G, I) from consecutive n-turns
/// 5. Detect β-bridges and assign strands (E) / bridges (B)
/// 6. Assign turns (T) for remaining n-turn residues
/// 7. Detect bends (S) from CA geometry
/// 8. Remaining residues are coil (-)
///
/// Priority: H > G > I > E > B > T > S > -
pub fn dssp_assign(chain : Chain) -> Array[DSSPEntry] {
  let residues = chain.get_residues()
  let n = residues.length()
  if n == 0 {
    return []
  }
  let backbones = extract_backbone(residues)
  let hbonds = dssp_find_hbonds(backbones)
  let donor_map = build_donor_map(hbonds)
  let acceptor_map = build_acceptor_map(hbonds)
  let dihedrals = dssp_calc_dihedrals(backbones)

  let ca_coords : Array[Vector3] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    ca_coords.push(backbones[i].ca)
  }
  let bends = dssp_find_bends(ca_coords)

  // Detect n-turns
  let turn3 : Array[Bool] = Array::make(n, false)
  let turn4 : Array[Bool] = Array::make(n, false)
  let turn5 : Array[Bool] = Array::make(n, false)
  for i = 0; i < n; i = i + 1 {
    if dssp_has_n_turn(donor_map, i, 3) {
      turn3[i] = true
    }
    if dssp_has_n_turn(donor_map, i, 4) {
      turn4[i] = true
    }
    if dssp_has_n_turn(donor_map, i, 5) {
      turn5[i] = true
    }
  }

  // Assign SS types (priority: H > G > I > E > B > T > S > -)
  let ss : Array[DSSPType] = Array::make(n, DSSP_C)

  // α-helix (H): 4 consecutive 4-turns
  for i = 0; i < n; i = i + 1 {
    if i >= 3 && turn4[i] && turn4[i - 1] && turn4[i - 2] && turn4[i - 3] {
      ss[i] = DSSP_H
    }
  }

  // 3₁₀-helix (G): 2 consecutive 3-turns, not already H
  for i = 0; i < n; i = i + 1 {
    if i >= 1 && turn3[i] && turn3[i - 1] && ss[i] == DSSP_C {
      ss[i] = DSSP_G
    }
  }

  // π-helix (I): 5 consecutive 5-turns, not already H or G
  for i = 0; i < n; i = i + 1 {
    if i >= 4 && turn5[i] && turn5[i - 1] && turn5[i - 2] && turn5[i - 3] && turn5[i - 4] && ss[i] == DSSP_C {
      ss[i] = DSSP_I
    }
  }

  // β-bridges and strands
  let bridges = dssp_find_bridges(donor_map, n)
  let _ = build_bridge_map(bridges)
  let in_bridge : Array[Bool] = Array::make(n, false)
  for b in bridges {
    in_bridge[b.i] = true
    in_bridge[b.j] = true
  }
  for i = 0; i < n; i = i + 1 {
    if in_bridge[i] && ss[i] == DSSP_C {
      let prev_in = i > 0 && in_bridge[i - 1]
      let next_in = i + 1 < n && in_bridge[i + 1]
      if prev_in || next_in {
        ss[i] = DSSP_E
      } else {
        ss[i] = DSSP_B
      }
    }
  }

  // Turn (T): has an n-turn but not assigned to H, G, I, E, B
  for i = 0; i < n; i = i + 1 {
    if ss[i] == DSSP_C && (turn3[i] || turn4[i] || turn5[i]) {
      ss[i] = DSSP_T
    }
  }

  // Bend (S): bend detected and not already assigned
  for i = 0; i < n; i = i + 1 {
    if ss[i] == DSSP_C && bends[i] {
      ss[i] = DSSP_S
    }
  }

  // Build DSSP entries
  let entries : Array[DSSPEntry] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let res = residues[i]
    let (phi, psi) = dihedrals[i]
    let donors : Array[(Int, Double)] = match donor_map.get(i) {
      Some(arr) => arr
      None => []
    }
    let acceptors : Array[(Int, Double)] = match acceptor_map.get(i) {
      Some(arr) => arr
      None => []
    }
    entries.push(
      DSSPEntry::new(
        resseq=res.resseq,
        resname=res.resname,
        chainid=res.chainid,
        ss=ss[i],
        phi=phi,
        psi=psi,
        hbond_donors=donors,
        hbond_acceptors=acceptors,
      ),
    )
  }
  entries
}

// ===== Summary and output =====

///|
/// Count residues by DSSP type.
pub fn dssp_summary(entries : Array[DSSPEntry]) -> Map[String, Int] {
  let result : Map[String, Int] = Map([], capacity=8)
  for entry in entries {
    let key = entry.ss.to_str()
    match result.get(key) {
      Some(count) => result[key] = count + 1
      None => result[key] = 1
    }
  }
  result
}

///|
/// Get the secondary structure as a single-letter string sequence.
pub fn dssp_ss_sequence(entries : Array[DSSPEntry]) -> String {
  let mut result = ""
  for entry in entries {
    result = result + entry.ss.to_str()
  }
  result
}

///|
/// Count the number of hydrogen bonds found.
pub fn dssp_count_hbonds(entries : Array[DSSPEntry]) -> Int {
  let mut count = 0
  for entry in entries {
    count = count + entry.hbond_donors.length()
  }
  count
}

///|
/// Convert DSSP entries to a summary text string.
pub fn dssp_to_string(entries : Array[DSSPEntry]) -> String {
  let mut output = "DSSP Secondary Structure Assignment\n"
  output = output + "  #  RESIDUE AA SS    PHI     PSI\n"
  for i = 0; i < entries.length(); i = i + 1 {
    let e = entries[i]
    let idx = dssp_pad_left((i + 1).to_string(), 5)
    let resseq = dssp_pad_left(e.resseq.to_string(), 5)
    let resname = e.resname
    let ss = e.ss.to_str()
    let phi = dssp_pad_left(e.phi.to_string(), 7)
    let psi = dssp_pad_left(e.psi.to_string(), 7)
    output = output + idx + resseq + "  " + resname + "   " + ss + "  " + phi + "  " + psi + "\n"
  }
  output
}

fn dssp_pad_left(s : String, width : Int) -> String {
  let len = s.length()
  if len >= width {
    return s
  }
  let mut result = ""
  for _ in 0..<(width - len) {
    result = result + " "
  }
  result + s
}

// ===== Structure generators for testing =====

///|
/// Create an ideal α-helix structure with n residues.
/// Uses standard α-helix backbone geometry:
///   rise = 1.5 Å/residue, turn = 100°/residue, radius = 2.3 Å
pub fn create_ideal_helix(n : Int) -> Structure {
  let residues : Array[Residue] = Array(capacity=n)
  let rise = 1.5
  let turn = 100.0
  let radius = 2.3
  let pi = 3.14159265358979323846
  for i = 0; i < n; i = i + 1 {
    let angle = turn * i.to_double() * pi / 180.0
    let z = rise * i.to_double()
    let ca = Vector3::new(radius * @math.cos(angle), radius * @math.sin(angle), z)
    let n_angle = angle - 0.5
    let n_z = z - 0.8
    let n_pos = Vector3::new(radius * @math.cos(n_angle) * 0.95, radius * @math.sin(n_angle) * 0.95, n_z)
    let c_angle = angle + 0.6
    let c_z = z + 0.5
    let c_pos = Vector3::new(radius * @math.cos(c_angle), radius * @math.sin(c_angle), c_z)
    let o_angle = c_angle + 2.5
    let o_z = c_z + 0.3
    let o_pos = Vector3::new(radius * @math.cos(o_angle) * 1.1, radius * @math.sin(o_angle) * 1.1, o_z)
    let resname = if i == 0 {
      "ALA"
    } else if i == 1 {
      "GLY"
    } else if i == 2 {
      "VAL"
    } else {
      "LEU"
    }
    residues.push(
      Residue::new(
        resname=resname, chainid='A', resseq=i + 1,
        atoms=[
          Atom::new(name="N", coord=n_pos, resname=resname, chainid='A', resseq=i + 1),
          Atom::new(name="CA", coord=ca, resname=resname, chainid='A', resseq=i + 1),
          Atom::new(name="C", coord=c_pos, resname=resname, chainid='A', resseq=i + 1),
          Atom::new(name="O", coord=o_pos, resname=resname, chainid='A', resseq=i + 1),
        ],
      ),
    )
  }
  let chain = Chain::new(id='A', residues=residues)
  let model = Model::new(id=0, chains=[chain])
  Structure::new(id="ideal_helix", models=[model])
}

///|
/// Create an ideal antiparallel β-hairpin structure.
/// Two strands of n_strand residues each, connected by a 2-residue turn.
pub fn create_beta_hairpin(n_strand : Int) -> Structure {
  let total = n_strand * 2 + 2
  let residues : Array[Residue] = Array(capacity=total)
  let spacing = 3.5
  for i = 0; i < total; i = i + 1 {
    let (ca, n_pos, c_pos, o_pos) = if i < n_strand {
      let x = spacing * i.to_double()
      (Vector3::new(x, 0.0, 0.0),
       Vector3::new(x - 1.2, 0.8, 0.0),
       Vector3::new(x + 1.2, 0.8, 0.0),
       Vector3::new(x + 1.0, 2.0, 0.0))
    } else if i < n_strand + 2 {
      let x = spacing * n_strand.to_double() + spacing * 0.5 * (i - n_strand).to_double()
      let y = 2.0 + (i - n_strand).to_double() * 1.5
      (Vector3::new(x, y, 0.0),
       Vector3::new(x - 1.0, y + 0.5, 0.0),
       Vector3::new(x + 1.0, y + 0.5, 0.0),
       Vector3::new(x + 0.5, y + 1.8, 0.0))
    } else {
      let j = i - n_strand - 2
      let x = spacing * (n_strand - 1 - j).to_double()
      (Vector3::new(x, 5.0, 0.0),
       Vector3::new(x + 1.2, 5.8, 0.0),
       Vector3::new(x - 1.2, 5.8, 0.0),
       Vector3::new(x - 1.0, 7.0, 0.0))
    }
    let resname = "ALA"
    residues.push(
      Residue::new(
        resname=resname, chainid='A', resseq=i + 1,
        atoms=[
          Atom::new(name="N", coord=n_pos, resname=resname, chainid='A', resseq=i + 1),
          Atom::new(name="CA", coord=ca, resname=resname, chainid='A', resseq=i + 1),
          Atom::new(name="C", coord=c_pos, resname=resname, chainid='A', resseq=i + 1),
          Atom::new(name="O", coord=o_pos, resname=resname, chainid='A', resseq=i + 1),
        ],
      ),
    )
  }
  let chain = Chain::new(id='A', residues=residues)
  let model = Model::new(id=0, chains=[chain])
  Structure::new(id="beta_hairpin", models=[model])
}

///|
/// Create a simple extended chain (all coil) with n residues.
pub fn create_extended_chain(n : Int) -> Structure {
  let residues : Array[Residue] = Array(capacity=n)
  let spacing = 3.8
  for i = 0; i < n; i = i + 1 {
    let x = spacing * i.to_double()
    let ca = Vector3::new(x, 0.0, 0.0)
    let n_pos = Vector3::new(x - 1.2, 0.8, 0.0)
    let c_pos = Vector3::new(x + 1.2, 0.8, 0.0)
    let o_pos = Vector3::new(x + 1.0, 2.0, 0.0)
    let resname = "ALA"
    residues.push(
      Residue::new(
        resname=resname, chainid='A', resseq=i + 1,
        atoms=[
          Atom::new(name="N", coord=n_pos, resname=resname, chainid='A', resseq=i + 1),
          Atom::new(name="CA", coord=ca, resname=resname, chainid='A', resseq=i + 1),
          Atom::new(name="C", coord=c_pos, resname=resname, chainid='A', resseq=i + 1),
          Atom::new(name="O", coord=o_pos, resname=resname, chainid='A', resseq=i + 1),
        ],
      ),
    )
  }
  let chain = Chain::new(id='A', residues=residues)
  let model = Model::new(id=0, chains=[chain])
  Structure::new(id="extended", models=[model])
}
