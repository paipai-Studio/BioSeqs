///|
/// structure_internal_coords.mbt — Backbone internal coordinates
/// forward ↔ Cartesian coordinate reconstruction.
///
/// Inspired by Biopython's Bio.PDB.internal_coords, focused on the
/// practical subset: backbone φ/ψ/ω → Cartesian coordinates.
///
/// Uses Engh & Huber standard bond lengths / angles and a
/// peptide-plane + chained-rotation reconstruction method.
///
/// NOTE: Vector3 uses .sub(a), .add(a), .scale(s), .dot(a), .cross(a),
/// .norm(), .normalize(), .negate() — no + - * operators.

// ============================================================================
// Standard geometry (Engh & Huber 1991)
// ============================================================================

pub const BOND_N_CA : Double = 1.458 // Å

///|
pub const BOND_CA_C : Double = 1.523 // Å

///|
pub const BOND_C_N : Double = 1.329 // Å (peptide bond)

///|
pub const BOND_C_O : Double = 1.231 // Å

///|
pub const ANGLE_N_CA_C : Double = 1.9408 // 111.2° rad

///|
pub const ANGLE_CA_C_N : Double = 2.0280 // 116.2°

///|
pub const ANGLE_C_N_CA : Double = 2.1240 // 121.7°

///|
pub const ANGLE_C_C_O : Double = 2.1065 // 120.7° (C=O in peptide plane)

///|
pub const DEG2RAD : Double = 0.017453292519943295

// ============================================================================
// Data structures
// ============================================================================

///|
/// Backbone dihedral angles for one residue, degrees.
/// phi   = C(i-1)-N(i)-CA(i)-C(i)      (None for i=0)
/// psi   = N(i)-CA(i)-C(i)-N(i+1)      (None for last residue)
/// omega = CA(i)-C(i)-N(i+1)-CA(i+1)   (peptide bond, ≈180° trans)
pub struct BackboneAngles {
  phi : Double?
  psi : Double?
  omega : Double?
}

///|
/// Per-residue backbone angles for a full chain.
pub struct BackboneAngleChain {
  angles : Array[BackboneAngles]
}

///|
pub fn BackboneAngleChain::new(
  angles : Array[BackboneAngles],
) -> BackboneAngleChain {
  BackboneAngleChain::{ angles, }
}

///|
pub fn BackboneAngleChain::len(self : BackboneAngleChain) -> Int32 {
  self.angles.length()
}

// ============================================================================
// Forward: extract backbone angles from an existing AtomArray
// ============================================================================

///|
/// Extract all backbone dihedral angles from a single-chain AtomArray.
pub fn AtomArray::get_backbone_angles(
  self : AtomArray,
) -> Array[BackboneAngles] {
  let res_starts = self.get_residue_starts()
  let n_res = res_starts.length()
  let out : Array[BackboneAngles] = Array(capacity=n_res)
  let mut i = 0
  while i < n_res {
    let phi = self.backbone_phi(i)
    let psi = self.backbone_psi(i)
    let omega = if i < n_res - 1 { self.backbone_omega(i) } else { None }
    out.push(BackboneAngles::{ phi, psi, omega, })
    i = i + 1
  }
  out
}

// ============================================================================
// Inverse geometry helpers
// ============================================================================

///|
/// Canonical placement of the first 3 backbone atoms.
/// N0 at origin, N0-CA0 along +x, CA0-C0 in xy-plane below x-axis
/// at angle (π - ANGLE_N_CA_C) from CA0→N0.
fn initial_backbone_triad() -> (Vector3, Vector3, Vector3) {
  let n0 = Vector3::new(0.0, 0.0, 0.0)
  let ca0 = Vector3::new(BOND_N_CA, 0.0, 0.0)
  let t = @math.PI - ANGLE_N_CA_C
  let c0_local = Vector3::new(
    BOND_CA_C * @math.cos(t),
    -BOND_CA_C * @math.sin(t),
    0.0,
  )
  let c0 = ca0.add(c0_local)
  (n0, ca0, c0)
}

///|
/// Build an orthonormal right-handed frame from two non-collinear vectors.
/// Returns (origin, e1, e2, e3) where e1=normalize(v1),
/// e2 perpendicular to e1 in the plane of v1 & v2, e3 = e1×e2.
fn build_frame(
  origin : Vector3,
  v1 : Vector3,
  v2 : Vector3,
) -> (Vector3, Vector3, Vector3, Vector3) {
  let e1 = v1.normalize()
  let e3 = e1.cross(v2).normalize()
  let e2 = e3.cross(e1).normalize()
  (origin, e1, e2, e3)
}

///|
/// Express global point p in local frame (origin, e1, e2, e3).
fn to_local(
  p : Vector3,
  origin : Vector3,
  e1 : Vector3,
  e2 : Vector3,
  e3 : Vector3,
) -> Vector3 {
  let d = p.sub(origin)
  let x = d.dot(e1)
  let y = d.dot(e2)
  let z = d.dot(e3)
  Vector3::new(x, y, z)
}

///|
/// Map a local frame point back to global coordinates.
fn to_global(
  lcl : Vector3,
  origin : Vector3,
  e1 : Vector3,
  e2 : Vector3,
  e3 : Vector3,
) -> Vector3 {
  Vector3::new(
    origin.x + lcl.x * e1.x + lcl.y * e2.x + lcl.z * e3.x,
    origin.y + lcl.x * e1.y + lcl.y * e2.y + lcl.z * e3.y,
    origin.z + lcl.x * e1.z + lcl.y * e2.z + lcl.z * e3.z,
  )
}

///|
/// Given 3 previous backbone atoms (C_prev, N_cur, CA_prev) and the
/// peptide bond angle ω (C_prev-N_cur-CA_cur-CA_prev ≈ 180°),
/// place CA_cur in the peptide plane.
///
/// Strategy: build a local frame at N_cur from (C_prev→N_cur) and
/// (N_cur→CA_prev). In this frame the peptide plane is z=0, so
/// CA_cur is placed at (BOND_N_CA * cos α, BOND_N_CA * sin α, 0)
/// where α is determined by φ + ω geometry.
///
/// The `phi_deg` argument is the dihedral C_prev-N_cur-CA_cur-C_cur.
/// We use it as a rotation out of peptide plane for CA_cur.
fn place_ca_cur(
  c_prev : Vector3,
  n_cur : Vector3,
  ca_prev : Vector3,
  phi_deg : Double,
) -> Vector3 {
  // Frame at N_cur
  let v_nc = n_cur.sub(c_prev) // C_prev → N_cur
  let v_ncaprev = ca_prev.sub(n_cur) // N_cur → CA_prev
  let (_, e1, e2, e3) = build_frame(n_cur, v_nc, v_ncaprev)
  // e1 = C_prev→N_cur direction (peptide plane normal e3)
  // In the peptide plane, N_cur→CA_cur makes angle with C_prev→N_cur
  // determined by peptide bond geometry. For simplicity we align it
  // with the peptide bond continuation direction and let φ rotate it.
  // Local direction in plane = (cos φ, sin φ, 0)
  let phi = phi_deg * DEG2RAD
  let lcl = Vector3::new(@math.cos(phi), @math.sin(phi), 0.0)
  to_global(lcl.scale(BOND_N_CA), n_cur, e1, e2, e3)
}

///|
/// Place C_cur given N_cur, CA_cur and ψ(N_cur-CA_cur-C_cur-N_next).
///
/// Build frame at CA_cur from (N_cur→CA_cur) and peptide plane normal
/// (derived from C_prev, N_cur, CA_prev). Then C_cur local position
/// comes from angle (π - ANGLE_N_CA_C) around N-CA axis + ψ rotation.
fn place_c_cur(
  n_cur : Vector3,
  ca_cur : Vector3,
  c_prev : Vector3,
  ca_prev : Vector3,
  psi_deg : Double,
) -> Vector3 {
  // Peptide plane normal from previous peptide bond
  let v_caprev_cprev = c_prev.sub(ca_prev)
  let v_cprev_ncur = n_cur.sub(c_prev)
  let n_peptide = v_caprev_cprev
    .normalize()
    .cross(v_cprev_ncur.normalize())
    .normalize()

  // Frame at CA_cur: e1 = N_cur→CA_cur, e2 from peptide plane, e3 = e1×e2
  let v_nca = ca_cur.sub(n_cur)
  let e1 = v_nca.normalize()
  let e2_part = n_peptide.cross(e1)
  let e2 = if e2_part.norm() > 1.0e-10 {
    e2_part.normalize()
  } else {
    Vector3::unit_y()
  }
  let e3 = e1.cross(e2).normalize()

  // C_cur in local frame: angle = (π - ANGLE_N_CA_C) around N-CA,
  // then rotated around N-CA by ψ
  let lcl = Vector3::new(
    @math.cos(@math.PI - ANGLE_N_CA_C),
    @math.sin(@math.PI - ANGLE_N_CA_C),
    0.0,
  )
  let rot = RotationMatrix3::rotation_axis_angle(e1, psi_deg * DEG2RAD)
  let rotated = rot.transform(lcl)

  to_global(rotated.scale(BOND_CA_C), ca_cur, e1, e2, e3)
}

///|
/// Place carbonyl O(i) in the peptide plane C(i)-N(i+1).
/// The peptide plane is defined by (C(i), N(i+1), CA(i+1)) — or if that
/// isn't available, by (C(i), N(i+1), CA(i)) — O sits ~120° from
/// N(i+1)-C(i) direction, on the opposite side from N(i+1).
fn place_carbonyl(c_i : Vector3, n_next : Vector3, ca_cur : Vector3) -> Vector3 {
  let v_c_nn = n_next.sub(c_i)
  let v_c_ca = ca_cur.sub(c_i)
  let e1 = v_c_nn.normalize()
  let n_plane = e1.cross(v_c_ca.normalize()).normalize()
  let e2 = n_plane.cross(e1).normalize()
  // O direction: angle t = π - ANGLE_C_C_O from C→N_next, in opposite direction
  let t = @math.PI - ANGLE_C_C_O
  let lcl = Vector3::new(-@math.cos(t), @math.sin(t), 0.0)
  to_global(lcl.scale(BOND_C_O), c_i, e1, e2, n_plane)
}

///|
/// Place carbonyl O for the LAST residue (no N(i+1) available).
/// O sits in the N-CA-C plane on the "outside" of the triangle.
fn place_carbonyl_terminal(
  n_i : Vector3,
  ca_i : Vector3,
  c_i : Vector3,
) -> Vector3 {
  let v_c_n = n_i.sub(c_i)
  let v_c_ca = ca_i.sub(c_i)
  let e1 = v_c_n.normalize()
  let n_plane = e1.cross(v_c_ca.normalize()).normalize()
  let e2 = n_plane.cross(e1).normalize()
  let t = @math.PI - ANGLE_C_C_O
  let lcl = Vector3::new(-@math.cos(t), @math.sin(t), 0.0)
  to_global(lcl.scale(BOND_C_O), c_i, e1, e2, n_plane)
}

// ============================================================================
// Peptide bond N placement given ω
// ============================================================================

///|
/// Place N(i) given C(i-1), CA(i-1), and ω (CA(i-1)-C(i-1)-N(i)-CA(i) ≈ 180°).
///
/// Strategy: start from initial guess N_guess = C + (BOND_C_N * e1)
/// where e1 = normalize(C_prev→CA_prev). Then rotate around C_prev→CA_prev
/// by (ω - 180°) to enforce the peptide bond dihedral.
fn place_n_cur(
  c_prev : Vector3,
  ca_prev : Vector3,
  omega_deg : Double,
) -> Vector3 {
  // Initial guess: N along C→CA direction (ω=180° case)
  let axis = ca_prev.sub(c_prev).normalize()
  let n_guess = c_prev.add(axis.scale(BOND_C_N))
  // Rotate N_guess around this axis by (ω - 180°)
  let offset = (omega_deg - 180.0) * DEG2RAD
  let rot = RotationMatrix3::rotation_axis_angle(axis, offset)
  let lcl = n_guess.sub(c_prev)
  let rotated = rot.transform(lcl)
  c_prev.add(rotated)
}

// ============================================================================
// Public API — backbone builder
// ============================================================================

///|
/// Reconstruct a single-chain backbone AtomArray from per-residue
/// φ/ψ/ω angles (degrees). Peptide bonds default to ω = 180° (trans)
/// if not specified. Residues without φ/ψ use 0°.
///
/// Atom order per residue: [N, CA, C, O].
pub fn BackboneAngleChain::build_backbone(
  self : BackboneAngleChain,
  res_names~ : Array[String] = [],
  chain_id~ : Char = 'A',
  start_res_id~ : Int32 = 1,
) -> AtomArray {
  let n_res = self.angles.length()
  if n_res == 0 {
    return AtomArray::new()
  }
  let rn = if res_names.length() == n_res {
    res_names
  } else {
    let out : Array[String] = Array(capacity=n_res)
    let mut i = 0
    while i < n_res {
      out.push("GLY")
      i = i + 1
    }
    out
  }

  let coords : Array[Vector3] = Array(capacity=n_res * 4)
  let names : Array[String] = Array(capacity=n_res * 4)
  let res_names_out : Array[String] = Array(capacity=n_res * 4)
  let chain_ids_arr : Array[Char] = Array(capacity=n_res * 4)
  let res_ids : Array[Int32] = Array(capacity=n_res * 4)
  let elements : Array[String] = Array(capacity=n_res * 4)
  let b_factors : Array[Double] = Array(capacity=n_res * 4)
  let occupancies : Array[Double] = Array(capacity=n_res * 4)
  let altlocs : Array[Char] = Array(capacity=n_res * 4)

  // First residue in canonical frame
  let (n0, ca0, c0) = initial_backbone_triad()

  // We'll build residue i by first placing its N,CA,C, then O comes
  // after we know N of i+1
  // Approach: push (N,CA,C) immediately, defer O until next iteration

  // --- Residue 0: N, CA, C (O deferred until we know N1) ---
  coords.push(n0)
  names.push("N")
  res_names_out.push(rn[0])
  chain_ids_arr.push(chain_id)
  res_ids.push(start_res_id)
  elements.push("N")
  b_factors.push(0.0)
  occupancies.push(1.0)
  altlocs.push(' ')
  coords.push(ca0)
  names.push("CA")
  res_names_out.push(rn[0])
  chain_ids_arr.push(chain_id)
  res_ids.push(start_res_id)
  elements.push("C")
  b_factors.push(0.0)
  occupancies.push(1.0)
  altlocs.push(' ')
  coords.push(c0)
  names.push("C")
  res_names_out.push(rn[0])
  chain_ids_arr.push(chain_id)
  res_ids.push(start_res_id)
  elements.push("C")
  b_factors.push(0.0)
  occupancies.push(1.0)
  altlocs.push(' ')

  if n_res == 1 {
    let o0 = place_carbonyl_terminal(n0, ca0, c0)
    coords.push(o0)
    names.push("O")
    res_names_out.push(rn[0])
    chain_ids_arr.push(chain_id)
    res_ids.push(start_res_id)
    elements.push("O")
    b_factors.push(0.0)
    occupancies.push(1.0)
    altlocs.push(' ')
    return AtomArray::{
      coords,
      names,
      res_names: res_names_out,
      chain_ids: chain_ids_arr,
      res_ids,
      elements,
      b_factors,
      occupancies,
      altlocs,
    }
  }

  // --- Residues 1..n_res-1 ---
  let mut c_prev = c0
  let mut ca_prev = ca0
  let mut res_id_cur = start_res_id + 1
  let mut i = 1
  while i < n_res {
    // Peptide bond ω(i-1) from previous residue's omega field
    let omega_prev = match self.angles[i - 1].omega {
      Some(v) => v
      None => 180.0
    }
    let phi_i = match self.angles[i].phi {
      Some(v) => v
      None => 0.0
    }
    let psi_i = match self.angles[i].psi {
      Some(v) => v
      None => 0.0
    }

    // 1. Place N(i) using peptide bond geometry
    let n_cur = place_n_cur(c_prev, ca_prev, omega_prev)

    // 2. Place CA(i) — use φ(i) as the rotation of N(i)→CA(i) in peptide plane
    let ca_cur = place_ca_cur(c_prev, n_cur, ca_prev, phi_i)

    // 3. Place C(i) — use ψ(i) relative to N(i)-CA(i)
    let c_cur = place_c_cur(n_cur, ca_cur, c_prev, ca_prev, psi_i)

    // 4. NOW place carbonyl of residue i-1 (we know N(i))
    let o_prev = place_carbonyl(c_prev, n_cur, ca_prev)
    coords.push(o_prev)
    names.push("O")
    res_names_out.push(rn[i - 1])
    chain_ids_arr.push(chain_id)
    res_ids.push(start_res_id + i - 1)
    elements.push("O")
    b_factors.push(0.0)
    occupancies.push(1.0)
    altlocs.push(' ')

    // 5. Push current residue N, CA, C (O deferred)
    coords.push(n_cur)
    names.push("N")
    res_names_out.push(rn[i])
    chain_ids_arr.push(chain_id)
    res_ids.push(res_id_cur)
    elements.push("N")
    b_factors.push(0.0)
    occupancies.push(1.0)
    altlocs.push(' ')
    coords.push(ca_cur)
    names.push("CA")
    res_names_out.push(rn[i])
    chain_ids_arr.push(chain_id)
    res_ids.push(res_id_cur)
    elements.push("C")
    b_factors.push(0.0)
    occupancies.push(1.0)
    altlocs.push(' ')
    coords.push(c_cur)
    names.push("C")
    res_names_out.push(rn[i])
    chain_ids_arr.push(chain_id)
    res_ids.push(res_id_cur)
    elements.push("C")
    b_factors.push(0.0)
    occupancies.push(1.0)
    altlocs.push(' ')

    // Advance state
    c_prev = c_cur
    ca_prev = ca_cur
    res_id_cur = res_id_cur + 1
    i = i + 1
  }

  // --- Last residue O ---
  let last_idx = n_res - 1
  let o_last = place_carbonyl_terminal(
    coords[last_idx * 4],
    coords[last_idx * 4 + 1],
    coords[last_idx * 4 + 2],
  )
  coords.push(o_last)
  names.push("O")
  res_names_out.push(rn[last_idx])
  chain_ids_arr.push(chain_id)
  res_ids.push(start_res_id + last_idx)
  elements.push("O")
  b_factors.push(0.0)
  occupancies.push(1.0)
  altlocs.push(' ')

  AtomArray::{
    coords,
    names,
    res_names: res_names_out,
    chain_ids: chain_ids_arr,
    res_ids,
    elements,
    b_factors,
    occupancies,
    altlocs,
  }
}

// ============================================================================
// Utility — build a default trans-peptide chain from φ/ψ arrays
// ============================================================================

///|
/// Build a BackboneAngleChain from parallel φ and ψ arrays (degrees).
/// All peptide bonds default to trans (ω = 180°).
pub fn BackboneAngleChain::from_phi_psi(
  phi : Array[Double?],
  psi : Array[Double?],
) -> BackboneAngleChain {
  let n = phi.length()
  let out : Array[BackboneAngles] = Array(capacity=n)
  let mut i = 0
  while i < n {
    let omega : Double? = if i < n - 1 { Some(180.0) } else { None }
    out.push(BackboneAngles::{ phi: phi[i], psi: psi[i], omega, })
    i = i + 1
  }
  BackboneAngleChain::{ angles: out, }
}

///|
pub fn BackboneAngleChain::uniform(
  n_res : Int32,
  phi : Double,
  psi : Double,
) -> BackboneAngleChain {
  let out : Array[BackboneAngles] = Array(capacity=n_res)
  let mut i = 0
  while i < n_res {
    let p : Double? = if i == 0 { None } else { Some(phi) }
    let ps : Double? = if i == n_res - 1 { None } else { Some(psi) }
    let om : Double? = if i < n_res - 1 { Some(180.0) } else { None }
    out.push(BackboneAngles::{ phi: p, psi: ps, omega: om, })
    i = i + 1
  }
  BackboneAngleChain::{ angles: out, }
}

// ============================================================================
// Sidechain standard geometry
// ============================================================================

///|
pub const BOND_CA_CB : Double = 1.540

///|
pub const BOND_CB_CG : Double = 1.520

///|
pub const BOND_CG_CD : Double = 1.520

///|
pub const BOND_CD_CE : Double = 1.520

///|
pub const BOND_N_CA_CB_DEG : Double = 110.0

///|
pub const BOND_CB_CA_C_DEG : Double = 110.0

///|
pub const BOND_N_CA_CB_RAD : Double = 1.9199

///|
pub const BOND_CB_CA_C_RAD : Double = 1.9199

///|
pub const TETRAHEDRAL_DEG : Double = 109.5

///|
pub const TETRAHEDRAL_RAD : Double = 1.9199

// ============================================================================
// Sidechain data structures
// ============================================================================

///|
/// Side-chain χ¹..χ⁴ dihedral angles for one residue, degrees.
/// Up to 4 angles; unused slots are None.
pub struct SidechainAngles {
  chi1 : Double?
  chi2 : Double?
  chi3 : Double?
  chi4 : Double?
}

///|
/// Full internal-coordinate description for one residue:
/// backbone (φ/ψ/ω) + sidechain (χ¹..χ⁴).
pub struct ResidueAngles {
  backbone : BackboneAngles
  sidechain : SidechainAngles
}

///|
/// Full chain: backbone + sidechain angles for every residue.
pub struct ChainAngleChain {
  residues : Array[ResidueAngles]
}

// ============================================================================
// Sidechain forward: extract χ angles from an existing AtomArray
// ============================================================================

///|
/// Extract all available χ dihedral angles for residue `ri`.
pub fn AtomArray::get_sidechain_angles(
  self : AtomArray,
  ri : Int32,
) -> SidechainAngles {
  let chi1 = self.dihedral_side_chain(ri, 1)
  let chi2 = self.dihedral_side_chain(ri, 2)
  let chi3 = self.dihedral_side_chain(ri, 3)
  let chi4 = self.dihedral_side_chain(ri, 4)
  SidechainAngles::{ chi1, chi2, chi3, chi4, }
}

///|
/// Extract full chain-level ResidueAngles (backbone + sidechain).
pub fn AtomArray::get_all_angles(self : AtomArray) -> Array[ResidueAngles] {
  let res_starts = self.get_residue_starts()
  let n = res_starts.length()
  let out : Array[ResidueAngles] = Array(capacity=n)
  let mut i = 0
  while i < n {
    let phi = self.backbone_phi(i)
    let psi = self.backbone_psi(i)
    let omega : Double? = if i < n - 1 { self.backbone_omega(i) } else { None }
    let bb = BackboneAngles::{ phi, psi, omega, }
    let sc = self.get_sidechain_angles(i)
    out.push(ResidueAngles::{ backbone: bb, sidechain: sc, })
    i = i + 1
  }
  out
}

// ============================================================================
// Sidechain inverse geometry helpers
// ============================================================================

///|
/// Unwrap an optional Double with default fallback.
/// Avoids verbose match syntax inside builders.
fn unwrap_double(opt : Double?, default : Double) -> Double {
  match opt {
    Some(v) => v
    None => default
  }
}

///|
/// Generic sidechain atom placement.
///
///   pivot     = 旋转轴起点（上上级，已放好）
///   connect   = 旋转轴终点（直接连接目标的那个原子，已放好）
///   ref_atom  = 参考原子，用于定义 χ=0° 时的初始平面
///   bond_len  = connect → target 距离 Å
///   bond_angle = pivot-connect-target 角度 rad（默认四面体 TETRAHEDRAL_RAD）
///   chi_deg   = 绕 pivot→connect 轴旋转角度 degrees
///
/// 做法：在 connect 处建局部框架：
///   e1 = normalize(connect - pivot)
///   e2 = 在 (ref_atom - connect) × e1 平面内垂直 e1
///   e3 = e1 × e2
/// 目标初方向在 (cos(π-angle), sin(π-angle), 0) 局部坐标，
/// 再绕 e1 旋转 chi_deg。
fn place_side(
  pivot : Vector3,
  connect : Vector3,
  ref_atom : Vector3,
  bond_len : Double,
  bond_angle : Double,
  chi_deg : Double,
) -> Vector3 {
  let v_connect_pivot = connect.sub(pivot)
  let e1 = v_connect_pivot.normalize() // 旋转轴 pivot→connect
  let v_ref_connect = ref_atom.sub(connect)
  let e2_part = v_ref_connect.normalize().cross(e1)
  let e2 = if e2_part.norm() > 1.0e-10 {
    e2_part.normalize()
  } else {
    Vector3::unit_y()
  }
  let e3 = e1.cross(e2).normalize()

  // 目标初方向：在局部坐标中离 connect 距离 bond_len，
  // 与 e1 夹角 = bond_angle
  let t = @math.PI - bond_angle
  let lcl = Vector3::new(@math.cos(t), @math.sin(t), 0.0)
  // 绕 e1 旋转 chi
  let rot = RotationMatrix3::rotation_axis_angle(e1, chi_deg * DEG2RAD)
  let rotated = rot.transform(lcl.scale(bond_len))

  Vector3::new(
    connect.x + rotated.x * e1.x + rotated.y * e2.x + rotated.z * e3.x,
    connect.y + rotated.x * e1.y + rotated.y * e2.y + rotated.z * e3.y,
    connect.z + rotated.x * e1.z + rotated.y * e2.z + rotated.z * e3.z,
  )
}

///|
/// Place CB given N, CA, C (backbone). CB sits in the N-CA-C plane,
/// on the opposite side of CA from C, with bond angle N-CA-CB = 110°.
/// This works for all non-GLY residues.
fn place_cb(n : Vector3, ca : Vector3, c : Vector3) -> Vector3 {
  let v_can = n.sub(ca).normalize()
  let v_cac = c.sub(ca).normalize()
  // 平分 N-CA 和 C-CA 的方向，取反向 = CB 方向
  let mid = v_can.add(v_cac).normalize()
  let cb_dir = mid.scale(-1.0) // 反方向 = N-CA-C 的另一边
  // 但我们还需要精确的 N-CA-CB 键角 110°，以及 C-CA-CB 110°
  // 在平面内旋转 mid 方向以得到精确 110°
  // 简化：直接用 bond angle 在平面内旋转
  let n_plane = v_can.cross(v_cac).normalize() // N-CA-C 平面法向
  // 期望 CB 在平面内，使 N-CA-CB = C-CA-CB = 110°
  // 实际上两者互斥（加起来等于 N-CA-C），所以取 N-CA-CB = 110°
  // 角度 N-CA-C：
  let angle_nca_c = v_can.angle(v_cac)
  // 平分角 + (angle_nca_c - 110°)/2 旋转
  let offset = angle_nca_c - BOND_N_CA_CB_RAD
  let rot = RotationMatrix3::rotation_axis_angle(n_plane, offset / 2.0)
  let cb_dir_fine = rot.transform(v_can).normalize()
  ca.add(cb_dir_fine.scale(BOND_CA_CB))
}

// ============================================================================
// Sidechain inverse — per-AA atom-adding steps
// ============================================================================

///|
/// Internal accumulator: tracks placed backbone + sidechain atoms
/// as (name, res_idx, coord) tuples so we can chain geometric references.
struct AtomPool {
  names : Array[String]
  res_idxs : Array[Int32]
  coords : Array[Vector3]
}

///|
fn new_atom_pool() -> AtomPool {
  AtomPool::{ names: [], res_idxs: [], coords: [], }
}

///|
fn AtomPool::add(
  self : AtomPool,
  name : String,
  res_idx : Int32,
  coord : Vector3,
) -> Unit {
  self.names.push(name)
  self.res_idxs.push(res_idx)
  self.coords.push(coord)
}

///|
/// Find the coordinate of atom `name` within a specific residue.
fn AtomPool::get(self : AtomPool, res_idx : Int32, name : String) -> Vector3? {
  let mut i = 0
  while i < self.names.length() {
    if self.res_idxs[i] == res_idx && self.names[i] == name {
      return Some(self.coords[i])
    }
    i = i + 1
  }
  None
}

///|
/// Unwrap version — returns zero vector if missing (caller should ensure
/// atom is present; this simplifies MoonBit Option handling).
fn AtomPool::get_u(self : AtomPool, res_idx : Int32, name : String) -> Vector3 {
  let mut i = 0
  while i < self.names.length() {
    if self.res_idxs[i] == res_idx && self.names[i] == name {
      return self.coords[i]
    }
    i = i + 1
  }
  Vector3::new(0.0, 0.0, 0.0)
}

// ============================================================================
// Main sidechain builder — dispatch per residue type
// ============================================================================

///|
/// Full reconstruction: backbone + sidechain → AtomArray.
/// Takes an array of (res_name, ResidueAngles) per residue.
pub fn ChainAngleChain::build_full(
  self : ChainAngleChain,
  res_names : Array[String],
  chain_id~ : Char = 'A',
  start_res_id~ : Int32 = 1,
) -> AtomArray {
  let n_res = self.residues.length()
  if n_res == 0 {
    return AtomArray::new()
  }

  // 1. 用 BackboneAngleChain 构建 backbone-only AtomArray
  let bb_arr : Array[BackboneAngles] = Array(capacity=n_res)
  let mut i = 0
  while i < n_res {
    bb_arr.push(self.residues[i].backbone)
    i = i + 1
  }
  let bb_chain = BackboneAngleChain::{ angles: bb_arr, }
  let bb_only = bb_chain.build_backbone(res_names~, chain_id~, start_res_id~)

  // 2. 把 backbone 原子装进 AtomPool（N/CA/C/O 各 res）
  let pool = new_atom_pool()
  let mut resi = 0
  while resi < n_res {
    let base = resi * 4 // N, CA, C, O
    pool.add("N", resi, bb_only.get_coord(base))
    pool.add("CA", resi, bb_only.get_coord(base + 1))
    pool.add("C", resi, bb_only.get_coord(base + 2))
    pool.add("O", resi, bb_only.get_coord(base + 3))
    resi = resi + 1
  }

  // 3. 为每个 residue 放 sidechain 原子（向 pool 追加）
  let mut ri = 0
  while ri < n_res {
    let rn = res_names[ri]
    let sc = self.residues[ri].sidechain
    build_residue_sidechain(pool, ri, rn, sc)
    ri = ri + 1
  }

  // 4. 把 pool 输出成完整 AtomArray
  return pool_to_atomarray(pool, res_names, chain_id, start_res_id)
}

// ============================================================================
// Sidechain per-residue dispatch
// ============================================================================

///|
fn build_residue_sidechain(
  pool : AtomPool,
  ri : Int32,
  res_name : String,
  sc : SidechainAngles,
) -> Unit {
  match res_name {
    "ALA" => place_ala_sidechain(pool, ri)
    "GLY" => ()
    "VAL" => place_val_sidechain(pool, ri, sc)
    "LEU" => place_leu_sidechain(pool, ri, sc)
    "ILE" => place_ile_sidechain(pool, ri, sc)
    "SER" => place_ser_sidechain(pool, ri, sc)
    "THR" => place_thr_sidechain(pool, ri, sc)
    "CYS" => place_cys_sidechain(pool, ri, sc)
    "MET" => place_met_sidechain(pool, ri, sc)
    "ASP" => place_asp_sidechain(pool, ri, sc)
    "GLU" => place_glu_sidechain(pool, ri, sc)
    "ASN" => place_asn_sidechain(pool, ri, sc)
    "GLN" => place_gln_sidechain(pool, ri, sc)
    "LYS" => place_lys_sidechain(pool, ri, sc)
    "ARG" => place_arg_sidechain(pool, ri, sc)
    "PRO" => place_pro_sidechain(pool, ri, sc)
    "PHE" => place_phe_sidechain(pool, ri, sc)
    "TYR" => place_tyr_sidechain(pool, ri, sc)
    "TRP" => place_trp_sidechain(pool, ri, sc)
    "HIS" => place_his_sidechain(pool, ri, sc)
    _ => ()
  }
}

// ============================================================================
// Individual AA sidechain builders
// ============================================================================

///|
/// ALA: just CB (one hydrogen replaced by methyl — we model CB only)
fn place_ala_sidechain(pool : AtomPool, ri : Int32) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
}

///|
/// VAL: CB + CG1 + CG2 (χ¹ 绕 CA-CB)
fn place_val_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let cg1 = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1)
  pool.add("CG1", ri, cg1)
  // CG2: 绕 CA-CB 轴再 +120°
  let cg2 = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1 + 120.0)
  pool.add("CG2", ri, cg2)
}

///|
/// LEU: CB + CG + CD1 + CD2 (χ¹ + χ²)
fn place_leu_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let cg = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1)
  pool.add("CG", ri, cg)
  let chi2 = unwrap_double(sc.chi2, 180.0)
  let cd1 = place_side(cb, cg, ca, BOND_CG_CD, TETRAHEDRAL_RAD, chi2)
  pool.add("CD1", ri, cd1)
  let cd2 = place_side(cb, cg, ca, BOND_CG_CD, TETRAHEDRAL_RAD, chi2 + 120.0)
  pool.add("CD2", ri, cd2)
}

///|
/// ILE: CB + CG1 + CG2 + CD1 (特殊：χ¹ 绕 CA-CB 指向 CG1)
fn place_ile_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let cg1 = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1)
  pool.add("CG1", ri, cg1)
  let cg2 = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1 - 120.0)
  pool.add("CG2", ri, cg2)
  let chi2 = unwrap_double(sc.chi2, 180.0)
  let cd1 = place_side(cb, cg1, ca, BOND_CG_CD, TETRAHEDRAL_RAD, chi2)
  pool.add("CD1", ri, cd1)
}

///|
/// SER: CB + OG (χ¹)
fn place_ser_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let og = place_side(ca, cb, n, 1.41, BOND_N_CA_CB_RAD, chi1)
  pool.add("OG", ri, og)
}

///|
/// THR: CB + OG1 + CG2 (特殊)
fn place_thr_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let og1 = place_side(ca, cb, n, 1.43, BOND_N_CA_CB_RAD, chi1)
  pool.add("OG1", ri, og1)
  let cg2 = place_side(ca, cb, n, 1.51, TETRAHEDRAL_RAD, chi1 + 120.0)
  pool.add("CG2", ri, cg2)
}

///|
/// CYS: CB + SG (χ¹)
fn place_cys_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let sg = place_side(ca, cb, n, 1.81, TETRAHEDRAL_RAD, chi1)
  pool.add("SG", ri, sg)
}

///|
/// MET: CB + CG + SD + CE (χ¹ + χ² + χ³)
fn place_met_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let cg = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1)
  pool.add("CG", ri, cg)
  let chi2 = unwrap_double(sc.chi2, 180.0)
  let sd = place_side(cb, cg, ca, 1.79, TETRAHEDRAL_RAD, chi2)
  pool.add("SD", ri, sd)
  let chi3 = unwrap_double(sc.chi3, 180.0)
  let ce = place_side(cg, sd, cb, 1.78, TETRAHEDRAL_RAD, chi3)
  pool.add("CE", ri, ce)
}

///|
/// ASP: CB + CG + OD1 + OD2 (χ¹ + χ², 羧基 OD1/OD2 约 120°)
fn place_asp_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let cg = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1)
  pool.add("CG", ri, cg)
  let chi2 = unwrap_double(sc.chi2, 180.0)
  let od1 = place_side(cb, cg, ca, 1.25, 2.0944, chi2) // 120°
  pool.add("OD1", ri, od1)
  let od2 = place_side(cb, cg, ca, 1.25, 2.0944, chi2 + 180.0) // 对侧
  pool.add("OD2", ri, od2)
}

///|
/// GLU: CB + CG + CD + OE1 + OE2 (χ¹ + χ² + χ³)
fn place_glu_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let cg = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1)
  pool.add("CG", ri, cg)
  let chi2 = unwrap_double(sc.chi2, 180.0)
  let cd = place_side(cb, cg, ca, BOND_CG_CD, TETRAHEDRAL_RAD, chi2)
  pool.add("CD", ri, cd)
  let chi3 = unwrap_double(sc.chi3, 180.0)
  let oe1 = place_side(cg, cd, cb, 1.25, 2.0944, chi3)
  pool.add("OE1", ri, oe1)
  let oe2 = place_side(cg, cd, cb, 1.25, 2.0944, chi3 + 180.0)
  pool.add("OE2", ri, oe2)
}

///|
/// ASN: CB + CG + OD1 + ND2 (χ¹ + χ²)
fn place_asn_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let cg = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1)
  pool.add("CG", ri, cg)
  let chi2 = unwrap_double(sc.chi2, 180.0)
  let od1 = place_side(cb, cg, ca, 1.23, 2.0944, chi2)
  pool.add("OD1", ri, od1)
  let nd2 = place_side(cb, cg, ca, 1.33, 2.0944, chi2 + 180.0)
  pool.add("ND2", ri, nd2)
}

///|
/// GLN: CB + CG + CD + OE1 + NE2 (χ¹ + χ² + χ³)
fn place_gln_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let cg = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1)
  pool.add("CG", ri, cg)
  let chi2 = unwrap_double(sc.chi2, 180.0)
  let cd = place_side(cb, cg, ca, BOND_CG_CD, TETRAHEDRAL_RAD, chi2)
  pool.add("CD", ri, cd)
  let chi3 = unwrap_double(sc.chi3, 180.0)
  let oe1 = place_side(cg, cd, cb, 1.23, 2.0944, chi3)
  pool.add("OE1", ri, oe1)
  let ne2 = place_side(cg, cd, cb, 1.33, 2.0944, chi3 + 180.0)
  pool.add("NE2", ri, ne2)
}

///|
/// LYS: CB + CG + CD + CE + NZ (χ¹ + χ² + χ³ + χ⁴)
fn place_lys_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let cg = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1)
  pool.add("CG", ri, cg)
  let chi2 = unwrap_double(sc.chi2, 180.0)
  let cd = place_side(cb, cg, ca, BOND_CG_CD, TETRAHEDRAL_RAD, chi2)
  pool.add("CD", ri, cd)
  let chi3 = unwrap_double(sc.chi3, 180.0)
  let ce = place_side(cg, cd, cb, BOND_CD_CE, TETRAHEDRAL_RAD, chi3)
  pool.add("CE", ri, ce)
  let chi4 = unwrap_double(sc.chi4, 180.0)
  let nz = place_side(cd, ce, cg, 1.48, TETRAHEDRAL_RAD, chi4)
  pool.add("NZ", ri, nz)
}

///|
/// ARG: CB + CG + CD + NE + CZ + NH1 + NH2 (χ¹..χ⁴)
fn place_arg_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let cg = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1)
  pool.add("CG", ri, cg)
  let chi2 = unwrap_double(sc.chi2, 180.0)
  let cd = place_side(cb, cg, ca, BOND_CG_CD, TETRAHEDRAL_RAD, chi2)
  pool.add("CD", ri, cd)
  let chi3 = unwrap_double(sc.chi3, 180.0)
  let ne = place_side(cg, cd, cb, 1.46, TETRAHEDRAL_RAD, chi3)
  pool.add("NE", ri, ne)
  let chi4 = unwrap_double(sc.chi4, 180.0)
  let cz = place_side(cd, ne, cg, 1.33, 2.0944, chi4)
  pool.add("CZ", ri, cz)
  let nh1 = place_side(ne, cz, cd, 1.33, 2.0944, chi4 - 120.0)
  pool.add("NH1", ri, nh1)
  let nh2 = place_side(ne, cz, cd, 1.33, 2.0944, chi4 + 120.0)
  pool.add("NH2", ri, nh2)
}

///|
/// PRO: CB + CG + CD (5 元环闭合, χ¹)
fn place_pro_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let cg = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1)
  pool.add("CG", ri, cg)
  // CD: 环闭合到 N，用 bond length ≈ 1.47，键角用环约束近似
  let cd = place_side(cb, cg, ca, 1.50, 1.8326, chi1) // 105°
  pool.add("CD", ri, cd)
}

///|
/// PHE: CB + CG + CD1 + CD2 + CE1 + CE2 + CZ (芳香环, χ¹+χ²)
fn place_phe_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let cg = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1)
  pool.add("CG", ri, cg)
  // 芳香环 6 元环：χ² = 0°（planar）时 CD1 在 CG 的一侧
  let chi2 = unwrap_double(sc.chi2, 90.0)
  let cd1 = place_side(cb, cg, ca, 1.39, 2.0944, chi2)
  pool.add("CD1", ri, cd1)
  let cd2 = place_side(cb, cg, ca, 1.39, 2.0944, chi2 + 180.0)
  pool.add("CD2", ri, cd2)
  // CE1/CE2/CZ 用简单延伸近似（芳香环精确闭合工作量大）
  let ce1 = place_side(cg, cd1, cb, 1.39, 2.0944, chi2)
  pool.add("CE1", ri, ce1)
  let ce2 = place_side(cg, cd2, cb, 1.39, 2.0944, chi2 + 180.0)
  pool.add("CE2", ri, ce2)
  let cz = place_side(cd1, ce1, cg, 1.39, 2.0944, chi2 + 180.0)
  pool.add("CZ", ri, cz)
}

///|
/// TYR: PHE + OH
fn place_tyr_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  place_phe_sidechain(pool, ri, sc)
  let cz = pool.get_u(ri, "CZ")
  let ce1 = pool.get_u(ri, "CE1")
  let ce2 = pool.get_u(ri, "CE2")
  let oh = place_side(ce1, cz, ce2, 1.38, 2.0944, 0.0)
  pool.add("OH", ri, oh)
}

///|
/// HIS: CB + CG + ND1 + CD2 + CE1 (咪唑环)
fn place_his_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let cg = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1)
  pool.add("CG", ri, cg)
  let chi2 = unwrap_double(sc.chi2, 90.0)
  let nd1 = place_side(cb, cg, ca, 1.37, 2.0944, chi2)
  pool.add("ND1", ri, nd1)
  let cd2 = place_side(cb, cg, ca, 1.35, 2.0944, chi2 + 180.0)
  pool.add("CD2", ri, cd2)
  let ce1 = place_side(cg, nd1, cb, 1.32, 2.0944, chi2 + 180.0)
  pool.add("CE1", ri, ce1)
}

///|
/// TRP: 吲哚环 — 近似（最复杂芳香环，先做框架）
fn place_trp_sidechain(
  pool : AtomPool,
  ri : Int32,
  sc : SidechainAngles,
) -> Unit {
  let n = pool.get_u(ri, "N")
  let ca = pool.get_u(ri, "CA")
  let c = pool.get_u(ri, "C")
  let cb = place_cb(n, ca, c)
  pool.add("CB", ri, cb)
  let chi1 = unwrap_double(sc.chi1, 180.0)
  let cg = place_side(ca, cb, n, BOND_CB_CG, TETRAHEDRAL_RAD, chi1)
  pool.add("CG", ri, cg)
  let chi2 = unwrap_double(sc.chi2, 90.0)
  let cd1 = place_side(cb, cg, ca, 1.37, 2.0944, chi2)
  pool.add("CD1", ri, cd1)
  let cd2 = place_side(cb, cg, ca, 1.43, 2.0944, chi2 + 180.0)
  pool.add("CD2", ri, cd2)
  let ne1 = place_side(cg, cd1, cb, 1.37, 2.0944, chi2 + 180.0)
  pool.add("NE1", ri, ne1)
  let ce2 = place_side(cd1, ne1, cg, 1.39, 2.0944, chi2)
  pool.add("CE2", ri, ce2)
  let ce3 = place_side(cg, cd2, cb, 1.39, 2.0944, chi2)
  pool.add("CE3", ri, ce3)
  let cz2 = place_side(ne1, ce2, cd1, 1.39, 2.0944, chi2)
  pool.add("CZ2", ri, cz2)
  let cz3 = place_side(cd2, ce3, cg, 1.39, 2.0944, chi2 + 180.0)
  pool.add("CZ3", ri, cz3)
  let ch2 = place_side(ce2, cz2, ne1, 1.38, 2.0944, chi2 + 180.0)
  pool.add("CH2", ri, ch2)
}

// ============================================================================
// AtomPool → AtomArray output
// ============================================================================

///|
fn pool_to_atomarray(
  pool : AtomPool,
  res_names : Array[String],
  chain_id : Char,
  start_res_id : Int32,
) -> AtomArray {
  let n = pool.names.length()
  let coords : Array[Vector3] = Array(capacity=n)
  let names : Array[String] = Array(capacity=n)
  let res_names_out : Array[String] = Array(capacity=n)
  let chain_ids : Array[Char] = Array(capacity=n)
  let res_ids : Array[Int32] = Array(capacity=n)
  let elements : Array[String] = Array(capacity=n)
  let b_factors : Array[Double] = Array(capacity=n)
  let occupancies : Array[Double] = Array(capacity=n)
  let altlocs : Array[Char] = Array(capacity=n)
  let mut i = 0
  while i < n {
    coords.push(pool.coords[i])
    names.push(pool.names[i])
    let ri = pool.res_idxs[i]
    let rn = if ri < res_names.length() { res_names[ri] } else { "UNK" }
    res_names_out.push(rn)
    chain_ids.push(chain_id)
    res_ids.push(start_res_id + ri)
    elements.push(infer_element_simple(pool.names[i]))
    b_factors.push(0.0)
    occupancies.push(1.0)
    altlocs.push(' ')
    i = i + 1
  }
  AtomArray::{
    coords,
    names,
    res_names: res_names_out,
    chain_ids,
    res_ids,
    elements,
    b_factors,
    occupancies,
    altlocs,
  }
}

///|
/// Very simple PDB atom name → element for our synthesized atoms.
fn infer_element_simple(atom_name : String) -> String {
  let first = atom_name[0]
  match first {
    'C' => "C"
    'N' => "N"
    'O' => "O"
    'S' => "S"
    _ => "C" // default
  }
}

// ============================================================================
// Backward-compatible: ChainAngleChain helpers
// ============================================================================

///|
/// Build a ChainAngleChain from arrays of (backbone angles, sidechain angles).
pub fn ChainAngleChain::from_arrays(
  bb : Array[BackboneAngles],
  sc : Array[SidechainAngles],
) -> ChainAngleChain {
  let n = bb.length()
  let out : Array[ResidueAngles] = Array(capacity=n)
  let mut i = 0
  while i < n {
    out.push(ResidueAngles::{ backbone: bb[i], sidechain: sc[i], })
    i = i + 1
  }
  ChainAngleChain::{ residues: out, }
}

///|
/// Extract full forward angles from an AtomArray → ChainAngleChain.
pub fn AtomArray::get_chain_angle_chain(self : AtomArray) -> ChainAngleChain {
  let arr = self.get_all_angles()
  ChainAngleChain::{ residues: arr, }
}
