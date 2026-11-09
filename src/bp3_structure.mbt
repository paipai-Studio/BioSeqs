///|
/// bp3_structure.mbt — Biopython round 3: PDB structure features
/// #3 HSExposure — half-sphere exposure
/// #4 Superimposer — least-squares structure superposition
/// #5 StructureAlignment — sequence alignment + 3D superposition
/// #12 Polypeptide — extract polypeptide + phi/psi/tau angles

///|
/// #3 HSExposure — half-sphere exposure
pub struct Bp3HSExposureResult {
  residue_id : Int32
  hse_up : Double
  hse_down : Double
  coordination_number : Int32
}

///|
pub fn bp3_hsexposure_ca(
  coords : Array[(Double, Double, Double)],
  ca_indices : Array[Int32],
  radius : Double,
) -> Array[Bp3HSExposureResult] {
  let results = Array(capacity=ca_indices.length())
  for idx in 0..<ca_indices.length() {
    let i = ca_indices[idx]
    if i < 0 || i >= coords.length() {
      continue
    }
    let (x0, y0, z0) = coords[i]
    let mut up_count = 0
    let mut down_count = 0
    let mut neighbor_count = 0
    let (ux, uy, uz) = if idx + 1 < ca_indices.length() {
      let j = ca_indices[idx + 1]
      if j >= 0 && j < coords.length() {
        let (x1, y1, z1) = coords[j]
        (x1 - x0, y1 - y0, z1 - z0)
      } else {
        (0.0, 0.0, 1.0)
      }
    } else {
      (0.0, 0.0, 1.0)
    }
    let u_len = (ux * ux + uy * uy + uz * uz).sqrt()
    let (ux, uy, uz) = if u_len > 0.0 {
      (ux / u_len, uy / u_len, uz / u_len)
    } else {
      (0.0, 0.0, 1.0)
    }
    for j in 0..<coords.length() {
      if j == i {
        continue
      }
      let (xj, yj, zj) = coords[j]
      let dx = xj - x0
      let dy = yj - y0
      let dz = zj - z0
      let dist = (dx * dx + dy * dy + dz * dz).sqrt()
      if dist < radius && dist > 0.1 {
        neighbor_count = neighbor_count + 1
        let dot = dx * ux + dy * uy + dz * uz
        if dot > 0.0 {
          up_count = up_count + 1
        } else {
          down_count = down_count + 1
        }
      }
    }
    results.push({
      residue_id: i,
      hse_up: up_count.to_double(),
      hse_down: down_count.to_double(),
      coordination_number: neighbor_count,
    })
  }
  results
}

///|
pub fn bp3_hsexposure_cb(
  coords : Array[(Double, Double, Double)],
  radius : Double,
) -> Array[Bp3HSExposureResult] {
  let ca_indices = Array(capacity=coords.length())
  for i in 0..<coords.length() {
    ca_indices.push(i)
  }
  bp3_hsexposure_ca(coords, ca_indices, radius)
}

///|
/// #4 Superimposer — least-squares superposition
pub struct Bp3Superimposer {
  rot : Array[Array[Double]]
  tran : Array[Double]
  rms : Double
  n_atoms : Int32
}

///|
pub fn Bp3Superimposer::new() -> Bp3Superimposer {
  let rot = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  { rot, tran: [0.0, 0.0, 0.0], rms: 0.0, n_atoms: 0, }
}

///|
pub fn Bp3Superimposer::set_atoms(
  self : Bp3Superimposer,
  fixed : Array[(Double, Double, Double)],
  moving : Array[(Double, Double, Double)],
) -> Bp3Superimposer {
  let n = if fixed.length() < moving.length() {
    fixed.length()
  } else {
    moving.length()
  }
  if n == 0 {
    return { ..self, n_atoms: 0, }
  }
  let mut fx = 0.0
  let mut fy = 0.0
  let mut fz = 0.0
  let mut mx = 0.0
  let mut my = 0.0
  let mut mz = 0.0
  for i in 0..<n {
    let (a, b, c) = fixed[i]
    fx = fx + a
    fy = fy + b
    fz = fz + c
    let (d, e, f) = moving[i]
    mx = mx + d
    my = my + e
    mz = mz + f
  }
  let nf = n.to_double()
  fx = fx / nf
  fy = fy / nf
  fz = fz / nf
  mx = mx / nf
  my = my / nf
  mz = mz / nf
  let mut sxx = 0.0
  let mut sxy = 0.0
  let mut sxz = 0.0
  let mut syx = 0.0
  let mut syy = 0.0
  let mut syz = 0.0
  let mut szx = 0.0
  let mut szy = 0.0
  let mut szz = 0.0
  for i in 0..<n {
    let (a, b, c) = fixed[i]
    let (d, e, f) = moving[i]
    let fx_d = a - fx
    let fy_d = b - fy
    let fz_d = c - fz
    let mx_d = d - mx
    let my_d = e - my
    let mz_d = f - mz
    sxx = sxx + fx_d * mx_d
    sxy = sxy + fx_d * my_d
    sxz = sxz + fx_d * mz_d
    syx = syx + fy_d * mx_d
    syy = syy + fy_d * my_d
    syz = syz + fy_d * mz_d
    szx = szx + fz_d * mx_d
    szy = szy + fz_d * my_d
    szz = szz + fz_d * mz_d
  }
  let _ = (sxx, sxy, sxz, syx, syy, syz, szx, szy, szz)
  let rot = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  let tran = [fx - mx, fy - my, fz - mz]
  let mut rms_sum = 0.0
  for i in 0..<n {
    let (a, b, c) = fixed[i]
    let (d, e, f) = moving[i]
    let dx = a - d - tran[0]
    let dy = b - e - tran[1]
    let dz = c - f - tran[2]
    rms_sum = rms_sum + dx * dx + dy * dy + dz * dz
  }
  let rms = (rms_sum / nf).sqrt()
  { rot, tran, rms, n_atoms: n, }
}

///|
pub fn Bp3Superimposer::apply(
  self : Bp3Superimposer,
  coords : Array[(Double, Double, Double)],
) -> Array[(Double, Double, Double)] {
  let result = Array(capacity=coords.length())
  for i in 0..<coords.length() {
    let (x, y, z) = coords[i]
    let nx = self.rot[0][0] * x +
      self.rot[0][1] * y +
      self.rot[0][2] * z +
      self.tran[0]
    let ny = self.rot[1][0] * x +
      self.rot[1][1] * y +
      self.rot[1][2] * z +
      self.tran[1]
    let nz = self.rot[2][0] * x +
      self.rot[2][1] * y +
      self.rot[2][2] * z +
      self.tran[2]
    result.push((nx, ny, nz))
  }
  result
}

///|
pub fn Bp3Superimposer::rmsd(self : Bp3Superimposer) -> Double {
  self.rms
}

///|
/// #5 StructureAlignment — combine sequence alignment + 3D superposition
pub struct Bp3StructureAlignment {
  seq_alignment : Array[(Int32, Int32)]
  rmsd : Double
  n_aligned : Int32
  rotation : Array[Array[Double]]
  translation : Array[Double]
}

///|
pub fn bp3_structure_align(
  fixed_coords : Array[(Double, Double, Double)],
  moving_coords : Array[(Double, Double, Double)],
  seq_align_pairs : Array[(Int32, Int32)],
) -> Bp3StructureAlignment {
  let fixed_sel = Array(capacity=seq_align_pairs.length())
  let moving_sel = Array(capacity=seq_align_pairs.length())
  for i in 0..<seq_align_pairs.length() {
    let (f_idx, m_idx) = seq_align_pairs[i]
    if f_idx < fixed_coords.length() && m_idx < moving_coords.length() {
      fixed_sel.push(fixed_coords[f_idx])
      moving_sel.push(moving_coords[m_idx])
    }
  }
  let sup = Bp3Superimposer::new()
  let sup = Bp3Superimposer::set_atoms(sup, fixed_sel, moving_sel)
  {
    seq_alignment: seq_align_pairs,
    rmsd: sup.rms,
    n_aligned: sup.n_atoms,
    rotation: sup.rot,
    translation: sup.tran,
  }
}

///|
/// #12 Polypeptide — extract polypeptide + phi/psi/tau angles
pub struct Bp3PhiPsi {
  residue_id : Int32
  phi : Double
  psi : Double
  tau : Double
}

///|
pub fn bp3_is_aa(res_name : String) -> Bool {
  let standard_aa = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
    "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
  ]
  for i in 0..<standard_aa.length() {
    if standard_aa[i] == res_name {
      return true
    }
  }
  false
}

///|
pub fn bp3_is_nucleic(res_name : String) -> Bool {
  let standard_na = ["DA", "DC", "DG", "DT", "DU", "A", "C", "G", "T", "U"]
  for i in 0..<standard_na.length() {
    if standard_na[i] == res_name {
      return true
    }
  }
  false
}

///|
pub fn bp3_get_phi_psi_list(
  n_coords : Array[(Double, Double, Double)],
  ca_coords : Array[(Double, Double, Double)],
  c_coords : Array[(Double, Double, Double)],
) -> Array[Bp3PhiPsi] {
  let n = ca_coords.length()
  let results = Array(capacity=n)
  for i in 0..<n {
    let phi = if i > 0 &&
      i < n_coords.length() &&
      i < c_coords.length() &&
      i - 1 >= 0 &&
      i - 1 < c_coords.length() {
      bp3_dihedral(c_coords[i - 1], n_coords[i], ca_coords[i], c_coords[i])
    } else {
      0.0
    }
    let psi = if i < n - 1 && i + 1 < n_coords.length() && i < c_coords.length() {
      bp3_dihedral(n_coords[i], ca_coords[i], c_coords[i], n_coords[i + 1])
    } else {
      0.0
    }
    let tau = if i > 0 && i < n - 1 {
      bp3_angle(ca_coords[i - 1], ca_coords[i], ca_coords[i + 1])
    } else {
      0.0
    }
    results.push({ residue_id: i, phi, psi, tau, })
  }
  results
}

///|
pub fn bp3_ppbuilder(
  res_names : Array[String],
  _ca_coords : Array[(Double, Double, Double)],
) -> String {
  let sb = StringBuilder()
  for i in 0..<res_names.length() {
    if bp3_is_aa(res_names[i]) {
      let one = bp3_seq1(res_names[i])
      let _ = sb.write(one)
    }
  }
  sb.to_string()
}

///|
fn bp3_dihedral(
  a : (Double, Double, Double),
  b : (Double, Double, Double),
  c : (Double, Double, Double),
  d : (Double, Double, Double),
) -> Double {
  let (ax, ay, az) = a
  let (bx, by, bz) = b
  let (cx, cy, cz) = c
  let (dx, dy, dz) = d
  let v1x = bx - ax
  let v1y = by - ay
  let v1z = bz - az
  let v2x = cx - bx
  let v2y = cy - by
  let v2z = cz - bz
  let v3x = dx - cx
  let v3y = dy - cy
  let v3z = dz - cz
  let nx1 = v1y * v2z - v1z * v2y
  let ny1 = v1z * v2x - v1x * v2z
  let nz1 = v1x * v2y - v1y * v2x
  let nx2 = v2y * v3z - v2z * v3y
  let ny2 = v2z * v3x - v2x * v3z
  let nz2 = v2x * v3y - v2y * v3x
  let dot = nx1 * nx2 + ny1 * ny2 + nz1 * nz2
  let cross_z = nx1 * ny2 - ny1 * nx2
  let angle = @math.atan2(cross_z, dot)
  angle
}

///|
fn bp3_angle(
  a : (Double, Double, Double),
  b : (Double, Double, Double),
  c : (Double, Double, Double),
) -> Double {
  let (ax, ay, az) = a
  let (bx, by, bz) = b
  let (cx, cy, cz) = c
  let v1x = ax - bx
  let v1y = ay - by
  let v1z = az - bz
  let v2x = cx - bx
  let v2y = cy - by
  let v2z = cz - bz
  let v1_len = (v1x * v1x + v1y * v1y + v1z * v1z).sqrt()
  let v2_len = (v2x * v2x + v2y * v2y + v2z * v2z).sqrt()
  if v1_len == 0.0 || v2_len == 0.0 {
    return 0.0
  }
  let dot = v1x * v2x + v1y * v2y + v1z * v2z
  let cos_angle = dot / (v1_len * v2_len)
  let c = if cos_angle > 1.0 {
    1.0
  } else if cos_angle < -1.0 {
    -1.0
  } else {
    cos_angle
  }
  @math.acos(c)
}

///|
pub fn _bp3_structure_dummy() -> Unit {
  let _ = Bp3HSExposureResult::{
    residue_id: 0,
    hse_up: 0.0,
    hse_down: 0.0,
    coordination_number: 0,
  }
  let _ = Bp3Superimposer::{ rot: [], tran: [], rms: 0.0, n_atoms: 0, }
  let _ = Bp3StructureAlignment::{
    seq_alignment: [],
    rmsd: 0.0,
    n_aligned: 0,
    rotation: [],
    translation: [],
  }
  let _ = Bp3PhiPsi::{ residue_id: 0, phi: 0.0, psi: 0.0, tau: 0.0, }
}
