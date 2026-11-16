/// System/Topology object for MD simulations.
///
/// OpenMM equivalent: `System` class
///
/// Unified container that holds particles (with mass), forces,
/// constraints, and periodic box information.

///|
/// A particle in the system
pub struct OmmParticle {
  mass : Double
  charge : Double
  sigma : Double // LJ sigma
  epsilon : Double // LJ epsilon
}

///|
pub fn OmmParticle::new(
  mass~ : Double = 1.0,
  charge~ : Double = 0.0,
  sigma~ : Double = 0.0,
  epsilon~ : Double = 0.0,
) -> OmmParticle {
  { mass, charge, sigma, epsilon, }
}

///|
/// Bond constraint specification
pub struct OmmBondConstraint {
  particle1 : Int32
  particle2 : Int32
  distance : Double
}

///|
/// The MD System
pub struct OmmSystem {
  particles : Array[OmmParticle]
  constraints : Array[OmmBondConstraint]
  mut box : Array[Double] // [lx, ly, lz], zeros = non-periodic
  bonds : Array[(Int32, Int32, Double, Double)] // (i, j, k, r0)
  angles : Array[(Int32, Int32, Int32, Double, Double)]
  dihedrals : Array[(Int32, Int32, Int32, Int32, Double, Double, Int32)]
}

///|
pub fn OmmSystem::new() -> OmmSystem {
  {
    particles: [],
    constraints: [],
    box: [0.0, 0.0, 0.0],
    bonds: [],
    angles: [],
    dihedrals: [],
  }
}

///|
pub fn OmmSystem::num_particles(self : OmmSystem) -> Int32 {
  self.particles.length()
}

///|
pub fn OmmSystem::num_constraints(self : OmmSystem) -> Int32 {
  self.constraints.length()
}

///|
pub fn OmmSystem::num_bonds(self : OmmSystem) -> Int32 {
  self.bonds.length()
}

///|
pub fn OmmSystem::num_angles(self : OmmSystem) -> Int32 {
  self.angles.length()
}

///|
pub fn OmmSystem::num_dihedrals(self : OmmSystem) -> Int32 {
  self.dihedrals.length()
}

///|
pub fn OmmSystem::add_particle(
  self : OmmSystem,
  particle~ : OmmParticle,
) -> Int32 {
  let idx = self.particles.length()
  self.particles.push(particle)
  idx
}

///|
pub fn OmmSystem::add_constraint(
  self : OmmSystem,
  p1 : Int32,
  p2 : Int32,
  distance : Double,
) -> Unit {
  self.constraints.push({ particle1: p1, particle2: p2, distance, })
}

///|
pub fn OmmSystem::add_bond(
  self : OmmSystem,
  i : Int32,
  j : Int32,
  k : Double,
  r0 : Double,
) -> Unit {
  self.bonds.push((i, j, k, r0))
}

///|
pub fn OmmSystem::add_angle(
  self : OmmSystem,
  i : Int32,
  j : Int32,
  k : Int32,
  k_theta : Double,
  theta0 : Double,
) -> Unit {
  self.angles.push((i, j, k, k_theta, theta0))
}

///|
pub fn OmmSystem::add_dihedral(
  self : OmmSystem,
  i : Int32,
  j : Int32,
  k : Int32,
  l : Int32,
  k_phi : Double,
  phi0 : Double,
  n : Int32,
) -> Unit {
  self.dihedrals.push((i, j, k, l, k_phi, phi0, n))
}

///|
pub fn OmmSystem::set_periodic_box(
  self : OmmSystem,
  lx : Double,
  ly : Double,
  lz : Double,
) -> Unit {
  self.box = [lx, ly, lz]
}

///|
pub fn OmmSystem::is_periodic(self : OmmSystem) -> Bool {
  self.box[0] > 0.0 && self.box[1] > 0.0 && self.box[2] > 0.0
}

///|
pub fn OmmSystem::get_masses(self : OmmSystem) -> Array[Double] {
  let n = self.particles.length()
  let masses : Array[Double] = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    masses[i] = self.particles[i].mass
  }
  masses
}

///|
pub fn OmmSystem::get_charges(self : OmmSystem) -> Array[Double] {
  let n = self.particles.length()
  let charges : Array[Double] = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    charges[i] = self.particles[i].charge
  }
  charges
}

///|
/// Compute total potential energy for the system at given positions.
pub fn OmmSystem::potential_energy(
  self : OmmSystem,
  positions : Array[Array[Double]],
) -> Double {
  let mut energy = 0.0
  // Bond energy
  for b in self.bonds {
    let dx = positions[b.1][0] - positions[b.0][0]
    let dy = positions[b.1][1] - positions[b.0][1]
    let dz = positions[b.1][2] - positions[b.0][2]
    let r = (dx * dx + dy * dy + dz * dz).sqrt()
    energy = energy + 0.5 * b.2 * (r - b.3) * (r - b.3)
  }
  // Angle energy
  for a in self.angles {
    let v1x = positions[a.0][0] - positions[a.1][0]
    let v1y = positions[a.0][1] - positions[a.1][1]
    let v1z = positions[a.0][2] - positions[a.1][2]
    let v2x = positions[a.2][0] - positions[a.1][0]
    let v2y = positions[a.2][1] - positions[a.1][1]
    let v2z = positions[a.2][2] - positions[a.1][2]
    let v1m = (v1x * v1x + v1y * v1y + v1z * v1z).sqrt()
    let v2m = (v2x * v2x + v2y * v2y + v2z * v2z).sqrt()
    if v1m > 0.0 && v2m > 0.0 {
      let ct = (v1x * v2x + v1y * v2y + v1z * v2z) / (v1m * v2m)
      let c = if ct > 1.0 { 1.0 } else if ct < -1.0 { -1.0 } else { ct }
      let theta = @math.acos(c)
      energy = energy + 0.5 * a.3 * (theta - a.4) * (theta - a.4)
    }
  }
  // Dihedral energy: periodic torsion E = k (1 + cos(n*phi - phi0)).
  // phi is the i-j-k-l dihedral in radians; phi0 and k use the same units.
  for di in self.dihedrals {
    let (i0, j0, k0, l0, k_phi, phi0, n) = di
    let b1x = positions[j0][0] - positions[i0][0]
    let b1y = positions[j0][1] - positions[i0][1]
    let b1z = positions[j0][2] - positions[i0][2]
    let b2x = positions[k0][0] - positions[j0][0]
    let b2y = positions[k0][1] - positions[j0][1]
    let b2z = positions[k0][2] - positions[j0][2]
    let b3x = positions[l0][0] - positions[k0][0]
    let b3y = positions[l0][1] - positions[k0][1]
    let b3z = positions[l0][2] - positions[k0][2]
    let m1 = (
      b1y * b2z - b1z * b2y,
      b1z * b2x - b1x * b2z,
      b1x * b2y - b1y * b2x,
    )
    let m2 = (
      b2y * b3z - b2z * b3y,
      b2z * b3x - b2x * b3z,
      b2x * b3y - b2y * b3x,
    )
    let m1_len = @math.pow(m1.0 * m1.0 + m1.1 * m1.1 + m1.2 * m1.2, 0.5)
    let m2_len = @math.pow(m2.0 * m2.0 + m2.1 * m2.1 + m2.2 * m2.2, 0.5)
    let b2_len = @math.pow(b2x * b2x + b2y * b2y + b2z * b2z, 0.5)
    if m1_len > 1.0e-12 && m2_len > 1.0e-12 && b2_len > 1.0e-12 {
      let cos_phi = (m1.0 * m2.0 + m1.1 * m2.1 + m1.2 * m2.2) /
        (m1_len * m2_len)
      let c = if cos_phi > 1.0 {
        1.0
      } else if cos_phi < -1.0 {
        -1.0
      } else {
        cos_phi
      }
      // Sign from the mixed product (m1 x b2) . m2 to place phi in (-pi, pi].
      let cross_x = m1.1 * b2z - m1.2 * b2y
      let cross_y = m1.2 * b2x - m1.0 * b2z
      let cross_z = m1.0 * b2y - m1.1 * b2x
      let sgn = if cross_x * m2.0 + cross_y * m2.1 + cross_z * m2.2 >= 0.0 {
        1.0
      } else {
        -1.0
      }
      let phi = sgn * @math.acos(c)
      energy = energy + k_phi * (1.0 + @math.cos(n.to_double() * phi - phi0))
    }
  }
  // LJ energy
  let n = self.particles.length()
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let sigma = (self.particles[i].sigma + self.particles[j].sigma) * 0.5
      let epsilon = (self.particles[i].epsilon * self.particles[j].epsilon).sqrt()
      if epsilon == 0.0 {
        continue
      }
      let dx = positions[j][0] - positions[i][0]
      let dy = positions[j][1] - positions[i][1]
      let dz = positions[j][2] - positions[i][2]
      let r_sq = dx * dx + dy * dy + dz * dz
      if r_sq == 0.0 {
        continue
      }
      let sr6 = sigma *
        sigma /
        r_sq *
        (sigma * sigma / r_sq) *
        (sigma * sigma / r_sq)
      energy = energy + 4.0 * epsilon * (sr6 * sr6 - sr6)
    }
  }
  energy
}

///|
/// Compute forces (negative gradient) via finite difference.
pub fn OmmSystem::forces(
  self : OmmSystem,
  positions : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = positions.length()
  let forces : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  for i = 0; i < n; i = i + 1 {
    forces[i] = [0.0, 0.0, 0.0]
  }
  let delta = 0.0001
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      let pos_plus : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
      let pos_minus : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
      for j = 0; j < n; j = j + 1 {
        pos_plus[j] = [positions[j][0], positions[j][1], positions[j][2]]
        pos_minus[j] = [positions[j][0], positions[j][1], positions[j][2]]
      }
      pos_plus[i][d] = pos_plus[i][d] + delta
      pos_minus[i][d] = pos_minus[i][d] - delta
      let e_plus = self.potential_energy(pos_plus)
      let e_minus = self.potential_energy(pos_minus)
      forces[i][d] = -(e_plus - e_minus) / (2.0 * delta)
    }
  }
  forces
}

///|
/// Lennard-Jones + Coulomb 成对能量（非周期边界，不做 PME）。
pub fn OmmSystem::potential_energy_nonbonded(
  self : OmmSystem,
  positions : Array[Array[Double]],
) -> (Double, Double) {
  let n = self.particles.length()
  let mut lj = 0.0
  let mut coul = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let sigma = (self.particles[i].sigma + self.particles[j].sigma) * 0.5
      let epsilon = (self.particles[i].epsilon * self.particles[j].epsilon).sqrt()
      let qi = self.particles[i].charge
      let qj = self.particles[j].charge
      if epsilon == 0.0 && qi * qj == 0.0 {
        continue
      }
      let dx = positions[j][0] - positions[i][0]
      let dy = positions[j][1] - positions[i][1]
      let dz = positions[j][2] - positions[i][2]
      let r2 = dx * dx + dy * dy + dz * dz
      let r = if r2 > 1.0e-12 { @math.pow(r2, 0.5) } else { 1.0e-6 }
      if epsilon > 0.0 && sigma > 0.0 {
        let sr6 = (sigma / r) * (sigma / r) * (sigma / r) * (sigma / r) * (sigma / r) * (sigma / r)
        lj = lj + 4.0 * epsilon * (sr6 * sr6 - sr6)
      }
      let qq = qi * qj
      if qq != 0.0 {
        coul = coul + 138.935485 * qq / r
      }
    }
  }
  (lj, coul)
}

///|
/// 能量分解：返回 (bond, angle, dihedral, lj, coul) 各分量（kJ/mol）。
pub fn OmmSystem::energy_decomposed(
  self : OmmSystem,
  positions : Array[Array[Double]],
) -> (Double, Double, Double, Double, Double) {
  let mut bond = 0.0
  for b in self.bonds {
    let dx = positions[b.1][0] - positions[b.0][0]
    let dy = positions[b.1][1] - positions[b.0][1]
    let dz = positions[b.1][2] - positions[b.0][2]
    let r = (dx * dx + dy * dy + dz * dz).sqrt()
    bond = bond + 0.5 * b.2 * (r - b.3) * (r - b.3)
  }
  let mut ang = 0.0
  for a in self.angles {
    let v1x = positions[a.0][0] - positions[a.1][0]
    let v1y = positions[a.0][1] - positions[a.1][1]
    let v1z = positions[a.0][2] - positions[a.1][2]
    let v2x = positions[a.2][0] - positions[a.1][0]
    let v2y = positions[a.2][1] - positions[a.1][1]
    let v2z = positions[a.2][2] - positions[a.1][2]
    let v1m = (v1x * v1x + v1y * v1y + v1z * v1z).sqrt()
    let v2m = (v2x * v2x + v2y * v2y + v2z * v2z).sqrt()
    if v1m > 0.0 && v2m > 0.0 {
      let ct = (v1x * v2x + v1y * v2y + v1z * v2z) / (v1m * v2m)
      let c = if ct > 1.0 { 1.0 } else if ct < -1.0 { -1.0 } else { ct }
      let theta = @math.acos(c)
      ang = ang + 0.5 * a.3 * (theta - a.4) * (theta - a.4)
    }
  }
  let mut dih = 0.0
  for di in self.dihedrals {
    let (i0, j0, k0, l0, k_phi, phi0, n) = di
    let b1x = positions[j0][0] - positions[i0][0]
    let b1y = positions[j0][1] - positions[i0][1]
    let b1z = positions[j0][2] - positions[i0][2]
    let b2x = positions[k0][0] - positions[j0][0]
    let b2y = positions[k0][1] - positions[j0][1]
    let b2z = positions[k0][2] - positions[j0][2]
    let b3x = positions[l0][0] - positions[k0][0]
    let b3y = positions[l0][1] - positions[k0][1]
    let b3z = positions[l0][2] - positions[k0][2]
    let m1 = (b1y * b2z - b1z * b2y, b1z * b2x - b1x * b2z, b1x * b2y - b1y * b2x)
    let m2 = (b2y * b3z - b2z * b3y, b2z * b3x - b2x * b3z, b2x * b3y - b2y * b3x)
    let m1_len = @math.pow(m1.0 * m1.0 + m1.1 * m1.1 + m1.2 * m1.2, 0.5)
    let m2_len = @math.pow(m2.0 * m2.0 + m2.1 * m2.1 + m2.2 * m2.2, 0.5)
    let b2_len = @math.pow(b2x * b2x + b2y * b2y + b2z * b2z, 0.5)
    if m1_len > 1.0e-12 && m2_len > 1.0e-12 && b2_len > 1.0e-12 {
      let cos_phi = (m1.0 * m2.0 + m1.1 * m2.1 + m1.2 * m2.2) / (m1_len * m2_len)
      let c = if cos_phi > 1.0 { 1.0 } else if cos_phi < -1.0 { -1.0 } else { cos_phi }
      let cross_x = m1.1 * b2z - m1.2 * b2y
      let cross_y = m1.2 * b2x - m1.0 * b2z
      let cross_z = m1.0 * b2y - m1.1 * b2x
      let sgn = if cross_x * m2.0 + cross_y * m2.1 + cross_z * m2.2 >= 0.0 { 1.0 } else { -1.0 }
      let phi = sgn * @math.acos(c)
      dih = dih + k_phi * (1.0 + @math.cos(n.to_double() * phi - phi0))
    }
  }
  let (lj, coul) = self.potential_energy_nonbonded(positions)
  (bond, ang, dih, lj, coul)
}

// ============================================================================
// A9: 约束违反率监控
// ============================================================================

pub struct OmmConstraintViolation {
  constraint_idx : Int32
  particle1 : Int32
  particle2 : Int32
  target_distance : Double
  actual_distance : Double
  absolute_violation : Double
  relative_violation : Double  // |actual - target| / target
}

///|
/// 检查所有约束的违反情况。
/// 返回 (max_violation, mean_violation, violations_list)。
pub fn OmmSystem::check_constraints(
  self : OmmSystem,
  positions : Array[Array[Double]],
) -> (Double, Double, Array[OmmConstraintViolation]) {
  let n = self.constraints.length()
  let out : Array[OmmConstraintViolation] = Array(capacity=n)
  let mut max_v = 0.0
  let mut sum_v = 0.0
  let mut i = 0
  while i < n {
    let c = self.constraints[i]
    let p1 = positions[c.particle1]
    let p2 = positions[c.particle2]
    let dx = p2[0] - p1[0]
    let dy = p2[1] - p1[1]
    let dz = p2[2] - p1[2]
    let r = (dx * dx + dy * dy + dz * dz).sqrt()
    let abs_v = (r - c.distance).abs()
    let rel_v = if c.distance > 1.0e-12 { abs_v / c.distance } else { 0.0 }
    if abs_v > max_v { max_v = abs_v }
    sum_v = sum_v + abs_v
    out.push({
      constraint_idx: i,
      particle1: c.particle1,
      particle2: c.particle2,
      target_distance: c.distance,
      actual_distance: r,
      absolute_violation: abs_v,
      relative_violation: rel_v,
    })
    i = i + 1
  }
  let mean_v = if n > 0 { sum_v / n.to_double() } else { 0.0 }
  (max_v, mean_v, out)
}

// ============================================================================
// B2: 轻量级 PDB → OmmSystem 参数化（按元素/残基查表）
// ============================================================================
//
// 不依赖外部 XML 力场文件。
// 从 PDB ATOM 记录推断 element/mass/sigma/epsilon/charge，
// 并根据残基名识别标准键合拓扑（水、单原子离子、简单分子）。
//
// 对于**完整蛋白质+配体**模拟，仍需 AMBER prmtop 或 OpenMM XML 力场。
// 本函数适合：水盒子、离子溶液、简单溶剂化体系、粗粒度快速测试。

fn omm_element_params(element : String) -> (Double, Double, Double, Double) {
  match element {
    "H"  => (1.008,  0.0,    0.0,    0.0)
    "C"  => (12.011, 3.4,    0.258,  0.0)
    "N"  => (14.007, 3.25,   0.163,  0.0)
    "O"  => (15.999, 3.166,  0.155,  0.0)
    "F"  => (18.998, 2.9,    0.08,   0.0)
    "P"  => (30.974, 3.742,  0.377,  0.0)
    "S"  => (32.06,  3.55,   0.35,   0.0)
    "Cl" => (35.453, 3.5,    0.53,   -1.0)
    "Na" => (22.99,  3.33,   0.0115, 1.0)
    "K"  => (39.1,   3.74,   0.0008, 1.0)
    "Ca" => (40.08,  3.362,  0.55,   2.0)
    "Mg" => (24.31,  2.86,   0.87,   2.0)
    "I"  => (126.9,  3.9,    0.40,   -1.0)
    _    => (1.0, 0.0, 0.0, 0.0)
  }
}

///|
/// 从 PDB 文本直接构建 OmmSystem。
/// 自动从原子名推断 element/mass/LJ 参数，从残基名补氢键拓扑。
///
/// 注意：这是轻量级参数化，**不适合真实蛋白质模拟**。
/// 真正的力场参数化请用 AMBER prmtop 或 OpenMM XML 力场。
pub fn omm_system_from_pdb_fast(pdb_text : String) -> OmmSystem {
  let system = OmmSystem::new()
  let lines_sv = pdb_text.split("\n").to_array()
  let atoms_info : Array[(String, String, Double, Double, Double, Double)] = Array(capacity=128)
  let mut i : Int32 = 0
  while i < lines_sv.length() {
    let ln = lines_sv[i].to_owned()
    if ln.length() >= 54 {
      let rec = ln.substring(start=0, end=4)
      if rec == "ATOM" || ln.substring(start=0, end=6) == "HETATM" {
        let name = ln.substring(start=12, end=16).trim().to_owned()
        let resname = ln.substring(start=17, end=20).trim().to_owned()
        let mut elem_upper = "C"
        let nc = name.length()
        let mut j = 0
        while j < nc {
          let ci = name.unsafe_get(j)
          if (ci >= 65 && ci <= 90) || (ci >= 97 && ci <= 122) {
            let upp_i = if ci >= 97 { ci - 32 } else { ci }
            elem_upper = String::make(1, upp_i.unsafe_to_char())
            break
          }
          j = j + 1
        }
        let is_water = resname == "HOH" || resname == "WAT" || resname == "TIP3"
        // omm_element_params 返回 (mass, sigma, epsilon, charge)
        let (mass, sig, eps, chg) = match elem_upper {
          "H" => if is_water { (1.008, 0.0, 0.0, 0.417) } else { omm_element_params("H") }
          "O" => if is_water { (15.999, 3.166, 0.155, -0.834) } else { omm_element_params("O") }
          _   => omm_element_params(elem_upper)
        }
        atoms_info.push((elem_upper, resname, mass, chg, sig, eps))
      }
    }
    i = i + 1
  }
  let n = atoms_info.length()
  i = 0
  while i < n {
    let (_, _, mass, chg, sig, eps) = atoms_info[i]
    ignore(system.add_particle(particle=OmmParticle::new(mass=mass, charge=chg, sigma=sig, epsilon=eps)))
    i = i + 1
  }
  i = 0
  while i + 2 < n {
    let r1 = atoms_info[i].1
    if r1 == "HOH" || r1 == "WAT" || r1 == "TIP3" {
      ignore(system.add_constraint(i, i + 1, 0.9572))
      ignore(system.add_constraint(i, i + 2, 0.9572))
      ignore(system.add_bond(i, i + 1, 450.0, 0.9572))
      ignore(system.add_bond(i, i + 2, 450.0, 0.9572))
      ignore(system.add_angle(i + 1, i, i + 2, 55.0, 104.52 * 3.141592653589793 / 180.0))
    }
    i = i + 3
  }
  system
}

// ============================================================================
// B2 完整版：PDB + OpenMM XML 力场 → OmmSystem 一键参数化
// ============================================================================

///|
/// 从 PDB ATOM 名称推断元素。处理双字母元素：Na/Cl/Fe/Ca/Mg/Zn/Br。
/// 返回标准格式：首字母大写 + 次字母小写（如 "Cl"，不是 "CL"）。
/// PDB 原子名可能全大写（如 "NA"、"CL"），统一处理。
fn omm_pdb_element_from_atom_name(name : String) -> String {
  let nc : Int32 = name.length()
  let mut j : Int32 = 0
  // 找第一个字母（ASCII）
  while j < nc {
    let ci : UInt16 = name.unsafe_get(j)
    if (ci >= 65 && ci <= 90) || (ci >= 97 && ci <= 122) {
      let upp_ci = if ci >= 97 { ci - 32 } else { ci }   // 首字母大写
      let fc : Char = upp_ci.unsafe_to_char()
      let upp_result = String::make(1, fc)
      // 看有没有第二个字母（PDB 可能是小写或大写）
      if j + 1 < nc {
        let c2 : UInt16 = name.unsafe_get(j + 1)
        let is_letter = (c2 >= 65 && c2 <= 90) || (c2 >= 97 && c2 <= 122)
        if is_letter {
          // 次字母统一为小写（标准元素符号：Na, Cl, Fe...）
          let low_c2 : UInt16 = if c2 >= 97 { c2 } else { c2 + 32 }
          let sc : Char = low_c2.unsafe_to_char()
          let two = upp_result + String::make(1, sc)
          if two == "Na" || two == "Cl" || two == "Fe" || two == "Ca" || two == "Mg" || two == "Zn" || two == "Br" {
            return two
          }
        }
      }
      return upp_result
    }
    j = j + 1
  }
  "C"
}

///|
/// **B2 一键入口**：PDB 文本 + OpenMM XML 力场 → 参数化 OmmSystem。
///
/// - `pdb_text`: 标准 PDB 格式（ATOM / HETATM 记录）
/// - `ff_content`: OpenMM XML 力场文件内容（如 amber14-all.xml）
///
/// 工作流：PDB 解析 → 从残基模板匹配 atom type → 无匹配则从元素回退 → 组装 bond/angle/torsion → 非键参数赋值。
/// 当 XML 缺少某个键/角/二面角参数时使用合理默认值（bond: 300 kcal/mol/A^2, 1.5 A; angle: 50 kcal/mol/rad^2, 109.5°; torsion: 跳过）。
/// 当 XML 也找不到元素模板时，非键参数回退到 omm_element_params 查表。
pub fn omm_system_from_pdb_ff(
  pdb_text : String,
  ff_content : String,
) -> OmmSystem {
  // 1. 解析 PDB → 提取每个原子的 (name, resname, element)
  let atoms : Array[OmmXmlAtomSpec] = Array(capacity=256)
  let lines_sv = pdb_text.split("\n").to_array()
  let mut i : Int32 = 0
  while i < lines_sv.length() {
    let ln = lines_sv[i].to_owned()
    let rec = if ln.length() >= 6 { ln.substring(start=0, end=6) } else { "" }
    if rec == "ATOM  " || rec == "HETATM" || (ln.length() >= 4 && ln.substring(start=0, end=4) == "ATOM") {
      let name = if ln.length() >= 16 { ln.substring(start=12, end=16).trim().to_owned() } else { "X" }
      let resname = if ln.length() >= 20 { ln.substring(start=17, end=20).trim().to_owned() } else { "UNK" }
      let elem = omm_pdb_element_from_atom_name(name)
      atoms.push(OmmXmlAtomSpec::{ name, resname, element: elem })
    }
    i = i + 1
  }
  // 2. 解析 XML 力场
  let ff = omm_parse_xml_ff(ff_content)
  // 3. 构建 OmmSystem（XML 力场参数化）
  omm_xml_ff_build_system(ff, atoms)
}

///|
/// B2 便捷入口：内置极简 XML 力场（覆盖 H/C/N/O/P/S/Na/Cl 的基本参数），不需要外部文件。
/// 适合快速测试或简单体系。
pub fn omm_system_from_pdb_ff_builtin(pdb_text : String) -> OmmSystem {
  let builtin_ff =
    "<?xml version=\"1.0\"?>\n" +
    "<ForceField name=\"builtin_basic\">\n" +
    "  <Type name=\"H\" element=\"H\" mass=\"1.008\"/>\n" +
    "  <Type name=\"C\" element=\"C\" mass=\"12.011\"/>\n" +
    "  <Type name=\"N\" element=\"N\" mass=\"14.007\"/>\n" +
    "  <Type name=\"O\" element=\"O\" mass=\"15.999\"/>\n" +
    "  <Type name=\"P\" element=\"P\" mass=\"30.974\"/>\n" +
    "  <Type name=\"S\" element=\"S\" mass=\"32.06\"/>\n" +
    "  <Type name=\"Na\" element=\"Na\" mass=\"22.99\"/>\n" +
    "  <Type name=\"Cl\" element=\"Cl\" mass=\"35.453\"/>\n" +
    "  <Residue name=\"HOH\">\n" +
    "    <Atom name=\"O\" type=\"O\"/>\n" +
    "    <Atom name=\"H1\" type=\"H\"/>\n" +
    "    <Atom name=\"H2\" type=\"H\"/>\n" +
    "    <Bond from=\"0\" to=\"1\"/>\n" +
    "    <Bond from=\"0\" to=\"2\"/>\n" +
    "  </Residue>\n" +
    "  <Atom class=\"A\" type=\"H\" charge=\"0.000\" sigma=\"0.00000\" epsilon=\"0.000\"/>\n" +
    "  <Atom class=\"A\" type=\"C\" charge=\"0.000\" sigma=\"0.33997\" epsilon=\"0.35982\"/>\n" +
    "  <Atom class=\"A\" type=\"N\" charge=\"0.000\" sigma=\"0.32499\" epsilon=\"0.19246\"/>\n" +
    "  <Atom class=\"A\" type=\"O\" charge=\"0.000\" sigma=\"0.31656\" epsilon=\"0.15500\"/>\n" +
    "  <Atom class=\"A\" type=\"P\" charge=\"0.000\" sigma=\"0.37417\" epsilon=\"0.46900\"/>\n" +
    "  <Atom class=\"A\" type=\"S\" charge=\"0.000\" sigma=\"0.35500\" epsilon=\"0.25000\"/>\n" +
    "  <Atom class=\"A\" type=\"Na\" charge=\"1.000\" sigma=\"0.33284\" epsilon=\"0.00310\"/>\n" +
    "  <Atom class=\"A\" type=\"Cl\" charge=\"-1.000\" sigma=\"0.44010\" epsilon=\"0.10000\"/>\n" +
    "</ForceField>\n"
  omm_system_from_pdb_ff(pdb_text, builtin_ff)
}
