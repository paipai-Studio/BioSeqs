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
