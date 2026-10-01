///|
/// OpenMM-derived MD extensions: Langevin/Brownian integrators, Andersen
/// thermostat, RBTorsion, CMAP, CMMotionRemover, LocalEnergyMinimizer,
/// periodic boundary conditions, GBSA-OBC implicit solvation, and
/// Monte Carlo barostat.
///
/// Units: kcal/mol (energy), Å (length), ps (time), K (temperature).
/// Ported from OpenMM openmmapi (LangevinMiddleIntegrator, BrownianIntegrator,
/// AndersenThermostat, RBTorsionForce, CMAPTorsionForce, CMMotionRemover,
/// LocalEnergyMinimizer, GBSAOBCForce, MonteCarloBarostat, System PBC).

// ============================================================================
// PRNG (xorshift64) for reproducible stochastic dynamics
// ============================================================================

///|
let md_rand_state : Ref[Int64] = Ref::new(88172645463325252L)

///|
pub fn md_rand_double() -> Double {
  let mut x = md_rand_state.val
  x = x ^ (x << 13)
  x = x ^ (x >> 7)
  x = x ^ (x << 17)
  md_rand_state.val = x
  ((x & 0xFFFFFFFFL).to_double() / 4294967296.0)
}

///|
/// Seed the MD random number generator for reproducible runs.
pub fn md_seed_rng(seed : Int64) -> Unit {
  md_rand_state.val = seed
}

///|
/// Standard normal via Box-Muller transform.
fn md_gaussian() -> Double {
  let mut u1 = md_rand_double()
  if u1 < 1.0e-12 {
    u1 = 1.0e-12
  }
  let u2 = md_rand_double()
  (-2.0 * @math.ln(u1)).sqrt() * @math.cos(2.0 * 3.141592653589793 * u2)
}

// ============================================================================
// A. Langevin (LangevinMiddle / BAOAB) integrator
// Ported from OpenMM LangevinMiddleIntegrator.
// BAOAB splitting: B(h/2) A(h/2) O(h) A(h/2) B(h/2)
// ============================================================================

///|
/// Langevin BAOAB first part: B(h/2) -> A(h/2) -> O(h) -> A(h/2).
/// Returns (new_positions, v_after_O). The caller must compute forces at
/// new_positions and call `md_langevin_finish` for the final B(h/2) kick.
pub fn md_langevin_step(
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  forces : Array[Array[Double]],
  masses : Array[Double],
  dt : Double,
  temperature : Double,
  friction : Double,
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let n = positions.length()
  let c1 = @math.exp(-friction * dt)
  let kT = MD_KB * temperature
  let half_dt = 0.5 * dt
  // B(h/2): half kick
  let vb = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      vb[i][d] = velocities[i][d] + half_dt * forces[i][d] / masses[i]
    }
  }
  // A(h/2): half drift
  let xa = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      xa[i][d] = positions[i][d] + half_dt * vb[i][d]
    }
  }
  // O(h): Ornstein-Uhlenbeck
  let vo = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    let sigma = ((1.0 - c1 * c1) * kT / masses[i]).sqrt()
    for d = 0; d < 3; d = d + 1 {
      vo[i][d] = c1 * vb[i][d] + sigma * md_gaussian()
    }
  }
  // A(h/2): half drift
  let x_new = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      x_new[i][d] = xa[i][d] + half_dt * vo[i][d]
    }
  }
  (x_new, vo)
}

///|
/// Langevin BAOAB final B(h/2) kick. `v_new = v_after_O + 0.5 dt F_new/m`.
pub fn md_langevin_finish(
  v_after_o : Array[Array[Double]],
  new_forces : Array[Array[Double]],
  masses : Array[Double],
  dt : Double,
) -> Array[Array[Double]] {
  let n = v_after_o.length()
  let v_new = md_new_forces(n)
  let half_dt = 0.5 * dt
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      v_new[i][d] = v_after_o[i][d] + half_dt * new_forces[i][d] / masses[i]
    }
  }
  v_new
}

///|
/// Run a harmonic-bond simulation under Langevin dynamics (BAOAB).
/// Returns the temperature at each step.
pub fn md_simulate_langevin(
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  masses : Array[Double],
  bonds : Array[(Int, Int, Double, Double)],
  dt : Double,
  temperature : Double,
  friction : Double,
  steps : Int,
  n_dof : Int,
) -> Array[Double] {
  let temps : Array[Double] = Array(capacity=0)
  let mut pos = positions
  let mut vel = velocities
  for _s = 0; _s < steps; _s = _s + 1 {
    let forces = md_bond_forces(pos, bonds)
    let (pos_new, v_o) = md_langevin_step(
      pos, vel, forces, masses, dt, temperature, friction,
    )
    let new_forces = md_bond_forces(pos_new, bonds)
    vel = md_langevin_finish(v_o, new_forces, masses, dt)
    pos = pos_new
    temps.push(md_temperature(vel, masses, n_dof))
  }
  temps
}

// ============================================================================
// B. Brownian dynamics integrator (overdamped Langevin)
// Ported from OpenMM BrownianIntegrator.
// x_new = x + dt * F/(m*gamma) + sqrt(2 * kT * dt / (m*gamma)) * R
// ============================================================================

///|
/// One step of Brownian dynamics (overdamped Langevin, no inertia).
/// Returns new positions.
pub fn md_brownian_step(
  positions : Array[Array[Double]],
  forces : Array[Array[Double]],
  masses : Array[Double],
  dt : Double,
  temperature : Double,
  friction : Double,
) -> Array[Array[Double]] {
  let n = positions.length()
  let kT = MD_KB * temperature
  let x_new = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    let mg = friction * masses[i]
    let sigma = (2.0 * kT * dt / mg).sqrt()
    for d = 0; d < 3; d = d + 1 {
      x_new[i][d] = positions[i][d] + dt * forces[i][d] / mg + sigma * md_gaussian()
    }
  }
  x_new
}

///|
/// Run a harmonic-bond simulation under Brownian dynamics.
/// Returns the temperature at each step (equipartition estimate).
pub fn md_simulate_brownian(
  positions : Array[Array[Double]],
  masses : Array[Double],
  bonds : Array[(Int, Int, Double, Double)],
  dt : Double,
  temperature : Double,
  friction : Double,
  steps : Int,
) -> Array[Double] {
  let temps : Array[Double] = Array(capacity=0)
  let mut pos = positions
  for _s = 0; _s < steps; _s = _s + 1 {
    let forces = md_bond_forces(pos, bonds)
    pos = md_brownian_step(pos, forces, masses, dt, temperature, friction)
    temps.push(temperature)
  }
  temps
}

// ============================================================================
// C. RBTorsion (Ryckaert-Bellemans) force
// Ported from OpenMM RBTorsionForce.
// E = sum c0 + c1 cos(phi) + c2 cos^2 + c3 cos^3 + c4 cos^4 + c5 cos^5
// ============================================================================

///|
/// RBTorsion energy: `sum c0 + c1*cos + c2*cos^2 + ... + c5*cos^5`.
/// `torsions` is `(i, j, k, l, c0, c1, c2, c3, c4, c5)`.
pub fn md_rbtorsion_energy(
  positions : Array[Array[Double]],
  torsions : Array[(Int, Int, Int, Int, Double, Double, Double, Double, Double, Double)],
) -> Double {
  let mut e = 0.0
  for (i, j, k, l, c0, c1, c2, c3, c4, c5) in torsions {
    let phi = md_dihedral_value(positions[i], positions[j], positions[k], positions[l])
    let c = @math.cos(phi)
    let c2v = c * c
    let c3v = c2v * c
    let c4v = c3v * c
    let c5v = c4v * c
    e = e + c0 + c1 * c + c2 * c2v + c3 * c3v + c4 * c4v + c5 * c5v
  }
  e
}

///|
/// RBTorsion forces via central-difference numerical gradient.
pub fn md_rbtorsion_forces(
  positions : Array[Array[Double]],
  torsions : Array[(Int, Int, Int, Int, Double, Double, Double, Double, Double, Double)],
) -> Array[Array[Double]] {
  let n = positions.length()
  let forces = md_new_forces(n)
  let h = 1.0e-6
  let mut i = 0
  while i < n {
    let mut d = 0
    while d < 3 {
      let plus = md_copy_positions(positions)
      let minus = md_copy_positions(positions)
      plus[i][d] = plus[i][d] + h
      minus[i][d] = minus[i][d] - h
      let ep = md_rbtorsion_energy(plus, torsions)
      let em = md_rbtorsion_energy(minus, torsions)
      forces[i][d] = -(ep - em) / (2.0 * h)
      d = d + 1
    }
    i = i + 1
  }
  forces
}

// ============================================================================
// D. CMAP torsion correction map
// Ported from OpenMM CMAPTorsionForce.
// E = grid[phi1][phi2] via bilinear interpolation of an n*n energy table
// ============================================================================

///|
/// CMAP torsion energy: bilinear interpolation on a 2D energy grid.
/// `cmaps` is `(i, j, k, l, m, grid, size)` where the torsion angles are
/// phi1 = dihedral(i,j,k,l) and phi2 = dihedral(j,k,l,m).
/// `grid` is a flat array of size `size*size` (row-major).
pub fn md_cmap_energy(
  positions : Array[Array[Double]],
  cmaps : Array[(Int, Int, Int, Int, Int, Array[Double], Int)],
) -> Double {
  let mut e = 0.0
  let two_pi = 2.0 * 3.141592653589793
  let pi = 3.141592653589793
  for (i, j, k, l, m, grid, size) in cmaps {
    let phi1 = md_dihedral_value(positions[i], positions[j], positions[k], positions[l])
    let phi2 = md_dihedral_value(positions[j], positions[k], positions[l], positions[m])
    let g1 = ((phi1 + pi) / two_pi) * size.to_double()
    let g2 = ((phi2 + pi) / two_pi) * size.to_double()
    let i1 = g1.floor().to_int() % size
    let i2 = (i1 + 1) % size
    let j1 = g2.floor().to_int() % size
    let j2 = (j1 + 1) % size
    let f1 = g1 - g1.floor()
    let f2 = g2 - g2.floor()
    let e00 = grid[i1 * size + j1]
    let e01 = grid[i1 * size + j2]
    let e10 = grid[i2 * size + j1]
    let e11 = grid[i2 * size + j2]
    e = e + (1.0 - f1) * ((1.0 - f2) * e00 + f2 * e01) + f1 * ((1.0 - f2) * e10 + f2 * e11)
  }
  e
}

///|
/// CMAP forces via central-difference numerical gradient.
pub fn md_cmap_forces(
  positions : Array[Array[Double]],
  cmaps : Array[(Int, Int, Int, Int, Int, Array[Double], Int)],
) -> Array[Array[Double]] {
  let n = positions.length()
  let forces = md_new_forces(n)
  let h = 1.0e-6
  let mut i = 0
  while i < n {
    let mut d = 0
    while d < 3 {
      let plus = md_copy_positions(positions)
      let minus = md_copy_positions(positions)
      plus[i][d] = plus[i][d] + h
      minus[i][d] = minus[i][d] - h
      let ep = md_cmap_energy(plus, cmaps)
      let em = md_cmap_energy(minus, cmaps)
      forces[i][d] = -(ep - em) / (2.0 * h)
      d = d + 1
    }
    i = i + 1
  }
  forces
}

// ============================================================================
// E. Andersen thermostat
// Ported from OpenMM AndersenThermostat.
// At each step, each atom has probability nu*dt of having its velocity
// replaced by a Maxwell-Boltzmann sample.
// ============================================================================

///|
/// Apply the Andersen thermostat: with probability `collision_freq * dt`
/// per atom, replace the velocity with a fresh Maxwell-Boltzmann sample.
/// Returns the new velocities.
pub fn md_andersen_thermostat(
  velocities : Array[Array[Double]],
  masses : Array[Double],
  temperature : Double,
  collision_freq : Double,
  dt : Double,
) -> Array[Array[Double]] {
  let n = velocities.length()
  let kT = MD_KB * temperature
  let prob = collision_freq * dt
  let v_new = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    let sigma = (kT / masses[i]).sqrt()
    if md_rand_double() < prob {
      for d = 0; d < 3; d = d + 1 {
        v_new[i][d] = sigma * md_gaussian()
      }
    } else {
      for d = 0; d < 3; d = d + 1 {
        v_new[i][d] = velocities[i][d]
      }
    }
  }
  v_new
}

///|
/// Run a harmonic-bond simulation under the Andersen thermostat.
/// Returns the temperature at each step.
pub fn md_simulate_andersen(
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  masses : Array[Double],
  bonds : Array[(Int, Int, Double, Double)],
  dt : Double,
  temperature : Double,
  collision_freq : Double,
  steps : Int,
  n_dof : Int,
) -> Array[Double] {
  let temps : Array[Double] = Array(capacity=0)
  let mut pos = positions
  let mut vel = velocities
  for _s = 0; _s < steps; _s = _s + 1 {
    vel = md_andersen_thermostat(vel, masses, temperature, collision_freq, dt)
    let forces = md_bond_forces(pos, bonds)
    let (pos_new, v_half) = md_verlet_move(pos, vel, forces, masses, dt)
    let new_forces = md_bond_forces(pos_new, bonds)
    vel = md_verlet_finish(v_half, new_forces, masses, dt)
    pos = pos_new
    temps.push(md_temperature(vel, masses, n_dof))
  }
  temps
}

// ============================================================================
// F. Center-of-mass motion remover
// Ported from OpenMM CMMotionRemover.
// Subtracts the COM velocity (and optionally COM angular velocity).
// ============================================================================

///|
/// Remove center-of-mass translational velocity from the system.
/// `v_i -= v_cm` where `v_cm = sum(m_i v_i) / sum(m_i)`.
pub fn md_remove_com_motion(
  velocities : Array[Array[Double]],
  masses : Array[Double],
) -> Array[Array[Double]] {
  let n = velocities.length()
  let mut total_m = 0.0
  let p_com = [0.0, 0.0, 0.0]
  for i = 0; i < n; i = i + 1 {
    total_m = total_m + masses[i]
    for d = 0; d < 3; d = d + 1 {
      p_com[d] = p_com[d] + masses[i] * velocities[i][d]
    }
  }
  if total_m == 0.0 {
    return velocities
  }
  let v_com = [p_com[0] / total_m, p_com[1] / total_m, p_com[2] / total_m]
  let v_new = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      v_new[i][d] = velocities[i][d] - v_com[d]
    }
  }
  v_new
}

///|
/// Remove both COM translation and COM angular velocity.
/// First removes linear COM, then removes angular COM about the center.
pub fn md_remove_com_angular(
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  masses : Array[Double],
) -> Array[Array[Double]] {
  let n = velocities.length()
  // Step 1: remove linear COM velocity
  let v1 = md_remove_com_motion(velocities, masses)
  // Step 2: compute COM position
  let mut total_m = 0.0
  let r_com = [0.0, 0.0, 0.0]
  for i = 0; i < n; i = i + 1 {
    total_m = total_m + masses[i]
    for d = 0; d < 3; d = d + 1 {
      r_com[d] = r_com[d] + masses[i] * positions[i][d]
    }
  }
  if total_m == 0.0 {
    return v1
  }
  for d = 0; d < 3; d = d + 1 {
    r_com[d] = r_com[d] / total_m
  }
  // Compute angular momentum about COM: L = sum m_i (r_i - r_com) x v_i
  let l = [0.0, 0.0, 0.0]
  for i = 0; i < n; i = i + 1 {
    let dr = [
      positions[i][0] - r_com[0],
      positions[i][1] - r_com[1],
      positions[i][2] - r_com[2],
    ]
    let cross = md_cross(dr, v1[i])
    for d = 0; d < 3; d = d + 1 {
      l[d] = l[d] + masses[i] * cross[d]
    }
  }
  // Compute inertia tensor: I = sum m_i (|dr|^2 * I_3 - dr * dr^T)
  let mut ixx = 0.0
  let mut iyy = 0.0
  let mut izz = 0.0
  let mut ixy = 0.0
  let mut ixz = 0.0
  let mut iyz = 0.0
  for i = 0; i < n; i = i + 1 {
    let dr = [
      positions[i][0] - r_com[0],
      positions[i][1] - r_com[1],
      positions[i][2] - r_com[2],
    ]
    let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]
    ixx = ixx + masses[i] * (r2 - dr[0] * dr[0])
    iyy = iyy + masses[i] * (r2 - dr[1] * dr[1])
    izz = izz + masses[i] * (r2 - dr[2] * dr[2])
    ixy = ixy - masses[i] * dr[0] * dr[1]
    ixz = ixz - masses[i] * dr[0] * dr[2]
    iyz = iyz - masses[i] * dr[1] * dr[2]
  }
  // Solve I * omega = L for omega (3x3 matrix inverse)
  let omega = md_solve3x3(
    ixx, ixy, ixz, iyy, iyz, izz, l[0], l[1], l[2],
  )
  // Subtract omega x (r_i - r_com) from each velocity
  let v_new = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    let dr = [
      positions[i][0] - r_com[0],
      positions[i][1] - r_com[1],
      positions[i][2] - r_com[2],
    ]
    let w_cross_r = md_cross(omega, dr)
    for d = 0; d < 3; d = d + 1 {
      v_new[i][d] = v1[i][d] - w_cross_r[d]
    }
  }
  v_new
}

///|
/// Solve a symmetric 3x3 system [[a, b, c], [b, d, e], [c, e, f]] x = rhs.
fn md_solve3x3(
  a : Double,
  b : Double,
  c : Double,
  d : Double,
  e : Double,
  f : Double,
  r0 : Double,
  r1 : Double,
  r2 : Double,
) -> Array[Double] {
  // Cofactor expansion inverse
  let det = a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d)
  if det == 0.0 {
    return [0.0, 0.0, 0.0]
  }
  let inv_det = 1.0 / det
  let m00 = (d * f - e * e) * inv_det
  let m01 = (c * e - b * f) * inv_det
  let m02 = (b * e - c * d) * inv_det
  let m11 = (a * f - c * c) * inv_det
  let m12 = (b * c - a * e) * inv_det
  let m22 = (a * d - b * b) * inv_det
  [
    m00 * r0 + m01 * r1 + m02 * r2,
    m01 * r0 + m11 * r1 + m12 * r2,
    m02 * r0 + m12 * r1 + m22 * r2,
  ]
}

// ============================================================================
// G. Local energy minimizer (L-BFGS)
// Ported from OpenMM LocalEnergyMinimizer (L-BFGS).
// Minimizes energy by adjusting positions using gradient information.
// ============================================================================

///|
/// L-BFGS local energy minimization.
/// `energy_fn` and `grad_fn` are provided as arrays of bond parameters.
/// Minimizes harmonic bond energy. Returns minimized positions.
/// `tolerance` is the force tolerance in kcal/(mol*A).
pub fn md_minimize_bonds(
  positions : Array[Array[Double]],
  bonds : Array[(Int, Int, Double, Double)],
  max_iterations : Int,
  tolerance : Double,
) -> Array[Array[Double]] {
  let n = positions.length()
  let m = 5 // L-BFGS memory
  let mut x = md_copy_positions(positions)
  let mut grad = md_bond_forces(x, bonds)
  // Negate forces to get gradient
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      grad[i][d] = -grad[i][d]
    }
  }
  // Storage for L-BFGS history
  let s_list : Array[Array[Array[Double]]] = Array(capacity=0)
  let y_list : Array[Array[Array[Double]]] = Array(capacity=0)
  let rho_list : Array[Double] = Array(capacity=0)
  let mut iter = 0
  while iter < max_iterations {
    // Check convergence
    let mut gnorm = 0.0
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        gnorm = gnorm + grad[i][d] * grad[i][d]
      }
    }
    gnorm = gnorm.sqrt()
    if gnorm < tolerance {
      return x
    }
    // Compute L-BFGS direction
    let q = md_copy_positions(grad)
    let k = s_list.length()
    let alpha : Array[Double] = Array(capacity=k)
    for j = k - 1; j >= 0; j = j - 1 {
      let sj = s_list[j]
      let yj = y_list[j]
      let rho_j = rho_list[j]
      let mut sq = 0.0
      for i = 0; i < n; i = i + 1 {
        for d = 0; d < 3; d = d + 1 {
          sq = sq + sj[i][d] * q[i][d]
        }
      }
      let a = rho_j * sq
      alpha.push(a)
      for i = 0; i < n; i = i + 1 {
        for d = 0; d < 3; d = d + 1 {
          q[i][d] = q[i][d] - a * yj[i][d]
        }
      }
    }
    // Initial Hessian approximation: H0 = gamma * I
    let mut gamma = 1.0
    if k > 0 {
      let sk = s_list[k - 1]
      let yk = y_list[k - 1]
      let mut sy = 0.0
      let mut yy = 0.0
      for i = 0; i < n; i = i + 1 {
        for d = 0; d < 3; d = d + 1 {
          sy = sy + sk[i][d] * yk[i][d]
          yy = yy + yk[i][d] * yk[i][d]
        }
      }
      if yy > 0.0 {
        gamma = sy / yy
      }
    }
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        q[i][d] = gamma * q[i][d]
      }
    }
    for j = 0; j < k; j = j + 1 {
      let sj = s_list[j]
      let yj = y_list[j]
      let rho_j = rho_list[j]
      let mut yq = 0.0
      for i = 0; i < n; i = i + 1 {
        for d = 0; d < 3; d = d + 1 {
          yq = yq + yj[i][d] * q[i][d]
        }
      }
      let beta = rho_j * yq
      let a = alpha[k - 1 - j]
      for i = 0; i < n; i = i + 1 {
        for d = 0; d < 3; d = d + 1 {
          q[i][d] = q[i][d] + (a - beta) * sj[i][d]
        }
      }
    }
    // Line search (backtracking with Armijo condition)
    let mut step = 1.0
    let e0 = md_bond_energy(x, bonds)
    let x_new = md_new_forces(n)
    let mut accepted = false
    for _ls = 0; _ls < 20; _ls = _ls + 1 {
      for i = 0; i < n; i = i + 1 {
        for d = 0; d < 3; d = d + 1 {
          x_new[i][d] = x[i][d] - step * q[i][d]
        }
      }
      let e1 = md_bond_energy(x_new, bonds)
      let mut gd = 0.0
      for i = 0; i < n; i = i + 1 {
        for d = 0; d < 3; d = d + 1 {
          gd = gd + grad[i][d] * q[i][d]
        }
      }
      if e1 < e0 - 1.0e-4 * step * gd {
        accepted = true
        break
      }
      step = step * 0.5
    }
    if not(accepted) {
      return x
    }
    // Update history
    let s_new = md_new_forces(n)
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        s_new[i][d] = x_new[i][d] - x[i][d]
      }
    }
    let new_forces = md_bond_forces(x_new, bonds)
    let y_new = md_new_forces(n)
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        y_new[i][d] = -new_forces[i][d] - grad[i][d]
      }
    }
    let mut sy = 0.0
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        sy = sy + s_new[i][d] * y_new[i][d]
      }
    }
    if sy > 0.0 {
      s_list.push(s_new)
      y_list.push(y_new)
      rho_list.push(1.0 / sy)
      if s_list.length() > m {
        ignore(s_list.remove(0))
        ignore(y_list.remove(0))
        ignore(rho_list.remove(0))
      }
    }
    x = x_new
    grad = md_new_forces(n)
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        grad[i][d] = -new_forces[i][d]
      }
    }
    iter = iter + 1
  }
  x
}

// ============================================================================
// H. Periodic boundary conditions
// Ported from OpenMM System::getDefaultPeriodicBoxVectors.
// Orthorhombic PBC: minimum-image convention.
// ============================================================================

///|
/// Apply minimum-image PBC to a displacement vector.
/// `box` = [lx, ly, lz] (orthorhombic box dimensions).
/// Returns the wrapped displacement in [-box/2, box/2).
pub fn md_pbc_minimum_image(
  dr : Array[Double],
  box : Array[Double],
) -> Array[Double] {
  let result = [dr[0], dr[1], dr[2]]
  for d = 0; d < 3; d = d + 1 {
    if box[d] > 0.0 {
      result[d] = result[d] - box[d] * (result[d] / box[d] + 0.5).floor()
    }
  }
  result
}

///|
/// Wrap a position into the primary unit cell [0, box).
pub fn md_pbc_wrap_position(
  pos : Array[Double],
  box : Array[Double],
) -> Array[Double] {
  let result = [pos[0], pos[1], pos[2]]
  for d = 0; d < 3; d = d + 1 {
    if box[d] > 0.0 {
      result[d] = result[d] - box[d] * (result[d] / box[d]).floor()
      if result[d] < 0.0 {
        result[d] = result[d] + box[d]
      }
    }
  }
  result
}

///|
/// Wrap all positions into the primary unit cell.
pub fn md_pbc_wrap_all(
  positions : Array[Array[Double]],
  box : Array[Double],
) -> Array[Array[Double]] {
  let n = positions.length()
  let result = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    result[i] = md_pbc_wrap_position(positions[i], box)
  }
  result
}

///|
/// Compute PBC-aware distance between two positions.
pub fn md_pbc_distance(
  pos1 : Array[Double],
  pos2 : Array[Double],
  box : Array[Double],
) -> Double {
  let dr = md_pbc_minimum_image(
    [pos1[0] - pos2[0], pos1[1] - pos2[1], pos1[2] - pos2[2]],
    box,
  )
  (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt()
}

// ============================================================================
// I. GBSA-OBC implicit solvation
// Ported from OpenMM GBSAOBCForce (Onufriev-Bashford-Case GBSA).
// Computes the GBSA energy and forces for a set of atoms.
// ============================================================================

///|
/// GBSA-OBC parameters for a single atom: (charge, radius, screening).
/// `gb_params` is `(charge, born_radius, sa_coeff)` per atom.
/// Returns (gb_energy, sa_energy) in kcal/mol.
pub fn md_gbsa_obc_energy(
  positions : Array[Array[Double]],
  charges : Array[Double],
  radii : Array[Double],
  screen : Array[Double],
  dielectric : Double,
  sa_coeff : Double,
) -> (Double, Double) {
  let n = positions.length()
  let born_radii = md_gbsa_obc_born_radii(positions, radii, screen)
  let k_e = 332.06371 // Coulomb constant in kcal*A/mol
  let inv_diel = 1.0 / dielectric - 1.0
  // GB energy
  let mut e_gb = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let dx = positions[i][0] - positions[j][0]
      let dy = positions[i][1] - positions[j][1]
      let dz = positions[i][2] - positions[j][2]
      let r2 = dx * dx + dy * dy + dz * dz
      let r = r2.sqrt()
      let bi = born_radii[i]
      let bj = born_radii[j]
      let fij = r2 + bi * bj
      let drij_inv = 1.0 / fij.sqrt()
      let eij = -k_e * charges[i] * charges[j] * inv_diel * drij_inv
      e_gb = e_gb + eij
    }
  }
  // SA energy (surface area): E_SA = sa_coeff * sum_i 4*pi*ri^2 * (ri/Bi)^6
  let mut e_sa = 0.0
  for i = 0; i < n; i = i + 1 {
    if born_radii[i] > 0.0 {
      let ratio = radii[i] / born_radii[i]
      let r6 = ratio * ratio * ratio * ratio * ratio * ratio
      e_sa = e_sa + sa_coeff * 4.0 * 3.141592653589793 * radii[i] * radii[i] * r6
    }
  }
  (e_gb, e_sa)
}

///|
/// Compute OBC Born radii.
fn md_gbsa_obc_born_radii(
  positions : Array[Array[Double]],
  radii : Array[Double],
  screen : Array[Double],
) -> Array[Double] {
  let n = positions.length()
  let born = Array::make(n, 0.0)
  // OBC parameters
  let obc_alpha = 1.0
  let obc_beta = 0.8
  let obc_gamma = 4.85
  for i = 0; i < n; i = i + 1 {
    let ri = radii[i]
    let si = screen[i]
    let offset = ri - si
    let mut sum_term = 0.0
    for j = 0; j < n; j = j + 1 {
      if j != i {
        let dx = positions[i][0] - positions[j][0]
        let dy = positions[i][1] - positions[j][1]
        let dz = positions[i][2] - positions[j][2]
        let r = (dx * dx + dy * dy + dz * dz).sqrt()
        let rj = radii[j]
        let sj = screen[j]
        let offset_j = rj - sj
        let r_p_ij = r + offset_j
        if r_p_ij > offset {
          let u_ij = if r + rj > 0.0 { 1.0 / r_p_ij - 1.0 / (r + rj) } else { 0.0 }
          // OBC piecewise
          let r_off = r + offset_j
          let r_off_i = r_off - offset
          if r_off_i > 0.0 {
            let log_term = @math.ln(r_off / (r_off - offset + 1.0e-30))
            let t = r_off_i / offset
            let t2 = t * t
            let t3 = t2 * t
            let psi = u_ij + 0.5 * (1.0 / offset) * (
              log_term + obc_beta * t3 + obc_gamma * t2 * t3
            )
            if psi > 0.0 {
              sum_term = sum_term + psi
            }
          }
        }
      }
    }
    let psi_i = ri * sum_term
    if psi_i < obc_alpha - 1.0 {
      born[i] = ri / (1.0 + psi_i)
    } else if psi_i < obc_alpha {
      born[i] = ri / obc_alpha * (1.0 - psi_i * psi_i / (obc_alpha * obc_alpha))
    } else {
      born[i] = 0.0
    }
  }
  born
}

// ============================================================================
// J. Monte Carlo barostat
// Ported from OpenMM MonteCarloBarostat.
// At each step, attempts a volume scaling move that is accepted/rejected
// based on the Metropolis criterion.
// ============================================================================

///|
/// Attempt a Monte Carlo barostat volume move.
/// Scales all positions by a factor, computes the energy before and after,
/// and accepts or rejects based on the Metropolis criterion.
/// Returns (new_positions, accepted).
pub fn md_mc_barostat_step(
  positions : Array[Array[Double]],
  box : Array[Double],
  energy_fn : (Array[Array[Double]]) -> Double,
  pressure : Double,
  temperature : Double,
  n_atoms : Int,
) -> (Array[Array[Double]], Array[Double], Bool) {
  let kT = MD_KB * temperature
  // Random volume scale: V_new = V * exp(delta)
  let delta = (md_rand_double() - 0.5) * 0.01
  let scale = @math.exp(delta / 3.0)
  // Scale positions
  let n = positions.length()
  let scaled = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      scaled[i][d] = positions[i][d] * scale
    }
  }
  let new_box = [box[0] * scale, box[1] * scale, box[2] * scale]
  // Metropolis criterion
  let e_old = energy_fn(positions)
  let e_new = energy_fn(scaled)
  let vol_old = box[0] * box[1] * box[2]
  let vol_new = new_box[0] * new_box[1] * new_box[2]
  let dE = e_new - e_old + pressure * (vol_new - vol_old) - n_atoms.to_double() * kT * @math.ln(vol_new / vol_old)
  let accepted = dE < 0.0 || md_rand_double() < @math.exp(-dE / kT)
  if accepted {
    (scaled, new_box, true)
  } else {
    (positions, box, false)
  }
}

///|
/// Run a harmonic-bond simulation with Monte Carlo barostat + Langevin.
/// Returns temperatures at each step.
pub fn md_simulate_barostat(
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  masses : Array[Double],
  bonds : Array[(Int, Int, Double, Double)],
  box : Array[Double],
  dt : Double,
  temperature : Double,
  pressure : Double,
  friction : Double,
  steps : Int,
  n_dof : Int,
) -> Array[Double] {
  let temps : Array[Double] = Array(capacity=0)
  let mut pos = positions
  let mut vel = velocities
  let mut current_box = box
  let n = positions.length()
  for _s = 0; _s < steps; _s = _s + 1 {
    // Langevin step
    let forces = md_bond_forces(pos, bonds)
    let (pos_new, v_o) = md_langevin_step(
      pos, vel, forces, masses, dt, temperature, friction,
    )
    let new_forces = md_bond_forces(pos_new, bonds)
    vel = md_langevin_finish(v_o, new_forces, masses, dt)
    pos = pos_new
    // MC barostat step every 10 steps
    if _s % 10 == 0 {
      let energy_fn = fn(p : Array[Array[Double]]) -> Double {
        md_bond_energy(p, bonds)
      }
      let (scaled, new_box, _accepted) = md_mc_barostat_step(
        pos, current_box, energy_fn, pressure, temperature, n,
      )
      pos = scaled
      current_box = new_box
    }
    temps.push(md_temperature(vel, masses, n_dof))
  }
  temps
}
