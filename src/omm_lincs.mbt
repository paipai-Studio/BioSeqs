/// omm_lincs.mbt — LINCS constraint solver.
///
/// OpenMM equivalent: `CustomIntegrator` constraint handling / GROMACS
/// LINCS (Hess et al. 1997).
///
/// Restores holonomic distance constraints r_ij = d after an unconstrained
/// integration step, for both positions and velocities. The coupled
/// constraint system
///
///   A lambda = sigma,  A = B M^-1 B^T
///
/// (B the constraint Jacobian, sigma the violations) is solved exactly by
/// Gaussian elimination with partial pivoting; the nonlinearity is handled
/// by re-linearizing the bond geometry for a few outer iterations, which
/// is the LINCS update scheme. The correction
///
///   dx_i = +lambda_c * r_hat_c / m_i,   dx_j = -lambda_c * r_hat_c / m_j
///
/// is mass-weighted and preserves the center of mass. Velocities lose
/// their radial component along constrained bonds through the same solve
/// applied to B v.
///
/// Units: A (length), ps (time), amu (mass).

///|
/// One distance constraint (0-based particle indices + target length).
pub(all) struct OmmLincsConstraint {
  i : Int32
  j : Int32
  d : Double
}

///|
/// Solve A x = b by Gaussian elimination with partial pivoting (A is
/// s x s, symmetric positive-ish for valid constraint sets).
fn omm_lincs_solve_linear(
  a : Array[Array[Double]],
  b : Array[Double],
) -> Array[Double] {
  let s = b.length()
  let x : Array[Double] = Array::make(s, 0.0)
  if s == 0 {
    return x
  }
  // Forward elimination on an augmented copy.
  for col = 0; col < s; col = col + 1 {
    // Pivot search.
    let mut piv = col
    let mut best = a[col][col].abs()
    for r = col + 1; r < s; r = r + 1 {
      if a[r][col].abs() > best {
        best = a[r][col].abs()
        piv = r
      }
    }
    if best < 1.0e-14 {
      continue
    }
    if piv != col {
      let tmp = a[col]
      a[col] = a[piv]
      a[piv] = tmp
      let tb = b[col]
      b[col] = b[piv]
      b[piv] = tb
    }
    let pivot = a[col][col]
    for r = col + 1; r < s; r = r + 1 {
      let factor = a[r][col] / pivot
      if factor != 0.0 {
        for k = col; k < s; k = k + 1 {
          a[r][k] = a[r][k] - factor * a[col][k]
        }
        b[r] = b[r] - factor * b[col]
      }
    }
  }
  // Back substitution.
  for r = s - 1; r >= 0; r = r - 1 {
    let mut acc = b[r]
    for k = r + 1; k < s; k = k + 1 {
      acc = acc - a[r][k] * x[k]
    }
    x[r] = if a[r][r].abs() > 1.0e-14 { acc / a[r][r] } else { 0.0 }
  }
  x
}

///|
/// Project positions and velocities onto all constraints.
///
/// `constraints` entries are (i, j, d). After an unconstrained step the
/// bond lengths are restored to `d` (positions) and the radial relative
/// velocities along constrained bonds are removed (velocities). The
/// projection is mass-weighted and preserves the center of mass.
///
/// Returns (projected_positions, projected_velocities); the input arrays
/// are not modified.
pub fn omm_lincs_apply(
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  constraints : Array[OmmLincsConstraint],
  masses : Array[Double],
  n_iter~ : Int32 = 3,
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let n = positions.length()
  let s = constraints.length()
  let pos : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  let vel : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  for i = 0; i < n; i = i + 1 {
    pos[i] = [positions[i][0], positions[i][1], positions[i][2]]
    vel[i] = [velocities[i][0], velocities[i][1], velocities[i][2]]
  }
  if s == 0 || n == 0 {
    return (pos, vel)
  }
  fn mass_of(idx : Int32) -> Double {
    if idx < masses.length() && idx >= 0 {
      masses[idx]
    } else {
      1.0
    }
  }

  // Bond unit vectors (r_hat points from i to j).
  let ux : Array[Double] = Array::make(s, 0.0)
  let uy : Array[Double] = Array::make(s, 0.0)
  let uz : Array[Double] = Array::make(s, 0.0)
  fn refresh() -> Unit {
    for c = 0; c < s; c = c + 1 {
      let i = constraints[c].i
      let j = constraints[c].j
      let rx = pos[j][0] - pos[i][0]
      let ry = pos[j][1] - pos[i][1]
      let rz = pos[j][2] - pos[i][2]
      let len = @math.pow(rx * rx + ry * ry + rz * rz, 0.5)
      if len > 1.0e-12 {
        ux[c] = rx / len
        uy[c] = ry / len
        uz[c] = rz / len
      } else {
        ux[c] = 1.0
        uy[c] = 0.0
        uz[c] = 0.0
      }
    }
  }

  // Coupling matrix A_cd = B_c M^-1 B_d^T. B_c acts as -r_c on i_c and
  // +r_c on j_c; entries only couple constraints that share an atom.
  fn build_a() -> Array[Array[Double]] {
    let a : Array[Array[Double]] = Array::make(s, [])
    for c = 0; c < s; c = c + 1 {
      let row : Array[Double] = Array::make(s, 0.0)
      for d2 = 0; d2 < s; d2 = d2 + 1 {
        let mut acc = 0.0
        // shared atom i_c
        if constraints[c].i == constraints[d2].i {
          acc = acc +
            (ux[c] * ux[d2] + uy[c] * uy[d2] + uz[c] * uz[d2]) /
            mass_of(constraints[c].i)
        }
        if constraints[c].i == constraints[d2].j {
          acc = acc -
            (ux[c] * ux[d2] + uy[c] * uy[d2] + uz[c] * uz[d2]) /
            mass_of(constraints[c].i)
        }
        // shared atom j_c
        if constraints[c].j == constraints[d2].i {
          acc = acc -
            (ux[c] * ux[d2] + uy[c] * uy[d2] + uz[c] * uz[d2]) /
            mass_of(constraints[c].j)
        }
        if constraints[c].j == constraints[d2].j {
          acc = acc +
            (ux[c] * ux[d2] + uy[c] * uy[d2] + uz[c] * uz[d2]) /
            mass_of(constraints[c].j)
        }
        row[d2] = acc
      }
      a[c] = row
    }
    a
  }

  fn apply_lambda(lam : Array[Double], into_positions : Bool) -> Unit {
    for c = 0; c < s; c = c + 1 {
      let i = constraints[c].i
      let j = constraints[c].j
      let gi = lam[c] / mass_of(i)
      let gj = -lam[c] / mass_of(j)
      if into_positions {
        pos[i][0] = pos[i][0] + gi * ux[c]
        pos[i][1] = pos[i][1] + gi * uy[c]
        pos[i][2] = pos[i][2] + gi * uz[c]
        pos[j][0] = pos[j][0] + gj * ux[c]
        pos[j][1] = pos[j][1] + gj * uy[c]
        pos[j][2] = pos[j][2] + gj * uz[c]
      } else {
        vel[i][0] = vel[i][0] + gi * ux[c]
        vel[i][1] = vel[i][1] + gi * uy[c]
        vel[i][2] = vel[i][2] + gi * uz[c]
        vel[j][0] = vel[j][0] + gj * ux[c]
        vel[j][1] = vel[j][1] + gj * uy[c]
        vel[j][2] = vel[j][2] + gj * uz[c]
      }
    }
  }

  // Positions: re-linearize and solve exactly each outer iteration.
  for _it = 0; _it < n_iter; _it = _it + 1 {
    refresh()
    let sigma : Array[Double] = Array::make(s, 0.0)
    for c = 0; c < s; c = c + 1 {
      let i = constraints[c].i
      let j = constraints[c].j
      let dx = pos[j][0] - pos[i][0]
      let dy = pos[j][1] - pos[i][1]
      let dz = pos[j][2] - pos[i][2]
      let len = @math.pow(dx * dx + dy * dy + dz * dz, 0.5)
      sigma[c] = len - constraints[c].d
    }
    // Correction: x -= M^-1 B^T lambda with A lambda = sigma, matching the
    // apply sign convention (lambda positive shrinks overlong bonds).
    let a = build_a()
    let lam = omm_lincs_solve_linear(a, sigma)
    apply_lambda(lam, true)
  }
  // Velocities: enforce B v = 0 with the same coupled solve.
  for _it = 0; _it < n_iter; _it = _it + 1 {
    refresh()
    let sigma : Array[Double] = Array::make(s, 0.0)
    for c = 0; c < s; c = c + 1 {
      let i = constraints[c].i
      let j = constraints[c].j
      sigma[c] = (vel[j][0] - vel[i][0]) * ux[c] +
        (vel[j][1] - vel[i][1]) * uy[c] +
        (vel[j][2] - vel[i][2]) * uz[c]
    }
    let a = build_a()
    let lam = omm_lincs_solve_linear(a, sigma)
    apply_lambda(lam, false)
  }
  (pos, vel)
}
