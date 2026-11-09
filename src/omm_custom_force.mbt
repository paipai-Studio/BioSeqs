/// omm_custom_force.mbt — user-defined force terms with math expressions.
///
/// OpenMM equivalents: `CustomBondForce`, `CustomAngleForce`,
/// `CustomTorsionForce`, `CustomNonbondedForce`.
///
/// Each force stores an expression over one geometric variable:
/// - bond / nonbonded: `r` (distance in A)
/// - angle: `theta` (radians)
/// - torsion: `phi` (radians)
///
/// Expressions support the operator set and math functions of the custom
/// integrator engine (sqrt/exp/log/sin/cos/tan/abs/floor/ceil/square/step/
/// min/max/select). Global parameters and per-term parameters (added with
/// `add_*`) are visible as variables inside the expression; per-term
/// parameters override globals on name clash.
///
/// Energy is evaluated exactly; forces use a finite-difference derivative
/// of the energy with respect to the geometric variable (bond/nonbonded)
/// or Cartesian coordinates (angle/torsion). Units: kcal/mol, A, radians.

///|
/// One user-defined bond term: particle pair + per-term parameters.
pub(all) struct OmmCustomBondTerm {
  i : Int32
  j : Int32
  params : Array[(String, Double)]
}

///|
/// CustomBondForce: energy = sum of `expression(r)` over bond terms.
pub struct OmmCustomBondForce {
  expression : String
  globals : Array[(String, Double)]
  terms : Array[OmmCustomBondTerm]
}

///|
pub fn OmmCustomBondForce::new(expression~ : String) -> OmmCustomBondForce {
  { expression, globals: [], terms: [], }
}

///|
pub fn OmmCustomBondForce::add_global(
  self : OmmCustomBondForce,
  name : String,
  value : Double,
) -> Unit {
  self.globals.push((name, value))
}

///|
pub fn OmmCustomBondForce::add_bond(
  self : OmmCustomBondForce,
  i : Int32,
  j : Int32,
  params~ : Array[(String, Double)] = [],
) -> Unit {
  self.terms.push({ i, j, params, })
}

///|
/// Merge globals with per-term parameters (per-term wins on clash).
fn omm_custom_merged_params(
  globals : Array[(String, Double)],
  params : Array[(String, Double)],
) -> Array[(String, Double)] {
  let merged : Array[(String, Double)] = Array(capacity=0)
  for g in globals {
    let mut overridden = false
    for p in params {
      if p.0 == g.0 {
        merged.push(p)
        overridden = true
        break
      }
    }
    if not(overridden) {
      merged.push(g)
    }
  }
  for p in params {
    let mut already = false
    for g in globals {
      if g.0 == p.0 {
        already = true
        break
      }
    }
    if not(already) {
      merged.push(p)
    }
  }
  merged
}

///|
/// Distance between two points.
pub fn omm_pair_distance(a : Array[Double], b : Array[Double]) -> Double {
  let dx = b[0] - a[0]
  let dy = b[1] - a[1]
  let dz = b[2] - a[2]
  @math.pow(dx * dx + dy * dy + dz * dz, 0.5)
}

///|
/// Bond angle at vertex j (radians, in [0, pi]).
pub fn omm_angle_rad(
  pi : Array[Double],
  pj : Array[Double],
  pk : Array[Double],
) -> Double {
  let v1x = pi[0] - pj[0]
  let v1y = pi[1] - pj[1]
  let v1z = pi[2] - pj[2]
  let v2x = pk[0] - pj[0]
  let v2y = pk[1] - pj[1]
  let v2z = pk[2] - pj[2]
  let m1 = @math.pow(v1x * v1x + v1y * v1y + v1z * v1z, 0.5)
  let m2 = @math.pow(v2x * v2x + v2y * v2y + v2z * v2z, 0.5)
  if m1 <= 0.0 || m2 <= 0.0 {
    return 0.0
  }
  let ct = (v1x * v2x + v1y * v2y + v1z * v2z) / (m1 * m2)
  let c = if ct > 1.0 { 1.0 } else if ct < -1.0 { -1.0 } else { ct }
  @math.acos(c)
}

///|
/// Dihedral angle i-j-k-l (radians, in (-pi, pi]).
pub fn omm_dihedral_rad(
  p1 : Array[Double],
  p2 : Array[Double],
  p3 : Array[Double],
  p4 : Array[Double],
) -> Double {
  let b1x = p2[0] - p1[0]
  let b1y = p2[1] - p1[1]
  let b1z = p2[2] - p1[2]
  let b2x = p3[0] - p2[0]
  let b2y = p3[1] - p2[1]
  let b2z = p3[2] - p2[2]
  let b3x = p4[0] - p3[0]
  let b3y = p4[1] - p3[1]
  let b3z = p4[2] - p3[2]
  let m1x = b1y * b2z - b1z * b2y
  let m1y = b1z * b2x - b1x * b2z
  let m1z = b1x * b2y - b1y * b2x
  let m2x = b2y * b3z - b2z * b3y
  let m2y = b2z * b3x - b2x * b3z
  let m2z = b2x * b3y - b2y * b3x
  let m1_len = @math.pow(m1x * m1x + m1y * m1y + m1z * m1z, 0.5)
  let m2_len = @math.pow(m2x * m2x + m2y * m2y + m2z * m2z, 0.5)
  let b2_len = @math.pow(b2x * b2x + b2y * b2y + b2z * b2z, 0.5)
  if m1_len <= 1.0e-12 || m2_len <= 1.0e-12 || b2_len <= 1.0e-12 {
    return 0.0
  }
  let cos_phi = (m1x * m2x + m1y * m2y + m1z * m2z) / (m1_len * m2_len)
  let c = if cos_phi > 1.0 {
    1.0
  } else if cos_phi < -1.0 {
    -1.0
  } else {
    cos_phi
  }
  let cross_x = m1y * b2z - m1z * b2y
  let cross_y = m1z * b2x - m1x * b2z
  let cross_z = m1x * b2y - m1y * b2x
  let sgn = if cross_x * m2x + cross_y * m2y + cross_z * m2z >= 0.0 {
    1.0
  } else {
    -1.0
  }
  sgn * @math.acos(c)
}

///|
/// Evaluate the bond-force expression at a distance r. The geometric
/// variable is exposed as `r` (and as `x` for convenience).
fn omm_custom_bond_eval(
  f : OmmCustomBondForce,
  r : Double,
  params : Array[(String, Double)],
) -> Double {
  let merged = omm_custom_merged_params(f.globals, params)
  merged.push(("r", r))
  omm_eval_expr(f.expression, r, 0.0, 0.0, 1.0, 0.0, merged)
}

///|
/// Total CustomBondForce energy.
pub fn omm_custom_bond_energy(
  f : OmmCustomBondForce,
  positions : Array[Array[Double]],
) -> Double {
  let mut e = 0.0
  for t in f.terms {
    if t.i < positions.length() && t.j < positions.length() {
      e = e +
        omm_custom_bond_eval(
          f,
          omm_pair_distance(positions[t.i], positions[t.j]),
          t.params,
        )
    }
  }
  e
}

///|
/// Forces from a CustomBondForce (analytic chain rule via numeric dE/dr).
pub fn omm_custom_bond_forces(
  f : OmmCustomBondForce,
  positions : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = positions.length()
  let out = md_new_forces(n)
  let h = 1.0e-5
  for t in f.terms {
    if t.i >= n || t.j >= n {
      continue
    }
    let a = positions[t.i]
    let b = positions[t.j]
    let r0 = omm_pair_distance(a, b)
    if r0 <= 0.0 {
      continue
    }
    let dedr = (
        omm_custom_bond_eval(f, r0 + h, t.params) -
        omm_custom_bond_eval(f, r0 - h, t.params)
      ) /
      (2.0 * h)
    // r-hat from i to j; F_i = +dE/dr * r-hat, F_j = -dE/dr * r-hat.
    let ux = (b[0] - a[0]) / r0
    let uy = (b[1] - a[1]) / r0
    let uz = (b[2] - a[2]) / r0
    out[t.i][0] = out[t.i][0] + dedr * ux
    out[t.i][1] = out[t.i][1] + dedr * uy
    out[t.i][2] = out[t.i][2] + dedr * uz
    out[t.j][0] = out[t.j][0] - dedr * ux
    out[t.j][1] = out[t.j][1] - dedr * uy
    out[t.j][2] = out[t.j][2] - dedr * uz
  }
  out
}

///|
/// One user-defined angle term: particle triple (j is the vertex).
pub(all) struct OmmCustomAngleTerm {
  i : Int32
  j : Int32
  k : Int32
  params : Array[(String, Double)]
}

///|
/// CustomAngleForce: energy = sum of `expression(theta)` over angle terms.
pub struct OmmCustomAngleForce {
  expression : String
  globals : Array[(String, Double)]
  terms : Array[OmmCustomAngleTerm]
}

///|
pub fn OmmCustomAngleForce::new(expression~ : String) -> OmmCustomAngleForce {
  { expression, globals: [], terms: [], }
}

///|
pub fn OmmCustomAngleForce::add_global(
  self : OmmCustomAngleForce,
  name : String,
  value : Double,
) -> Unit {
  self.globals.push((name, value))
}

///|
pub fn OmmCustomAngleForce::add_angle(
  self : OmmCustomAngleForce,
  i : Int32,
  j : Int32,
  k : Int32,
  params~ : Array[(String, Double)] = [],
) -> Unit {
  self.terms.push({ i, j, k, params, })
}

///|
/// Evaluate the angle-force expression at a theta (radians). The geometric
/// variable is exposed as `theta`.
fn omm_custom_angle_eval(
  f : OmmCustomAngleForce,
  theta : Double,
  params : Array[(String, Double)],
) -> Double {
  let merged = omm_custom_merged_params(f.globals, params)
  merged.push(("theta", theta))
  omm_eval_expr(f.expression, 0.0, 0.0, 0.0, 1.0, 0.0, merged)
}

///|
/// Total CustomAngleForce energy.
pub fn omm_custom_angle_energy(
  f : OmmCustomAngleForce,
  positions : Array[Array[Double]],
) -> Double {
  let mut e = 0.0
  for t in f.terms {
    if t.i < positions.length() &&
      t.j < positions.length() &&
      t.k < positions.length() {
      e = e +
        omm_custom_angle_eval(
          f,
          omm_angle_rad(positions[t.i], positions[t.j], positions[t.k]),
          t.params,
        )
    }
  }
  e
}

///|
/// Cartesian forces from a CustomAngleForce (finite differences).
pub fn omm_custom_angle_forces(
  f : OmmCustomAngleForce,
  positions : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = positions.length()
  let out = md_new_forces(n)
  let h = 1.0e-5
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      let plus = Array::make(n, [0.0, 0.0, 0.0])
      let minus = Array::make(n, [0.0, 0.0, 0.0])
      for j = 0; j < n; j = j + 1 {
        plus[j] = [positions[j][0], positions[j][1], positions[j][2]]
        minus[j] = [positions[j][0], positions[j][1], positions[j][2]]
      }
      plus[i][d] = plus[i][d] + h
      minus[i][d] = minus[i][d] - h
      out[i][d] = -(omm_custom_angle_energy(f, plus) -
        omm_custom_angle_energy(f, minus)) /
        (2.0 * h)
    }
  }
  out
}

///|
/// One user-defined torsion term: particle quadruple i-j-k-l.
pub(all) struct OmmCustomTorsionTerm {
  i : Int32
  j : Int32
  k : Int32
  l : Int32
  params : Array[(String, Double)]
}

///|
/// CustomTorsionForce: energy = sum of `expression(phi)` over terms.
pub struct OmmCustomTorsionForce {
  expression : String
  globals : Array[(String, Double)]
  terms : Array[OmmCustomTorsionTerm]
}

///|
pub fn OmmCustomTorsionForce::new(
  expression~ : String,
) -> OmmCustomTorsionForce {
  { expression, globals: [], terms: [], }
}

///|
pub fn OmmCustomTorsionForce::add_global(
  self : OmmCustomTorsionForce,
  name : String,
  value : Double,
) -> Unit {
  self.globals.push((name, value))
}

///|
pub fn OmmCustomTorsionForce::add_torsion(
  self : OmmCustomTorsionForce,
  i : Int32,
  j : Int32,
  k : Int32,
  l : Int32,
  params~ : Array[(String, Double)] = [],
) -> Unit {
  self.terms.push({ i, j, k, l, params, })
}

///|
/// Total CustomTorsionForce energy.
pub fn omm_custom_torsion_energy(
  f : OmmCustomTorsionForce,
  positions : Array[Array[Double]],
) -> Double {
  let mut e = 0.0
  for t in f.terms {
    if t.i < positions.length() &&
      t.j < positions.length() &&
      t.k < positions.length() &&
      t.l < positions.length() {
      let phi = omm_dihedral_rad(
        positions[t.i],
        positions[t.j],
        positions[t.k],
        positions[t.l],
      )
      let merged = omm_custom_merged_params(f.globals, t.params)
      merged.push(("phi", phi))
      e = e + omm_eval_expr(f.expression, 0.0, 0.0, 0.0, 1.0, 0.0, merged)
    }
  }
  e
}

///|
/// Cartesian forces from a CustomTorsionForce (finite differences).
pub fn omm_custom_torsion_forces(
  f : OmmCustomTorsionForce,
  positions : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = positions.length()
  let out = md_new_forces(n)
  let h = 1.0e-5
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      let plus = Array::make(n, [0.0, 0.0, 0.0])
      let minus = Array::make(n, [0.0, 0.0, 0.0])
      for j = 0; j < n; j = j + 1 {
        plus[j] = [positions[j][0], positions[j][1], positions[j][2]]
        minus[j] = [positions[j][0], positions[j][1], positions[j][2]]
      }
      plus[i][d] = plus[i][d] + h
      minus[i][d] = minus[i][d] - h
      out[i][d] = -(omm_custom_torsion_energy(f, plus) -
        omm_custom_torsion_energy(f, minus)) /
        (2.0 * h)
    }
  }
  out
}

///|
/// One user-defined nonbonded pair: particle pair + per-pair parameters.
pub(all) struct OmmCustomNonbondedTerm {
  i : Int32
  j : Int32
  params : Array[(String, Double)]
}

///|
/// CustomNonbondedForce: energy = sum of `expression(r)` over pair terms.
pub struct OmmCustomNonbondedForce {
  expression : String
  globals : Array[(String, Double)]
  terms : Array[OmmCustomNonbondedTerm]
}

///|
pub fn OmmCustomNonbondedForce::new(
  expression~ : String,
) -> OmmCustomNonbondedForce {
  { expression, globals: [], terms: [], }
}

///|
pub fn OmmCustomNonbondedForce::add_global(
  self : OmmCustomNonbondedForce,
  name : String,
  value : Double,
) -> Unit {
  self.globals.push((name, value))
}

///|
pub fn OmmCustomNonbondedForce::add_pair(
  self : OmmCustomNonbondedForce,
  i : Int32,
  j : Int32,
  params~ : Array[(String, Double)] = [],
) -> Unit {
  self.terms.push({ i, j, params, })
}

///|
/// Total CustomNonbondedForce energy.
pub fn omm_custom_nonbonded_energy(
  f : OmmCustomNonbondedForce,
  positions : Array[Array[Double]],
) -> Double {
  let mut e = 0.0
  for t in f.terms {
    if t.i < positions.length() && t.j < positions.length() {
      let r = omm_pair_distance(positions[t.i], positions[t.j])
      let merged = omm_custom_merged_params(f.globals, t.params)
      merged.push(("r", r))
      e = e + omm_eval_expr(f.expression, r, 0.0, 0.0, 1.0, 0.0, merged)
    }
  }
  e
}

///|
/// Forces from a CustomNonbondedForce (analytic chain rule via dE/dr).
pub fn omm_custom_nonbonded_forces(
  f : OmmCustomNonbondedForce,
  positions : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = positions.length()
  let out = md_new_forces(n)
  let h = 1.0e-5
  for t in f.terms {
    if t.i >= n || t.j >= n {
      continue
    }
    let a = positions[t.i]
    let b = positions[t.j]
    let r0 = omm_pair_distance(a, b)
    if r0 <= 0.0 {
      continue
    }
    let merged = omm_custom_merged_params(f.globals, t.params)
    merged.push(("r", r0 + h))
    let e_plus = omm_eval_expr(f.expression, r0 + h, 0.0, 0.0, 1.0, 0.0, merged)
    merged[merged.length() - 1] = ("r", r0 - h)
    let e_minus = omm_eval_expr(
      f.expression,
      r0 - h,
      0.0,
      0.0,
      1.0,
      0.0,
      merged,
    )
    let dedr = (e_plus - e_minus) / (2.0 * h)
    let ux = (b[0] - a[0]) / r0
    let uy = (b[1] - a[1]) / r0
    let uz = (b[2] - a[2]) / r0
    out[t.i][0] = out[t.i][0] + dedr * ux
    out[t.i][1] = out[t.i][1] + dedr * uy
    out[t.i][2] = out[t.i][2] + dedr * uz
    out[t.j][0] = out[t.j][0] - dedr * ux
    out[t.j][1] = out[t.j][1] - dedr * uy
    out[t.j][2] = out[t.j][2] - dedr * uz
  }
  out
}
