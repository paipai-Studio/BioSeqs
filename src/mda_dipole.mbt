///|
/// Trajectory dipole moments (MDAnalysis `dipole moments`): the molecular
/// dipole vector of a charge distribution, per-frame dipoles over a
/// trajectory, unit conversion to Debye, magnitude, and the dipole
/// autocorrelation function.

///|
/// Dipole moment `mu = sum_i q_i r_i` of a set of point charges
/// (MDAnalysis dipole moments), returned as `[mux, muy, muz]` in units of
/// `e * Angstrom`. Raises when `coords` and `charges` have different
/// lengths.
pub fn mda_dipole(
  coords : Array[Array[Double]],
  charges : Array[Double],
) -> Array[Double] raise {
  if coords.length() != charges.length() {
    raise SeqError("mda_dipole: coords and charges must have the same length")
  }
  let mut mx = 0.0
  let mut my = 0.0
  let mut mz = 0.0
  for i = 0; i < coords.length(); i = i + 1 {
    let q = charges[i]
    let r = coords[i]
    mx = mx + q * r[0]
    my = my + q * r[1]
    mz = mz + q * r[2]
  }
  [mx, my, mz]
}

///|
/// Per-frame dipole moments over a trajectory (MDAnalysis dipole
/// moments): one `[mux, muy, muz]` row (units `e * Angstrom`) per frame,
/// each obtained with `mda_dipole`. Raises when a frame length differs
/// from the charge array length.
pub fn mda_dipole_traj(
  traj : Array[Array[Array[Double]]],
  charges : Array[Double],
) -> Array[Array[Double]] raise {
  let out : Array[Array[Double]] = Array(capacity=traj.length())
  for frame in traj {
    out.push(mda_dipole(frame, charges))
  }
  out
}

///|
/// Convert a dipole from `e * Angstrom` to Debye
/// (`1 e*A = 4.80320425 D`).
pub fn mda_dipole_debye(mu : Array[Double]) -> Array[Double] {
  [mu[0] * 4.80320425, mu[1] * 4.80320425, mu[2] * 4.80320425]
}

///|
/// Magnitude `|mu| = sqrt(mu . mu)` of a dipole vector.
pub fn mda_dipole_magnitude(mu : Array[Double]) -> Double {
  (mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2]).sqrt()
}

///|
/// Normalized dipole autocorrelation function
/// `C(t) = <mu(tau) . mu(tau + t)>_tau / <mu(tau) . mu(tau)>_tau`,
/// averaged over all valid time origins `tau` for each lag
/// `t = 0 .. n-1`. `C(0) = 1` for any nonzero dipole trajectory; when the
/// denominator vanishes (all dipoles zero) an all-zero array is returned.
pub fn mda_dipole_correlation(dipoles : Array[Array[Double]]) -> Array[Double] {
  let n = dipoles.length()
  let c = Array::make(n, 0.0)
  if n == 0 {
    return c
  }
  let mut denom = 0.0
  for mu in dipoles {
    denom = denom + mda_vec3_dot(mu, mu)
  }
  denom = denom / n.to_double()
  if denom == 0.0 {
    return c
  }
  for t = 0; t < n; t = t + 1 {
    let m = n - t
    let mut s = 0.0
    for tau = 0; tau < m; tau = tau + 1 {
      s = s + mda_vec3_dot(dipoles[tau], dipoles[tau + t])
    }
    c[t] = s / m.to_double() / denom
  }
  c
}
