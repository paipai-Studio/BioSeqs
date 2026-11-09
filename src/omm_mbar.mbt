/// omm_mbar.mbt — free-energy machinery: soft-core potentials, MBAR, and
/// lambda-schedule helpers.
///
/// OpenMM equivalents: soft-core interactions (CustomNonbondedForce free
/// energy recipes), `barostat`-free alchemical analysis (pymbar-style
/// MBAR), and adaptive lambda scheduling for TI/FEP.
///
/// Units: kcal/mol, A, e. MBAR works in reduced (dimensionless) energies:
/// pass u = beta * U (beta already folded in).

///|
/// Soft-core Lennard-Jones energy for an alchemical pair (Beutler-style).
///
///   r_eff^6 = r^6 + alpha_sc * sigma^6 * lambda^2
///   E = 4*eps*lambda * ((sigma/r_eff)^12 - (sigma/r_eff)^6)
///
/// lambda = 0 switches the interaction off smoothly; with alpha_sc -> 0
/// and lambda = 1 the standard LJ energy is recovered. The singularity at
/// r -> 0 is removed (finite energy).
pub fn omm_softcore_lj(
  sigma : Double,
  epsilon : Double,
  r : Double,
  lambda : Double,
  alpha_sc~ : Double = 0.5,
) -> Double {
  if lambda <= 0.0 {
    return 0.0
  }
  let r6 = @math.pow(r, 6.0)
  let s6 = @math.pow(sigma, 6.0)
  let r_eff6 = r6 + alpha_sc * s6 * lambda * lambda
  let sr6 = s6 / r_eff6
  4.0 * epsilon * lambda * (sr6 * sr6 - sr6)
}

///|
/// Soft-core Coulomb energy for an alchemical pair.
///
///   r_eff^2 = r^2 + alpha_sc * lambda^2
///   E = 332.06371 * q1*q2 * lambda / r_eff   (kcal/mol, A, e)
///
/// lambda = 0 switches the interaction off; lambda = 1 with small
/// alpha_sc approaches the bare Coulomb law.
pub fn omm_softcore_coulomb(
  q1q2 : Double,
  r : Double,
  lambda : Double,
  alpha_sc~ : Double = 0.5,
) -> Double {
  if lambda <= 0.0 {
    return 0.0
  }
  let r_eff = @math.pow(r * r + alpha_sc * lambda * lambda, 0.5)
  332.06371 * q1q2 * lambda / r_eff
}

///|
/// Log-sum-exp of an array (numerically stable).
fn omm_logsumexp(xs : Array[Double]) -> Double {
  if xs.length() == 0 {
    return -1.0e18
  }
  let mut mx = xs[0]
  for x in xs {
    if x > mx {
      mx = x
    }
  }
  if mx < -5.0e17 {
    return mx
  }
  let mut s = 0.0
  for x in xs {
    s = s + @math.exp(x - mx)
  }
  mx + @math.ln(s)
}

///|
/// Solve MBAR (multistate Bennett acceptance ratio) by log-domain
/// self-consistent iteration.
///
/// `u_nk[k][n][j]` is the reduced (beta-folded) potential of sample n from
/// state k evaluated at state j. Returns the dimensionless free energies
/// log_f with the convention log_f[0] = 0; free-energy differences are
/// Delta f_j = log_f[j] - log_f[i].
///
/// Iteration (Shirts & Chodera):
///   logC_j = logsumexp_{k,n} ( logF_k - u_kn[j] )
///   logF_j = -logC_j + logC_0          (anchor f_0 = 0)
pub fn omm_mbar_solve(
  u_nk : Array[Array[Array[Double]]],
  max_iter~ : Int32 = 1000,
  tol~ : Double = 1.0e-10,
) -> Array[Double] {
  let k_states = u_nk.length()
  let log_f : Array[Double] = Array::make(k_states, 0.0)
  if k_states == 0 {
    return log_f
  }
  for k = 0; k < k_states; k = k + 1 {
    log_f[k] = -@math.ln(u_nk[k].length().max(1).to_double())
  }
  let mut iter = 0
  while iter < max_iter {
    // Normalization of the mixture at every sample:
    // log c(x_kn) = logsumexp_j ( logF_j - u_kn[j] ).
    let log_c : Array[Array[Double]] = Array::make(k_states, [])
    for k = 0; k < k_states; k = k + 1 {
      let rows : Array[Double] = Array::make(u_nk[k].length(), 0.0)
      for n = 0; n < u_nk[k].length(); n = n + 1 {
        let terms : Array[Double] = Array(capacity=0)
        for j = 0; j < k_states; j = j + 1 {
          terms.push(log_f[j] - u_nk[k][n][j])
        }
        rows[n] = omm_logsumexp(terms)
      }
      log_c[k] = rows
    }
    // Self-consistent update from the likelihood stationarity condition
    // N_i = sum_{k,n} exp(logF_i - u_kn[i]) / c(x_kn):
    //   logF_i = log N_i - logsumexp_{k,n} ( -u_kn[i] - logc(x_kn) ).
    let new_f : Array[Double] = Array::make(k_states, 0.0)
    for i = 0; i < k_states; i = i + 1 {
      let terms : Array[Double] = Array(capacity=0)
      for k = 0; k < k_states; k = k + 1 {
        for n = 0; n < u_nk[k].length(); n = n + 1 {
          terms.push(-u_nk[k][n][i] - log_c[k][n])
        }
      }
      new_f[i] = @math.ln(u_nk[i].length().max(1).to_double()) -
        omm_logsumexp(terms)
    }
    // Anchor log f_0 = 0 and fold the overall constant in.
    let anchor = new_f[0]
    let mut delta = 0.0
    for j = 0; j < k_states; j = j + 1 {
      let v = new_f[j] - anchor
      let d = (v - log_f[j]).abs()
      if d > delta {
        delta = d
      }
      log_f[j] = v
    }
    if delta < tol {
      break
    }
    iter = iter + 1
  }
  log_f
}

///|
/// Uniform lambda schedule with `n` points from lo to hi (inclusive).
pub fn omm_lambda_linear(
  n~ : Int32,
  lo~ : Double = 0.0,
  hi~ : Double = 1.0,
) -> Array[Double] {
  let out : Array[Double] = Array(capacity=0)
  if n <= 1 {
    out.push(lo)
    return out
  }
  for i = 0; i < n; i = i + 1 {
    out.push(lo + (hi - lo) * i.to_double() / (n - 1).to_double())
  }
  out
}

///|
/// Rebalance a lambda schedule so that each interval covers the same
/// estimated free-energy change.
///
/// `gradients[i]` is the estimated <dU/dlambda> at lambdas[i]. The
/// cumulative free-energy profile is integrated by the trapezoid rule and
/// re-sampled at equal increments of `target` (at least 2 segments).
/// Endpoints are preserved (0 and 1).
pub fn omm_lambda_rebalance(
  lambdas : Array[Double],
  gradients : Array[Double],
  target~ : Double,
  lo~ : Double = 0.0,
  hi~ : Double = 1.0,
) -> Array[Double] {
  let n = lambdas.length()
  if n < 2 {
    return omm_lambda_linear(n=2, lo~, hi~)
  }
  // Cumulative free energy F_i (F_0 = 0).
  let cum : Array[Double] = Array::make(n, 0.0)
  for i = 1; i < n; i = i + 1 {
    let dl = (lambdas[i] - lambdas[i - 1]).abs()
    let g = 0.5 * (gradients[i].abs() + gradients[i - 1].abs())
    cum[i] = cum[i - 1] + g * dl
  }
  let total = cum[n - 1]
  if total <= 0.0 {
    return omm_lambda_linear(n~, lo~, hi~)
  }
  let segs = (total / target.max(1.0e-12)).ceil().to_int().max(2)
  let out : Array[Double] = Array(capacity=0)
  out.push(lo)
  let span = hi - lo
  for s = 1; s < segs; s = s + 1 {
    // Target cumulative value at fraction s/segs of the total dF.
    let f_target = total * s.to_double() / segs.to_double()
    // Invert the cumulative curve by linear interpolation.
    let mut i = 1
    while i < n && cum[i] < f_target {
      i = i + 1
    }
    let lam = if i >= n {
      hi
    } else {
      let c0 = cum[i - 1]
      let c1 = cum[i]
      let frac = if c1 > c0 { (f_target - c0) / (c1 - c0) } else { 0.0 }
      lambdas[i - 1] + frac * (lambdas[i] - lambdas[i - 1])
    }
    out.push(lo + span * lam)
  }
  out.push(hi)
  out
}
