/// Free energy calculations: thermodynamic integration (TI) and FEP.
///
/// OpenMM equivalent: `AlchemicalCalculator`, `TIIntegrator`
///
/// Computes free energy differences via alchemical transformations
/// using a coupling parameter lambda that interpolates between states.

///|
/// Thermodynamic integration: dG = integral_0^1 <dU/dlambda> dlambda
///
/// energy_fn: function(positions, lambda) -> energy
/// positions: atomic positions (fixed for this calculation)
/// n_lambda: number of lambda windows
/// Returns free energy difference (G(1) - G(0))
pub fn omm_thermodynamic_integration(
  positions : Array[Array[Double]],
  energy_fn : (Array[Array[Double]], Double) -> Double,
  n_lambda~ : Int32 = 21,
) -> Double {
  if n_lambda < 2 {
    return 0.0
  }
  // Compute <dU/dlambda> at each lambda window via finite difference
  let dl = 1.0 / (n_lambda - 1).to_double()
  let delta = dl * 0.001 // finite difference step
  let dudl : Array[Double] = Array::make(n_lambda, 0.0)
  for i = 0; i < n_lambda; i = i + 1 {
    let lam = i.to_double() * dl
    let e_plus = energy_fn(positions, lam + delta)
    let e_minus = energy_fn(positions, lam - delta)
    dudl[i] = (e_plus - e_minus) / (2.0 * delta)
  }
  // Integrate using trapezoidal rule
  let mut integral = 0.0
  for i = 0; i < n_lambda - 1; i = i + 1 {
    integral = integral + (dudl[i] + dudl[i + 1]) * 0.5 * dl
  }
  integral
}

///|
/// Free Energy Perturbation (FEP): dG = -kT * ln(<exp(-beta*dU)>_A)
///
/// energies_A: energies of samples from state A
/// energies_B: energies of the same samples evaluated in state B
pub fn omm_fep(
  energies_A : Array[Double],
  energies_B : Array[Double],
  temperature~ : Double = 300.0,
) -> Double {
  let n = if energies_A.length() < energies_B.length() {
    energies_A.length()
  } else {
    energies_B.length()
  }
  if n == 0 {
    return 0.0
  }
  let kT = 0.001987204 * temperature // kcal/mol
  let beta = 1.0 / kT
  // dG = -kT * ln( <exp(-beta*(E_B - E_A))>_A )
  let mut sum = 0.0
  for i = 0; i < n; i = i + 1 {
    let dU = energies_B[i] - energies_A[i]
    sum = sum + @math.exp(-beta * dU)
  }
  -kT * @math.ln(sum / n.to_double())
}

///|
/// Bennett Acceptance Ratio (BAR) estimator.
/// Iteratively solves for dG that satisfies the BAR equation.
pub fn omm_bar(
  energies_A_in_A : Array[Double], // samples from A, energy in A
  energies_B_in_A : Array[Double], // samples from A, energy in B
  energies_A_in_B : Array[Double], // samples from B, energy in A
  energies_B_in_B : Array[Double], // samples from B, energy in B
  temperature~ : Double = 300.0,
) -> Double {
  let n_A = energies_A_in_A.length()
  let n_B = energies_A_in_B.length()
  if n_A == 0 || n_B == 0 {
    return 0.0
  }
  let kT = 0.001987204 * temperature
  let beta = 1.0 / kT
  // Initial guess: FEP estimate
  let mut dG = 0.0
  // Iterate BAR equation
  for iter = 0; iter < 100; iter = iter + 1 {
    let c_val = dG + kT * @math.ln(n_B.to_double() / n_A.to_double())
    // f(x) = 1/(1+exp(beta*x))
    let mut sum_A = 0.0
    for i = 0; i < n_A; i = i + 1 {
      let dU = energies_B_in_A[i] - energies_A_in_A[i]
      sum_A = sum_A + 1.0 / (1.0 + @math.exp(beta * (dU - c_val)))
    }
    let mut sum_B = 0.0
    for i = 0; i < n_B; i = i + 1 {
      let dU = energies_B_in_B[i] - energies_A_in_B[i]
      sum_B = sum_B + 1.0 / (1.0 + @math.exp(-beta * (dU - c_val)))
    }
    let ratio = sum_A / sum_B
    let dG_new = c_val - kT * @math.ln(ratio)
    let diff = dG_new - dG
    dG = dG_new
    if diff < 0.0 {
      if -diff < 0.0001 {
        break
      }
    } else if diff < 0.0001 {
      break
    }
  }
  dG
}

///|
/// Alchemical lambda schedule (soft-core to avoid singularities).
/// U(lambda) = (1-lambda)*U_A + lambda*U_B
/// With soft-core: r_eff = (r^6 + sigma^6*(1-lambda)*lambda)^(1/6)
pub fn omm_alchemical_energy(
  r : Double,
  epsilon_A : Double,
  sigma_A : Double,
  epsilon_B : Double,
  sigma_B : Double,
  lambda : Double,
  soft_core~ : Double = 0.5,
) -> Double {
  // Soft-core distance
  let r6 = r * r * r * r * r * r
  let sc = soft_core *
    sigma_A *
    sigma_A *
    sigma_A *
    sigma_A *
    sigma_A *
    sigma_A *
    (1.0 - lambda) *
    lambda
  let r_eff6 = r6 + sc
  let r_eff = r_eff6.pow(1.0 / 6.0)
  // LJ potential for state A and B
  let sr6_A = (sigma_A / r_eff).pow(6.0)
  let u_A = 4.0 * epsilon_A * (sr6_A * sr6_A - sr6_A)
  let sr6_B = (sigma_B / r_eff).pow(6.0)
  let u_B = 4.0 * epsilon_B * (sr6_B * sr6_B - sr6_B)
  (1.0 - lambda) * u_A + lambda * u_B
}
