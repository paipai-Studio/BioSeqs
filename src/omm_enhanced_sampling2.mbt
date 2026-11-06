/// omm_enhanced_sampling2.mbt — Steered MD + Metadynamics + REST2
/// Items 12-14: Advanced enhanced sampling methods

///|
// === Steered MD ===

///|
/// Steered MD: apply a moving harmonic restraint to pull along a CV.
/// The center moves at constant velocity: center(t) = center0 + v * t
pub fn omm_steered_md_force(
  cv_value : Double, // current collective variable value
  center0 : Double, // initial center
  velocity : Double, // pulling velocity (units/time)
  time : Double, // current simulation time
  k : Double, // spring constant
) -> Double {
  let center = center0 + velocity * time
  -k * (cv_value - center)
}

///|
/// Steered MD energy
pub fn omm_steered_md_energy(
  cv_value : Double,
  center0 : Double,
  velocity : Double,
  time : Double,
  k : Double,
) -> Double {
  let center = center0 + velocity * time
  0.5 * k * (cv_value - center) * (cv_value - center)
}

///|
/// Compute work accumulated during steered MD (Jarzynski equality).
pub fn omm_steered_md_work(
  cv_values : Array[Double], // CV at each step
  times : Array[Double], // time at each step
  center0 : Double,
  velocity : Double,
  k : Double,
) -> Double {
  let n = cv_values.length()
  if n < 2 {
    return 0.0
  }
  let mut work = 0.0
  for i = 1; i < n; i = i + 1 {
    let f_prev = omm_steered_md_force(
      cv_values[i - 1],
      center0,
      velocity,
      times[i - 1],
      k,
    )
    let dcv = cv_values[i] - cv_values[i - 1]
    work = work + f_prev * dcv
  }
  work
}

///|
// === Metadynamics ===

///|
/// Metadynamics bias potential: sum of Gaussian hills deposited along CV.
pub struct OmmMetaHill {
  center : Array[Double] // hill center in CV space
  height : Double // hill height
  width : Array[Double] // hill width (sigma) per CV
}

///|
/// Metadynamics state
pub struct OmmMetadynamics {
  hills : Array[OmmMetaHill]
  n_cv : Int32
  bias_factor : Double // well-tempered: deltaT/T
  current_bias : Double
}

///|
/// Compute metadynamics bias potential at a given CV point.
pub fn omm_metad_bias(
  meta : OmmMetadynamics,
  cv_values : Array[Double],
) -> Double {
  let mut bias = 0.0
  for hill in meta.hills {
    // Multi-dimensional Gaussian
    let mut arg2 = 0.0
    for d = 0; d < meta.n_cv; d = d + 1 {
      if d < cv_values.length() &&
        d < hill.center.length() &&
        d < hill.width.length() {
        let diff = cv_values[d] - hill.center[d]
        let sigma = hill.width[d]
        arg2 = arg2 + diff * diff / (2.0 * sigma * sigma)
      }
    }
    bias = bias + hill.height * @math.exp(-arg2)
  }
  bias
}

///|
/// Compute metadynamics force on a CV.
pub fn omm_metad_force(
  meta : OmmMetadynamics,
  cv_values : Array[Double],
  cv_index : Int32, // which CV to compute force for
) -> Double {
  let mut force = 0.0
  for hill in meta.hills {
    let mut arg2 = 0.0
    let mut diff_cv = 0.0
    let mut sigma_cv = 1.0
    for d = 0; d < meta.n_cv; d = d + 1 {
      if d < cv_values.length() &&
        d < hill.center.length() &&
        d < hill.width.length() {
        let diff = cv_values[d] - hill.center[d]
        let sigma = hill.width[d]
        arg2 = arg2 + diff * diff / (2.0 * sigma * sigma)
        if d == cv_index {
          diff_cv = diff
          sigma_cv = sigma
        }
      }
    }
    // dV/dcv = -height * exp(-arg2) * (-diff/sigma^2) = height * diff/sigma^2 * exp(-arg2)
    force = force +
      hill.height * diff_cv / (sigma_cv * sigma_cv) * @math.exp(-arg2)
  }
  // Force is -dV/dcv
  -force
}

///|
/// Deposit a new Gaussian hill (well-tempered metadynamics).
pub fn omm_metad_deposit_hill(
  meta : OmmMetadynamics,
  cv_values : Array[Double],
  hill_height : Double,
  hill_width : Array[Double],
) -> OmmMetadynamics {
  // Well-tempered: height = height0 * exp(-bias / (kB * deltaT))
  let current_bias = omm_metad_bias(meta, cv_values)
  let scaled_height = if meta.bias_factor > 0.0 {
    hill_height * @math.exp(-current_bias / meta.bias_factor)
  } else {
    hill_height
  }
  let new_hill = OmmMetaHill::{
    center: cv_values.copy(),
    height: scaled_height,
    width: hill_width,
  }
  let new_hills : Array[OmmMetaHill] = Array(capacity=meta.hills.length() + 1)
  for h in meta.hills {
    new_hills.push(h)
  }
  new_hills.push(new_hill)
  OmmMetadynamics::{
    hills: new_hills,
    n_cv: meta.n_cv,
    bias_factor: meta.bias_factor,
    current_bias: current_bias + scaled_height,
  }
}

///|
// === REST2 (Replica Exchange with Solute Tempering) ===

///|
/// REST2: scale solute-solute and solute-solvent interactions by a factor.
/// factor = (T_ref / T)^2 for solute-solute, (T_ref / T) for solute-solvent.
pub fn omm_rest2_scale_factors(
  temp_ref : Double,
  temp : Double,
) -> (Double, Double) {
  let ratio = temp_ref / temp
  (ratio * ratio, ratio) // (solute-solute scale, solute-solvent scale)
}

///|
/// REST2: scale energy components.
pub fn omm_rest2_scaled_energy(
  e_solute_solute : Double,
  e_solute_solvent : Double,
  e_solvent_solvent : Double,
  temp_ref : Double,
  temp : Double,
) -> Double {
  let (ss_scale, sv_scale) = omm_rest2_scale_factors(temp_ref, temp)
  ss_scale * e_solute_solute + sv_scale * e_solute_solvent + e_solvent_solvent
}

///|
/// REST2 replica exchange acceptance probability.
pub fn omm_rest2_swap_accept(
  e1_ss : Double, // replica 1 solute-solute
  e1_sv : Double, // replica 1 solute-solvent
  e1_vv : Double, // replica 1 solvent-solvent
  e2_ss : Double,
  e2_sv : Double,
  e2_vv : Double,
  temp1 : Double,
  temp2 : Double,
  temp_ref~ : Double = 300.0,
) -> Bool {
  // Scaled energies at each temperature
  let e1 = omm_rest2_scaled_energy(e1_ss, e1_sv, e1_vv, temp_ref, temp1)
  let e2 = omm_rest2_scaled_energy(e2_ss, e2_sv, e2_vv, temp_ref, temp2)
  // After swap
  let e1_swapped = omm_rest2_scaled_energy(e2_ss, e2_sv, e2_vv, temp_ref, temp1)
  let e2_swapped = omm_rest2_scaled_energy(e1_ss, e1_sv, e1_vv, temp_ref, temp2)
  let delta = e1_swapped + e2_swapped - (e1 + e2)
  let beta1 = 1.0 / (1.380649e-23 * temp1 + 0.0000000001)
  let beta2 = 1.0 / (1.380649e-23 * temp2 + 0.0000000001)
  let arg = (beta1 - beta2) * delta
  arg > 0.0
}

///|
pub fn omm_enhanced_sampling2_dummy() -> Unit {
  ignore(OmmMetaHill::{
    center: Array(capacity=0),
    height: 0.0,
    width: Array(capacity=0),
  })
  ignore(OmmMetadynamics::{
    hills: Array(capacity=0),
    n_cv: 0,
    bias_factor: 0.0,
    current_bias: 0.0,
  })
}
