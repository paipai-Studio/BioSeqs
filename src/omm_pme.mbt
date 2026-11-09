/// Particle Mesh Ewald (PME) for long-range electrostatics.
///
/// OpenMM equivalent: `NonbondedForce::PME`
///
/// Splits Coulomb into short-range (real-space, cutoff) + long-range
/// (reciprocal-space, FFT on grid). This implements the Ewald sum
/// with grid-based reciprocal sum (simplified SPME).

///|
/// PME parameters
pub struct OmmPMEParams {
  alpha : Double // Ewald splitting parameter
  cutoff : Double // real-space cutoff
  grid_x : Int32 // FFT grid dimensions
  grid_y : Int32
  grid_z : Int32
  box : Array[Double] // [lx, ly, lz]
}

///|
pub fn OmmPMEParams::new(
  cutoff : Double,
  grid_x~ : Int32 = 32,
  grid_y~ : Int32 = 32,
  grid_z~ : Int32 = 32,
  box~ : Array[Double] = [0.0, 0.0, 0.0],
) -> OmmPMEParams {
  // Optimal alpha: alpha = sqrt(pi) / cutoff (heuristic)
  let alpha = 3.141592653589793.sqrt() / cutoff
  { alpha, cutoff, grid_x, grid_y, grid_z, box, }
}

///|
/// Real-space Ewald sum (short-range part).
/// E_real = sum_{i<j} q_i*q_j * erfc(alpha*r) / r
pub fn omm_pme_real_space(
  positions : Array[Array[Double]],
  charges : Array[Double],
  params : OmmPMEParams,
) -> Double {
  let n = positions.length()
  let mut energy = 0.0
  let alpha = params.alpha
  let cutoff = params.cutoff
  let cutoff_sq = cutoff * cutoff
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let dx = positions[j][0] - positions[i][0]
      let dy = positions[j][1] - positions[i][1]
      let dz = positions[j][2] - positions[i][2]
      // PBC minimum image
      let mut rx = dx
      let mut ry = dy
      let mut rz = dz
      if params.box[0] > 0.0 {
        rx = rx - params.box[0] * (rx / params.box[0] + 0.5).floor()
        ry = ry - params.box[1] * (ry / params.box[1] + 0.5).floor()
        rz = rz - params.box[2] * (rz / params.box[2] + 0.5).floor()
      }
      let r_sq = rx * rx + ry * ry + rz * rz
      if r_sq > cutoff_sq || r_sq == 0.0 {
        continue
      }
      let r = r_sq.sqrt()
      // erfc(alpha*r) / r
      let arg = alpha * r
      let erfc_val = omm_erfc(arg)
      energy = energy + charges[i] * charges[j] * erfc_val / r
    }
  }
  energy
}

///|
/// Reciprocal-space Ewald sum (long-range part).
/// Uses direct summation (not FFT) for simplicity — O(N^2) but exact.
/// For production use, this would be replaced by FFT-based SPME.
pub fn omm_pme_reciprocal_space(
  positions : Array[Array[Double]],
  charges : Array[Double],
  params : OmmPMEParams,
) -> Double {
  let n = positions.length()
  let alpha = params.alpha
  let box = params.box
  if box[0] <= 0.0 {
    return 0.0
  }
  let mut energy = 0.0
  let k_max = 5 // reciprocal space vector cutoff
  let inv_vol = 1.0 / (box[0] * box[1] * box[2])
  let pi = 3.141592653589793
  let alpha_sq = alpha * alpha
  for kx = -k_max; kx <= k_max; kx = kx + 1 {
    for ky = -k_max; ky <= k_max; ky = ky + 1 {
      for kz = -k_max; kz <= k_max; kz = kz + 1 {
        if kx == 0 && ky == 0 && kz == 0 {
          continue
        }
        let kx_f = 2.0 * pi * kx.to_double() / box[0]
        let ky_f = 2.0 * pi * ky.to_double() / box[1]
        let kz_f = 2.0 * pi * kz.to_double() / box[2]
        let k_sq = kx_f * kx_f + ky_f * ky_f + kz_f * kz_f
        // Structure factor S(k)
        let mut s_re = 0.0
        let mut s_im = 0.0
        for i = 0; i < n; i = i + 1 {
          let phase = kx_f * positions[i][0] +
            ky_f * positions[i][1] +
            kz_f * positions[i][2]
          s_re = s_re + charges[i] * @math.cos(phase)
          s_im = s_im + charges[i] * @math.sin(phase)
        }
        let s_sq = s_re * s_re + s_im * s_im
        // Ewald factor: exp(-k^2/(4*alpha^2)) / k^2
        let factor = @math.exp(-k_sq / (4.0 * alpha_sq)) / k_sq
        energy = energy + s_sq * factor
      }
    }
  }
  energy * 2.0 * pi * inv_vol
}

///|
/// Self-interaction correction for Ewald sum.
/// E_self = -alpha/sqrt(pi) * sum(q_i^2)
pub fn omm_pme_self_energy(
  charges : Array[Double],
  params : OmmPMEParams,
) -> Double {
  let mut sum_q2 = 0.0
  for i = 0; i < charges.length(); i = i + 1 {
    sum_q2 = sum_q2 + charges[i] * charges[i]
  }
  -params.alpha / 3.141592653589793.sqrt() * sum_q2
}

///|
/// Full PME electrostatic energy.
pub fn omm_pme_energy(
  positions : Array[Array[Double]],
  charges : Array[Double],
  params : OmmPMEParams,
) -> Double {
  let e_real = omm_pme_real_space(positions, charges, params)
  let e_recip = omm_pme_reciprocal_space(positions, charges, params)
  let e_self = omm_pme_self_energy(charges, params)
  e_real + e_recip + e_self
}

///|
/// Complementary error function (erfc), Abramowitz & Stegun 7.1.26:
/// erfc(x) = t * exp(-x^2) * (a1 + a2*t + ... + a5*t^4), t = 1/(1+p*x).
fn omm_erfc(x : Double) -> Double {
  if x < 0.0 {
    return 2.0 - omm_erfc(-x)
  }
  if x > 10.0 {
    return 0.0
  }
  let t = 1.0 / (1.0 + 0.3275911 * x)
  let t2 = t * t
  let t3 = t2 * t
  let t4 = t3 * t
  let poly = 0.254829592 +
    -0.284496736 * t +
    1.421413741 * t2 +
    -1.453152027 * t3 +
    1.061405429 * t4
  poly * t * @math.exp(-x * x)
}
