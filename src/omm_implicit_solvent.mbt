/// omm_implicit_solvent.mbt — GBn/GBn2 + OBC2 implicit solvent models
/// Items 6-7: Advanced generalized Born implicit solvent

///|
/// GBn parameters (GBneck model)
pub struct OmmGBnParams {
  // Solute radii (Born radii offset parameters)
  radii : Array[Double]
  // Charges
  charges : Array[Double]
  // GBn-specific: neck scale factor and neck cut
  neck_scale : Double
  neck_cut : Double
  // SA surface tension (kcal/mol/Å²)
  surface_tension : Double
  // Solvent dielectric
  solute_dielectric : Double
  solvent_dielectric : Double
}

///|
/// Compute OBC2 (Onufriev-Bashford-Case 2) Born radii.
/// OBC2 uses a modified Still equation with empirical scaling.
pub fn omm_obc2_born_radii(
  coords : Array[(Double, Double, Double)],
  radii : Array[Double],
) -> Array[Double] {
  let n = coords.length()
  let born_radii : Array[Double] = Array(capacity=n)
  // OBC2 parameters
  let obc_alpha = 1.0
  let obc_beta = 0.8
  let obc_gamma = 4.85
  for i = 0; i < n; i = i + 1 {
    let ri = radii[i]
    // Compute offset and radius
    let mut offset = 0.0
    let mut radius = ri
    for j = 0; j < n; j = j + 1 {
      if i == j {
        continue
      }
      let dx = coords[i].0 - coords[j].0
      let dy = coords[i].1 - coords[j].1
      let dz = coords[i].2 - coords[j].2
      let r = @math.pow(dx * dx + dy * dy + dz * dz, 0.5)
      let rj = radii[j]
      let sj = rj - 0.09
      let sk = ri + 0.14
      if r > sk + sj {
        // Case 1: atoms far apart
        let th = 1.0 / (r - sj)
        offset = offset + th
        let th2 = th * th
        radius = radius + rj * rj * th2 / r
      } else if r > sk - sj && r > 0.0000000001 {
        // Case 2: overlapping region
        let log_arg = (r + sj - sk) / (r - sj + sk)
        if log_arg > 0.0 {
          let lr = @math.ln(log_arg)
          offset = offset + 2.0 * lr / (r + sj - sk)
          radius = radius + 2.0 * rj * rj * lr / ((r + sj - sk) * (r - sj + sk))
        }
      } else if r > 0.0000000001 {
        // Case 3: r < sk - sj
        let log_arg = (sk + sj - r) / (sk - sj + r)
        if log_arg > 0.0 {
          let lr = @math.ln(log_arg)
          offset = offset +
            (2.0 * lr + 4.0 * sj / (r + sk - sj)) / (r + sj - sk)
          radius = radius + 2.0 * rj * rj * lr / ((r + sj - sk) * (r - sj + sk))
        }
      }
    }
    // OBC2 empirical correction
    radius = ri * ri / (ri - offset * ri * ri)
    let mut reff = 1.0 /
      (
        1.0 / ri -
        obc_alpha * offset +
        obc_beta * radius +
        obc_gamma * radius * radius
      )
    if reff < 0.0 {
      reff = 0.5 * ri
    }
    born_radii.push(reff)
  }
  born_radii
}

///|
/// Compute GBn (GBneck) Born radii with neck correction.
pub fn omm_gbn_born_radii(
  coords : Array[(Double, Double, Double)],
  radii : Array[Double],
  params : OmmGBnParams,
) -> Array[Double] {
  let n = coords.length()
  // Start with OBC2 radii
  let born_radii = omm_obc2_born_radii(coords, radii)
  // Apply neck correction
  for i = 0; i < n; i = i + 1 {
    let mut neck_sum = 0.0
    for j = 0; j < n; j = j + 1 {
      if i == j {
        continue
      }
      let dx = coords[i].0 - coords[j].0
      let dy = coords[i].1 - coords[j].1
      let dz = coords[i].2 - coords[j].2
      let r = @math.pow(dx * dx + dy * dy + dz * dz, 0.5)
      if r < params.neck_cut {
        // Neck correction term (simplified)
        let ri = radii[i]
        let rj = radii[j]
        let u = (r - ri - rj) / params.neck_cut
        if u > -1.0 && u < 1.0 {
          // Smooth switching function
          let sw = 0.5 * (1.0 + @math.cos(3.141592653589793 * (u + 1.0) / 2.0))
          neck_sum = neck_sum +
            params.neck_scale * sw * ri * rj / (r * r + 0.0001)
        }
      }
    }
    born_radii[i] = born_radii[i] + neck_sum
  }
  born_radii
}

///|
/// Compute GBn electrostatic energy.
/// E_GB = -0.5 * (1/eps_in - 1/eps_out) * sum_i sum_j (q_i*q_j / f_ij)
/// where f_ij = sqrt(r_ij^2 + B_i*B_j*exp(-r_ij^2/(4*B_i*B_j)))
pub fn omm_gbn_energy(
  params : OmmGBnParams,
  coords : Array[(Double, Double, Double)],
  born_radii : Array[Double],
) -> Double {
  let n = coords.length()
  let prefactor = -0.5 *
    (1.0 / params.solute_dielectric - 1.0 / params.solvent_dielectric) *
    138.935485
  let mut energy = 0.0
  for i = 0; i < n; i = i + 1 {
    // Self energy
    energy = energy +
      prefactor * params.charges[i] * params.charges[i] / born_radii[i]
    // Pair interactions
    for j = i + 1; j < n; j = j + 1 {
      let dx = coords[i].0 - coords[j].0
      let dy = coords[i].1 - coords[j].1
      let dz = coords[i].2 - coords[j].2
      let r2 = dx * dx + dy * dy + dz * dz
      let bi = born_radii[i]
      let bj = born_radii[j]
      let f_ij = @math.pow(
        r2 + bi * bj * @math.exp(-r2 / (4.0 * bi * bj + 0.0000000001)),
        0.5,
      )
      energy = energy +
        2.0 * prefactor * params.charges[i] * params.charges[j] / f_ij
    }
  }
  energy
}

///|
/// Compute GBn2 energy (GBn + surface area term).
pub fn omm_gbn2_energy(
  params : OmmGBnParams,
  coords : Array[(Double, Double, Double)],
  born_radii : Array[Double],
) -> Double {
  let gb_energy = omm_gbn_energy(params, coords, born_radii)
  // Surface area correction
  let mut sa_energy = 0.0
  for i = 0; i < coords.length(); i = i + 1 {
    let ri = params.radii[i]
    let bi = born_radii[i]
    // SA term: surface_tension * 4*pi*ri^2 * (ri/bi)^6
    let ratio = ri / bi
    sa_energy = sa_energy +
      params.surface_tension *
      4.0 *
      3.141592653589793 *
      ri *
      ri *
      ratio *
      ratio *
      ratio *
      ratio *
      ratio *
      ratio
  }
  gb_energy + sa_energy
}

///|
/// Compute OBC2 total energy (convenience wrapper).
pub fn omm_obc2_energy(
  charges : Array[Double],
  radii : Array[Double],
  coords : Array[(Double, Double, Double)],
  solute_dielectric~ : Double = 1.0,
  solvent_dielectric~ : Double = 78.5,
) -> Double {
  let born_radii = omm_obc2_born_radii(coords, radii)
  let params = OmmGBnParams::{
    radii,
    charges,
    neck_scale: 0.0,
    neck_cut: 0.0,
    surface_tension: 0.0,
    solute_dielectric,
    solvent_dielectric,
  }
  omm_gbn_energy(params, coords, born_radii)
}

///|
pub fn omm_implicit_solvent_dummy() -> Unit {
  ignore(OmmGBnParams::{
    radii: Array(capacity=0),
    charges: Array(capacity=0),
    neck_scale: 0.0,
    neck_cut: 0.0,
    surface_tension: 0.0,
    solute_dielectric: 1.0,
    solvent_dielectric: 78.5,
  })
}
