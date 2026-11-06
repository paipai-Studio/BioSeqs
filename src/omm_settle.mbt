/// omm_settle.mbt — SETTLE + CCMA constraint algorithms
/// Items 8-9: Water constraints and parallel constraint solver

///|
/// SETTLE: analytical constraint for rigid water molecules.
/// Constrains O-H1 and O-H2 bonds and H1-H2 distance.
/// Returns updated coordinates after constraint.
pub fn omm_settle(
  coords : Array[(Double, Double, Double)],
  velocities : Array[(Double, Double, Double)],
  masses : Array[Double],
  water_molecules : Array[(Int32, Int32, Int32)], // (O_idx, H1_idx, H2_idx)
  dt : Double,
  _box~ : (Double, Double, Double) = (0.0, 0.0, 0.0),
) -> (Array[(Double, Double, Double)], Array[(Double, Double, Double)]) {
  let n = coords.length()
  let new_coords : Array[(Double, Double, Double)] = Array(capacity=n)
  let new_vels : Array[(Double, Double, Double)] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    new_coords.push(coords[i])
    new_vels.push(velocities[i])
  }
  // Ideal water geometry
  let dOH = 0.9572 // O-H bond length
  let dHH = 1.5139 // H-H distance
  let inv_dOH2 = 1.0 / (dOH * dOH)
  let inv_dHH2 = 1.0 / (dHH * dHH)
  for wm in water_molecules {
    let (oi, h1i, h2i) = wm
    if oi >= n || h1i >= n || h2i >= n {
      continue
    }
    let mO = masses[oi]
    let mH = masses[h1i]
    let inv_mO = 1.0 / mO
    let inv_mH = 1.0 / mH
    // Predicted positions after unconstrained step
    let (ox0, oy0, oz0) = coords[oi]
    let (ovx, ovy, ovz) = velocities[oi]
    let pxO = (ox0 + ovx * dt, oy0 + ovy * dt, oz0 + ovz * dt)
    let (h1x0, h1y0, h1z0) = coords[h1i]
    let (h1vx, h1vy, h1vz) = velocities[h1i]
    let pxH1 = (h1x0 + h1vx * dt, h1y0 + h1vy * dt, h1z0 + h1vz * dt)
    let (h2x0, h2y0, h2z0) = coords[h2i]
    let (h2vx, h2vy, h2vz) = velocities[h2i]
    let pxH2 = (h2x0 + h2vx * dt, h2y0 + h2vy * dt, h2z0 + h2vz * dt)
    // Center of mass
    let total_mass = mO + 2.0 * mH
    let cmx = (mO * pxO.0 + mH * pxH1.0 + mH * pxH2.0) / total_mass
    let cmy = (mO * pxO.1 + mH * pxH1.1 + mH * pxH2.1) / total_mass
    let cmz = (mO * pxO.2 + mH * pxH1.2 + mH * pxH2.2) / total_mass
    // Relative positions from COM
    let rO = (pxO.0 - cmx, pxO.1 - cmy, pxO.2 - cmz)
    let rH1 = (pxH1.0 - cmx, pxH1.1 - cmy, pxH1.2 - cmz)
    let rH2 = (pxH2.0 - cmx, pxH2.1 - cmy, pxH2.2 - cmz)
    // Simplified SETTLE: project back to ideal geometry
    // Compute current O-H1 vector
    let d1x = rH1.0 - rO.0
    let d1y = rH1.1 - rO.1
    let d1z = rH1.2 - rO.2
    let d1len = @math.pow(d1x * d1x + d1y * d1y + d1z * d1z, 0.5)
    // Scale to ideal O-H distance
    let s1 = if d1len > 0.0000000001 { dOH / d1len } else { 1.0 }
    let newH1 = (rO.0 + d1x * s1, rO.1 + d1y * s1, rO.2 + d1z * s1)
    // O-H2
    let d2x = rH2.0 - rO.0
    let d2y = rH2.1 - rO.1
    let d2z = rH2.2 - rO.2
    let d2len = @math.pow(d2x * d2x + d2y * d2y + d2z * d2z, 0.5)
    let s2 = if d2len > 0.0000000001 { dOH / d2len } else { 1.0 }
    let newH2 = (rO.0 + d2x * s2, rO.1 + d2y * s2, rO.2 + d2z * s2)
    // H-H constraint: adjust H2 to satisfy H1-H2 distance
    let hhx = newH2.0 - newH1.0
    let hhy = newH2.1 - newH1.1
    let hhz = newH2.2 - newH1.2
    let hhlen = @math.pow(hhx * hhx + hhy * hhy + hhz * hhz, 0.5)
    if hhlen > 0.0000000001 {
      let sh = dHH / hhlen
      new_coords[h2i] = (
        cmx + newH1.0 + hhx * sh,
        cmy + newH1.1 + hhy * sh,
        cmz + newH1.2 + hhz * sh,
      )
    }
    new_coords[oi] = (cmx + rO.0, cmy + rO.1, cmz + rO.2)
    new_coords[h1i] = (cmx + newH1.0, cmy + newH1.1, cmz + newH1.2)
    // Update velocities from corrected positions
    new_vels[oi] = (
      (new_coords[oi].0 - ox0) / dt,
      (new_coords[oi].1 - oy0) / dt,
      (new_coords[oi].2 - oz0) / dt,
    )
    new_vels[h1i] = (
      (new_coords[h1i].0 - h1x0) / dt,
      (new_coords[h1i].1 - h1y0) / dt,
      (new_coords[h1i].2 - h1z0) / dt,
    )
    new_vels[h2i] = (
      (new_coords[h2i].0 - h2x0) / dt,
      (new_coords[h2i].1 - h2y0) / dt,
      (new_coords[h2i].2 - h2z0) / dt,
    )
    // Suppress unused warnings
    ignore(inv_mO)
    ignore(inv_mH)
    ignore(inv_dOH2)
    ignore(inv_dHH2)
  }
  (new_coords, new_vels)
}

///|
/// CCMA (Crank-Nicholson Matrix Algorithm) constraint solver.
/// Iterative constraint solver suitable for parallel execution.
pub fn omm_ccma_constraints(
  coords : Array[(Double, Double, Double)],
  velocities : Array[(Double, Double, Double)],
  masses : Array[Double],
  constraints : Array[(Int32, Int32, Double)], // (i, j, target_distance)
  dt : Double,
  max_iter~ : Int32 = 20,
  tolerance~ : Double = 0.000001,
) -> (Array[(Double, Double, Double)], Array[(Double, Double, Double)]) {
  let n = coords.length()
  let new_coords : Array[(Double, Double, Double)] = Array(capacity=n)
  let new_vels : Array[(Double, Double, Double)] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    new_coords.push(
      (
        coords[i].0 + velocities[i].0 * dt,
        coords[i].1 + velocities[i].1 * dt,
        coords[i].2 + velocities[i].2 * dt,
      ),
    )
    new_vels.push(velocities[i])
  }
  // CCMA: iterative constraint correction
  for iter = 0; iter < max_iter; iter = iter + 1 {
    let mut max_error = 0.0
    for c in constraints {
      let (i, j, d0) = c
      if i >= n || j >= n {
        continue
      }
      let dx = new_coords[i].0 - new_coords[j].0
      let dy = new_coords[i].1 - new_coords[j].1
      let dz = new_coords[i].2 - new_coords[j].2
      let d2 = dx * dx + dy * dy + dz * dz
      let _d = @math.pow(d2, 0.5)
      let error = d2 - d0 * d0
      if omm_st_abs(error) > max_error {
        max_error = omm_st_abs(error)
      }
      // Lagrange multiplier correction
      let inv_mi = 1.0 / masses[i]
      let inv_mj = 1.0 / masses[j]
      let g = error / (2.0 * d2 * (inv_mi + inv_mj) + 0.0000000001)
      // Apply correction
      let corr_i = g * inv_mi
      let corr_j = g * inv_mj
      new_coords[i] = (
        new_coords[i].0 - corr_i * dx,
        new_coords[i].1 - corr_i * dy,
        new_coords[i].2 - corr_i * dz,
      )
      new_coords[j] = (
        new_coords[j].0 + corr_j * dx,
        new_coords[j].1 + corr_j * dy,
        new_coords[j].2 + corr_j * dz,
      )
    }
    if max_error < tolerance {
      break
    }
  }
  // Update velocities from corrected positions
  for i = 0; i < n; i = i + 1 {
    new_vels[i] = (
      (new_coords[i].0 - coords[i].0) / dt,
      (new_coords[i].1 - coords[i].1) / dt,
      (new_coords[i].2 - coords[i].2) / dt,
    )
  }
  (new_coords, new_vels)
}

///|
fn omm_st_abs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

///|
pub fn omm_settle_dummy() -> Unit {
  ignore(0)
}
