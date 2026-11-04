/// cobra_achr.mbt — ACHR Sampler (COBRApy sampling.achr)
/// Artificial Center Hit-and-Run sampling of the flux space.
/// Algorithm: Kaufman & Smith (1998); Schellenberger & Palsson (2012).

///|
/// ACHR sampling result.
pub struct CobraAchrResult {
  samples : Array[Array[Double]]
  n_samples : Int32
  n_warmup : Int32
}

///|
/// ACHR sampling: generate samples from the feasible flux space.
///
/// Parameters:
/// - n_samples: number of samples to generate
/// - thinning: keep every `thinning`-th sample
/// - n_warmup: warmup steps
pub fn cobra_achr_sample(
  model : FbaModel,
  n_samples~ : Int32 = 100,
  thinning~ : Int32 = 100,
  n_warmup~ : Int32 = 1000,
) -> CobraAchrResult raise FbaError {
  let n_r = model.reactions.length()
  let _m = model.metabolites.length()
  // Get flux bounds
  let lb : Array[Double] = Array(capacity=n_r)
  let ub : Array[Double] = Array(capacity=n_r)
  for j = 0; j < n_r; j = j + 1 {
    lb.push(model.reactions[j].lower_bound)
    ub.push(model.reactions[j].upper_bound)
  }
  // Find a center point (average of min/max from FVA)
  let center : Array[Double] = Array(capacity=n_r)
  for j = 0; j < n_r; j = j + 1 {
    center.push((lb[j] + ub[j]) / 2.0)
  }
  // Warmup: generate initial directions by solving FBA with random objectives
  let directions : Array[Array[Double]] = Array(capacity=0)
  let mut seed = 42
  for _w = 0; _w < n_warmup; _w = _w + 1 {
    // Random objective
    let reactions : Array[FbaReaction] = Array(capacity=0)
    for j = 0; j < n_r; j = j + 1 {
      let r = model.reactions[j]
      seed = (seed * 1103515245 + 12345) % 2147483647
      let rand_obj = (seed.to_double() / 2147483647.0 - 0.5) * 2.0
      reactions.push({
        id: r.id,
        name: r.name,
        lower_bound: r.lower_bound,
        upper_bound: r.upper_bound,
        objective: rand_obj,
        genes: r.genes,
      })
    }
    let rand_model = {
      metabolites: model.metabolites,
      reactions,
      stoichiometry: model.stoichiometry,
    }
    let sol = fba_solve(rand_model)
    if sol.status == "optimal" {
      directions.push(sol.fluxes)
    }
  }
  if directions.length() == 0 {
    // Fallback: use center
    directions.push(center)
  }
  // Hit-and-run sampling
  let samples : Array[Array[Double]] = Array(capacity=0)
  let mut current = center
  let total_steps = n_samples * thinning
  for step = 0; step < total_steps; step = step + 1 {
    // Pick a random direction from warmup
    seed = (seed * 1103515245 + 12345) % 2147483647
    let dir_idx = seed % directions.length()
    let dir = directions[dir_idx.abs()]
    // Compute step size: find max alpha such that current + alpha * dir stays in bounds
    let mut alpha_max = 1.0
    let mut alpha_min = -1.0
    for j = 0; j < n_r; j = j + 1 {
      if dir[j].abs() > 1.0e-9 {
        let hi = (ub[j] - current[j]) / dir[j]
        let lo = (lb[j] - current[j]) / dir[j]
        if hi > lo {
          if hi < alpha_max {
            alpha_max = hi
          }
          if lo > alpha_min {
            alpha_min = lo
          }
        } else {
          if lo < alpha_max {
            alpha_max = lo
          }
          if hi > alpha_min {
            alpha_min = hi
          }
        }
      }
    }
    // Random alpha in [alpha_min, alpha_max]
    seed = (seed * 1103515245 + 12345) % 2147483647
    let t = seed.to_double() / 2147483647.0
    let alpha = alpha_min + t * (alpha_max - alpha_min)
    // Update current
    let next : Array[Double] = Array(capacity=n_r)
    for j = 0; j < n_r; j = j + 1 {
      next.push(current[j] + alpha * dir[j])
    }
    current = next
    // Store sample
    if step % thinning == 0 {
      let sample : Array[Double] = Array(capacity=n_r)
      for j = 0; j < n_r; j = j + 1 {
        sample.push(current[j])
      }
      samples.push(sample)
    }
  }
  { samples, n_samples: samples.length(), n_warmup: directions.length(), }
}
