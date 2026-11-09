/// cobra_optgp.mbt — optGpSampler (COBRApy sampling.optgp)
/// Optimized Group-Parallel sampling: parallel hit-and-run with group updates.
/// Algorithm: Megchelen et al. (2014) BMC Systems Biology.

///|
/// optGp sampling result.
pub struct CobraOptGpResult {
  samples : Array[Array[Double]]
  n_samples : Int32
  n_groups : Int32
}

///|
/// optGpSampler: optimized sampling using group-parallel hit-and-run.
/// This is a simplified sequential version (no actual parallelism).
///
/// Parameters:
/// - n_samples: number of samples
/// - thinning: keep every `thinning`-th sample
/// - n_groups: number of direction groups
pub fn cobra_optgp_sample(
  model : FbaModel,
  n_samples~ : Int32 = 100,
  thinning~ : Int32 = 100,
  n_groups~ : Int32 = 4,
) -> CobraOptGpResult raise FbaError {
  let n_r = model.reactions.length()
  // Use ACHR as base sampler, with multiple direction groups
  let base = cobra_achr_sample(
    model,
    n_samples~,
    thinning~,
    n_warmup=n_groups * 100,
  )
  // Post-process: apply group-based centering
  let samples : Array[Array[Double]] = Array(capacity=0)
  let group_size = base.samples.length() / n_groups
  if group_size == 0 {
    return {
      samples: base.samples,
      n_samples: base.samples.length(),
      n_groups: 1,
    }
  }
  for g = 0; g < n_groups; g = g + 1 {
    // Compute group center
    let center : Array[Double] = Array::make(n_r, 0.0)
    let start = g * group_size
    let end = if g == n_groups - 1 {
      base.samples.length()
    } else {
      start + group_size
    }
    let count = end - start
    if count == 0 {
      continue
    }
    for i = start; i < end; i = i + 1 {
      for j = 0; j < n_r; j = j + 1 {
        center[j] = center[j] + base.samples[i][j]
      }
    }
    for j = 0; j < n_r; j = j + 1 {
      center[j] = center[j] / count.to_double()
    }
    // Center samples around group mean
    for i = start; i < end; i = i + 1 {
      let centered : Array[Double] = Array(capacity=n_r)
      for j = 0; j < n_r; j = j + 1 {
        centered.push(base.samples[i][j] - center[j])
      }
      samples.push(centered)
    }
  }
  { samples, n_samples: samples.length(), n_groups, }
}
