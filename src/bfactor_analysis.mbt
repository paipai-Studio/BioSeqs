///|
/// B-factor analysis — thermal mobility and flexibility assessment.
///
/// Biotite provides B-factor analysis for identifying flexible/rigid regions.

///|
/// B-factor statistics.
pub struct BFactorStats {
  mean : Double
  std_dev : Double
  min_val : Double
  max_val : Double
  median : Double
  q25 : Double // 25th percentile
  q75 : Double // 75th percentile
}

///|
/// Get mean B-factor.
pub fn BFactorStats::get_mean(self : BFactorStats) -> Double {
  self.mean
}

///|
/// Get standard deviation.
pub fn BFactorStats::get_std_dev(self : BFactorStats) -> Double {
  self.std_dev
}

///|
/// Calculate B-factor statistics.
pub fn bfactor_stats(bfactors : Array[Double]) -> BFactorStats {
  let n = bfactors.length()
  if n == 0 {
    return BFactorStats::{
      mean: 0.0,
      std_dev: 0.0,
      min_val: 0.0,
      max_val: 0.0,
      median: 0.0,
      q25: 0.0,
      q75: 0.0,
    }
  }
  let mut sum = 0.0
  let mut min_v = bfactors[0]
  let mut max_v = bfactors[0]
  for b in bfactors {
    sum = sum + b
    if b < min_v {
      min_v = b
    }
    if b > max_v {
      max_v = b
    }
  }
  let mean = sum / n.to_double()
  let mut var_sum = 0.0
  for b in bfactors {
    let diff = b - mean
    var_sum = var_sum + diff * diff
  }
  let std_dev = (var_sum / n.to_double()).sqrt()
  // Simple median (sort and pick middle)
  let sorted = bfactor_sort_doubles(bfactors)
  let median = if n % 2 == 1 {
    sorted[n / 2]
  } else {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  }
  let q25 = sorted[n / 4]
  let q75 = sorted[3 * n / 4]
  BFactorStats::{
    mean,
    std_dev,
    min_val: min_v,
    max_val: max_v,
    median,
    q25,
    q75,
  }
}

///|
/// Sort an array of doubles (insertion sort).
fn bfactor_sort_doubles(arr : Array[Double]) -> Array[Double] {
  let result = Array(capacity=arr.length())
  for v in arr {
    result.push(v)
  }
  for i in 1..<result.length() {
    let key = result[i]
    let mut j = i - 1
    while j >= 0 && result[j] > key {
      result[j + 1] = result[j]
      j = j - 1
    }
    result[j + 1] = key
  }
  result
}

///|
/// Normalize B-factors to Z-scores.
///
/// Z = (B - mean) / std_dev
pub fn normalize_bfactors(bfactors : Array[Double]) -> Array[Double] {
  let stats = bfactor_stats(bfactors)
  if stats.std_dev == 0.0 {
    return Array::make(bfactors.length(), 0.0)
  }
  let result : Array[Double] = []
  for b in bfactors {
    result.push((b - stats.mean) / stats.std_dev)
  }
  result
}

///|
/// Convert B-factors to RMS displacement (Å).
///
/// RMS displacement = sqrt(B / (8 * pi^2))
pub fn bfactor_to_rmsd(bfactors : Array[Double]) -> Array[Double] {
  let factor = 8.0 * 3.141592653589793 * 3.141592653589793
  let result : Array[Double] = []
  for b in bfactors {
    result.push((b / factor).sqrt())
  }
  result
}

///|
/// Identify high B-factor residues (flexible regions).
///
/// Returns indices of residues with B-factor > mean + threshold * std_dev.
pub fn find_flexible_residues(
  bfactors : Array[Double],
  threshold~ : Double = 1.5,
) -> Array[Int32] {
  let stats = bfactor_stats(bfactors)
  let cutoff = stats.mean + threshold * stats.std_dev
  let result : Array[Int32] = []
  for i in 0..<bfactors.length() {
    if bfactors[i] > cutoff {
      result.push(i)
    }
  }
  result
}

///|
/// Identify low B-factor residues (rigid core).
pub fn find_rigid_residues(
  bfactors : Array[Double],
  threshold~ : Double = 1.0,
) -> Array[Int32] {
  let stats = bfactor_stats(bfactors)
  let cutoff = stats.mean - threshold * stats.std_dev
  let result : Array[Int32] = []
  for i in 0..<bfactors.length() {
    if bfactors[i] < cutoff {
      result.push(i)
    }
  }
  result
}

///|
/// Calculate per-residue average B-factor.
///
/// Averages B-factors across all atoms in each residue.
pub fn residue_average_bfactor(
  residue_bfactors : Array[Array[Double]],
) -> Array[Double] {
  let result : Array[Double] = []
  for atom_bfactors in residue_bfactors {
    if atom_bfactors.length() == 0 {
      result.push(0.0)
    } else {
      let mut sum = 0.0
      for b in atom_bfactors {
        sum = sum + b
      }
      result.push(sum / atom_bfactors.length().to_double())
    }
  }
  result
}
