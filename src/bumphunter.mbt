///|
/// Bumphunter - Methylation/signal region (bump) detection.
///
/// References:
///   - Bioconductor bumphunter package
///   - Aryee MJ, Jaffe AE, Corrada-Bravo H, et al. (2014) Minfi: a
///     flexible and comprehensive Bioconductor package for the analysis
///     of Infinium DNA methylation microarrays. Bioinformatics 30(10):1363-9.
///   - Jaffe AE, Murakami P, Lee H, et al. (2012) Bump hunting to
///     identify differentially methylated regions in epigenetic
///     epidemiology studies. Int J Epidemiol 41(1):200-9.
///
/// The bumphunter algorithm identifies genomic regions ("bumps") where
/// a smoothed signal (e.g. methylation difference between cases and
/// controls) exceeds a threshold, then assesses significance via
/// label permutations.
///
/// All public identifiers are prefixed with `bump_` to avoid conflicts.

// ============================================================================
// Data structures
// ============================================================================

///|
/// A single genomic position with a measured value for a set of samples.
pub struct BumpPosition {
  chrom : String
  pos : Int
  values : Array[Double] // per-sample value (e.g. beta value)
}

///|
/// Create a genomic position.
pub fn BumpPosition::new(
  chrom : String,
  pos : Int,
  values : Array[Double],
) -> BumpPosition {
  { chrom, pos, values }
}

///|
/// Get chromosome.
pub fn BumpPosition::chrom(self : BumpPosition) -> String {
  self.chrom
}

///|
/// Get position.
pub fn BumpPosition::pos(self : BumpPosition) -> Int {
  self.pos
}

///|
/// Get values.
pub fn BumpPosition::values(self : BumpPosition) -> Array[Double] {
  self.values
}

///|
/// A detected bump: a contiguous region where the smoothed signal
/// exceeds the cutoff threshold.
pub struct BumpRegion {
  chrom : String
  start : Int
  end_ : Int
  value : Double // max |smoothed| within the bump
  area : Double // sum of smoothed values within the bump
  index_start : Int // index into the position array
  index_end : Int // index into the position array (inclusive)
}

///|
/// Create a bump region.
pub fn BumpRegion::new(
  chrom : String,
  start : Int,
  end_ : Int,
  value : Double,
  area : Double,
  index_start : Int,
  index_end : Int,
) -> BumpRegion {
  { chrom, start, end_, value, area, index_start, index_end }
}

///|
/// Get chromosome.
pub fn BumpRegion::chrom(self : BumpRegion) -> String {
  self.chrom
}

///|
/// Get start position.
pub fn BumpRegion::start(self : BumpRegion) -> Int {
  self.start
}

///|
/// Get end position.
pub fn BumpRegion::end_(self : BumpRegion) -> Int {
  self.end_
}

///|
/// Get peak value (max |smoothed|).
pub fn BumpRegion::value(self : BumpRegion) -> Double {
  self.value
}

///|
/// Get area (sum of smoothed values).
pub fn BumpRegion::area(self : BumpRegion) -> Double {
  self.area
}

///|
/// Get index start.
pub fn BumpRegion::index_start(self : BumpRegion) -> Int {
  self.index_start
}

///|
/// Get index end.
pub fn BumpRegion::index_end(self : BumpRegion) -> Int {
  self.index_end
}

///|
/// Length of the bump in genomic coordinates.
pub fn BumpRegion::length(self : BumpRegion) -> Int {
  self.end_ - self.start
}

///|
/// A significant bump with permutation-based statistics.
pub struct BumpResult {
  chrom : String
  start : Int
  end_ : Int
  value : Double
  area : Double
  p_value : Double // permutation-based p-value
  fdr : Double // false discovery rate
  index_start : Int
  index_end : Int
}

///|
/// Create a bump result.
pub fn BumpResult::new(
  chrom : String,
  start : Int,
  end_ : Int,
  value : Double,
  area : Double,
  p_value : Double,
  fdr : Double,
  index_start : Int,
  index_end : Int,
) -> BumpResult {
  {
    chrom,
    start,
    end_,
    value,
    area,
    p_value,
    fdr,
    index_start,
    index_end,
  }
}

///|
/// Get chromosome.
pub fn BumpResult::chrom(self : BumpResult) -> String {
  self.chrom
}

///|
/// Get start position.
pub fn BumpResult::start(self : BumpResult) -> Int {
  self.start
}

///|
/// Get end position.
pub fn BumpResult::end_(self : BumpResult) -> Int {
  self.end_
}

///|
/// Get peak value.
pub fn BumpResult::value(self : BumpResult) -> Double {
  self.value
}

///|
/// Get area.
pub fn BumpResult::area(self : BumpResult) -> Double {
  self.area
}

///|
/// Get p-value.
pub fn BumpResult::p_value(self : BumpResult) -> Double {
  self.p_value
}

///|
/// Get FDR.
pub fn BumpResult::fdr(self : BumpResult) -> Double {
  self.fdr
}

///|
/// Get index start.
pub fn BumpResult::index_start(self : BumpResult) -> Int {
  self.index_start
}

///|
/// Get index end.
pub fn BumpResult::index_end(self : BumpResult) -> Int {
  self.index_end
}

// ============================================================================
// Core algorithm
// ============================================================================

///|
/// Compute per-position t-statistic comparing two groups of samples.
/// Positive values mean group1 > group2.
/// Returns an array of t-statistics, one per position.
pub fn bump_compute_t_statistics_test(
  positions : Array[BumpPosition],
  group1_size : Int,
) -> Array[Double] {
  let n = positions.length()
  let n1 = group1_size
  let n2 = n1 // assume balanced; total = n1 + n2
  let stats : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    let vals = positions[i].values
    let g1 : Array[Double] = []
    let g2 : Array[Double] = []
    for j in 0..<vals.length() {
      if j < n1 {
        g1.push(vals[j])
      } else {
        g2.push(vals[j])
      }
    }
    let m1 = stat_mean(g1)
    let m2 = stat_mean(g2)
    let v1 = stat_variance(g1)
    let v2 = stat_variance(g2)
    // Add small epsilon to avoid division by zero when variance is 0
    let denom = (v1 / n1.to_double() + v2 / n2.to_double() + 1.0e-10).sqrt()
    stats[i] = (m1 - m2) / denom
  }
  stats
}

///|
/// Smooth a signal using a running mean with the given half-window size.
/// The smoothed value at position i is the mean of values in [i-k, i+k].
pub fn bump_smooth_test(stats : Array[Double], k : Int) -> Array[Double] {
  let n = stats.length()
  let smoothed : Array[Double] = Array::make(n, 0.0)
  if n == 0 {
    return smoothed
  }
  for i in 0..<n {
    let lo = if i - k < 0 { 0 } else { i - k }
    let hi = if i + k >= n { n - 1 } else { i + k }
    let mut sum = 0.0
    let mut cnt = 0
    for j in lo..=hi {
      sum = sum + stats[j]
      cnt = cnt + 1
    }
    smoothed[i] = sum / cnt.to_double()
  }
  smoothed
}

///|
/// Find candidate bumps: maximal contiguous runs where |smoothed| > cutoff.
/// Bumps cannot span chromosomes.
pub fn bump_find_candidates_test(
  positions : Array[BumpPosition],
  smoothed : Array[Double],
  cutoff : Double,
) -> Array[BumpRegion] {
  let n = positions.length()
  let bumps : Array[BumpRegion] = []
  if n == 0 {
    return bumps
  }
  let mut i = 0
  while i < n {
    if bump_abs(smoothed[i]) > cutoff {
      let chrom = positions[i].chrom
      let start_idx = i
      let mut max_val = bump_abs(smoothed[i])
      let mut area = 0.0
      let mut signed_max = smoothed[i]
      while i < n &&
        bump_abs(smoothed[i]) > cutoff &&
        positions[i].chrom == chrom {
        area = area + smoothed[i]
        if bump_abs(smoothed[i]) > max_val {
          max_val = bump_abs(smoothed[i])
          signed_max = smoothed[i]
        }
        i = i + 1
      }
      let end_idx = i - 1
      bumps.push(
        BumpRegion::new(
          chrom,
          positions[start_idx].pos,
          positions[end_idx].pos,
          signed_max,
          area,
          start_idx,
          end_idx,
        ),
      )
    } else {
      i = i + 1
    }
  }
  bumps
}

///|
/// Absolute value helper.
fn bump_abs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

///|
/// Permutation test: shuffle the sample labels and recompute bumps.
/// Returns the number of permutation bumps with value >= threshold
/// (used for p-value estimation).
fn bump_permutation_count(
  positions : Array[BumpPosition],
  group1_size : Int,
  k : Int,
  cutoff : Double,
  threshold : Double,
  seed : Int,
) -> Int {
  let n = positions.length()
  if n == 0 {
    return 0
  }
  // Shuffle the values within each position
  let shuffled : Array[BumpPosition] = []
  let mut rng_state = seed
  for i in 0..<n {
    let vals = positions[i].values
    let n_vals = vals.length()
    let n_vals_minus_1 = n_vals - 1
    // Fisher-Yates shuffle with simple LCG
    let shuffled_vals : Array[Double] = []
    for j in 0..<n_vals {
      shuffled_vals.push(vals[j])
    }
    for j in 0..<n_vals_minus_1 {
      rng_state = (rng_state * 1103515245 + 12345) & 0x7FFFFFFF
      let r = (rng_state % (n_vals - j)).to_double()
      let swap_idx = j + r.to_int()
      let tmp = shuffled_vals[j]
      shuffled_vals[j] = shuffled_vals[swap_idx]
      shuffled_vals[swap_idx] = tmp
    }
    shuffled.push(BumpPosition::new(positions[i].chrom, positions[i].pos, shuffled_vals))
  }
  let perm_stats = bump_compute_t_statistics_test(shuffled, group1_size)
  let perm_smoothed = bump_smooth_test(perm_stats, k)
  let perm_bumps = bump_find_candidates_test(positions, perm_smoothed, cutoff)
  // Count bumps with value >= threshold
  let mut count = 0
  for b in perm_bumps {
    if bump_abs(b.value) >= threshold {
      count = count + 1
    }
  }
  count
}

///|
/// Main bumphunter analysis.
///
/// Parameters:
/// - positions: genomic positions with per-sample values
/// - group1_size: number of samples in group 1 (first group1_size samples)
/// - cutoff: threshold for |smoothed| to define a candidate bump
/// - smooth_k: half-window size for smoothing
/// - n_permutations: number of permutations for significance testing
///
/// Returns significant bumps with p-values and FDR.
pub fn bump_hunt(
  positions : Array[BumpPosition],
  group1_size : Int,
  cutoff : Double,
  smooth_k : Int,
  n_permutations : Int,
) -> Array[BumpResult] {
  let n = positions.length()
  if n == 0 || group1_size <= 0 {
    return []
  }
  // 1. Compute per-position t-statistics
  let stats = bump_compute_t_statistics_test(positions, group1_size)
  // 2. Smooth
  let smoothed = bump_smooth_test(stats, smooth_k)
  // 3. Find candidate bumps
  let candidates = bump_find_candidates_test(positions, smoothed, cutoff)
  if candidates.length() == 0 {
    return []
  }
  // 4. Permutation testing for p-values
  let total_perm_bumps : Array[Int] = Array::make(candidates.length(), 0)
  let mut total_perm_count = 0
  for p in 0..<n_permutations {
    let seed = 42 + p * 1009
    let perm_count = bump_permutation_count(
      positions,
      group1_size,
      smooth_k,
      cutoff,
      0.0,
      seed,
    )
    total_perm_count = total_perm_count + perm_count
    for i in 0..<candidates.length() {
      let perm_bumps_for_threshold = bump_permutation_count(
        positions,
        group1_size,
        smooth_k,
        cutoff,
        bump_abs(candidates[i].value),
        seed,
      )
      total_perm_bumps[i] = total_perm_bumps[i] + perm_bumps_for_threshold
    }
  }
  // 5. Compute p-values and FDR
  let n_perm = n_permutations.to_double()
  let p_values : Array[Double] = Array::make(candidates.length(), 1.0)
  let fdrs : Array[Double] = Array::make(candidates.length(), 1.0)
  for i in 0..<candidates.length() {
    // p-value: fraction of permutations producing a bump >= this value
    p_values[i] = (total_perm_bumps[i].to_double() + 1.0) / (n_perm + 1.0)
    // FDR: expected false positives / observed positives
    let expected_fp = total_perm_count.to_double() / n_perm
    let observed = candidates.length().to_double()
    if observed > 0.0 {
      fdrs[i] = expected_fp / observed
    } else {
      fdrs[i] = 1.0
    }
    if fdrs[i] > 1.0 {
      fdrs[i] = 1.0
    }
  }
  // 6. Build results
  let results : Array[BumpResult] = []
  for i in 0..<candidates.length() {
    let c = candidates[i]
    results.push(
      BumpResult::new(
        c.chrom,
        c.start,
        c.end_,
        c.value,
        c.area,
        p_values[i],
        fdrs[i],
        c.index_start,
        c.index_end,
      ),
    )
  }
  results
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Generate sample methylation data with a known bump on chr1.
/// 8 samples (4 cases + 4 controls), 30 CpG sites on chr1.
/// Cases have elevated methylation at positions 10-20 (the bump).
pub fn bump_sample_data() -> Array[BumpPosition] {
  let positions : Array[BumpPosition] = []
  for i in 0..<30 {
    let chrom = "chr1"
    let pos = 1000 + i * 100
    let values : Array[Double] = []
    // 4 cases (group 1) + 4 controls (group 2)
    for s in 0..<8 {
      let base = 0.5
      let mut val = base
      // Add the bump: elevated methylation in cases at positions 10-20
      if s < 4 && i >= 10 && i <= 20 {
        val = base + 0.3
      }
      // Add per-sample noise (deterministic, varies by sample to ensure
      // non-zero within-group variance)
      val = val + ((i + s) % 3).to_double() * 0.01
      values.push(val)
    }
    positions.push(BumpPosition::new(chrom, pos, values))
  }
  positions
}

///|
/// Format a BumpResult as a readable string.
pub fn BumpResult::to_string(self : BumpResult) -> String {
  let sign = if self.value >= 0.0 { "+" } else { "-" }
  "chr=" + self.chrom + " start=" + self.start.to_string() + " end=" + self.end_.to_string() + " value=" + sign + bump_abs(self.value).to_string() + " p=" + self.p_value.to_string() + " fdr=" + self.fdr.to_string()
}

///|
/// Format all results as a table.
pub fn bump_results_to_string(results : Array[BumpResult]) -> String {
  let mut s = "Bumphunter Results (" + results.length().to_string() + " bumps)\n"
  s = s + "chrom\tstart\tend\tvalue\tp_value\tfdr\n"
  for r in results {
    s = s + r.chrom + "\t" + r.start.to_string() + "\t" + r.end_.to_string() + "\t" + r.value.to_string() + "\t" + r.p_value.to_string() + "\t" + r.fdr.to_string() + "\n"
  }
  s
}
