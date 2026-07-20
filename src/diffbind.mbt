///|
/// DiffBind - Differential Binding analysis for ChIP-seq data.
/// Inspired by Bioconductor DiffBind package.
///
/// Provides functions for peak consistency analysis, differential binding
/// site detection, and visualization data generation.

///|
/// Peak file format enumeration.
pub enum PeakFormat {
  BedFormat
  NarrowPeak
  BroadPeak
  MacsXls
} derive(Eq)

///|
/// Score type enumeration for peak scoring.
pub enum ScoreType {
  DBA_SCORE_TISSUE_MINUS
  DBA_SCORE_READS
  DBA_SCORE_RPKM
  DBA_SCORE_NORMALIZED
} derive(Eq)

///|
/// Genomic peak region.
pub struct DBAPeak {
  chr : String
  start : Int
  end : Int
  score : Double
  summit : Int
}

///|
/// Sample metadata and peaks.
pub struct DBASample {
  sample_id : String
  condition : String
  replicate : Int
  peak_file : String
  bam_file : String
  peaks : Array[DBAPeak]
}

///|
/// Configuration object for DiffBind analysis.
pub struct DBAConfig {
  samples : Array[DBASample]
  peaks : Array[Array[DBAPeak]]
  peak_format : PeakFormat
  score : ScoreType
  conditions : Array[String]
  replicates : Array[Int]
}

///|
/// Result of differential binding analysis.
pub struct DBAResult {
  peaks : Array[DBAPeak]
  pvalues : Array[Double]
  fdrs : Array[Double]
  fold_changes : Array[Double]
  contrast : (String, String)
  n_significant : Int
}

///|
/// Result of overlap analysis.
pub struct OverlapResult {
  consensus_peaks : Array[DBAPeak]
  overlap_matrix : Array[Array[Int]]
  sample_names : Array[String]
  min_overlap : Double
}

///|
/// Construct a new DBAPeak with default score and summit at midpoint.
pub fn DBAPeak::new(chr : String, start : Int, end : Int) -> DBAPeak {
  DBAPeak::{ chr, start, end, score: 0.0, summit: (start + end) / 2 }
}

///|
/// Construct a new DBASample with empty peaks and file paths.
pub fn DBASample::new(
  sample_id : String,
  condition : String,
  replicate : Int,
) -> DBASample {
  DBASample::{
    sample_id,
    condition,
    replicate,
    peak_file: "",
    bam_file: "",
    peaks: Array::new(),
  }
}

///|
/// Construct a new DBAConfig from an array of samples.
pub fn DBAConfig::new(samples : Array[DBASample]) -> DBAConfig {
  let conditions : Array[String] = Array::new()
  let replicates : Array[Int] = Array::new()
  let peaks : Array[Array[DBAPeak]] = Array::new()
  let mut i = 0
  while i < samples.length() {
    conditions.push(samples[i].condition)
    replicates.push(samples[i].replicate)
    peaks.push(samples[i].peaks)
    i = i + 1
  }
  DBAConfig::{
    samples,
    peaks,
    peak_format: PeakFormat::BedFormat,
    score: ScoreType::DBA_SCORE_READS,
    conditions,
    replicates,
  }
}

///|
/// String representation of PeakFormat.
pub fn PeakFormat::to_string(self : PeakFormat) -> String {
  match self {
    PeakFormat::BedFormat => "BED"
    PeakFormat::NarrowPeak => "NarrowPeak"
    PeakFormat::BroadPeak => "BroadPeak"
    PeakFormat::MacsXls => "MACS xls"
  }
}

///|
/// String representation of ScoreType.
pub fn ScoreType::to_string(self : ScoreType) -> String {
  match self {
    ScoreType::DBA_SCORE_TISSUE_MINUS => "TissueMinus"
    ScoreType::DBA_SCORE_READS => "Reads"
    ScoreType::DBA_SCORE_RPKM => "RPKM"
    ScoreType::DBA_SCORE_NORMALIZED => "Normalized"
  }
}

///|
/// String representation of DBAPeak.
pub fn DBAPeak::to_string(self : DBAPeak) -> String {
  self.chr +
  ":" +
  self.start.to_string() +
  "-" +
  self.end.to_string() +
  " (score=" +
  self.score.to_string() +
  ", summit=" +
  self.summit.to_string() +
  ")"
}

///|
/// String representation of DBAResult.
pub fn DBAResult::to_string(self : DBAResult) -> String {
  let (c1, c2) = self.contrast
  let mut result = "DBAResult: " + c1 + " vs " + c2 + "\n"
  result = result + "  Peaks: " + self.peaks.length().to_string() + "\n"
  result = result + "  Significant: " + self.n_significant.to_string() + "\n"
  result = result + "  Peak\tP-value\tFDR\tFoldChange\n"
  let mut i = 0
  while i < self.peaks.length() {
    result = result +
      self.peaks[i].to_string() +
      "\t" +
      self.pvalues[i].to_string() +
      "\t" +
      self.fdrs[i].to_string() +
      "\t" +
      self.fold_changes[i].to_string() +
      "\n"
    i = i + 1
  }
  result
}

///|
/// String representation of OverlapResult.
pub fn OverlapResult::to_string(self : OverlapResult) -> String {
  let mut result = "OverlapResult:\n"
  result = result +
    "  Consensus peaks: " +
    self.consensus_peaks.length().to_string() +
    "\n"
  result = result + "  Min overlap: " + self.min_overlap.to_string() + "\n"
  result = result + "  Sample names:"
  let mut i = 0
  while i < self.sample_names.length() {
    result = result + " " + self.sample_names[i]
    i = i + 1
  }
  result = result + "\n"
  result = result + "  Overlap matrix:\n"
  let mut r = 0
  while r < self.overlap_matrix.length() {
    result = result + "    ["
    let mut c = 0
    while c < self.overlap_matrix[r].length() {
      result = result + self.overlap_matrix[r][c].to_string()
      if c < self.overlap_matrix[r].length() - 1 {
        result = result + ", "
      }
      c = c + 1
    }
    result = result + "]\n"
    r = r + 1
  }
  result
}

///|
/// Square root via Newton's method.
pub fn diffbind_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x == 1.0 {
    return 1.0
  }
  let mut guess = x / 2.0
  if guess < 1.0 {
    guess = 1.0
  }
  let mut i = 0
  while i < 60 {
    let next = (guess + x / guess) / 2.0
    let diff = next - guess
    if diff < 0.0 {
      if -diff < 0.00000000000001 {
        return next
      }
    } else if diff < 0.00000000000001 {
      return next
    }
    guess = next
    i = i + 1
  }
  guess
}

///|
/// Natural logarithm via series expansion.
pub fn diffbind_ln(x : Double) -> Double {
  if x <= 0.0 {
    return -10000000000.0
  }
  if x == 1.0 {
    return 0.0
  }
  // Use ln(x) = 2 * sum_{k=0}^{inf} y^(2k+1)/(2k+1), where y = (x-1)/(x+1)
  // This converges for x > 0. For large x, reduce via ln(x) = ln(x/e^k) + k.
  let mut v = x
  let mut offset = 0.0
  while v > 2.0 {
    v = v / 2.718281828459045
    offset = offset + 1.0
  }
  while v < 0.5 {
    v = v * 2.718281828459045
    offset = offset - 1.0
  }
  let y = (v - 1.0) / (v + 1.0)
  let y2 = y * y
  let mut sum = y
  let mut term = y
  let mut i = 1
  while i <= 50 {
    term = term * y2
    let delta = term / (2 * i + 1).to_double()
    sum = sum + delta
    if delta < 0.0 {
      if -delta < 0.0000000000000001 {
        return 2.0 * sum + offset
      }
    } else if delta < 0.0000000000000001 {
      return 2.0 * sum + offset
    }
    i = i + 1
  }
  2.0 * sum + offset
}

///|
/// Exponential function via Taylor series.
pub fn diffbind_exp(x : Double) -> Double {
  if x == 0.0 {
    return 1.0
  }
  if x < 0.0 {
    return 1.0 / diffbind_exp(-x)
  }
  // Split x = k + r where k is integer part, r is fractional.
  // exp(x) = e^k * exp(r). Use e = 2.718281828459045.
  let mut k = 0
  let mut r = x
  while r >= 1.0 {
    r = r - 1.0
    k = k + 1
  }
  let mut sum = 1.0
  let mut term = 1.0
  let mut i = 1
  while i <= 30 {
    term = term * r / i.to_double()
    sum = sum + term
    if term < 0.0 {
      if -term < 0.0000000000000001 {
        break
      }
    } else if term < 0.0000000000000001 {
      break
    }
    i = i + 1
  }
  let mut e_pow_k = 1.0
  let mut j = 0
  while j < k {
    e_pow_k = e_pow_k * 2.718281828459045
    j = j + 1
  }
  e_pow_k * sum
}

///|
/// Base-2 logarithm.
pub fn diffbind_log2(x : Double) -> Double {
  diffbind_ln(x) / diffbind_ln(2.0)
}

///|
/// Base-10 logarithm.
pub fn diffbind_log10(x : Double) -> Double {
  diffbind_ln(x) / diffbind_ln(10.0)
}

///|
/// Standard normal CDF (Phi) via Abramowitz-Stegun approximation.
pub fn diffbind_normal_cdf(z : Double) -> Double {
  let abs_z = if z < 0.0 { -z } else { z }
  let t = 1.0 / (1.0 + 0.2316419 * abs_z)
  let d = 0.3989423 * diffbind_exp(-abs_z * abs_z / 2.0)
  let p = 1.0 -
    d *
    t *
    (
      0.3193815 +
      t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274)))
    )
  if z >= 0.0 {
    p
  } else {
    1.0 - p
  }
}

///|
/// Check if two intervals [start1, end1] and [start2, end2] overlap.
pub fn diffbind_is_overlapping(
  start1 : Int,
  end1 : Int,
  start2 : Int,
  end2 : Int,
) -> Bool {
  start1 <= end2 && start2 <= end1
}

///|
/// Sort peaks by chromosome then start position (bubble sort).
fn diffbind_sort_peaks(peaks : Array[DBAPeak]) -> Array[DBAPeak] {
  let sorted = peaks.copy()
  let n = sorted.length()
  let mut i = 0
  while i < n - 1 {
    let mut j = 0
    while j < n - i - 1 {
      let a = sorted[j]
      let b = sorted[j + 1]
      let swap = if a.chr != b.chr { a.chr > b.chr } else { a.start > b.start }
      if swap {
        let tmp = sorted[j]
        sorted[j] = sorted[j + 1]
        sorted[j + 1] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  sorted
}

///|
/// Sort an array of doubles (bubble sort), returning a new array.
fn diffbind_sort_doubles(arr : Array[Double]) -> Array[Double] {
  let sorted = arr.copy()
  let n = sorted.length()
  let mut i = 0
  while i < n - 1 {
    let mut j = 0
    while j < n - i - 1 {
      if sorted[j] > sorted[j + 1] {
        let tmp = sorted[j]
        sorted[j] = sorted[j + 1]
        sorted[j + 1] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  sorted
}

///|
/// Add a sample to the config, returning a new config.
pub fn diffbind_peakset(config : DBAConfig, sample : DBASample) -> DBAConfig {
  let new_samples = config.samples.copy()
  new_samples.push(sample)
  let new_peaks = config.peaks.copy()
  new_peaks.push(sample.peaks)
  let new_conditions = config.conditions.copy()
  new_conditions.push(sample.condition)
  let new_replicates = config.replicates.copy()
  new_replicates.push(sample.replicate)
  DBAConfig::{
    samples: new_samples,
    peaks: new_peaks,
    peak_format: config.peak_format,
    score: config.score,
    conditions: new_conditions,
    replicates: new_replicates,
  }
}

///|
/// Count peaks overlapping each query range across all samples.
/// Returns a Map from "chr:start-end" string to peak count.
pub fn diffbind_count(
  config : DBAConfig,
  ranges : Array[(String, Int, Int)],
) -> Map[String, Int] {
  let result : Map[String, Int] = Map([], capacity=ranges.length())
  let mut i = 0
  while i < ranges.length() {
    let (chr, start, end) = ranges[i]
    let key = chr + ":" + start.to_string() + "-" + end.to_string()
    let mut count = 0
    let mut s = 0
    while s < config.samples.length() {
      let sample = config.samples[s]
      let mut p = 0
      while p < sample.peaks.length() {
        let peak = sample.peaks[p]
        if peak.chr == chr &&
          diffbind_is_overlapping(peak.start, peak.end, start, end) {
          count = count + 1
        }
        p = p + 1
      }
      s = s + 1
    }
    result[key] = count
    i = i + 1
  }
  result
}

///|
/// Identify consensus peaks present in at least min_overlap fraction of samples.
pub fn diffbind_consensus_peaks(
  samples : Array[DBASample],
  min_overlap : Double,
) -> Array[DBAPeak] {
  if samples.length() == 0 {
    return Array::new()
  }
  let n_samples = samples.length()
  let threshold = min_overlap * n_samples.to_double()
  // Collect all peaks across samples.
  let all_peaks : Array[DBAPeak] = Array::new()
  let mut s = 0
  while s < n_samples {
    let mut p = 0
    while p < samples[s].peaks.length() {
      all_peaks.push(samples[s].peaks[p])
      p = p + 1
    }
    s = s + 1
  }
  if all_peaks.length() == 0 {
    return Array::new()
  }
  // Sort peaks by chromosome then start.
  let sorted = diffbind_sort_peaks(all_peaks)
  // Merge overlapping peaks within the same chromosome.
  let merged : Array[DBAPeak] = Array::new()
  merged.push(sorted[0])
  let mut i = 1
  while i < sorted.length() {
    let last = merged[merged.length() - 1]
    let current = sorted[i]
    if current.chr == last.chr && current.start <= last.end {
      let new_end = if current.end > last.end { current.end } else { last.end }
      let new_summit = if current.summit > last.summit {
        current.summit
      } else {
        last.summit
      }
      let new_score = if current.score > last.score {
        current.score
      } else {
        last.score
      }
      merged[merged.length() - 1] = DBAPeak::{
        chr: last.chr,
        start: last.start,
        end: new_end,
        score: new_score,
        summit: new_summit,
      }
    } else {
      merged.push(current)
    }
    i = i + 1
  }
  // For each merged peak, count how many samples contain an overlapping peak.
  let consensus : Array[DBAPeak] = Array::new()
  let mut m = 0
  while m < merged.length() {
    let mut sample_count = 0
    let mut s2 = 0
    while s2 < n_samples {
      let mut found = false
      let mut p2 = 0
      while p2 < samples[s2].peaks.length() {
        let peak = samples[s2].peaks[p2]
        if peak.chr == merged[m].chr &&
          diffbind_is_overlapping(
            peak.start,
            peak.end,
            merged[m].start,
            merged[m].end,
          ) {
          found = true
          break
        }
        p2 = p2 + 1
      }
      if found {
        sample_count = sample_count + 1
      }
      s2 = s2 + 1
    }
    if sample_count.to_double() >= threshold {
      consensus.push(merged[m])
    }
    m = m + 1
  }
  consensus
}

///|
/// Build overlap matrix and consensus peaks for a config.
pub fn diffbind_overlap(
  config : DBAConfig,
  min_overlap : Double,
) -> OverlapResult {
  let consensus = diffbind_consensus_peaks(config.samples, min_overlap)
  let sample_names : Array[String] = Array::new()
  let mut s = 0
  while s < config.samples.length() {
    sample_names.push(config.samples[s].sample_id)
    s = s + 1
  }
  let overlap_matrix : Array[Array[Int]] = Array::new()
  let mut i = 0
  while i < consensus.length() {
    let row : Array[Int] = Array::make(config.samples.length(), 0)
    let mut j = 0
    while j < config.samples.length() {
      let sample = config.samples[j]
      let mut found = 0
      let mut k = 0
      while k < sample.peaks.length() {
        let peak = sample.peaks[k]
        if peak.chr == consensus[i].chr &&
          diffbind_is_overlapping(
            peak.start,
            peak.end,
            consensus[i].start,
            consensus[i].end,
          ) {
          found = 1
          break
        }
        k = k + 1
      }
      row[j] = found
      j = j + 1
    }
    overlap_matrix.push(row)
    i = i + 1
  }
  OverlapResult::{
    consensus_peaks: consensus,
    overlap_matrix,
    sample_names,
    min_overlap,
  }
}

///|
/// Extract consensus peaks from an OverlapResult.
pub fn diffbind_find_consensus_peaks(
  overlap_result : OverlapResult,
) -> Array[DBAPeak] {
  overlap_result.consensus_peaks.copy()
}

///|
/// Merge two arrays of peaks, combining overlapping peaks on the same chromosome.
pub fn diffbind_merge_peaks(
  peaks1 : Array[DBAPeak],
  peaks2 : Array[DBAPeak],
) -> Array[DBAPeak] {
  let combined : Array[DBAPeak] = Array::new()
  let mut i = 0
  while i < peaks1.length() {
    combined.push(peaks1[i])
    i = i + 1
  }
  let mut j = 0
  while j < peaks2.length() {
    combined.push(peaks2[j])
    j = j + 1
  }
  if combined.length() == 0 {
    return Array::new()
  }
  let sorted = diffbind_sort_peaks(combined)
  let merged : Array[DBAPeak] = Array::new()
  merged.push(sorted[0])
  let mut k = 1
  while k < sorted.length() {
    let last = merged[merged.length() - 1]
    let current = sorted[k]
    if current.chr == last.chr && current.start <= last.end {
      let new_end = if current.end > last.end { current.end } else { last.end }
      let new_score = if current.score > last.score {
        current.score
      } else {
        last.score
      }
      merged[merged.length() - 1] = DBAPeak::{
        chr: last.chr,
        start: last.start,
        end: new_end,
        score: new_score,
        summit: last.summit,
      }
    } else {
      merged.push(current)
    }
    k = k + 1
  }
  merged
}

///|
/// Compute TMM normalization factors for a count matrix (peaks x samples).
/// Uses sample 0 as reference. Returns one factor per sample.
pub fn diffbind_compute_tmm(counts : Array[Array[Double]]) -> Array[Double] {
  let n_genes = counts.length()
  if n_genes == 0 {
    return Array::new()
  }
  let n_samples = counts[0].length()
  let factors : Array[Double] = Array::make(n_samples, 1.0)
  if n_samples <= 1 {
    return factors
  }
  // Compute library sizes.
  let lib_sizes : Array[Double] = Array::make(n_samples, 0.0)
  let mut j = 0
  while j < n_samples {
    let mut sum = 0.0
    let mut i = 0
    while i < n_genes {
      sum = sum + counts[i][j]
      i = i + 1
    }
    lib_sizes[j] = sum
    j = j + 1
  }
  let ref_lib = lib_sizes[0]
  if ref_lib <= 0.0 {
    return factors
  }
  // For each non-reference sample, compute TMM factor.
  let mut s = 1
  while s < n_samples {
    let lib_s = lib_sizes[s]
    if lib_s <= 0.0 {
      s = s + 1
      continue
    }
    // Collect (M, A) values for genes with positive counts in both samples.
    let m_values : Array[Double] = Array::new()
    let a_values : Array[Double] = Array::new()
    let weights : Array[Double] = Array::new()
    let mut g = 0
    while g < n_genes {
      let y_ref = counts[g][0]
      let y_s = counts[g][s]
      if y_ref > 0.0 && y_s > 0.0 {
        // Normalize by library size.
        let n_ref = y_ref / ref_lib
        let n_s = y_s / lib_s
        if n_ref > 0.0 && n_s > 0.0 {
          let m = diffbind_log2(n_s) - diffbind_log2(n_ref)
          let a = (diffbind_log2(n_s) + diffbind_log2(n_ref)) / 2.0
          m_values.push(m)
          a_values.push(a)
          // Asymptotic variance weight.
          let denom = n_ref + n_s
          let w = if denom > 0.0 { n_ref * n_s / denom } else { 0.0 }
          weights.push(w)
        }
      }
      g = g + 1
    }
    if m_values.length() == 0 {
      s = s + 1
      continue
    }
    // Sort M values and trim 30% from each end (15% top, 15% bottom by M).
    let m_sorted = diffbind_sort_doubles(m_values)
    let n_valid = m_sorted.length()
    let trim_count = n_valid * 30 / 100
    let lo = trim_count / 2
    let hi = n_valid - trim_count / 2 - 1
    let actual_lo = if lo > 0 { lo } else { 0 }
    let actual_hi = if hi < n_valid - 1 { hi } else { n_valid - 1 }
    // Compute weighted mean of M values within trimmed range.
    // Use the original (unsorted) m_values and weights, filtering by sorted M range.
    let m_low = m_sorted[actual_lo]
    let m_high = m_sorted[actual_hi]
    let mut weighted_sum = 0.0
    let mut weight_total = 0.0
    let mut idx = 0
    while idx < m_values.length() {
      let m = m_values[idx]
      if m >= m_low && m <= m_high {
        let w = if weights[idx] > 0.0 { weights[idx] } else { 1.0 }
        weighted_sum = weighted_sum + m * w
        weight_total = weight_total + w
      }
      idx = idx + 1
    }
    let mean_m = if weight_total > 0.0 {
      weighted_sum / weight_total
    } else {
      0.0
    }
    factors[s] = diffbind_exp(mean_m * diffbind_ln(2.0))
    s = s + 1
  }
  factors
}

///|
/// Two-sample negative binomial test.
/// Returns two-sided p-value for the difference between count1 and count2.
pub fn diffbind_negative_binomial_test(
  count1 : Double,
  count2 : Double,
  dispersion : Double,
) -> Double {
  let total = count1 + count2
  if total <= 0.0 {
    return 1.0
  }
  let mean = total / 2.0
  let phi = if dispersion > 0.0001 { dispersion } else { 0.0001 }
  // NB variance: mean + mean^2 / phi
  let variance = mean + mean * mean / phi
  let std_dev = diffbind_sqrt(variance * 2.0)
  if std_dev <= 0.0 {
    return 1.0
  }
  let z = (count1 - count2) / std_dev
  let abs_z = if z < 0.0 { -z } else { z }
  // Two-sided p-value.
  let p = 2.0 * (1.0 - diffbind_normal_cdf(abs_z))
  if p < 0.0 {
    0.0
  } else if p > 1.0 {
    1.0
  } else {
    p
  }
}

///|
/// Benjamini-Hochberg FDR correction.
pub fn diffbind_bh_fdr(pvalues : Array[Double]) -> Array[Double] {
  let n = pvalues.length()
  if n == 0 {
    return Array::new()
  }
  // Create (original_index, pvalue) pairs.
  let indexed : Array[(Int, Double)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((i, pvalues[i]))
    i = i + 1
  }
  // Sort by p-value ascending (bubble sort on the indexed array).
  let sorted = indexed.copy()
  let mut a = 0
  while a < n - 1 {
    let mut b = 0
    while b < n - a - 1 {
      if sorted[b].1 > sorted[b + 1].1 {
        let tmp = sorted[b]
        sorted[b] = sorted[b + 1]
        sorted[b + 1] = tmp
      }
      b = b + 1
    }
    a = a + 1
  }
  // Compute adjusted p-values, enforcing monotonicity from the largest.
  let adjusted : Array[Double] = Array::make(n, 1.0)
  let mut prev = 1.0
  let mut k = n - 1
  while k >= 0 {
    let rank = k + 1
    let raw = sorted[k].1 * n.to_double() / rank.to_double()
    let adj = if raw > prev { prev } else { raw }
    let capped = if adj > 1.0 { 1.0 } else { adj }
    adjusted[k] = capped
    prev = capped
    k = k - 1
  }
  // Map back to original ordering.
  let result : Array[Double] = Array::make(n, 1.0)
  let mut j = 0
  while j < n {
    result[sorted[j].0] = adjusted[j]
    j = j + 1
  }
  result
}

///|
/// Build a count matrix (peaks x samples) from consensus peaks.
fn diffbind_build_count_matrix(
  config : DBAConfig,
  consensus : Array[DBAPeak],
) -> Array[Array[Double]] {
  let n_peaks = consensus.length()
  let n_samples = config.samples.length()
  let counts : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_peaks {
    let row : Array[Double] = Array::make(n_samples, 0.0)
    let mut s = 0
    while s < n_samples {
      let sample = config.samples[s]
      let mut best_score = 0.0
      let mut found = false
      let mut p = 0
      while p < sample.peaks.length() {
        let peak = sample.peaks[p]
        if peak.chr == consensus[i].chr &&
          diffbind_is_overlapping(
            peak.start,
            peak.end,
            consensus[i].start,
            consensus[i].end,
          ) {
          if !found || peak.score > best_score {
            best_score = peak.score
            found = true
          }
        }
        p = p + 1
      }
      row[s] = if found { best_score } else { 0.0 }
      s = s + 1
    }
    counts.push(row)
    i = i + 1
  }
  counts
}

///|
/// Normalize the count matrix using the specified method.
/// Supported methods: "TMM", "RPKM", "CPM", "library" (default).
pub fn diffbind_normalize(
  config : DBAConfig,
  method : String,
) -> Array[Array[Double]] {
  let consensus = diffbind_consensus_peaks(config.samples, 0.5)
  let counts = diffbind_build_count_matrix(config, consensus)
  let n_peaks = counts.length()
  if n_peaks == 0 {
    return counts
  }
  let n_samples = counts[0].length()
  if n_samples == 0 {
    return counts
  }
  if method == "TMM" {
    let factors = diffbind_compute_tmm(counts)
    let mut i = 0
    while i < n_peaks {
      let mut s = 0
      while s < n_samples {
        let f = factors[s]
        counts[i][s] = if f > 0.0 { counts[i][s] / f } else { counts[i][s] }
        s = s + 1
      }
      i = i + 1
    }
  } else if method == "RPKM" {
    let lib_sizes : Array[Double] = Array::make(n_samples, 0.0)
    let mut s = 0
    while s < n_samples {
      let mut sum = 0.0
      let mut i = 0
      while i < n_peaks {
        sum = sum + counts[i][s]
        i = i + 1
      }
      lib_sizes[s] = sum
      s = s + 1
    }
    let mut i = 0
    while i < n_peaks {
      let length = (consensus[i].end - consensus[i].start + 1).to_double() /
        1000.0
      let mut s = 0
      while s < n_samples {
        if lib_sizes[s] > 0.0 && length > 0.0 {
          counts[i][s] = counts[i][s] / lib_sizes[s] * 1000000000.0 / length
        } else {
          counts[i][s] = 0.0
        }
        s = s + 1
      }
      i = i + 1
    }
  } else if method == "CPM" {
    let lib_sizes : Array[Double] = Array::make(n_samples, 0.0)
    let mut s = 0
    while s < n_samples {
      let mut sum = 0.0
      let mut i = 0
      while i < n_peaks {
        sum = sum + counts[i][s]
        i = i + 1
      }
      lib_sizes[s] = sum
      s = s + 1
    }
    let mut i = 0
    while i < n_peaks {
      let mut s = 0
      while s < n_samples {
        if lib_sizes[s] > 0.0 {
          counts[i][s] = counts[i][s] / lib_sizes[s] * 1000000.0
        } else {
          counts[i][s] = 0.0
        }
        s = s + 1
      }
      i = i + 1
    }
  } else {
    // "library" or default: divide by library size (in millions).
    let lib_sizes : Array[Double] = Array::make(n_samples, 0.0)
    let mut s = 0
    while s < n_samples {
      let mut sum = 0.0
      let mut i = 0
      while i < n_peaks {
        sum = sum + counts[i][s]
        i = i + 1
      }
      lib_sizes[s] = sum
      s = s + 1
    }
    let mut i = 0
    while i < n_peaks {
      let mut s = 0
      while s < n_samples {
        if lib_sizes[s] > 0.0 {
          counts[i][s] = counts[i][s] / lib_sizes[s]
        }
        s = s + 1
      }
      i = i + 1
    }
  }
  counts
}

///|
/// Estimate a global dispersion value from a count matrix.
fn diffbind_estimate_dispersion(
  counts : Array[Array[Double]],
  group1_idx : Array[Int],
  group2_idx : Array[Int],
) -> Double {
  let n_peaks = counts.length()
  if n_peaks == 0 {
    return 0.5
  }
  let mut total_dispersion = 0.0
  let mut n_valid = 0
  let mut i = 0
  while i < n_peaks {
    // Compute mean and variance across all samples in both groups.
    let mut sum = 0.0
    let mut n = 0
    let mut g = 0
    while g < group1_idx.length() {
      sum = sum + counts[i][group1_idx[g]]
      n = n + 1
      g = g + 1
    }
    g = 0
    while g < group2_idx.length() {
      sum = sum + counts[i][group2_idx[g]]
      n = n + 1
      g = g + 1
    }
    if n < 2 {
      i = i + 1
      continue
    }
    let mean = sum / n.to_double()
    if mean <= 0.0 {
      i = i + 1
      continue
    }
    let mut var_sum = 0.0
    g = 0
    while g < group1_idx.length() {
      let diff = counts[i][group1_idx[g]] - mean
      var_sum = var_sum + diff * diff
      g = g + 1
    }
    g = 0
    while g < group2_idx.length() {
      let diff = counts[i][group2_idx[g]] - mean
      var_sum = var_sum + diff * diff
      g = g + 1
    }
    let variance = var_sum / (n - 1).to_double()
    // NB: variance = mean + mean^2 / phi  =>  phi = mean^2 / (variance - mean)
    if variance > mean {
      let phi = mean * mean / (variance - mean)
      if phi > 0.0 {
        total_dispersion = total_dispersion + phi
        n_valid = n_valid + 1
      }
    }
    i = i + 1
  }
  if n_valid == 0 {
    0.5
  } else {
    let avg = total_dispersion / n_valid.to_double()
    if avg < 0.1 {
      0.1
    } else if avg > 10.0 {
      10.0
    } else {
      avg
    }
  }
}

///|
/// Perform differential binding analysis between two conditions.
/// Uses a negative binomial model similar to DESeq2/edgeR.
pub fn diffbind_analyze(
  config : DBAConfig,
  condition1 : String,
  condition2 : String,
) -> DBAResult {
  let empty_result = DBAResult::{
    peaks: Array::new(),
    pvalues: Array::new(),
    fdrs: Array::new(),
    fold_changes: Array::new(),
    contrast: (condition1, condition2),
    n_significant: 0,
  }
  if config.samples.length() == 0 {
    return empty_result
  }
  // Identify sample indices for each condition.
  let group1_idx : Array[Int] = Array::new()
  let group2_idx : Array[Int] = Array::new()
  let mut s = 0
  while s < config.samples.length() {
    if config.samples[s].condition == condition1 {
      group1_idx.push(s)
    } else if config.samples[s].condition == condition2 {
      group2_idx.push(s)
    }
    s = s + 1
  }
  if group1_idx.length() == 0 || group2_idx.length() == 0 {
    return empty_result
  }
  // Build consensus peaks across all samples.
  let consensus = diffbind_consensus_peaks(config.samples, 0.5)
  if consensus.length() == 0 {
    return empty_result
  }
  // Build count matrix.
  let counts = diffbind_build_count_matrix(config, consensus)
  // Apply TMM normalization.
  let factors = diffbind_compute_tmm(counts)
  let n_peaks = counts.length()
  // Estimate dispersion.
  let dispersion = diffbind_estimate_dispersion(counts, group1_idx, group2_idx)
  // Compute per-peak statistics.
  let pvalues : Array[Double] = Array::make(n_peaks, 1.0)
  let fold_changes : Array[Double] = Array::make(n_peaks, 0.0)
  let mut i = 0
  while i < n_peaks {
    let mut sum1 = 0.0
    let mut g1 = 0
    while g1 < group1_idx.length() {
      let f = factors[group1_idx[g1]]
      let v = if f > 0.0 {
        counts[i][group1_idx[g1]] / f
      } else {
        counts[i][group1_idx[g1]]
      }
      sum1 = sum1 + v
      g1 = g1 + 1
    }
    let mut sum2 = 0.0
    let mut g2 = 0
    while g2 < group2_idx.length() {
      let f = factors[group2_idx[g2]]
      let v = if f > 0.0 {
        counts[i][group2_idx[g2]] / f
      } else {
        counts[i][group2_idx[g2]]
      }
      sum2 = sum2 + v
      g2 = g2 + 1
    }
    // Add pseudo-count to avoid log(0).
    let pseudo1 = sum1 + 0.5
    let pseudo2 = sum2 + 0.5
    fold_changes[i] = diffbind_log2(pseudo1 / pseudo2)
    pvalues[i] = diffbind_negative_binomial_test(sum1, sum2, dispersion)
    i = i + 1
  }
  // BH-FDR correction.
  let fdrs = diffbind_bh_fdr(pvalues)
  // Count significant peaks (FDR < 0.05).
  let mut n_sig = 0
  let mut k = 0
  while k < fdrs.length() {
    if fdrs[k] < 0.05 {
      n_sig = n_sig + 1
    }
    k = k + 1
  }
  DBAResult::{
    peaks: consensus,
    pvalues,
    fdrs,
    fold_changes,
    contrast: (condition1, condition2),
    n_significant: n_sig,
  }
}

///|
/// Report significant peaks with FDR below the given threshold.
pub fn diffbind_report(
  result : DBAResult,
  threshold : Double,
) -> Array[DBAPeak] {
  let significant : Array[DBAPeak] = Array::new()
  let mut i = 0
  while i < result.peaks.length() {
    if result.fdrs[i] < threshold {
      significant.push(result.peaks[i])
    }
    i = i + 1
  }
  significant
}

///|
/// Generate volcano plot data: (fold_change, -log10(pvalue)) pairs.
pub fn diffbind_volcano_data(result : DBAResult) -> Array[(Double, Double)] {
  let data : Array[(Double, Double)] = Array::new()
  let mut i = 0
  while i < result.peaks.length() {
    let fc = result.fold_changes[i]
    let p = result.pvalues[i]
    let neg_log10_p = if p > 0.0 { -diffbind_log10(p) } else { 300.0 }
    data.push((fc, neg_log10_p))
    i = i + 1
  }
  data
}

///|
/// Generate MA plot data: (mean_expression_log2, log2_fold_change) pairs.
pub fn diffbind_ma_data(result : DBAResult) -> Array[(Double, Double)] {
  let data : Array[(Double, Double)] = Array::new()
  let mut i = 0
  while i < result.peaks.length() {
    let fc = result.fold_changes[i]
    // Approximate mean expression: derive from fold change and contrast.
    // Use the peak score as a proxy for base expression.
    let base_score = result.peaks[i].score
    let mean_log2 = if base_score > 0.0 {
      diffbind_log2(base_score)
    } else {
      0.0
    }
    data.push((mean_log2, fc))
    i = i + 1
  }
  data
}

///|
/// Create an example DBAConfig with two conditions and multiple replicates.
pub fn diffbind_create_example_config() -> DBAConfig {
  let samples : Array[DBASample] = Array::new()
  // Condition A: 2 replicates.
  samples.push(DBASample::{
    sample_id: "A1",
    condition: "ConditionA",
    replicate: 1,
    peak_file: "",
    bam_file: "",
    peaks: diffbind_create_example_peaks(),
  })
  samples.push(DBASample::{
    sample_id: "A2",
    condition: "ConditionA",
    replicate: 2,
    peak_file: "",
    bam_file: "",
    peaks: diffbind_create_example_peaks(),
  })
  // Condition B: 2 replicates with elevated scores to mimic differential binding.
  samples.push(DBASample::{
    sample_id: "B1",
    condition: "ConditionB",
    replicate: 1,
    peak_file: "",
    bam_file: "",
    peaks: diffbind_create_example_peaks_alt(),
  })
  samples.push(DBASample::{
    sample_id: "B2",
    condition: "ConditionB",
    replicate: 2,
    peak_file: "",
    bam_file: "",
    peaks: diffbind_create_example_peaks_alt(),
  })
  DBAConfig::new(samples)
}

///|
/// Create an example array of peaks (ConditionA baseline).
pub fn diffbind_create_example_peaks() -> Array[DBAPeak] {
  let peaks : Array[DBAPeak] = Array::new()
  peaks.push(DBAPeak::{
    chr: "chr1",
    start: 1000,
    end: 2000,
    score: 50.0,
    summit: 1500,
  })
  peaks.push(DBAPeak::{
    chr: "chr1",
    start: 3000,
    end: 4000,
    score: 80.0,
    summit: 3500,
  })
  peaks.push(DBAPeak::{
    chr: "chr2",
    start: 1500,
    end: 2500,
    score: 30.0,
    summit: 2000,
  })
  peaks.push(DBAPeak::{
    chr: "chr2",
    start: 5000,
    end: 6000,
    score: 70.0,
    summit: 5500,
  })
  peaks.push(DBAPeak::{
    chr: "chr3",
    start: 1000,
    end: 3000,
    score: 90.0,
    summit: 2000,
  })
  peaks
}

///|
/// Create an example array of peaks for ConditionB (with shifted scores).
pub fn diffbind_create_example_peaks_alt() -> Array[DBAPeak] {
  let peaks : Array[DBAPeak] = Array::new()
  peaks.push(DBAPeak::{
    chr: "chr1",
    start: 1000,
    end: 2000,
    score: 200.0,
    summit: 1500,
  })
  peaks.push(DBAPeak::{
    chr: "chr1",
    start: 3000,
    end: 4000,
    score: 75.0,
    summit: 3500,
  })
  peaks.push(DBAPeak::{
    chr: "chr2",
    start: 1500,
    end: 2500,
    score: 180.0,
    summit: 2000,
  })
  peaks.push(DBAPeak::{
    chr: "chr2",
    start: 5000,
    end: 6000,
    score: 65.0,
    summit: 5500,
  })
  peaks.push(DBAPeak::{
    chr: "chr3",
    start: 1000,
    end: 3000,
    score: 95.0,
    summit: 2000,
  })
  peaks
}

///|
/// Define a contrast for differential binding analysis.
/// Returns a tuple of (condition1, condition2) for use with diffbind_analyze.
pub fn diffbind_contrast(
  config : DBAConfig,
  condition1 : String,
  condition2 : String,
) -> (String, String) {
  let mut found1 = false
  let mut found2 = false
  let mut i = 0
  while i < config.samples.length() {
    if config.samples[i].condition == condition1 {
      found1 = true
    } else if config.samples[i].condition == condition2 {
      found2 = true
    }
    if found1 && found2 {
      break
    }
    i = i + 1
  }
  if !found1 {
    panic()
  }
  if !found2 {
    panic()
  }
  (condition1, condition2)
}

///|
/// Transpose a matrix (peaks x samples -> samples x peaks).
fn diffbind_transpose(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  if matrix.length() == 0 {
    return Array::new()
  }
  let n_rows = matrix.length()
  let n_cols = matrix[0].length()
  let transposed : Array[Array[Double]] = Array::make(n_cols, Array::make(n_rows, 0.0))
  let mut i = 0
  while i < n_rows {
    let mut j = 0
    while j < n_cols {
      transposed[j][i] = matrix[i][j]
      j = j + 1
    }
    i = i + 1
  }
  transposed
}

///|
/// Compute mean of an array.
fn diffbind_mean(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < arr.length() {
    sum = sum + arr[i]
    i = i + 1
  }
  sum / arr.length().to_double()
}

///|
/// Extract column from matrix.
fn diffbind_extract_column(matrix : Array[Array[Double]], col_idx : Int) -> Array[Double] {
  let n_rows = matrix.length()
  let col : Array[Double] = Array::make(n_rows, 0.0)
  let mut i = 0
  while i < n_rows {
    col[i] = matrix[i][col_idx]
    i = i + 1
  }
  col
}

///|
/// Compute sum of array elements.
fn diffbind_sum(arr : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < arr.length() {
    sum = sum + arr[i]
    i = i + 1
  }
  sum
}

///|
/// Compute covariance between two arrays.
fn diffbind_covariance(x : Array[Double], y : Array[Double]) -> Double {
  if x.length() != y.length() || x.length() < 2 {
    return 0.0
  }
  let mean_x = diffbind_mean(x)
  let mean_y = diffbind_mean(y)
  let mut sum = 0.0
  let mut i = 0
  while i < x.length() {
    sum = sum + (x[i] - mean_x) * (y[i] - mean_y)
    i = i + 1
  }
  sum / (x.length() - 1).to_double()
}

///|
/// Compute variance of an array.
fn diffbind_variance(arr : Array[Double]) -> Double {
  diffbind_covariance(arr, arr)
}

///|
/// Compute Pearson correlation between two arrays.
fn diffbind_correlation(x : Array[Double], y : Array[Double]) -> Double {
  let cov = diffbind_covariance(x, y)
  let std_x = diffbind_sqrt(diffbind_variance(x))
  let std_y = diffbind_sqrt(diffbind_variance(y))
  if std_x == 0.0 || std_y == 0.0 {
    return 0.0
  }
  cov / (std_x * std_y)
}

///|
/// Compute PCA on a matrix (samples x peaks).
/// Returns principal components (samples x components) and explained variance ratios.
fn diffbind_pca(matrix : Array[Array[Double]]) -> (Array[Array[Double]], Array[Double]) {
  if matrix.length() == 0 || matrix[0].length() == 0 {
    return (Array::new(), Array::new())
  }
  let n_samples = matrix.length()
  let n_peaks = matrix[0].length()
  let k = if n_samples < n_peaks { n_samples - 1 } else { 3 }
  
  let centered : Array[Array[Double]] = Array::make(n_samples, Array::make(n_peaks, 0.0))
  let mut p = 0
  while p < n_peaks {
    let peak_values : Array[Double] = Array::make(n_samples, 0.0)
    let mut s = 0
    while s < n_samples {
      peak_values[s] = matrix[s][p]
      s = s + 1
    }
    let mean_p = diffbind_mean(peak_values)
    s = 0
    while s < n_samples {
      centered[s][p] = matrix[s][p] - mean_p
      s = s + 1
    }
    p = p + 1
  }
  
  let cov_matrix : Array[Array[Double]] = Array::make(n_peaks, Array::make(n_peaks, 0.0))
  let mut i = 0
  while i < n_peaks {
    let mut j = i
    while j < n_peaks {
      let col_i = diffbind_extract_column(centered, i)
      let col_j = diffbind_extract_column(centered, j)
      let cov = diffbind_covariance(col_i, col_j)
      cov_matrix[i][j] = cov
      cov_matrix[j][i] = cov
      j = j + 1
    }
    i = i + 1
  }
  
  let eigvals : Array[Double] = Array::make(k, 0.0)
  let eigvecs : Array[Array[Double]] = Array::make(k, Array::make(n_peaks, 0.0))
  
  let mut comp = 0
  while comp < k {
    let init_val = 1.0 / diffbind_sqrt(n_peaks.to_double())
    let mut vec : Array[Double] = Array::make(n_peaks, init_val)
    
    let mut val = 0.0
    let mut iter = 0
    while iter < 100 {
      let new_vec : Array[Double] = Array::make(n_peaks, 0.0)
      let mut p2 = 0
      while p2 < n_peaks {
        let mut sum = 0.0
        let mut q = 0
        while q < n_peaks {
          sum = sum + cov_matrix[p2][q] * vec[q]
          q = q + 1
        }
        new_vec[p2] = sum
        p2 = p2 + 1
      }
      
      let mut sum_sq = 0.0
      let mut idx = 0
      while idx < new_vec.length() {
        sum_sq = sum_sq + new_vec[idx] * new_vec[idx]
        idx = idx + 1
      }
      
      let norm = diffbind_sqrt(sum_sq)
      
      let new_vec_norm : Array[Double] = Array::make(n_peaks, 0.0)
      if norm > 0.0 {
        let mut n_idx = 0
        while n_idx < new_vec.length() {
          new_vec_norm[n_idx] = new_vec[n_idx] / norm
          n_idx = n_idx + 1
        }
      }
      
      vec = new_vec_norm
      val = sum_sq
      iter = iter + 1
    }
    
    eigvals[comp] = val
    let mut v_idx = 0
    while v_idx < n_peaks {
      eigvecs[comp][v_idx] = vec[v_idx]
      v_idx = v_idx + 1
    }
    comp = comp + 1
  }
  
  let mut total_var = 0.0
  let mut r_idx = 0
  while r_idx < cov_matrix.length() {
    let mut c_idx = 0
    while c_idx < cov_matrix[r_idx].length() {
      total_var = total_var + cov_matrix[r_idx][c_idx]
      c_idx = c_idx + 1
    }
    r_idx = r_idx + 1
  }
  
  let explained_var : Array[Double] = Array::make(k, 0.0)
  let mut ev_idx = 0
  while ev_idx < k {
    explained_var[ev_idx] = if total_var > 0.0 { eigvals[ev_idx] / total_var } else { 0.0 }
    ev_idx = ev_idx + 1
  }
  
  let pc_scores : Array[Array[Double]] = Array::make(n_samples, Array::make(k, 0.0))
  let mut s2 = 0
  while s2 < n_samples {
    let mut comp2 = 0
    while comp2 < k {
      let mut sum = 0.0
      let mut p3 = 0
      while p3 < n_peaks {
        sum = sum + centered[s2][p3] * eigvecs[comp2][p3]
        p3 = p3 + 1
      }
      pc_scores[s2][comp2] = sum
      comp2 = comp2 + 1
    }
    s2 = s2 + 1
  }
  
  (pc_scores, explained_var)
}

///|
/// Generate PCA plot data for visualization.
pub fn diffbind_plot_pca(config : DBAConfig, method : String) -> String {
  let consensus = diffbind_consensus_peaks(config.samples, 0.5)
  if consensus.length() == 0 {
    return "PCA Plot: No consensus peaks found\n"
  }
  let counts = diffbind_normalize(config, method)
  let matrix = diffbind_transpose(counts)
  
  let (pc_scores, explained_var) = diffbind_pca(matrix)
  if pc_scores.length() == 0 || explained_var.length() < 2 {
    return "PCA Plot: Insufficient data for PCA\n"
  }
  
  let mut plot = "PCA Plot\n"
  plot = plot + "========\n\n"
  plot = plot + "PC1: " + (explained_var[0] * 100.0).to_string() + "% variance\n"
  plot = plot + "PC2: " + (explained_var[1] * 100.0).to_string() + "% variance\n\n"
  
  let mut min_x = pc_scores[0][0]
  let mut max_x = pc_scores[0][0]
  let mut min_y = pc_scores[0][1]
  let mut max_y = pc_scores[0][1]
  let mut i = 1
  while i < pc_scores.length() {
    if pc_scores[i][0] < min_x { min_x = pc_scores[i][0] }
    if pc_scores[i][0] > max_x { max_x = pc_scores[i][0] }
    if pc_scores[i][1] < min_y { min_y = pc_scores[i][1] }
    if pc_scores[i][1] > max_y { max_y = pc_scores[i][1] }
    i = i + 1
  }
  
  let padding = 0.1
  min_x = min_x - (max_x - min_x) * padding
  max_x = max_x + (max_x - min_x) * padding
  min_y = min_y - (max_y - min_y) * padding
  max_y = max_y + (max_y - min_y) * padding
  
  let width = 40
  let height = 20
  
  let canvas : Array[Array[String]] = Array::make(height, Array::make(width, " "))
  
  let mut s = 0
  while s < pc_scores.length() {
    let x = pc_scores[s][0]
    let y = pc_scores[s][1]
    let px = ((x - min_x) / (max_x - min_x) * (width - 1).to_double()).to_int()
    let py = height - 1 - ((y - min_y) / (max_y - min_y) * (height - 1).to_double()).to_int()
    if px >= 0 && px < width && py >= 0 && py < height {
      let sample = config.samples[s]
      canvas[py][px] = if sample.condition == config.conditions[0] { "A" } else { "B" }
    }
    s = s + 1
  }
  
  plot = plot + "Y\n"
  let mut row = height - 1
  while row >= 0 {
    plot = plot + "|"
    let mut col = 0
    while col < width {
      plot = plot + canvas[row][col]
      col = col + 1
    }
    plot = plot + "|\n"
    row = row - 1
  }
  plot = plot + "+" + "-".repeat(width) + "+\n"
  plot = plot + "X\n\n"
  
  let colors : Map[String, String] = Map([], capacity=2)
  let mut cond_idx = 0
  while cond_idx < config.conditions.length() {
    let cond = config.conditions[cond_idx]
    if colors.get(cond).is_none() {
      if colors.length() == 0 {
        colors[cond] = "A"
      } else {
        colors[cond] = "B"
      }
    }
    cond_idx = cond_idx + 1
  }
  
  plot = plot + "Legend:\n"
  let keys = colors.keys().collect()
  let mut k = 0
  while k < keys.length() {
    plot = plot + "  " + colors[keys[k]] + " = " + keys[k] + "\n"
    k = k + 1
  }
  
  plot
}

///|
/// Generate heatmap for differential binding analysis results.
pub fn diffbind_plot_heatmap(result : DBAResult, top_n : Int) -> String {
  if result.peaks.length() == 0 {
    return "Heatmap: No peaks in result\n"
  }
  
  let n = if top_n < result.peaks.length() { top_n } else { result.peaks.length() }
  
  let mut heatmap = "Differential Binding Heatmap\n"
  heatmap = heatmap + "============================\n\n"
  heatmap = heatmap + "Top " + n.to_string() + " peaks\n\n"
  
  let width = if n > 30 { 30 } else { n }
  let height = 10
  
  let mut min_fc = result.fold_changes[0]
  let mut max_fc = result.fold_changes[0]
  let mut i = 1
  while i < result.fold_changes.length() {
    let fc = result.fold_changes[i]
    if fc < min_fc { min_fc = fc }
    if fc > max_fc { max_fc = fc }
    i = i + 1
  }
  
  let mut row = 0
  while row < height {
    heatmap = heatmap + "|"
    let mut col = 0
    while col < width {
      let idx = col % result.fold_changes.length()
      let fc = result.fold_changes[idx]
      let norm = (fc - min_fc) / (max_fc - min_fc + 0.0001)
      let char = if norm < 0.2 { " " } else if norm < 0.4 { "." } else if norm < 0.6 { "+" } else if norm < 0.8 { "*" } else { "#" }
      heatmap = heatmap + char
      col = col + 1
    }
    heatmap = heatmap + "|\n"
    row = row + 1
  }
  
  heatmap = heatmap + "\nColor scale (log2 fold change):\n"
  heatmap = heatmap + "  " + min_fc.to_string() + " <-" + " ".repeat(8) + ".     +     *     #" + "-> " + max_fc.to_string() + "\n"
  
  heatmap
}

///|
/// Generate volcano plot for differential binding analysis results.
pub fn diffbind_plot_volcano(result : DBAResult) -> String {
  let data = diffbind_volcano_data(result)
  if data.length() == 0 {
    return "Volcano Plot: No data\n"
  }
  
  let mut min_fc = data[0].0
  let mut max_fc = data[0].0
  let mut max_neg_log_p = data[0].1
  let mut i = 1
  while i < data.length() {
    if data[i].0 < min_fc { min_fc = data[i].0 }
    if data[i].0 > max_fc { max_fc = data[i].0 }
    if data[i].1 > max_neg_log_p { max_neg_log_p = data[i].1 }
    i = i + 1
  }
  
  let padding = 0.1
  min_fc = min_fc - (max_fc - min_fc) * padding
  max_fc = max_fc + (max_fc - min_fc) * padding
  max_neg_log_p = max_neg_log_p + max_neg_log_p * padding
  
  let width = 50
  let height = 20
  
  let canvas : Array[Array[String]] = Array::make(height, Array::make(width, " "))
  
  let mut j = 0
  while j < data.length() {
    let (fc, neg_log_p) = data[j]
    let px = ((fc - min_fc) / (max_fc - min_fc) * (width - 1).to_double()).to_int()
    let py = height - 1 - (neg_log_p / max_neg_log_p * (height - 1).to_double()).to_int()
    if px >= 0 && px < width && py >= 0 && py < height {
      let idx = j
      let fdr = result.fdrs[idx]
      canvas[py][px] = if fdr < 0.05 { "*" } else { "." }
    }
    j = j + 1
  }
  
  let mid_x = ((0.0 - min_fc) / (max_fc - min_fc) * (width - 1).to_double()).to_int()
  let mut row = height - 1
  while row >= 0 {
    if row == height - 1 {
      canvas[row][mid_x] = "+"
    } else {
      canvas[row][mid_x] = "|"
    }
    row = row - 1
  }
  
  let mut plot = "Volcano Plot\n"
  plot = plot + "============\n\n"
  plot = plot + "Y: -log10(p-value)\n"
  plot = plot + "X: log2(fold change)\n\n"
  
  row = height - 1
  while row >= 0 {
    plot = plot + "|"
    let mut col = 0
    while col < width {
      plot = plot + canvas[row][col]
      col = col + 1
    }
    plot = plot + "|\n"
    row = row - 1
  }
  plot = plot + "+" + "-".repeat(mid_x) + "+" + "-".repeat(width - mid_x - 1) + "+\n"
  plot = plot + "X\n\n"
  
  plot = plot + "Legend:\n"
  plot = plot + "  * = Significant (FDR < 0.05)\n"
  plot = plot + "  . = Not significant\n\n"
  
  let mut sig_count = 0
  let mut up_count = 0
  let mut down_count = 0
  let mut k = 0
  while k < result.fdrs.length() {
    if result.fdrs[k] < 0.05 {
      sig_count = sig_count + 1
      if result.fold_changes[k] > 0 {
        up_count = up_count + 1
      } else {
        down_count = down_count + 1
      }
    }
    k = k + 1
  }
  
  plot = plot + "Summary:\n"
  plot = plot + "  Total peaks: " + result.peaks.length().to_string() + "\n"
  plot = plot + "  Significant (FDR < 0.05): " + sig_count.to_string() + "\n"
  plot = plot + "    Up-regulated: " + up_count.to_string() + "\n"
  plot = plot + "    Down-regulated: " + down_count.to_string() + "\n"
  
  plot
}

///|
/// Generate MA plot for differential binding analysis results.
pub fn diffbind_plot_ma(result : DBAResult) -> String {
  let data = diffbind_ma_data(result)
  if data.length() == 0 {
    return "MA Plot: No data\n"
  }
  
  let mut min_mean = data[0].0
  let mut max_mean = data[0].0
  let mut min_fc = data[0].1
  let mut max_fc = data[0].1
  let mut i = 1
  while i < data.length() {
    if data[i].0 < min_mean { min_mean = data[i].0 }
    if data[i].0 > max_mean { max_mean = data[i].0 }
    if data[i].1 < min_fc { min_fc = data[i].1 }
    if data[i].1 > max_fc { max_fc = data[i].1 }
    i = i + 1
  }
  
  let padding = 0.1
  min_mean = min_mean - (max_mean - min_mean) * padding
  max_mean = max_mean + (max_mean - min_mean) * padding
  min_fc = min_fc - (max_fc - min_fc) * padding
  max_fc = max_fc + (max_fc - min_fc) * padding
  
  let width = 50
  let height = 20
  
  let canvas : Array[Array[String]] = Array::make(height, Array::make(width, " "))
  
  let mut j = 0
  while j < data.length() {
    let (mean_log2, fc) = data[j]
    let px = ((mean_log2 - min_mean) / (max_mean - min_mean) * (width - 1).to_double()).to_int()
    let py = height - 1 - ((fc - min_fc) / (max_fc - min_fc) * (height - 1).to_double()).to_int()
    if px >= 0 && px < width && py >= 0 && py < height {
      let idx = j
      let fdr = result.fdrs[idx]
      canvas[py][px] = if fdr < 0.05 { "*" } else { "." }
    }
    j = j + 1
  }
  
  let mid_y = ((0.0 - min_fc) / (max_fc - min_fc) * (height - 1).to_double()).to_int()
  let mut row = height - 1
  while row >= 0 {
    if row == mid_y {
      let mut col = 0
      while col < width {
        canvas[row][col] = "-"
        col = col + 1
      }
    }
    row = row - 1
  }
  
  let mut plot = "MA Plot\n"
  plot = plot + "=======\n\n"
  plot = plot + "Y: log2(fold change)\n"
  plot = plot + "X: log2(mean expression)\n\n"
  
  row = height - 1
  while row >= 0 {
    plot = plot + "|"
    let mut col = 0
    while col < width {
      plot = plot + canvas[row][col]
      col = col + 1
    }
    plot = plot + "|\n"
    row = row - 1
  }
  plot = plot + "+" + "-".repeat(width) + "+\n"
  plot = plot + "X\n\n"
  
  plot = plot + "Legend:\n"
  plot = plot + "  * = Significant (FDR < 0.05)\n"
  plot = plot + "  . = Not significant\n"
  
  plot
}

///|
/// Create an example DBAResult for testing visualization functions.
pub fn diffbind_create_example_result() -> DBAResult {
  let peaks = diffbind_create_example_peaks()
  let n = peaks.length()
  let pvalues : Array[Double] = Array::make(n, 0.0)
  let fdrs : Array[Double] = Array::make(n, 0.0)
  let fold_changes : Array[Double] = Array::make(n, 0.0)
  // Example values for illustration.
  pvalues[0] = 0.001
  pvalues[1] = 0.45
  pvalues[2] = 0.0001
  pvalues[3] = 0.78
  pvalues[4] = 0.012
  fdrs[0] = 0.005
  fdrs[1] = 0.6
  fdrs[2] = 0.0005
  fdrs[3] = 0.85
  fdrs[4] = 0.06
  fold_changes[0] = 2.5
  fold_changes[1] = -0.3
  fold_changes[2] = 3.1
  fold_changes[3] = 0.1
  fold_changes[4] = 1.8
  let mut n_sig = 0
  let mut i = 0
  while i < n {
    if fdrs[i] < 0.05 {
      n_sig = n_sig + 1
    }
    i = i + 1
  }
  DBAResult::{
    peaks,
    pvalues,
    fdrs,
    fold_changes,
    contrast: ("ConditionA", "ConditionB"),
    n_significant: n_sig,
  }
}
