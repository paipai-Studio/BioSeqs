///|
/// Bioconductor EnrichedHeatmap - normalize genomic signals to target regions.
///
/// Mirrors the `normalizeToMatrix` function from the Bioconductor
/// `EnrichedHeatmap` package (Gu, Eils & Schlesner, 2018, BMC Genomics
/// 19:234). Given genomic signal regions (e.g. CpG methylation, histone
/// marks) and target regions (e.g. TSS, promoters), this module normalizes
/// the signal into a matrix with equal number of columns per target for
/// visualization as an "enriched heatmap."
///
/// The normalization process:
///   1. Extend each target upstream and downstream by `extend` bp.
///   2. Split the flanking regions into windows of `w` bp.
///   3. For each window, compute the mean signal of overlapping signal
///      regions (using one of several mean modes).
///   4. Optionally smooth each row and clip extreme values.
///
/// Mean modes:
///   * `absolute` — unweighted mean of overlapping signal values.
///   * `weighted` — overlap-width-weighted mean.
///   * `w0` — like weighted, but denominator includes non-covered bp.
///   * `coverage` — weighted sum divided by total window width.

// ============================================================================
// Data structures
// ============================================================================

///|
/// A genomic interval with an associated signal value.
pub struct GenomicSignal {
  chr : String
  start : Int // 1-based start
  end_ : Int // 1-based end (inclusive)
  value : Double
}

///|
pub fn GenomicSignal::new(
  chr : String,
  start : Int,
  end_ : Int,
  value : Double,
) -> GenomicSignal {
  GenomicSignal::{ chr, start, end_, value }
}

///|
/// A target genomic region (e.g. TSS, promoter).
pub struct TargetRegion {
  chr : String
  start : Int // 1-based start
  end_ : Int // 1-based end (inclusive)
  strand : String // "+" or "-"
}

///|
pub fn TargetRegion::new(
  chr : String,
  start : Int,
  end_ : Int,
  strand : String,
) -> TargetRegion {
  TargetRegion::{ chr, start, end_, strand }
}

///|
/// Mean computation mode for signal normalization.
pub enum MeanMode {
  Absolute
  Weighted
  W0
  Coverage
}

///|
pub fn mean_mode_absolute() -> MeanMode {
  MeanMode::Absolute
}

///|
pub fn mean_mode_weighted() -> MeanMode {
  MeanMode::Weighted
}

///|
pub fn mean_mode_w0() -> MeanMode {
  MeanMode::W0
}

///|
pub fn mean_mode_coverage() -> MeanMode {
  MeanMode::Coverage
}

///|
/// Configuration for signal normalization.
pub struct EnrichedHeatmapConfig {
  /// Upstream extension in bp.
  extendUp : Int
  /// Downstream extension in bp.
  extendDown : Int
  /// Window size in bp for splitting flanking regions.
  w : Int
  /// Number of windows for the target body (0 = no target body).
  k : Int
  /// Mean computation mode.
  meanMode : MeanMode
  /// Whether to include target body columns.
  includeTarget : Bool
  /// Whether to apply row-wise smoothing.
  smooth : Bool
  /// Smoothing window size (number of columns).
  smoothWindow : Int
  /// Background value for windows with no signal overlap.
  background : Double
  /// Whether to flip upstream to downstream for single-point targets.
  flipUpstream : Bool
  /// Lower percentile for clipping (0.0 = no clipping).
  keepLow : Double
  /// Upper percentile for clipping (1.0 = no clipping).
  keepHigh : Double
}

///|
pub fn EnrichedHeatmapConfig::new(
  extendUp? : Int = 5000,
  extendDown? : Int = 5000,
  w? : Int = 100,
  k? : Int = 0,
  meanMode? : MeanMode = Absolute,
  includeTarget? : Bool = false,
  smooth? : Bool = false,
  smoothWindow? : Int = 5,
  background? : Double = 0.0,
  flipUpstream? : Bool = false,
  keepLow? : Double = 0.0,
  keepHigh? : Double = 1.0,
) -> EnrichedHeatmapConfig {
  EnrichedHeatmapConfig::{
    extendUp,
    extendDown,
    w,
    k,
    meanMode,
    includeTarget,
    smooth,
    smoothWindow,
    background,
    flipUpstream,
    keepLow,
    keepHigh,
  }
}

///|
/// Result of signal normalization: a matrix and metadata.
pub struct NormalizedMatrix {
  /// n_targets × n_cols signal matrix.
  matrix : Array[Array[Double]]
  /// Column indices for upstream region.
  upstreamIndex : Array[Int]
  /// Column indices for target body.
  targetIndex : Array[Int]
  /// Column indices for downstream region.
  downstreamIndex : Array[Int]
  /// Extension in bp [up, down].
  extend : (Int, Int)
  /// Whether smoothing was applied.
  smooth : Bool
  /// Background value used.
  background : Double
  /// Number of rows (targets).
  nRows : Int
  /// Number of columns (windows).
  nCols : Int
  /// Signal name (for labeling).
  signalName : String
  /// Target name (for labeling).
  targetName : String
  /// Indices of rows where smoothing failed.
  failedRows : Array[Int]
}

// ============================================================================
// Window creation
// ============================================================================

///|
/// A genomic window with a reference to its parent target.
struct Window {
  chr : String
  start : Int
  end_ : Int
  targetIndex : Int // index of the parent target region
  windowIndex : Int // window index within the target (0-based)
  isTarget : Bool // true for target body windows, false for flanking
}

///|
/// Create windows for upstream, target body, and downstream of each target.
fn make_windows(
  targets : Array[TargetRegion],
  config : EnrichedHeatmapConfig,
) -> Array[Window] {
  let windows : Array[Window] = []
  let n_up = config.extendUp / config.w
  let n_down = config.extendDown / config.w
  let mut ti = 0
  while ti < targets.length() {
    let t = targets[ti]
    let is_minus = t.strand == "-"
    let tStart = t.start
    let tEnd = t.end_
    // Upstream windows.
    let mut wi = 0
    while wi < n_up {
      let w_start = if is_minus {
        tEnd + 1 + wi * config.w
      } else {
        tStart - config.extendUp + wi * config.w
      }
      let w_end = w_start + config.w - 1
      // For minus strand, reverse the window index.
      let actual_wi = if is_minus { n_up - 1 - wi } else { wi }
      windows.push(
        Window::{
          chr: t.chr,
          start: w_start,
          end_: w_end,
          targetIndex: ti,
          windowIndex: actual_wi,
          isTarget: false,
        },
      )
      wi = wi + 1
    }
    // Target body windows (if includeTarget and k > 0).
    if config.includeTarget && config.k > 0 && tEnd >= tStart {
      let target_width = tEnd - tStart + 1
      let kw = if target_width >= config.k {
        target_width / config.k
      } else {
        1
      }
      let mut ki = 0
      while ki < config.k {
        let w_start = tStart + ki * kw
        let w_end = if ki == config.k - 1 {
          tEnd
        } else {
          w_start + kw - 1
        }
        let actual_ki = if is_minus { config.k - 1 - ki } else { ki }
        windows.push(
          Window::{
            chr: t.chr,
            start: w_start,
            end_: w_end,
            targetIndex: ti,
            windowIndex: n_up + actual_ki,
            isTarget: true,
          },
        )
        ki = ki + 1
      }
    }
    // Downstream windows.
    let mut di = 0
    while di < n_down {
      let w_start = if is_minus {
        tStart - config.extendDown + di * config.w
      } else {
        tEnd + 1 + di * config.w
      }
      let w_end = w_start + config.w - 1
      let actual_di = if is_minus { n_down - 1 - di } else { di }
      let base_idx = n_up + (if config.includeTarget && config.k > 0 { config.k } else { 0 })
      windows.push(
        Window::{
          chr: t.chr,
          start: w_start,
          end_: w_end,
          targetIndex: ti,
          windowIndex: base_idx + actual_di,
          isTarget: false,
        },
      )
      di = di + 1
    }
    ti = ti + 1
  }
  windows
}

// ============================================================================
// Overlap detection and mean computation
// ============================================================================

///|
/// Check if two intervals overlap (1-based, inclusive).
fn overlaps(s1 : Int, e1 : Int, s2 : Int, e2 : Int) -> Bool {
  s1 <= e2 && s2 <= e1
}

///|
/// Compute overlap width (1-based inclusive).
fn overlap_width(s1 : Int, e1 : Int, s2 : Int, e2 : Int) -> Int {
  let lo = if s1 > s2 { s1 } else { s2 }
  let hi = if e1 < e2 { e1 } else { e2 }
  if hi >= lo {
    hi - lo + 1
  } else {
    0
  }
}

///|
/// Compute the mean signal for a single window based on the mean mode.
fn compute_window_mean(
  window : Window,
  signals : Array[GenomicSignal],
  mode : MeanMode,
  background : Double,
) -> Double {
  let w_width = window.end_ - window.start + 1
  if w_width <= 0 {
    return background
  }
  // Find overlapping signals.
  let values : Array[Double] = []
  let widths : Array[Int] = []
  let mut total_overlap = 0
  let mut si = 0
  while si < signals.length() {
    let s = signals[si]
    if s.chr == window.chr && overlaps(s.start, s.end_, window.start, window.end_) {
      let ow = overlap_width(s.start, s.end_, window.start, window.end_)
      values.push(s.value)
      widths.push(ow)
      total_overlap = total_overlap + ow
    }
    si = si + 1
  }
  if values.length() == 0 {
    return background
  }
  match mode {
    Absolute => {
      // Unweighted mean of overlapping values.
      let mut sum = 0.0
      let mut i = 0
      while i < values.length() {
        sum = sum + values[i]
        i = i + 1
      }
      sum / values.length().to_double()
    }
    Weighted => {
      // Overlap-width-weighted mean.
      let mut wsum = 0.0
      let mut wtot = 0.0
      let mut i = 0
      while i < values.length() {
        let w = widths[i].to_double()
        wsum = wsum + values[i] * w
        wtot = wtot + w
        i = i + 1
      }
      if wtot > 0.0 { wsum / wtot } else { background }
    }
    W0 => {
      // Like weighted, but denominator includes non-covered bp.
      let mut wsum = 0.0
      let mut wtot = 0.0
      let mut i = 0
      while i < values.length() {
        let w = widths[i].to_double()
        wsum = wsum + values[i] * w
        wtot = wtot + w
        i = i + 1
      }
      let non_overlap = w_width.to_double() - total_overlap.to_double()
      let denom = wtot + non_overlap
      if denom > 0.0 { wsum / denom } else { background }
    }
    Coverage => {
      // Weighted sum divided by total window width.
      let mut wsum = 0.0
      let mut i = 0
      while i < values.length() {
        let w = widths[i].to_double()
        wsum = wsum + values[i] * w
        i = i + 1
      }
      let denom = w_width.to_double()
      if denom > 0.0 { wsum / denom } else { background }
    }
  }
}

// ============================================================================
// Smoothing
// ============================================================================

///|
/// Apply moving-average smoothing to a row of the matrix.
fn smooth_row(row : Array[Double], window : Int) -> Array[Double] {
  let n = row.length()
  if n == 0 || window <= 1 {
    return row
  }
  let half = window / 2
  let result : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let lo = if i - half > 0 { i - half } else { 0 }
    let hi = if i + half < n - 1 { i + half } else { n - 1 }
    let mut sum = 0.0
    let mut count = 0
    let mut j = lo
    while j <= hi {
      sum = sum + row[j]
      count = count + 1
      j = j + 1
    }
    result[i] = if count > 0 { sum / count.to_double() } else { row[i] }
    i = i + 1
  }
  result
}

// ============================================================================
// Percentile clipping
// ============================================================================

///|
/// Compute the p-th percentile of an array (0 <= p <= 1).
fn percentile(arr : Array[Double], p : Double) -> Double {
  let sorted = eh_sort_doubles(arr)
  let n = sorted.length()
  if n == 0 {
    return 0.0
  }
  let idx = p * (n - 1).to_double()
  let lo = idx.floor().to_int()
  let hi = idx.ceil().to_int()
  if lo == hi {
    return sorted[lo]
  }
  let frac = idx - lo.to_double()
  sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

///|
/// Sort doubles ascending (returns new array).
fn eh_sort_doubles(arr : Array[Double]) -> Array[Double] {
  let result : Array[Double] = []
  let mut i = 0
  while i < arr.length() {
    result.push(arr[i])
    i = i + 1
  }
  let n = result.length()
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n - 1 - i {
      if result[j] > result[j + 1] {
        let tmp = result[j]
        result[j] = result[j + 1]
        result[j + 1] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Clip matrix values to [q_low, q_high] percentiles.
fn clip_matrix(
  matrix : Array[Array[Double]],
  qLow : Double,
  qHigh : Double,
) -> Array[Array[Double]] {
  if qLow <= 0.0 && qHigh >= 1.0 {
    return matrix
  }
  // Collect all values.
  let all_vals : Array[Double] = []
  let mut i = 0
  while i < matrix.length() {
    let mut j = 0
    while j < matrix[i].length() {
      all_vals.push(matrix[i][j])
      j = j + 1
    }
    i = i + 1
  }
  let lo = percentile(all_vals, qLow)
  let hi = percentile(all_vals, qHigh)
  let result : Array[Array[Double]] = []
  let mut i = 0
  while i < matrix.length() {
    let row : Array[Double] = []
    let mut j = 0
    while j < matrix[i].length() {
      let v = matrix[i][j]
      if v < lo {
        row.push(lo)
      } else if v > hi {
        row.push(hi)
      } else {
        row.push(v)
      }
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  result
}

// ============================================================================
// Main normalization function
// ============================================================================

///|
/// Normalize genomic signals to target regions.
///
/// Parameters:
///   - `signals`: array of GenomicSignal (signal regions with values).
///   - `targets`: array of TargetRegion (target regions to normalize to).
///   - `config`: normalization configuration.
pub fn normalize_to_matrix(
  signals : Array[GenomicSignal],
  targets : Array[TargetRegion],
  config : EnrichedHeatmapConfig,
) -> NormalizedMatrix {
  let n_targets = targets.length()
  if n_targets == 0 {
    return NormalizedMatrix::{
      matrix: [],
      upstreamIndex: [],
      targetIndex: [],
      downstreamIndex: [],
      extend: (config.extendUp, config.extendDown),
      smooth: config.smooth,
      background: config.background,
      nRows: 0,
      nCols: 0,
      signalName: "",
      targetName: "",
      failedRows: [],
    }
  }

  // Create windows for all targets.
  let windows = make_windows(targets, config)
  let n_up = config.extendUp / config.w
  let n_down = config.extendDown / config.w
  let n_target = if config.includeTarget && config.k > 0 { config.k } else { 0 }
  let n_cols = n_up + n_target + n_down

  // Compute column indices.
  let upstreamIndex : Array[Int] = []
  let targetIndex : Array[Int] = []
  let downstreamIndex : Array[Int] = []
  let mut ci = 0
  while ci < n_up {
    upstreamIndex.push(ci)
    ci = ci + 1
  }
  let mut ci = 0
  while ci < n_target {
    targetIndex.push(n_up + ci)
    ci = ci + 1
  }
  let mut ci = 0
  while ci < n_down {
    downstreamIndex.push(n_up + n_target + ci)
    ci = ci + 1
  }

  // Build the matrix: for each target, compute mean signal per window.
  let matrix : Array[Array[Double]] = Array::make(n_targets, [])
  let mut ri = 0
  while ri < n_targets {
    matrix[ri] = Array::make(n_cols, config.background)
    ri = ri + 1
  }

  let mut wi = 0
  while wi < windows.length() {
    let win = windows[wi]
    let mean_val = compute_window_mean(win, signals, config.meanMode, config.background)
    matrix[win.targetIndex][win.windowIndex] = mean_val
    wi = wi + 1
  }

  // Apply smoothing if requested.
  let failedRows : Array[Int] = []
  let final_matrix = if config.smooth {
    let smoothed : Array[Array[Double]] = []
    let mut ri = 0
    while ri < n_targets {
      smoothed.push(smooth_row(matrix[ri], config.smoothWindow))
      ri = ri + 1
    }
    smoothed
  } else {
    matrix
  }

  // Apply percentile clipping.
  let clipped = clip_matrix(final_matrix, config.keepLow, config.keepHigh)

  NormalizedMatrix::{
    matrix: clipped,
    upstreamIndex,
    targetIndex,
    downstreamIndex,
    extend: (config.extendUp, config.extendDown),
    smooth: config.smooth,
    background: config.background,
    nRows: n_targets,
    nCols: n_cols,
    signalName: "",
    targetName: "",
    failedRows,
  }
}

// ============================================================================
// Enrichment profile computation
// ============================================================================

///|
/// Compute the mean enrichment profile (column means of the normalized matrix).
pub fn enrichment_profile(mat : NormalizedMatrix) -> Array[Double] {
  let n_cols = mat.nCols
  let n_rows = mat.nRows
  if n_cols == 0 || n_rows == 0 {
    return []
  }
  let profile : Array[Double] = Array::make(n_cols, 0.0)
  let mut ci = 0
  while ci < n_cols {
    let mut sum = 0.0
    let mut count = 0
    let mut ri = 0
    while ri < n_rows {
      sum = sum + mat.matrix[ri][ci]
      count = count + 1
      ri = ri + 1
    }
    profile[ci] = if count > 0 { sum / count.to_double() } else { 0.0 }
    ci = ci + 1
  }
  profile
}

///|
/// Compute row means of the normalized matrix (per-target average signal).
pub fn eh_row_means(mat : NormalizedMatrix) -> Array[Double] {
  let means : Array[Double] = []
  let mut ri = 0
  while ri < mat.nRows {
    let mut sum = 0.0
    let mut count = 0
    let mut ci = 0
    while ci < mat.nCols {
      sum = sum + mat.matrix[ri][ci]
      count = count + 1
      ci = ci + 1
    }
    means.push(if count > 0 { sum / count.to_double() } else { 0.0 })
    ri = ri + 1
  }
  means
}

// ============================================================================
// Summary and sample data
// ============================================================================

///|
/// Generate a summary string of the normalized matrix.
pub fn enriched_heatmap_summary(mat : NormalizedMatrix) -> String {
  let buf = StringBuilder::new()
  buf.write_string("EnrichedHeatmap Normalized Matrix Summary\n")
  buf.write_string("  Rows (targets): " + mat.nRows.to_string() + "\n")
  buf.write_string("  Columns (windows): " + mat.nCols.to_string() + "\n")
  buf.write_string("  Upstream windows: " + mat.upstreamIndex.length().to_string() + "\n")
  buf.write_string("  Target body windows: " + mat.targetIndex.length().to_string() + "\n")
  buf.write_string("  Downstream windows: " + mat.downstreamIndex.length().to_string() + "\n")
  buf.write_string("  Extension: " + mat.extend.0.to_string() + " up, " + mat.extend.1.to_string() + " down\n")
  buf.write_string("  Smoothed: " + mat.smooth.to_string() + "\n")
  buf.write_string("  Background: " + mat.background.to_string() + "\n")
  buf.write_string("  Failed rows: " + mat.failedRows.length().to_string() + "\n")
  let profile = enrichment_profile(mat)
  if profile.length() > 0 {
    let mut minVal = profile[0]
    let mut maxVal = profile[0]
    let mut sum = 0.0
    let mut i = 0
    while i < profile.length() {
      if profile[i] < minVal {
        minVal = profile[i]
      }
      if profile[i] > maxVal {
        maxVal = profile[i]
      }
      sum = sum + profile[i]
      i = i + 1
    }
    let mean = sum / profile.length().to_double()
    buf.write_string("  Enrichment profile: min=" + minVal.to_string() + " max=" + maxVal.to_string() + " mean=" + mean.to_string() + "\n")
  }
  buf.to_string()
}

///|
/// Generate sample signal and target data for demos.
/// Returns (signals, targets) with 10 targets (TSS points) and
/// 30 signal regions with values around the targets.
pub fn enriched_heatmap_sample() -> (Array[GenomicSignal], Array[TargetRegion]) {
  let signals : Array[GenomicSignal] = []
  let targets : Array[TargetRegion] = []
  // Create 10 target TSS points on chr1.
  let mut i = 0
  while i < 10 {
    let tss = 10000 + i * 5000
    targets.push(TargetRegion::new("chr1", tss, tss, "+"))
    i = i + 1
  }
  // Create signal regions with high values near each TSS.
  let mut i = 0
  while i < 10 {
    let tss = 10000 + i * 5000
    // Upstream signals (increasing towards TSS).
    signals.push(GenomicSignal::new("chr1", tss - 4000, tss - 3901, 0.1))
    signals.push(GenomicSignal::new("chr1", tss - 3000, tss - 2901, 0.3))
    signals.push(GenomicSignal::new("chr1", tss - 2000, tss - 1901, 0.6))
    signals.push(GenomicSignal::new("chr1", tss - 1000, tss - 901, 0.8))
    // At TSS: peak.
    signals.push(GenomicSignal::new("chr1", tss - 200, tss + 200, 1.0))
    // Downstream signals (decreasing from TSS).
    signals.push(GenomicSignal::new("chr1", tss + 901, tss + 1000, 0.7))
    signals.push(GenomicSignal::new("chr1", tss + 1901, tss + 2000, 0.4))
    signals.push(GenomicSignal::new("chr1", tss + 2901, tss + 3000, 0.2))
    i = i + 1
  }
  (signals, targets)
}
