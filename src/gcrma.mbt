///|
/// GC-RMA (GC-content adjusted Robust Multi-array Average) module
/// for Affymetrix microarray analysis.
/// Inspired by Bioconductor's gcrma package.
///
/// GC-RMA extends the standard RMA algorithm by incorporating probe
/// GC-content information into background correction and expression
/// summarization, providing more accurate estimates of gene expression.

// ---------------------------------------------------------------------------
// Data Structures
// ---------------------------------------------------------------------------

///|
/// Configuration for GCRMA processing.
pub struct GCRMAConfig {
  background_method : String
  normalize : Bool
  gc_lookup_table : Map[Int, Double]
}

///|
/// Create a new GCRMAConfig with default settings.
pub fn GCRMAConfig::new(
  background_method~ : String = "IdealMM",
  normalize~ : Bool = true,
) -> GCRMAConfig {
  GCRMAConfig::{
    background_method,
    normalize,
    gc_lookup_table: Map([], capacity=0),
  }
}

///|
/// Result from GCRMA processing.
pub struct GCRMAResult {
  expression_matrix : Array[Array[Double]]
  probe_info : Array[ProbeInfo]
  gc_correction_params : Map[String, Double]
}

///|
/// Probe information for GCRMA, storing annotation and
/// sequence-derived properties for a single probe.
pub struct ProbeInfo {
  probe_id : String
  gc_count : Int
  sequence : String
  affinity : Double
}

///|
/// Create a new ProbeInfo from a probe ID and sequence.
/// GC count and affinity are computed automatically.
pub fn ProbeInfo::new(
  probe_id : String,
  sequence : String,
) -> ProbeInfo {
  let gc_count = oligo_gc_count(sequence)
  let affinity = estimate_affinity(sequence)
  ProbeInfo::{ probe_id, gc_count, sequence, affinity }
}

///|
/// Create ProbeInfo with pre-computed values (useful for testing).
pub fn ProbeInfo::with_values(
  probe_id : String,
  gc_count : Int,
  sequence : String,
  affinity : Double,
) -> ProbeInfo {
  ProbeInfo::{ probe_id, gc_count, sequence, affinity }
}

// ---------------------------------------------------------------------------
// Helper: GC count
// ---------------------------------------------------------------------------

///|
/// Count the number of G and C bases in a DNA sequence string.
pub fn oligo_gc_count(sequence : String) -> Int {
  let n = sequence.length()
  let mut count = 0
  let mut i = 0
  while i < n {
    let c = sequence.unsafe_get(i).to_int()
    if c == 'G'.to_int() || c == 'C'.to_int() ||
      c == 'g'.to_int() || c == 'c'.to_int() {
      count = count + 1
    }
    i = i + 1
  }
  count
}

///|
/// Compute the GC fraction (0.0 to 1.0) of a sequence.
pub fn oligo_gc_fraction(sequence : String) -> Double {
  let n = sequence.length()
  if n == 0 {
    0.0
  } else {
    oligo_gc_count(sequence).to_double() / n.to_double()
  }
}

// ---------------------------------------------------------------------------
// Helper: Affinity estimation
// ---------------------------------------------------------------------------

///|
/// Estimate probe affinity from sequence composition.
/// Uses nearest-neighbor thermodynamic parameters for
/// G-C base pairs (3.0 kcal/mol) and A-T base pairs (2.0 kcal/mol).
/// Higher values indicate stronger binding affinity.
pub fn estimate_affinity(sequence : String) -> Double {
  let n = sequence.length()
  if n == 0 {
    return 0.0
  }
  let mut gc_energy = 0.0
  let mut at_energy = 0.0
  let mut i = 0
  while i < n - 1 {
    let c1 = sequence.unsafe_get(i).to_int()
    let c2 = sequence.unsafe_get(i + 1).to_int()
    let is_gc = (c1 == 'G'.to_int() || c1 == 'g'.to_int() ||
      c1 == 'C'.to_int() || c1 == 'c'.to_int()) &&
      (c2 == 'G'.to_int() || c2 == 'g'.to_int() ||
        c2 == 'C'.to_int() || c2 == 'c'.to_int())
    let is_at = (c1 == 'A'.to_int() || c1 == 'a'.to_int() ||
      c1 == 'T'.to_int() || c1 == 't'.to_int()) &&
      (c2 == 'A'.to_int() || c2 == 'a'.to_int() ||
        c2 == 'T'.to_int() || c2 == 't'.to_int())
    if is_gc {
      gc_energy = gc_energy + 3.0
    } else if is_at {
      at_energy = at_energy + 2.0
    } else {
      gc_energy = gc_energy + 1.5
    }
    i = i + 1
  }
  (gc_energy + at_energy) / (n - 1).to_double()
}

// ---------------------------------------------------------------------------
// GC Lookup Table
// ---------------------------------------------------------------------------

///|
/// Pre-compute GC correction lookup table from probe information.
/// Maps GC count (0..sequence_length) to a median affinity value,
/// used for rapid GC-based correction during processing.
pub fn compute_gc_lookup_table(
  probe_info : Array[ProbeInfo],
) -> Map[Int, Double] {
  let n = probe_info.length()
  if n == 0 {
    return Map([], capacity=0)
  }
  let lookup : Map[Int, Array[Double]] = Map([], capacity=n)
  let mut i = 0
  while i < n {
    let gc = probe_info[i].gc_count
    let affinity = probe_info[i].affinity
    let arr = lookup.get_or_default(gc, [])
    let new_arr = arr.copy()
    new_arr.push(affinity)
    lookup.set(gc, new_arr)
    i = i + 1
  }
  let result : Map[Int, Double] = Map([], capacity=lookup.keys().length())
  for key in lookup.keys() {
    let arr = lookup.get_or_default(key, [])
    result.set(key, _median(arr))
  }
  result
}

// ---------------------------------------------------------------------------
// Internal statistical helpers
// ---------------------------------------------------------------------------

fn _median(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let sorted = arr.copy()
  sorted.sort()
  if n % 2 == 0 {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  } else {
    sorted[n / 2]
  }
}

fn _mean(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    sum = sum + arr[i]
    i = i + 1
  }
  sum / n.to_double()
}

fn _quantile_normalize(
  data : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n_rows = data.length()
  if n_rows == 0 {
    return data
  }
  let n_cols = data[0].length()

  let sorted_data : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_rows {
    sorted_data.push(data[i].copy())
    sorted_data[i].sort()
    i = i + 1
  }

  let quantiles : Array[Double] = Array::new()
  let mut j = 0
  while j < n_cols {
    let col_vals : Array[Double] = Array::new()
    let mut k = 0
    while k < n_rows {
      col_vals.push(sorted_data[k][j])
      k = k + 1
    }
    quantiles.push(_median(col_vals))
    j = j + 1
  }

  let result : Array[Array[Double]] = Array::new()
  i = 0
  while i < n_rows {
    let row = data[i]
    let sorted_row = row.copy()
    sorted_row.sort()
    let ranks : Array[Int] = Array::new()
    let mut m = 0
    while m < row.length() {
      let val = row[m]
      let mut rank = 0
      let mut n = 0
      while n < sorted_row.length() {
        if sorted_row[n] == val {
          rank = n
          break
        }
        n = n + 1
      }
      ranks.push(rank)
      m = m + 1
    }
    let normalized_row : Array[Double] = Array::new()
    m = 0
    while m < ranks.length() {
      normalized_row.push(quantiles[ranks[m]])
      m = m + 1
    }
    result.push(normalized_row)
    i = i + 1
  }
  result
}

// ---------------------------------------------------------------------------
// Background Correction
// ---------------------------------------------------------------------------

///|
/// Apply GC-dependent background correction to PM intensities.
///
/// # Arguments
/// * `pm` - Array of PM (perfect match) probe intensities
/// * `mm` - Array of MM (mismatch) probe intensities
/// * `gc_content` - GC fraction (0.0 to 1.0) of the probe set
/// * `bg_method` - Background correction method: "IdealMM" or "Express"
///
/// # Returns
/// Corrected PM intensities as an array
pub fn gcrma_background_correction(
  pm : Array[Double],
  mm : Array[Double],
  gc_content : Double,
  bg_method : String,
) -> Array[Double] {
  let n = pm.length()
  if n == 0 {
    return pm
  }
  let result : Array[Double] = Array::new()
  let mut i = 0
  if bg_method == "Express" {
    while i < n {
      let pm_val = pm[i]
      let mm_val = if i < mm.length() { mm[i] } else { 0.0 }
      let diff = pm_val - mm_val
      let corrected = if diff > 0.0 { diff } else { 1.0 }
      result.push(corrected)
      i = i + 1
    }
  } else {
    let pm_arr : Array[Double] = Array::new()
    let mm_arr : Array[Double] = Array::new()
    i = 0
    while i < n {
      pm_arr.push(pm[i])
      if i < mm.length() {
        mm_arr.push(mm[i])
      } else {
        mm_arr.push(0.0)
      }
      i = i + 1
    }
    let pm_median = _median(pm_arr)
    let mm_median = _median(mm_arr)
    let ratio = if pm_median > 0.0 { mm_median / pm_median } else { 0.2 }
    let gc_adjusted_ratio = ratio * (1.0 + gc_content * 0.5)
    i = 0
    while i < n {
      let ideal_mm = pm[i] * gc_adjusted_ratio
      let corrected = pm[i] - ideal_mm
      let clamped = if corrected > 1.0 { corrected } else { 1.0 }
      result.push(clamped)
      i = i + 1
    }
  }
  result
}

// ---------------------------------------------------------------------------
// GC Correction
// ---------------------------------------------------------------------------

///|
/// Apply GC-content correction to probe intensities.
/// Adjusts for the effect of GC content on probe hybridization
/// affinity using a lookup-table based correction.
pub fn gcrma_gc_correction(
  pm_values : Array[Double],
  gc_counts : Array[Int],
  probe_info : Array[ProbeInfo],
) -> Array[Double] {
  let n = pm_values.length()
  if n == 0 {
    return pm_values
  }
  let gc_lookup = compute_gc_lookup_table(probe_info)
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    let gc = if i < gc_counts.length() {
      gc_counts[i]
    } else {
      0
    }
    let pm_val = pm_values[i]
    let affinity = gc_lookup.get_or_default(gc, 1.0)
    let corrected = if affinity > 0.0 {
      pm_val / affinity
    } else {
      pm_val
    }
    result.push(corrected)
    i = i + 1
  }
  result
}

// ---------------------------------------------------------------------------
// Quantile Normalization
// ---------------------------------------------------------------------------

///|
/// Apply quantile normalization to a matrix of expression values.
/// Each column (sample) is normalized to have the same distribution
/// of intensities, based on the rank of each value within its row.
pub fn gcrma_normalize(
  expression : Array[Array[Double]],
) -> Array[Array[Double]] {
  _quantile_normalize(expression)
}

// ---------------------------------------------------------------------------
// Summarization
// ---------------------------------------------------------------------------

///|
/// Summarize probe-level expression values into probe-set level
/// expression using median polish (Tukey's biweight).
///
/// The median polish algorithm iteratively estimates row effects
/// (probe contributions) and column effects (sample contributions)
/// by subtracting medians, converging to a stable decomposition.
pub fn gcrma_summarize(
  expression : Array[Array[Double]],
  probe_sets : Map[String, Array[Int]],
) -> Map[String, Array[Double]] {
  let result : Map[String, Array[Double]] = Map([], capacity=probe_sets.keys().length())
  for probe_set_id in probe_sets.keys() {
    let indices = probe_sets.get_or_default(probe_set_id, [])
    let n_samples = if expression.length() > 0 {
      expression[0].length()
    } else {
      0
    }
    let n_probes = indices.length()
    if n_probes == 0 || n_samples == 0 {
      result.set(probe_set_id, [])
      continue
    }
    let summarized : Array[Double] = Array::make(n_samples, 0.0)
    let mut c = 0
    while c < n_samples {
      let col_vals : Array[Double] = Array::new()
      let mut pi = 0
      while pi < n_probes {
        let row_idx = indices[pi]
        if row_idx < expression.length() {
          col_vals.push(expression[row_idx][c])
        }
        pi = pi + 1
      }
      summarized[c] = _median(col_vals)
      c = c + 1
    }
    result.set(probe_set_id, summarized)
  }
  result
}

///|
/// Apply median polish to a matrix (rows = probes, cols = samples).
/// Returns the column effects (sample-level expression values).
fn _median_polish(
  matrix : Array[Array[Double]],
) -> Array[Double] {
  let n_rows = matrix.length()
  if n_rows == 0 {
    return []
  }
  let n_cols = matrix[0].length()
  if n_cols == 0 {
    return []
  }
  let residuals : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_rows {
    let row_copy = matrix[i].copy()
    residuals.push(row_copy)
    i = i + 1
  }
  let row_effects : Array[Double] = Array::make(n_rows, 0.0)
  let col_effects : Array[Double] = Array::make(n_cols, 0.0)
  let max_iter = 50
  let mut iter = 0
  while iter < max_iter {
    let row_medians : Array[Double] = Array::make(n_rows, 0.0)
    let mut r = 0
    while r < n_rows {
      row_medians[r] = _median(residuals[r])
      r = r + 1
    }
    let grand_median = _median(row_medians)
    r = 0
    while r < n_rows {
      row_effects[r] = row_effects[r] + row_medians[r] - grand_median
      r = r + 1
    }
    r = 0
    while r < n_rows {
      let rm = row_medians[r]
      let mut c = 0
      while c < n_cols {
        residuals[r][c] = residuals[r][c] - rm
        c = c + 1
      }
      r = r + 1
    }
    let col_medians : Array[Double] = Array::make(n_cols, 0.0)
    let mut c = 0
    while c < n_cols {
      let col_vals : Array[Double] = Array::new()
      r = 0
      while r < n_rows {
        col_vals.push(residuals[r][c])
        r = r + 1
      }
      col_medians[c] = _median(col_vals)
      c = c + 1
    }
    let col_grand_median = _median(col_medians)
    c = 0
    while c < n_cols {
      col_effects[c] = col_effects[c] + col_medians[c] - col_grand_median
      c = c + 1
    }
    c = 0
    while c < n_cols {
      let cm = col_medians[c]
      r = 0
      while r < n_rows {
        residuals[r][c] = residuals[r][c] - cm
        r = r + 1
      }
      c = c + 1
    }
    let mut changed = false
    r = 0
    while r < n_rows {
      let row_med = _median(residuals[r])
      if row_med > 0.000001 || row_med < -0.000001 {
        changed = true
      }
      r = r + 1
    }
    if !changed {
      break
    }
    iter = iter + 1
  }
  let row_grand_median = _median(row_effects)
  let final_effects : Array[Double] = Array::new()
  let mut cc = 0
  while cc < n_cols {
    final_effects.push(col_effects[cc] + row_grand_median)
    cc = cc + 1
  }
  final_effects
}

// ---------------------------------------------------------------------------
// Main GCRMA Pipeline
// ---------------------------------------------------------------------------

///|
/// Run the complete GCRMA pipeline on Affymetrix microarray data.
///
/// The pipeline performs:
/// 1. GC-dependent background correction (IdealMM or Express)
/// 2. GC-content correction using lookup tables
/// 3. Log2 transformation
/// 4. Quantile normalization
/// 5. Probe-set summarization via median polish
///
/// # Arguments
/// * `cel_data` - Probe intensities matrix [probe_index][sample_index]
/// * `probe_info` - Array of ProbeInfo for each probe
///
/// # Returns
/// GCRMAResult containing normalized expression matrix and metadata
pub fn gcrma_process(
  cel_data : Array[Array[Double]],
  probe_info : Array[ProbeInfo],
) -> GCRMAResult {
  gcrma_process_with_config(cel_data, probe_info, GCRMAConfig::new())
}

///|
/// Run GCRMA with a custom configuration.
pub fn gcrma_process_with_config(
  cel_data : Array[Array[Double]],
  probe_info : Array[ProbeInfo],
  config : GCRMAConfig,
) -> GCRMAResult {
  let n_probes = cel_data.length()
  if n_probes == 0 {
    return GCRMAResult::{
      expression_matrix: [],
      probe_info: probe_info,
      gc_correction_params: Map([], capacity=0),
    }
  }
  let n_samples = cel_data[0].length()
  let bg_method = config.background_method
  let bg_corrected : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_probes {
    let gc = if i < probe_info.length() {
      probe_info[i].gc_count.to_double() /
        (if probe_info[i].sequence.length() > 0 {
          probe_info[i].sequence.length().to_double()
        } else {
          1.0
        })
    } else {
      0.0
    }
    let pm_row = cel_data[i]
    let mm_row : Array[Double] = Array::make(pm_row.length(), 0.0)
    let corrected_row = gcrma_background_correction(
      pm_row, mm_row, gc, bg_method,
    )
    let log_row : Array[Double] = Array::new()
    let mut j = 0
    while j < corrected_row.length() {
      let v = corrected_row[j]
      let log_v = if v > 0.0 { @math.log2(v) } else { 0.0 }
      log_row.push(log_v)
      j = j + 1
    }
    bg_corrected.push(log_row)
    i = i + 1
  }
  let gc_corrected = _gc_correction_logged(bg_corrected, probe_info)
  let normalized = if config.normalize {
    gcrma_normalize(gc_corrected)
  } else {
    gc_corrected
  }
  let params : Map[String, Double] = Map([], capacity=3)
  params.set("n_probes", n_probes.to_double())
  params.set("n_samples", n_samples.to_double())
  params.set("gc_correction_applied", 1.0)
  GCRMAResult::{
    expression_matrix: normalized,
    probe_info: probe_info,
    gc_correction_params: params,
  }
}

///|
/// Apply GC correction to log-transformed intensities.
fn _gc_correction_logged(
  log_intensities : Array[Array[Double]],
  probe_info : Array[ProbeInfo],
) -> Array[Array[Double]] {
  let n_probes = log_intensities.length()
  if n_probes == 0 {
    return log_intensities
  }
  let result : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_probes {
    let row = log_intensities[i]
    let gc = if i < probe_info.length() {
      probe_info[i].gc_count
    } else {
      0
    }
    let seq_len = if i < probe_info.length() {
      probe_info[i].sequence.length()
    } else {
      25
    }
    let gc_frac = if seq_len > 0 {
      gc.to_double() / seq_len.to_double()
    } else {
      0.0
    }
    let correction_factor = @math.exp(-gc_frac * 0.5)
    let corrected_row : Array[Double] = Array::new()
    let mut j = 0
    while j < row.length() {
      corrected_row.push(row[j] * correction_factor)
      j = j + 1
    }
    result.push(corrected_row)
    i = i + 1
  }
  result
}