///|
/// MSnbase - Base structures and processing for mass spectrometry data.
/// Inspired by Bioconductor MSnbase package.
///
/// Provides core data structures for mass spectrometry and proteomics:
/// - MSnSet: Quantitative feature x sample data with feature annotations
/// - Spectrum: MS1 and MS2 spectra representation
/// - Chromatogram: Extracted ion chromatograms (XIC)
/// - Processing steps tracking
/// - Spectral processing (smoothing, baseline correction, centroiding)
/// - Peak picking and feature grouping
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/MSnbase.html
///   - Gatto L, Lilley KS (2012) MSnbase - an R/Bioconductor package for
///     isobaric tagged mass spectrometry data visualization, processing and
///     quantitation. Bioinformatics 28:288-289.

// ============================================================================
// Enums and basic types
// ============================================================================

///|
/// Mass spectrum acquisition level (MS1 = full scan, MS2 = MS/MS, MSn = higher).
pub enum MSLevel {
  MS1
  MS2
  MSn(Int)
} derive(Eq, Debug)

///|
/// Create MS1 level.
pub fn ms_level_ms1() -> MSLevel { MSLevel::MS1 }

///|
/// Create MS2 level.
pub fn ms_level_ms2() -> MSLevel { MSLevel::MS2 }

///|
/// Create MSn level with specific n.
pub fn ms_level_msn(n : Int) -> MSLevel { MSLevel::MSn(n) }

///|
/// Ionization mode.
pub enum Polarity {
  Positive
  Negative
  Unknown
} derive(Eq, Debug)

///|
/// Create positive polarity.
pub fn polarity_positive() -> Polarity { Polarity::Positive }

///|
/// Create negative polarity.
pub fn polarity_negative() -> Polarity { Polarity::Negative }

///|
/// Spectrum processing step (for tracking).
pub struct ProcessingStep {
  description : String
  timestamp : String  // simplified timestamp as string
}

///|
/// Create a processing step record.
pub fn ProcessingStep::new(
  description : String,
  timestamp : String,
) -> ProcessingStep {
  ProcessingStep::{ description, timestamp }
}

// ============================================================================
// Spectrum data structure
// ============================================================================

///|
/// A mass spectrum (MS1, MS2, or MSn).
/// Stores m/z values, intensity values, and metadata.
pub struct Spectrum {
  mz : Array[Double]
  intensity : Array[Double]
  ms_level : MSLevel
  polarity : Polarity
  rt : Double              // retention time in seconds
  scan_index : Int
  precursor_mz : Double    // for MS2
  precursor_charge : Int   // for MS2
  tic : Double             // total ion current
  peaks_count : Int
  mut centroided : Bool
  processing_log : Array[ProcessingStep]
}

///|
/// Create a Spectrum from m/z and intensity arrays.
pub fn Spectrum::new(
  mz : Array[Double],
  intensity : Array[Double],
  ms_level : MSLevel,
  polarity : Polarity,
  rt : Double,
) -> Spectrum {
  let n = mz.length()
  if n != intensity.length() {
    // Return empty spectrum if mismatch
    return Spectrum::empty()
  }
  let mut tic = 0.0
  for i in 0..<n {
    tic = tic + intensity[i]
  }
  Spectrum::{
    mz,
    intensity,
    ms_level,
    polarity,
    rt,
    scan_index: 0,
    precursor_mz: 0.0,
    precursor_charge: 0,
    tic,
    peaks_count: n,
    centroided: false,
    processing_log: [],
  }
}

///|
/// Create an empty spectrum.
pub fn Spectrum::empty() -> Spectrum {
  Spectrum::{
    mz: [],
    intensity: [],
    ms_level: MSLevel::MS1,
    polarity: Polarity::Unknown,
    rt: 0.0,
    scan_index: 0,
    precursor_mz: 0.0,
    precursor_charge: 0,
    tic: 0.0,
    peaks_count: 0,
    centroided: false,
    processing_log: [],
  }
}

///|
/// Create a copy with precursor info (for MS2).
pub fn Spectrum::with_precursor(
  self : Spectrum,
  precursor_mz : Double,
  precursor_charge : Int,
) -> Spectrum {
  Spectrum::{
    mz: self.mz,
    intensity: self.intensity,
    ms_level: self.ms_level,
    polarity: self.polarity,
    rt: self.rt,
    scan_index: self.scan_index,
    precursor_mz,
    precursor_charge,
    tic: self.tic,
    peaks_count: self.peaks_count,
    centroided: self.centroided,
    processing_log: self.processing_log,
  }
}

///|
/// Get m/z array.
pub fn Spectrum::mz(self : Spectrum) -> Array[Double] { self.mz }

///|
/// Get intensity array.
pub fn Spectrum::intensity(self : Spectrum) -> Array[Double] { self.intensity }

///|
/// Get number of peaks.
pub fn Spectrum::peaks_count(self : Spectrum) -> Int { self.peaks_count }

///|
/// Get retention time.
pub fn Spectrum::rt(self : Spectrum) -> Double { self.rt }

///|
/// Get total ion current.
pub fn Spectrum::tic(self : Spectrum) -> Double { self.tic }

///|
/// Get the base peak (highest intensity) info: (mz, intensity).
pub fn Spectrum::base_peak(self : Spectrum) -> (Double, Double) {
  if self.peaks_count == 0 { return (0.0, 0.0) }
  let mut max_i = 0
  let mut max_v = self.intensity[0]
  for i in 1..<self.peaks_count {
    if self.intensity[i] > max_v {
      max_v = self.intensity[i]
      max_i = i
    }
  }
  (self.mz[max_i], max_v)
}

///|
/// Find the closest peak to a given m/z within ppm tolerance.
/// Returns index of matched peak or -1.
pub fn Spectrum::find_peak(
  self : Spectrum,
  target_mz : Double,
  ppm~ : Double = 10.0,
) -> Int {
  let n = self.peaks_count
  if n == 0 { return -1 }
  let tolerance = target_mz * ppm / 1000000.0
  let mut best_idx = -1
  let mut best_diff = 1000000.0
  for i in 0..<n {
    let diff = mbn_abs(self.mz[i] - target_mz)
    if diff <= tolerance && diff < best_diff {
      best_diff = diff
      best_idx = i
    }
  }
  best_idx
}

///|
/// Extract a subset of the spectrum within an m/z range.
pub fn Spectrum::filter_mz_range(
  self : Spectrum,
  mz_min : Double,
  mz_max : Double,
) -> Spectrum {
  let new_mz : Array[Double] = []
  let new_int : Array[Double] = []
  for i in 0..<self.peaks_count {
    if self.mz[i] >= mz_min && self.mz[i] <= mz_max {
      new_mz.push(self.mz[i])
      new_int.push(self.intensity[i])
    }
  }
  Spectrum::new(new_mz, new_int, self.ms_level, self.polarity, self.rt)
}

///|
/// Normalize spectrum intensity by TIC (so sum = 1.0).
pub fn Spectrum::normalize_tic(self : Spectrum) -> Spectrum {
  if self.tic <= 0.0 { return self }
  let new_int : Array[Double] = []
  for i in 0..<self.peaks_count {
    new_int.push(self.intensity[i] / self.tic)
  }
  Spectrum::new(self.mz, new_int, self.ms_level, self.polarity, self.rt)
}

// ============================================================================
// Chromatogram data structure
// ============================================================================

///|
/// An extracted ion chromatogram (XIC): intensity vs retention time.
pub struct Chromatogram {
  rt : Array[Double]         // retention times in seconds
  intensity : Array[Double]  // intensities
  mz_target : Double         // target m/z (0 if full range)
  ppm_tolerance : Double
  n_points : Int
  acquisition_mode : String  // "MRM", "SIM", "XIC", etc.
}

///|
/// Create a Chromatogram.
pub fn Chromatogram::new(
  rt : Array[Double],
  intensity : Array[Double],
  mz_target~ : Double = 0.0,
  ppm_tolerance~ : Double = 10.0,
  acquisition_mode~ : String = "XIC",
) -> Chromatogram {
  let n = rt.length()
  Chromatogram::{
    rt,
    intensity,
    mz_target,
    ppm_tolerance,
    n_points: if n == intensity.length() { n } else { 0 },
    acquisition_mode,
  }
}

///|
/// Get n_points.
pub fn Chromatogram::n_points(self : Chromatogram) -> Int { self.n_points }

///|
/// Get total signal (AUC) using trapezoid rule.
pub fn Chromatogram::total_auc(self : Chromatogram) -> Double {
  if self.n_points < 2 { return 0.0 }
  let mut auc = 0.0
  for i in 1..<self.n_points {
    let dt = self.rt[i] - self.rt[i - 1]
    let avg_int = (self.intensity[i] + self.intensity[i - 1]) / 2.0
    auc = auc + dt * avg_int
  }
  auc
}

///|
/// Find apex (maximum intensity) retention time and intensity.
pub fn Chromatogram::apex(self : Chromatogram) -> (Double, Double) {
  if self.n_points == 0 { return (0.0, 0.0) }
  let mut max_idx = 0
  let mut max_v = self.intensity[0]
  for i in 1..<self.n_points {
    if self.intensity[i] > max_v {
      max_v = self.intensity[i]
      max_idx = i
    }
  }
  (self.rt[max_idx], max_v)
}

///|
/// Full width at half maximum (FWHM) in seconds.
pub fn Chromatogram::fwhm(self : Chromatogram) -> Double {
  if self.n_points < 3 { return 0.0 }
  let (_apex_rt, apex_int) = self.apex()
  if apex_int <= 0.0 { return 0.0 }
  let half_max = apex_int / 2.0
  // Find index of apex
  let mut apex_idx = 0
  for i in 0..<self.n_points {
    if self.intensity[i] == apex_int { apex_idx = i; break }
  }
  // Walk left
  let mut left_idx = apex_idx
  while left_idx > 0 && self.intensity[left_idx] > half_max {
    left_idx = left_idx - 1
  }
  // Walk right
  let mut right_idx = apex_idx
  while right_idx < self.n_points - 1 && self.intensity[right_idx] > half_max {
    right_idx = right_idx + 1
  }
  self.rt[right_idx] - self.rt[left_idx]
}

// ============================================================================
// MSnSet - quantitative feature data (core data structure)
// ============================================================================

///|
/// A single feature (row) annotation in MSnSet.
/// Represents a peptide/protein/feature with metadata.
pub struct MSnFeatureData {
  feature_name : String
  protein_accessions : Array[String]
  sequence : String        // peptide sequence if applicable
  modifications : Array[String]
  charge : Int
  mz : Double
  rt : Double
  additional : Map[String, String]
}

///|
/// Create a feature annotation record.
pub fn MSnFeatureData::new(
  feature_name : String,
  protein_accessions : Array[String],
  sequence : String,
  charge : Int,
  mz : Double,
  rt : Double,
) -> MSnFeatureData {
  MSnFeatureData::{
    feature_name,
    protein_accessions,
    sequence,
    modifications: [],
    charge,
    mz,
    rt,
    additional: Map([], capacity=0),
  }
}

///|
/// A single sample (column) annotation.
pub struct MSnSampleData {
  sample_name : String
  group : String
  subject : String
  replicate : Int
  run : String
  acquisition_date : String
  additional : Map[String, String]
}

///|
/// Create a sample annotation record.
pub fn MSnSampleData::new(
  sample_name : String,
  group : String,
  subject : String,
  replicate : Int,
  run : String,
) -> MSnSampleData {
  MSnSampleData::{
    sample_name,
    group,
    subject,
    replicate,
    run,
    acquisition_date: "",
    additional: Map([], capacity=0),
  }
}

///|
/// MSnSet: the main container for quantitative MS-based proteomics data.
/// Stores:
/// - exprs: expression (intensity) matrix, rows=features, cols=samples
/// - feature_data: row annotations
/// - sample_data: column annotations
/// - processing: processing log
pub struct MSnSet {
  exprs : Array[Array[Double]]    // rows=features, cols=samples
  feature_data : Array[MSnFeatureData]
  sample_data : Array[MSnSampleData]
  n_features : Int
  n_samples : Int
  processing_log : Array[ProcessingStep]
  mut is_log2 : Bool
}

///|
/// Create an MSnSet from data matrix, feature annotations, and sample annotations.
pub fn MSnSet::new(
  exprs : Array[Array[Double]],
  feature_data : Array[MSnFeatureData],
  sample_data : Array[MSnSampleData],
) -> MSnSet {
  let n_f = exprs.length()
  let n_s = if n_f > 0 { exprs[0].length() } else { 0 }
  MSnSet::{
    exprs,
    feature_data,
    sample_data,
    n_features: n_f,
    n_samples: n_s,
    processing_log: [],
    is_log2: false,
  }
}

///|
/// Quick create an MSnSet from names and matrix only.
pub fn MSnSet::from_names(
  exprs : Array[Array[Double]],
  feature_names : Array[String],
  sample_names : Array[String],
) -> MSnSet {
  let n_f = exprs.length()
  let n_s = if n_f > 0 { exprs[0].length() } else { 0 }
  let fd : Array[MSnFeatureData] = []
  for i in 0..<n_f {
    let name = if i < feature_names.length() { feature_names[i] } else { "feat_\{i}" }
    fd.push(MSnFeatureData::new(name, [], "", 0, 0.0, 0.0))
  }
  let sd : Array[MSnSampleData] = []
  for j in 0..<n_s {
    let name = if j < sample_names.length() { sample_names[j] } else { "samp_\{j}" }
    sd.push(MSnSampleData::new(name, "", "", 0, ""))
  }
  MSnSet::new(exprs, fd, sd)
}

///|
/// Get n_features.
pub fn MSnSet::n_features(self : MSnSet) -> Int { self.n_features }

///|
/// Get n_samples.
pub fn MSnSet::n_samples(self : MSnSet) -> Int { self.n_samples }

///|
/// Get expression matrix (intensities).
pub fn MSnSet::exprs(self : MSnSet) -> Array[Array[Double]] { self.exprs }

///|
/// Get feature names.
pub fn MSnSet::feature_names(self : MSnSet) -> Array[String] {
  self.feature_data.map(fn(f) { f.feature_name })
}

///|
/// Get sample names.
pub fn MSnSet::sample_names(self : MSnSet) -> Array[String] {
  self.sample_data.map(fn(s) { s.sample_name })
}

///|
/// Get expression for a specific feature across samples.
pub fn MSnSet::get_feature(self : MSnSet, feature_name : String) -> Array[Double] {
  let mut idx = -1
  for i in 0..<self.n_features {
    if self.feature_data[i].feature_name == feature_name { idx = i; break }
  }
  if idx < 0 { return [] }
  self.exprs[idx].copy()
}

///|
/// Get expression for a specific sample across features.
pub fn MSnSet::get_sample(self : MSnSet, sample_name : String) -> Array[Double] {
  let mut idx = -1
  for j in 0..<self.n_samples {
    if self.sample_data[j].sample_name == sample_name { idx = j; break }
  }
  if idx < 0 { return [] }
  let result : Array[Double] = []
  for i in 0..<self.n_features {
    result.push(self.exprs[i][idx])
  }
  result
}

// ============================================================================
// Processing / transformation utilities
// ============================================================================

///|
/// Absolute value helper.
fn mbn_abs(x : Double) -> Double {
  if x < 0.0 { -x } else { x }
}

///|
/// Log2 transform intensities in MSnSet (with offset for zeros).
pub fn MSnSet::log2_transform(
  self : MSnSet,
  offset~ : Double = 1.0,
) -> MSnSet {
  if self.is_log2 { return self }
  let n_f = self.n_features
  let n_s = self.n_samples
  let new_exprs : Array[Array[Double]] = []
  for i in 0..<n_f {
    let row : Array[Double] = []
    for j in 0..<n_s {
      let v = self.exprs[i][j] + offset
      if v <= 0.0 {
        row.push(0.0)
      } else {
        row.push(@math.log2(v))
      }
    }
    new_exprs.push(row)
  }
  let result = MSnSet::new(new_exprs, self.feature_data, self.sample_data)
  result.is_log2 = true
  result
}

///|
/// Missing value imputation:
/// method = "mean" (replace with row mean), "median" (row median),
///          "zero" (set to 0), "min" (row min / 2)
pub fn MSnSet::impute_missing(
  self : MSnSet,
  method~ : String = "mean",
) -> MSnSet {
  let n_f = self.n_features
  let n_s = self.n_samples
  let new_exprs : Array[Array[Double]] = []
  for i in 0..<n_f {
    let row : Array[Double] = []
    for v in self.exprs[i] { row.push(v) }
    // Collect non-NA values (treat NaN and negative as missing)
    let non_missing : Array[Double] = []
    for j in 0..<n_s {
      let v = self.exprs[i][j]
      if !v.is_nan() && v > 0.0 {
        non_missing.push(v)
      }
    }
    if non_missing.length() > 0 {
      let replace_val = if method == "median" {
        // Simple median
        let sorted = mbn_sort_dbl(non_missing)
        let k = sorted.length()
        if k % 2 == 0 {
          (sorted[k / 2 - 1] + sorted[k / 2]) / 2.0
        } else {
          sorted[k / 2]
        }
      } else if method == "min" {
        let sorted = mbn_sort_dbl(non_missing)
        sorted[0] / 2.0
      } else if method == "zero" {
        0.0
      } else {
        // mean
        let mut s = 0.0
        for v in non_missing { s = s + v }
        s / non_missing.length().to_double()
      }
      for j in 0..<n_s {
        let v = self.exprs[i][j]
        if v.is_nan() || v <= 0.0 {
          row[j] = replace_val
        }
      }
    }
    new_exprs.push(row)
  }
  MSnSet::new(new_exprs, self.feature_data, self.sample_data)
}

///|
/// Normalize each sample (column) by scaling to sum = target (default 1e6).
/// Common for proteomics TMT data (sum normalization).
pub fn MSnSet::normalize_sum(
  self : MSnSet,
  target~ : Double = 1000000.0,
) -> MSnSet {
  let n_f = self.n_features
  let n_s = self.n_samples
  // First fill rows with copies of original data
  let new_exprs : Array[Array[Double]] = []
  for i in 0..<n_f {
    let row : Array[Double] = []
    for _j in 0..<n_s { row.push(0.0) }
    new_exprs.push(row)
  }
  for j in 0..<n_s {
    let mut col_sum = 0.0
    for i in 0..<n_f {
      if self.exprs[i][j] > 0.0 {
        col_sum = col_sum + self.exprs[i][j]
      }
    }
    if col_sum > 0.0 {
      let scale = target / col_sum
      for i in 0..<n_f {
        new_exprs[i][j] = if self.exprs[i][j] > 0.0 {
          self.exprs[i][j] * scale
        } else {
          self.exprs[i][j]
        }
      }
    } else {
      for i in 0..<n_f { new_exprs[i][j] = self.exprs[i][j] }
    }
  }
  MSnSet::new(new_exprs, self.feature_data, self.sample_data)
}

///|
/// Median centering normalization: subtract per-sample median from each feature.
pub fn MSnSet::normalize_median_center(self : MSnSet) -> MSnSet {
  let n_f = self.n_features
  let n_s = self.n_samples
  let new_exprs : Array[Array[Double]] = []
  for i in 0..<n_f {
    let row : Array[Double] = []
    for _j in 0..<n_s { row.push(0.0) }
    new_exprs.push(row)
  }
  for j in 0..<n_s {
    let col_vals : Array[Double] = []
    for i in 0..<n_f {
      if self.exprs[i][j] > 0.0 {
        col_vals.push(self.exprs[i][j])
      }
    }
    let median_j = if col_vals.length() > 0 {
      let sorted = mbn_sort_dbl(col_vals)
      let k = sorted.length()
      if k % 2 == 0 {
        (sorted[k / 2 - 1] + sorted[k / 2]) / 2.0
      } else {
        sorted[k / 2]
      }
    } else {
      0.0
    }
    for i in 0..<n_f {
      if self.exprs[i][j] > 0.0 {
        new_exprs[i][j] = self.exprs[i][j] - median_j
      } else {
        new_exprs[i][j] = self.exprs[i][j]
      }
    }
  }
  MSnSet::new(new_exprs, self.feature_data, self.sample_data)
}

///|
/// Summarize features (peptides) into protein-level using
/// specified summarization method ("sum", "mean", "median").
/// Returns new MSnSet with one row per unique protein accession.
pub fn MSnSet::summarize_to_proteins(
  self : MSnSet,
  method~ : String = "sum",
) -> MSnSet {
  // Map from protein accession -> list of row indices
  let protein_rows : Map[String, Array[Int]] = Map([], capacity=0)
  for i in 0..<self.n_features {
    let prots = self.feature_data[i].protein_accessions
    if prots.length() == 0 {
      // Use feature name as pseudo-protein
      let name = self.feature_data[i].feature_name
      if protein_rows.contains(name) {
        let existing = protein_rows.get(name)
        match existing {
          Some(arr) => {
            let new_arr : Array[Int] = []
            for x in arr { new_arr.push(x) }
            new_arr.push(i)
            protein_rows.set(name, new_arr)
          }
          None => protein_rows.set(name, [i])
        }
      } else {
        protein_rows.set(name, [i])
      }
    } else {
      for p in prots {
        if protein_rows.contains(p) {
          let existing = protein_rows.get(p)
          match existing {
            Some(arr) => {
              let new_arr : Array[Int] = []
              for x in arr { new_arr.push(x) }
              new_arr.push(i)
              protein_rows.set(p, new_arr)
            }
            None => protein_rows.set(p, [i])
          }
        } else {
          protein_rows.set(p, [i])
        }
      }
    }
  }
  let all_prots = protein_rows.keys().to_array()
  let n_prot = all_prots.length()
  let n_s = self.n_samples

  let new_exprs : Array[Array[Double]] = []
  for _i in 0..<n_prot {
    let row : Array[Double] = []
    for _j in 0..<n_s { row.push(0.0) }
    new_exprs.push(row)
  }
  let new_fd : Array[MSnFeatureData] = []
  for p_idx in 0..<n_prot {
    let pname = all_prots[p_idx]
    let row_opt = protein_rows.get(pname)
    let rows = match row_opt {
      Some(r) => r
      None => []
    }
    // Summarize each sample column
    let prot_row : Array[Double] = []
    for _j in 0..<n_s { prot_row.push(0.0) }
    for j in 0..<n_s {
      let vals : Array[Double] = []
      for r in rows {
        if r < self.n_features && self.exprs[r][j] > 0.0 {
          vals.push(self.exprs[r][j])
        }
      }
      if vals.length() > 0 {
        if method == "median" {
          let sorted = mbn_sort_dbl(vals)
          let k = sorted.length()
          prot_row[j] = if k % 2 == 0 {
            (sorted[k / 2 - 1] + sorted[k / 2]) / 2.0
          } else {
            sorted[k / 2]
          }
        } else if method == "mean" {
          let mut s = 0.0
          for v in vals { s = s + v }
          prot_row[j] = s / vals.length().to_double()
        } else {
          // sum
          let mut s = 0.0
          for v in vals { s = s + v }
          prot_row[j] = s
        }
      }
    }
    new_exprs[p_idx] = prot_row
    new_fd.push(MSnFeatureData::new(
      pname, [pname], "", 0, 0.0, 0.0,
    ))
  }
  MSnSet::new(new_exprs, new_fd, self.sample_data)
}

///|
/// Simple Savitzky-Golay-like smoothing via moving average for spectra.
/// Window size is the full window (will be made odd). Returns smoothed Spectrum.
pub fn Spectrum::smooth_moving_average(
  self : Spectrum,
  half_window~ : Int = 2,
) -> Spectrum {
  let n = self.peaks_count
  if n == 0 { return self }
  let hw = if half_window < 1 { 1 } else { half_window }
  let new_int : Array[Double] = []
  for i in 0..<n {
    let mut sum = 0.0
    let mut count = 0
    let lo = if i - hw > 0 { i - hw } else { 0 }
    let hi = if i + hw < n { i + hw } else { n - 1 }
    for k in lo..<(hi + 1) {
      sum = sum + self.intensity[k]
      count = count + 1
    }
    new_int.push(if count > 0 { sum / count.to_double() } else { 0.0 })
  }
  Spectrum::new(self.mz, new_int, self.ms_level, self.polarity, self.rt)
}

///|
/// Simple baseline correction using asymmetric least squares-style simplified:
/// subtract minimum from a sliding window.
pub fn Spectrum::baseline_correct_minwin(
  self : Spectrum,
  half_window~ : Int = 20,
) -> Spectrum {
  let n = self.peaks_count
  if n == 0 { return self }
  let hw = if half_window < 1 { 1 } else { half_window }
  // Compute baseline estimate per point = min in window
  let baseline : Array[Double] = []
  for i in 0..<n {
    let lo = if i - hw > 0 { i - hw } else { 0 }
    let hi = if i + hw < n { i + hw } else { n - 1 }
    let mut min_v = self.intensity[lo]
    for k in lo..<(hi + 1) {
      if self.intensity[k] < min_v { min_v = self.intensity[k] }
    }
    baseline.push(min_v)
  }
  let new_int : Array[Double] = []
  for i in 0..<n {
    let v = self.intensity[i] - baseline[i]
    new_int.push(if v > 0.0 { v } else { 0.0 })
  }
  Spectrum::new(self.mz, new_int, self.ms_level, self.polarity, self.rt)
}

///|
/// Simple peak picking / centroiding: identify local maxima above noise threshold.
/// Returns a new centroided Spectrum with only peak locations.
pub fn Spectrum::centroid_simple(
  self : Spectrum,
  snr_threshold~ : Double = 3.0,
) -> Spectrum {
  let n = self.peaks_count
  if n < 3 { return self }
  // Estimate noise as median of small differences
  let noise_estimate = 0.0
  // Compute median intensity / 50 as simple noise floor
  let sorted_int = mbn_sort_dbl(self.intensity)
  let noise_estimate = sorted_int[n / 4]
  let noise_estimate = if noise_estimate <= 0.0 { 1.0 } else { noise_estimate }
  let threshold = noise_estimate * snr_threshold

  let new_mz : Array[Double] = []
  let new_int : Array[Double] = []
  for i in 1..<(n - 1) {
    let cur = self.intensity[i]
    if cur >= threshold &&
       cur > self.intensity[i - 1] &&
       cur > self.intensity[i + 1] {
      // Simple 3-point parabolic apex refinement
      let a = self.intensity[i - 1]
      let b = cur
      let c = self.intensity[i + 1]
      let denom = 2.0 * (2.0 * b - a - c)
      if mbn_abs(denom) > 0.0000000001 {
        let delta = (a - c) / denom
        let mz_shift = delta * (self.mz[i] - self.mz[i - 1])
        new_mz.push(self.mz[i] + mz_shift)
      } else {
        new_mz.push(self.mz[i])
      }
      new_int.push(cur)
    }
  }
  let result = Spectrum::new(new_mz, new_int, self.ms_level, self.polarity, self.rt)
  result.centroided = true
  result
}

// ============================================================================
// Helper: sort doubles
// ============================================================================

fn mbn_sort_dbl(arr : Array[Double]) -> Array[Double] {
  let result = arr.copy()
  let n = result.length()
  for i in 0..<n {
    let start = i + 1
    for j in start..<n {
      if result[i] > result[j] {
        let t = result[i]
        result[i] = result[j]
        result[j] = t
      }
    }
  }
  result
}

// ============================================================================
// Quality assessment helpers
// ============================================================================

///|
/// Compute per-sample QC metrics from MSnSet:
/// Returns array of (sample_name, n_features_detected, total_intensity, median_intensity, coef_of_variation)
pub fn MSnSet::sample_qc(self : MSnSet) -> Array[(String, Int, Double, Double, Double)] {
  let n_s = self.n_samples
  let result : Array[(String, Int, Double, Double, Double)] = []
  for j in 0..<n_s {
    let mut detected = 0
    let mut total = 0.0
    let vals : Array[Double] = []
    for i in 0..<self.n_features {
      let v = self.exprs[i][j]
      if v > 0.0 {
        detected = detected + 1
        total = total + v
        vals.push(v)
      }
    }
    let k = vals.length()
    let median_v = if k > 0 {
      let sorted = mbn_sort_dbl(vals)
      if k % 2 == 0 {
        (sorted[k / 2 - 1] + sorted[k / 2]) / 2.0
      } else {
        sorted[k / 2]
      }
    } else {
      0.0
    }
    // Coefficient of variation (sd / mean)
    let mut cv = 0.0
    if k > 1 {
      let mut mean_v = 0.0
      for v in vals { mean_v = mean_v + v }
      mean_v = mean_v / k.to_double()
      let mut sum_sq = 0.0
      for v in vals {
        let d = v - mean_v
        sum_sq = sum_sq + d * d
      }
      let sd = (sum_sq / (k - 1).to_double()).sqrt()
      if mean_v > 0.0 { cv = sd / mean_v }
    }
    let sname = if j < self.sample_data.length() {
      self.sample_data[j].sample_name
    } else {
      "sample_\{j}"
    }
    result.push((sname, detected, total, median_v, cv))
  }
  result
}
