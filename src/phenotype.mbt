///|
/// Bio.phenotype module - Phenotype microarray analysis.
/// Inspired by Biopython's Bio.phenotype package.
/// Provides plate/well record data structures, CSV/JSON parsing,
/// growth curve fitting (logistic, Gompertz), and parameter extraction.

// ===== Core data structures =====

///|
/// WellRecord - a single well's time-series measurements.
pub struct WellRecord {
  /// Well identifier (e.g., "A01")
  well_id : String
  /// Time points (hours)
  times : Array[Double]
  /// Signal measurements (e.g., OD)
  signals : Array[Double]
  /// Plate identifier this well belongs to
  plate_id : String
}

///|
/// PlateRecord - a 96-well plate of phenotype data.
pub struct PlateRecord {
  /// Plate identifier
  plate_id : String
  /// Wells (8 rows A-H x 12 cols 01-12 = up to 96)
  wells : Array[WellRecord]
  /// Time points shared across the plate (hours)
  times : Array[Double]
}

///|
/// Fitted parameters for a growth curve.
pub struct PhenFitParams {
  /// Asymptote (carrying capacity, max signal)
  asymptote : Double
  /// Inflection point (time of fastest growth)
  inflection : Double
  /// Maximum growth slope
  slope : Double
  /// Lag time (start of growth phase)
  lag_time : Double
  /// Plateau signal value (final)
  plateau : Double
  /// Initial signal value
  initial : Double
  /// Model type used ("logistic" or "gompertz")
  model : String
  /// Residual sum of squares
  rss : Double
  /// R-squared value
  r_squared : Double
}

///|
/// Result of control subtraction.
pub struct PhenControlSubtracted {
  /// Well IDs
  well_ids : Array[String]
  /// Time points
  times : Array[Double]
  /// Subtracted signals per well
  signals : Array[Array[Double]]
}

// ===== WellRecord methods =====

///|
/// Create a new WellRecord.
pub fn WellRecord::new(well_id : String, plate_id : String, times : Array[Double], signals : Array[Double]) -> WellRecord {
  WellRecord::{ well_id, plate_id, times, signals }
}

///|
/// Get number of time points.
pub fn WellRecord::phen_n_points(self : WellRecord) -> Int {
  self.times.length()
}

///|
/// Get signal at a time index.
pub fn WellRecord::phen_get_signal(self : WellRecord, idx : Int) -> Double {
  self.signals[idx]
}

///|
/// Get time at an index.
pub fn WellRecord::phen_get_time(self : WellRecord, idx : Int) -> Double {
  self.times[idx]
}

///|
/// Get maximum signal value.
pub fn WellRecord::phen_max_signal(self : WellRecord) -> Double {
  let mut max_val = if self.signals.length() > 0 { self.signals[0] } else { 0.0 }
  let mut i = 1
  while i < self.signals.length() {
    if self.signals[i] > max_val {
      max_val = self.signals[i]
    }
    i = i + 1
  }
  max_val
}

///|
/// Get minimum signal value.
pub fn WellRecord::phen_min_signal(self : WellRecord) -> Double {
  let mut min_val = if self.signals.length() > 0 { self.signals[0] } else { 0.0 }
  let mut i = 1
  while i < self.signals.length() {
    if self.signals[i] < min_val {
      min_val = self.signals[i]
    }
    i = i + 1
  }
  min_val
}

///|
/// Get the area under the curve (trapezoidal rule).
pub fn WellRecord::phen_area(self : WellRecord) -> Double {
  let n = self.times.length()
  if n < 2 {
    return 0.0
  }
  let mut area = 0.0
  let mut i = 1
  while i < n {
    let dt = self.times[i] - self.times[i - 1]
    let avg = (self.signals[i] + self.signals[i - 1]) / 2.0
    area = area + dt * avg
    i = i + 1
  }
  area
}

///|
/// Compute the maximum slope (numerical derivative).
pub fn WellRecord::phen_max_slope(self : WellRecord) -> Double {
  let n = self.times.length()
  if n < 2 {
    return 0.0
  }
  let mut max_slope = 0.0
  let mut i = 1
  while i < n {
    let dt = self.times[i] - self.times[i - 1]
    let local_slope = if dt > 0.0 { (self.signals[i] - self.signals[i - 1]) / dt } else { 0.0 }
    if local_slope > max_slope {
      max_slope = local_slope
    }
    i = i + 1
  }
  max_slope
}

// ===== PlateRecord methods =====

///|
/// Create a new PlateRecord.
pub fn PlateRecord::new(plate_id : String, times : Array[Double]) -> PlateRecord {
  PlateRecord::{ plate_id, wells: Array::new(), times }
}

///|
/// Add a well to the plate.
pub fn PlateRecord::phen_add_well(self : PlateRecord, well : WellRecord) -> PlateRecord {
  let new_wells = self.wells.copy()
  new_wells.push(well)
  PlateRecord::{ plate_id: self.plate_id, wells: new_wells, times: self.times }
}

///|
/// Get number of wells.
pub fn PlateRecord::phen_n_wells(self : PlateRecord) -> Int {
  self.wells.length()
}

///|
/// Get a well by well_id (returns Option).
pub fn PlateRecord::phen_get_well(self : PlateRecord, well_id : String) -> WellRecord? {
  let mut found : WellRecord? = None
  let mut i = 0
  while i < self.wells.length() {
    if self.wells[i].well_id == well_id {
      found = Some(self.wells[i])
      i = self.wells.length()
    }
    i = i + 1
  }
  found
}

///|
/// Get a well at an index.
pub fn PlateRecord::phen_get_well_at(self : PlateRecord, idx : Int) -> WellRecord {
  self.wells[idx]
}

///|
/// Number of time points in the plate.
pub fn PlateRecord::phen_n_times(self : PlateRecord) -> Int {
  self.times.length()
}

///|
/// Get mean signal across all wells at a time index.
pub fn PlateRecord::phen_mean_signal_at(self : PlateRecord, time_idx : Int) -> Double {
  if self.wells.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < self.wells.length() {
    if time_idx < self.wells[i].signals.length() {
      sum = sum + self.wells[i].signals[time_idx]
    }
    i = i + 1
  }
  sum / self.wells.length().to_double()
}

// ===== Curve fitting =====

///|
/// Estimate initial parameters from raw data.
fn phen_estimate_initial(well : WellRecord) -> (Double, Double, Double, Double) {
  let n = well.signals.length()
  if n == 0 {
    return (0.0, 0.0, 0.0, 0.0)
  }
  let initial = well.signals[0]
  let plateau = well.phen_max_signal()
  let slope = well.phen_max_slope()
  // Inflection estimate: time of max derivative
  let mut inflection = 0.0
  let mut max_slope_val = 0.0
  let mut i = 1
  while i < n {
    let dt = well.times[i] - well.times[i - 1]
    let local_slope = if dt > 0.0 { (well.signals[i] - well.signals[i - 1]) / dt } else { 0.0 }
    if local_slope > max_slope_val {
      max_slope_val = local_slope
      inflection = well.times[i]
    }
    i = i + 1
  }
  (initial, plateau, slope, inflection)
}

///|
/// Logistic model: L / (1 + exp(-k * (t - t0))) + offset.
/// L = asymptote, k = growth rate, t0 = inflection, offset = initial.
fn phen_logistic_model(t : Double, asymptote : Double, slope : Double, inflection : Double, initial : Double) -> Double {
  // f(t) = initial + (asymptote - initial) / (1 + exp(-slope * (t - inflection)))
  let exponent = -slope * (t - inflection)
  let denom = 1.0 + phen_safe_exp(exponent)
  initial + (asymptote - initial) / denom
}

///|
/// Gompertz model: initial + (asymptote - initial) * exp(-exp(-slope * (t - inflection))).
fn phen_gompertz_model(t : Double, asymptote : Double, slope : Double, inflection : Double, initial : Double) -> Double {
  let inner = -slope * (t - inflection)
  let outer = -phen_safe_exp(inner)
  initial + (asymptote - initial) * phen_safe_exp(outer)
}

///|
/// Safe exponential that won't overflow.
fn phen_safe_exp(x : Double) -> Double {
  if x > 700.0 {
    return 1.0e304
  }
  if x < -700.0 {
    return 0.0
  }
  @math.exp(x)
}

///|
/// Compute logistic curve fit using simplified grid search around estimates.
pub fn phen_fit_logistic(well : WellRecord) -> PhenFitParams {
  let (initial, plateau, slope_est, inflection_est) = phen_estimate_initial(well)
  let asymptote = plateau
  // Refine inflection and slope via coarse grid search
  let mut best_rss = 1.0e300
  let mut best_slope = slope_est
  let mut best_inflection = inflection_est
  let slope_candidates = [slope_est * 0.5, slope_est, slope_est * 1.5, slope_est * 2.0, 0.1, 0.5, 1.0]
  let infl_candidates = [inflection_est - 5.0, inflection_est - 2.0, inflection_est, inflection_est + 2.0, inflection_est + 5.0, 0.0, 5.0, 10.0, 15.0]
  let mut si = 0
  while si < slope_candidates.length() {
    let mut ii = 0
    while ii < infl_candidates.length() {
      let s = slope_candidates[si]
      let infl = infl_candidates[ii]
      let mut rss = 0.0
      let mut k = 0
      while k < well.times.length() {
        let pred = phen_logistic_model(well.times[k], asymptote, s, infl, initial)
        let diff = well.signals[k] - pred
        rss = rss + diff * diff
        k = k + 1
      }
      if rss < best_rss {
        best_rss = rss
        best_slope = s
        best_inflection = infl
      }
      ii = ii + 1
    }
    si = si + 1
  }
  // Lag time estimate: inflection - 1/slope
  let lag_time = if best_slope > 0.0 { best_inflection - 1.0 / best_slope } else { 0.0 }
  let lag_clamped = if lag_time < 0.0 { 0.0 } else { lag_time }
  // R-squared
  let mean_y = phen_mean(well.signals)
  let ss_tot = phen_sum_sq(well.signals, mean_y)
  let r_sq = if ss_tot > 0.0 { 1.0 - best_rss / ss_tot } else { 0.0 }
  let r_sq_clamped = if r_sq < 0.0 { 0.0 } else { r_sq }
  PhenFitParams::{
    asymptote,
    inflection: best_inflection,
    slope: best_slope,
    lag_time: lag_clamped,
    plateau,
    initial,
    model: "logistic",
    rss: best_rss,
    r_squared: r_sq_clamped
  }
}

///|
/// Compute Gompertz curve fit using simplified grid search.
pub fn phen_fit_gompertz(well : WellRecord) -> PhenFitParams {
  let (initial, plateau, slope_est, inflection_est) = phen_estimate_initial(well)
  let asymptote = plateau
  let mut best_rss = 1.0e300
  let mut best_slope = slope_est
  let mut best_inflection = inflection_est
  let slope_candidates = [slope_est * 0.5, slope_est, slope_est * 1.5, slope_est * 2.0, 0.1, 0.5, 1.0]
  let infl_candidates = [inflection_est - 5.0, inflection_est - 2.0, inflection_est, inflection_est + 2.0, inflection_est + 5.0, 0.0, 5.0, 10.0, 15.0]
  let mut si = 0
  while si < slope_candidates.length() {
    let mut ii = 0
    while ii < infl_candidates.length() {
      let s = slope_candidates[si]
      let infl = infl_candidates[ii]
      let mut rss = 0.0
      let mut k = 0
      while k < well.times.length() {
        let pred = phen_gompertz_model(well.times[k], asymptote, s, infl, initial)
        let diff = well.signals[k] - pred
        rss = rss + diff * diff
        k = k + 1
      }
      if rss < best_rss {
        best_rss = rss
        best_slope = s
        best_inflection = infl
      }
      ii = ii + 1
    }
    si = si + 1
  }
  let lag_time = if best_slope > 0.0 { best_inflection - 1.0 / best_slope } else { 0.0 }
  let lag_clamped = if lag_time < 0.0 { 0.0 } else { lag_time }
  let mean_y = phen_mean(well.signals)
  let ss_tot = phen_sum_sq(well.signals, mean_y)
  let r_sq = if ss_tot > 0.0 { 1.0 - best_rss / ss_tot } else { 0.0 }
  let r_sq_clamped = if r_sq < 0.0 { 0.0 } else { r_sq }
  PhenFitParams::{
    asymptote,
    inflection: best_inflection,
    slope: best_slope,
    lag_time: lag_clamped,
    plateau,
    initial,
    model: "gompertz",
    rss: best_rss,
    r_squared: r_sq_clamped
  }
}

///|
/// Compute the mean of an array.
fn phen_mean(arr : Array[Double]) -> Double {
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
/// Compute sum of squared deviations from a value.
fn phen_sum_sq(arr : Array[Double], center : Double) -> Double {
  let mut ss = 0.0
  let mut i = 0
  while i < arr.length() {
    let d = arr[i] - center
    ss = ss + d * d
    i = i + 1
  }
  ss
}

///|
/// Compute predicted values from a fitted model.
pub fn phen_predict(fit : PhenFitParams, times : Array[Double]) -> Array[Double] {
  let preds : Array[Double] = Array::new()
  let mut i = 0
  while i < times.length() {
    let pred = if fit.model == "logistic" {
      phen_logistic_model(times[i], fit.asymptote, fit.slope, fit.inflection, fit.initial)
    } else if fit.model == "gompertz" {
      phen_gompertz_model(times[i], fit.asymptote, fit.slope, fit.inflection, fit.initial)
    } else {
      fit.asymptote
    }
    preds.push(pred)
    i = i + 1
  }
  preds
}

// ===== Control subtraction =====

///|
/// Subtract control well signals from sample wells.
pub fn phen_subtract_control(plate : PlateRecord, control_well_id : String) -> PhenControlSubtracted {
  let control_opt = plate.phen_get_well(control_well_id)
  let well_ids : Array[String] = Array::new()
  let signals : Array[Array[Double]] = Array::new()
  let mut control_signals : Array[Double] = Array::new()
  if control_opt is Some(_) {
    control_signals = control_opt.unwrap().signals
  }
  let mut i = 0
  while i < plate.wells.length() {
    let w = plate.wells[i]
    if w.well_id != control_well_id {
      well_ids.push(w.well_id)
      let sub : Array[Double] = Array::new()
      let mut j = 0
      while j < w.signals.length() {
        let ctrl = if j < control_signals.length() { control_signals[j] } else { 0.0 }
        let v = w.signals[j] - ctrl
        sub.push(v)
        j = j + 1
      }
      signals.push(sub)
    }
    i = i + 1
  }
  PhenControlSubtracted::{ well_ids, times: plate.times, signals }
}

// ===== CSV / JSON parsing =====

///|
/// Parse a CSV string into a PlateRecord.
/// Expected format (Biopython-like):
///   Plate,Well,Time,Signal
///   P1,A01,0.0,0.05
///   P1,A01,1.0,0.06
///   ...
pub fn phen_parse_csv(csv_text : String) -> PlateRecord {
  let lines = csv_text.split("\n").to_array()
  // Group wells
  let well_map : Map[String, WellRecord] = Map([], capacity=96)
  let times_collected : Array[Double] = Array::new()
  let plate_id_var : Array[String] = Array::new()
  // Skip header line; parse rest
  let mut idx = 0
  while idx < lines.length() {
    let line = lines[idx].trim()
    if line.length() == 0 || idx == 0 {
      idx = idx + 1
      continue
    }
    let parts = line.split(",").to_array()
    if parts.length() >= 4 {
      let pid = parts[0].trim().to_owned()
      let wid = parts[1].trim().to_owned()
      let t = phen_parse_double(parts[2].trim().to_owned())
      let s = phen_parse_double(parts[3].trim().to_owned())
      if plate_id_var.length() == 0 {
        plate_id_var.push(pid)
      }
      let existing = well_map.get(wid)
      if existing is Some(_) {
        let w = existing.unwrap()
        let new_times = w.times.copy()
        new_times.push(t)
        let new_signals = w.signals.copy()
        new_signals.push(s)
        well_map.set(wid, WellRecord::{ well_id: wid, plate_id: pid, times: new_times, signals: new_signals })
      } else {
        let new_times : Array[Double] = Array::new()
        new_times.push(t)
        let new_signals : Array[Double] = Array::new()
        new_signals.push(s)
        well_map.set(wid, WellRecord::{ well_id: wid, plate_id: pid, times: new_times, signals: new_signals })
      }
      times_collected.push(t)
    }
    idx = idx + 1
  }
  let pid = if plate_id_var.length() > 0 { plate_id_var[0] } else { "plate" }
  // Deduplicate time points
  let unique_times = phen_dedup_sort(times_collected)
  // Build plate
  let plate = PlateRecord::new(pid, unique_times)
  // Add wells in alphabetical order
  let well_ids = well_map.keys().collect()
  let well_ids_sorted = phen_sort_well_ids(well_ids)
  let mut p = plate
  let mut wi = 0
  while wi < well_ids_sorted.length() {
    let w_opt = well_map.get(well_ids_sorted[wi])
    if w_opt is Some(_) {
      p = p.phen_add_well(w_opt.unwrap())
    }
    wi = wi + 1
  }
  p
}

///|
/// Parse a JSON-like string into a PlateRecord.
/// Expected simplified format:
///   {"plate":"P1","times":[0,1,2],"wells":[{"id":"A01","signals":[0.05,0.06,0.1]}, ...]}
pub fn phen_parse_json(json_text : String) -> PlateRecord {
  // Very lightweight JSON parser for the specific subset we expect.
  let trimmed = json_text.trim().to_owned()
  let plate_id = phen_json_extract_string(trimmed, "plate")
  let times = phen_json_extract_double_array(trimmed, "times")
  let wells = phen_json_extract_wells(trimmed)
  let mut p = PlateRecord::new(plate_id, times)
  let mut i = 0
  while i < wells.length() {
    p = p.phen_add_well(wells[i])
    i = i + 1
  }
  p
}

///|
/// Extract a string value for a key from simplified JSON.
fn phen_json_extract_string(json : String, key : String) -> String {
  let pattern = "\"" + key + "\":\""
  let start = json.find(pattern)
  if start is Some(_) {
    let s = start.unwrap() + pattern.length()
    let rest = json[s:].to_owned()
    let end = rest.find("\"")
    if end is Some(_) {
      return rest[0:end.unwrap()].to_owned()
    }
  }
  "plate"
}

///|
/// Extract a numeric array for a key from simplified JSON.
fn phen_json_extract_double_array(json : String, key : String) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let pattern = "\"" + key + "\":["
  let start = json.find(pattern)
  if start is None {
    return result
  }
  let s = start.unwrap() + pattern.length()
  let rest = json[s:].to_owned()
  let end = rest.find("]")
  if end is None {
    return result
  }
  let arr_str = rest[0:end.unwrap()].to_owned()
  let parts = arr_str.split(",").to_array()
  let mut i = 0
  while i < parts.length() {
    let v = phen_parse_double(parts[i].trim().to_owned())
    result.push(v)
    i = i + 1
  }
  result
}

///|
/// Extract wells array from simplified JSON.
fn phen_json_extract_wells(json : String) -> Array[WellRecord] {
  let wells : Array[WellRecord] = Array::new()
  let pattern = "\"wells\":["
  let start = json.find(pattern)
  if start is None {
    return wells
  }
  let mut pos = start.unwrap() + pattern.length()
  // Find each well object { ... }
  let mut depth = 1
  while pos < json.length() && depth > 0 {
    let ch = json[pos]
    if ch == '{' {
      // Find the matching close brace
      let obj_start = pos
      let mut obj_depth = 1
      pos = pos + 1
      while pos < json.length() && obj_depth > 0 {
        let c2 = json[pos]
        if c2 == '{' {
          obj_depth = obj_depth + 1
        }
        if c2 == '}' {
          obj_depth = obj_depth - 1
        }
        if obj_depth > 0 {
          pos = pos + 1
        }
      }
      let obj_str = json[obj_start:pos + 1].to_owned()
      let wid = phen_json_extract_string(obj_str, "id")
      let sigs = phen_json_extract_double_array(obj_str, "signals")
      let pid = phen_json_extract_string(json, "plate")
      let times = phen_json_extract_double_array(json, "times")
      wells.push(WellRecord::{ well_id: wid, plate_id: pid, times, signals: sigs })
    }
    if json[pos] == ']' {
      depth = depth - 1
    }
    pos = pos + 1
  }
  wells
}

///|
/// Parse a double from a string (handles simple decimals and integers).
fn phen_parse_double(s : String) -> Double {
  let trimmed = s.trim().to_owned()
  if trimmed.length() == 0 {
    return 0.0
  }
  let neg = trimmed[0] == '-'
  let digits : String = if neg { trimmed[1:].to_owned() } else { trimmed }
  let mut int_part = 0.0
  let mut frac_part = 0.0
  let mut frac_div = 1.0
  let mut in_frac = false
  let mut i = 0
  while i < digits.length() {
    let c = digits[i]
    if c >= '0' && c <= '9' {
      let d = (c.to_int() - '0'.to_int()).to_double()
      if in_frac {
        frac_div = frac_div * 10.0
        frac_part = frac_part + d / frac_div
      } else {
        int_part = int_part * 10.0 + d
      }
    } else if c == '.' {
      in_frac = true
    } else {
      // Stop at first non-numeric (e.g. 'e' for exponent); treat as 0 contribution
      ignore(c)
    }
    i = i + 1
  }
  let v = int_part + frac_part
  if neg { -v } else { v }
}

///|
/// Deduplicate and sort a list of doubles.
fn phen_dedup_sort(arr : Array[Double]) -> Array[Double] {
  let sorted = arr.copy()
  // Simple insertion sort
  let mut i = 1
  while i < sorted.length() {
    let key = sorted[i]
    let mut j = i - 1
    while j >= 0 && sorted[j] > key {
      sorted[j + 1] = sorted[j]
      j = j - 1
    }
    sorted[j + 1] = key
    i = i + 1
  }
  let result : Array[Double] = Array::new()
  let mut k = 0
  while k < sorted.length() {
    if result.length() == 0 || (result[result.length() - 1] - sorted[k]).abs() > 1.0e-9 {
      result.push(sorted[k])
    }
    k = k + 1
  }
  result
}

///|
/// Sort well IDs alphabetically (e.g., A01, A02, ..., B01).
fn phen_sort_well_ids(ids : Array[String]) -> Array[String] {
  let sorted = ids.copy()
  let mut i = 1
  while i < sorted.length() {
    let key = sorted[i]
    let mut j = i - 1
    while j >= 0 && sorted[j] > key {
      sorted[j + 1] = sorted[j]
      j = j - 1
    }
    sorted[j + 1] = key
    i = i + 1
  }
  sorted
}

// ===== Example data generators =====

///|
/// Create a small example plate with growth curves.
pub fn phen_create_example_plate() -> PlateRecord {
  let times = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0]
  let mut plate = PlateRecord::new("ExamplePlate", times)
  // Well A01 - control
  let a01_signals = [0.05, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
  plate = plate.phen_add_well(WellRecord::new("A01", "ExamplePlate", times, a01_signals))
  // Well A02 - growing culture
  let a02_signals = [0.05, 0.06, 0.08, 0.12, 0.2, 0.35, 0.55, 0.75, 0.9, 0.97, 0.99, 1.0, 1.0]
  plate = plate.phen_add_well(WellRecord::new("A02", "ExamplePlate", times, a02_signals))
  // Well A03 - another growing culture
  let a03_signals = [0.04, 0.05, 0.06, 0.08, 0.12, 0.2, 0.32, 0.5, 0.7, 0.85, 0.93, 0.97, 0.99]
  plate = plate.phen_add_well(WellRecord::new("A03", "ExamplePlate", times, a03_signals))
  // Well B01 - slow growth
  let b01_signals = [0.05, 0.05, 0.05, 0.06, 0.07, 0.09, 0.12, 0.17, 0.25, 0.35, 0.48, 0.6, 0.72]
  plate = plate.phen_add_well(WellRecord::new("B01", "ExamplePlate", times, b01_signals))
  plate
}

///|
/// Create a single example well with a logistic-like growth curve.
pub fn phen_create_example_well() -> WellRecord {
  let times = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
  let signals = [0.05, 0.06, 0.08, 0.13, 0.25, 0.5, 0.78, 0.93, 0.98, 1.0, 1.0]
  WellRecord::new("A01", "ExamplePlate", times, signals)
}

///|
/// Example CSV string for testing.
pub fn phen_example_csv() -> String {
  "Plate,Well,Time,Signal\nP1,A01,0.0,0.05\nP1,A01,2.0,0.06\nP1,A01,4.0,0.08\nP1,A02,0.0,0.05\nP1,A02,2.0,0.1\nP1,A02,4.0,0.2\n"
}
