///|
/// nucleR module, inspired by Bioconductor nucleR package.
/// Provides functions for nucleosome positioning analysis from sequencing data.
/// Supports nucleosome calling, occupancy calculation, and dynamic nucleosome analysis.

/// Nucleosome positioning result
pub struct NucPosition {
  /// Chromosome
  chrom : String
  /// Start position
  start : Int
  /// End position
  end : Int
  /// Center position
  center : Double
  /// Width
  width : Double
  /// Occupancy score
  occupancy : Double
  /// Positioning score
  positioning_score : Double
  /// Is well-positioned
  well_positioned : Bool
}

/// Create new nucleosome position
pub fn NucPosition::new(chrom : String, start : Int, end : Int, center : Double,
                        width : Double, occupancy : Double,
                        positioning_score : Double, well_positioned : Bool) -> NucPosition {
  NucPosition::{ chrom, start, end, center, width, occupancy, positioning_score, well_positioned }
}

/// Get chromosome
pub fn NucPosition::nuc_chrom(self : NucPosition) -> String {
  self.chrom
}

/// Get start position
pub fn NucPosition::nuc_start(self : NucPosition) -> Int {
  self.start
}

/// Get end position
pub fn NucPosition::nuc_end(self : NucPosition) -> Int {
  self.end
}

/// Get center position
pub fn NucPosition::nuc_center(self : NucPosition) -> Double {
  self.center
}

/// Get width
pub fn NucPosition::nuc_width(self : NucPosition) -> Double {
  self.width
}

/// Get occupancy
pub fn NucPosition::nuc_occupancy(self : NucPosition) -> Double {
  self.occupancy
}

/// Get positioning score
pub fn NucPosition::nuc_positioning_score(self : NucPosition) -> Double {
  self.positioning_score
}

/// Check if well-positioned
pub fn NucPosition::nuc_is_well_positioned(self : NucPosition) -> Bool {
  self.well_positioned
}

/// Nucleosome call result
pub struct NucCallResult {
  /// Nucleosome positions
  nucleosomes : Array[NucPosition]
  /// Number of nucleosomes
  n_nucleosomes : Int
  /// Mean nucleosome spacing
  mean_spacing : Double
  /// Mean occupancy
  mean_occupancy : Double
  /// Fraction well-positioned
  frac_well_positioned : Double
}

/// Create new call result
pub fn NucCallResult::new(nucleosomes : Array[NucPosition], n_nucleosomes : Int,
                          mean_spacing : Double, mean_occupancy : Double,
                          frac_well_positioned : Double) -> NucCallResult {
  NucCallResult::{ nucleosomes, n_nucleosomes, mean_spacing, mean_occupancy, frac_well_positioned }
}

/// Get number of nucleosomes
pub fn NucCallResult::nuc_n_count(self : NucCallResult) -> Int {
  self.n_nucleosomes
}

/// Get nucleosome at index
pub fn NucCallResult::nuc_get(self : NucCallResult, idx : Int) -> NucPosition {
  self.nucleosomes[idx]
}

/// Get mean spacing
pub fn NucCallResult::nuc_mean_spacing(self : NucCallResult) -> Double {
  self.mean_spacing
}

/// Get mean occupancy
pub fn NucCallResult::nuc_mean_occupancy(self : NucCallResult) -> Double {
  self.mean_occupancy
}

/// Get fraction well-positioned
pub fn NucCallResult::nuc_frac_well_positioned(self : NucCallResult) -> Double {
  self.frac_well_positioned
}

/// Filter by occupancy threshold
pub fn NucCallResult::nuc_filter_by_occupancy(self : NucCallResult, threshold : Double) -> NucCallResult {
  let filtered : Array[NucPosition] = Array::new()

  let mut i = 0
  while i < self.nucleosomes.length() {
    if self.nucleosomes[i].occupancy >= threshold {
      filtered.push(self.nucleosomes[i])
    }
    i = i + 1
  }

  NucCallResult::new(filtered, filtered.length(), self.mean_spacing,
                      self.mean_occupancy, self.frac_well_positioned)
}

/// Dynamic nucleosome result for comparison
pub struct NucDynamicResult {
  /// Shared nucleosomes
  shared : Array[NucPosition]
  /// Gained nucleosomes
  gained : Array[NucPosition]
  /// Lost nucleosomes
  lost : Array[NucPosition]
  /// Fraction changed
  frac_changed : Double
  /// Direction of change (+1 gain, -1 loss)
  direction : Int
}

/// Create new dynamic result
pub fn NucDynamicResult::new(shared : Array[NucPosition], gained : Array[NucPosition],
                              lost : Array[NucPosition], frac_changed : Double, direction : Int) -> NucDynamicResult {
  NucDynamicResult::{ shared, gained, lost, frac_changed, direction }
}

/// Get number of shared nucleosomes
pub fn NucDynamicResult::nuc_n_shared(self : NucDynamicResult) -> Int {
  self.shared.length()
}

/// Get number of gained nucleosomes
pub fn NucDynamicResult::nuc_n_gained(self : NucDynamicResult) -> Int {
  self.gained.length()
}

/// Get number of lost nucleosomes
pub fn NucDynamicResult::nuc_n_lost(self : NucDynamicResult) -> Int {
  self.lost.length()
}

/// Get fraction changed
pub fn NucDynamicResult::nuc_frac_changed(self : NucDynamicResult) -> Double {
  self.frac_changed
}

/// Get direction
pub fn NucDynamicResult::nuc_direction(self : NucDynamicResult) -> Int {
  self.direction
}

/// Smoothing methods
pub enum NucSmoothMethod {
  /// Moving average
  MovingAverage
  /// Gaussian kernel
  Gaussian
  /// Savitzky-Golay
  SavitzkyGolay
}

/// Nucleosome calling parameters
pub struct NucCallParams {
  /// Smoothing window size
  smooth_window : Int
  /// Peak detection threshold
  peak_threshold : Double
  /// Minimum distance between nucleosomes
  min_distance : Int
  /// Maximum distance between nucleosomes
  max_distance : Int
  /// Occupancy threshold
  occupancy_threshold : Double
  /// Width threshold for well-positioned
  width_threshold : Double
}

/// Create default parameters
pub fn NucCallParams::new() -> NucCallParams {
  NucCallParams::{ smooth_window: 100, peak_threshold: 0.5, min_distance: 100, max_distance: 300,
    occupancy_threshold: 0.3, width_threshold: 60.0 }
}

/// Set smoothing window
pub fn NucCallParams::nuc_set_smooth_window(self : NucCallParams, window : Int) -> NucCallParams {
  NucCallParams::{ smooth_window: window, peak_threshold: self.peak_threshold, min_distance: self.min_distance, max_distance: self.max_distance, occupancy_threshold: self.occupancy_threshold, width_threshold: self.width_threshold }
}

/// Set peak threshold
pub fn NucCallParams::nuc_set_peak_threshold(self : NucCallParams, threshold : Double) -> NucCallParams {
  NucCallParams::{ smooth_window: self.smooth_window, peak_threshold: threshold, min_distance: self.min_distance, max_distance: self.max_distance, occupancy_threshold: self.occupancy_threshold, width_threshold: self.width_threshold }
}

/// Set minimum distance
pub fn NucCallParams::nuc_set_min_distance(self : NucCallParams, dist : Int) -> NucCallParams {
  NucCallParams::{ smooth_window: self.smooth_window, peak_threshold: self.peak_threshold, min_distance: dist, max_distance: self.max_distance, occupancy_threshold: self.occupancy_threshold, width_threshold: self.width_threshold }
}

/// Set occupancy threshold
pub fn NucCallParams::nuc_set_occupancy_threshold(self : NucCallParams, threshold : Double) -> NucCallParams {
  NucCallParams::{ smooth_window: self.smooth_window, peak_threshold: self.peak_threshold, min_distance: self.min_distance, max_distance: self.max_distance, occupancy_threshold: threshold, width_threshold: self.width_threshold }
}

/// Signal track for nucleosome analysis
pub struct NucSignalTrack {
  /// Chromosome
  chrom : String
  /// Signal values
  signal : Array[Double]
  /// Start position
  start : Int
  /// Step size
  step : Int
  /// Length
  length : Int
}

/// Create new signal track
pub fn NucSignalTrack::new(chrom : String, signal : Array[Double], start : Int, step : Int) -> NucSignalTrack {
  NucSignalTrack::{ chrom, signal, start, step, length: signal.length() }
}

/// Get signal value at position
pub fn NucSignalTrack::nuc_get_signal(self : NucSignalTrack, pos : Int) -> Double {
  let idx = (pos - self.start) / self.step
  if idx >= 0 && idx < self.length {
    self.signal[idx]
  } else {
    0.0
  }
}

/// Get position for index
pub fn NucSignalTrack::nuc_get_position(self : NucSignalTrack, idx : Int) -> Int {
  self.start + idx * self.step
}

/// Get length
pub fn NucSignalTrack::nuc_get_length(self : NucSignalTrack) -> Int {
  self.length
}

/// Get chromosome
pub fn NucSignalTrack::nuc_get_chrom(self : NucSignalTrack) -> String {
  self.chrom
}

/// Smooth signal track
pub fn nuc_smooth_signal(track : NucSignalTrack, window : Int) -> NucSignalTrack {
  let n = track.signal.length()
  let smoothed : Array[Double] = Array::make(n, 0.0)

  let half_window = window / 2

  let mut i = 0
  while i < n {
    let start = if i - half_window > 0 { i - half_window } else { 0 }
    let end = if i + half_window < n { i + half_window } else { n - 1 }

    let mut sum = 0.0
    let mut count = 0
    let mut j = start
    while j <= end {
      sum = sum + track.signal[j]
      count = count + 1
      j = j + 1
    }

    smoothed[i] = if count > 0 { sum / count.to_double() } else { 0.0 }
    i = i + 1
  }

  NucSignalTrack::new(track.chrom, smoothed, track.start, track.step)
}

/// Find peaks in signal
pub fn nuc_find_peaks(track : NucSignalTrack, threshold : Double, min_distance : Int) -> Array[Int] {
  let n = track.signal.length()
  let peaks : Array[Int] = Array::new()

  let mut i = 1
  while i < n - 1 {
    let val = track.signal[i]
    let prev = track.signal[i - 1]
    let next = track.signal[i + 1]

    // Peak detection
    if val > threshold && val >= prev && val >= next {
      let pos = track.start + i * track.step

      // Check minimum distance from last peak
      if peaks.length() == 0 || (pos - track.start) - peaks[peaks.length() - 1] >= min_distance {
        peaks.push(pos - track.start)
      } else if val > track.signal[peaks[peaks.length() - 1]] {
        // Replace if higher
        peaks[peaks.length() - 1] = i
      }
    }
    i = i + 1
  }

  peaks
}

/// Calculate occupancy from signal
pub fn nuc_calculate_occupancy(track : NucSignalTrack,
                                start_pos : Int, end_pos : Int) -> Double {
  let n = track.signal.length()
  let start_idx = if start_pos >= track.start { (start_pos - track.start) / track.step } else { 0 }
  let end_idx = if end_pos > track.start {
    let idx = (end_pos - track.start) / track.step
    if idx < n { idx } else { n }
  } else { 0 }

  if start_idx >= end_idx {
    return 0.0
  }

  let mut sum = 0.0
  let mut max_val = 0.0
  let mut i = start_idx
  while i < end_idx {
    sum = sum + track.signal[i]
    if track.signal[i] > max_val {
      max_val = track.signal[i]
    }
    i = i + 1
  }

  let n_positions = (end_idx - start_idx).to_double()
  let mean = if n_positions > 0.0 { sum / n_positions } else { 0.0 }

  // Normalize by maximum possible occupancy
  if max_val > 0.0 {
    mean / max_val
  } else {
    0.0
  }
}

/// Call nucleosomes from signal
pub fn nuc_call_nucleosomes(track : NucSignalTrack, params : NucCallParams) -> NucCallResult {
  // Smooth the signal
  let smoothed = nuc_smooth_signal(track, params.smooth_window)

  // Find peaks
  let peak_indices = nuc_find_peaks(smoothed, params.peak_threshold, params.min_distance)

  // Convert peaks to nucleosome positions
  let nucleosomes : Array[NucPosition] = Array::new()

  let mut i = 0
  while i < peak_indices.length() {
    let idx = peak_indices[i]
    let center_pos = track.start + idx * track.step

    // Estimate width from peak shape
    let width = estimate_peak_width(track, idx, params.occupancy_threshold)

    // Occupancy
    let start_pos = center_pos - (width / 2.0).to_int()
    let end_pos = center_pos + (width / 2.0).to_int()
    let occupancy = nuc_calculate_occupancy(track, start_pos, end_pos)

    // Positioning score (narrow = well-positioned)
    let positioning_score = if width > 0.0 { 1.0 / width } else { 0.0 }

    let well_positioned = width < params.width_threshold

    nucleosomes.push(NucPosition::new(
      track.chrom,
      start_pos,
      end_pos,
      center_pos.to_double(),
      width,
      occupancy,
      positioning_score,
      well_positioned
    ))

    i = i + 1
  }

  // Calculate statistics
  let n_nucleosomes = nucleosomes.length()
  let mean_spacing = calculate_mean_spacing(nucleosomes)
  let mean_occupancy = calculate_mean_occupancy(nucleosomes)
  let frac_well_positioned = calculate_frac_well_positioned(nucleosomes)

  NucCallResult::new(nucleosomes, n_nucleosomes, mean_spacing,
                      mean_occupancy, frac_well_positioned)
}

/// Estimate peak width
fn estimate_peak_width(track : NucSignalTrack, center_idx : Int, threshold : Double) -> Double {
  let n = track.signal.length()
  let max_val = track.signal[center_idx]
  let half_max = max_val * threshold

  // Find left boundary
  let mut left = center_idx
  while left > 0 && track.signal[left] > half_max {
    left = left - 1
  }

  // Find right boundary
  let mut right = center_idx
  while right < n - 1 && track.signal[right] > half_max {
    right = right + 1
  }

  ((right - left) * track.step).to_double()
}

/// Calculate mean spacing between nucleosomes
fn calculate_mean_spacing(nucleosomes : Array[NucPosition]) -> Double {
  if nucleosomes.length() < 2 {
    return 0.0
  }

  // Sort by position
  let sorted : Array[NucPosition] = Array::new()
  let mut i = 0
  while i < nucleosomes.length() {
    sorted.push(nucleosomes[i])
    i = i + 1
  }

  // Simple bubble sort
  i = 0
  while i < sorted.length() {
    let mut j = i + 1
    while j < sorted.length() {
      if sorted[j].center < sorted[i].center {
        let tmp = sorted[i]
        sorted[i] = sorted[j]
        sorted[j] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }

  let mut total = 0.0
  let mut count = 0
  i = 1
  while i < sorted.length() {
    total = total + (sorted[i].center - sorted[i - 1].center)
    count = count + 1
    i = i + 1
  }

  if count > 0 { total / count.to_double() } else { 0.0 }
}

/// Calculate mean occupancy
fn calculate_mean_occupancy(nucleosomes : Array[NucPosition]) -> Double {
  if nucleosomes.length() == 0 {
    return 0.0
  }

  let mut sum = 0.0
  let mut i = 0
  while i < nucleosomes.length() {
    sum = sum + nucleosomes[i].occupancy
    i = i + 1
  }

  sum / nucleosomes.length().to_double()
}

/// Calculate fraction well-positioned
fn calculate_frac_well_positioned(nucleosomes : Array[NucPosition]) -> Double {
  if nucleosomes.length() == 0 {
    return 0.0
  }

  let mut count = 0
  let mut i = 0
  while i < nucleosomes.length() {
    if nucleosomes[i].well_positioned {
      count = count + 1
    }
    i = i + 1
  }

  count.to_double() / nucleosomes.length().to_double()
}

/// Compare nucleosome positioning between two conditions
pub fn nuc_compare_positioning(result1 : NucCallResult, result2 : NucCallResult,
                                tolerance : Int) -> NucDynamicResult {
  let shared : Array[NucPosition] = Array::new()
  let gained : Array[NucPosition] = Array::new()
  let lost : Array[NucPosition] = Array::new()

  // Find shared nucleosomes (within tolerance)
  let mut i = 0
  while i < result1.nucleosomes.length() {
    let nuc1 = result1.nucleosomes[i]
    let mut found = false

    let mut j = 0
    while j < result2.nucleosomes.length() {
      let nuc2 = result2.nucleosomes[j]
      if nuc1.chrom == nuc2.chrom &&
         (nuc1.center - nuc2.center).abs() <= tolerance.to_double() {
        shared.push(nuc1)
        found = true
        break
      }
      j = j + 1
    }

    if !found {
      lost.push(nuc1)
    }

    i = i + 1
  }

  // Find gained nucleosomes
  let mut k = 0
  while k < result2.nucleosomes.length() {
    let nuc2 = result2.nucleosomes[k]
    let mut found = false

    let mut m = 0
    while m < result1.nucleosomes.length() {
      let nuc1 = result1.nucleosomes[m]
      if nuc1.chrom == nuc2.chrom &&
         (nuc1.center - nuc2.center).abs() <= tolerance.to_double() {
        found = true
        break
      }
      m = m + 1
    }

    if !found {
      gained.push(nuc2)
    }

    k = k + 1
  }

  let total = result1.n_nucleosomes + result2.n_nucleosomes
  let frac_changed = if total > 0 {
    (gained.length() + lost.length()).to_double() / total.to_double()
  } else { 0.0 }

  let direction = if gained.length() > lost.length() { 1 }
                   else if lost.length() > gained.length() { -1 }
                   else { 0 }

  NucDynamicResult::new(shared, gained, lost, frac_changed, direction)
}

/// Create example signal track for testing
pub fn nuc_create_example_track() -> NucSignalTrack {
  let length = 1000
  let signal : Array[Double] = Array::make(length, 0.0)
  let step = 1
  let start = 0

  // Create synthetic nucleosome signal
  let positions = [50, 150, 250, 350, 450, 550, 650, 750, 850]
  let amplitudes = [0.8, 0.9, 0.7, 0.95, 0.85, 0.75, 0.92, 0.78, 0.88]
  let widths = [30.0, 28.0, 35.0, 25.0, 32.0, 38.0, 26.0, 34.0, 30.0]

  let mut i = 0
  while i < positions.length() {
    let pos = positions[i]
    let amp = amplitudes[i]
    let width = widths[i]

    let mut j = 0
    while j < length {
      let dist = (j.to_double() - pos.to_double()) / width * 3.0
      if dist.abs() < 3.0 {
        // Gaussian-like shape
        let val = amp * @math.exp(-dist * dist)
        if val > signal[j] {
          signal[j] = val
        }
      }
      j = j + 1
    }

    i = i + 1
  }

  // Add some noise
  let mut k = 0
  while k < length {
    let noise = 0.05 * @math.sin(k.to_double() * 0.1)
    signal[k] = signal[k] + noise
    k = k + 1
  }

  NucSignalTrack::new("chr1", signal, start, step)
}

/// Create example nucleosome call result
pub fn nuc_create_example_result() -> NucCallResult {
  let track = nuc_create_example_track()
  let params = NucCallParams::new()
  nuc_call_nucleosomes(track, params)
}