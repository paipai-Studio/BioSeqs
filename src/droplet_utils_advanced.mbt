///|
/// Advanced, portable implementation of the core DropletUtils 1.33.0
/// emptyDrops workflow. Matrices use the Bioconductor feature x barcode
/// convention.

///|
pub suberror DropletUtilsAdvancedError {
  DropletUtilsAdvancedError(String)
}

///|
pub struct DropletCountMatrix {
  counts : Array[Array[Int]]
  feature_names : Array[String]
  barcodes : Array[String]
} derive(Debug)
pub extend DropletCountMatrix with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct DropletUtilsConfig {
  lower : Int
  iterations : Int
  fdr_threshold : Double
  retain : Double
  by_rank : Int
  ignore : Int
  test_ambient : Bool
  alpha : Double
  estimate_alpha : Bool
  rank_exclude : Int
  rank_window : Double
  gradient_threshold : Double
  round_counts : Bool
  seed : Int
} derive(Debug)
pub extend DropletUtilsConfig with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct DropletAmbientProfile {
  raw_counts : Array[Int]
  proportions : Array[Double]
  assumed_empty : Array[Bool]
  lower : Int
  total : Int
} derive(Debug)
pub extend DropletAmbientProfile with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct DropletBarcodeRanks {
  barcodes : Array[String]
  ranks : Array[Double]
  totals : Array[Int]
  knee : Double
  inflection : Double
} derive(Debug)
pub extend DropletBarcodeRanks with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct DropletEmptyDropsResult {
  barcodes : Array[String]
  totals : Array[Int]
  log_probabilities : Array[Double?]
  p_values : Array[Double?]
  limited : Array[Bool?]
  fdr : Array[Double?]
  calls : Array[Bool]
  always_retained : Array[Bool]
  ambient : DropletAmbientProfile
  lower : Int
  retain : Double
  alpha : Double
  alpha_estimated : Bool
  iterations : Int
  fdr_threshold : Double
} derive(Debug)
pub extend DropletEmptyDropsResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct DropletUtilsSingleCellExperimentOutput {
  experiment : SingleCellExperiment
  result : DropletEmptyDropsResult
}

///|
fn dua_fail(message : String) -> Unit raise DropletUtilsAdvancedError {
  raise DropletUtilsAdvancedError::DropletUtilsAdvancedError(message)
}

///|
fn dua_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn dua_copy_int_matrix(values : Array[Array[Int]]) -> Array[Array[Int]] {
  let output : Array[Array[Int]] = []
  for row in values {
    output.push(row.copy())
  }
  output
}

///|
fn dua_copy_double_matrix(
  values : Array[Array[Double]],
) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in values {
    output.push(row.copy())
  }
  output
}

///|
fn dua_default_names(prefix : String, count : Int) -> Array[String] {
  let output : Array[String] = []
  for index in 0..<count {
    output.push(prefix + (index + 1).to_string())
  }
  output
}

///|
fn dua_prepare_names(
  supplied : Array[String],
  count : Int,
  prefix : String,
) -> Array[String] raise DropletUtilsAdvancedError {
  let names = if supplied.length() == 0 {
    dua_default_names(prefix, count)
  } else {
    if supplied.length() != count {
      dua_fail(
        "DropletUtils " + prefix + " names must match the matrix dimension",
      )
    }
    supplied.copy()
  }
  let seen : Map[String, Bool] = Map([])
  for name in names {
    if name == "" {
      dua_fail("DropletUtils " + prefix + " names cannot be empty")
    }
    if seen.contains(name) {
      dua_fail("DropletUtils " + prefix + " names must be unique")
    }
    seen[name] = true
  }
  names
}

///|
pub fn DropletCountMatrix::create(
  counts : Array[Array[Int]],
  feature_names? : Array[String] = [],
  barcodes? : Array[String] = [],
) -> DropletCountMatrix raise DropletUtilsAdvancedError {
  if counts.length() == 0 {
    dua_fail("DropletUtils count matrix must contain at least one feature")
  }
  if counts[0].length() == 0 {
    dua_fail("DropletUtils count matrix must contain at least one barcode")
  }
  let barcode_count = counts[0].length()
  let copied : Array[Array[Int]] = []
  for row in counts {
    if row.length() != barcode_count {
      dua_fail("DropletUtils count matrix must be rectangular")
    }
    let output_row : Array[Int] = []
    for value in row {
      if value < 0 {
        dua_fail("DropletUtils counts must be non-negative integers")
      }
      output_row.push(value)
    }
    copied.push(output_row)
  }
  DropletCountMatrix::{
    counts: copied,
    feature_names: dua_prepare_names(feature_names, counts.length(), "feature"),
    barcodes: dua_prepare_names(barcodes, barcode_count, "barcode"),
  }
}

///|
pub fn DropletCountMatrix::n_features(self : DropletCountMatrix) -> Int {
  self.counts.length()
}

///|
pub fn DropletCountMatrix::n_barcodes(self : DropletCountMatrix) -> Int {
  if self.counts.length() == 0 {
    0
  } else {
    self.counts[0].length()
  }
}

///|
pub fn DropletCountMatrix::copy_counts(
  self : DropletCountMatrix,
) -> Array[Array[Int]] {
  dua_copy_int_matrix(self.counts)
}

///|
pub fn DropletCountMatrix::barcode_totals(
  self : DropletCountMatrix,
) -> Array[Int] {
  let totals = Array::make(self.n_barcodes(), 0)
  for feature in 0..<self.n_features() {
    for barcode in 0..<self.n_barcodes() {
      totals[barcode] = totals[barcode] + self.counts[feature][barcode]
    }
  }
  totals
}

///|
pub fn DropletCountMatrix::feature_totals(
  self : DropletCountMatrix,
) -> Array[Int] {
  let totals = Array::make(self.n_features(), 0)
  for feature in 0..<self.n_features() {
    let mut total = 0
    for barcode in 0..<self.n_barcodes() {
      total = total + self.counts[feature][barcode]
    }
    totals[feature] = total
  }
  totals
}

///|
pub fn DropletUtilsConfig::create(
  lower? : Int = 100,
  iterations? : Int = 10000,
  fdr_threshold? : Double = 0.001,
  retain? : Double = -1.0,
  by_rank? : Int = 0,
  ignore? : Int = -1,
  test_ambient? : Bool = false,
  alpha? : Double = -1.0,
  estimate_alpha? : Bool = false,
  rank_exclude? : Int = 50,
  rank_window? : Double = 1.0,
  gradient_threshold? : Double = -1.0,
  round_counts? : Bool = true,
  seed? : Int = 1,
) -> DropletUtilsConfig raise DropletUtilsAdvancedError {
  if lower < 0 {
    dua_fail("DropletUtils lower threshold must be non-negative")
  }
  if iterations <= 0 {
    dua_fail("DropletUtils iterations must be positive")
  }
  if !dua_is_finite(fdr_threshold) ||
    fdr_threshold <= 0.0 ||
    fdr_threshold > 1.0 {
    dua_fail("DropletUtils FDR threshold must be in (0, 1]")
  }
  if !dua_is_finite(retain) || (retain != -1.0 && retain <= 0.0) {
    dua_fail("DropletUtils retain must be -1 for automatic or positive")
  }
  if by_rank < 0 {
    dua_fail("DropletUtils by_rank must be non-negative")
  }
  if ignore < -1 {
    dua_fail("DropletUtils ignore must be -1 or non-negative")
  }
  if !dua_is_finite(alpha) || alpha == 0.0 {
    dua_fail(
      "DropletUtils alpha must be negative for multinomial or positive for Dirichlet-multinomial",
    )
  }
  if rank_exclude < 0 {
    dua_fail("DropletUtils rank exclusion must be non-negative")
  }
  if !dua_is_finite(rank_window) || rank_window <= 0.0 {
    dua_fail("DropletUtils rank window must be positive")
  }
  if !dua_is_finite(gradient_threshold) {
    dua_fail("DropletUtils gradient threshold must be finite")
  }
  if seed <= 0 {
    dua_fail("DropletUtils random seed must be positive")
  }
  DropletUtilsConfig::{
    lower,
    iterations,
    fdr_threshold,
    retain,
    by_rank,
    ignore,
    test_ambient,
    alpha,
    estimate_alpha,
    rank_exclude,
    rank_window,
    gradient_threshold,
    round_counts,
    seed,
  }
}

///|
pub fn DropletUtilsConfig::default() -> DropletUtilsConfig {
  DropletUtilsConfig::{
    lower: 100,
    iterations: 10000,
    fdr_threshold: 0.001,
    retain: -1.0,
    by_rank: 0,
    ignore: -1,
    test_ambient: false,
    alpha: -1.0,
    estimate_alpha: false,
    rank_exclude: 50,
    rank_window: 1.0,
    gradient_threshold: -1.0,
    round_counts: true,
    seed: 1,
  }
}

///|
fn dua_log_gamma(value : Double) -> Double {
  let coefficients = [
    0.9999999999998099, 676.5203681218851, -1259.1392167224028, 771.3234287776531,
    -176.6150291621406, 12.507343278686905, -0.13857109526572012, 9.984369578019572e-6,
    1.5056327351493116e-7,
  ]
  if value < 0.5 {
    let pi = 3.141592653589793
    return @math.ln(pi) -
      @math.ln(@math.sin(pi * value).abs()) -
      dua_log_gamma(1.0 - value)
  }
  let shifted = value - 1.0
  let mut series = coefficients[0]
  for index in 1..<coefficients.length() {
    series = series + coefficients[index] / (shifted + index.to_double())
  }
  let t = shifted + 7.5
  0.9189385332046727 + (shifted + 0.5) * @math.ln(t) - t + @math.ln(series)
}

///|
fn dua_log_factorial(value : Int) -> Double {
  dua_log_gamma(value.to_double() + 1.0)
}

///|
pub fn droplet_utils_good_turing(
  counts : Array[Int],
) -> Array[Double] raise DropletUtilsAdvancedError {
  if counts.length() == 0 {
    dua_fail("DropletUtils Good-Turing input cannot be empty")
  }
  let positive : Array[Int] = []
  let mut zero_count = 0
  let mut grand_total = 0
  for count in counts {
    if count < 0 {
      dua_fail("DropletUtils Good-Turing counts must be non-negative")
    }
    grand_total = grand_total + count
    if count == 0 {
      zero_count = zero_count + 1
    } else {
      positive.push(count)
    }
  }
  if grand_total == 0 {
    dua_fail("DropletUtils Good-Turing requires at least one count")
  }
  positive.sort_by(fn(left : Int, right : Int) -> Int { left - right })
  let observed : Array[Int] = []
  let frequencies : Array[Int] = []
  for count in positive {
    if observed.length() == 0 || observed[observed.length() - 1] != count {
      observed.push(count)
      frequencies.push(1)
    } else {
      let last = frequencies.length() - 1
      frequencies[last] = frequencies[last] + 1
    }
  }
  let point_count = observed.length()
  let log_observed : Array[Double] = []
  let log_z : Array[Double] = []
  for index in 0..<point_count {
    log_observed.push(@math.ln(observed[index].to_double()))
    let previous = if index == 0 { 0 } else { observed[index - 1] }
    let width = if index + 1 == point_count {
      2 * (observed[index] - previous)
    } else {
      observed[index + 1] - previous
    }
    log_z.push(
      @math.ln(2.0 * frequencies[index].to_double()) -
      @math.ln(width.to_double()),
    )
  }
  let mut slope = -1.0
  if point_count > 1 {
    let mut mean_x = 0.0
    let mut mean_y = 0.0
    for index in 0..<point_count {
      mean_x = mean_x + log_observed[index]
      mean_y = mean_y + log_z[index]
    }
    mean_x = mean_x / point_count.to_double()
    mean_y = mean_y / point_count.to_double()
    let mut cross = 0.0
    let mut squares = 0.0
    for index in 0..<point_count {
      let dx = log_observed[index] - mean_x
      cross = cross + dx * (log_z[index] - mean_y)
      squares = squares + dx * dx
    }
    if squares > 1.0e-15 {
      slope = cross / squares
    }
  }
  let p_zero = if observed[0] == 1 {
    frequencies[0].to_double() / grand_total.to_double()
  } else {
    0.0
  }
  let adjusted : Array[Double] = Array::make(point_count, 0.0)
  let mut use_smoothed = false
  let mut adjusted_total = 0.0
  for index in 0..<point_count {
    let next_count = observed[index] + 1
    let smoothed = next_count.to_double() *
      @math.exp(
        slope * (@math.ln(next_count.to_double()) - log_observed[index]),
      )
    if index + 1 == point_count || observed[index + 1] != next_count {
      use_smoothed = true
    }
    if !use_smoothed {
      let raw = next_count.to_double() *
        frequencies[index + 1].to_double() /
        frequencies[index].to_double()
      let uncertainty = 1.96 *
        raw *
        (1.0 / frequencies[index + 1].to_double() +
        1.0 / frequencies[index].to_double()).sqrt()
      if (raw - smoothed).abs() <= uncertainty {
        use_smoothed = true
      } else {
        adjusted[index] = raw
      }
    }
    if use_smoothed {
      adjusted[index] = smoothed
    }
    adjusted_total = adjusted_total +
      adjusted[index] * frequencies[index].to_double()
  }
  let scale = if adjusted_total > 0.0 {
    (1.0 - p_zero) / adjusted_total
  } else {
    0.0
  }
  let lookup : Map[Int, Double] = Map([])
  for index in 0..<point_count {
    lookup[observed[index]] = adjusted[index] * scale
  }
  let output = Array::make(counts.length(), 0.0)
  let zero_probability = if zero_count > 0 {
    p_zero / zero_count.to_double()
  } else {
    0.0
  }
  for index in 0..<counts.length() {
    output[index] = if counts[index] == 0 {
      zero_probability
    } else {
      lookup[counts[index]]
    }
  }
  let mut still_zero = 0
  for value in output {
    if value <= 0.0 {
      still_zero = still_zero + 1
    }
  }
  if still_zero > 0 {
    let pseudo_probability = (1.0 / grand_total.to_double()).min(0.5)
    for index in 0..<output.length() {
      if output[index] <= 0.0 {
        output[index] = pseudo_probability / still_zero.to_double()
      } else {
        output[index] = output[index] * (1.0 - pseudo_probability)
      }
    }
  }
  let mut total_probability = 0.0
  for value in output {
    total_probability = total_probability + value
  }
  if total_probability <= 0.0 {
    dua_fail("DropletUtils Good-Turing produced an invalid profile")
  }
  for index in 0..<output.length() {
    output[index] = output[index] / total_probability
  }
  output
}

///|
fn dua_empty_mask(
  totals : Array[Int],
  lower : Int,
  by_rank : Int,
  known_empty : Array[Bool],
) -> (Array[Bool], Int) raise DropletUtilsAdvancedError {
  if known_empty.length() > 0 {
    if known_empty.length() != totals.length() {
      dua_fail(
        "DropletUtils known-empty mask must match the number of barcodes",
      )
    }
    let output = known_empty.copy()
    let mut found = false
    for value in output {
      if value {
        found = true
      }
    }
    if !found {
      dua_fail("DropletUtils known-empty mask must select at least one barcode")
    }
    return (output, lower)
  }
  if by_rank > 0 {
    if by_rank >= totals.length() {
      dua_fail("DropletUtils by_rank must be smaller than barcode count")
    }
    let order : Array[Int] = []
    for index in 0..<totals.length() {
      order.push(index)
    }
    order.sort_by(fn(left : Int, right : Int) -> Int {
      if totals[left] > totals[right] {
        -1
      } else if totals[left] < totals[right] {
        1
      } else {
        left - right
      }
    })
    let output = Array::make(totals.length(), false)
    let mut actual_lower = 0
    for rank in by_rank..<order.length() {
      let index = order[rank]
      output[index] = true
      if totals[index] > actual_lower {
        actual_lower = totals[index]
      }
    }
    return (output, actual_lower)
  }
  let output = Array::make(totals.length(), false)
  for index in 0..<totals.length() {
    output[index] = totals[index] <= lower
  }
  (output, lower)
}

///|
pub fn droplet_utils_ambient_profile(
  data : DropletCountMatrix,
  lower? : Int = 100,
  by_rank? : Int = 0,
  known_empty? : Array[Bool] = [],
) -> DropletAmbientProfile raise DropletUtilsAdvancedError {
  if lower < 0 || by_rank < 0 {
    dua_fail("DropletUtils ambient thresholds must be non-negative")
  }
  let totals = data.barcode_totals()
  let (empty, actual_lower) = dua_empty_mask(
    totals, lower, by_rank, known_empty,
  )
  let raw = Array::make(data.n_features(), 0)
  let all_feature_totals = data.feature_totals()
  let mut ambient_total = 0
  for feature in 0..<data.n_features() {
    let mut value = 0
    for barcode in 0..<data.n_barcodes() {
      if empty[barcode] {
        value = value + data.counts[feature][barcode]
      }
    }
    raw[feature] = value
    ambient_total = ambient_total + value
  }
  if ambient_total == 0 {
    dua_fail("DropletUtils has no counts for estimating the ambient profile")
  }
  let active_indices : Array[Int] = []
  let active_counts : Array[Int] = []
  for feature in 0..<data.n_features() {
    if all_feature_totals[feature] > 0 {
      active_indices.push(feature)
      active_counts.push(raw[feature])
    }
  }
  let active_probabilities = droplet_utils_good_turing(active_counts)
  let proportions = Array::make(data.n_features(), 0.0)
  for index in 0..<active_indices.length() {
    proportions[active_indices[index]] = active_probabilities[index]
  }
  DropletAmbientProfile::{
    raw_counts: raw,
    proportions,
    assumed_empty: empty,
    lower: actual_lower,
    total: ambient_total,
  }
}

///|
pub fn DropletAmbientProfile::n_empty(self : DropletAmbientProfile) -> Int {
  let mut count = 0
  for value in self.assumed_empty {
    if value {
      count = count + 1
    }
  }
  count
}

///|
fn dua_log10(value : Double) -> Double {
  @math.ln(value) / 2.302585092994046
}

///|
fn dua_pow10(value : Double) -> Double {
  @math.exp(value * 2.302585092994046)
}

///|
fn dua_interpolate_curve(
  target : Double,
  cumulative : Array[Double],
  step : Array[Double],
  x : Array[Double],
  y : Array[Double],
) -> (Double, Double) {
  let mut left = 0
  while left + 1 < cumulative.length() && cumulative[left + 1] < target {
    left = left + 1
  }
  if left + 1 >= cumulative.length() {
    return (x[x.length() - 1], y[y.length() - 1])
  }
  let proportion = if step[left + 1] > 0.0 {
    (target - cumulative[left]) / step[left + 1]
  } else {
    0.0
  }
  (
    x[left] + proportion * (x[left + 1] - x[left]),
    y[left] + proportion * (y[left + 1] - y[left]),
  )
}

///|
pub fn droplet_utils_barcode_ranks(
  data : DropletCountMatrix,
  lower? : Int = 100,
  exclude_from? : Int = 50,
  window? : Double = 1.0,
  gradient_threshold? : Double = -1.0,
) -> DropletBarcodeRanks raise DropletUtilsAdvancedError {
  if lower < 0 || exclude_from < 0 {
    dua_fail("DropletUtils barcode-rank thresholds must be non-negative")
  }
  if !dua_is_finite(window) || window <= 0.0 {
    dua_fail("DropletUtils barcode-rank window must be positive")
  }
  if !dua_is_finite(gradient_threshold) {
    dua_fail("DropletUtils barcode-rank gradient must be finite")
  }
  let totals = data.barcode_totals()
  let order : Array[Int] = []
  for index in 0..<totals.length() {
    order.push(index)
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    if totals[left] > totals[right] {
      -1
    } else if totals[left] < totals[right] {
      1
    } else {
      left - right
    }
  })
  let ranks = Array::make(totals.length(), 0.0)
  let run_totals : Array[Int] = []
  let run_ranks : Array[Double] = []
  let mut position = 0
  while position < order.length() {
    let start = position
    let current_total = totals[order[position]]
    while position < order.length() && totals[order[position]] == current_total {
      position = position + 1
    }
    let length = position - start
    let mid_rank = start.to_double() + (length + 1).to_double() / 2.0
    run_totals.push(current_total)
    run_ranks.push(mid_rank)
    for sorted_position in start..<position {
      ranks[order[sorted_position]] = mid_rank
    }
  }
  let x : Array[Double] = []
  let y : Array[Double] = []
  for run in 0..<run_totals.length() {
    if run_totals[run] > lower && run_ranks[run] > exclude_from.to_double() {
      x.push(dua_log10(run_ranks[run]))
      y.push(dua_log10(run_totals[run].to_double()))
    }
  }
  if x.length() < 2 {
    dua_fail(
      "DropletUtils has insufficient unique points for knee/inflection detection",
    )
  }
  let step = Array::make(x.length(), 0.0)
  let cumulative = Array::make(x.length(), 0.0)
  for index in 1..<x.length() {
    let dx = x[index] - x[index - 1]
    let dy = y[index] - y[index - 1]
    step[index] = (dx * dx + dy * dy).sqrt()
    cumulative[index] = cumulative[index - 1] + step[index]
  }
  let scan : Array[Int] = []
  for index in 0..<x.length() {
    if cumulative[index] + window <= cumulative[cumulative.length() - 1] {
      scan.push(index)
    }
  }
  let mut knee = dua_pow10(y[y.length() - 1])
  let mut inflection = knee
  if scan.length() > 0 {
    let left_x : Array[Double] = []
    let right_x : Array[Double] = []
    let right_y : Array[Double] = []
    let mid_x : Array[Double] = []
    let mid_y : Array[Double] = []
    let gaps : Array[Double] = []
    let gradients : Array[Double] = []
    let mid_above : Array[Bool] = []
    for scan_index in scan {
      let lx = x[scan_index]
      let ly = y[scan_index]
      let (rx, ry) = dua_interpolate_curve(
        cumulative[scan_index] + window,
        cumulative,
        step,
        x,
        y,
      )
      let (mx, my) = dua_interpolate_curve(
        cumulative[scan_index] + window / 2.0,
        cumulative,
        step,
        x,
        y,
      )
      let gradient = (ry - ly) / (rx - lx)
      let intercept = ry - gradient * rx
      let dx = lx - rx
      let dy = ly - ry
      left_x.push(lx)
      right_x.push(rx)
      right_y.push(ry)
      mid_x.push(mx)
      mid_y.push(my)
      gaps.push((dx * dx + dy * dy).sqrt())
      gradients.push(gradient)
      mid_above.push(my > gradient * mx + intercept)
    }
    let mut first_elbow = -1
    for index in 0..<scan.length() {
      if !mid_above[index] &&
        gradients[index] < gradient_threshold &&
        first_elbow < 0 {
        first_elbow = index
      }
    }
    let mut inflection_limit = scan.length()
    if first_elbow >= 0 {
      inflection_limit = 0
      while inflection_limit < left_x.length() &&
            left_x[inflection_limit] <= mid_x[first_elbow] {
        inflection_limit = inflection_limit + 1
      }
      if inflection_limit == 0 {
        inflection_limit = 1
      }
    }
    let mut inflection_index = 0
    for index in 1..<inflection_limit {
      if gradients[index] < gradients[inflection_index] {
        inflection_index = index
      }
    }
    let knee_limit = if first_elbow >= 0 {
      first_elbow + 1
    } else {
      scan.length()
    }
    let mut knee_index = -1
    for index in 0..<knee_limit {
      if mid_above[index] && (knee_index < 0 || gaps[index] < gaps[knee_index]) {
        knee_index = index
      }
    }
    if knee_index < 0 {
      knee_index = inflection_index
    }
    knee = dua_pow10(mid_y[knee_index])
    inflection = dua_pow10(mid_y[inflection_index])
  }
  DropletBarcodeRanks::{
    barcodes: data.barcodes.copy(),
    ranks,
    totals,
    knee,
    inflection,
  }
}

///|
fn dua_validate_probability_input(
  counts : Array[Int],
  proportions : Array[Double],
  alpha : Double,
) -> Unit raise DropletUtilsAdvancedError {
  if counts.length() == 0 || counts.length() != proportions.length() {
    dua_fail(
      "DropletUtils counts and ambient proportions must have equal non-zero length",
    )
  }
  if !dua_is_finite(alpha) || alpha == 0.0 {
    dua_fail("DropletUtils alpha must be negative or positive")
  }
  let mut total_probability = 0.0
  for index in 0..<counts.length() {
    if counts[index] < 0 {
      dua_fail("DropletUtils probability counts must be non-negative")
    }
    let probability = proportions[index]
    if !dua_is_finite(probability) || probability < 0.0 {
      dua_fail(
        "DropletUtils ambient proportions must be finite and non-negative",
      )
    }
    if probability == 0.0 && counts[index] > 0 {
      dua_fail(
        "DropletUtils positive counts cannot have zero ambient probability",
      )
    }
    total_probability = total_probability + probability
  }
  if (total_probability - 1.0).abs() > 1.0e-8 {
    dua_fail("DropletUtils ambient proportions must sum to one")
  }
}

///|
fn dua_probability_data_component(
  counts : Array[Int],
  proportions : Array[Double],
  alpha : Double,
) -> Double {
  let mut output = 0.0
  for index in 0..<counts.length() {
    let count = counts[index]
    if count > 0 {
      if alpha < 0.0 {
        output = output +
          count.to_double() * @math.ln(proportions[index]) -
          dua_log_factorial(count)
      } else {
        let scaled = alpha * proportions[index]
        output = output +
          dua_log_gamma(scaled + count.to_double()) -
          dua_log_factorial(count) -
          dua_log_gamma(scaled)
      }
    }
  }
  output
}

///|
pub fn droplet_utils_log_probability(
  counts : Array[Int],
  proportions : Array[Double],
  alpha? : Double = -1.0,
) -> Double raise DropletUtilsAdvancedError {
  dua_validate_probability_input(counts, proportions, alpha)
  let mut total = 0
  for count in counts {
    total = total + count
  }
  let data_component = dua_probability_data_component(
    counts, proportions, alpha,
  )
  if alpha < 0.0 {
    data_component + dua_log_factorial(total)
  } else {
    data_component +
    dua_log_factorial(total) +
    dua_log_gamma(alpha) -
    dua_log_gamma(total.to_double() + alpha)
  }
}

///|
priv struct DropletUtilsRng {
  mut state : Int
}

///|
fn DropletUtilsRng::new(seed : Int) -> DropletUtilsRng {
  let normalized = (seed - 1) % 2147483646 + 1
  DropletUtilsRng::{ state: normalized, }
}

///|
fn DropletUtilsRng::uniform(self : DropletUtilsRng) -> Double {
  let quotient = self.state / 44488
  let remainder = self.state % 44488
  let candidate = 48271 * remainder - 3399 * quotient
  self.state = if candidate > 0 { candidate } else { candidate + 2147483647 }
  self.state.to_double() / 2147483647.0
}

///|
fn DropletUtilsRng::normal(self : DropletUtilsRng) -> Double {
  let first = self.uniform().max(1.0e-15)
  let second = self.uniform()
  (-2.0 * @math.ln(first)).sqrt() * @math.cos(6.283185307179586 * second)
}

///|
fn DropletUtilsRng::gamma(self : DropletUtilsRng, shape : Double) -> Double {
  if shape < 1.0 {
    let uniform = self.uniform().max(1.0e-15)
    return self.gamma(shape + 1.0) * @math.exp(@math.ln(uniform) / shape)
  }
  let d = shape - 1.0 / 3.0
  let c = 1.0 / (9.0 * d).sqrt()
  while true {
    let normal = self.normal()
    let base = 1.0 + c * normal
    if base > 0.0 {
      let volume = base * base * base
      let uniform = self.uniform().max(1.0e-15)
      if uniform < 1.0 - 0.0331 * normal * normal * normal * normal {
        return d * volume
      }
      if @math.ln(uniform) <
        0.5 * normal * normal + d * (1.0 - volume + @math.ln(volume)) {
        return d * volume
      }
    }
  }
  shape
}

///|
fn DropletUtilsRng::categorical(
  self : DropletUtilsRng,
  weights : Array[Double],
) -> Int {
  let mut total = 0.0
  for weight in weights {
    total = total + weight
  }
  if total <= 0.0 {
    return 0
  }
  let target = self.uniform() * total
  let mut cumulative = 0.0
  for index in 0..<weights.length() {
    cumulative = cumulative + weights[index]
    if target <= cumulative {
      return index
    }
  }
  weights.length() - 1
}

///|
fn dua_monte_carlo_extreme(
  totals : Array[Int],
  observed : Array[Double],
  ambient : Array[Double],
  iterations : Int,
  alpha : Double,
  seed : Int,
) -> Array[Int] {
  let order : Array[Int] = []
  for index in 0..<totals.length() {
    order.push(index)
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    if totals[left] < totals[right] {
      -1
    } else if totals[left] > totals[right] {
      1
    } else if observed[left] < observed[right] {
      -1
    } else if observed[left] > observed[right] {
      1
    } else {
      left - right
    }
  })
  let extreme = Array::make(totals.length(), 0)
  let rng = DropletUtilsRng::new(seed)
  for _ in 0..<iterations {
    let weights = if alpha < 0.0 {
      ambient.copy()
    } else {
      let sampled = Array::make(ambient.length(), 0.0)
      let mut sampled_total = 0.0
      for feature in 0..<ambient.length() {
        if ambient[feature] > 0.0 {
          sampled[feature] = rng.gamma(ambient[feature] * alpha)
          sampled_total = sampled_total + sampled[feature]
        }
      }
      if sampled_total <= 0.0 {
        ambient.copy()
      } else {
        sampled
      }
    }
    let tracker = Array::make(ambient.length(), 0)
    let mut current_total = 0
    let mut current_probability = 0.0
    let mut position = 0
    while position < order.length() {
      let target_total = totals[order[position]]
      while current_total < target_total {
        let feature = rng.categorical(weights)
        let previous = tracker[feature]
        if alpha < 0.0 {
          current_probability = current_probability + @math.ln(ambient[feature])
        } else {
          current_probability = current_probability +
            @math.ln(alpha * ambient[feature] + previous.to_double())
        }
        tracker[feature] = previous + 1
        current_probability = current_probability -
          @math.ln(tracker[feature].to_double())
        current_total = current_total + 1
      }
      let mut end = position + 1
      while end < order.length() && totals[order[end]] == target_total {
        end = end + 1
      }
      for grouped in position..<end {
        let index = order[grouped]
        if current_probability <= observed[index] + 1.0e-12 {
          extreme[index] = extreme[index] + 1
        }
      }
      position = end
    }
  }
  extreme
}

///|
pub fn droplet_utils_monte_carlo_p_value(
  counts : Array[Int],
  proportions : Array[Double],
  iterations? : Int = 1000,
  alpha? : Double = -1.0,
  seed? : Int = 1,
) -> (Double, Bool) raise DropletUtilsAdvancedError {
  if iterations <= 0 || seed <= 0 {
    dua_fail("DropletUtils Monte Carlo iterations and seed must be positive")
  }
  dua_validate_probability_input(counts, proportions, alpha)
  let mut total = 0
  for value in counts {
    total = total + value
  }
  if total <= 0 {
    dua_fail("DropletUtils Monte Carlo count vector must have positive total")
  }
  let observed = dua_probability_data_component(counts, proportions, alpha)
  let extreme = dua_monte_carlo_extreme(
      [total],
      [observed],
      proportions,
      iterations,
      alpha,
      seed,
    )[0]
  ((extreme + 1).to_double() / (iterations + 1).to_double(), extreme == 0)
}

///|
fn dua_alpha_log_likelihood(
  data : DropletCountMatrix,
  empty : Array[Bool],
  proportions : Array[Double],
  alpha : Double,
) -> Double {
  let mut output = 0.0
  for barcode in 0..<data.n_barcodes() {
    if empty[barcode] {
      let mut total = 0
      for feature in 0..<data.n_features() {
        total = total + data.counts[feature][barcode]
      }
      output = output +
        dua_log_gamma(alpha) -
        dua_log_gamma(total.to_double() + alpha)
      for feature in 0..<data.n_features() {
        let count = data.counts[feature][barcode]
        if count > 0 && proportions[feature] > 0.0 {
          let scaled = alpha * proportions[feature]
          output = output +
            dua_log_gamma(count.to_double() + scaled) -
            dua_log_gamma(scaled)
        }
      }
    }
  }
  output
}

///|
pub fn droplet_utils_estimate_alpha(
  data : DropletCountMatrix,
  ambient : DropletAmbientProfile,
) -> Double raise DropletUtilsAdvancedError {
  if ambient.assumed_empty.length() != data.n_barcodes() ||
    ambient.proportions.length() != data.n_features() {
    dua_fail("DropletUtils ambient profile dimensions do not match counts")
  }
  let ratio = 0.6180339887498949
  let mut left = @math.ln(0.01)
  let mut right = @math.ln(10000.0)
  let mut first = right - ratio * (right - left)
  let mut second = left + ratio * (right - left)
  let mut first_value = dua_alpha_log_likelihood(
    data,
    ambient.assumed_empty,
    ambient.proportions,
    @math.exp(first),
  )
  let mut second_value = dua_alpha_log_likelihood(
    data,
    ambient.assumed_empty,
    ambient.proportions,
    @math.exp(second),
  )
  for _ in 0..<80 {
    if first_value < second_value {
      left = first
      first = second
      first_value = second_value
      second = left + ratio * (right - left)
      second_value = dua_alpha_log_likelihood(
        data,
        ambient.assumed_empty,
        ambient.proportions,
        @math.exp(second),
      )
    } else {
      right = second
      second = first
      second_value = first_value
      first = right - ratio * (right - left)
      first_value = dua_alpha_log_likelihood(
        data,
        ambient.assumed_empty,
        ambient.proportions,
        @math.exp(first),
      )
    }
  }
  @math.exp((left + right) / 2.0)
}

///|
fn dua_bh_adjust(values : Array[Double?]) -> Array[Double?] {
  let indexed : Array[(Double, Int)] = []
  for index in 0..<values.length() {
    match values[index] {
      Some(value) => indexed.push((value.max(0.0).min(1.0), index))
      None => ()
    }
  }
  indexed.sort_by(fn(left : (Double, Int), right : (Double, Int)) -> Int {
    if left.0 < right.0 {
      -1
    } else if left.0 > right.0 {
      1
    } else {
      left.1 - right.1
    }
  })
  let output : Array[Double?] = Array::make(values.length(), None)
  let count = indexed.length()
  let mut running = 1.0
  let mut position = count - 1
  while position >= 0 {
    let candidate = (indexed[position].0 *
    count.to_double() /
    (position + 1).to_double()).min(1.0)
    running = running.min(candidate)
    output[indexed[position].1] = Some(running)
    position = position - 1
  }
  output
}

///|
fn dua_barcode_counts(data : DropletCountMatrix, barcode : Int) -> Array[Int] {
  let output : Array[Int] = []
  for feature in 0..<data.n_features() {
    output.push(data.counts[feature][barcode])
  }
  output
}

///|
pub fn droplet_utils_empty_drops(
  data : DropletCountMatrix,
  known_empty? : Array[Bool] = [],
  config? : DropletUtilsConfig = DropletUtilsConfig::default(),
) -> DropletEmptyDropsResult raise DropletUtilsAdvancedError {
  let ambient = droplet_utils_ambient_profile(
    data,
    lower=config.lower,
    by_rank=config.by_rank,
    known_empty~,
  )
  let alpha = if config.estimate_alpha {
    droplet_utils_estimate_alpha(data, ambient)
  } else {
    config.alpha
  }
  let totals = data.barcode_totals()
  let tested_indices : Array[Int] = []
  let tested_totals : Array[Int] = []
  let observed_components : Array[Double] = []
  let log_probabilities : Array[Double?] = Array::make(data.n_barcodes(), None)
  for barcode in 0..<data.n_barcodes() {
    let eligible = if config.test_ambient {
      totals[barcode] > 0
    } else {
      !ambient.assumed_empty[barcode]
    }
    if eligible && (config.ignore < 0 || totals[barcode] > config.ignore) {
      let barcode_counts = dua_barcode_counts(data, barcode)
      dua_validate_probability_input(barcode_counts, ambient.proportions, alpha)
      let component = dua_probability_data_component(
        barcode_counts,
        ambient.proportions,
        alpha,
      )
      let rest = if alpha < 0.0 {
        dua_log_factorial(totals[barcode])
      } else {
        dua_log_factorial(totals[barcode]) +
        dua_log_gamma(alpha) -
        dua_log_gamma(totals[barcode].to_double() + alpha)
      }
      tested_indices.push(barcode)
      tested_totals.push(totals[barcode])
      observed_components.push(component)
      log_probabilities[barcode] = Some(component + rest)
    }
  }
  let extreme = dua_monte_carlo_extreme(
    tested_totals,
    observed_components,
    ambient.proportions,
    config.iterations,
    alpha,
    config.seed,
  )
  let p_values : Array[Double?] = Array::make(data.n_barcodes(), None)
  let limited : Array[Bool?] = Array::make(data.n_barcodes(), None)
  for tested in 0..<tested_indices.length() {
    let barcode = tested_indices[tested]
    p_values[barcode] = Some(
      (extreme[tested] + 1).to_double() / (config.iterations + 1).to_double(),
    )
    limited[barcode] = Some(extreme[tested] == 0)
  }
  let retain = if config.retain > 0.0 {
    config.retain
  } else {
    droplet_utils_barcode_ranks(
      data,
      lower=ambient.lower,
      exclude_from=config.rank_exclude,
      window=config.rank_window,
      gradient_threshold=config.gradient_threshold,
    ).knee
  }
  let always_retained = Array::make(data.n_barcodes(), false)
  let correction_values = p_values.copy()
  for barcode in 0..<data.n_barcodes() {
    if totals[barcode].to_double() >= retain {
      always_retained[barcode] = true
      correction_values[barcode] = Some(0.0)
    }
  }
  if config.test_ambient {
    for barcode in 0..<data.n_barcodes() {
      if ambient.assumed_empty[barcode] {
        correction_values[barcode] = None
      }
    }
  }
  let fdr = dua_bh_adjust(correction_values)
  let calls = Array::make(data.n_barcodes(), false)
  for barcode in 0..<data.n_barcodes() {
    match fdr[barcode] {
      Some(value) => if value <= config.fdr_threshold { calls[barcode] = true }
      None => ()
    }
  }
  DropletEmptyDropsResult::{
    barcodes: data.barcodes.copy(),
    totals,
    log_probabilities,
    p_values,
    limited,
    fdr,
    calls,
    always_retained,
    ambient,
    lower: ambient.lower,
    retain,
    alpha,
    alpha_estimated: config.estimate_alpha,
    iterations: config.iterations,
    fdr_threshold: config.fdr_threshold,
  }
}

///|
pub fn DropletEmptyDropsResult::n_called(self : DropletEmptyDropsResult) -> Int {
  let mut count = 0
  for call in self.calls {
    if call {
      count = count + 1
    }
  }
  count
}

///|
pub fn DropletEmptyDropsResult::n_tested(self : DropletEmptyDropsResult) -> Int {
  let mut count = 0
  for value in self.p_values {
    match value {
      Some(_) => count = count + 1
      None => ()
    }
  }
  count
}

///|
pub fn DropletEmptyDropsResult::index_of(
  self : DropletEmptyDropsResult,
  barcode : String,
) -> Int? {
  for index in 0..<self.barcodes.length() {
    if self.barcodes[index] == barcode {
      return Some(index)
    }
  }
  None
}

///|
pub fn DropletEmptyDropsResult::called_barcodes(
  self : DropletEmptyDropsResult,
) -> Array[String] {
  let output : Array[String] = []
  for index in 0..<self.barcodes.length() {
    if self.calls[index] {
      output.push(self.barcodes[index])
    }
  }
  output
}

///|
pub fn DropletEmptyDropsResult::summary(
  self : DropletEmptyDropsResult,
) -> String {
  "DropletUtils emptyDrops: " +
  self.barcodes.length().to_string() +
  " barcodes, " +
  self.ambient.n_empty().to_string() +
  " assumed empty, " +
  self.n_tested().to_string() +
  " tested, " +
  self.n_called().to_string() +
  " called; lower=" +
  self.lower.to_string() +
  ", retain=" +
  self.retain.to_string() +
  ", iterations=" +
  self.iterations.to_string()
}

///|
pub fn DropletEmptyDropsResult::filter(
  self : DropletEmptyDropsResult,
  data : DropletCountMatrix,
) -> DropletCountMatrix raise DropletUtilsAdvancedError {
  if self.barcodes.length() != data.n_barcodes() {
    dua_fail("DropletUtils result does not match the count matrix")
  }
  for index in 0..<self.barcodes.length() {
    if self.barcodes[index] != data.barcodes[index] {
      dua_fail("DropletUtils result barcodes do not match the count matrix")
    }
  }
  let selected_names : Array[String] = []
  for barcode in 0..<data.n_barcodes() {
    if self.calls[barcode] {
      selected_names.push(data.barcodes[barcode])
    }
  }
  if selected_names.length() == 0 {
    dua_fail("DropletUtils result contains no called barcodes")
  }
  let selected_counts : Array[Array[Int]] = []
  for feature in 0..<data.n_features() {
    let row : Array[Int] = []
    for barcode in 0..<data.n_barcodes() {
      if self.calls[barcode] {
        row.push(data.counts[feature][barcode])
      }
    }
    selected_counts.push(row)
  }
  DropletCountMatrix::create(
    selected_counts,
    feature_names=data.feature_names,
    barcodes=selected_names,
  )
}

///|
fn dua_copy_sce(source : SingleCellExperiment) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  for key in source.assays.keys() {
    assays[key] = dua_copy_double_matrix(source.assays[key])
  }
  let row_data : Map[String, Array[String]] = Map([])
  for key in source.row_data.keys() {
    row_data[key] = source.row_data[key].copy()
  }
  let col_data : Map[String, Array[String]] = Map([])
  for key in source.col_data.keys() {
    col_data[key] = source.col_data[key].copy()
  }
  let reduced_dims : Map[String, Array[Array[Double]]] = Map([])
  for key in source.reduced_dims.keys() {
    reduced_dims[key] = dua_copy_double_matrix(source.reduced_dims[key])
  }
  let alternatives : Map[String, SingleCellExperiment] = Map([])
  for key in source.alternative_experiments.keys() {
    alternatives[key] = source.alternative_experiments[key]
  }
  let metadata : Map[String, String] = Map([])
  for key in source.metadata.keys() {
    metadata[key] = source.metadata[key]
  }
  SingleCellExperiment::{
    assays,
    row_data,
    col_data,
    reduced_dims,
    alternative_experiments: alternatives,
    metadata,
    row_names: source.row_names.copy(),
    col_names: source.col_names.copy(),
  }
}

///|
fn dua_option_double_string(value : Double?) -> String {
  match value {
    Some(number) => number.to_string()
    None => "NA"
  }
}

///|
pub fn droplet_utils_empty_drops_sce(
  experiment : SingleCellExperiment,
  assay_name? : String = "counts",
  output_prefix? : String = "emptyDrops",
  known_empty_column? : String = "",
  config? : DropletUtilsConfig = DropletUtilsConfig::default(),
) -> DropletUtilsSingleCellExperimentOutput raise DropletUtilsAdvancedError {
  if assay_name == "" || output_prefix == "" {
    dua_fail("DropletUtils assay name and output prefix cannot be empty")
  }
  if !experiment.assays.contains(assay_name) {
    dua_fail("DropletUtils assay not found: " + assay_name)
  }
  let assay = experiment.assays[assay_name]
  if assay.length() == 0 || assay[0].length() == 0 {
    dua_fail("DropletUtils assay cannot be empty")
  }
  if experiment.row_names.length() != assay.length() ||
    experiment.col_names.length() != assay[0].length() {
    dua_fail(
      "DropletUtils SingleCellExperiment names must match assay dimensions",
    )
  }
  let counts : Array[Array[Int]] = []
  for row in assay {
    if row.length() != experiment.col_names.length() {
      dua_fail("DropletUtils SingleCellExperiment assay must be rectangular")
    }
    let output_row : Array[Int] = []
    for value in row {
      if !dua_is_finite(value) || value < 0.0 {
        dua_fail(
          "DropletUtils SingleCellExperiment counts must be finite and non-negative",
        )
      }
      if !config.round_counts && (value - value.round()).abs() > 1.0e-8 {
        dua_fail(
          "DropletUtils SingleCellExperiment counts must be integer-valued when rounding is disabled",
        )
      }
      output_row.push(value.round().to_int())
    }
    counts.push(output_row)
  }
  let data = DropletCountMatrix::create(
    counts,
    feature_names=experiment.row_names,
    barcodes=experiment.col_names,
  )
  let known_empty : Array[Bool] = if known_empty_column == "" {
    []
  } else {
    if !experiment.col_data.contains(known_empty_column) {
      dua_fail(
        "DropletUtils known-empty column not found: " + known_empty_column,
      )
    }
    let source = experiment.col_data[known_empty_column]
    if source.length() != data.n_barcodes() {
      dua_fail("DropletUtils known-empty column must match the barcode count")
    }
    let parsed : Array[Bool] = []
    for value in source {
      if value == "true" || value == "empty" || value == "1" {
        parsed.push(true)
      } else if value == "false" || value == "cell" || value == "0" {
        parsed.push(false)
      } else {
        dua_fail(
          "DropletUtils known-empty column contains invalid value: " + value,
        )
      }
    }
    parsed
  }
  let result = droplet_utils_empty_drops(data, known_empty~, config~)
  let output = dua_copy_sce(experiment)
  let total_values : Array[String] = []
  let log_values : Array[String] = []
  let p_values : Array[String] = []
  let limited_values : Array[String] = []
  let fdr_values : Array[String] = []
  let class_values : Array[String] = []
  let retained_values : Array[String] = []
  for barcode in 0..<data.n_barcodes() {
    total_values.push(result.totals[barcode].to_string())
    log_values.push(dua_option_double_string(result.log_probabilities[barcode]))
    p_values.push(dua_option_double_string(result.p_values[barcode]))
    limited_values.push(
      match result.limited[barcode] {
        Some(value) => value.to_string()
        None => "NA"
      },
    )
    fdr_values.push(dua_option_double_string(result.fdr[barcode]))
    class_values.push(if result.calls[barcode] { "cell" } else { "empty" })
    retained_values.push(result.always_retained[barcode].to_string())
  }
  output.col_data[output_prefix + ".total"] = total_values
  output.col_data[output_prefix + ".logProb"] = log_values
  output.col_data[output_prefix + ".pValue"] = p_values
  output.col_data[output_prefix + ".limited"] = limited_values
  output.col_data[output_prefix + ".fdr"] = fdr_values
  output.col_data[output_prefix + ".class"] = class_values
  output.col_data[output_prefix + ".retained"] = retained_values
  let ambient_values : Array[String] = []
  let ambient_count_values : Array[String] = []
  for feature in 0..<data.n_features() {
    ambient_values.push(result.ambient.proportions[feature].to_string())
    ambient_count_values.push(result.ambient.raw_counts[feature].to_string())
  }
  output.row_data[output_prefix + ".ambient"] = ambient_values
  output.row_data[output_prefix + ".ambientCount"] = ambient_count_values
  output.metadata[output_prefix + ".lower"] = result.lower.to_string()
  output.metadata[output_prefix + ".retain"] = result.retain.to_string()
  output.metadata[output_prefix + ".iterations"] = result.iterations.to_string()
  output.metadata[output_prefix + ".alpha"] = result.alpha.to_string()
  output.metadata[output_prefix + ".called"] = result.n_called().to_string()
  DropletUtilsSingleCellExperimentOutput::{ experiment: output, result, }
}

///|
pub fn droplet_utils_advanced_example() -> DropletCountMatrix raise DropletUtilsAdvancedError {
  let feature_count = 12
  let barcode_count = 64
  let counts : Array[Array[Int]] = []
  for _ in 0..<feature_count {
    counts.push(Array::make(barcode_count, 0))
  }
  for barcode in 0..<40 {
    let total = 7 + barcode % 8
    for molecule in 0..<total {
      let feature = (molecule * 3 + barcode * 5) % 4
      counts[feature][barcode] = counts[feature][barcode] + 1
    }
  }
  for barcode in 40..<52 {
    counts[0][barcode] = 12 + barcode % 5
    counts[1][barcode] = 9 + barcode % 4
    counts[2][barcode] = 7 + barcode % 3
    counts[3][barcode] = 5 + barcode % 2
    counts[4][barcode] = 34 + barcode % 7
    counts[5][barcode] = 26 + barcode % 6
    counts[6][barcode] = 8 + barcode % 4
    counts[7][barcode] = 3
  }
  for barcode in 52..<64 {
    counts[0][barcode] = 13 + barcode % 5
    counts[1][barcode] = 8 + barcode % 4
    counts[2][barcode] = 6 + barcode % 3
    counts[3][barcode] = 5 + barcode % 2
    counts[8][barcode] = 35 + barcode % 7
    counts[9][barcode] = 25 + barcode % 6
    counts[10][barcode] = 9 + barcode % 4
    counts[11][barcode] = 3
  }
  let features : Array[String] = []
  for feature in 0..<feature_count {
    features.push("gene_" + (feature + 1).to_string())
  }
  let barcodes : Array[String] = []
  for barcode in 0..<barcode_count {
    barcodes.push("barcode_" + (barcode + 1).to_string())
  }
  DropletCountMatrix::create(counts, feature_names=features, barcodes~)
}
