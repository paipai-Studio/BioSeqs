// Compositional differential-abundance inference inspired by Bioconductor
// ALDEx2. Count matrices use feature x sample orientation throughout.

///|
pub suberror Aldex2Error {
  Aldex2Error(String)
}

///|
pub(all) enum Aldex2Denominator {
  Aldex2All
  Aldex2Median
  Aldex2Iqlr
  Aldex2Zero
  Aldex2Lvha
  Aldex2User
} derive(Eq, Debug)

///|
pub struct Aldex2Config {
  mc_samples : Int
  prior : Double
  denominator : Aldex2Denominator
  denominator_indices : Array[Int]
  seed : Int
  paired : Bool
  interval_level : Double
} derive(Eq, Debug)

///|
pub fn Aldex2Config::create(
  mc_samples? : Int = 128,
  prior? : Double = 0.5,
  denominator? : Aldex2Denominator = Aldex2All,
  denominator_indices? : Array[Int] = [],
  seed? : Int = 42,
  paired? : Bool = false,
  interval_level? : Double = 0.95,
) -> Aldex2Config raise Aldex2Error {
  if mc_samples < 1 {
    raise Aldex2Error("ALDEx2 Monte Carlo sample count must be positive")
  }
  if !aldex2_is_finite(prior) || prior <= 0.0 {
    raise Aldex2Error("ALDEx2 Dirichlet prior must be finite and positive")
  }
  if seed <= 0 {
    raise Aldex2Error("ALDEx2 random seed must be positive")
  }
  if !aldex2_is_finite(interval_level) ||
    interval_level <= 0.0 ||
    interval_level >= 1.0 {
    raise Aldex2Error("ALDEx2 interval level must be finite and in (0, 1)")
  }
  if denominator == Aldex2User && denominator_indices.length() == 0 {
    raise Aldex2Error("ALDEx2 user denominator must contain feature indices")
  }
  if denominator != Aldex2User && denominator_indices.length() > 0 {
    raise Aldex2Error(
      "ALDEx2 denominator indices require the user denominator mode",
    )
  }
  Aldex2Config::{
    mc_samples,
    prior,
    denominator,
    denominator_indices: denominator_indices.copy(),
    seed,
    paired,
    interval_level,
  }
}

///|
pub fn Aldex2Config::default() -> Aldex2Config {
  Aldex2Config::{
    mc_samples: 128,
    prior: 0.5,
    denominator: Aldex2All,
    denominator_indices: [],
    seed: 42,
    paired: false,
    interval_level: 0.95,
  }
}

///|
pub struct Aldex2Clr {
  counts : Array[Array[Double]]
  feature_names : Array[String]
  sample_names : Array[String]
  conditions : Array[String]
  group_names : Array[String]
  group_indices : Array[Array[Int]]
  original_feature_indices : Array[Int]
  denominator_by_sample : Array[Array[Int]]
  dirichlet : Array[Array[Array[Double]]]
  analysis : Array[Array[Array[Double]]]
  scale_samples : Array[Array[Double]]
  original_feature_count : Int
  config : Aldex2Config
}

///|
pub struct Aldex2TestFeature {
  feature_name : String
  welch_p_value : Double
  welch_adjusted : Double
  wilcoxon_p_value : Double
  wilcoxon_adjusted : Double
} derive(Eq, Debug)

///|
pub struct Aldex2TestResult {
  features : Array[Aldex2TestFeature]
  group_a : String
  group_b : String
  paired : Bool
  mc_samples : Int
}

///|
pub struct Aldex2EffectFeature {
  feature_name : String
  rab_all : Double
  rab_group_a : Double
  rab_group_b : Double
  diff_between : Double
  diff_within : Double
  effect : Double
  effect_low : Double
  effect_high : Double
  overlap : Double
} derive(Eq, Debug)

///|
pub struct Aldex2EffectResult {
  features : Array[Aldex2EffectFeature]
  group_a : String
  group_b : String
  paired : Bool
  interval_level : Double
}

///|
pub struct Aldex2Result {
  clr : Aldex2Clr
  tests : Aldex2TestResult
  effects : Aldex2EffectResult
}

///|
pub struct Aldex2Summary {
  feature_count : Int
  sample_count : Int
  removed_feature_count : Int
  significant_welch_count : Int
  significant_wilcoxon_count : Int
  large_effect_count : Int
  group_a : String
  group_b : String
  mc_samples : Int
  denominator : String
  scale_model : Bool
} derive(Eq, Debug)

///|
pub struct Aldex2SummarizedExperimentOutput {
  experiment : SummarizedExperiment
  result : Aldex2Result
}

///|
fn aldex2_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn aldex2_log2(value : Double) -> Double {
  @math.ln(value) / @math.ln(2.0)
}

///|
fn aldex2_mean(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let mut total = 0.0
  for value in values {
    total = total + value
  }
  total / values.length().to_double()
}

///|
fn aldex2_variance(values : Array[Double]) -> Double {
  if values.length() < 2 {
    return 0.0
  }
  let center = aldex2_mean(values)
  let mut total = 0.0
  for value in values {
    let delta = value - center
    total = total + delta * delta
  }
  total / (values.length() - 1).to_double()
}

///|
fn aldex2_quantile(values : Array[Double], probability : Double) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let sorted = values.copy()
  sorted.sort_by(fn(left : Double, right : Double) -> Int {
    if left < right {
      -1
    } else if left > right {
      1
    } else {
      0
    }
  })
  if sorted.length() == 1 {
    return sorted[0]
  }
  let position = probability.max(0.0).min(1.0) *
    (sorted.length() - 1).to_double()
  let lower = @math.floor(position).to_int()
  let upper = (lower + 1).min(sorted.length() - 1)
  let fraction = position - lower.to_double()
  sorted[lower] * (1.0 - fraction) + sorted[upper] * fraction
}

///|
fn aldex2_median(values : Array[Double]) -> Double {
  aldex2_quantile(values, 0.5)
}

///|
fn aldex2_default_names(prefix : String, count : Int) -> Array[String] {
  let names : Array[String] = []
  for index in 0..<count {
    names.push(prefix + (index + 1).to_string())
  }
  names
}

///|
fn aldex2_validate_names(
  names : Array[String],
  expected : Int,
  label : String,
) -> Unit raise Aldex2Error {
  if names.length() != expected {
    raise Aldex2Error("ALDEx2 " + label + " name count does not match data")
  }
  for index in 0..<names.length() {
    if names[index].trim().length() == 0 {
      raise Aldex2Error("ALDEx2 " + label + " names must not be empty")
    }
    for prior in 0..<index {
      if names[prior] == names[index] {
        raise Aldex2Error(
          "ALDEx2 duplicate " + label + " name: " + names[index],
        )
      }
    }
  }
}

///|
fn aldex2_groups(
  conditions : Array[String],
) -> (Array[String], Array[Array[Int]]) raise Aldex2Error {
  let names : Array[String] = []
  let indices : Array[Array[Int]] = []
  for sample in 0..<conditions.length() {
    let condition = conditions[sample]
    if condition.trim().length() == 0 {
      raise Aldex2Error("ALDEx2 condition labels must not be empty")
    }
    let mut found = -1
    for index in 0..<names.length() {
      if names[index] == condition {
        found = index
      }
    }
    if found < 0 {
      names.push(condition)
      indices.push([sample])
    } else {
      indices[found].push(sample)
    }
  }
  (names, indices)
}

///|
fn aldex2_all_indices(count : Int) -> Array[Int] {
  let output : Array[Int] = []
  for index in 0..<count {
    output.push(index)
  }
  output
}

///|
fn aldex2_contains(values : Array[Int], target : Int) -> Bool {
  for value in values {
    if value == target {
      return true
    }
  }
  false
}

///|
fn aldex2_intersection(left : Array[Int], right : Array[Int]) -> Array[Int] {
  let output : Array[Int] = []
  for value in left {
    if aldex2_contains(right, value) {
      output.push(value)
    }
  }
  output
}

///|
fn aldex2_lcg(state : Int) -> (Int, Double) {
  let normalized = (state - 1) % 2147483646 + 1
  let quotient = normalized / 44488
  let remainder = normalized % 44488
  let candidate = 48271 * remainder - 3399 * quotient
  let next = if candidate > 0 { candidate } else { candidate + 2147483647 }
  (next, next.to_double() / 2147483647.0)
}

///|
fn aldex2_normal(state : Int) -> (Int, Double) {
  let (state1, raw1) = aldex2_lcg(state)
  let (state2, raw2) = aldex2_lcg(state1)
  let first = raw1.max(1.0e-15)
  let value = (-2.0 * @math.ln(first)).sqrt() *
    @math.cos(6.283185307179586 * raw2)
  (state2, value)
}

///|
fn aldex2_gamma(shape : Double, state : Int) -> (Int, Double) {
  if shape < 1.0 {
    let (state1, uniform) = aldex2_lcg(state)
    let (state2, sampled) = aldex2_gamma(shape + 1.0, state1)
    return (state2, sampled * @math.exp(@math.ln(uniform.max(1.0e-15)) / shape))
  }
  let d = shape - 1.0 / 3.0
  let c = 1.0 / (9.0 * d).sqrt()
  let mut current = state
  while true {
    let (state1, normal) = aldex2_normal(current)
    current = state1
    let base = 1.0 + c * normal
    if base > 0.0 {
      let volume = base * base * base
      let (state2, uniform) = aldex2_lcg(current)
      current = state2
      if uniform < 1.0 - 0.0331 * normal * normal * normal * normal {
        return (current, d * volume)
      }
      if @math.ln(uniform.max(1.0e-15)) <
        0.5 * normal * normal + d * (1.0 - volume + @math.ln(volume)) {
        return (current, d * volume)
      }
    }
  }
  (current, shape)
}

///|
fn aldex2_pseudoclr(
  counts : Array[Array[Double]],
  prior : Double,
) -> Array[Array[Double]] {
  let feature_count = counts.length()
  let sample_count = counts[0].length()
  let output : Array[Array[Double]] = []
  for sample in 0..<sample_count {
    let values : Array[Double] = []
    for feature in 0..<feature_count {
      values.push(aldex2_log2(counts[feature][sample] + prior))
    }
    let center = aldex2_mean(values)
    for feature in 0..<feature_count {
      values[feature] = values[feature] - center
    }
    output.push(values)
  }
  output
}

///|
fn aldex2_variances_for_samples(
  transformed : Array[Array[Double]],
  samples : Array[Int],
) -> Array[Double] {
  let feature_count = transformed[0].length()
  let output : Array[Double] = []
  for feature in 0..<feature_count {
    let values : Array[Double] = []
    for sample in samples {
      values.push(transformed[sample][feature])
    }
    output.push(aldex2_variance(values))
  }
  output
}

///|
fn aldex2_iqlr_indices(
  counts : Array[Array[Double]],
  group_indices : Array[Array[Int]],
  prior : Double,
) -> Array[Int] raise Aldex2Error {
  let transformed = aldex2_pseudoclr(counts, prior)
  let sets : Array[Array[Int]] = []
  for samples in group_indices {
    let variances = aldex2_variances_for_samples(transformed, samples)
    let lower = aldex2_quantile(variances, 0.25)
    let upper = aldex2_quantile(variances, 0.75)
    let selected : Array[Int] = []
    for feature in 0..<variances.length() {
      if variances[feature] > lower && variances[feature] < upper {
        selected.push(feature)
      }
    }
    sets.push(selected)
  }
  let all_samples = aldex2_all_indices(transformed.length())
  let all_variances = aldex2_variances_for_samples(transformed, all_samples)
  let all_lower = aldex2_quantile(all_variances, 0.25)
  let all_upper = aldex2_quantile(all_variances, 0.75)
  let all_selected : Array[Int] = []
  for feature in 0..<all_variances.length() {
    if all_variances[feature] > all_lower && all_variances[feature] < all_upper {
      all_selected.push(feature)
    }
  }
  sets.push(all_selected)
  let mut intersection = sets[0].copy()
  for index in 1..<sets.length() {
    intersection = aldex2_intersection(intersection, sets[index])
  }
  if intersection.length() < 2 {
    raise Aldex2Error(
      "ALDEx2 IQLR denominator has fewer than two invariant features",
    )
  }
  intersection
}

///|
fn aldex2_lvha_indices(
  counts : Array[Array[Double]],
  group_indices : Array[Array[Int]],
  prior : Double,
) -> Array[Int] raise Aldex2Error {
  let transformed = aldex2_pseudoclr(counts, prior)
  let sets : Array[Array[Int]] = []
  for samples in group_indices {
    let variances = aldex2_variances_for_samples(transformed, samples)
    let abundances : Array[Double] = []
    for feature in 0..<counts.length() {
      let values : Array[Double] = []
      for sample in samples {
        values.push(transformed[sample][feature])
      }
      abundances.push(aldex2_mean(values))
    }
    let variance_cutoff = aldex2_quantile(variances, 0.25)
    let abundance_cutoff = aldex2_quantile(abundances, 0.75)
    let selected : Array[Int] = []
    for feature in 0..<counts.length() {
      if variances[feature] <= variance_cutoff &&
        abundances[feature] >= abundance_cutoff {
        selected.push(feature)
      }
    }
    sets.push(selected)
  }
  let mut intersection = sets[0].copy()
  for index in 1..<sets.length() {
    intersection = aldex2_intersection(intersection, sets[index])
  }
  if intersection.length() < 2 {
    raise Aldex2Error(
      "ALDEx2 LVHA denominator has fewer than two shared features",
    )
  }
  intersection
}

///|
fn aldex2_denominators(
  counts : Array[Array[Double]],
  groups : Array[Array[Int]],
  conditions : Array[String],
  config : Aldex2Config,
) -> Array[Array[Int]] raise Aldex2Error {
  let feature_count = counts.length()
  let shared = match config.denominator {
    Aldex2Iqlr => aldex2_iqlr_indices(counts, groups, config.prior)
    Aldex2Lvha => aldex2_lvha_indices(counts, groups, config.prior)
    Aldex2User => {
      let selected : Array[Int] = []
      for index in config.denominator_indices {
        if index < 0 || index >= feature_count {
          raise Aldex2Error(
            "ALDEx2 user denominator feature index is out of bounds",
          )
        }
        if aldex2_contains(selected, index) {
          raise Aldex2Error(
            "ALDEx2 user denominator contains duplicate feature indices",
          )
        }
        selected.push(index)
      }
      selected
    }
    _ => aldex2_all_indices(feature_count)
  }
  let by_group : Array[Array[Int]] = []
  if config.denominator == Aldex2Zero {
    for samples in groups {
      let selected : Array[Int] = []
      for feature in 0..<feature_count {
        let mut total = 0.0
        for sample in samples {
          total = total + counts[feature][sample]
        }
        if total > 0.0 {
          selected.push(feature)
        }
      }
      if selected.length() == 0 {
        raise Aldex2Error(
          "ALDEx2 zero denominator has no observed features in a condition",
        )
      }
      by_group.push(selected)
    }
  } else {
    for _ in groups {
      by_group.push(shared.copy())
    }
  }
  let by_sample : Array[Array[Int]] = []
  for sample in 0..<conditions.length() {
    let mut group = -1
    for group_index in 0..<groups.length() {
      if aldex2_contains(groups[group_index], sample) {
        group = group_index
      }
    }
    if group < 0 {
      raise Aldex2Error("ALDEx2 internal condition index is invalid")
    }
    by_sample.push(by_group[group].copy())
  }
  by_sample
}

///|
fn aldex2_validate_scale_samples(
  scale_samples : Array[Array[Double]],
  sample_count : Int,
  mc_samples : Int,
) -> Unit raise Aldex2Error {
  if scale_samples.length() == 0 {
    return
  }
  if scale_samples.length() != sample_count {
    raise Aldex2Error("ALDEx2 scale sample rows must match samples")
  }
  for sample in 0..<scale_samples.length() {
    if scale_samples[sample].length() != mc_samples {
      raise Aldex2Error(
        "ALDEx2 scale sample columns must match Monte Carlo samples",
      )
    }
    for value in scale_samples[sample] {
      if !aldex2_is_finite(value) || value <= 0.0 {
        raise Aldex2Error("ALDEx2 scale samples must be finite and positive")
      }
    }
  }
}

///|
pub fn aldex2_clr(
  counts : Array[Array[Double]],
  conditions : Array[String],
  config? : Aldex2Config = Aldex2Config::default(),
  feature_names? : Array[String] = [],
  sample_names? : Array[String] = [],
  scale_samples? : Array[Array[Double]] = [],
) -> Aldex2Clr raise Aldex2Error {
  if counts.length() == 0 {
    raise Aldex2Error("ALDEx2 count matrix must contain features")
  }
  let sample_count = counts[0].length()
  if sample_count == 0 {
    raise Aldex2Error("ALDEx2 count matrix must contain samples")
  }
  if conditions.length() != sample_count {
    raise Aldex2Error("ALDEx2 condition count must match samples")
  }
  let supplied_features = if feature_names.length() == 0 {
    aldex2_default_names("feature_", counts.length())
  } else {
    feature_names.copy()
  }
  let supplied_samples = if sample_names.length() == 0 {
    aldex2_default_names("sample_", sample_count)
  } else {
    sample_names.copy()
  }
  aldex2_validate_names(supplied_features, counts.length(), "feature")
  aldex2_validate_names(supplied_samples, sample_count, "sample")
  let library_sizes = Array::make(sample_count, 0.0)
  let filtered_counts : Array[Array[Double]] = []
  let filtered_names : Array[String] = []
  let original_indices : Array[Int] = []
  for feature in 0..<counts.length() {
    if counts[feature].length() != sample_count {
      raise Aldex2Error("ALDEx2 count matrix must not be ragged")
    }
    let mut row_total = 0.0
    let copied : Array[Double] = []
    for sample in 0..<sample_count {
      let value = counts[feature][sample]
      if !aldex2_is_finite(value) {
        raise Aldex2Error("ALDEx2 counts must be finite")
      }
      if value < 0.0 {
        raise Aldex2Error("ALDEx2 counts must be non-negative")
      }
      if (value - value.round()).abs() > 1.0e-9 {
        raise Aldex2Error("ALDEx2 counts must be integers")
      }
      copied.push(value)
      row_total = row_total + value
      library_sizes[sample] = library_sizes[sample] + value
    }
    if row_total > 0.0 {
      filtered_counts.push(copied)
      filtered_names.push(supplied_features[feature])
      original_indices.push(feature)
    }
  }
  if filtered_counts.length() == 0 {
    raise Aldex2Error("ALDEx2 count matrix contains only zero features")
  }
  for sample in 0..<sample_count {
    if library_sizes[sample] <= 0.0 {
      raise Aldex2Error(
        "ALDEx2 sample has an empty library: " + supplied_samples[sample],
      )
    }
  }
  let (group_names, group_indices) = aldex2_groups(conditions)
  if group_names.length() < 2 {
    raise Aldex2Error("ALDEx2 requires at least two conditions")
  }
  if config.paired {
    if group_names.length() != 2 {
      raise Aldex2Error("ALDEx2 paired analysis requires two conditions")
    }
    if group_indices[0].length() != group_indices[1].length() {
      raise Aldex2Error(
        "ALDEx2 paired conditions must have equal sample counts",
      )
    }
  }
  aldex2_validate_scale_samples(scale_samples, sample_count, config.mc_samples)
  let denominator_config = if config.denominator == Aldex2User {
    let mapped_indices : Array[Int] = []
    for requested in config.denominator_indices {
      if requested < 0 || requested >= counts.length() {
        raise Aldex2Error(
          "ALDEx2 user denominator feature index is out of bounds",
        )
      }
      let mut filtered_index = -1
      for index in 0..<original_indices.length() {
        if original_indices[index] == requested {
          filtered_index = index
        }
      }
      if filtered_index < 0 {
        raise Aldex2Error(
          "ALDEx2 user denominator references an all-zero feature",
        )
      }
      mapped_indices.push(filtered_index)
    }
    Aldex2Config::{
      mc_samples: config.mc_samples,
      prior: config.prior,
      denominator: config.denominator,
      denominator_indices: mapped_indices,
      seed: config.seed,
      paired: config.paired,
      interval_level: config.interval_level,
    }
  } else {
    config
  }
  let denominators = aldex2_denominators(
    filtered_counts, group_indices, conditions, denominator_config,
  )
  let feature_count = filtered_counts.length()
  let dirichlet : Array[Array[Array[Double]]] = []
  let analysis : Array[Array[Array[Double]]] = []
  let mut state = config.seed
  for sample in 0..<sample_count {
    let sample_dirichlet : Array[Array[Double]] = []
    let sample_analysis : Array[Array[Double]] = []
    for _ in 0..<feature_count {
      sample_dirichlet.push(Array::make(config.mc_samples, 0.0))
      sample_analysis.push(Array::make(config.mc_samples, 0.0))
    }
    for instance in 0..<config.mc_samples {
      let gamma_values = Array::make(feature_count, 0.0)
      let mut gamma_total = 0.0
      for feature in 0..<feature_count {
        let shape = filtered_counts[feature][sample] + config.prior
        let (next_state, sampled) = aldex2_gamma(shape, state)
        state = next_state
        gamma_values[feature] = sampled.max(1.0e-300)
        gamma_total = gamma_total + gamma_values[feature]
      }
      let log_values = Array::make(feature_count, 0.0)
      for feature in 0..<feature_count {
        let proportion = gamma_values[feature] / gamma_total
        sample_dirichlet[feature][instance] = proportion
        log_values[feature] = aldex2_log2(proportion)
      }
      if scale_samples.length() > 0 {
        let scale = aldex2_log2(scale_samples[sample][instance])
        for feature in 0..<feature_count {
          sample_analysis[feature][instance] = log_values[feature] + scale
        }
      } else {
        let denominator_values : Array[Double] = []
        for feature in denominators[sample] {
          denominator_values.push(log_values[feature])
        }
        let center = if config.denominator == Aldex2Median {
          aldex2_median(denominator_values)
        } else {
          aldex2_mean(denominator_values)
        }
        for feature in 0..<feature_count {
          sample_analysis[feature][instance] = log_values[feature] - center
        }
      }
    }
    dirichlet.push(sample_dirichlet)
    analysis.push(sample_analysis)
  }
  Aldex2Clr::{
    counts: filtered_counts,
    feature_names: filtered_names,
    sample_names: supplied_samples,
    conditions: conditions.copy(),
    group_names,
    group_indices,
    original_feature_indices: original_indices,
    denominator_by_sample: denominators,
    dirichlet,
    analysis,
    scale_samples: {
      let copied : Array[Array[Double]] = []
      for row in scale_samples {
        copied.push(row.copy())
      }
      copied
    },
    original_feature_count: counts.length(),
    config,
  }
}

///|
pub fn Aldex2Clr::feature_count(self : Aldex2Clr) -> Int {
  self.feature_names.length()
}

///|
pub fn Aldex2Clr::sample_count(self : Aldex2Clr) -> Int {
  self.sample_names.length()
}

///|
pub fn Aldex2Clr::mc_sample_count(self : Aldex2Clr) -> Int {
  self.config.mc_samples
}

///|
pub fn Aldex2Clr::removed_feature_count(self : Aldex2Clr) -> Int {
  self.original_feature_count - self.feature_count()
}

///|
pub fn Aldex2Clr::sample_instance(
  self : Aldex2Clr,
  sample : Int,
  instance : Int,
) -> Array[Double]? {
  if sample < 0 ||
    sample >= self.sample_count() ||
    instance < 0 ||
    instance >= self.config.mc_samples {
    return None
  }
  let output : Array[Double] = []
  for feature in 0..<self.feature_count() {
    output.push(self.analysis[sample][feature][instance])
  }
  Some(output)
}

///|
pub fn Aldex2Clr::feature_instance(
  self : Aldex2Clr,
  feature : Int,
  instance : Int,
) -> Array[Double]? {
  if feature < 0 ||
    feature >= self.feature_count() ||
    instance < 0 ||
    instance >= self.config.mc_samples {
    return None
  }
  let output : Array[Double] = []
  for sample in 0..<self.sample_count() {
    output.push(self.analysis[sample][feature][instance])
  }
  Some(output)
}

///|
pub fn Aldex2Clr::expected_clr(
  self : Aldex2Clr,
  feature : Int,
  sample : Int,
) -> Double? {
  if feature < 0 ||
    feature >= self.feature_count() ||
    sample < 0 ||
    sample >= self.sample_count() {
    return None
  }
  Some(aldex2_median(self.analysis[sample][feature]))
}

///|
pub fn Aldex2Clr::expected_proportion(
  self : Aldex2Clr,
  feature : Int,
  sample : Int,
) -> Double? {
  if feature < 0 ||
    feature >= self.feature_count() ||
    sample < 0 ||
    sample >= self.sample_count() {
    return None
  }
  Some(aldex2_median(self.dirichlet[sample][feature]))
}

///|
fn aldex2_welch(
  group_a : Array[Double],
  group_b : Array[Double],
) -> (Double, Double) {
  if group_a.length() < 2 || group_b.length() < 2 {
    return (0.0, 1.0)
  }
  let mean_a = aldex2_mean(group_a)
  let mean_b = aldex2_mean(group_b)
  let variance_a = aldex2_variance(group_a)
  let variance_b = aldex2_variance(group_b)
  let component_a = variance_a / group_a.length().to_double()
  let component_b = variance_b / group_b.length().to_double()
  let standard_error = (component_a + component_b).sqrt()
  if standard_error <= 1.0e-15 {
    if (mean_a - mean_b).abs() <= 1.0e-15 {
      return (0.0, 1.0)
    }
    return (if mean_b > mean_a { 1.0e12 } else { -1.0e12 }, 0.0)
  }
  let statistic = (mean_b - mean_a) / standard_error
  let denominator = component_a *
    component_a /
    (group_a.length() - 1).to_double() +
    component_b * component_b / (group_b.length() - 1).to_double()
  let degrees = if denominator > 0.0 {
    (component_a + component_b) * (component_a + component_b) / denominator
  } else {
    (group_a.length() + group_b.length() - 2).to_double()
  }
  let p_value = (2.0 * (1.0 - stat_t_cdf(statistic.abs(), degrees)))
    .max(0.0)
    .min(1.0)
  (statistic, p_value)
}

///|
fn aldex2_paired_t(
  group_a : Array[Double],
  group_b : Array[Double],
) -> (Double, Double) {
  if group_a.length() != group_b.length() || group_a.length() < 2 {
    return (0.0, 1.0)
  }
  let differences : Array[Double] = []
  for index in 0..<group_a.length() {
    differences.push(group_b[index] - group_a[index])
  }
  let center = aldex2_mean(differences)
  let variance = aldex2_variance(differences)
  let standard_error = (variance / differences.length().to_double()).sqrt()
  if standard_error <= 1.0e-15 {
    if center.abs() <= 1.0e-15 {
      return (0.0, 1.0)
    }
    return (if center > 0.0 { 1.0e12 } else { -1.0e12 }, 0.0)
  }
  let statistic = center / standard_error
  let p_value = (2.0 *
    (1.0 - stat_t_cdf(statistic.abs(), (differences.length() - 1).to_double())))
    .max(0.0)
    .min(1.0)
  (statistic, p_value)
}

///|
fn aldex2_bh_adjust(p_values : Array[Double]) -> Array[Double] {
  let count = p_values.length()
  if count == 0 {
    return []
  }
  let sorted : Array[(Double, Int)] = []
  for index in 0..<count {
    sorted.push((p_values[index].max(0.0).min(1.0), index))
  }
  sorted.sort_by(fn(left : (Double, Int), right : (Double, Int)) -> Int {
    if left.0 < right.0 {
      -1
    } else if left.0 > right.0 {
      1
    } else {
      left.1 - right.1
    }
  })
  let adjusted = Array::make(count, 1.0)
  let mut running_minimum = 1.0
  let mut rank = count - 1
  while rank >= 0 {
    let candidate = (sorted[rank].0 * count.to_double() / (rank + 1).to_double()).min(
      1.0,
    )
    running_minimum = running_minimum.min(candidate)
    adjusted[sorted[rank].1] = running_minimum
    rank = rank - 1
  }
  adjusted
}

///|
pub fn aldex2_ttest(clr : Aldex2Clr) -> Aldex2TestResult raise Aldex2Error {
  if clr.group_names.length() != 2 {
    raise Aldex2Error("ALDEx2 t-test requires exactly two conditions")
  }
  if clr.group_indices[0].length() < 2 || clr.group_indices[1].length() < 2 {
    raise Aldex2Error(
      "ALDEx2 t-test requires at least two replicates per condition",
    )
  }
  let feature_count = clr.feature_count()
  let mc_samples = clr.config.mc_samples
  let welch : Array[Array[Double]] = []
  let wilcoxon : Array[Array[Double]] = []
  let welch_adjusted : Array[Array[Double]] = []
  let wilcoxon_adjusted : Array[Array[Double]] = []
  for _ in 0..<feature_count {
    welch.push(Array::make(mc_samples, 1.0))
    wilcoxon.push(Array::make(mc_samples, 1.0))
    welch_adjusted.push(Array::make(mc_samples, 1.0))
    wilcoxon_adjusted.push(Array::make(mc_samples, 1.0))
  }
  for instance in 0..<mc_samples {
    let welch_instance = Array::make(feature_count, 1.0)
    let wilcoxon_instance = Array::make(feature_count, 1.0)
    for feature in 0..<feature_count {
      let group_a : Array[Double] = []
      let group_b : Array[Double] = []
      for sample in clr.group_indices[0] {
        group_a.push(clr.analysis[sample][feature][instance])
      }
      for sample in clr.group_indices[1] {
        group_b.push(clr.analysis[sample][feature][instance])
      }
      let (_, welch_p) = if clr.config.paired {
        aldex2_paired_t(group_a, group_b)
      } else {
        aldex2_welch(group_a, group_b)
      }
      let (_, wilcoxon_p) = if clr.config.paired {
        stat_wilcoxon_signed_rank(group_a, group_b)
      } else {
        stat_mann_whitney(group_a, group_b)
      }
      welch_instance[feature] = welch_p
      wilcoxon_instance[feature] = wilcoxon_p
    }
    let adjusted_welch = aldex2_bh_adjust(welch_instance)
    let adjusted_wilcoxon = aldex2_bh_adjust(wilcoxon_instance)
    for feature in 0..<feature_count {
      welch[feature][instance] = welch_instance[feature]
      wilcoxon[feature][instance] = wilcoxon_instance[feature]
      welch_adjusted[feature][instance] = adjusted_welch[feature]
      wilcoxon_adjusted[feature][instance] = adjusted_wilcoxon[feature]
    }
  }
  let features : Array[Aldex2TestFeature] = []
  for feature in 0..<feature_count {
    features.push(Aldex2TestFeature::{
      feature_name: clr.feature_names[feature],
      welch_p_value: aldex2_mean(welch[feature]),
      welch_adjusted: aldex2_mean(welch_adjusted[feature]),
      wilcoxon_p_value: aldex2_mean(wilcoxon[feature]),
      wilcoxon_adjusted: aldex2_mean(wilcoxon_adjusted[feature]),
    })
  }
  Aldex2TestResult::{
    features,
    group_a: clr.group_names[0],
    group_b: clr.group_names[1],
    paired: clr.config.paired,
    mc_samples,
  }
}

///|
fn aldex2_pooled_feature(
  clr : Aldex2Clr,
  feature : Int,
  samples : Array[Int],
) -> Array[Double] {
  let output : Array[Double] = []
  for sample in samples {
    for instance in 0..<clr.config.mc_samples {
      output.push(clr.analysis[sample][feature][instance])
    }
  }
  output
}

///|
fn aldex2_unpaired_effect_distribution(
  clr : Aldex2Clr,
  feature : Int,
) -> (Array[Double], Array[Double]) {
  let between : Array[Double] = []
  let within : Array[Double] = []
  for instance in 0..<clr.config.mc_samples {
    let group_a : Array[Double] = []
    let group_b : Array[Double] = []
    for sample in clr.group_indices[0] {
      group_a.push(clr.analysis[sample][feature][instance])
    }
    for sample in clr.group_indices[1] {
      group_b.push(clr.analysis[sample][feature][instance])
    }
    let within_a : Array[Double] = []
    let within_b : Array[Double] = []
    for left in 0..<group_a.length() {
      for right in (left + 1)..<group_a.length() {
        within_a.push((group_a[left] - group_a[right]).abs())
      }
    }
    for left in 0..<group_b.length() {
      for right in (left + 1)..<group_b.length() {
        within_b.push((group_b[left] - group_b[right]).abs())
      }
    }
    let default_a = if within_a.length() == 0 {
      0.0
    } else {
      aldex2_median(within_a)
    }
    let default_b = if within_b.length() == 0 {
      0.0
    } else {
      aldex2_median(within_b)
    }
    let comparisons = group_a.length() * group_b.length()
    let mut comparison = 0
    for value_a in group_a {
      for value_b in group_b {
        between.push(value_b - value_a)
        let spread_a = if within_a.length() == 0 {
          default_a
        } else {
          within_a[comparison % within_a.length()]
        }
        let spread_b = if within_b.length() == 0 {
          default_b
        } else {
          within_b[comparison % within_b.length()]
        }
        within.push(spread_a.max(spread_b))
        comparison = comparison + 1
      }
    }
    ignore(comparisons)
  }
  (between, within)
}

///|
fn aldex2_paired_effect_distribution(
  clr : Aldex2Clr,
  feature : Int,
) -> (Array[Double], Array[Double]) {
  let between : Array[Double] = []
  for instance in 0..<clr.config.mc_samples {
    for pair in 0..<clr.group_indices[0].length() {
      let sample_a = clr.group_indices[0][pair]
      let sample_b = clr.group_indices[1][pair]
      between.push(
        clr.analysis[sample_b][feature][instance] -
        clr.analysis[sample_a][feature][instance],
      )
    }
  }
  let center = aldex2_median(between)
  let deviations : Array[Double] = []
  for value in between {
    deviations.push((value - center).abs())
  }
  let mad = 1.428 * aldex2_median(deviations)
  (between, Array::make(between.length(), mad))
}

///|
pub fn aldex2_effect(clr : Aldex2Clr) -> Aldex2EffectResult raise Aldex2Error {
  if clr.group_names.length() != 2 {
    raise Aldex2Error(
      "ALDEx2 effect estimation requires exactly two conditions",
    )
  }
  if clr.group_indices[0].length() < 2 || clr.group_indices[1].length() < 2 {
    raise Aldex2Error(
      "ALDEx2 effect estimation requires at least two replicates per condition",
    )
  }
  let all_samples = aldex2_all_indices(clr.sample_count())
  let alpha = (1.0 - clr.config.interval_level) / 2.0
  let features : Array[Aldex2EffectFeature] = []
  for feature in 0..<clr.feature_count() {
    let pooled_all = aldex2_pooled_feature(clr, feature, all_samples)
    let pooled_a = aldex2_pooled_feature(clr, feature, clr.group_indices[0])
    let pooled_b = aldex2_pooled_feature(clr, feature, clr.group_indices[1])
    let (between, within) = if clr.config.paired {
      aldex2_paired_effect_distribution(clr, feature)
    } else {
      aldex2_unpaired_effect_distribution(clr, feature)
    }
    let effect_distribution : Array[Double] = []
    let mut negative = 0
    let mut positive = 0
    for index in 0..<between.length() {
      let spread = within[index]
      let effect = if spread > 1.0e-12 {
        between[index] / spread
      } else if between[index].abs() <= 1.0e-12 {
        0.0
      } else {
        between[index] / 1.0e-12
      }
      effect_distribution.push(effect)
      if effect < 0.0 {
        negative = negative + 1
      } else if effect > 0.0 {
        positive = positive + 1
      }
    }
    let overlap_denominator = (negative + positive).to_double() + 1.0
    let overlap = if overlap_denominator > 0.0 {
      (negative.min(positive).to_double() + 0.5) / overlap_denominator
    } else {
      0.5
    }
    features.push(Aldex2EffectFeature::{
      feature_name: clr.feature_names[feature],
      rab_all: aldex2_median(pooled_all),
      rab_group_a: aldex2_median(pooled_a),
      rab_group_b: aldex2_median(pooled_b),
      diff_between: aldex2_median(between),
      diff_within: aldex2_median(within),
      effect: aldex2_median(effect_distribution),
      effect_low: aldex2_quantile(effect_distribution, alpha),
      effect_high: aldex2_quantile(effect_distribution, 1.0 - alpha),
      overlap,
    })
  }
  Aldex2EffectResult::{
    features,
    group_a: clr.group_names[0],
    group_b: clr.group_names[1],
    paired: clr.config.paired,
    interval_level: clr.config.interval_level,
  }
}

///|
pub fn aldex2(
  counts : Array[Array[Double]],
  conditions : Array[String],
  config? : Aldex2Config = Aldex2Config::default(),
  feature_names? : Array[String] = [],
  sample_names? : Array[String] = [],
  scale_samples? : Array[Array[Double]] = [],
) -> Aldex2Result raise Aldex2Error {
  let clr = aldex2_clr(
    counts,
    conditions,
    config~,
    feature_names~,
    sample_names~,
    scale_samples~,
  )
  let tests = aldex2_ttest(clr)
  let effects = aldex2_effect(clr)
  Aldex2Result::{ clr, tests, effects }
}

///|
pub fn Aldex2TestResult::find_feature(
  self : Aldex2TestResult,
  name : String,
) -> Aldex2TestFeature? {
  for feature in self.features {
    if feature.feature_name == name {
      return Some(feature)
    }
  }
  None
}

///|
pub fn Aldex2EffectResult::find_feature(
  self : Aldex2EffectResult,
  name : String,
) -> Aldex2EffectFeature? {
  for feature in self.features {
    if feature.feature_name == name {
      return Some(feature)
    }
  }
  None
}

///|
pub fn Aldex2Result::select(
  self : Aldex2Result,
  maximum_adjusted? : Double = 0.05,
  minimum_effect? : Double = 1.0,
) -> Array[Aldex2EffectFeature] {
  let selected : Array[Aldex2EffectFeature] = []
  for index in 0..<self.effects.features.length() {
    if self.tests.features[index].welch_adjusted <= maximum_adjusted &&
      self.effects.features[index].effect.abs() >= minimum_effect {
      selected.push(self.effects.features[index])
    }
  }
  selected
}

///|
pub fn Aldex2Result::ranked(self : Aldex2Result) -> Array[Aldex2EffectFeature] {
  let ranked = self.effects.features.copy()
  ranked.sort_by(fn(
    left : Aldex2EffectFeature,
    right : Aldex2EffectFeature,
  ) -> Int {
    if left.effect.abs() > right.effect.abs() {
      -1
    } else if left.effect.abs() < right.effect.abs() {
      1
    } else if left.feature_name < right.feature_name {
      -1
    } else if left.feature_name > right.feature_name {
      1
    } else {
      0
    }
  })
  ranked
}

///|
fn aldex2_denominator_name(value : Aldex2Denominator) -> String {
  match value {
    Aldex2All => "all"
    Aldex2Median => "median"
    Aldex2Iqlr => "iqlr"
    Aldex2Zero => "zero"
    Aldex2Lvha => "lvha"
    Aldex2User => "user"
  }
}

///|
pub fn Aldex2Result::summary(
  self : Aldex2Result,
  maximum_adjusted? : Double = 0.05,
  minimum_effect? : Double = 1.0,
) -> Aldex2Summary {
  let mut welch_count = 0
  let mut wilcoxon_count = 0
  let mut effect_count = 0
  for index in 0..<self.effects.features.length() {
    if self.tests.features[index].welch_adjusted <= maximum_adjusted {
      welch_count = welch_count + 1
    }
    if self.tests.features[index].wilcoxon_adjusted <= maximum_adjusted {
      wilcoxon_count = wilcoxon_count + 1
    }
    if self.effects.features[index].effect.abs() >= minimum_effect {
      effect_count = effect_count + 1
    }
  }
  Aldex2Summary::{
    feature_count: self.clr.feature_count(),
    sample_count: self.clr.sample_count(),
    removed_feature_count: self.clr.removed_feature_count(),
    significant_welch_count: welch_count,
    significant_wilcoxon_count: wilcoxon_count,
    large_effect_count: effect_count,
    group_a: self.tests.group_a,
    group_b: self.tests.group_b,
    mc_samples: self.clr.config.mc_samples,
    denominator: aldex2_denominator_name(self.clr.config.denominator),
    scale_model: self.clr.scale_samples.length() > 0,
  }
}

///|
pub fn Aldex2Clr::expected_distance(self : Aldex2Clr) -> Array[Array[Double]] {
  let sample_count = self.sample_count()
  let distances : Array[Array[Double]] = []
  for _ in 0..<sample_count {
    distances.push(Array::make(sample_count, 0.0))
  }
  for left in 0..<sample_count {
    for right in (left + 1)..<sample_count {
      let posterior : Array[Double] = []
      for instance in 0..<self.config.mc_samples {
        let mut squared = 0.0
        for feature in 0..<self.feature_count() {
          let delta = self.analysis[left][feature][instance] -
            self.analysis[right][feature][instance]
          squared = squared + delta * delta
        }
        posterior.push(squared.sqrt())
      }
      let expected = aldex2_median(posterior)
      distances[left][right] = expected
      distances[right][left] = expected
    }
  }
  distances
}

///|
pub fn Aldex2Result::to_tsv(self : Aldex2Result) -> String {
  let buffer = StringBuilder::new()
  buffer.write_string(
    "feature\trab.all\trab.group.a\trab.group.b\tdiff.btw\tdiff.win\teffect\teffect.low\teffect.high\toverlap\twe.ep\twe.eBH\twi.ep\twi.eBH\n",
  )
  for index in 0..<self.effects.features.length() {
    let effect = self.effects.features[index]
    let test_result = self.tests.features[index]
    buffer.write_string(
      effect.feature_name +
      "\t" +
      effect.rab_all.to_string() +
      "\t" +
      effect.rab_group_a.to_string() +
      "\t" +
      effect.rab_group_b.to_string() +
      "\t" +
      effect.diff_between.to_string() +
      "\t" +
      effect.diff_within.to_string() +
      "\t" +
      effect.effect.to_string() +
      "\t" +
      effect.effect_low.to_string() +
      "\t" +
      effect.effect_high.to_string() +
      "\t" +
      effect.overlap.to_string() +
      "\t" +
      test_result.welch_p_value.to_string() +
      "\t" +
      test_result.welch_adjusted.to_string() +
      "\t" +
      test_result.wilcoxon_p_value.to_string() +
      "\t" +
      test_result.wilcoxon_adjusted.to_string() +
      "\n",
    )
  }
  buffer.to_string()
}

///|
fn aldex2_copy_experiment(
  experiment : SummarizedExperiment,
) -> SummarizedExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  for key in experiment.assays.keys() {
    let matrix : Array[Array[Double]] = []
    for row in experiment.assays[key] {
      matrix.push(row.copy())
    }
    assays[key] = matrix
  }
  let metadata : Map[String, String] = Map([])
  for key in experiment.metadata.keys() {
    metadata[key] = experiment.metadata[key]
  }
  SummarizedExperiment::{
    assays,
    row_ranges: experiment.row_ranges.copy(),
    col_data: experiment.col_data.copy(),
    metadata,
  }
}

///|
pub fn aldex2_summarized_experiment(
  experiment : SummarizedExperiment,
  conditions : Array[String],
  config? : Aldex2Config = Aldex2Config::default(),
  assay_name? : String = "counts",
  feature_names? : Array[String] = [],
  sample_names? : Array[String] = [],
  scale_samples? : Array[Array[Double]] = [],
) -> Aldex2SummarizedExperimentOutput raise Aldex2Error {
  let counts = match se_assay(experiment, assay_name) {
    Some(value) => value
    None =>
      raise Aldex2Error(
        "ALDEx2 SummarizedExperiment assay not found: " + assay_name,
      )
  }
  let result = aldex2(
    counts,
    conditions,
    config~,
    feature_names~,
    sample_names~,
    scale_samples~,
  )
  let enriched = aldex2_copy_experiment(experiment)
  let effect = Array::make(result.clr.original_feature_count, [0.0])
  let overlap = Array::make(result.clr.original_feature_count, [1.0])
  let welch = Array::make(result.clr.original_feature_count, [1.0])
  let wilcoxon = Array::make(result.clr.original_feature_count, [1.0])
  for index in 0..<result.clr.feature_count() {
    let original = result.clr.original_feature_indices[index]
    effect[original] = [result.effects.features[index].effect]
    overlap[original] = [result.effects.features[index].overlap]
    welch[original] = [result.tests.features[index].welch_adjusted]
    wilcoxon[original] = [result.tests.features[index].wilcoxon_adjusted]
  }
  enriched.assays["aldex2_effect"] = effect
  enriched.assays["aldex2_overlap"] = overlap
  enriched.assays["aldex2_welch_bh"] = welch
  enriched.assays["aldex2_wilcoxon_bh"] = wilcoxon
  enriched.metadata["aldex2_mc_samples"] = config.mc_samples.to_string()
  enriched.metadata["aldex2_denominator"] = aldex2_denominator_name(
    config.denominator,
  )
  enriched.metadata["aldex2_groups"] = result.tests.group_a +
    "," +
    result.tests.group_b
  Aldex2SummarizedExperimentOutput::{ experiment: enriched, result }
}

///|
pub fn aldex2_example_data() -> (
  Array[Array[Double]],
  Array[String],
  Array[String],
  Array[String],
) {
  let counts = [
    [120.0, 132.0, 115.0, 126.0, 480.0, 510.0, 495.0, 525.0],
    [360.0, 340.0, 375.0, 355.0, 88.0, 95.0, 82.0, 90.0],
    [210.0, 220.0, 205.0, 215.0, 225.0, 218.0, 230.0, 222.0],
    [75.0, 82.0, 78.0, 80.0, 92.0, 86.0, 89.0, 94.0],
    [42.0, 38.0, 45.0, 40.0, 44.0, 41.0, 43.0, 39.0],
    [18.0, 24.0, 20.0, 22.0, 19.0, 21.0, 23.0, 20.0],
    [6.0, 8.0, 5.0, 7.0, 26.0, 30.0, 28.0, 32.0],
    [14.0, 12.0, 16.0, 13.0, 15.0, 14.0, 13.0, 16.0],
    [90.0, 96.0, 88.0, 93.0, 97.0, 91.0, 95.0, 92.0],
    [33.0, 36.0, 31.0, 35.0, 34.0, 32.0, 37.0, 33.0],
    [55.0, 51.0, 58.0, 54.0, 57.0, 53.0, 56.0, 52.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
  ]
  let conditions = [
    "control", "control", "control", "control", "treated", "treated", "treated",
    "treated",
  ]
  let features = [
    "increased", "decreased", "stable_high", "mild", "stable_1", "stable_2", "rare_increased",
    "stable_3", "stable_4", "stable_5", "stable_6", "all_zero",
  ]
  let samples = ["C1", "C2", "C3", "C4", "T1", "T2", "T3", "T4"]
  (counts, conditions, features, samples)
}
