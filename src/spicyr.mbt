// Differential cell-type colocalization inspired by Bioconductor spicyR.
//
// The spatial layer computes cross-L(r) curves for ordered cell-type pairs
// within each image. The modelling layer fits weighted linear models or
// subject-level random-intercept models through variancePartition/dream.

///|
pub suberror SpicyError {
  SpicyError(String)
}

///|
pub enum SpicyModelKind {
  WeightedLinear
  WeightedRandomIntercept
} derive(Eq, Debug)
pub extend SpicyModelKind with Eq::{equal, not_equal}
pub extend SpicyModelKind with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct SpicyCell {
  image_id : String
  cell_type : String
  x : Double
  y : Double
} derive(Eq, Debug)
pub extend SpicyCell with Eq::{equal, not_equal}
pub extend SpicyCell with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn SpicyCell::create(
  image_id : String,
  cell_type : String,
  x : Double,
  y : Double,
) -> SpicyCell raise SpicyError {
  if image_id == "" {
    raise SpicyError("spicyR cell image ID must not be empty")
  }
  if cell_type == "" {
    raise SpicyError("spicyR cell type must not be empty")
  }
  if !spicy_is_finite(x) || !spicy_is_finite(y) {
    raise SpicyError("spicyR cell coordinates must be finite")
  }
  SpicyCell::{ image_id, cell_type, x, y, }
}

///|
pub struct SpicyImageMetadata {
  image_id : String
  condition : String
  subject : String
  covariates : Map[String, Double]
} derive(Debug)
pub extend SpicyImageMetadata with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn SpicyImageMetadata::create(
  image_id : String,
  condition : String,
  subject? : String = "",
  covariates? : Map[String, Double] = Map([]),
) -> SpicyImageMetadata raise SpicyError {
  if image_id == "" {
    raise SpicyError("spicyR metadata image ID must not be empty")
  }
  if condition == "" {
    raise SpicyError("spicyR condition must not be empty")
  }
  let copied : Map[String, Double] = Map([], capacity=covariates.length())
  for name in covariates.keys() {
    if name == "" {
      raise SpicyError("spicyR covariate names must not be empty")
    }
    let value = covariates[name]
    if !spicy_is_finite(value) {
      raise SpicyError("spicyR covariates must be finite")
    }
    copied[name] = value
  }
  SpicyImageMetadata::{ image_id, condition, subject, covariates: copied, }
}

///|
pub struct SpicyPair {
  from_cell_type : String
  to_cell_type : String
} derive(Eq, Debug)
pub extend SpicyPair with Eq::{equal, not_equal}
pub extend SpicyPair with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn SpicyPair::create(
  from_cell_type : String,
  to_cell_type : String,
) -> SpicyPair raise SpicyError {
  if from_cell_type == "" || to_cell_type == "" {
    raise SpicyError("spicyR pair cell types must not be empty")
  }
  SpicyPair::{ from_cell_type, to_cell_type, }
}

///|
pub fn SpicyPair::label(self : SpicyPair) -> String {
  self.from_cell_type + "__" + self.to_cell_type
}

///|
pub struct SpicyConfig {
  radii : Array[Double]
  edge_correction : Bool
  include_zero_cells : Bool
  use_weights : Bool
  weight_factor : Double
  window_padding : Double
  fdr_threshold : Double
  reference_condition : String
  max_iterations : Int
  tolerance : Double
} derive(Debug)
pub extend SpicyConfig with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn SpicyConfig::create(
  radii? : Array[Double] = [20.0, 50.0, 100.0],
  edge_correction? : Bool = true,
  include_zero_cells? : Bool = false,
  use_weights? : Bool = true,
  weight_factor? : Double = 1.0,
  window_padding? : Double = 0.0,
  fdr_threshold? : Double = 0.05,
  reference_condition? : String = "",
  max_iterations? : Int = 100,
  tolerance? : Double = 0.005,
) -> SpicyConfig raise SpicyError {
  if radii.length() == 0 {
    raise SpicyError("spicyR radii must not be empty")
  }
  let copied = radii.copy()
  let mut previous = 0.0
  for index in 0..<copied.length() {
    let radius = copied[index]
    if !spicy_is_finite(radius) || radius <= 0.0 {
      raise SpicyError("spicyR radii must be finite and positive")
    }
    if index > 0 && radius <= previous {
      raise SpicyError("spicyR radii must be strictly increasing")
    }
    previous = radius
  }
  if !spicy_is_finite(weight_factor) || weight_factor < 0.0 {
    raise SpicyError("spicyR weight factor must be finite and non-negative")
  }
  if !spicy_is_finite(window_padding) || window_padding < 0.0 {
    raise SpicyError("spicyR window padding must be finite and non-negative")
  }
  if !spicy_is_finite(fdr_threshold) ||
    fdr_threshold <= 0.0 ||
    fdr_threshold > 1.0 {
    raise SpicyError("spicyR FDR threshold must be in (0, 1]")
  }
  if max_iterations <= 0 {
    raise SpicyError("spicyR max_iterations must be positive")
  }
  if !spicy_is_finite(tolerance) || tolerance <= 0.0 || tolerance >= 1.0 {
    raise SpicyError("spicyR tolerance must be in (0, 1)")
  }
  SpicyConfig::{
    radii: copied,
    edge_correction,
    include_zero_cells,
    use_weights,
    weight_factor,
    window_padding,
    fdr_threshold,
    reference_condition,
    max_iterations,
    tolerance,
  }
}

///|
pub fn SpicyConfig::default() -> SpicyConfig {
  SpicyConfig::{
    radii: [20.0, 50.0, 100.0],
    edge_correction: true,
    include_zero_cells: false,
    use_weights: true,
    weight_factor: 1.0,
    window_padding: 0.0,
    fdr_threshold: 0.05,
    reference_condition: "",
    max_iterations: 100,
    tolerance: 0.005,
  }
}

///|
pub struct SpicyAssociation {
  image_id : String
  from_cell_type : String
  to_cell_type : String
  radii : Array[Double]
  l_values : Array[Double]
  statistic : Double
  from_count : Int
  to_count : Int
  observed : Bool
  available : Bool
  weight : Double
} derive(Debug)
pub extend SpicyAssociation with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn SpicyAssociation::create(
  image_id : String,
  from_cell_type : String,
  to_cell_type : String,
  statistic : Double,
  from_count? : Int = 1,
  to_count? : Int = 1,
  observed? : Bool = true,
  available? : Bool = true,
) -> SpicyAssociation raise SpicyError {
  if image_id == "" {
    raise SpicyError("spicyR association image ID must not be empty")
  }
  if from_cell_type == "" || to_cell_type == "" {
    raise SpicyError("spicyR association cell types must not be empty")
  }
  if !spicy_is_finite(statistic) {
    raise SpicyError("spicyR association statistic must be finite")
  }
  if from_count < 0 || to_count < 0 {
    raise SpicyError("spicyR association cell counts must be non-negative")
  }
  SpicyAssociation::{
    image_id,
    from_cell_type,
    to_cell_type,
    radii: [],
    l_values: [],
    statistic,
    from_count,
    to_count,
    observed,
    available,
    weight: 1.0,
  }
}

///|
pub struct SpicyConditionResult {
  condition_level : String
  estimate : Double
  standard_error : Double
  statistic : Double
  degrees_of_freedom : Double
  p_value : Double
  adjusted_p_value : Double
  converged : Bool
} derive(Debug)
pub extend SpicyConditionResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct SpicyPairResult {
  pair : SpicyPair
  model : SpicyModelKind
  intercept : Double
  conditions : Array[SpicyConditionResult]
  random_intercept_variance : Double
  residual_variance : Double
  associations : Array[SpicyAssociation]
} derive(Debug)
pub extend SpicyPairResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct SpicyResult {
  reference_condition : String
  condition_levels : Array[String]
  image_ids : Array[String]
  pairs : Array[SpicyPairResult]
  config : SpicyConfig
} derive(Debug)
pub extend SpicyResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct SpicySpatialExperimentOutput {
  experiment : SpatialExperiment
  result : SpicyResult
}

///|
priv struct SpicyWindow {
  min_x : Double
  max_x : Double
  min_y : Double
  max_y : Double
  area : Double
}

///|
fn spicy_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn spicy_find_string(values : Array[String], target : String) -> Int {
  for index in 0..<values.length() {
    if values[index] == target {
      return index
    }
  }
  -1
}

///|
fn spicy_copy_string_map(source : Map[String, String]) -> Map[String, String] {
  let output : Map[String, String] = Map([], capacity=source.length())
  for key in source.keys() {
    output[key] = source[key]
  }
  output
}

///|
fn spicy_copy_config(config : SpicyConfig) -> SpicyConfig {
  SpicyConfig::{
    radii: config.radii.copy(),
    edge_correction: config.edge_correction,
    include_zero_cells: config.include_zero_cells,
    use_weights: config.use_weights,
    weight_factor: config.weight_factor,
    window_padding: config.window_padding,
    fdr_threshold: config.fdr_threshold,
    reference_condition: config.reference_condition,
    max_iterations: config.max_iterations,
    tolerance: config.tolerance,
  }
}

///|
fn spicy_copy_association(
  value : SpicyAssociation,
  weight : Double,
) -> SpicyAssociation {
  SpicyAssociation::{
    image_id: value.image_id,
    from_cell_type: value.from_cell_type,
    to_cell_type: value.to_cell_type,
    radii: value.radii.copy(),
    l_values: value.l_values.copy(),
    statistic: value.statistic,
    from_count: value.from_count,
    to_count: value.to_count,
    observed: value.observed,
    available: value.available,
    weight,
  }
}

///|
fn spicy_image_cells(
  cells : Array[SpicyCell],
  image_id : String,
) -> Array[SpicyCell] {
  let output : Array[SpicyCell] = []
  for cell in cells {
    if cell.image_id == image_id {
      output.push(cell)
    }
  }
  output
}

///|
fn spicy_window(
  cells : Array[SpicyCell],
  padding : Double,
) -> SpicyWindow raise SpicyError {
  if cells.length() == 0 {
    raise SpicyError("spicyR image must contain cells")
  }
  let mut min_x = cells[0].x
  let mut max_x = cells[0].x
  let mut min_y = cells[0].y
  let mut max_y = cells[0].y
  for index in 1..<cells.length() {
    min_x = min_x.min(cells[index].x)
    max_x = max_x.max(cells[index].x)
    min_y = min_y.min(cells[index].y)
    max_y = max_y.max(cells[index].y)
  }
  min_x = min_x - padding
  max_x = max_x + padding
  min_y = min_y - padding
  max_y = max_y + padding
  let width = max_x - min_x
  let height = max_y - min_y
  if width <= 0.0 || height <= 0.0 {
    raise SpicyError("spicyR image window must have positive width and height")
  }
  SpicyWindow::{ min_x, max_x, min_y, max_y, area: width * height, }
}

///|
fn spicy_effective_radii(
  config : SpicyConfig,
  window : SpicyWindow,
) -> Array[Double] {
  let width = window.max_x - window.min_x
  let height = window.max_y - window.min_y
  let maximum = width.min(height) / 2.01
  let output : Array[Double] = []
  for radius in config.radii {
    let effective = radius.min(maximum)
    if effective > 0.0 {
      if output.length() == 0 ||
        (output[output.length() - 1] - effective).abs() > 1.0e-12 {
        output.push(effective)
      }
    }
  }
  output
}

///|
fn spicy_circle_vertical_length(
  relative_x : Double,
  center_y : Double,
  radius : Double,
  window : SpicyWindow,
) -> Double {
  let squared = radius * radius - relative_x * relative_x
  if squared <= 0.0 {
    return 0.0
  }
  let half_height = squared.sqrt()
  let lower = (center_y - half_height).max(window.min_y)
  let upper = (center_y + half_height).min(window.max_y)
  (upper - lower).max(0.0)
}

///|
fn spicy_circle_window_fraction(
  x : Double,
  y : Double,
  radius : Double,
  window : SpicyWindow,
) -> Double {
  let left = (-radius).max(window.min_x - x)
  let right = radius.min(window.max_x - x)
  if right <= left {
    return 1.0e-8
  }
  // Simpson integration gives a deterministic rectangle/disc intersection.
  let slices = 96
  let step = (right - left) / slices.to_double()
  let mut total = 0.0
  for index in 0..<=slices {
    let relative_x = left + index.to_double() * step
    let height = spicy_circle_vertical_length(relative_x, y, radius, window)
    let coefficient = if index == 0 || index == slices {
      1.0
    } else if index % 2 == 0 {
      2.0
    } else {
      4.0
    }
    total = total + coefficient * height
  }
  let intersection = total * step / 3.0
  (intersection / (3.141592653589793 * radius * radius)).max(1.0e-8).min(1.0)
}

///|
fn spicy_pair_indices(
  cells : Array[SpicyCell],
  cell_type : String,
) -> Array[Int] {
  let output : Array[Int] = []
  for index in 0..<cells.length() {
    if cells[index].cell_type == cell_type {
      output.push(index)
    }
  }
  output
}

///|
fn spicy_pair_association_for_image(
  image : SpicyImageMetadata,
  cells : Array[SpicyCell],
  pair : SpicyPair,
  config : SpicyConfig,
) -> SpicyAssociation raise SpicyError {
  let window = spicy_window(cells, config.window_padding)
  let radii = spicy_effective_radii(config, window)
  if radii.length() == 0 {
    raise SpicyError("spicyR effective radii must not be empty")
  }
  let from_indices = spicy_pair_indices(cells, pair.from_cell_type)
  let to_indices = spicy_pair_indices(cells, pair.to_cell_type)
  let same_type = pair.from_cell_type == pair.to_cell_type
  let present = from_indices.length() > 0 && to_indices.length() > 0
  let observed = present && (!same_type || from_indices.length() > 1)
  let available = present || config.include_zero_cells
  let l_values = Array::make(radii.length(), 0.0)
  if observed {
    let denominator = from_indices.length().to_double() *
      to_indices.length().to_double()
    for radius_index in 0..<radii.length() {
      let radius = radii[radius_index]
      let edge_weights = Array::make(from_indices.length(), 1.0)
      if config.edge_correction {
        for source_index in 0..<from_indices.length() {
          let source = cells[from_indices[source_index]]
          edge_weights[source_index] = 1.0 /
            spicy_circle_window_fraction(source.x, source.y, radius, window)
        }
      }
      let mut weighted_pairs = 0.0
      for source_index in 0..<from_indices.length() {
        let source_cell_index = from_indices[source_index]
        let source = cells[source_cell_index]
        for target_cell_index in to_indices {
          if source_cell_index != target_cell_index {
            let target = cells[target_cell_index]
            let delta_x = source.x - target.x
            let delta_y = source.y - target.y
            let distance_squared = delta_x * delta_x + delta_y * delta_y
            if distance_squared <= radius * radius + 1.0e-12 {
              weighted_pairs = weighted_pairs + edge_weights[source_index]
            }
          }
        }
      }
      let k_value = window.area * weighted_pairs / denominator
      l_values[radius_index] = (k_value.max(0.0) / 3.141592653589793).sqrt()
    }
  }
  let mut statistic = 0.0
  for index in 0..<radii.length() {
    statistic = statistic + l_values[index] - radii[index]
  }
  SpicyAssociation::{
    image_id: image.image_id,
    from_cell_type: pair.from_cell_type,
    to_cell_type: pair.to_cell_type,
    radii,
    l_values,
    statistic,
    from_count: from_indices.length(),
    to_count: to_indices.length(),
    observed,
    available,
    weight: 1.0,
  }
}

///|
fn spicy_validate_metadata(
  cells : Array[SpicyCell],
  metadata : Array[SpicyImageMetadata],
) -> Unit raise SpicyError {
  if cells.length() == 0 {
    raise SpicyError("spicyR cells must not be empty")
  }
  if metadata.length() < 3 {
    raise SpicyError("spicyR requires at least three images")
  }
  let seen : Map[String, Bool] = Map([], capacity=metadata.length())
  for image in metadata {
    if seen.contains(image.image_id) {
      raise SpicyError("spicyR image metadata IDs must be unique")
    }
    seen[image.image_id] = true
    let mut found = false
    for cell in cells {
      if cell.image_id == image.image_id {
        found = true
        break
      }
    }
    if !found {
      raise SpicyError("spicyR metadata image has no cells: " + image.image_id)
    }
  }
  for cell in cells {
    if !seen.contains(cell.image_id) {
      raise SpicyError("spicyR cell references unknown image: " + cell.image_id)
    }
  }
}

///|
fn spicy_cell_types(cells : Array[SpicyCell]) -> Array[String] {
  let output : Array[String] = []
  for cell in cells {
    if !output.contains(cell.cell_type) {
      output.push(cell.cell_type)
    }
  }
  output
}

///|
fn spicy_resolve_pairs(
  cells : Array[SpicyCell],
  requested : Array[SpicyPair],
) -> Array[SpicyPair] raise SpicyError {
  let cell_types = spicy_cell_types(cells)
  let output : Array[SpicyPair] = []
  if requested.length() == 0 {
    for to_type in cell_types {
      for from_type in cell_types {
        output.push(SpicyPair::{
          from_cell_type: from_type,
          to_cell_type: to_type,
        })
      }
    }
    return output
  }
  let labels : Array[String] = []
  for pair in requested {
    if !cell_types.contains(pair.from_cell_type) ||
      !cell_types.contains(pair.to_cell_type) {
      raise SpicyError("spicyR pair cell type is absent from the data")
    }
    let label = pair.label()
    if labels.contains(label) {
      raise SpicyError("spicyR pairs must not contain duplicates")
    }
    labels.push(label)
    output.push(pair)
  }
  output
}

///|
pub fn spicy_pairwise(
  cells : Array[SpicyCell],
  metadata : Array[SpicyImageMetadata],
  pair : SpicyPair,
  config? : SpicyConfig = SpicyConfig::default(),
) -> Array[SpicyAssociation] raise SpicyError {
  spicy_validate_metadata(cells, metadata)
  let types = spicy_cell_types(cells)
  if !types.contains(pair.from_cell_type) || !types.contains(pair.to_cell_type) {
    raise SpicyError("spicyR pair cell type is absent from the data")
  }
  let output : Array[SpicyAssociation] = []
  for image in metadata {
    output.push(
      spicy_pair_association_for_image(
        image,
        spicy_image_cells(cells, image.image_id),
        pair,
        config,
      ),
    )
  }
  output
}

///|
fn spicy_pair_key(association : SpicyAssociation) -> String {
  association.from_cell_type + "__" + association.to_cell_type
}

///|
fn spicy_effective_count(association : SpicyAssociation) -> Double {
  let from_count = association.from_count.to_double()
  let to_count = association.to_count.to_double()
  if association.from_cell_type == association.to_cell_type {
    (from_count * (from_count - 1.0)).max(1.0)
  } else if from_count + to_count <= 0.0 {
    1.0
  } else {
    (2.0 * from_count * to_count / (from_count + to_count)).max(1.0)
  }
}

///|
fn spicy_weight_associations(
  associations : Array[SpicyAssociation],
  config : SpicyConfig,
) -> Array[SpicyAssociation] {
  if !config.use_weights {
    let output : Array[SpicyAssociation] = []
    for value in associations {
      output.push(spicy_copy_association(value, 1.0))
    }
    return output
  }
  let raw = Array::make(associations.length(), 1.0)
  let mut total = 0.0
  let mut count = 0
  for index in 0..<associations.length() {
    if associations[index].available {
      raw[index] = spicy_effective_count(associations[index])
      total = total + raw[index]
      count = count + 1
    }
  }
  let mean = if count == 0 { 1.0 } else { total / count.to_double() }
  let output : Array[SpicyAssociation] = []
  for index in 0..<associations.length() {
    let weight = if associations[index].available {
      @math.pow((raw[index] / mean).max(1.0e-8), config.weight_factor)
    } else {
      1.0
    }
    output.push(spicy_copy_association(associations[index], weight))
  }
  output
}

///|
fn spicy_metadata_for_image(
  metadata : Array[SpicyImageMetadata],
  image_id : String,
) -> SpicyImageMetadata? {
  for value in metadata {
    if value.image_id == image_id {
      return Some(value)
    }
  }
  None
}

///|
fn spicy_condition_levels(
  metadata : Array[SpicyImageMetadata],
) -> Array[String] {
  let output : Array[String] = []
  for image in metadata {
    if !output.contains(image.condition) {
      output.push(image.condition)
    }
  }
  output
}

///|
fn spicy_model_kind(metadata : Array[SpicyImageMetadata]) -> SpicyModelKind {
  let subjects : Array[String] = []
  let mut all_present = true
  for image in metadata {
    if image.subject == "" {
      all_present = false
    } else if !subjects.contains(image.subject) {
      subjects.push(image.subject)
    }
  }
  if all_present &&
    subjects.length() >= 2 &&
    subjects.length() < metadata.length() {
    WeightedRandomIntercept
  } else {
    WeightedLinear
  }
}

///|
fn spicy_build_design(
  metadata : Array[SpicyImageMetadata],
  covariate_names : Array[String],
  reference : String,
  model : SpicyModelKind,
) -> VariancePartitionDesign raise SpicyError {
  let sample_names : Array[String] = []
  let conditions : Array[String] = []
  for image in metadata {
    sample_names.push(image.image_id)
    conditions.push(image.condition)
  }
  let fixed : Array[VariancePartitionFixedEffect] = []
  let condition_effect = vp_categorical_effect(
    "condition",
    conditions,
    reference~,
  ) catch {
    VariancePartitionError(message) =>
      raise SpicyError("spicyR condition design failed: " + message)
  }
  fixed.push(condition_effect)
  let seen_covariates : Array[String] = []
  for name in covariate_names {
    if name == "" || seen_covariates.contains(name) {
      raise SpicyError("spicyR covariate names must be non-empty and unique")
    }
    seen_covariates.push(name)
    let values : Array[Double] = []
    for image in metadata {
      match image.covariates.get(name) {
        Some(value) => values.push(value)
        None =>
          raise SpicyError(
            "spicyR covariate '" +
            name +
            "' is missing for image " +
            image.image_id,
          )
      }
    }
    let effect = vp_numeric_effect(name, values) catch {
      VariancePartitionError(message) =>
        raise SpicyError("spicyR covariate design failed: " + message)
    }
    fixed.push(effect)
  }
  let random : Array[VariancePartitionRandomEffect] = []
  if model == WeightedRandomIntercept {
    let subjects : Array[String] = []
    for image in metadata {
      subjects.push(image.subject)
    }
    let effect = vp_random_effect("subject", subjects) catch {
      VariancePartitionError(message) =>
        raise SpicyError("spicyR subject design failed: " + message)
    }
    random.push(effect)
  }
  vp_design(sample_names, fixed, random) catch {
    VariancePartitionError(message) =>
      raise SpicyError("spicyR model design failed: " + message)
  }
}

///|
fn spicy_bh_adjust(p_values : Array[Double]) -> Array[Double] {
  let count = p_values.length()
  let order : Array[Int] = []
  for index in 0..<count {
    order.push(index)
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    if p_values[left] < p_values[right] {
      -1
    } else if p_values[left] > p_values[right] {
      1
    } else if left < right {
      -1
    } else if left > right {
      1
    } else {
      0
    }
  })
  let adjusted = Array::make(count, 1.0)
  let mut running = 1.0
  let mut rank = count
  while rank > 0 {
    let original = order[rank - 1]
    let candidate = (p_values[original] * count.to_double() / rank.to_double()).min(
      1.0,
    )
    running = running.min(candidate)
    adjusted[original] = running.max(p_values[original])
    rank = rank - 1
  }
  adjusted
}

///|
fn spicy_fit_pair(
  pair : SpicyPair,
  associations : Array[SpicyAssociation],
  metadata : Array[SpicyImageMetadata],
  covariate_names : Array[String],
  condition_levels : Array[String],
  reference : String,
  config : SpicyConfig,
) -> SpicyPairResult raise SpicyError {
  let weighted = spicy_weight_associations(associations, config)
  let selected_metadata : Array[SpicyImageMetadata] = []
  let values : Array[Double] = []
  let weights : Array[Double] = []
  let selected_associations : Array[SpicyAssociation] = []
  for association in weighted {
    if association.available {
      let image = match
        spicy_metadata_for_image(metadata, association.image_id) {
        Some(value) => value
        None =>
          raise SpicyError(
            "spicyR association references unknown image: " +
            association.image_id,
          )
      }
      selected_metadata.push(image)
      values.push(association.statistic)
      weights.push(association.weight)
      selected_associations.push(association)
    }
  }
  if selected_associations.length() < 3 {
    raise SpicyError(
      "spicyR pair " + pair.label() + " has fewer than three usable images",
    )
  }
  let model = spicy_model_kind(selected_metadata)
  let design = spicy_build_design(
    selected_metadata, covariate_names, reference, model,
  )
  let fit = vp_fit_gene(
    values,
    design,
    gene_name=pair.label(),
    weights~,
    reml=true,
    max_iterations=config.max_iterations,
    tolerance=config.tolerance,
  ) catch {
    VariancePartitionError(message) =>
      raise SpicyError(
        "spicyR model fit failed for " + pair.label() + ": " + message,
      )
  }
  let intercept_index = spicy_find_string(fit.coefficient_names, "(Intercept)")
  let intercept = if intercept_index < 0 {
    0.0
  } else {
    fit.coefficients[intercept_index]
  }
  let condition_results : Array[SpicyConditionResult] = []
  for level in condition_levels {
    if level != reference {
      let coefficient_name = "condition:" + level
      let contrast = vp_coefficient_contrast(
        design,
        level + "-" + reference,
        coefficient_name,
      ) catch {
        VariancePartitionError(message) =>
          raise SpicyError(
            "spicyR condition contrast failed for " +
            pair.label() +
            ": " +
            message,
          )
      }
      let result = dream(
        [values],
        design,
        contrast,
        gene_names=[pair.label()],
        weights=[weights],
        reml=true,
        max_iterations=config.max_iterations,
        tolerance=config.tolerance,
      ) catch {
        VariancePartitionError(message) =>
          raise SpicyError(
            "spicyR condition test failed for " + pair.label() + ": " + message,
          )
      }
      let tested = result.genes[0]
      condition_results.push(SpicyConditionResult::{
        condition_level: level,
        estimate: tested.estimate,
        standard_error: tested.standard_error,
        statistic: tested.statistic,
        degrees_of_freedom: tested.degrees_of_freedom,
        p_value: tested.p_value,
        adjusted_p_value: tested.p_value,
        converged: tested.converged,
      })
    }
  }
  let random_variance = if fit.random_variances.length() == 0 {
    0.0
  } else {
    fit.random_variances[0]
  }
  SpicyPairResult::{
    pair,
    model,
    intercept,
    conditions: condition_results,
    random_intercept_variance: random_variance,
    residual_variance: fit.residual_variance,
    associations: weighted,
  }
}

///|
fn spicy_adjust_pair_results(
  pairs : Array[SpicyPairResult],
  condition_levels : Array[String],
  reference : String,
) -> Array[SpicyPairResult] {
  let adjusted_by_level : Map[String, Array[Double]] = Map([])
  for level in condition_levels {
    if level != reference {
      let p_values : Array[Double] = []
      for pair in pairs {
        let mut found = 1.0
        for condition in pair.conditions {
          if condition.condition_level == level {
            found = condition.p_value
            break
          }
        }
        p_values.push(found)
      }
      adjusted_by_level[level] = spicy_bh_adjust(p_values)
    }
  }
  let output : Array[SpicyPairResult] = []
  for pair_index in 0..<pairs.length() {
    let pair = pairs[pair_index]
    let conditions : Array[SpicyConditionResult] = []
    for condition in pair.conditions {
      conditions.push(SpicyConditionResult::{
        condition_level: condition.condition_level,
        estimate: condition.estimate,
        standard_error: condition.standard_error,
        statistic: condition.statistic,
        degrees_of_freedom: condition.degrees_of_freedom,
        p_value: condition.p_value,
        adjusted_p_value: adjusted_by_level[condition.condition_level][pair_index],
        converged: condition.converged,
      })
    }
    output.push(SpicyPairResult::{
      pair: pair.pair,
      model: pair.model,
      intercept: pair.intercept,
      conditions,
      random_intercept_variance: pair.random_intercept_variance,
      residual_variance: pair.residual_variance,
      associations: pair.associations,
    })
  }
  output
}

///|
pub fn spicy_fit_associations(
  associations : Array[SpicyAssociation],
  metadata : Array[SpicyImageMetadata],
  covariate_names? : Array[String] = [],
  config? : SpicyConfig = SpicyConfig::default(),
) -> SpicyResult raise SpicyError {
  if associations.length() == 0 {
    raise SpicyError("spicyR associations must not be empty")
  }
  if metadata.length() < 3 {
    raise SpicyError("spicyR requires at least three images")
  }
  let metadata_ids : Array[String] = []
  for image in metadata {
    if metadata_ids.contains(image.image_id) {
      raise SpicyError("spicyR image metadata IDs must be unique")
    }
    metadata_ids.push(image.image_id)
  }
  let condition_levels = spicy_condition_levels(metadata)
  if condition_levels.length() < 2 {
    raise SpicyError("spicyR requires at least two condition levels")
  }
  let reference = if config.reference_condition == "" {
    condition_levels[0]
  } else {
    config.reference_condition
  }
  if !condition_levels.contains(reference) {
    raise SpicyError("spicyR reference condition is absent from metadata")
  }
  let pair_labels : Array[String] = []
  let pairs : Array[SpicyPair] = []
  for association in associations {
    if !metadata_ids.contains(association.image_id) {
      raise SpicyError(
        "spicyR association references unknown image: " + association.image_id,
      )
    }
    let key = spicy_pair_key(association)
    if !pair_labels.contains(key) {
      pair_labels.push(key)
      pairs.push(SpicyPair::{
        from_cell_type: association.from_cell_type,
        to_cell_type: association.to_cell_type,
      })
    }
  }
  let raw_results : Array[SpicyPairResult] = []
  for pair_index in 0..<pairs.length() {
    let pair = pairs[pair_index]
    let selected : Array[SpicyAssociation] = []
    let seen_images : Array[String] = []
    for association in associations {
      if spicy_pair_key(association) == pair_labels[pair_index] {
        if seen_images.contains(association.image_id) {
          raise SpicyError(
            "spicyR pair has duplicate association for image " +
            association.image_id,
          )
        }
        seen_images.push(association.image_id)
        selected.push(association)
      }
    }
    raw_results.push(
      spicy_fit_pair(
        pair, selected, metadata, covariate_names, condition_levels, reference, config,
      ),
    )
  }
  SpicyResult::{
    reference_condition: reference,
    condition_levels: condition_levels.copy(),
    image_ids: metadata_ids,
    pairs: spicy_adjust_pair_results(raw_results, condition_levels, reference),
    config: spicy_copy_config(config),
  }
}

///|
pub fn spicyr(
  cells : Array[SpicyCell],
  metadata : Array[SpicyImageMetadata],
  pairs? : Array[SpicyPair] = [],
  covariate_names? : Array[String] = [],
  config? : SpicyConfig = SpicyConfig::default(),
) -> SpicyResult raise SpicyError {
  spicy_validate_metadata(cells, metadata)
  let selected_pairs = spicy_resolve_pairs(cells, pairs)
  let associations : Array[SpicyAssociation] = []
  for pair in selected_pairs {
    for image in metadata {
      associations.push(
        spicy_pair_association_for_image(
          image,
          spicy_image_cells(cells, image.image_id),
          pair,
          config,
        ),
      )
    }
  }
  spicy_fit_associations(associations, metadata, covariate_names~, config~)
}

///|
pub fn SpicyPairResult::condition(
  self : SpicyPairResult,
  level : String,
) -> SpicyConditionResult? {
  for value in self.conditions {
    if value.condition_level == level {
      return Some(value)
    }
  }
  None
}

///|
pub fn SpicyPairResult::association(
  self : SpicyPairResult,
  image_id : String,
) -> SpicyAssociation? {
  for value in self.associations {
    if value.image_id == image_id {
      return Some(value)
    }
  }
  None
}

///|
pub fn SpicyResult::pair(
  self : SpicyResult,
  from_cell_type : String,
  to_cell_type : String,
) -> SpicyPairResult? {
  for pair in self.pairs {
    if pair.pair.from_cell_type == from_cell_type &&
      pair.pair.to_cell_type == to_cell_type {
      return Some(pair)
    }
  }
  None
}

///|
pub fn SpicyResult::n_pairs(self : SpicyResult) -> Int {
  self.pairs.length()
}

///|
pub fn SpicyResult::n_images(self : SpicyResult) -> Int {
  self.image_ids.length()
}

///|
pub fn SpicyResult::significant(
  self : SpicyResult,
  condition_level : String,
  threshold? : Double = -1.0,
) -> Array[SpicyPairResult] raise SpicyError {
  let selected_threshold = if threshold < 0.0 {
    self.config.fdr_threshold
  } else {
    threshold
  }
  if !spicy_is_finite(selected_threshold) ||
    selected_threshold < 0.0 ||
    selected_threshold > 1.0 {
    raise SpicyError("spicyR significance threshold must be in [0, 1]")
  }
  if condition_level == self.reference_condition ||
    !self.condition_levels.contains(condition_level) {
    raise SpicyError("spicyR requested condition contrast is absent")
  }
  let output : Array[SpicyPairResult] = []
  for pair in self.pairs {
    match pair.condition(condition_level) {
      Some(value) =>
        if value.adjusted_p_value <= selected_threshold {
          output.push(pair)
        }
      None => ()
    }
  }
  output
}

///|
pub fn SpicyResult::top(
  self : SpicyResult,
  condition_level : String,
  limit? : Int = 10,
) -> Array[SpicyPairResult] raise SpicyError {
  if limit < 0 {
    raise SpicyError("spicyR top-pair limit must be non-negative")
  }
  if condition_level == self.reference_condition ||
    !self.condition_levels.contains(condition_level) {
    raise SpicyError("spicyR requested condition contrast is absent")
  }
  let output = self.pairs.copy()
  output.sort_by(fn(left : SpicyPairResult, right : SpicyPairResult) -> Int {
    let left_p = match left.condition(condition_level) {
      Some(value) => value.p_value
      None => 1.0
    }
    let right_p = match right.condition(condition_level) {
      Some(value) => value.p_value
      None => 1.0
    }
    if left_p < right_p {
      -1
    } else if left_p > right_p {
      1
    } else {
      left.pair.label().compare(right.pair.label())
    }
  })
  let count = limit.min(output.length())
  let selected : Array[SpicyPairResult] = []
  for index in 0..<count {
    selected.push(output[index])
  }
  selected
}

///|
pub fn SpicyResult::summary(self : SpicyResult) -> String {
  let mut significant = 0
  for pair in self.pairs {
    for condition in pair.conditions {
      if condition.adjusted_p_value <= self.config.fdr_threshold {
        significant = significant + 1
      }
    }
  }
  "spicyR(images=" +
  self.n_images().to_string() +
  ", pairs=" +
  self.n_pairs().to_string() +
  ", reference=" +
  self.reference_condition +
  ", significant_contrasts=" +
  significant.to_string() +
  ")"
}

///|
fn spicy_copy_spatial_experiment(
  experiment : SpatialExperiment,
) -> SpatialExperiment {
  let output = SpatialExperiment::new()
  for assay_name in experiment.assay.keys() {
    let matrix : Array[Array[Double]] = []
    for row in experiment.assay[assay_name] {
      matrix.push(row.copy())
    }
    output.assay[assay_name] = matrix
  }
  for row in experiment.row_data {
    output.row_data.push(spicy_copy_string_map(row))
  }
  for column in experiment.col_data {
    output.col_data.push(spicy_copy_string_map(column))
  }
  for coordinate in experiment.spatial_coords {
    output.spatial_coords.push(coordinate)
  }
  for image in experiment.images {
    output.images.push(image)
  }
  for key in experiment.metadata.keys() {
    output.metadata[key] = experiment.metadata[key]
  }
  output
}

///|
fn spicy_required_column(
  column : Map[String, String],
  key : String,
  row : Int,
) -> String raise SpicyError {
  match column.get(key) {
    Some(value) =>
      if value == "" {
        raise SpicyError(
          "spicyR SpatialExperiment column '" +
          key +
          "' is empty at cell " +
          row.to_string(),
        )
      } else {
        value
      }
    None =>
      raise SpicyError(
        "spicyR SpatialExperiment is missing column '" + key + "'",
      )
  }
}

///|
pub fn spicyr_spatial_experiment(
  experiment : SpatialExperiment,
  image_id_key? : String = "imageID",
  cell_type_key? : String = "cellType",
  condition_key? : String = "condition",
  subject_key? : String = "",
  covariate_keys? : Array[String] = [],
  pairs? : Array[SpicyPair] = [],
  config? : SpicyConfig = SpicyConfig::default(),
) -> SpicySpatialExperimentOutput raise SpicyError {
  if experiment.col_data.length() == 0 {
    raise SpicyError("spicyR SpatialExperiment col_data must not be empty")
  }
  if experiment.col_data.length() != experiment.spatial_coords.length() {
    raise SpicyError(
      "spicyR SpatialExperiment coordinates must match col_data rows",
    )
  }
  if image_id_key == "" || cell_type_key == "" || condition_key == "" {
    raise SpicyError("spicyR SpatialExperiment key names must not be empty")
  }
  let cells : Array[SpicyCell] = []
  let metadata : Array[SpicyImageMetadata] = []
  let metadata_indices : Map[String, Int] = Map([])
  for row in 0..<experiment.col_data.length() {
    let column = experiment.col_data[row]
    let image_id = spicy_required_column(column, image_id_key, row)
    let cell_type = spicy_required_column(column, cell_type_key, row)
    let condition = spicy_required_column(column, condition_key, row)
    let subject = if subject_key == "" {
      ""
    } else {
      spicy_required_column(column, subject_key, row)
    }
    let covariates : Map[String, Double] = Map(
      [],
      capacity=covariate_keys.length(),
    )
    for key in covariate_keys {
      let text = spicy_required_column(column, key, row)
      let value = match parse_double(text) {
        Some(parsed) => parsed
        None =>
          raise SpicyError(
            "spicyR SpatialExperiment covariate '" + key + "' must be numeric",
          )
      }
      if !spicy_is_finite(value) {
        raise SpicyError("spicyR SpatialExperiment covariates must be finite")
      }
      covariates[key] = value
    }
    cells.push(
      SpicyCell::create(
        image_id,
        cell_type,
        experiment.spatial_coords[row].x,
        experiment.spatial_coords[row].y,
      ),
    )
    match metadata_indices.get(image_id) {
      Some(index) => {
        let existing = metadata[index]
        if existing.condition != condition || existing.subject != subject {
          raise SpicyError(
            "spicyR image metadata must be constant within each image",
          )
        }
        for key in covariate_keys {
          match existing.covariates.get(key) {
            Some(value) =>
              if value != covariates[key] {
                raise SpicyError(
                  "spicyR covariates must be constant within each image",
                )
              }
            None =>
              raise SpicyError(
                "spicyR covariate is missing from existing image metadata",
              )
          }
        }
      }
      None => {
        metadata_indices[image_id] = metadata.length()
        metadata.push(
          SpicyImageMetadata::create(image_id, condition, subject~, covariates~),
        )
      }
    }
  }
  let result = spicyr(
    cells,
    metadata,
    pairs~,
    covariate_names=covariate_keys,
    config~,
  )
  let enriched = spicy_copy_spatial_experiment(experiment)
  enriched.metadata["spicyr_images"] = result.n_images().to_string()
  enriched.metadata["spicyr_pairs"] = result.n_pairs().to_string()
  enriched.metadata["spicyr_reference"] = result.reference_condition
  let mut significant = 0
  for pair in result.pairs {
    for condition in pair.conditions {
      if condition.adjusted_p_value <= config.fdr_threshold {
        significant = significant + 1
      }
    }
  }
  enriched.metadata["spicyr_significant_contrasts"] = significant.to_string()
  SpicySpatialExperimentOutput::{ experiment: enriched, result, }
}

///|
pub fn spicyr_example_data() -> (
  Array[SpicyCell],
  Array[SpicyImageMetadata],
  Array[SpicyPair],
) {
  let cells : Array[SpicyCell] = []
  let metadata : Array[SpicyImageMetadata] = []
  let anchors = [
    (2.0, 2.0),
    (2.4, 2.0),
    (2.0, 2.4),
    (2.4, 2.4),
    (1.8, 2.2),
    (2.6, 2.2),
    (0.0, 0.0),
    (10.0, 10.0),
  ]
  for subject_index in 0..<6 {
    let condition = if subject_index < 3 { "control" } else { "treated" }
    let subject = "subject_" + (subject_index + 1).to_string()
    for field in 0..<2 {
      let image_id = subject + "_field_" + (field + 1).to_string()
      let covariates : Map[String, Double] = Map([
        ("age", 40.0 + subject_index.to_double() * 3.0),
      ])
      metadata.push(SpicyImageMetadata::{
        image_id,
        condition,
        subject,
        covariates,
      })
      let jitter = subject_index.to_double() * 0.025 + field.to_double() * 0.015
      for point in anchors {
        cells.push(SpicyCell::{
          image_id,
          cell_type: "T_cell",
          x: point.0 + jitter,
          y: point.1,
        })
        let shift = if condition == "control" { 4.0 } else { 0.35 }
        cells.push(SpicyCell::{
          image_id,
          cell_type: "Tumour",
          x: (point.0 + shift).min(10.0) + jitter,
          y: (point.1 + shift).min(10.0),
        })
      }
    }
  }
  let pairs = [SpicyPair::{ from_cell_type: "T_cell", to_cell_type: "Tumour", }]
  (cells, metadata, pairs)
}
