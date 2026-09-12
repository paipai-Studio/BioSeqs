// Advanced multi-sample, multi-cluster differential-state analysis inspired
// by Bioconductor muscat 1.27.4.
//
// Matrices follow the Bioconductor gene x cell convention. Pseudobulks are
// represented as cluster x gene x sample arrays. The implementation is
// dependency-free: it uses a negative-binomial IRLS model with empirical
// Bayes dispersion shrinkage instead of dispatching to edgeR or DESeq2.

///|
pub suberror MuscatAdvancedError {
  MuscatAdvancedError(String)
}

///|
pub enum MuscatAdvancedAggregation {
  SumCounts
  MeanExpression
  MedianExpression
  ProportionDetected
  NumberDetected
} derive(Eq, Debug)

///|
pub fn muscat_sum_counts() -> MuscatAdvancedAggregation {
  SumCounts
}

///|
pub fn muscat_mean_expression() -> MuscatAdvancedAggregation {
  MeanExpression
}

///|
pub fn muscat_median_expression() -> MuscatAdvancedAggregation {
  MedianExpression
}

///|
pub fn muscat_proportion_detected() -> MuscatAdvancedAggregation {
  ProportionDetected
}

///|
pub fn muscat_number_detected() -> MuscatAdvancedAggregation {
  NumberDetected
}

///|
pub struct MuscatAdvancedConfig {
  min_cells : Int
  min_count : Double
  min_samples : Int
  max_iterations : Int
  tolerance : Double
  minimum_dispersion : Double
  maximum_dispersion : Double
  dispersion_prior_df : Double
  ridge : Double
  fdr_threshold : Double
  lfc_threshold : Double
  detection_filter : Double
} derive(Debug)

///|
pub fn MuscatAdvancedConfig::create(
  min_cells? : Int = 10,
  min_count? : Double = 1.0,
  min_samples? : Int = 2,
  max_iterations? : Int = 100,
  tolerance? : Double = 1.0e-8,
  minimum_dispersion? : Double = 1.0e-8,
  maximum_dispersion? : Double = 100.0,
  dispersion_prior_df? : Double = 10.0,
  ridge? : Double = 1.0e-8,
  fdr_threshold? : Double = 0.05,
  lfc_threshold? : Double = 0.0,
  detection_filter? : Double = 0.9,
) -> MuscatAdvancedConfig raise MuscatAdvancedError {
  if min_cells < 0 {
    raise MuscatAdvancedError("min_cells must be non-negative")
  }
  if !muscat_adv_finite(min_count) || min_count < 0.0 {
    raise MuscatAdvancedError("min_count must be finite and non-negative")
  }
  if min_samples < 1 {
    raise MuscatAdvancedError("min_samples must be positive")
  }
  if max_iterations < 1 {
    raise MuscatAdvancedError("max_iterations must be positive")
  }
  if !muscat_adv_finite(tolerance) || tolerance <= 0.0 {
    raise MuscatAdvancedError("tolerance must be finite and positive")
  }
  if !muscat_adv_finite(minimum_dispersion) || minimum_dispersion <= 0.0 {
    raise MuscatAdvancedError("minimum_dispersion must be finite and positive")
  }
  if !muscat_adv_finite(maximum_dispersion) ||
    maximum_dispersion < minimum_dispersion {
    raise MuscatAdvancedError(
      "maximum_dispersion must be finite and no smaller than the minimum",
    )
  }
  if !muscat_adv_finite(dispersion_prior_df) || dispersion_prior_df < 0.0 {
    raise MuscatAdvancedError(
      "dispersion_prior_df must be finite and non-negative",
    )
  }
  if !muscat_adv_finite(ridge) || ridge < 0.0 {
    raise MuscatAdvancedError("ridge must be finite and non-negative")
  }
  if !muscat_adv_finite(fdr_threshold) ||
    fdr_threshold <= 0.0 ||
    fdr_threshold > 1.0 {
    raise MuscatAdvancedError("fdr_threshold must be in (0, 1]")
  }
  if !muscat_adv_finite(lfc_threshold) || lfc_threshold < 0.0 {
    raise MuscatAdvancedError("lfc_threshold must be finite and non-negative")
  }
  if !muscat_adv_finite(detection_filter) ||
    detection_filter <= 0.0 ||
    detection_filter > 1.0 {
    raise MuscatAdvancedError("detection_filter must be in (0, 1]")
  }
  MuscatAdvancedConfig::{
    min_cells,
    min_count,
    min_samples,
    max_iterations,
    tolerance,
    minimum_dispersion,
    maximum_dispersion,
    dispersion_prior_df,
    ridge,
    fdr_threshold,
    lfc_threshold,
    detection_filter,
  }
}

///|
pub fn MuscatAdvancedConfig::default() -> MuscatAdvancedConfig {
  MuscatAdvancedConfig::{
    min_cells: 10,
    min_count: 1.0,
    min_samples: 2,
    max_iterations: 100,
    tolerance: 1.0e-8,
    minimum_dispersion: 1.0e-8,
    maximum_dispersion: 100.0,
    dispersion_prior_df: 10.0,
    ridge: 1.0e-8,
    fdr_threshold: 0.05,
    lfc_threshold: 0.0,
    detection_filter: 0.9,
  }
}

///|
pub struct MuscatAdvancedData {
  counts : Array[Array[Double]]
  gene_names : Array[String]
  cell_names : Array[String]
  sample_ids : Array[String]
  cluster_ids : Array[String]
  group_ids : Array[String]
  sample_names : Array[String]
  sample_groups : Array[String]
  cluster_names : Array[String]
  group_names : Array[String]
} derive(Debug)

///|
pub struct MuscatAdvancedPseudoBulk {
  values : Array[Array[Array[Double]]]
  gene_names : Array[String]
  sample_names : Array[String]
  sample_groups : Array[String]
  cluster_names : Array[String]
  cell_counts : Array[Array[Int]]
  library_sizes : Array[Array[Double]]
  aggregation : MuscatAdvancedAggregation
  scaled_cpm : Bool
} derive(Debug)

///|
pub struct MuscatAdvancedDesign {
  matrix : Array[Array[Double]]
  sample_names : Array[String]
  coefficient_names : Array[String]
} derive(Debug)

///|
pub struct MuscatAdvancedContrast {
  values : Array[Double]
  name : String
} derive(Debug)

///|
pub struct MuscatAdvancedGeneResult {
  gene_name : String
  cluster_name : String
  contrast_name : String
  mode : String
  tested : Bool
  log_fc : Double
  log_cpm : Double
  statistic : Double
  p_value : Double
  mut local_fdr : Double
  mut global_fdr : Double
  raw_dispersion : Double
  dispersion : Double
  mean_count : Double
  n_samples : Int
  converged : Bool
} derive(Debug)

///|
pub struct MuscatAdvancedResults {
  results : Array[MuscatAdvancedGeneResult]
  gene_names : Array[String]
  cluster_names : Array[String]
  contrast_names : Array[String]
  mode : String
  fdr_threshold : Double
} derive(Debug)

///|
pub struct MuscatStagewiseGeneResult {
  gene_name : String
  cluster_name : String
  contrast_name : String
  screen_p_value : Double
  mut screen_fdr : Double
  ds_p_value : Double
  dd_p_value : Double
  mut ds_confirmation : Double
  mut dd_confirmation : Double
  mut classification : String
} derive(Debug)

///|
pub struct MuscatStagewiseResults {
  results : Array[MuscatStagewiseGeneResult]
  gene_names : Array[String]
  cluster_names : Array[String]
  contrast_names : Array[String]
  alpha : Double
} derive(Debug)

///|
pub struct MuscatAdvancedSCEOutput {
  experiment : SingleCellExperiment
  sum_pseudobulk : MuscatAdvancedPseudoBulk
  detection_pseudobulk : MuscatAdvancedPseudoBulk
  ds_result : MuscatAdvancedResults
  dd_result : MuscatAdvancedResults
  stagewise_result : MuscatStagewiseResults
}

///|
priv struct MuscatAdvCoreFit {
  coefficients : Array[Double]
  covariance : Array[Array[Double]]
  fitted_values : Array[Double]
  converged : Bool
  fitted : Bool
}

///|
fn muscat_adv_finite(value : Double) -> Bool {
  value == value && value.abs() < 1.0e300
}

///|
fn muscat_adv_clamp(value : Double, lower : Double, upper : Double) -> Double {
  if value < lower {
    lower
  } else if value > upper {
    upper
  } else {
    value
  }
}

///|
fn muscat_adv_zero_matrix(rows : Int, columns : Int) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for _ in 0..<rows {
    output.push(Array::make(columns, 0.0))
  }
  output
}

///|
fn muscat_adv_copy_matrix(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in matrix {
    output.push(row.copy())
  }
  output
}

///|
fn muscat_adv_levels(
  values : Array[String],
  label : String,
) -> Array[String] raise MuscatAdvancedError {
  let output : Array[String] = []
  for value in values {
    if value.trim().length() == 0 {
      raise MuscatAdvancedError(label + " must not contain empty values")
    }
    if !output.contains(value) {
      output.push(value)
    }
  }
  output
}

///|
fn muscat_adv_names(
  supplied : Array[String],
  length : Int,
  prefix : String,
) -> Array[String] raise MuscatAdvancedError {
  if supplied.length() != 0 && supplied.length() != length {
    raise MuscatAdvancedError(prefix + " names do not match the data dimension")
  }
  let output : Array[String] = []
  for index in 0..<length {
    let name = if supplied.length() == 0 {
      prefix + (index + 1).to_string()
    } else {
      supplied[index]
    }
    if name.trim().length() == 0 {
      raise MuscatAdvancedError(prefix + " names must not be empty")
    }
    if output.contains(name) {
      raise MuscatAdvancedError(prefix + " names must be unique: " + name)
    }
    output.push(name)
  }
  output
}

///|
pub fn MuscatAdvancedData::create(
  counts : Array[Array[Double]],
  sample_ids : Array[String],
  cluster_ids : Array[String],
  group_ids : Array[String],
  gene_names? : Array[String] = [],
  cell_names? : Array[String] = [],
) -> MuscatAdvancedData raise MuscatAdvancedError {
  if counts.length() == 0 {
    raise MuscatAdvancedError("counts must contain at least one gene")
  }
  let cells = counts[0].length()
  if cells == 0 {
    raise MuscatAdvancedError("counts must contain at least one cell")
  }
  for row in counts {
    if row.length() != cells {
      raise MuscatAdvancedError("counts must be rectangular")
    }
    for value in row {
      if !muscat_adv_finite(value) ||
        value < 0.0 ||
        (value - value.round()).abs() > 1.0e-8 {
        raise MuscatAdvancedError(
          "counts must contain finite non-negative integer values",
        )
      }
    }
  }
  if sample_ids.length() != cells ||
    cluster_ids.length() != cells ||
    group_ids.length() != cells {
    raise MuscatAdvancedError(
      "sample, cluster and group identifiers must match the number of cells",
    )
  }
  let sample_names = muscat_adv_levels(sample_ids, "sample identifiers")
  let cluster_names = muscat_adv_levels(cluster_ids, "cluster identifiers")
  let group_names = muscat_adv_levels(group_ids, "group identifiers")
  if sample_names.length() < 2 {
    raise MuscatAdvancedError("at least two biological samples are required")
  }
  if group_names.length() < 2 {
    raise MuscatAdvancedError("at least two experimental groups are required")
  }
  let sample_groups : Array[String] = []
  for sample in sample_names {
    let mut assigned = ""
    for cell in 0..<cells {
      if sample_ids[cell] == sample {
        if assigned == "" {
          assigned = group_ids[cell]
        } else if assigned != group_ids[cell] {
          raise MuscatAdvancedError(
            "each sample must belong to exactly one experimental group",
          )
        }
      }
    }
    sample_groups.push(assigned)
  }
  MuscatAdvancedData::{
    counts: muscat_adv_copy_matrix(counts),
    gene_names: muscat_adv_names(gene_names, counts.length(), "Gene"),
    cell_names: muscat_adv_names(cell_names, cells, "Cell"),
    sample_ids: sample_ids.copy(),
    cluster_ids: cluster_ids.copy(),
    group_ids: group_ids.copy(),
    sample_names,
    sample_groups,
    cluster_names,
    group_names,
  }
}

///|
fn muscat_adv_median(values : Array[Double]) -> Double {
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
  if sorted.length() % 2 == 0 {
    (sorted[sorted.length() / 2 - 1] + sorted[sorted.length() / 2]) / 2.0
  } else {
    sorted[sorted.length() / 2]
  }
}

///|
pub fn muscat_aggregate_advanced(
  data : MuscatAdvancedData,
  aggregation? : MuscatAdvancedAggregation = SumCounts,
  scale_cpm? : Bool = false,
) -> MuscatAdvancedPseudoBulk raise MuscatAdvancedError {
  if scale_cpm && aggregation == NumberDetected {
    raise MuscatAdvancedError(
      "number-detected pseudobulks cannot be CPM scaled",
    )
  }
  let clusters = data.cluster_names.length()
  let samples = data.sample_names.length()
  let genes = data.gene_names.length()
  let values : Array[Array[Array[Double]]] = []
  let cell_counts : Array[Array[Int]] = []
  let library_sizes : Array[Array[Double]] = []
  for cluster in 0..<clusters {
    let cluster_values = muscat_adv_zero_matrix(genes, samples)
    let cluster_cells = Array::make(samples, 0)
    let cluster_libraries = Array::make(samples, 0.0)
    for sample in 0..<samples {
      let selected : Array[Int] = []
      for cell in 0..<data.cell_names.length() {
        if data.cluster_ids[cell] == data.cluster_names[cluster] &&
          data.sample_ids[cell] == data.sample_names[sample] {
          selected.push(cell)
        }
      }
      cluster_cells[sample] = selected.length()
      let mut library_size = 0.0
      for gene in 0..<genes {
        let observations : Array[Double] = []
        let mut detected = 0
        let mut sum = 0.0
        for cell in selected {
          let value = data.counts[gene][cell]
          observations.push(value)
          sum = sum + value
          if value > 0.0 {
            detected = detected + 1
          }
        }
        library_size = library_size + sum
        cluster_values[gene][sample] = match aggregation {
          SumCounts => sum
          MeanExpression =>
            if selected.length() == 0 {
              0.0
            } else {
              sum / selected.length().to_double()
            }
          MedianExpression => muscat_adv_median(observations)
          ProportionDetected =>
            if selected.length() == 0 {
              0.0
            } else {
              detected.to_double() / selected.length().to_double()
            }
          NumberDetected => detected.to_double()
        }
      }
      cluster_libraries[sample] = library_size
    }
    if scale_cpm {
      for sample in 0..<samples {
        let library_size = cluster_libraries[sample]
        if library_size > 0.0 {
          for gene in 0..<genes {
            cluster_values[gene][sample] = cluster_values[gene][sample] /
              library_size *
              1.0e6
          }
        }
      }
    }
    values.push(cluster_values)
    cell_counts.push(cluster_cells)
    library_sizes.push(cluster_libraries)
  }
  MuscatAdvancedPseudoBulk::{
    values,
    gene_names: data.gene_names.copy(),
    sample_names: data.sample_names.copy(),
    sample_groups: data.sample_groups.copy(),
    cluster_names: data.cluster_names.copy(),
    cell_counts,
    library_sizes,
    aggregation,
    scaled_cpm: scale_cpm,
  }
}

///|
fn muscat_adv_matrix_rank(matrix : Array[Array[Double]]) -> Int {
  if matrix.length() == 0 || matrix[0].length() == 0 {
    return 0
  }
  let work = muscat_adv_copy_matrix(matrix)
  let rows = work.length()
  let columns = work[0].length()
  let mut rank = 0
  let mut column = 0
  while column < columns && rank < rows {
    let mut pivot = rank
    let mut maximum = work[pivot][column].abs()
    for candidate in (rank + 1)..<rows {
      if work[candidate][column].abs() > maximum {
        pivot = candidate
        maximum = work[candidate][column].abs()
      }
    }
    if maximum > 1.0e-10 {
      if pivot != rank {
        let temporary = work[rank]
        work[rank] = work[pivot]
        work[pivot] = temporary
      }
      let divisor = work[rank][column]
      for inner in column..<columns {
        work[rank][inner] = work[rank][inner] / divisor
      }
      for row in 0..<rows {
        if row != rank {
          let factor = work[row][column]
          for inner in column..<columns {
            work[row][inner] = work[row][inner] - factor * work[rank][inner]
          }
        }
      }
      rank = rank + 1
    }
    column = column + 1
  }
  rank
}

///|
pub fn MuscatAdvancedDesign::create(
  matrix : Array[Array[Double]],
  sample_names : Array[String],
  coefficient_names : Array[String],
) -> MuscatAdvancedDesign raise MuscatAdvancedError {
  if matrix.length() == 0 || matrix.length() != sample_names.length() {
    raise MuscatAdvancedError(
      "design rows and sample names must be non-empty and have equal length",
    )
  }
  if matrix[0].length() == 0 || matrix[0].length() != coefficient_names.length() {
    raise MuscatAdvancedError(
      "design columns and coefficient names must be non-empty and have equal length",
    )
  }
  ignore(muscat_adv_names(sample_names, sample_names.length(), "Sample"))
  ignore(
    muscat_adv_names(
      coefficient_names,
      coefficient_names.length(),
      "Coefficient",
    ),
  )
  for row in matrix {
    if row.length() != coefficient_names.length() {
      raise MuscatAdvancedError("design matrix must be rectangular")
    }
    for value in row {
      if !muscat_adv_finite(value) {
        raise MuscatAdvancedError("design matrix must contain finite values")
      }
    }
  }
  if matrix.length() <= coefficient_names.length() {
    raise MuscatAdvancedError("design must have more samples than coefficients")
  }
  if muscat_adv_matrix_rank(matrix) != coefficient_names.length() {
    raise MuscatAdvancedError("design matrix must have full column rank")
  }
  MuscatAdvancedDesign::{
    matrix: muscat_adv_copy_matrix(matrix),
    sample_names: sample_names.copy(),
    coefficient_names: coefficient_names.copy(),
  }
}

///|
pub fn muscat_group_design_advanced(
  pseudobulk : MuscatAdvancedPseudoBulk,
  reference? : String = "",
) -> MuscatAdvancedDesign raise MuscatAdvancedError {
  let groups = muscat_adv_levels(
    pseudobulk.sample_groups,
    "pseudobulk sample groups",
  )
  let selected_reference = if reference == "" { groups[0] } else { reference }
  if !groups.contains(selected_reference) {
    raise MuscatAdvancedError("reference group is not present")
  }
  let encoded : Array[String] = []
  for group in groups {
    if group != selected_reference {
      encoded.push(group)
    }
  }
  let matrix : Array[Array[Double]] = []
  for group in pseudobulk.sample_groups {
    let row = Array::make(encoded.length() + 1, 0.0)
    row[0] = 1.0
    for index in 0..<encoded.length() {
      if group == encoded[index] {
        row[index + 1] = 1.0
      }
    }
    matrix.push(row)
  }
  let coefficient_names : Array[String] = ["(Intercept)"]
  for group in encoded {
    coefficient_names.push(group + "-" + selected_reference)
  }
  MuscatAdvancedDesign::create(
    matrix,
    pseudobulk.sample_names,
    coefficient_names,
  )
}

///|
pub fn muscat_contrast_advanced(
  design : MuscatAdvancedDesign,
  values : Array[Double],
  name : String,
) -> MuscatAdvancedContrast raise MuscatAdvancedError {
  if values.length() != design.coefficient_names.length() {
    raise MuscatAdvancedError(
      "contrast length must match the number of design coefficients",
    )
  }
  if name.trim().length() == 0 {
    raise MuscatAdvancedError("contrast name must not be empty")
  }
  let mut magnitude = 0.0
  for value in values {
    if !muscat_adv_finite(value) {
      raise MuscatAdvancedError("contrast must contain finite values")
    }
    magnitude = magnitude + value.abs()
  }
  if magnitude == 0.0 {
    raise MuscatAdvancedError("contrast must contain a non-zero coefficient")
  }
  MuscatAdvancedContrast::{ values: values.copy(), name, }
}

///|
pub fn muscat_default_contrasts_advanced(
  design : MuscatAdvancedDesign,
) -> Array[MuscatAdvancedContrast] raise MuscatAdvancedError {
  if design.coefficient_names.length() < 2 {
    raise MuscatAdvancedError(
      "default contrasts require at least one non-intercept coefficient",
    )
  }
  let output : Array[MuscatAdvancedContrast] = []
  for coefficient in 1..<design.coefficient_names.length() {
    let values = Array::make(design.coefficient_names.length(), 0.0)
    values[coefficient] = 1.0
    output.push(MuscatAdvancedContrast::{
      values,
      name: design.coefficient_names[coefficient],
    })
  }
  output
}

///|
fn muscat_adv_cholesky(matrix : Array[Array[Double]]) -> Array[Array[Double]]? {
  let size = matrix.length()
  if size == 0 {
    return None
  }
  let lower = muscat_adv_zero_matrix(size, size)
  for row in 0..<size {
    if matrix[row].length() != size {
      return None
    }
    for column in 0..<(row + 1) {
      let mut value = matrix[row][column]
      for inner in 0..<column {
        value = value - lower[row][inner] * lower[column][inner]
      }
      if row == column {
        if !muscat_adv_finite(value) || value <= 1.0e-14 {
          return None
        }
        lower[row][column] = value.sqrt()
      } else {
        lower[row][column] = value / lower[column][column]
      }
    }
  }
  Some(lower)
}

///|
fn muscat_adv_solve(
  lower : Array[Array[Double]],
  right_hand_side : Array[Double],
) -> Array[Double] {
  let size = lower.length()
  let forward = Array::make(size, 0.0)
  for row in 0..<size {
    let mut value = right_hand_side[row]
    for column in 0..<row {
      value = value - lower[row][column] * forward[column]
    }
    forward[row] = value / lower[row][row]
  }
  let output = Array::make(size, 0.0)
  let mut row = size - 1
  while row >= 0 {
    let mut value = forward[row]
    for column in (row + 1)..<size {
      value = value - lower[column][row] * output[column]
    }
    output[row] = value / lower[row][row]
    row = row - 1
  }
  output
}

///|
fn muscat_adv_inverse(lower : Array[Array[Double]]) -> Array[Array[Double]] {
  let size = lower.length()
  let output = muscat_adv_zero_matrix(size, size)
  for column in 0..<size {
    let basis = Array::make(size, 0.0)
    basis[column] = 1.0
    let solution = muscat_adv_solve(lower, basis)
    for row in 0..<size {
      output[row][column] = solution[row]
    }
  }
  output
}

///|
fn muscat_adv_regularized_cholesky(
  matrix : Array[Array[Double]],
  ridge : Double,
) -> Array[Array[Double]]? {
  let mut added = ridge
  for attempt in 0..<14 {
    let candidate = muscat_adv_copy_matrix(matrix)
    for index in 0..<candidate.length() {
      candidate[index][index] = candidate[index][index] + added
    }
    match muscat_adv_cholesky(candidate) {
      Some(lower) => return Some(lower)
      None =>
        added = if attempt == 0 && added == 0.0 {
          1.0e-10
        } else {
          added.max(1.0e-10) * 10.0
        }
    }
  }
  None
}

///|
fn muscat_adv_fit_gene(
  response : Array[Double],
  design : Array[Array[Double]],
  offsets : Array[Double],
  dispersion : Double,
  config : MuscatAdvancedConfig,
) -> MuscatAdvCoreFit {
  let observations = response.length()
  let parameters = design[0].length()
  let coefficients = Array::make(parameters, 0.0)
  let mut total_response = 0.0
  let mut total_exposure = 0.0
  for observation in 0..<observations {
    total_response = total_response + response[observation]
    total_exposure = total_exposure + @math.exp(offsets[observation])
  }
  if total_response == 0.0 {
    return MuscatAdvCoreFit::{
      coefficients,
      covariance: muscat_adv_zero_matrix(parameters, parameters),
      fitted_values: Array::make(observations, 0.0),
      converged: true,
      fitted: false,
    }
  }
  coefficients[0] = @math.ln((total_response + 0.5) / (total_exposure + 0.5))
  let fitted_values = Array::make(observations, 0.0)
  let mut converged = false
  for _ in 0..<config.max_iterations {
    let information = muscat_adv_zero_matrix(parameters, parameters)
    let right_hand_side = Array::make(parameters, 0.0)
    for observation in 0..<observations {
      let mut linear = offsets[observation]
      for parameter in 0..<parameters {
        linear = linear +
          design[observation][parameter] * coefficients[parameter]
      }
      let mean = @math.exp(muscat_adv_clamp(linear, -30.0, 30.0)).max(1.0e-12)
      fitted_values[observation] = mean
      let weight = mean / (1.0 + dispersion * mean)
      let working = linear +
        (response[observation] - mean) / mean -
        offsets[observation]
      for left in 0..<parameters {
        right_hand_side[left] = right_hand_side[left] +
          design[observation][left] * weight * working
        for right in 0..<parameters {
          information[left][right] = information[left][right] +
            design[observation][left] * weight * design[observation][right]
        }
      }
    }
    let lower = match
      muscat_adv_regularized_cholesky(information, config.ridge) {
      Some(value) => value
      None => break
    }
    let updated = muscat_adv_solve(lower, right_hand_side)
    let mut maximum_change = 0.0
    for parameter in 0..<parameters {
      maximum_change = maximum_change.max(
        (updated[parameter] - coefficients[parameter]).abs(),
      )
      coefficients[parameter] = muscat_adv_clamp(
        updated[parameter],
        -30.0,
        30.0,
      )
    }
    if maximum_change < config.tolerance {
      converged = true
      break
    }
  }
  let information = muscat_adv_zero_matrix(parameters, parameters)
  for observation in 0..<observations {
    let mut linear = offsets[observation]
    for parameter in 0..<parameters {
      linear = linear + design[observation][parameter] * coefficients[parameter]
    }
    let mean = @math.exp(muscat_adv_clamp(linear, -30.0, 30.0)).max(1.0e-12)
    fitted_values[observation] = mean
    let weight = mean / (1.0 + dispersion * mean)
    for left in 0..<parameters {
      for right in 0..<parameters {
        information[left][right] = information[left][right] +
          design[observation][left] * weight * design[observation][right]
      }
    }
  }
  let covariance = match
    muscat_adv_regularized_cholesky(information, config.ridge) {
    Some(lower) => muscat_adv_inverse(lower)
    None => muscat_adv_zero_matrix(parameters, parameters)
  }
  MuscatAdvCoreFit::{
    coefficients,
    covariance,
    fitted_values,
    converged,
    fitted: true,
  }
}

///|
fn muscat_adv_raw_dispersion(
  response : Array[Double],
  fitted : Array[Double],
  residual_df : Int,
  config : MuscatAdvancedConfig,
) -> Double {
  if residual_df <= 0 {
    return config.minimum_dispersion
  }
  let mut numerator = 0.0
  for observation in 0..<response.length() {
    let mean = fitted[observation].max(1.0e-12)
    let residual = response[observation] - mean
    numerator = numerator + (residual * residual - mean) / (mean * mean)
  }
  muscat_adv_clamp(
    numerator.max(0.0) / residual_df.to_double(),
    config.minimum_dispersion,
    config.maximum_dispersion,
  )
}

///|
fn muscat_adv_normal_cdf(value : Double) -> Double {
  let absolute = value.abs()
  let t = 1.0 / (1.0 + 0.2316419 * absolute)
  let density = 0.3989422804014327 * @math.exp(-0.5 * absolute * absolute)
  let polynomial = t *
    (
      0.319381530 +
      t *
      (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429)))
    )
  let upper = density * polynomial
  if value >= 0.0 {
    1.0 - upper
  } else {
    upper
  }
}

///|
fn muscat_adv_bh(values : Array[Double]) -> Array[Double] {
  if values.length() == 0 {
    return []
  }
  let order : Array[Int] = []
  for index in 0..<values.length() {
    order.push(index)
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    if values[left] < values[right] {
      -1
    } else if values[left] > values[right] {
      1
    } else {
      left - right
    }
  })
  let output = Array::make(values.length(), 1.0)
  let mut running = 1.0
  let mut position = values.length() - 1
  while position >= 0 {
    let adjusted = values[order[position]] *
      values.length().to_double() /
      (position + 1).to_double()
    running = running.min(adjusted).min(1.0)
    output[order[position]] = running
    position = position - 1
  }
  output
}

///|
fn muscat_adv_validate_fit_inputs(
  pseudobulk : MuscatAdvancedPseudoBulk,
  design : MuscatAdvancedDesign,
  contrasts : Array[MuscatAdvancedContrast],
) -> Unit raise MuscatAdvancedError {
  if pseudobulk.sample_names != design.sample_names {
    raise MuscatAdvancedError(
      "design sample names must exactly match pseudobulk sample order",
    )
  }
  if contrasts.length() == 0 {
    raise MuscatAdvancedError("at least one contrast is required")
  }
  let names : Array[String] = []
  for contrast in contrasts {
    if contrast.values.length() != design.coefficient_names.length() {
      raise MuscatAdvancedError(
        "all contrasts must match the design coefficient count",
      )
    }
    if contrast.name.trim().length() == 0 || names.contains(contrast.name) {
      raise MuscatAdvancedError("contrast names must be non-empty and unique")
    }
    names.push(contrast.name)
    let mut magnitude = 0.0
    for value in contrast.values {
      if !muscat_adv_finite(value) {
        raise MuscatAdvancedError("contrasts must contain finite values")
      }
      magnitude = magnitude + value.abs()
    }
    if magnitude == 0.0 {
      raise MuscatAdvancedError("contrasts must be non-zero")
    }
  }
}

///|
fn muscat_adv_offsets(
  pseudobulk : MuscatAdvancedPseudoBulk,
  cluster : Int,
  active : Array[Int],
  detection : Bool,
  config : MuscatAdvancedConfig,
) -> Array[Double] {
  let exposures = Array::make(active.length(), 1.0)
  if detection {
    let eligible : Array[Int] = []
    for gene in 0..<pseudobulk.gene_names.length() {
      let fractions : Array[Double] = []
      for sample in active {
        let cells = pseudobulk.cell_counts[cluster][sample]
        fractions.push(
          if cells == 0 {
            0.0
          } else {
            pseudobulk.values[cluster][gene][sample] / cells.to_double()
          },
        )
      }
      if muscat_adv_median(fractions) < config.detection_filter {
        eligible.push(gene)
      }
    }
    for position in 0..<active.length() {
      let sample = active[position]
      let cells = pseudobulk.cell_counts[cluster][sample].to_double()
      let mut fraction = 0.0
      if eligible.length() > 0 {
        for gene in eligible {
          fraction = fraction +
            pseudobulk.values[cluster][gene][sample] / cells.max(1.0)
        }
        fraction = fraction / eligible.length().to_double()
      } else {
        fraction = 1.0
      }
      exposures[position] = (cells * fraction.max(1.0e-6)).max(1.0e-8)
    }
  } else {
    for position in 0..<active.length() {
      exposures[position] = pseudobulk.library_sizes[cluster][active[position]].max(
        1.0,
      )
    }
    let mut log_sum = 0.0
    for exposure in exposures {
      log_sum = log_sum + @math.ln(exposure)
    }
    let geometric = @math.exp(log_sum / exposures.length().to_double())
    for position in 0..<exposures.length() {
      exposures[position] = exposures[position] / geometric.max(1.0e-12)
    }
  }
  exposures.map(fn(exposure) { @math.ln(exposure.max(1.0e-12)) })
}

///|
fn muscat_adv_log_cpm(
  response : Array[Double],
  exposures : Array[Double],
) -> Double {
  let mut normalized = 0.0
  for index in 0..<response.length() {
    normalized = normalized +
      response[index] / @math.exp(exposures[index]).max(1.0e-12)
  }
  @math.log2(normalized / response.length().to_double() + 0.5)
}

///|
fn muscat_adv_test(
  pseudobulk : MuscatAdvancedPseudoBulk,
  design : MuscatAdvancedDesign,
  contrasts : Array[MuscatAdvancedContrast],
  detection : Bool,
  config : MuscatAdvancedConfig,
) -> MuscatAdvancedResults raise MuscatAdvancedError {
  muscat_adv_validate_fit_inputs(pseudobulk, design, contrasts)
  if pseudobulk.scaled_cpm {
    raise MuscatAdvancedError(
      "differential testing requires unscaled pseudobulk values",
    )
  }
  if detection {
    if pseudobulk.aggregation != NumberDetected {
      raise MuscatAdvancedError(
        "differential detection requires number-detected pseudobulks",
      )
    }
  } else if pseudobulk.aggregation != SumCounts {
    raise MuscatAdvancedError(
      "differential state testing requires sum-count pseudobulks",
    )
  }
  let results : Array[MuscatAdvancedGeneResult] = []
  let mode = if detection { "DD" } else { "DS" }
  let parameters = design.coefficient_names.length()
  let mut any_tested = false
  for cluster in 0..<pseudobulk.cluster_names.length() {
    let active : Array[Int] = []
    for sample in 0..<pseudobulk.sample_names.length() {
      if pseudobulk.cell_counts[cluster][sample] >= config.min_cells &&
        pseudobulk.library_sizes[cluster][sample] > 0.0 {
        active.push(sample)
      }
    }
    let active_design : Array[Array[Double]] = []
    for sample in active {
      active_design.push(design.matrix[sample].copy())
    }
    let design_valid = active.length() > parameters &&
      muscat_adv_matrix_rank(active_design) == parameters
    let offsets = if design_valid {
      muscat_adv_offsets(pseudobulk, cluster, active, detection, config)
    } else {
      []
    }
    let first_fits : Array[MuscatAdvCoreFit] = []
    let raw_dispersions = Array::make(
      pseudobulk.gene_names.length(),
      config.minimum_dispersion,
    )
    let filtered = Array::make(pseudobulk.gene_names.length(), true)
    let trend_candidates : Array[Double] = []
    for gene in 0..<pseudobulk.gene_names.length() {
      let response : Array[Double] = []
      let detection_fractions : Array[Double] = []
      let mut passing = 0
      for sample in active {
        let value = pseudobulk.values[cluster][gene][sample]
        response.push(value)
        if detection {
          detection_fractions.push(
            value / pseudobulk.cell_counts[cluster][sample].to_double().max(1.0),
          )
        }
        if value >= config.min_count {
          passing = passing + 1
        }
      }
      let ubiquitous_detection = detection &&
        muscat_adv_median(detection_fractions) >= config.detection_filter
      filtered[gene] = !design_valid ||
        passing < config.min_samples ||
        ubiquitous_detection
      if filtered[gene] {
        first_fits.push(MuscatAdvCoreFit::{
          coefficients: Array::make(parameters, 0.0),
          covariance: muscat_adv_zero_matrix(parameters, parameters),
          fitted_values: Array::make(active.length(), 0.0),
          converged: false,
          fitted: false,
        })
      } else {
        let fit = muscat_adv_fit_gene(
          response,
          active_design,
          offsets,
          config.minimum_dispersion,
          config,
        )
        first_fits.push(fit)
        let raw = if fit.fitted {
          muscat_adv_raw_dispersion(
            response,
            fit.fitted_values,
            active.length() - parameters,
            config,
          )
        } else {
          config.minimum_dispersion
        }
        raw_dispersions[gene] = raw
        if fit.fitted {
          trend_candidates.push(raw)
        }
      }
    }
    let trend = if trend_candidates.length() == 0 {
      config.minimum_dispersion
    } else {
      muscat_adv_median(trend_candidates)
    }
    for gene in 0..<pseudobulk.gene_names.length() {
      let response : Array[Double] = []
      let mut mean_count = 0.0
      for sample in active {
        let value = pseudobulk.values[cluster][gene][sample]
        response.push(value)
        mean_count = mean_count + value
      }
      if active.length() > 0 {
        mean_count = mean_count / active.length().to_double()
      }
      let residual_df = (active.length() - parameters).max(0).to_double()
      let dispersion = if filtered[gene] {
        config.minimum_dispersion
      } else {
        muscat_adv_clamp(
          (
            residual_df * raw_dispersions[gene] +
            config.dispersion_prior_df * trend
          ) /
          (residual_df + config.dispersion_prior_df).max(1.0),
          config.minimum_dispersion,
          config.maximum_dispersion,
        )
      }
      let fit = if filtered[gene] {
        first_fits[gene]
      } else {
        muscat_adv_fit_gene(
          response, active_design, offsets, dispersion, config,
        )
      }
      for contrast in contrasts {
        let mut effect = 0.0
        let mut variance = 0.0
        if fit.fitted {
          for left in 0..<parameters {
            effect = effect + contrast.values[left] * fit.coefficients[left]
            for right in 0..<parameters {
              variance = variance +
                contrast.values[left] *
                fit.covariance[left][right] *
                contrast.values[right]
            }
          }
        }
        let threshold = config.lfc_threshold * @math.ln(2.0)
        let adjusted_effect = if effect > threshold {
          effect - threshold
        } else if effect < -threshold {
          effect + threshold
        } else {
          0.0
        }
        let tested = !filtered[gene] &&
          fit.fitted &&
          fit.converged &&
          muscat_adv_finite(variance) &&
          variance > 0.0
        let statistic = if tested {
          adjusted_effect * adjusted_effect / variance
        } else {
          0.0
        }
        let p_value = if tested {
          muscat_adv_clamp(
            2.0 * (1.0 - muscat_adv_normal_cdf(statistic.sqrt())),
            0.0,
            1.0,
          )
        } else {
          1.0
        }
        if tested {
          any_tested = true
        }
        results.push(MuscatAdvancedGeneResult::{
          gene_name: pseudobulk.gene_names[gene],
          cluster_name: pseudobulk.cluster_names[cluster],
          contrast_name: contrast.name,
          mode,
          tested,
          log_fc: if fit.fitted {
            effect / @math.ln(2.0)
          } else {
            0.0
          },
          log_cpm: if active.length() > 0 {
            muscat_adv_log_cpm(response, offsets)
          } else {
            0.0
          },
          statistic,
          p_value,
          local_fdr: 1.0,
          global_fdr: 1.0,
          raw_dispersion: raw_dispersions[gene],
          dispersion,
          mean_count,
          n_samples: active.length(),
          converged: fit.converged,
        })
      }
    }
  }
  if !any_tested {
    raise MuscatAdvancedError(
      "no genes could be tested; relax filtering or provide more replicated samples",
    )
  }
  for contrast in contrasts {
    for cluster in pseudobulk.cluster_names {
      let indices : Array[Int] = []
      let p_values : Array[Double] = []
      for index in 0..<results.length() {
        let result = results[index]
        if result.tested &&
          result.contrast_name == contrast.name &&
          result.cluster_name == cluster {
          indices.push(index)
          p_values.push(result.p_value)
        }
      }
      let adjusted = muscat_adv_bh(p_values)
      for position in 0..<indices.length() {
        results[indices[position]].local_fdr = adjusted[position]
      }
    }
    let indices : Array[Int] = []
    let p_values : Array[Double] = []
    for index in 0..<results.length() {
      let result = results[index]
      if result.tested && result.contrast_name == contrast.name {
        indices.push(index)
        p_values.push(result.p_value)
      }
    }
    let adjusted = muscat_adv_bh(p_values)
    for position in 0..<indices.length() {
      results[indices[position]].global_fdr = adjusted[position]
    }
  }
  MuscatAdvancedResults::{
    results,
    gene_names: pseudobulk.gene_names.copy(),
    cluster_names: pseudobulk.cluster_names.copy(),
    contrast_names: contrasts.map(fn(contrast) { contrast.name }),
    mode,
    fdr_threshold: config.fdr_threshold,
  }
}

///|
pub fn muscat_pbds_advanced(
  pseudobulk : MuscatAdvancedPseudoBulk,
  design : MuscatAdvancedDesign,
  contrasts : Array[MuscatAdvancedContrast],
  config? : MuscatAdvancedConfig = MuscatAdvancedConfig::default(),
) -> MuscatAdvancedResults raise MuscatAdvancedError {
  muscat_adv_test(pseudobulk, design, contrasts, false, config)
}

///|
pub fn muscat_pbdd_advanced(
  pseudobulk : MuscatAdvancedPseudoBulk,
  design : MuscatAdvancedDesign,
  contrasts : Array[MuscatAdvancedContrast],
  config? : MuscatAdvancedConfig = MuscatAdvancedConfig::default(),
) -> MuscatAdvancedResults raise MuscatAdvancedError {
  muscat_adv_test(pseudobulk, design, contrasts, true, config)
}

///|
pub fn MuscatAdvancedResults::get(
  self : MuscatAdvancedResults,
  gene_name : String,
  cluster_name : String,
  contrast_name : String,
) -> MuscatAdvancedGeneResult? {
  for result in self.results {
    if result.gene_name == gene_name &&
      result.cluster_name == cluster_name &&
      result.contrast_name == contrast_name {
      return Some(result)
    }
  }
  None
}

///|
pub fn MuscatAdvancedResults::n_tested(self : MuscatAdvancedResults) -> Int {
  let mut count = 0
  for result in self.results {
    if result.tested {
      count = count + 1
    }
  }
  count
}

///|
pub fn MuscatAdvancedResults::n_significant(
  self : MuscatAdvancedResults,
  global? : Bool = false,
) -> Int {
  let mut count = 0
  for result in self.results {
    let adjusted = if global { result.global_fdr } else { result.local_fdr }
    if result.tested && adjusted <= self.fdr_threshold {
      count = count + 1
    }
  }
  count
}

///|
pub fn MuscatAdvancedResults::summary(self : MuscatAdvancedResults) -> String {
  "muscat advanced " +
  self.mode +
  ": " +
  self.gene_names.length().to_string() +
  " genes x " +
  self.cluster_names.length().to_string() +
  " clusters x " +
  self.contrast_names.length().to_string() +
  " contrasts; " +
  self.n_tested().to_string() +
  " tested, " +
  self.n_significant().to_string() +
  " locally significant"
}

///|
pub fn muscat_stagewise_ds_dd(
  ds_result : MuscatAdvancedResults,
  dd_result : MuscatAdvancedResults,
  alpha? : Double = 0.05,
) -> MuscatStagewiseResults raise MuscatAdvancedError {
  if !muscat_adv_finite(alpha) || alpha <= 0.0 || alpha > 1.0 {
    raise MuscatAdvancedError("stagewise alpha must be in (0, 1]")
  }
  if ds_result.mode != "DS" || dd_result.mode != "DD" {
    raise MuscatAdvancedError("stagewise testing requires DS and DD results")
  }
  if ds_result.gene_names != dd_result.gene_names ||
    ds_result.cluster_names != dd_result.cluster_names ||
    ds_result.contrast_names != dd_result.contrast_names {
    raise MuscatAdvancedError(
      "DS and DD results must have matching genes, clusters and contrasts",
    )
  }
  let output : Array[MuscatStagewiseGeneResult] = []
  for contrast in ds_result.contrast_names {
    for cluster in ds_result.cluster_names {
      for gene in ds_result.gene_names {
        let ds = match ds_result.get(gene, cluster, contrast) {
          Some(value) => value
          None => raise MuscatAdvancedError("missing DS result")
        }
        let dd = match dd_result.get(gene, cluster, contrast) {
          Some(value) => value
          None => raise MuscatAdvancedError("missing DD result")
        }
        let ds_p = if ds.tested { ds.p_value } else { 1.0 }
        let dd_p = if dd.tested { dd.p_value } else { 1.0 }
        let screen = if ds_p == 0.0 || dd_p == 0.0 {
          0.0
        } else {
          2.0 / (1.0 / ds_p + 1.0 / dd_p)
        }
        output.push(MuscatStagewiseGeneResult::{
          gene_name: gene,
          cluster_name: cluster,
          contrast_name: contrast,
          screen_p_value: screen,
          screen_fdr: 1.0,
          ds_p_value: ds_p,
          dd_p_value: dd_p,
          ds_confirmation: 1.0,
          dd_confirmation: 1.0,
          classification: "none",
        })
      }
    }
  }
  for contrast in ds_result.contrast_names {
    let indices : Array[Int] = []
    let p_values : Array[Double] = []
    for index in 0..<output.length() {
      if output[index].contrast_name == contrast {
        indices.push(index)
        p_values.push(output[index].screen_p_value)
      }
    }
    let adjusted = muscat_adv_bh(p_values)
    let mut selected = 0
    for position in 0..<indices.length() {
      output[indices[position]].screen_fdr = adjusted[position]
      if adjusted[position] <= alpha {
        selected = selected + 1
      }
    }
    let selection_factor = if selected == 0 {
      indices.length().to_double()
    } else {
      indices.length().to_double() / selected.to_double()
    }
    for index in indices {
      if output[index].screen_fdr <= alpha {
        let ds_p = output[index].ds_p_value
        let dd_p = output[index].dd_p_value
        let ds_holm = if ds_p <= dd_p {
          (2.0 * ds_p).min(1.0)
        } else {
          dd_p.max(ds_p).min(1.0)
        }
        let dd_holm = if dd_p <= ds_p {
          (2.0 * dd_p).min(1.0)
        } else {
          ds_p.max(dd_p).min(1.0)
        }
        output[index].ds_confirmation = (ds_holm * selection_factor)
          .min(1.0)
          .max(output[index].screen_fdr)
        output[index].dd_confirmation = (dd_holm * selection_factor)
          .min(1.0)
          .max(output[index].screen_fdr)
        let ds_hit = output[index].ds_confirmation <= alpha
        let dd_hit = output[index].dd_confirmation <= alpha
        output[index].classification = if ds_hit && dd_hit {
          "both"
        } else if ds_hit {
          "DS"
        } else if dd_hit {
          "DD"
        } else {
          "screen_only"
        }
      }
    }
  }
  MuscatStagewiseResults::{
    results: output,
    gene_names: ds_result.gene_names.copy(),
    cluster_names: ds_result.cluster_names.copy(),
    contrast_names: ds_result.contrast_names.copy(),
    alpha,
  }
}

///|
pub fn MuscatStagewiseResults::get(
  self : MuscatStagewiseResults,
  gene_name : String,
  cluster_name : String,
  contrast_name : String,
) -> MuscatStagewiseGeneResult? {
  for result in self.results {
    if result.gene_name == gene_name &&
      result.cluster_name == cluster_name &&
      result.contrast_name == contrast_name {
      return Some(result)
    }
  }
  None
}

///|
pub fn MuscatStagewiseResults::n_confirmed(
  self : MuscatStagewiseResults,
) -> Int {
  let mut count = 0
  for result in self.results {
    if result.classification == "DS" ||
      result.classification == "DD" ||
      result.classification == "both" {
      count = count + 1
    }
  }
  count
}

///|
fn muscat_adv_copy_string_map(
  source : Map[String, Array[String]],
) -> Map[String, Array[String]] {
  let output : Map[String, Array[String]] = Map([])
  for key in source.keys() {
    output[key] = source[key].copy()
  }
  output
}

///|
fn muscat_adv_copy_sce(source : SingleCellExperiment) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  for key in source.assays.keys() {
    assays[key] = muscat_adv_copy_matrix(source.assays[key])
  }
  let reduced_dims : Map[String, Array[Array[Double]]] = Map([])
  for key in source.reduced_dims.keys() {
    reduced_dims[key] = muscat_adv_copy_matrix(source.reduced_dims[key])
  }
  let alternatives : Map[String, SingleCellExperiment] = Map([])
  for key in source.alternative_experiments.keys() {
    alternatives[key] = muscat_adv_copy_sce(source.alternative_experiments[key])
  }
  let metadata : Map[String, String] = Map([])
  for key in source.metadata.keys() {
    metadata[key] = source.metadata[key]
  }
  SingleCellExperiment::{
    assays,
    row_data: muscat_adv_copy_string_map(source.row_data),
    col_data: muscat_adv_copy_string_map(source.col_data),
    reduced_dims,
    alternative_experiments: alternatives,
    metadata,
    row_names: source.row_names.copy(),
    col_names: source.col_names.copy(),
  }
}

///|
fn muscat_adv_result_string(value : MuscatAdvancedGeneResult?) -> String {
  match value {
    Some(result) => result.global_fdr.to_string()
    None => "NA"
  }
}

///|
pub fn muscat_advanced_sce(
  experiment : SingleCellExperiment,
  sample_column : String,
  cluster_column : String,
  group_column : String,
  assay_name? : String = "counts",
  reference? : String = "",
  output_prefix? : String = "muscat",
  config? : MuscatAdvancedConfig = MuscatAdvancedConfig::default(),
) -> MuscatAdvancedSCEOutput raise MuscatAdvancedError {
  if output_prefix.trim().length() == 0 {
    raise MuscatAdvancedError("output_prefix must not be empty")
  }
  let counts = match experiment.assays.get(assay_name) {
    Some(value) => value
    None => raise MuscatAdvancedError("SCE assay not found: " + assay_name)
  }
  let sample_ids = match experiment.col_data.get(sample_column) {
    Some(value) => value
    None =>
      raise MuscatAdvancedError("SCE column data not found: " + sample_column)
  }
  let cluster_ids = match experiment.col_data.get(cluster_column) {
    Some(value) => value
    None =>
      raise MuscatAdvancedError("SCE column data not found: " + cluster_column)
  }
  let group_ids = match experiment.col_data.get(group_column) {
    Some(value) => value
    None =>
      raise MuscatAdvancedError("SCE column data not found: " + group_column)
  }
  let data = MuscatAdvancedData::create(
    counts,
    sample_ids,
    cluster_ids,
    group_ids,
    gene_names=experiment.row_names,
    cell_names=experiment.col_names,
  )
  let sum_pseudobulk = muscat_aggregate_advanced(data, aggregation=SumCounts)
  let detection_pseudobulk = muscat_aggregate_advanced(
    data,
    aggregation=NumberDetected,
  )
  let design = muscat_group_design_advanced(sum_pseudobulk, reference~)
  let contrasts = muscat_default_contrasts_advanced(design)
  let ds_result = muscat_pbds_advanced(
    sum_pseudobulk,
    design,
    contrasts,
    config~,
  )
  let dd_result = muscat_pbdd_advanced(
    detection_pseudobulk,
    design,
    contrasts,
    config~,
  )
  let stagewise_result = muscat_stagewise_ds_dd(
    ds_result,
    dd_result,
    alpha=config.fdr_threshold,
  )
  let enriched = muscat_adv_copy_sce(experiment)
  let contrast = contrasts[0].name
  for cluster in data.cluster_names {
    let ds_log_fc : Array[String] = []
    let ds_fdr : Array[String] = []
    let dd_log_fc : Array[String] = []
    let dd_fdr : Array[String] = []
    let stage_class : Array[String] = []
    for gene in data.gene_names {
      let ds = ds_result.get(gene, cluster, contrast)
      let dd = dd_result.get(gene, cluster, contrast)
      ds_log_fc.push(
        match ds {
          Some(result) => result.log_fc.to_string()
          None => "NA"
        },
      )
      ds_fdr.push(muscat_adv_result_string(ds))
      dd_log_fc.push(
        match dd {
          Some(result) => result.log_fc.to_string()
          None => "NA"
        },
      )
      dd_fdr.push(muscat_adv_result_string(dd))
      stage_class.push(
        match stagewise_result.get(gene, cluster, contrast) {
          Some(result) => result.classification
          None => "none"
        },
      )
    }
    let stem = output_prefix + "." + cluster
    enriched.row_data[stem + ".dsLogFC"] = ds_log_fc
    enriched.row_data[stem + ".dsFdr"] = ds_fdr
    enriched.row_data[stem + ".ddLogFC"] = dd_log_fc
    enriched.row_data[stem + ".ddFdr"] = dd_fdr
    enriched.row_data[stem + ".stageClass"] = stage_class
  }
  enriched.col_data[output_prefix + ".sample_id"] = sample_ids.copy()
  enriched.col_data[output_prefix + ".cluster_id"] = cluster_ids.copy()
  enriched.col_data[output_prefix + ".group_id"] = group_ids.copy()
  enriched.metadata[output_prefix + ".version"] = "muscat 1.27.4"
  enriched.metadata[output_prefix + ".assay"] = assay_name
  enriched.metadata[output_prefix + ".reference"] = if reference == "" {
    data.group_names[0]
  } else {
    reference
  }
  enriched.metadata[output_prefix + ".contrast"] = contrast
  enriched.metadata[output_prefix + ".clusters"] = data.cluster_names.join(",")
  enriched.metadata[output_prefix + ".dsTested"] = ds_result
    .n_tested()
    .to_string()
  enriched.metadata[output_prefix + ".ddTested"] = dd_result
    .n_tested()
    .to_string()
  enriched.metadata[output_prefix + ".stageConfirmed"] = stagewise_result
    .n_confirmed()
    .to_string()
  MuscatAdvancedSCEOutput::{
    experiment: enriched,
    sum_pseudobulk,
    detection_pseudobulk,
    ds_result,
    dd_result,
    stagewise_result,
  }
}

///|
pub fn muscat_advanced_example() -> MuscatAdvancedData raise MuscatAdvancedError {
  let genes = [
    "stable", "ds_a", "dd_b", "both_a", "sparse", "zero", "cdr1", "cdr2", "cdr3",
    "cdr4",
  ]
  let samples = ["C1", "C2", "C3", "C4", "T1", "T2", "T3", "T4"]
  let cells_per_stratum = 10
  let counts = muscat_adv_zero_matrix(
    genes.length(),
    samples.length() * 2 * cells_per_stratum,
  )
  let sample_ids : Array[String] = []
  let cluster_ids : Array[String] = []
  let group_ids : Array[String] = []
  let cell_names : Array[String] = []
  let mut cell = 0
  for sample_index in 0..<samples.length() {
    let treated = sample_index >= 4
    for cluster_index in 0..<2 {
      for replicate in 0..<cells_per_stratum {
        sample_ids.push(samples[sample_index])
        cluster_ids.push(if cluster_index == 0 { "A" } else { "B" })
        group_ids.push(if treated { "stim" } else { "ctrl" })
        cell_names.push("cell" + (cell + 1).to_string())
        counts[0][cell] = 100.0 + (replicate % 3).to_double()
        counts[1][cell] = if cluster_index == 0 && treated {
          18.0 + (replicate % 3).to_double()
        } else {
          3.0 + (replicate % 3).to_double()
        }
        counts[2][cell] = if cluster_index == 1 {
          if treated {
            10.0
          } else if replicate == 0 {
            100.0
          } else {
            0.0
          }
        } else {
          2.0
        }
        counts[3][cell] = if cluster_index == 0 {
          if treated {
            10.0
          } else if replicate == 0 {
            30.0
          } else {
            0.0
          }
        } else {
          2.0
        }
        counts[4][cell] = if replicate == 0 { 1.0 } else { 0.0 }
        for background in 0..<4 {
          counts[background + 6][cell] = if replicate == background + 1 {
            2.0
          } else {
            0.0
          }
        }
        cell = cell + 1
      }
    }
  }
  MuscatAdvancedData::create(
    counts,
    sample_ids,
    cluster_ids,
    group_ids,
    gene_names=genes,
    cell_names~,
  )
}
