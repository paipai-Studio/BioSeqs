///|
/// Cohort-scale pseudobulk differential expression inspired by Bioconductor
/// dreamlet. The mixed-model solver is provided by variance_partition.mbt;
/// this module owns aggregation, filtering, TMM normalization, voom-style
/// precision weights, per-cell-type orchestration, and study-wide FDR.

///|
pub suberror DreamletError {
  DreamletError(String)
}

///|
pub enum DreamletEffectKind {
  Numeric
  Categorical
  Random
} derive(Debug)

///|
pub struct DreamletEffectSpec {
  name : String
  kind : DreamletEffectKind
  reference : String
} derive(Debug)

///|
pub struct DreamletModelSpec {
  effects : Array[DreamletEffectSpec]
  intercept : Bool
} derive(Debug)

///|
pub struct DreamletPseudobulkAssay {
  cluster_id : String
  sample_names : Array[String]
  counts : Array[Array[Double]]
  cell_counts : Array[Int]
  library_sizes : Array[Double]
  metadata : Map[String, Array[String]]
} derive(Debug)

///|
pub struct DreamletPseudobulk {
  gene_names : Array[String]
  sample_names : Array[String]
  assays : Array[DreamletPseudobulkAssay]
} derive(Debug)

///|
pub struct DreamletProcessConfig {
  min_cells : Int
  min_count : Double
  min_samples : Int
  min_prop : Double
  min_total_count : Double
  prior_count : Double
  span : Double
  logratio_trim : Double
  sum_trim : Double
  use_poisson_initial_weights : Bool
  rescale_initial_weights : Bool
  max_iterations : Int
  tolerance : Double
} derive(Debug)

///|
pub struct DreamletTmmResult {
  factors : Array[Double]
  reference_sample : Int
} derive(Debug)

///|
pub struct DreamletExclusion {
  cluster_id : String
  reason : String
} derive(Debug)

///|
pub struct DreamletProcessedAssay {
  cluster_id : String
  sample_names : Array[String]
  gene_names : Array[String]
  counts : Array[Array[Double]]
  cell_counts : Array[Int]
  library_sizes : Array[Double]
  norm_factors : Array[Double]
  effective_library_sizes : Array[Double]
  normalized_cpm : Array[Array[Double]]
  log_cpm : Array[Array[Double]]
  initial_weights : Array[Array[Double]]
  precision_weights : Array[Array[Double]]
  design : VariancePartitionDesign
  metadata : Map[String, Array[String]]
  dropped_sample_names : Array[String]
  dropped_gene_names : Array[String]
  dropped_effects : Array[String]
  tmm_reference_sample : String
  voom_x : Array[Double]
  voom_y : Array[Double]
  voom_fitted : Array[Double]
} derive(Debug)

///|
pub struct DreamletProcessedData {
  assays : Array[DreamletProcessedAssay]
  exclusions : Array[DreamletExclusion]
  config : DreamletProcessConfig
} derive(Debug)

///|
pub struct DreamletAssayResult {
  cluster_id : String
  result : DreamResult
} derive(Debug)

///|
pub struct DreamletGeneResult {
  cluster_id : String
  gene_name : String
  estimate : Double
  standard_error : Double
  statistic : Double
  degrees_of_freedom : Double
  p_value : Double
  within_cluster_fdr : Double
  study_wide_fdr : Double
  converged : Bool
  random_variances : Array[Double]
  residual_variance : Double
} derive(Debug)

///|
pub struct DreamletResult {
  contrast_name : String
  coefficient_name : String
  assays : Array[DreamletAssayResult]
  genes : Array[DreamletGeneResult]
  exclusions : Array[DreamletExclusion]
} derive(Debug)

///|
fn dreamlet_fail(message : String) -> Unit raise DreamletError {
  raise DreamletError(message)
}

///|
fn dreamlet_is_finite(value : Double) -> Bool {
  value == value && value.abs() < 1.0e300
}

///|
fn dreamlet_log2(value : Double) -> Double {
  @math.ln(value) / @math.ln(2.0)
}

///|
fn dreamlet_copy_matrix(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in matrix {
    output.push(row.copy())
  }
  output
}

///|
fn dreamlet_copy_metadata(
  metadata : Map[String, Array[String]],
) -> Map[String, Array[String]] {
  let keys = metadata.keys().collect()
  let output : Map[String, Array[String]] = Map([], capacity=keys.length())
  for key in keys {
    output[key] = metadata[key].copy()
  }
  output
}

///|
fn dreamlet_mean(values : Array[Double]) -> Double {
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
fn dreamlet_median(values : Array[Double]) -> Double {
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
  let middle = sorted.length() / 2
  if sorted.length() % 2 == 1 {
    sorted[middle]
  } else {
    (sorted[middle - 1] + sorted[middle]) / 2.0
  }
}

///|
fn dreamlet_quantile(values : Array[Double], probability : Double) -> Double {
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
  let location = probability * (sorted.length() - 1).to_double()
  let lower = location.floor().to_int()
  let upper = (lower + 1).min(sorted.length() - 1)
  let fraction = location - lower.to_double()
  sorted[lower] + fraction * (sorted[upper] - sorted[lower])
}

///|
fn dreamlet_validate_names(
  values : Array[String],
  kind : String,
) -> Unit raise DreamletError {
  let seen : Map[String, Bool] = Map([], capacity=values.length())
  for value in values {
    if value == "" {
      dreamlet_fail("dreamlet " + kind + " names must not be empty")
    }
    if seen.contains(value) {
      dreamlet_fail("dreamlet " + kind + " names must be unique")
    }
    seen[value] = true
  }
}

///|
fn dreamlet_validate_count_matrix(
  counts : Array[Array[Double]],
  observations : Int,
  require_integer : Bool,
) -> Unit raise DreamletError {
  if counts.length() == 0 || observations == 0 {
    dreamlet_fail("dreamlet count matrix must not be empty")
  }
  for row in counts {
    if row.length() != observations {
      dreamlet_fail(
        "dreamlet count matrix must be rectangular and match the observations",
      )
    }
    for value in row {
      if !dreamlet_is_finite(value) || value < 0.0 {
        dreamlet_fail("dreamlet counts must be finite and non-negative")
      }
      if require_integer && (value - value.round()).abs() > 1.0e-8 {
        dreamlet_fail("dreamlet raw counts must contain integer values")
      }
    }
  }
}

///|
pub fn dreamlet_numeric_effect(name : String) -> DreamletEffectSpec {
  DreamletEffectSpec::{ name, kind: DreamletEffectKind::Numeric, reference: "" }
}

///|
pub fn dreamlet_categorical_effect(
  name : String,
  reference? : String = "",
) -> DreamletEffectSpec {
  DreamletEffectSpec::{ name, kind: DreamletEffectKind::Categorical, reference }
}

///|
pub fn dreamlet_random_effect(name : String) -> DreamletEffectSpec {
  DreamletEffectSpec::{ name, kind: DreamletEffectKind::Random, reference: "" }
}

///|
pub fn dreamlet_model(
  effects : Array[DreamletEffectSpec],
  intercept? : Bool = true,
) -> DreamletModelSpec raise DreamletError {
  let seen : Map[String, Bool] = Map([], capacity=effects.length())
  let mut fixed_effects = 0
  for effect in effects {
    if effect.name == "" {
      dreamlet_fail("dreamlet effect names must not be empty")
    }
    if seen.contains(effect.name) {
      dreamlet_fail("dreamlet effect names must be unique")
    }
    seen[effect.name] = true
    match effect.kind {
      DreamletEffectKind::Categorical => fixed_effects = fixed_effects + 1
      DreamletEffectKind::Numeric => {
        fixed_effects = fixed_effects + 1
        if effect.reference != "" {
          dreamlet_fail(
            "dreamlet reference levels are only valid for categorical effects",
          )
        }
      }
      DreamletEffectKind::Random =>
        if effect.reference != "" {
          dreamlet_fail(
            "dreamlet reference levels are only valid for categorical effects",
          )
        }
    }
  }
  if !intercept && fixed_effects == 0 {
    dreamlet_fail("dreamlet model must contain at least one fixed coefficient")
  }
  DreamletModelSpec::{ effects: effects.copy(), intercept }
}

///|
pub fn DreamletProcessConfig::create(
  min_cells? : Int = 5,
  min_count? : Double = 5.0,
  min_samples? : Int = 4,
  min_prop? : Double = 0.4,
  min_total_count? : Double = 15.0,
  prior_count? : Double = 0.5,
  span? : Double = 0.5,
  logratio_trim? : Double = 0.3,
  sum_trim? : Double = 0.05,
  use_poisson_initial_weights? : Bool = true,
  rescale_initial_weights? : Bool = false,
  max_iterations? : Int = 100,
  tolerance? : Double = 0.005,
) -> DreamletProcessConfig raise DreamletError {
  if min_cells < 1 {
    dreamlet_fail("dreamlet min_cells must be positive")
  }
  if !dreamlet_is_finite(min_count) || min_count < 0.0 {
    dreamlet_fail("dreamlet min_count must be finite and non-negative")
  }
  if min_samples < 3 {
    dreamlet_fail("dreamlet min_samples must be at least three")
  }
  if !dreamlet_is_finite(min_prop) || min_prop <= 0.0 || min_prop > 1.0 {
    dreamlet_fail("dreamlet min_prop must be in (0, 1]")
  }
  if !dreamlet_is_finite(min_total_count) || min_total_count < 0.0 {
    dreamlet_fail("dreamlet min_total_count must be finite and non-negative")
  }
  if !dreamlet_is_finite(prior_count) || prior_count <= 0.0 {
    dreamlet_fail("dreamlet prior_count must be finite and positive")
  }
  if !dreamlet_is_finite(span) || span <= 0.0 || span > 1.0 {
    dreamlet_fail("dreamlet LOWESS span must be in (0, 1]")
  }
  if !dreamlet_is_finite(logratio_trim) ||
    logratio_trim < 0.0 ||
    logratio_trim >= 0.5 {
    dreamlet_fail("dreamlet TMM log-ratio trim must be in [0, 0.5)")
  }
  if !dreamlet_is_finite(sum_trim) || sum_trim < 0.0 || sum_trim >= 0.5 {
    dreamlet_fail("dreamlet TMM abundance trim must be in [0, 0.5)")
  }
  if max_iterations <= 0 {
    dreamlet_fail("dreamlet max_iterations must be positive")
  }
  if !dreamlet_is_finite(tolerance) || tolerance <= 0.0 {
    dreamlet_fail("dreamlet tolerance must be finite and positive")
  }
  DreamletProcessConfig::{
    min_cells,
    min_count,
    min_samples,
    min_prop,
    min_total_count,
    prior_count,
    span,
    logratio_trim,
    sum_trim,
    use_poisson_initial_weights,
    rescale_initial_weights,
    max_iterations,
    tolerance,
  }
}

///|
pub fn DreamletProcessConfig::default() -> DreamletProcessConfig {
  DreamletProcessConfig::{
    min_cells: 5,
    min_count: 5.0,
    min_samples: 4,
    min_prop: 0.4,
    min_total_count: 15.0,
    prior_count: 0.5,
    span: 0.5,
    logratio_trim: 0.3,
    sum_trim: 0.05,
    use_poisson_initial_weights: true,
    rescale_initial_weights: false,
    max_iterations: 100,
    tolerance: 0.005,
  }
}

///|
pub fn dreamlet_aggregate_to_pseudobulk(
  counts : Array[Array[Double]],
  gene_names : Array[String],
  cell_names : Array[String],
  sample_ids : Array[String],
  cluster_ids : Array[String],
  cell_metadata? : Map[String, Array[String]] = Map([], capacity=0),
) -> DreamletPseudobulk raise DreamletError {
  let cells = cell_names.length()
  if gene_names.length() != counts.length() {
    dreamlet_fail("dreamlet gene names must match count matrix rows")
  }
  dreamlet_validate_names(gene_names, "gene")
  dreamlet_validate_names(cell_names, "cell")
  if sample_ids.length() != cells || cluster_ids.length() != cells {
    dreamlet_fail(
      "dreamlet sample and cluster identifiers must match count matrix columns",
    )
  }
  dreamlet_validate_count_matrix(counts, cells, true)
  let sample_names : Array[String] = []
  let cluster_names : Array[String] = []
  let sample_index : Map[String, Int] = Map([], capacity=cells)
  let cluster_index : Map[String, Int] = Map([], capacity=cells)
  for cell in 0..<cells {
    let sample = sample_ids[cell]
    let cluster = cluster_ids[cell]
    if sample == "" || cluster == "" {
      dreamlet_fail("dreamlet sample and cluster identifiers must not be empty")
    }
    if !sample_index.contains(sample) {
      sample_index[sample] = sample_names.length()
      sample_names.push(sample)
    }
    if !cluster_index.contains(cluster) {
      cluster_index[cluster] = cluster_names.length()
      cluster_names.push(cluster)
    }
  }
  let metadata_keys = cell_metadata.keys().collect()
  let sample_metadata : Map[String, Array[String]] = Map(
    [],
    capacity=metadata_keys.length(),
  )
  for key in metadata_keys {
    if key == "" {
      dreamlet_fail("dreamlet metadata names must not be empty")
    }
    let values = cell_metadata[key]
    if values.length() != cells {
      dreamlet_fail("dreamlet cell metadata must match count matrix columns")
    }
    let collapsed = Array::make(sample_names.length(), "")
    for cell in 0..<cells {
      if values[cell] == "" {
        dreamlet_fail("dreamlet metadata values must not be empty")
      }
      let sample = sample_index[sample_ids[cell]]
      if collapsed[sample] == "" {
        collapsed[sample] = values[cell]
      } else if collapsed[sample] != values[cell] {
        dreamlet_fail(
          "dreamlet metadata '" +
          key +
          "' is not constant within sample '" +
          sample_ids[cell] +
          "'",
        )
      }
    }
    sample_metadata[key] = collapsed
  }
  let assay_counts : Array[Array[Array[Double]]] = []
  let assay_cells : Array[Array[Int]] = []
  for _ in cluster_names {
    let matrix : Array[Array[Double]] = []
    for _ in gene_names {
      matrix.push(Array::make(sample_names.length(), 0.0))
    }
    assay_counts.push(matrix)
    assay_cells.push(Array::make(sample_names.length(), 0))
  }
  for cell in 0..<cells {
    let cluster = cluster_index[cluster_ids[cell]]
    let sample = sample_index[sample_ids[cell]]
    assay_cells[cluster][sample] = assay_cells[cluster][sample] + 1
    for gene in 0..<counts.length() {
      assay_counts[cluster][gene][sample] = assay_counts[cluster][gene][sample] +
        counts[gene][cell]
    }
  }
  let assays : Array[DreamletPseudobulkAssay] = []
  for cluster in 0..<cluster_names.length() {
    let library_sizes = Array::make(sample_names.length(), 0.0)
    for sample in 0..<sample_names.length() {
      let mut total = 0.0
      for gene in 0..<gene_names.length() {
        total = total + assay_counts[cluster][gene][sample]
      }
      library_sizes[sample] = total
    }
    assays.push(DreamletPseudobulkAssay::{
      cluster_id: cluster_names[cluster],
      sample_names: sample_names.copy(),
      counts: dreamlet_copy_matrix(assay_counts[cluster]),
      cell_counts: assay_cells[cluster].copy(),
      library_sizes,
      metadata: dreamlet_copy_metadata(sample_metadata),
    })
  }
  DreamletPseudobulk::{ gene_names: gene_names.copy(), sample_names, assays }
}

///|
pub fn dreamlet_aggregate_sce(
  sce : SingleCellExperiment,
  sample_id : String,
  cluster_id : String,
  assay? : String = "counts",
  metadata_fields? : Array[String] = [],
) -> DreamletPseudobulk raise DreamletError {
  if sample_id == "" || cluster_id == "" || assay == "" {
    dreamlet_fail("dreamlet SCE field and assay names must not be empty")
  }
  let counts = sce_get_assay(sce, assay)
  if counts.length() == 0 {
    dreamlet_fail("dreamlet SCE assay '" + assay + "' is absent or empty")
  }
  let samples = sce_get_col_data(sce, sample_id)
  let clusters = sce_get_col_data(sce, cluster_id)
  if samples.length() == 0 {
    dreamlet_fail("dreamlet SCE sample field '" + sample_id + "' is absent")
  }
  if clusters.length() == 0 {
    dreamlet_fail("dreamlet SCE cluster field '" + cluster_id + "' is absent")
  }
  let metadata : Map[String, Array[String]] = Map(
    [],
    capacity=metadata_fields.length(),
  )
  let seen : Map[String, Bool] = Map([], capacity=metadata_fields.length())
  for field in metadata_fields {
    if field == "" || seen.contains(field) {
      dreamlet_fail(
        "dreamlet SCE metadata field names must be non-empty and unique",
      )
    }
    seen[field] = true
    let values = sce_get_col_data(sce, field)
    if values.length() == 0 {
      dreamlet_fail("dreamlet SCE metadata field '" + field + "' is absent")
    }
    metadata[field] = values.copy()
  }
  dreamlet_aggregate_to_pseudobulk(
    counts,
    sce.row_names,
    sce.col_names,
    samples,
    clusters,
    cell_metadata=metadata,
  )
}

///|
fn dreamlet_tmm_rank(values : Array[Double], index : Int) -> Double {
  let mut less = 0
  let mut equal = 0
  for value in values {
    if value < values[index] {
      less = less + 1
    } else if value == values[index] {
      equal = equal + 1
    }
  }
  less.to_double() + (equal.to_double() + 1.0) / 2.0
}

///|
fn dreamlet_tmm_pair(
  counts : Array[Array[Double]],
  library_sizes : Array[Double],
  sample : Int,
  reference : Int,
  logratio_trim : Double,
  sum_trim : Double,
) -> Double {
  if sample == reference {
    return 1.0
  }
  let m_values : Array[Double] = []
  let a_values : Array[Double] = []
  let variances : Array[Double] = []
  let observed_library = library_sizes[sample]
  let reference_library = library_sizes[reference]
  for gene in 0..<counts.length() {
    let observed = counts[gene][sample]
    let reference_count = counts[gene][reference]
    if observed > 0.0 && reference_count > 0.0 {
      let observed_rate = observed / observed_library
      let reference_rate = reference_count / reference_library
      let m_value = dreamlet_log2(observed_rate / reference_rate)
      let a_value = 0.5 * dreamlet_log2(observed_rate * reference_rate)
      let variance = (observed_library - observed) /
        (observed_library * observed) +
        (reference_library - reference_count) /
        (reference_library * reference_count)
      if dreamlet_is_finite(m_value) &&
        dreamlet_is_finite(a_value) &&
        dreamlet_is_finite(variance) &&
        variance > 0.0 {
        m_values.push(m_value)
        a_values.push(a_value)
        variances.push(variance)
      }
    }
  }
  let count = m_values.length()
  if count == 0 {
    return 1.0
  }
  let lower_m = (count.to_double() * logratio_trim).floor().to_int() + 1
  let upper_m = count + 1 - lower_m
  let lower_a = (count.to_double() * sum_trim).floor().to_int() + 1
  let upper_a = count + 1 - lower_a
  let mut weighted_sum = 0.0
  let mut weight_sum = 0.0
  for index in 0..<count {
    let rank_m = dreamlet_tmm_rank(m_values, index)
    let rank_a = dreamlet_tmm_rank(a_values, index)
    if rank_m >= lower_m.to_double() &&
      rank_m <= upper_m.to_double() &&
      rank_a >= lower_a.to_double() &&
      rank_a <= upper_a.to_double() {
      let weight = 1.0 / variances[index]
      weighted_sum = weighted_sum + weight * m_values[index]
      weight_sum = weight_sum + weight
    }
  }
  if weight_sum <= 0.0 {
    1.0
  } else {
    @math.exp(weighted_sum / weight_sum * @math.ln(2.0))
  }
}

///|
pub fn dreamlet_tmm(
  counts : Array[Array[Double]],
  logratio_trim? : Double = 0.3,
  sum_trim? : Double = 0.05,
) -> DreamletTmmResult raise DreamletError {
  if counts.length() == 0 || counts[0].length() == 0 {
    dreamlet_fail("dreamlet TMM requires a non-empty count matrix")
  }
  let samples = counts[0].length()
  dreamlet_validate_count_matrix(counts, samples, false)
  if logratio_trim < 0.0 || logratio_trim >= 0.5 {
    dreamlet_fail("dreamlet TMM log-ratio trim must be in [0, 0.5)")
  }
  if sum_trim < 0.0 || sum_trim >= 0.5 {
    dreamlet_fail("dreamlet TMM abundance trim must be in [0, 0.5)")
  }
  let library_sizes = Array::make(samples, 0.0)
  for sample in 0..<samples {
    for gene in 0..<counts.length() {
      library_sizes[sample] = library_sizes[sample] + counts[gene][sample]
    }
    if library_sizes[sample] <= 0.0 {
      dreamlet_fail("dreamlet TMM requires positive library sizes")
    }
  }
  let upper_quartiles = Array::make(samples, 0.0)
  for sample in 0..<samples {
    let normalized : Array[Double] = []
    for gene in 0..<counts.length() {
      let mut all_zero = true
      for other in 0..<samples {
        if counts[gene][other] > 0.0 {
          all_zero = false
          break
        }
      }
      if !all_zero {
        normalized.push(counts[gene][sample] / library_sizes[sample])
      }
    }
    upper_quartiles[sample] = dreamlet_quantile(normalized, 0.75)
  }
  let mean_upper = dreamlet_mean(upper_quartiles)
  let median_upper = dreamlet_median(upper_quartiles)
  let mut reference = 0
  if median_upper > 1.0e-20 {
    let mut best = (upper_quartiles[0] - mean_upper).abs()
    for sample in 1..<samples {
      let distance = (upper_quartiles[sample] - mean_upper).abs()
      if distance < best {
        best = distance
        reference = sample
      }
    }
  } else {
    let mut best_mass = -1.0
    for sample in 0..<samples {
      let mut mass = 0.0
      for gene in 0..<counts.length() {
        mass = mass + counts[gene][sample].sqrt()
      }
      if mass > best_mass {
        best_mass = mass
        reference = sample
      }
    }
  }
  let factors = Array::make(samples, 1.0)
  for sample in 0..<samples {
    factors[sample] = dreamlet_tmm_pair(
      counts, library_sizes, sample, reference, logratio_trim, sum_trim,
    )
  }
  let mut mean_log = 0.0
  for factor in factors {
    mean_log = mean_log + @math.ln(factor)
  }
  mean_log = mean_log / samples.to_double()
  let scale = @math.exp(mean_log)
  for sample in 0..<samples {
    factors[sample] = factors[sample] / scale
  }
  DreamletTmmResult::{ factors, reference_sample: reference }
}

///|
pub fn dreamlet_compute_norm_counts(
  counts : Array[Array[Double]],
  effective_library_sizes : Array[Double],
) -> Array[Array[Double]] raise DreamletError {
  if effective_library_sizes.length() == 0 {
    dreamlet_fail("dreamlet normalized counts require library sizes")
  }
  dreamlet_validate_count_matrix(
    counts,
    effective_library_sizes.length(),
    false,
  )
  for size in effective_library_sizes {
    if !dreamlet_is_finite(size) || size <= 0.0 {
      dreamlet_fail(
        "dreamlet effective library sizes must be finite and positive",
      )
    }
  }
  let output : Array[Array[Double]] = []
  for row in counts {
    let normalized = Array::make(row.length(), 0.0)
    for sample in 0..<row.length() {
      normalized[sample] = row[sample] *
        1000000.0 /
        effective_library_sizes[sample]
    }
    output.push(normalized)
  }
  output
}

///|
pub fn dreamlet_compute_log_cpm(
  counts : Array[Array[Double]],
  effective_library_sizes : Array[Double],
  prior_count? : Double = 0.5,
) -> Array[Array[Double]] raise DreamletError {
  if !dreamlet_is_finite(prior_count) || prior_count <= 0.0 {
    dreamlet_fail("dreamlet log CPM prior count must be finite and positive")
  }
  if effective_library_sizes.length() == 0 {
    dreamlet_fail("dreamlet log CPM requires library sizes")
  }
  dreamlet_validate_count_matrix(
    counts,
    effective_library_sizes.length(),
    false,
  )
  for size in effective_library_sizes {
    if !dreamlet_is_finite(size) || size <= 0.0 {
      dreamlet_fail(
        "dreamlet effective library sizes must be finite and positive",
      )
    }
  }
  let output : Array[Array[Double]] = []
  for row in counts {
    let transformed = Array::make(row.length(), 0.0)
    for sample in 0..<row.length() {
      transformed[sample] = dreamlet_log2(
        (row[sample] + prior_count) *
        1000000.0 /
        (effective_library_sizes[sample] + 1.0),
      )
    }
    output.push(transformed)
  }
  output
}

///|
fn dreamlet_subset_metadata(
  metadata : Map[String, Array[String]],
  indices : Array[Int],
) -> Map[String, Array[String]] {
  let keys = metadata.keys().collect()
  let output : Map[String, Array[String]] = Map([], capacity=keys.length())
  for key in keys {
    let values : Array[String] = []
    for index in indices {
      values.push(metadata[key][index])
    }
    output[key] = values
  }
  output
}

///|
fn dreamlet_unique_count(values : Array[String]) -> Int {
  let seen : Map[String, Bool] = Map([], capacity=values.length())
  for value in values {
    seen[value] = true
  }
  seen.length()
}

///|
fn dreamlet_build_design(
  sample_names : Array[String],
  metadata : Map[String, Array[String]],
  model : DreamletModelSpec,
) -> (VariancePartitionDesign, Array[String]) raise DreamletError {
  let fixed : Array[VariancePartitionFixedEffect] = []
  let random : Array[VariancePartitionRandomEffect] = []
  let dropped : Array[String] = []
  for spec in model.effects {
    match metadata.get(spec.name) {
      None =>
        dreamlet_fail(
          "dreamlet model metadata field '" + spec.name + "' is absent",
        )
      Some(values) => {
        if values.length() != sample_names.length() {
          dreamlet_fail(
            "dreamlet model metadata does not match retained samples",
          )
        }
        match spec.kind {
          DreamletEffectKind::Numeric => {
            let numeric : Array[Double] = []
            for value in values {
              match parse_double(value) {
                Some(parsed) =>
                  if dreamlet_is_finite(parsed) {
                    numeric.push(parsed)
                  } else {
                    dreamlet_fail(
                      "dreamlet numeric metadata '" +
                      spec.name +
                      "' must be finite",
                    )
                  }
                None =>
                  dreamlet_fail(
                    "dreamlet numeric metadata '" +
                    spec.name +
                    "' contains an invalid value",
                  )
              }
            }
            let mut varying = false
            for index in 1..<numeric.length() {
              if (numeric[index] - numeric[0]).abs() > 1.0e-12 {
                varying = true
                break
              }
            }
            if varying {
              let effect = vp_numeric_effect(spec.name, numeric) catch {
                VariancePartitionError(message) =>
                  raise DreamletError("variancePartition: " + message)
              }
              fixed.push(effect)
            } else {
              dropped.push(spec.name + ": constant numeric effect")
            }
          }
          DreamletEffectKind::Categorical => {
            let levels = dreamlet_unique_count(values)
            if levels < 2 {
              dropped.push(spec.name + ": constant categorical effect")
            } else {
              if spec.reference != "" && !values.contains(spec.reference) {
                dreamlet_fail(
                  "dreamlet reference level '" +
                  spec.reference +
                  "' for '" +
                  spec.name +
                  "' is absent after filtering",
                )
              }
              let effect = vp_categorical_effect(
                spec.name,
                values,
                reference=spec.reference,
              ) catch {
                VariancePartitionError(message) =>
                  raise DreamletError("variancePartition: " + message)
              }
              fixed.push(effect)
            }
          }
          DreamletEffectKind::Random => {
            let levels = dreamlet_unique_count(values)
            if levels < 2 {
              dropped.push(spec.name + ": constant random effect")
            } else if levels >= values.length() {
              dropped.push(
                spec.name + ": random effect has no replicated levels",
              )
            } else {
              let effect = vp_random_effect(spec.name, values) catch {
                VariancePartitionError(message) =>
                  raise DreamletError("variancePartition: " + message)
              }
              random.push(effect)
            }
          }
        }
      }
    }
  }
  let design = vp_design(sample_names, fixed, random, intercept=model.intercept) catch {
    VariancePartitionError(message) =>
      raise DreamletError("variancePartition: " + message)
  }
  (design, dropped)
}

///|
fn dreamlet_initial_weights(
  counts : Array[Array[Double]],
  prior_count : Double,
  poisson : Bool,
) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in counts {
    let weights = Array::make(row.length(), 1.0)
    if poisson {
      let mut total = 0.0
      for sample in 0..<row.length() {
        weights[sample] = row[sample] + prior_count
        total = total + weights[sample]
      }
      let mean = total / row.length().to_double()
      for sample in 0..<row.length() {
        weights[sample] = weights[sample] / mean
      }
    }
    output.push(weights)
  }
  output
}

///|
fn dreamlet_process_one(
  source : DreamletPseudobulkAssay,
  gene_names : Array[String],
  model : DreamletModelSpec,
  config : DreamletProcessConfig,
) -> DreamletProcessedAssay raise DreamletError {
  let retained_samples : Array[Int] = []
  let dropped_samples : Array[String] = []
  for sample in 0..<source.sample_names.length() {
    if source.cell_counts[sample] >= config.min_cells &&
      source.library_sizes[sample] > 0.0 {
      retained_samples.push(sample)
    } else {
      dropped_samples.push(source.sample_names[sample])
    }
  }
  if retained_samples.length() < config.min_samples {
    dreamlet_fail(
      "fewer than " +
      config.min_samples.to_string() +
      " samples pass cell-count and library-size filters",
    )
  }
  let sample_names : Array[String] = []
  let cell_counts : Array[Int] = []
  let library_sizes : Array[Double] = []
  for index in retained_samples {
    sample_names.push(source.sample_names[index])
    cell_counts.push(source.cell_counts[index])
    library_sizes.push(source.library_sizes[index])
  }
  let retained_counts : Array[Array[Double]] = []
  for row in source.counts {
    let selected : Array[Double] = []
    for index in retained_samples {
      selected.push(row[index])
    }
    retained_counts.push(selected)
  }
  let metadata = dreamlet_subset_metadata(source.metadata, retained_samples)
  let tmm = dreamlet_tmm(
    retained_counts,
    logratio_trim=config.logratio_trim,
    sum_trim=config.sum_trim,
  )
  let effective_sizes = Array::make(sample_names.length(), 0.0)
  for sample in 0..<sample_names.length() {
    effective_sizes[sample] = library_sizes[sample] * tmm.factors[sample]
  }
  let median_library = dreamlet_median(library_sizes)
  if median_library <= 0.0 {
    dreamlet_fail("retained pseudobulk libraries have zero median size")
  }
  let cpm_cutoff = config.min_count / median_library * 1000000.0
  let required_samples = (config.min_prop * sample_names.length().to_double())
    .ceil()
    .to_int()
    .max(1)
  let retained_genes : Array[Int] = []
  let dropped_genes : Array[String] = []
  for gene in 0..<retained_counts.length() {
    let mut total = 0.0
    let mut expressed = 0
    for sample in 0..<sample_names.length() {
      total = total + retained_counts[gene][sample]
      let cpm = retained_counts[gene][sample] *
        1000000.0 /
        effective_sizes[sample]
      if cpm >= cpm_cutoff {
        expressed = expressed + 1
      }
    }
    if total >= config.min_total_count && expressed >= required_samples {
      retained_genes.push(gene)
    } else {
      dropped_genes.push(gene_names[gene])
    }
  }
  if retained_genes.length() < 2 {
    dreamlet_fail("fewer than two genes pass expression filters")
  }
  let filtered_names : Array[String] = []
  let filtered_counts : Array[Array[Double]] = []
  for gene in retained_genes {
    filtered_names.push(gene_names[gene])
    filtered_counts.push(retained_counts[gene].copy())
  }
  let normalized = dreamlet_compute_norm_counts(
    filtered_counts, effective_sizes,
  )
  let log_cpm = dreamlet_compute_log_cpm(
    filtered_counts,
    effective_sizes,
    prior_count=config.prior_count,
  )
  let (design, dropped_effects) = dreamlet_build_design(
    sample_names, metadata, model,
  )
  let initial = dreamlet_initial_weights(
    filtered_counts,
    config.prior_count,
    config.use_poisson_initial_weights,
  )
  let fits : Array[VariancePartitionFit] = []
  let trend_points : Array[(Double, Double)] = []
  let mean_log_library : Double = {
    let values : Array[Double] = []
    for size in effective_sizes {
      values.push(dreamlet_log2(size + 1.0))
    }
    dreamlet_mean(values)
  }
  for gene in 0..<filtered_names.length() {
    let fit = vp_fit_gene(
      log_cpm[gene],
      gene_name=filtered_names[gene],
      design,
      weights=initial[gene],
      reml=true,
      max_iterations=config.max_iterations,
      tolerance=config.tolerance,
    ) catch {
      VariancePartitionError(message) =>
        raise DreamletError(
          "initial precision model failed for gene '" +
          filtered_names[gene] +
          "': " +
          message,
        )
    }
    fits.push(fit)
    let x = dreamlet_mean(log_cpm[gene]) +
      mean_log_library -
      dreamlet_log2(1000000.0)
    let y = fit.residual_variance.max(0.0).sqrt().sqrt()
    trend_points.push((x, y))
  }
  trend_points.sort_by(fn(
    left : (Double, Double),
    right : (Double, Double),
  ) -> Int {
    if left.0 < right.0 {
      -1
    } else if left.0 > right.0 {
      1
    } else {
      0
    }
  })
  let trend_x : Array[Double] = []
  let trend_y : Array[Double] = []
  for point in trend_points {
    trend_x.push(point.0)
    trend_y.push(point.1)
  }
  let trend_fitted = lowess(
    trend_x,
    trend_y,
    f=config.span,
    n_steps=0,
    delta=0.0,
  )
  let precision : Array[Array[Double]] = []
  for gene in 0..<filtered_names.length() {
    let weights = Array::make(sample_names.length(), 1.0)
    for sample in 0..<sample_names.length() {
      let fitted_cpm = @math.exp(
        fits[gene].fitted_values[sample] * @math.ln(2.0),
      )
      let fitted_count = 1.0e-6 * fitted_cpm * (effective_sizes[sample] + 1.0)
      let fitted_log_count = dreamlet_log2(fitted_count.max(1.0e-12))
      let predicted = lowess_predict(
        fitted_log_count,
        trend_x,
        trend_y,
        trend_fitted,
        f=config.span,
      ).max(0.01)
      let predicted_variance = predicted * predicted * predicted * predicted
      let mut weight = 1.0 / predicted_variance
      if config.rescale_initial_weights {
        weight = weight * initial[gene][sample]
      }
      weights[sample] = weight.min(1.0e8).max(1.0e-8)
    }
    precision.push(weights)
  }
  DreamletProcessedAssay::{
    cluster_id: source.cluster_id,
    sample_names,
    gene_names: filtered_names,
    counts: filtered_counts,
    cell_counts,
    library_sizes,
    norm_factors: tmm.factors,
    effective_library_sizes: effective_sizes,
    normalized_cpm: normalized,
    log_cpm,
    initial_weights: initial,
    precision_weights: precision,
    design,
    metadata,
    dropped_sample_names: dropped_samples,
    dropped_gene_names: dropped_genes,
    dropped_effects,
    tmm_reference_sample: source.sample_names[retained_samples[tmm.reference_sample]],
    voom_x: trend_x,
    voom_y: trend_y,
    voom_fitted: trend_fitted,
  }
}

///|
pub fn dreamlet_process_assays(
  pseudobulk : DreamletPseudobulk,
  model : DreamletModelSpec,
  config? : DreamletProcessConfig = DreamletProcessConfig::default(),
) -> DreamletProcessedData raise DreamletError {
  if pseudobulk.assays.length() == 0 {
    dreamlet_fail("dreamlet pseudobulk object contains no assays")
  }
  if pseudobulk.gene_names.length() == 0 {
    dreamlet_fail("dreamlet pseudobulk object contains no genes")
  }
  ignore(
    DreamletProcessConfig::create(
      min_cells=config.min_cells,
      min_count=config.min_count,
      min_samples=config.min_samples,
      min_prop=config.min_prop,
      min_total_count=config.min_total_count,
      prior_count=config.prior_count,
      span=config.span,
      logratio_trim=config.logratio_trim,
      sum_trim=config.sum_trim,
      use_poisson_initial_weights=config.use_poisson_initial_weights,
      rescale_initial_weights=config.rescale_initial_weights,
      max_iterations=config.max_iterations,
      tolerance=config.tolerance,
    ),
  )
  let assays : Array[DreamletProcessedAssay] = []
  let exclusions : Array[DreamletExclusion] = []
  for source in pseudobulk.assays {
    assays.push(
      dreamlet_process_one(source, pseudobulk.gene_names, model, config),
    ) catch {
      DreamletError(message) =>
        exclusions.push(DreamletExclusion::{
          cluster_id: source.cluster_id,
          reason: message,
        })
    }
  }
  if assays.length() == 0 {
    let reason = if exclusions.length() > 0 {
      ": " + exclusions[0].reason
    } else {
      ""
    }
    dreamlet_fail("dreamlet retained no assays" + reason)
  }
  DreamletProcessedData::{ assays, exclusions, config }
}

///|
fn dreamlet_bh_adjust(p_values : Array[Double]) -> Array[Double] {
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
    } else {
      left - right
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
pub fn dreamlet(
  processed : DreamletProcessedData,
  coefficient_name : String,
  contrast_name? : String = "",
  reml? : Bool = true,
  max_iterations? : Int = 100,
  tolerance? : Double = 0.005,
) -> DreamletResult raise DreamletError {
  if coefficient_name == "" {
    dreamlet_fail("dreamlet coefficient name must not be empty")
  }
  if max_iterations <= 0 || !dreamlet_is_finite(tolerance) || tolerance <= 0.0 {
    dreamlet_fail("dreamlet fit controls are invalid")
  }
  let resolved_name = if contrast_name == "" {
    coefficient_name
  } else {
    contrast_name
  }
  let assay_results : Array[DreamletAssayResult] = []
  let exclusions = processed.exclusions.copy()
  let flattened : Array[DreamletGeneResult] = []
  for assay in processed.assays {
    try {
      let contrast = vp_coefficient_contrast(
        assay.design,
        resolved_name,
        coefficient_name,
      )
      let result = dream(
        assay.log_cpm,
        assay.design,
        contrast,
        gene_names=assay.gene_names,
        weights=assay.precision_weights,
        reml~,
        max_iterations~,
        tolerance~,
      )
      assay_results.push(DreamletAssayResult::{
        cluster_id: assay.cluster_id,
        result,
      })
      for gene in result.genes {
        flattened.push(DreamletGeneResult::{
          cluster_id: assay.cluster_id,
          gene_name: gene.gene_name,
          estimate: gene.estimate,
          standard_error: gene.standard_error,
          statistic: gene.statistic,
          degrees_of_freedom: gene.degrees_of_freedom,
          p_value: gene.p_value,
          within_cluster_fdr: gene.adjusted_p_value,
          study_wide_fdr: 1.0,
          converged: gene.converged,
          random_variances: gene.random_variances.copy(),
          residual_variance: gene.residual_variance,
        })
      }
    } catch {
      VariancePartitionError(message) =>
        exclusions.push(DreamletExclusion::{
          cluster_id: assay.cluster_id,
          reason: "variancePartition: " + message,
        })
    }
  }
  if assay_results.length() == 0 {
    dreamlet_fail(
      "dreamlet coefficient '" +
      coefficient_name +
      "' could not be fitted in any assay",
    )
  }
  let p_values : Array[Double] = []
  for gene in flattened {
    p_values.push(gene.p_value)
  }
  let adjusted = dreamlet_bh_adjust(p_values)
  let genes : Array[DreamletGeneResult] = []
  for index in 0..<flattened.length() {
    let gene = flattened[index]
    genes.push(DreamletGeneResult::{
      cluster_id: gene.cluster_id,
      gene_name: gene.gene_name,
      estimate: gene.estimate,
      standard_error: gene.standard_error,
      statistic: gene.statistic,
      degrees_of_freedom: gene.degrees_of_freedom,
      p_value: gene.p_value,
      within_cluster_fdr: gene.within_cluster_fdr,
      study_wide_fdr: adjusted[index],
      converged: gene.converged,
      random_variances: gene.random_variances.copy(),
      residual_variance: gene.residual_variance,
    })
  }
  DreamletResult::{
    contrast_name: resolved_name,
    coefficient_name,
    assays: assay_results,
    genes,
    exclusions,
  }
}

///|
pub fn DreamletPseudobulk::assay(
  self : DreamletPseudobulk,
  cluster_id : String,
) -> DreamletPseudobulkAssay? {
  for assay in self.assays {
    if assay.cluster_id == cluster_id {
      return Some(assay)
    }
  }
  None
}

///|
pub fn DreamletProcessedData::assay(
  self : DreamletProcessedData,
  cluster_id : String,
) -> DreamletProcessedAssay? {
  for assay in self.assays {
    if assay.cluster_id == cluster_id {
      return Some(assay)
    }
  }
  None
}

///|
pub fn DreamletResult::assay(
  self : DreamletResult,
  cluster_id : String,
) -> DreamResult? {
  for assay in self.assays {
    if assay.cluster_id == cluster_id {
      return Some(assay.result)
    }
  }
  None
}

///|
pub fn DreamletResult::top_table(
  self : DreamletResult,
  maximum_fdr? : Double = 1.0,
  cluster_id? : String = "",
) -> Array[DreamletGeneResult] raise DreamletError {
  if !dreamlet_is_finite(maximum_fdr) || maximum_fdr < 0.0 || maximum_fdr > 1.0 {
    dreamlet_fail("dreamlet maximum FDR must be between zero and one")
  }
  let selected : Array[DreamletGeneResult] = []
  for gene in self.genes {
    if gene.study_wide_fdr <= maximum_fdr &&
      (cluster_id == "" || gene.cluster_id == cluster_id) {
      selected.push(gene)
    }
  }
  selected.sort_by(fn(
    left : DreamletGeneResult,
    right : DreamletGeneResult,
  ) -> Int {
    if left.p_value < right.p_value {
      -1
    } else if left.p_value > right.p_value {
      1
    } else if left.cluster_id < right.cluster_id {
      -1
    } else if left.cluster_id > right.cluster_id {
      1
    } else if left.gene_name < right.gene_name {
      -1
    } else if left.gene_name > right.gene_name {
      1
    } else {
      0
    }
  })
  selected
}

///|
pub fn DreamletPseudobulk::summary(self : DreamletPseudobulk) -> String {
  let mut combinations = 0
  for assay in self.assays {
    for cells in assay.cell_counts {
      if cells > 0 {
        combinations = combinations + 1
      }
    }
  }
  "dreamlet Pseudobulk Summary\n" +
  "Genes: " +
  self.gene_names.length().to_string() +
  "\nSamples: " +
  self.sample_names.length().to_string() +
  "\nCell types: " +
  self.assays.length().to_string() +
  "\nObserved sample-cell type combinations: " +
  combinations.to_string() +
  "\n"
}

///|
pub fn DreamletProcessedData::summary(self : DreamletProcessedData) -> String {
  let mut genes = 0
  for assay in self.assays {
    genes = genes + assay.gene_names.length()
  }
  "dreamlet Processed Data Summary\n" +
  "Retained cell types: " +
  self.assays.length().to_string() +
  "\nExcluded cell types: " +
  self.exclusions.length().to_string() +
  "\nGene-cell type models: " +
  genes.to_string() +
  "\n"
}

///|
pub fn DreamletResult::summary(self : DreamletResult) -> String {
  let mut significant = 0
  let mut converged = 0
  for gene in self.genes {
    if gene.study_wide_fdr <= 0.05 {
      significant = significant + 1
    }
    if gene.converged {
      converged = converged + 1
    }
  }
  "dreamlet Differential Expression Summary\n" +
  "Contrast: " +
  self.contrast_name +
  "\nCell types: " +
  self.assays.length().to_string() +
  "\nTests: " +
  self.genes.length().to_string() +
  "\nConverged: " +
  converged.to_string() +
  "\nStudy-wide FDR <= 0.05: " +
  significant.to_string() +
  "\n"
}
