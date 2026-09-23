// Ambient RNA decontamination for single-cell count matrices.
//
// This module follows the model introduced by Bioconductor decontX:
// each cell is a mixture of its cluster-native multinomial distribution
// and a contamination distribution derived from the other populations.

///|
pub suberror DecontXError {
  DecontXError(String)
}

///|
/// Parameters for deterministic expectation-maximization inference.
pub struct DecontXConfig {
  max_iterations : Int
  tolerance : Double
  contamination_prior_alpha : Double
  contamination_prior_beta : Double
  profile_pseudocount : Double
  initial_contamination : Double
  minimum_probability : Double
} derive(Eq, Debug)
pub extend DecontXConfig with Eq::{equal, not_equal}
pub extend DecontXConfig with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn DecontXConfig::create(
  max_iterations? : Int = 100,
  tolerance? : Double = 1.0e-6,
  contamination_prior_alpha? : Double = 1.0,
  contamination_prior_beta? : Double = 9.0,
  profile_pseudocount? : Double = 0.5,
  initial_contamination? : Double = 0.1,
  minimum_probability? : Double = 1.0e-12,
) -> DecontXConfig raise DecontXError {
  if max_iterations < 1 {
    raise DecontXError::DecontXError("max_iterations must be positive")
  }
  if tolerance <= 0.0 {
    raise DecontXError::DecontXError("tolerance must be positive")
  }
  if contamination_prior_alpha <= 0.0 || contamination_prior_beta <= 0.0 {
    raise DecontXError::DecontXError("Beta prior parameters must be positive")
  }
  if profile_pseudocount <= 0.0 {
    raise DecontXError::DecontXError("profile_pseudocount must be positive")
  }
  if initial_contamination < 0.0 || initial_contamination >= 1.0 {
    raise DecontXError::DecontXError("initial_contamination must be in [0, 1)")
  }
  if minimum_probability <= 0.0 || minimum_probability >= 1.0 {
    raise DecontXError::DecontXError("minimum_probability must be in (0, 1)")
  }
  DecontXConfig::{
    max_iterations,
    tolerance,
    contamination_prior_alpha,
    contamination_prior_beta,
    profile_pseudocount,
    initial_contamination,
    minimum_probability,
  }
}

///|
pub fn DecontXConfig::default() -> DecontXConfig {
  DecontXConfig::{
    max_iterations: 100,
    tolerance: 1.0e-6,
    contamination_prior_alpha: 1.0,
    contamination_prior_beta: 9.0,
    profile_pseudocount: 0.5,
    initial_contamination: 0.1,
    minimum_probability: 1.0e-12,
  }
}

///|
pub struct DecontXDiagnostics {
  iterations : Int
  converged : Bool
  initial_log_likelihood : Double
  final_log_likelihood : Double
  log_likelihoods : Array[Double]
} derive(Eq, Debug)
pub extend DecontXDiagnostics with Eq::{equal, not_equal}
pub extend DecontXDiagnostics with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct DecontXCellEstimate {
  cell_name : String
  cluster : String
  contamination : Double
  library_size : Double
  native_counts : Double
  contaminant_counts : Double
} derive(Eq, Debug)
pub extend DecontXCellEstimate with Eq::{equal, not_equal}
pub extend DecontXCellEstimate with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct DecontXClusterEstimate {
  cluster : String
  n_cells : Int
  mean_contamination : Double
  total_counts : Double
  native_counts : Double
  contaminant_counts : Double
} derive(Eq, Debug)
pub extend DecontXClusterEstimate with Eq::{equal, not_equal}
pub extend DecontXClusterEstimate with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct DecontXResult {
  corrected_counts : Array[Array[Double]]
  contamination_counts : Array[Array[Double]]
  contamination : Array[Double]
  native_profiles : Array[Array[Double]]
  contaminant_profiles : Array[Array[Double]]
  cluster_names : Array[String]
  cluster_indices : Array[Int]
  gene_names : Array[String]
  cell_names : Array[String]
  used_background : Bool
  diagnostics : DecontXDiagnostics
} derive(Eq, Debug)
pub extend DecontXResult with Eq::{equal, not_equal}
pub extend DecontXResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Result of adding decontamination outputs to a SingleCellExperiment copy.
pub struct DecontXSCEOutput {
  experiment : SingleCellExperiment
  result : DecontXResult
}

///|
fn decontx_zero_matrix(rows : Int, columns : Int) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  let mut row = 0
  while row < rows {
    result.push(Array::make(columns, 0.0))
    row = row + 1
  }
  result
}

///|
fn decontx_copy_matrix(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  let mut row = 0
  while row < matrix.length() {
    result.push(matrix[row].copy())
    row = row + 1
  }
  result
}

///|
fn decontx_validate_matrix(
  matrix : Array[Array[Double]],
  label : String,
  expected_rows? : Int = -1,
) -> (Int, Int) raise DecontXError {
  if matrix.length() == 0 {
    raise DecontXError::DecontXError(label + " must contain at least one gene")
  }
  if expected_rows >= 0 && matrix.length() != expected_rows {
    raise DecontXError::DecontXError(
      label + " must have the same genes as counts",
    )
  }
  let columns = matrix[0].length()
  if columns == 0 {
    raise DecontXError::DecontXError(label + " must contain at least one cell")
  }
  let mut row = 0
  while row < matrix.length() {
    if matrix[row].length() != columns {
      raise DecontXError::DecontXError(label + " must be rectangular")
    }
    let mut column = 0
    while column < columns {
      let value = matrix[row][column]
      if value != value || value < 0.0 || value.abs() > 1.0e300 {
        raise DecontXError::DecontXError(
          label + " must contain finite non-negative counts",
        )
      }
      column = column + 1
    }
    row = row + 1
  }
  (matrix.length(), columns)
}

///|
fn decontx_prepare_names(
  supplied : Array[String],
  length : Int,
  prefix : String,
) -> Array[String] raise DecontXError {
  if supplied.length() != 0 && supplied.length() != length {
    raise DecontXError::DecontXError(
      prefix + " names must match the matrix dimension",
    )
  }
  let names : Array[String] = Array::new()
  let mut index = 0
  while index < length {
    let name = if supplied.length() == 0 {
      prefix + (index + 1).to_string()
    } else {
      supplied[index]
    }
    if name.trim().length() == 0 {
      raise DecontXError::DecontXError(prefix + " names cannot be empty")
    }
    let mut previous = 0
    while previous < names.length() {
      if names[previous] == name {
        raise DecontXError::DecontXError(
          prefix + " names must be unique: " + name,
        )
      }
      previous = previous + 1
    }
    names.push(name)
    index = index + 1
  }
  names
}

///|
fn decontx_encode_clusters(
  labels : Array[String],
  n_cells : Int,
) -> (Array[String], Array[Int]) raise DecontXError {
  if labels.length() != n_cells {
    raise DecontXError::DecontXError(
      "cluster labels must match the number of cells",
    )
  }
  let names : Array[String] = Array::new()
  let indices : Array[Int] = Array::new()
  let mut cell = 0
  while cell < labels.length() {
    let label = labels[cell]
    if label.trim().length() == 0 {
      raise DecontXError::DecontXError("cluster labels cannot be empty")
    }
    let mut found = -1
    let mut cluster = 0
    while cluster < names.length() {
      if names[cluster] == label {
        found = cluster
        break
      }
      cluster = cluster + 1
    }
    if found < 0 {
      names.push(label)
      found = names.length() - 1
    }
    indices.push(found)
    cell = cell + 1
  }
  (names, indices)
}

///|
fn decontx_normalize_profiles(
  sums : Array[Array[Double]],
  pseudocount : Double,
) -> Array[Array[Double]] {
  let profiles = decontx_zero_matrix(
    sums.length(),
    if sums.length() > 0 {
      sums[0].length()
    } else {
      0
    },
  )
  let mut cluster = 0
  while cluster < sums.length() {
    let mut total = pseudocount * sums[cluster].length().to_double()
    let mut gene = 0
    while gene < sums[cluster].length() {
      total = total + sums[cluster][gene]
      gene = gene + 1
    }
    gene = 0
    while gene < sums[cluster].length() {
      profiles[cluster][gene] = (sums[cluster][gene] + pseudocount) / total
      gene = gene + 1
    }
    cluster = cluster + 1
  }
  profiles
}

///|
fn decontx_initial_native_sums(
  counts : Array[Array[Double]],
  cluster_indices : Array[Int],
  n_clusters : Int,
) -> Array[Array[Double]] {
  let sums = decontx_zero_matrix(n_clusters, counts.length())
  let n_cells = if counts.length() > 0 { counts[0].length() } else { 0 }
  let mut cell = 0
  while cell < n_cells {
    let cluster = cluster_indices[cell]
    let mut gene = 0
    while gene < counts.length() {
      sums[cluster][gene] = sums[cluster][gene] + counts[gene][cell]
      gene = gene + 1
    }
    cell = cell + 1
  }
  sums
}

///|
fn decontx_background_profile(
  background : Array[Array[Double]],
  n_genes : Int,
  pseudocount : Double,
) -> Array[Double] raise DecontXError {
  let (_, n_droplets) = decontx_validate_matrix(
    background,
    "background",
    expected_rows=n_genes,
  )
  let profile = Array::make(n_genes, pseudocount)
  let mut total = pseudocount * n_genes.to_double()
  let mut gene = 0
  while gene < n_genes {
    let mut droplet = 0
    while droplet < n_droplets {
      profile[gene] = profile[gene] + background[gene][droplet]
      total = total + background[gene][droplet]
      droplet = droplet + 1
    }
    gene = gene + 1
  }
  gene = 0
  while gene < n_genes {
    profile[gene] = profile[gene] / total
    gene = gene + 1
  }
  profile
}

///|
fn decontx_contaminant_profiles(
  native_sums : Array[Array[Double]],
  pseudocount : Double,
  background_profile : Array[Double]?,
) -> Array[Array[Double]] {
  let n_clusters = native_sums.length()
  let n_genes = if n_clusters > 0 { native_sums[0].length() } else { 0 }
  let profiles = decontx_zero_matrix(n_clusters, n_genes)
  match background_profile {
    Some(background) => {
      let mut cluster = 0
      while cluster < n_clusters {
        let mut gene = 0
        while gene < n_genes {
          profiles[cluster][gene] = background[gene]
          gene = gene + 1
        }
        cluster = cluster + 1
      }
    }
    None => {
      let mut cluster = 0
      while cluster < n_clusters {
        let mut total = pseudocount * n_genes.to_double()
        let mut gene = 0
        while gene < n_genes {
          let mut other = 0
          let mut value = pseudocount
          while other < n_clusters {
            if other != cluster {
              value = value + native_sums[other][gene]
            }
            other = other + 1
          }
          profiles[cluster][gene] = value
          total = total + value - pseudocount
          gene = gene + 1
        }
        gene = 0
        while gene < n_genes {
          profiles[cluster][gene] = profiles[cluster][gene] / total
          gene = gene + 1
        }
        cluster = cluster + 1
      }
    }
  }
  profiles
}

///|
fn decontx_expectation(
  counts : Array[Array[Double]],
  cluster_indices : Array[Int],
  contamination : Array[Double],
  native_profiles : Array[Array[Double]],
  contaminant_profiles : Array[Array[Double]],
  minimum_probability : Double,
) -> (
  Array[Array[Double]],
  Array[Array[Double]],
  Array[Array[Double]],
  Array[Double],
  Double,
) {
  let n_genes = counts.length()
  let n_cells = if n_genes > 0 { counts[0].length() } else { 0 }
  let n_clusters = native_profiles.length()
  let corrected = decontx_zero_matrix(n_genes, n_cells)
  let contaminant = decontx_zero_matrix(n_genes, n_cells)
  let native_sums = decontx_zero_matrix(n_clusters, n_genes)
  let contaminant_totals = Array::make(n_cells, 0.0)
  let mut log_likelihood = 0.0
  let mut cell = 0
  while cell < n_cells {
    let cluster = cluster_indices[cell]
    let rho = contamination[cell]
    let mut gene = 0
    while gene < n_genes {
      let count = counts[gene][cell]
      let native_weight = (1.0 - rho) * native_profiles[cluster][gene]
      let contaminant_weight = rho * contaminant_profiles[cluster][gene]
      let mixture = native_weight + contaminant_weight
      let denominator = if mixture > minimum_probability {
        mixture
      } else {
        minimum_probability
      }
      let native_probability = native_weight / denominator
      let expected_native = count * native_probability
      let expected_contaminant = count - expected_native
      corrected[gene][cell] = expected_native
      contaminant[gene][cell] = expected_contaminant
      native_sums[cluster][gene] = native_sums[cluster][gene] + expected_native
      contaminant_totals[cell] = contaminant_totals[cell] + expected_contaminant
      if count > 0.0 {
        log_likelihood = log_likelihood + count * @math.ln(denominator)
      }
      gene = gene + 1
    }
    cell = cell + 1
  }
  (corrected, contaminant, native_sums, contaminant_totals, log_likelihood)
}

///|
fn decontx_library_sizes(counts : Array[Array[Double]]) -> Array[Double] {
  let n_cells = if counts.length() > 0 { counts[0].length() } else { 0 }
  let sizes = Array::make(n_cells, 0.0)
  let mut gene = 0
  while gene < counts.length() {
    let mut cell = 0
    while cell < n_cells {
      sizes[cell] = sizes[cell] + counts[gene][cell]
      cell = cell + 1
    }
    gene = gene + 1
  }
  sizes
}

///|
/// Estimate and remove ambient RNA contamination from a gene-by-cell matrix.
///
/// With no background matrix, each cluster's contaminant profile is estimated
/// from the native profiles of all other clusters. A supplied background
/// matrix is interpreted as genes-by-empty-droplets and fixes a common ambient
/// profile during inference.
pub fn decontx(
  counts : Array[Array[Double]],
  gene_names : Array[String],
  cell_names : Array[String],
  cluster_labels : Array[String],
  config? : DecontXConfig = DecontXConfig::default(),
  background? : Array[Array[Double]] = [],
) -> DecontXResult raise DecontXError {
  let (n_genes, n_cells) = decontx_validate_matrix(counts, "counts")
  let genes = decontx_prepare_names(gene_names, n_genes, "Gene")
  let cells = decontx_prepare_names(cell_names, n_cells, "Cell")
  let (cluster_names, cluster_indices) = decontx_encode_clusters(
    cluster_labels, n_cells,
  )
  let used_background = background.length() > 0
  if cluster_names.length() < 2 && !used_background {
    raise DecontXError::DecontXError(
      "at least two clusters are required without a background matrix",
    )
  }
  let background_profile = if used_background {
    Some(
      decontx_background_profile(
        background,
        n_genes,
        config.profile_pseudocount,
      ),
    )
  } else {
    None
  }
  let initial_sums = decontx_initial_native_sums(
    counts,
    cluster_indices,
    cluster_names.length(),
  )
  let mut native_profiles = decontx_normalize_profiles(
    initial_sums,
    config.profile_pseudocount,
  )
  let mut contaminant_profiles = decontx_contaminant_profiles(
    initial_sums,
    config.profile_pseudocount,
    background_profile,
  )
  let contamination = Array::make(n_cells, config.initial_contamination)
  let library_sizes = decontx_library_sizes(counts)
  let log_likelihoods : Array[Double] = Array::new()
  let mut iterations = 0
  let mut converged = false
  while iterations < config.max_iterations && !converged {
    let (_, _, native_sums, contaminant_totals, log_likelihood) = decontx_expectation(
      counts,
      cluster_indices,
      contamination,
      native_profiles,
      contaminant_profiles,
      config.minimum_probability,
    )
    log_likelihoods.push(log_likelihood)
    let next_native_profiles = decontx_normalize_profiles(
      native_sums,
      config.profile_pseudocount,
    )
    let next_contaminant_profiles = decontx_contaminant_profiles(
      native_sums,
      config.profile_pseudocount,
      background_profile,
    )
    let next_contamination = Array::make(n_cells, 0.0)
    let mut max_change = 0.0
    let mut cell = 0
    while cell < n_cells {
      let posterior = (
          contaminant_totals[cell] + config.contamination_prior_alpha
        ) /
        (
          library_sizes[cell] +
          config.contamination_prior_alpha +
          config.contamination_prior_beta
        )
      let bounded = if posterior < config.minimum_probability {
        config.minimum_probability
      } else if posterior > 1.0 - config.minimum_probability {
        1.0 - config.minimum_probability
      } else {
        posterior
      }
      next_contamination[cell] = bounded
      let change = (bounded - contamination[cell]).abs()
      if change > max_change {
        max_change = change
      }
      cell = cell + 1
    }
    let mut cluster = 0
    while cluster < native_profiles.length() {
      let mut gene = 0
      while gene < n_genes {
        let native_change = (next_native_profiles[cluster][gene] -
        native_profiles[cluster][gene]).abs()
        if native_change > max_change {
          max_change = native_change
        }
        let contaminant_change = (next_contaminant_profiles[cluster][gene] -
        contaminant_profiles[cluster][gene]).abs()
        if contaminant_change > max_change {
          max_change = contaminant_change
        }
        gene = gene + 1
      }
      cluster = cluster + 1
    }
    cell = 0
    while cell < n_cells {
      contamination[cell] = next_contamination[cell]
      cell = cell + 1
    }
    native_profiles = next_native_profiles
    contaminant_profiles = next_contaminant_profiles
    iterations = iterations + 1
    converged = max_change <= config.tolerance
  }
  let (corrected, contaminant, _, _, final_log_likelihood) = decontx_expectation(
    counts,
    cluster_indices,
    contamination,
    native_profiles,
    contaminant_profiles,
    config.minimum_probability,
  )
  let initial_log_likelihood = if log_likelihoods.length() > 0 {
    log_likelihoods[0]
  } else {
    final_log_likelihood
  }
  DecontXResult::{
    corrected_counts: corrected,
    contamination_counts: contaminant,
    contamination,
    native_profiles,
    contaminant_profiles,
    cluster_names,
    cluster_indices,
    gene_names: genes,
    cell_names: cells,
    used_background,
    diagnostics: DecontXDiagnostics::{
      iterations,
      converged,
      initial_log_likelihood,
      final_log_likelihood,
      log_likelihoods,
    },
  }
}

///|
/// Infer initial populations with deterministic k-means on log-normalized
/// cell profiles before running decontamination.
pub fn decontx_auto(
  counts : Array[Array[Double]],
  gene_names : Array[String],
  cell_names : Array[String],
  n_clusters : Int,
  config? : DecontXConfig = DecontXConfig::default(),
  background? : Array[Array[Double]] = [],
) -> DecontXResult raise DecontXError {
  let (n_genes, n_cells) = decontx_validate_matrix(counts, "counts")
  if n_clusters < 2 || n_clusters > n_cells {
    raise DecontXError::DecontXError(
      "n_clusters must be between 2 and the number of cells",
    )
  }
  let library_sizes = decontx_library_sizes(counts)
  let points : Array[Point] = Array::new()
  let mut cell = 0
  while cell < n_cells {
    let point : Array[Double] = Array::new()
    let scale = if library_sizes[cell] > 0.0 {
      10000.0 / library_sizes[cell]
    } else {
      1.0
    }
    let mut gene = 0
    while gene < n_genes {
      point.push(@math.ln(1.0 + counts[gene][cell] * scale))
      gene = gene + 1
    }
    points.push(point)
    cell = cell + 1
  }
  let model = KMeans::new(n_clusters, 100, 1.0e-8).fit(points)
  let labels : Array[String] = Array::new()
  let assigned = model.get_labels()
  cell = 0
  while cell < assigned.length() {
    labels.push("cluster_" + (assigned[cell] + 1).to_string())
    cell = cell + 1
  }
  decontx(counts, gene_names, cell_names, labels, config~, background~)
}

///|
pub fn DecontXResult::n_genes(self : DecontXResult) -> Int {
  self.corrected_counts.length()
}

///|
pub fn DecontXResult::n_cells(self : DecontXResult) -> Int {
  if self.corrected_counts.length() > 0 {
    self.corrected_counts[0].length()
  } else {
    0
  }
}

///|
pub fn DecontXResult::n_clusters(self : DecontXResult) -> Int {
  self.cluster_names.length()
}

///|
pub fn DecontXResult::mean_contamination(self : DecontXResult) -> Double {
  if self.contamination.length() == 0 {
    return 0.0
  }
  let mut total = 0.0
  let mut cell = 0
  while cell < self.contamination.length() {
    total = total + self.contamination[cell]
    cell = cell + 1
  }
  total / self.contamination.length().to_double()
}

///|
pub fn DecontXResult::gene_index(self : DecontXResult, name : String) -> Int? {
  let mut gene = 0
  while gene < self.gene_names.length() {
    if self.gene_names[gene] == name {
      return Some(gene)
    }
    gene = gene + 1
  }
  None
}

///|
pub fn DecontXResult::cell_index(self : DecontXResult, name : String) -> Int? {
  let mut cell = 0
  while cell < self.cell_names.length() {
    if self.cell_names[cell] == name {
      return Some(cell)
    }
    cell = cell + 1
  }
  None
}

///|
pub fn DecontXResult::cell_estimate(
  self : DecontXResult,
  index : Int,
) -> DecontXCellEstimate? {
  if index < 0 || index >= self.n_cells() {
    return None
  }
  let mut native = 0.0
  let mut contaminant = 0.0
  let mut gene = 0
  while gene < self.n_genes() {
    native = native + self.corrected_counts[gene][index]
    contaminant = contaminant + self.contamination_counts[gene][index]
    gene = gene + 1
  }
  Some(DecontXCellEstimate::{
    cell_name: self.cell_names[index],
    cluster: self.cluster_names[self.cluster_indices[index]],
    contamination: self.contamination[index],
    library_size: native + contaminant,
    native_counts: native,
    contaminant_counts: contaminant,
  })
}

///|
pub fn DecontXResult::cluster_estimates(
  self : DecontXResult,
) -> Array[DecontXClusterEstimate] {
  let result : Array[DecontXClusterEstimate] = Array::new()
  let mut cluster = 0
  while cluster < self.cluster_names.length() {
    let mut n_cells = 0
    let mut contamination_sum = 0.0
    let mut native = 0.0
    let mut contaminant = 0.0
    let mut cell = 0
    while cell < self.n_cells() {
      if self.cluster_indices[cell] == cluster {
        n_cells = n_cells + 1
        contamination_sum = contamination_sum + self.contamination[cell]
        let mut gene = 0
        while gene < self.n_genes() {
          native = native + self.corrected_counts[gene][cell]
          contaminant = contaminant + self.contamination_counts[gene][cell]
          gene = gene + 1
        }
      }
      cell = cell + 1
    }
    result.push(DecontXClusterEstimate::{
      cluster: self.cluster_names[cluster],
      n_cells,
      mean_contamination: if n_cells > 0 {
        contamination_sum / n_cells.to_double()
      } else {
        0.0
      },
      total_counts: native + contaminant,
      native_counts: native,
      contaminant_counts: contaminant,
    })
    cluster = cluster + 1
  }
  result
}

///|
pub fn DecontXResult::most_contaminated_cells(
  self : DecontXResult,
  n : Int,
) -> Array[DecontXCellEstimate] {
  if n <= 0 {
    return []
  }
  let indices : Array[Int] = Array::new()
  let mut cell = 0
  while cell < self.n_cells() {
    indices.push(cell)
    cell = cell + 1
  }
  let mut position = 0
  while position < indices.length() {
    let mut best = position
    let mut candidate = position + 1
    while candidate < indices.length() {
      if self.contamination[indices[candidate]] >
        self.contamination[indices[best]] {
        best = candidate
      }
      candidate = candidate + 1
    }
    if best != position {
      let temporary = indices[position]
      indices[position] = indices[best]
      indices[best] = temporary
    }
    position = position + 1
  }
  let result : Array[DecontXCellEstimate] = Array::new()
  let limit = if n < indices.length() { n } else { indices.length() }
  position = 0
  while position < limit {
    match self.cell_estimate(indices[position]) {
      Some(estimate) => result.push(estimate)
      None => ()
    }
    position = position + 1
  }
  result
}

///|
pub fn DecontXResult::summary(self : DecontXResult) -> String {
  "DecontX: " +
  self.n_genes().to_string() +
  " genes x " +
  self.n_cells().to_string() +
  " cells, " +
  self.n_clusters().to_string() +
  " clusters; mean contamination=" +
  self.mean_contamination().to_string() +
  ", iterations=" +
  self.diagnostics.iterations.to_string() +
  ", converged=" +
  self.diagnostics.converged.to_string() +
  ", background=" +
  self.used_background.to_string()
}

///|
fn decontx_copy_sce(sce : SingleCellExperiment) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  let assay_names = sce.assays.keys().collect()
  let mut index = 0
  while index < assay_names.length() {
    assays[assay_names[index]] = sce.assays[assay_names[index]]
    index = index + 1
  }
  let row_data : Map[String, Array[String]] = Map([])
  let row_keys = sce.row_data.keys().collect()
  index = 0
  while index < row_keys.length() {
    row_data[row_keys[index]] = sce.row_data[row_keys[index]]
    index = index + 1
  }
  let col_data : Map[String, Array[String]] = Map([])
  let col_keys = sce.col_data.keys().collect()
  index = 0
  while index < col_keys.length() {
    col_data[col_keys[index]] = sce.col_data[col_keys[index]]
    index = index + 1
  }
  let reduced_dims : Map[String, Array[Array[Double]]] = Map([])
  let reduced_keys = sce.reduced_dims.keys().collect()
  index = 0
  while index < reduced_keys.length() {
    reduced_dims[reduced_keys[index]] = sce.reduced_dims[reduced_keys[index]]
    index = index + 1
  }
  let alternative_experiments : Map[String, SingleCellExperiment] = Map([])
  let alternative_keys = sce.alternative_experiments.keys().collect()
  index = 0
  while index < alternative_keys.length() {
    alternative_experiments[alternative_keys[index]] = sce.alternative_experiments[alternative_keys[index]]
    index = index + 1
  }
  let metadata : Map[String, String] = Map([])
  let metadata_keys = sce.metadata.keys().collect()
  index = 0
  while index < metadata_keys.length() {
    metadata[metadata_keys[index]] = sce.metadata[metadata_keys[index]]
    index = index + 1
  }
  SingleCellExperiment::{
    assays,
    row_data,
    col_data,
    reduced_dims,
    alternative_experiments,
    metadata,
    row_names: sce.row_names.copy(),
    col_names: sce.col_names.copy(),
  }
}

///|
/// Run decontamination on an SCE assay and return an immutable enriched copy.
pub fn decontx_sce(
  sce : SingleCellExperiment,
  cluster_column? : String = "cluster",
  input_assay? : String = "counts",
  output_assay? : String = "decontXcounts",
  config? : DecontXConfig = DecontXConfig::default(),
  background? : Array[Array[Double]] = [],
) -> DecontXSCEOutput raise DecontXError {
  let counts = sce_get_assay(sce, input_assay)
  if counts.length() == 0 {
    raise DecontXError::DecontXError(
      "SingleCellExperiment assay not found: " + input_assay,
    )
  }
  let clusters = sce_get_col_data(sce, cluster_column)
  if clusters.length() == 0 {
    raise DecontXError::DecontXError(
      "SingleCellExperiment column data not found: " + cluster_column,
    )
  }
  let result = decontx(
    counts,
    sce.row_names,
    sce.col_names,
    clusters,
    config~,
    background~,
  )
  let enriched = decontx_copy_sce(sce)
  enriched.assays[output_assay] = decontx_copy_matrix(result.corrected_counts)
  let contamination_strings : Array[String] = Array::new()
  let mut cell = 0
  while cell < result.contamination.length() {
    contamination_strings.push(result.contamination[cell].to_string())
    cell = cell + 1
  }
  enriched.col_data["decontX_contamination"] = contamination_strings
  enriched.col_data["decontX_cluster"] = clusters.copy()
  enriched.metadata["decontX_iterations"] = result.diagnostics.iterations.to_string()
  enriched.metadata["decontX_converged"] = result.diagnostics.converged.to_string()
  enriched.metadata["decontX_background"] = result.used_background.to_string()
  DecontXSCEOutput::{ experiment: enriched, result, }
}

///|
/// Small two-population matrix with cross-population marker contamination.
pub fn decontx_example_data() -> (
  Array[Array[Double]],
  Array[String],
  Array[String],
  Array[String],
  Array[Array[Double]],
) {
  let counts = [
    [100.0, 90.0, 110.0, 95.0, 8.0, 5.0, 10.0, 6.0],
    [80.0, 70.0, 85.0, 75.0, 4.0, 3.0, 5.0, 4.0],
    [8.0, 6.0, 10.0, 7.0, 100.0, 90.0, 110.0, 95.0],
    [4.0, 3.0, 5.0, 4.0, 80.0, 70.0, 85.0, 75.0],
    [20.0, 22.0, 19.0, 21.0, 20.0, 22.0, 19.0, 21.0],
    [2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0],
  ]
  let gene_names = [
    "T_marker_1", "T_marker_2", "B_marker_1", "B_marker_2", "housekeeping", "low_expression",
  ]
  let cell_names = ["T1", "T2", "T3", "T4", "B1", "B2", "B3", "B4"]
  let clusters = ["T", "T", "T", "T", "B", "B", "B", "B"]
  let background = [
    [45.0, 50.0],
    [35.0, 40.0],
    [48.0, 52.0],
    [37.0, 42.0],
    [18.0, 20.0],
    [2.0, 2.0],
  ]
  (counts, gene_names, cell_names, clusters, background)
}
