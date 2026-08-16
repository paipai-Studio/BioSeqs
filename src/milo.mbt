///|
/// Graph-neighbourhood differential abundance analysis inspired by
/// Bioconductor miloR.
///
/// The implementation covers exact KNN graph construction, refined
/// neighbourhood sampling, neighbourhood-by-sample counting, negative
/// binomial GLM testing and graph-aware weighted FDR correction.

///|
pub suberror MiloError {
  MiloError(String)
}

///|
pub enum MiloSpatialWeighting {
  KDistance
  NeighbourDistance
  MaxDistance
  GraphOverlap
  NoSpatialWeighting
} derive(Eq, Debug)

///|
pub struct MiloExperiment {
  coordinates : Array[Array[Double]]
  cell_names : Array[String]
  dimensions : Int
  k : Int
  knn_indices : Array[Array[Int]]
  knn_distances : Array[Array[Double]]
  graph_neighbors : Array[Array[Int]]
  neighborhood_indices : Array[Int]
  neighborhoods : Array[Array[Int]]
  sample_names : Array[String]
  neighborhood_counts : Array[Array[Int]]
} derive(Debug)

///|
pub struct MiloDesign {
  matrix : Array[Array[Double]]
  levels : Array[String]
  reference : String
  coefficient_name : String
} derive(Debug)

///|
pub struct MiloDAResult {
  neighborhood : Int
  index_cell : Int
  tested : Bool
  log_fc : Double
  log_cpm : Double
  statistic : Double
  p_value : Double
  fdr : Double
  spatial_fdr : Double
  dispersion : Double
  mean_count : Double
} derive(Debug)

///|
pub struct MiloSampleData {
  coordinates : Array[Array[Double]]
  cell_names : Array[String]
  sample_ids : Array[String]
  sample_names : Array[String]
  sample_conditions : Array[String]
} derive(Debug)

///|
priv struct MiloFit {
  log_fc : Double
  log_cpm : Double
  statistic : Double
  p_value : Double
  dispersion : Double
  mean_count : Double
}

///|
pub fn milo_k_distance_weighting() -> MiloSpatialWeighting {
  KDistance
}

///|
pub fn milo_neighbour_distance_weighting() -> MiloSpatialWeighting {
  NeighbourDistance
}

///|
pub fn milo_max_distance_weighting() -> MiloSpatialWeighting {
  MaxDistance
}

///|
pub fn milo_graph_overlap_weighting() -> MiloSpatialWeighting {
  GraphOverlap
}

///|
pub fn milo_no_spatial_weighting() -> MiloSpatialWeighting {
  NoSpatialWeighting
}

///|
fn milo_is_finite(value : Double) -> Bool {
  value == value && value.abs() < 1.0e300
}

///|
fn milo_copy_double_matrix(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in matrix {
    output.push(row.copy())
  }
  output
}

///|
fn milo_copy_int_matrix(matrix : Array[Array[Int]]) -> Array[Array[Int]] {
  let output : Array[Array[Int]] = []
  for row in matrix {
    output.push(row.copy())
  }
  output
}

///|
fn milo_validate_coordinates(
  coordinates : Array[Array[Double]],
  dimensions : Int,
) -> (Int, Int) raise MiloError {
  let cells = coordinates.length()
  if cells < 2 {
    raise MiloError("miloR graph construction requires at least two cells")
  }
  let available_dimensions = coordinates[0].length()
  if available_dimensions == 0 {
    raise MiloError("miloR coordinates must contain at least one dimension")
  }
  for row in coordinates {
    if row.length() != available_dimensions {
      raise MiloError("miloR coordinates must be rectangular")
    }
    for value in row {
      if !milo_is_finite(value) {
        raise MiloError("miloR coordinates must contain finite values")
      }
    }
  }
  let used_dimensions = if dimensions == 0 {
    available_dimensions
  } else {
    dimensions
  }
  if used_dimensions <= 0 || used_dimensions > available_dimensions {
    raise MiloError(
      "miloR dimensions must be positive and no larger than the coordinate rank",
    )
  }
  (cells, used_dimensions)
}

///|
fn milo_validate_names(
  names : Array[String],
  cells : Int,
) -> Array[String] raise MiloError {
  if names.length() == 0 {
    let generated : Array[String] = []
    for index in 0..<cells {
      generated.push("Cell" + (index + 1).to_string())
    }
    return generated
  }
  if names.length() != cells {
    raise MiloError("miloR cell names must match the number of cells")
  }
  let seen : Map[String, Bool] = Map([], capacity=names.length())
  for name in names {
    if name == "" {
      raise MiloError("miloR cell names must not be empty")
    }
    if seen.contains(name) {
      raise MiloError("miloR cell names must be unique")
    }
    seen[name] = true
  }
  names.copy()
}

///|
fn milo_distance(
  left : Array[Double],
  right : Array[Double],
  dimensions : Int,
) -> Double {
  let mut total = 0.0
  for dimension in 0..<dimensions {
    let difference = left[dimension] - right[dimension]
    total = total + difference * difference
  }
  total.sqrt()
}

///|
fn milo_sort_ints(values : Array[Int]) -> Unit {
  values.sort_by(fn(left : Int, right : Int) -> Int {
    if left < right {
      -1
    } else if left > right {
      1
    } else {
      0
    }
  })
}

///|
fn milo_build_knn(
  coordinates : Array[Array[Double]],
  dimensions : Int,
  k : Int,
) -> (Array[Array[Int]], Array[Array[Double]]) {
  let cells = coordinates.length()
  let indices : Array[Array[Int]] = []
  let distances : Array[Array[Double]] = []
  for cell in 0..<cells {
    let candidates : Array[(Double, Int)] = []
    for other in 0..<cells {
      if other != cell {
        candidates.push(
          (
            milo_distance(coordinates[cell], coordinates[other], dimensions),
            other,
          ),
        )
      }
    }
    candidates.sort_by(fn(left : (Double, Int), right : (Double, Int)) -> Int {
      if left.0 < right.0 {
        -1
      } else if left.0 > right.0 {
        1
      } else if left.1 < right.1 {
        -1
      } else if left.1 > right.1 {
        1
      } else {
        0
      }
    })
    let cell_indices : Array[Int] = []
    let cell_distances : Array[Double] = []
    for neighbor in 0..<k {
      cell_indices.push(candidates[neighbor].1)
      cell_distances.push(candidates[neighbor].0)
    }
    indices.push(cell_indices)
    distances.push(cell_distances)
  }
  (indices, distances)
}

///|
fn milo_symmetrize_knn(knn_indices : Array[Array[Int]]) -> Array[Array[Int]] {
  let neighbors : Array[Array[Int]] = []
  for _ in 0..<knn_indices.length() {
    neighbors.push([])
  }
  for cell in 0..<knn_indices.length() {
    for other in knn_indices[cell] {
      if !neighbors[cell].contains(other) {
        neighbors[cell].push(other)
      }
      if !neighbors[other].contains(cell) {
        neighbors[other].push(cell)
      }
    }
  }
  for row in neighbors {
    milo_sort_ints(row)
  }
  neighbors
}

///|
/// Construct an exact, undirected KNN graph from cells by reduced dimensions.
///
/// `coordinates` uses cells in rows and dimensions in columns. `dimensions=0`
/// uses every available dimension. Cell names are generated when omitted.
pub fn milo_build_graph(
  coordinates : Array[Array[Double]],
  cell_names? : Array[String] = [],
  k? : Int = 10,
  dimensions? : Int = 0,
) -> MiloExperiment raise MiloError {
  let (cells, used_dimensions) = milo_validate_coordinates(
    coordinates, dimensions,
  )
  if k <= 0 || k >= cells {
    raise MiloError("miloR k must be positive and smaller than the cell count")
  }
  let validated_names = milo_validate_names(cell_names, cells)
  let copied_coordinates = milo_copy_double_matrix(coordinates)
  let (knn_indices, knn_distances) = milo_build_knn(
    copied_coordinates, used_dimensions, k,
  )
  let graph_neighbors = milo_symmetrize_knn(knn_indices)
  MiloExperiment::{
    coordinates: copied_coordinates,
    cell_names: validated_names,
    dimensions: used_dimensions,
    k,
    knn_indices,
    knn_distances,
    graph_neighbors,
    neighborhood_indices: [],
    neighborhoods: [],
    sample_names: [],
    neighborhood_counts: [],
  }
}

///|
/// Construct a Milo graph from a `SingleCellExperiment` reduced dimension.
pub fn milo_from_single_cell_experiment(
  sce : SingleCellExperiment,
  reduced_dim? : String = "PCA",
  k? : Int = 10,
  dimensions? : Int = 0,
) -> MiloExperiment raise MiloError {
  let coordinates = sce_get_reduced_dim(sce, reduced_dim)
  if coordinates.length() == 0 {
    raise MiloError(
      "miloR reduced dimension '" + reduced_dim + "' is missing or empty",
    )
  }
  milo_build_graph(coordinates, cell_names=sce.col_names, k~, dimensions~)
}

///|
fn milo_rng_next(state : Int) -> Int {
  let high = state / 44488
  let low = state % 44488
  let candidate = 48271 * low - 3399 * high
  if candidate > 0 {
    candidate
  } else {
    candidate + 2147483647
  }
}

///|
fn milo_sample_vertices(cells : Int, count : Int, seed : Int) -> Array[Int] {
  let sampled : Array[Int] = []
  let used = Array::make(cells, false)
  let normalized_seed = if seed == 0 {
    1
  } else {
    let positive = seed.abs() % 2147483646
    positive + 1
  }
  let mut state = normalized_seed
  while sampled.length() < count {
    state = milo_rng_next(state)
    let candidate = state % cells
    if !used[candidate] {
      used[candidate] = true
      sampled.push(candidate)
    }
  }
  sampled
}

///|
fn milo_refine_vertex(experiment : MiloExperiment, vertex : Int) -> Int {
  let median_profile = Array::make(experiment.dimensions, 0.0)
  for dimension in 0..<experiment.dimensions {
    let values : Array[Double] = []
    for neighbor in experiment.knn_indices[vertex] {
      values.push(experiment.coordinates[neighbor][dimension])
    }
    median_profile[dimension] = stat_median(values)
  }
  let mut nearest = 0
  let mut nearest_distance = milo_distance(
    median_profile,
    experiment.coordinates[0],
    experiment.dimensions,
  )
  for cell in 1..<experiment.coordinates.length() {
    let distance = milo_distance(
      median_profile,
      experiment.coordinates[cell],
      experiment.dimensions,
    )
    if distance < nearest_distance {
      nearest = cell
      nearest_distance = distance
    }
  }
  nearest
}

///|
/// Sample representative graph vertices and define their order-one
/// neighbourhoods. Refined sampling uses the component-wise median of each
/// seed's KNN profile and selects the nearest observed cell.
pub fn MiloExperiment::make_neighborhoods(
  self : MiloExperiment,
  proportion? : Double = 0.1,
  refined? : Bool = true,
  seed? : Int = 42,
) -> MiloExperiment raise MiloError {
  if !milo_is_finite(proportion) || proportion <= 0.0 || proportion >= 1.0 {
    raise MiloError("miloR sampling proportion must be in (0, 1)")
  }
  let cells = self.coordinates.length()
  let requested = (proportion * cells.to_double()).to_int()
  let sample_count = requested.max(1).min(cells)
  let initial = milo_sample_vertices(cells, sample_count, seed)
  let indices : Array[Int] = []
  for vertex in initial {
    let selected = if refined {
      milo_refine_vertex(self, vertex)
    } else {
      vertex
    }
    if !indices.contains(selected) {
      indices.push(selected)
    }
  }
  let neighborhoods : Array[Array[Int]] = []
  for index_cell in indices {
    let members : Array[Int] = [index_cell]
    for neighbor in self.graph_neighbors[index_cell] {
      if !members.contains(neighbor) {
        members.push(neighbor)
      }
    }
    milo_sort_ints(members)
    neighborhoods.push(members)
  }
  MiloExperiment::{
    coordinates: milo_copy_double_matrix(self.coordinates),
    cell_names: self.cell_names.copy(),
    dimensions: self.dimensions,
    k: self.k,
    knn_indices: milo_copy_int_matrix(self.knn_indices),
    knn_distances: milo_copy_double_matrix(self.knn_distances),
    graph_neighbors: milo_copy_int_matrix(self.graph_neighbors),
    neighborhood_indices: indices,
    neighborhoods,
    sample_names: [],
    neighborhood_counts: [],
  }
}

///|
/// Count cells from each experimental sample in every neighbourhood.
pub fn MiloExperiment::count_cells(
  self : MiloExperiment,
  sample_ids : Array[String],
) -> MiloExperiment raise MiloError {
  if self.neighborhoods.length() == 0 {
    raise MiloError(
      "miloR neighborhoods are missing; run make_neighborhoods first",
    )
  }
  if sample_ids.length() != self.coordinates.length() {
    raise MiloError("miloR sample IDs must match the number of cells")
  }
  let sample_names : Array[String] = []
  let sample_lookup : Map[String, Int] = Map([], capacity=sample_ids.length())
  for sample in sample_ids {
    if sample == "" {
      raise MiloError("miloR sample IDs must not be empty")
    }
    if !sample_lookup.contains(sample) {
      sample_lookup[sample] = sample_names.length()
      sample_names.push(sample)
    }
  }
  let counts : Array[Array[Int]] = []
  for neighborhood in self.neighborhoods {
    let row = Array::make(sample_names.length(), 0)
    for cell in neighborhood {
      let sample_index = sample_lookup[sample_ids[cell]]
      row[sample_index] = row[sample_index] + 1
    }
    counts.push(row)
  }
  MiloExperiment::{
    coordinates: milo_copy_double_matrix(self.coordinates),
    cell_names: self.cell_names.copy(),
    dimensions: self.dimensions,
    k: self.k,
    knn_indices: milo_copy_int_matrix(self.knn_indices),
    knn_distances: milo_copy_double_matrix(self.knn_distances),
    graph_neighbors: milo_copy_int_matrix(self.graph_neighbors),
    neighborhood_indices: self.neighborhood_indices.copy(),
    neighborhoods: milo_copy_int_matrix(self.neighborhoods),
    sample_names,
    neighborhood_counts: counts,
  }
}

///|
/// Compute feature means across the cells in each neighbourhood.
///
/// The input and output use features in rows. Input columns are cells and
/// output columns are neighbourhoods.
pub fn MiloExperiment::neighborhood_expression(
  self : MiloExperiment,
  expression : Array[Array[Double]],
) -> Array[Array[Double]] raise MiloError {
  if self.neighborhoods.length() == 0 {
    raise MiloError("miloR neighborhoods are missing")
  }
  let cells = self.coordinates.length()
  let output : Array[Array[Double]] = []
  for feature in expression {
    if feature.length() != cells {
      raise MiloError("miloR expression columns must match the number of cells")
    }
    let row = Array::make(self.neighborhoods.length(), 0.0)
    for neighborhood_index in 0..<self.neighborhoods.length() {
      let neighborhood = self.neighborhoods[neighborhood_index]
      let mut total = 0.0
      for cell in neighborhood {
        let value = feature[cell]
        if !milo_is_finite(value) {
          raise MiloError("miloR expression values must be finite")
        }
        total = total + value
      }
      row[neighborhood_index] = total / neighborhood.length().to_double()
    }
    output.push(row)
  }
  output
}

///|
/// Return pairwise cell-overlap counts between neighbourhoods.
pub fn MiloExperiment::neighborhood_overlap_matrix(
  self : MiloExperiment,
) -> Array[Array[Int]] {
  let count = self.neighborhoods.length()
  let overlaps : Array[Array[Int]] = []
  for _ in 0..<count {
    overlaps.push(Array::make(count, 0))
  }
  for left in 0..<count {
    for right in left..<count {
      let mut overlap = 0
      for cell in self.neighborhoods[left] {
        if self.neighborhoods[right].contains(cell) {
          overlap = overlap + 1
        }
      }
      overlaps[left][right] = overlap
      overlaps[right][left] = overlap
    }
  }
  overlaps
}

///|
fn milo_neighborhood_pairwise_mean(
  experiment : MiloExperiment,
  members : Array[Int],
) -> Double {
  if members.length() < 2 {
    return 0.0
  }
  let mut total = 0.0
  let mut pairs = 0
  for left in 0..<(members.length() - 1) {
    for right in (left + 1)..<members.length() {
      total = total +
        milo_distance(
          experiment.coordinates[members[left]],
          experiment.coordinates[members[right]],
          experiment.dimensions,
        )
      pairs = pairs + 1
    }
  }
  total / pairs.to_double()
}

///|
fn milo_spatial_connectivity(
  experiment : MiloExperiment,
  weighting : MiloSpatialWeighting,
) -> Array[Double] raise MiloError {
  let count = experiment.neighborhoods.length()
  if count == 0 {
    raise MiloError("miloR neighborhoods are missing")
  }
  let connectivity = Array::make(count, 1.0)
  match weighting {
    KDistance =>
      for neighborhood in 0..<count {
        let index_cell = experiment.neighborhood_indices[neighborhood]
        connectivity[neighborhood] = experiment.knn_distances[index_cell][experiment.k -
          1]
      }
    NeighbourDistance =>
      for neighborhood in 0..<count {
        connectivity[neighborhood] = milo_neighborhood_pairwise_mean(
          experiment,
          experiment.neighborhoods[neighborhood],
        )
      }
    MaxDistance =>
      for neighborhood in 0..<count {
        let index_cell = experiment.neighborhood_indices[neighborhood]
        let mut maximum = 0.0
        for cell in experiment.neighborhoods[neighborhood] {
          let distance = milo_distance(
            experiment.coordinates[index_cell],
            experiment.coordinates[cell],
            experiment.dimensions,
          )
          if distance > maximum {
            maximum = distance
          }
        }
        connectivity[neighborhood] = maximum
      }
    GraphOverlap => {
      let overlaps = experiment.neighborhood_overlap_matrix()
      for neighborhood in 0..<count {
        let mut total = 0.0
        for other in 0..<count {
          if neighborhood != other {
            total = total + overlaps[neighborhood][other].to_double()
          }
        }
        connectivity[neighborhood] = total
      }
    }
    NoSpatialWeighting => ()
  }
  connectivity
}

///|
/// Apply miloR/cydar-style frequency-weighted BH correction.
pub fn milo_graph_spatial_fdr(
  experiment : MiloExperiment,
  p_values : Array[Double],
  weighting? : MiloSpatialWeighting = KDistance,
) -> Array[Double] raise MiloError {
  let count = experiment.neighborhoods.length()
  if p_values.length() != count {
    raise MiloError("miloR p-values must match the neighborhood count")
  }
  for p_value in p_values {
    if !milo_is_finite(p_value) || p_value < 0.0 || p_value > 1.0 {
      raise MiloError("miloR p-values must be finite values in [0, 1]")
    }
  }
  if weighting is NoSpatialWeighting {
    return Array::make(count, @double.not_a_number)
  }
  let connectivity = milo_spatial_connectivity(experiment, weighting)
  let weights = Array::make(count, 1.0)
  for index in 0..<count {
    let value = connectivity[index]
    weights[index] = if milo_is_finite(value) && value > 1.0e-12 {
      1.0 / value
    } else {
      1.0
    }
  }
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
  let mut total_weight = 0.0
  for weight in weights {
    total_weight = total_weight + weight
  }
  let raw = Array::make(count, 1.0)
  let mut cumulative_weight = 0.0
  for rank in 0..<count {
    let index = order[rank]
    cumulative_weight = cumulative_weight + weights[index]
    raw[rank] = (total_weight * p_values[index] / cumulative_weight).min(1.0)
  }
  let adjusted = Array::make(count, 1.0)
  let mut cumulative_minimum = 1.0
  let mut rank = count - 1
  while rank >= 0 {
    cumulative_minimum = cumulative_minimum.min(raw[rank])
    adjusted[order[rank]] = cumulative_minimum
    rank = rank - 1
  }
  adjusted
}

///|
/// Build an intercept-plus-binary-condition design matrix.
pub fn milo_binary_design(
  conditions : Array[String],
  reference? : String = "",
) -> MiloDesign raise MiloError {
  if conditions.length() == 0 {
    raise MiloError("miloR binary design requires at least one sample")
  }
  let levels : Array[String] = []
  for condition in conditions {
    if condition == "" {
      raise MiloError("miloR condition names must not be empty")
    }
    if !levels.contains(condition) {
      levels.push(condition)
    }
  }
  if levels.length() != 2 {
    raise MiloError("miloR binary design requires exactly two conditions")
  }
  let selected_reference = if reference == "" { levels[0] } else { reference }
  if !levels.contains(selected_reference) {
    raise MiloError("miloR reference condition is not present")
  }
  let comparison = if levels[0] == selected_reference {
    levels[1]
  } else {
    levels[0]
  }
  let matrix : Array[Array[Double]] = []
  for condition in conditions {
    matrix.push([1.0, if condition == comparison { 1.0 } else { 0.0 }])
  }
  MiloDesign::{
    matrix,
    levels: [selected_reference, comparison],
    reference: selected_reference,
    coefficient_name: comparison + "_vs_" + selected_reference,
  }
}

///|
fn milo_matrix_inverse(matrix : Array[Array[Double]]) -> Array[Array[Double]]? {
  let size = matrix.length()
  if size == 0 {
    return None
  }
  let augmented : Array[Array[Double]] = []
  for row in 0..<size {
    if matrix[row].length() != size {
      return None
    }
    let output = Array::make(2 * size, 0.0)
    for column in 0..<size {
      output[column] = matrix[row][column]
    }
    output[size + row] = 1.0
    augmented.push(output)
  }
  for pivot_column in 0..<size {
    let mut pivot_row = pivot_column
    let mut pivot_value = augmented[pivot_row][pivot_column].abs()
    for candidate in (pivot_column + 1)..<size {
      let candidate_value = augmented[candidate][pivot_column].abs()
      if candidate_value > pivot_value {
        pivot_row = candidate
        pivot_value = candidate_value
      }
    }
    if pivot_value < 1.0e-12 {
      return None
    }
    if pivot_row != pivot_column {
      let temporary = augmented[pivot_column]
      augmented[pivot_column] = augmented[pivot_row]
      augmented[pivot_row] = temporary
    }
    let divisor = augmented[pivot_column][pivot_column]
    for column in 0..<(2 * size) {
      augmented[pivot_column][column] = augmented[pivot_column][column] /
        divisor
    }
    for row in 0..<size {
      if row != pivot_column {
        let factor = augmented[row][pivot_column]
        for column in 0..<(2 * size) {
          augmented[row][column] = augmented[row][column] -
            factor * augmented[pivot_column][column]
        }
      }
    }
  }
  let inverse : Array[Array[Double]] = []
  for row in 0..<size {
    let output = Array::make(size, 0.0)
    for column in 0..<size {
      output[column] = augmented[row][size + column]
    }
    inverse.push(output)
  }
  Some(inverse)
}

///|
fn milo_matrix_vector_product(
  matrix : Array[Array[Double]],
  vector : Array[Double],
) -> Array[Double] {
  let output = Array::make(matrix.length(), 0.0)
  for row in 0..<matrix.length() {
    for column in 0..<vector.length() {
      output[row] = output[row] + matrix[row][column] * vector[column]
    }
  }
  output
}

///|
fn milo_design_information(
  design : Array[Array[Double]],
) -> Array[Array[Double]] {
  let parameters = design[0].length()
  let information : Array[Array[Double]] = []
  for _ in 0..<parameters {
    information.push(Array::make(parameters, 0.0))
  }
  for row in design {
    for left in 0..<parameters {
      for right in 0..<parameters {
        information[left][right] = information[left][right] +
          row[left] * row[right]
      }
    }
  }
  information
}

///|
fn milo_validate_design(
  design : Array[Array[Double]],
  samples : Int,
  contrast_column : Int,
) -> Int raise MiloError {
  if design.length() != samples {
    raise MiloError("miloR design rows must match neighborhood count columns")
  }
  if samples < 2 {
    raise MiloError("miloR testing requires at least two samples")
  }
  let parameters = design[0].length()
  if parameters == 0 || contrast_column < 0 || contrast_column >= parameters {
    raise MiloError("miloR contrast column is out of bounds")
  }
  if samples <= parameters {
    raise MiloError("miloR testing requires residual degrees of freedom")
  }
  for row in design {
    if row.length() != parameters {
      raise MiloError("miloR design matrix must be rectangular")
    }
    for value in row {
      if !milo_is_finite(value) {
        raise MiloError("miloR design matrix must contain finite values")
      }
    }
  }
  if milo_matrix_inverse(milo_design_information(design)) is None {
    raise MiloError("miloR design matrix must have full column rank")
  }
  parameters
}

///|
fn milo_library_offsets(
  counts : Array[Array[Int]],
  cell_sizes : Array[Double],
) -> (Array[Double], Array[Double]) raise MiloError {
  let samples = counts[0].length()
  let library_sizes = Array::make(samples, 0.0)
  if cell_sizes.length() > 0 {
    if cell_sizes.length() != samples {
      raise MiloError("miloR cell sizes must match the sample count")
    }
    for sample in 0..<samples {
      if !milo_is_finite(cell_sizes[sample]) || cell_sizes[sample] <= 0.0 {
        raise MiloError("miloR cell sizes must be finite and positive")
      }
      library_sizes[sample] = cell_sizes[sample]
    }
  } else {
    for row in counts {
      for sample in 0..<samples {
        library_sizes[sample] = library_sizes[sample] + row[sample].to_double()
      }
    }
    for size in library_sizes {
      if size <= 0.0 {
        raise MiloError("miloR neighborhood library sizes must be positive")
      }
    }
  }
  let mut mean_log = 0.0
  for size in library_sizes {
    mean_log = mean_log + @math.ln(size)
  }
  mean_log = mean_log / samples.to_double()
  let reference = @math.exp(mean_log)
  let offsets = Array::make(samples, 0.0)
  let scales = Array::make(samples, 1.0)
  for sample in 0..<samples {
    scales[sample] = library_sizes[sample] / reference
    offsets[sample] = @math.ln(scales[sample])
  }
  (offsets, scales)
}

///|
fn milo_raw_dispersion(counts : Array[Int], scales : Array[Double]) -> Double {
  let normalized = Array::make(counts.length(), 0.0)
  for sample in 0..<counts.length() {
    normalized[sample] = counts[sample].to_double() / scales[sample]
  }
  let mean = stat_mean(normalized)
  if mean <= 1.0e-12 {
    return 1.0
  }
  let variance = stat_variance(normalized)
  ((variance - mean) / (mean * mean)).max(1.0e-4).min(100.0)
}

///|
fn milo_fit_nb_glm(
  counts : Array[Int],
  design : Array[Array[Double]],
  offsets : Array[Double],
  library_sizes : Array[Double],
  contrast_column : Int,
  dispersion : Double,
) -> MiloFit {
  let samples = counts.length()
  let parameters = design[0].length()
  let coefficients = Array::make(parameters, 0.0)
  let mut count_total = 0.0
  let mut exposure_total = 0.0
  for sample in 0..<samples {
    count_total = count_total + counts[sample].to_double()
    exposure_total = exposure_total + @math.exp(offsets[sample])
  }
  coefficients[0] = @math.ln((count_total + 0.5) / exposure_total)
  let fitted = Array::make(samples, 0.0)
  let mut covariance : Array[Array[Double]]? = None
  for _ in 0..<80 {
    let weights = Array::make(samples, 0.0)
    let working = Array::make(samples, 0.0)
    for sample in 0..<samples {
      let mut linear = 0.0
      for parameter in 0..<parameters {
        linear = linear + design[sample][parameter] * coefficients[parameter]
      }
      let eta = (offsets[sample] + linear).max(-30.0).min(30.0)
      let mean = @math.exp(eta).max(1.0e-8)
      fitted[sample] = mean
      weights[sample] = mean / (1.0 + dispersion * mean)
      working[sample] = linear + (counts[sample].to_double() - mean) / mean
    }
    let information : Array[Array[Double]] = []
    for _ in 0..<parameters {
      information.push(Array::make(parameters, 0.0))
    }
    let score = Array::make(parameters, 0.0)
    for sample in 0..<samples {
      for left in 0..<parameters {
        score[left] = score[left] +
          design[sample][left] * weights[sample] * working[sample]
        for right in 0..<parameters {
          information[left][right] = information[left][right] +
            design[sample][left] * weights[sample] * design[sample][right]
        }
      }
    }
    for parameter in 0..<parameters {
      information[parameter][parameter] = information[parameter][parameter] +
        1.0e-10
    }
    match milo_matrix_inverse(information) {
      Some(inverse) => {
        let updated = milo_matrix_vector_product(inverse, score)
        let mut largest_change = 0.0
        for parameter in 0..<parameters {
          largest_change = largest_change.max(
            (updated[parameter] - coefficients[parameter]).abs(),
          )
          coefficients[parameter] = updated[parameter]
        }
        covariance = Some(inverse)
        if largest_change < 1.0e-8 {
          break
        }
      }
      None => break
    }
  }
  let standard_error = match covariance {
    Some(matrix) =>
      matrix[contrast_column][contrast_column].abs().sqrt().max(1.0e-12)
    None => 1.0e300
  }
  let statistic = coefficients[contrast_column] / standard_error
  let p_value = if milo_is_finite(statistic) {
    (2.0 * (1.0 - stat_normal_cdf(statistic.abs()))).max(0.0).min(1.0)
  } else {
    1.0
  }
  let mut total_library = 0.0
  for size in library_sizes {
    total_library = total_library + size
  }
  let log_cpm = @math.ln(count_total / total_library * 1000000.0 + 0.5) /
    @math.ln(2.0)
  MiloFit::{
    log_fc: coefficients[contrast_column] / @math.ln(2.0),
    log_cpm,
    statistic,
    p_value,
    dispersion,
    mean_count: count_total / samples.to_double(),
  }
}

///|
/// Fit a fixed-effect negative binomial GLM to every neighbourhood.
///
/// The design matrix is sample by coefficient and must follow
/// `experiment.sample_names`. The coefficient selected by `contrast_column`
/// is reported as log2 fold change. Raw dispersions are method-of-moments
/// estimates shrunk towards their median before IRLS fitting.
pub fn MiloExperiment::test_neighborhoods(
  self : MiloExperiment,
  design : Array[Array[Double]],
  contrast_column : Int,
  min_mean? : Double = 0.0,
  weighting? : MiloSpatialWeighting = KDistance,
  cell_sizes? : Array[Double] = [],
) -> Array[MiloDAResult] raise MiloError {
  if self.neighborhood_counts.length() == 0 {
    raise MiloError(
      "miloR neighborhood counts are missing; run count_cells first",
    )
  }
  if !milo_is_finite(min_mean) || min_mean < 0.0 {
    raise MiloError("miloR minimum mean count must be finite and non-negative")
  }
  let samples = self.sample_names.length()
  let parameters = milo_validate_design(design, samples, contrast_column)
  ignore(parameters)
  for row in self.neighborhood_counts {
    if row.length() != samples {
      raise MiloError("miloR neighborhood count matrix is malformed")
    }
    for value in row {
      if value < 0 {
        raise MiloError("miloR neighborhood counts must be non-negative")
      }
    }
  }
  let (offsets, scales) = milo_library_offsets(
    self.neighborhood_counts,
    cell_sizes,
  )
  let library_sizes = Array::make(samples, 0.0)
  let mut geometric_reference = 0.0
  for scale in scales {
    geometric_reference = geometric_reference + @math.ln(scale)
  }
  geometric_reference = @math.exp(geometric_reference / samples.to_double())
  for sample in 0..<samples {
    library_sizes[sample] = scales[sample] / geometric_reference
  }
  if cell_sizes.length() > 0 {
    for sample in 0..<samples {
      library_sizes[sample] = cell_sizes[sample]
    }
  } else {
    for sample in 0..<samples {
      library_sizes[sample] = 0.0
      for row in self.neighborhood_counts {
        library_sizes[sample] = library_sizes[sample] + row[sample].to_double()
      }
    }
  }
  let raw_dispersions : Array[Double] = []
  for row in self.neighborhood_counts {
    raw_dispersions.push(milo_raw_dispersion(row, scales))
  }
  let common_dispersion = stat_median(raw_dispersions).max(1.0e-4)
  let residual_degrees = (samples - design[0].length()).max(1).to_double()
  let prior_degrees = 5.0
  let fits : Array[MiloFit] = []
  let tested = Array::make(self.neighborhood_counts.length(), false)
  let p_values = Array::make(self.neighborhood_counts.length(), 1.0)
  for neighborhood in 0..<self.neighborhood_counts.length() {
    let row = self.neighborhood_counts[neighborhood]
    let mut total = 0.0
    for value in row {
      total = total + value.to_double()
    }
    let mean = total / samples.to_double()
    let shrunk_dispersion = (
        residual_degrees * raw_dispersions[neighborhood] +
        prior_degrees * common_dispersion
      ) /
      (residual_degrees + prior_degrees)
    if mean >= min_mean {
      tested[neighborhood] = true
      let fit = milo_fit_nb_glm(
        row, design, offsets, library_sizes, contrast_column, shrunk_dispersion,
      )
      p_values[neighborhood] = fit.p_value
      fits.push(fit)
    } else {
      fits.push(MiloFit::{
        log_fc: 0.0,
        log_cpm: 0.0,
        statistic: 0.0,
        p_value: 1.0,
        dispersion: shrunk_dispersion,
        mean_count: mean,
      })
    }
  }
  let fdr = stat_bh(p_values)
  let spatial_fdr = milo_graph_spatial_fdr(self, p_values, weighting~)
  let results : Array[MiloDAResult] = []
  for neighborhood in 0..<fits.length() {
    let fit = fits[neighborhood]
    results.push(MiloDAResult::{
      neighborhood,
      index_cell: self.neighborhood_indices[neighborhood],
      tested: tested[neighborhood],
      log_fc: fit.log_fc,
      log_cpm: fit.log_cpm,
      statistic: fit.statistic,
      p_value: fit.p_value,
      fdr: fdr[neighborhood],
      spatial_fdr: spatial_fdr[neighborhood],
      dispersion: fit.dispersion,
      mean_count: fit.mean_count,
    })
  }
  results
}

///|
pub fn MiloExperiment::summary(self : MiloExperiment) -> String {
  "Milo(cells=" +
  self.coordinates.length().to_string() +
  ", dimensions=" +
  self.dimensions.to_string() +
  ", k=" +
  self.k.to_string() +
  ", neighborhoods=" +
  self.neighborhoods.length().to_string() +
  ", samples=" +
  self.sample_names.length().to_string() +
  ")"
}

///|
/// Deterministic two-state, eight-sample data for examples and tests.
///
/// Condition A is enriched around the left state and condition B around the
/// right state. Every sample contributes four cells.
pub fn milo_sample_data() -> MiloSampleData {
  let coordinates : Array[Array[Double]] = []
  let cell_names : Array[String] = []
  let sample_ids : Array[String] = []
  let sample_names : Array[String] = []
  let sample_conditions : Array[String] = []
  let mut cell_number = 1
  for sample in 0..<8 {
    let condition = if sample < 4 { "A" } else { "B" }
    let sample_name = condition + (sample % 4 + 1).to_string()
    sample_names.push(sample_name)
    sample_conditions.push(condition)
    let left_count = if condition == "A" { 3 } else { 1 }
    let right_count = 4 - left_count
    for local in 0..<left_count {
      coordinates.push([
        (sample % 4).to_double() * 0.03 + local.to_double() * 0.08,
        local.to_double() * 0.15,
      ])
      cell_names.push("Cell" + cell_number.to_string())
      sample_ids.push(sample_name)
      cell_number = cell_number + 1
    }
    for local in 0..<right_count {
      coordinates.push([
        5.0 + (sample % 4).to_double() * 0.03 + local.to_double() * 0.08,
        local.to_double() * 0.15,
      ])
      cell_names.push("Cell" + cell_number.to_string())
      sample_ids.push(sample_name)
      cell_number = cell_number + 1
    }
  }
  MiloSampleData::{
    coordinates,
    cell_names,
    sample_ids,
    sample_names,
    sample_conditions,
  }
}
