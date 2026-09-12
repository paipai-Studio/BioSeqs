// Local indicators of spatial association inspired by Bioconductor lisaClust.
//
// Curves are generated independently per image. Each row describes one source
// cell across every target cell type and radius, then deterministic multi-start
// k-means identifies tissue regions with similar local spatial organisation.

///|
pub suberror LisaError {
  LisaError(String)
}

///|
pub enum LisaCurveKind {
  LisaStandardizedK
  LisaCenteredL
} derive(Eq, Debug)

///|
pub fn lisa_k_curve() -> LisaCurveKind {
  LisaStandardizedK
}

///|
pub fn lisa_l_curve() -> LisaCurveKind {
  LisaCenteredL
}

///|
pub enum LisaWindowKind {
  LisaRectangle
  LisaConvexHull
} derive(Eq, Debug)

///|
pub fn lisa_rectangle_window() -> LisaWindowKind {
  LisaRectangle
}

///|
pub fn lisa_convex_window() -> LisaWindowKind {
  LisaConvexHull
}

///|
fn lisa_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
pub struct LisaCell {
  cell_id : String
  image_id : String
  cell_type : String
  x : Double
  y : Double
} derive(Eq, Debug)

///|
pub fn LisaCell::create(
  cell_id : String,
  image_id : String,
  cell_type : String,
  x : Double,
  y : Double,
) -> LisaCell raise LisaError {
  if cell_id == "" {
    raise LisaError("lisaClust cell ID must not be empty")
  }
  if image_id == "" {
    raise LisaError("lisaClust image ID must not be empty")
  }
  if cell_type == "" {
    raise LisaError("lisaClust cell type must not be empty")
  }
  if !lisa_is_finite(x) || !lisa_is_finite(y) {
    raise LisaError("lisaClust coordinates must be finite")
  }
  LisaCell::{ cell_id, image_id, cell_type, x, y, }
}

///|
pub struct LisaConfig {
  radii : Array[Double]
  curve_kind : LisaCurveKind
  window_kind : LisaWindowKind
  bandwidth : Double
  min_density : Double
  edge_correction : Bool
  window_padding : Double
  edge_samples : Int
  n_clusters : Int
  n_starts : Int
  max_iterations : Int
  tolerance : Double
  seed : Int
  region_prefix : String
} derive(Debug)

///|
pub fn LisaConfig::create(
  radii? : Array[Double] = [20.0, 50.0, 100.0],
  curve_kind? : LisaCurveKind = LisaStandardizedK,
  window_kind? : LisaWindowKind = LisaConvexHull,
  bandwidth? : Double = 100000.0,
  min_density? : Double = 0.05,
  edge_correction? : Bool = true,
  window_padding? : Double = 0.01,
  edge_samples? : Int = 360,
  n_clusters? : Int = 2,
  n_starts? : Int = 10,
  max_iterations? : Int = 100,
  tolerance? : Double = 1.0e-8,
  seed? : Int = 51773,
  region_prefix? : String = "region",
) -> LisaConfig raise LisaError {
  if radii.length() == 0 {
    raise LisaError("lisaClust radii must not be empty")
  }
  let checked_radii = radii.copy()
  let mut previous = 0.0
  for index in 0..<checked_radii.length() {
    let radius = checked_radii[index]
    if !lisa_is_finite(radius) || radius <= 0.0 {
      raise LisaError("lisaClust radii must be finite and positive")
    }
    if index > 0 && radius <= previous {
      raise LisaError("lisaClust radii must be strictly increasing")
    }
    previous = radius
  }
  if !lisa_is_finite(bandwidth) || bandwidth <= 0.0 {
    raise LisaError("lisaClust bandwidth must be finite and positive")
  }
  if !lisa_is_finite(min_density) || min_density <= 0.0 {
    raise LisaError("lisaClust minimum density must be finite and positive")
  }
  if !lisa_is_finite(window_padding) || window_padding < 0.0 {
    raise LisaError("lisaClust window padding must be finite and non-negative")
  }
  if edge_samples < 32 {
    raise LisaError("lisaClust edge integration requires at least 32 samples")
  }
  if n_clusters < 1 {
    raise LisaError("lisaClust cluster count must be positive")
  }
  if n_starts < 1 || max_iterations < 1 {
    raise LisaError("lisaClust k-means controls must be positive")
  }
  if !lisa_is_finite(tolerance) || tolerance <= 0.0 || tolerance >= 1.0 {
    raise LisaError("lisaClust tolerance must be finite and in (0, 1)")
  }
  if seed < 0 {
    raise LisaError("lisaClust seed must be non-negative")
  }
  if region_prefix == "" {
    raise LisaError("lisaClust region prefix must not be empty")
  }
  LisaConfig::{
    radii: checked_radii,
    curve_kind,
    window_kind,
    bandwidth,
    min_density,
    edge_correction,
    window_padding,
    edge_samples,
    n_clusters,
    n_starts,
    max_iterations,
    tolerance,
    seed,
    region_prefix,
  }
}

///|
pub fn LisaConfig::default() -> LisaConfig {
  LisaConfig::{
    radii: [20.0, 50.0, 100.0],
    curve_kind: LisaStandardizedK,
    window_kind: LisaConvexHull,
    bandwidth: 100000.0,
    min_density: 0.05,
    edge_correction: true,
    window_padding: 0.01,
    edge_samples: 360,
    n_clusters: 2,
    n_starts: 10,
    max_iterations: 100,
    tolerance: 1.0e-8,
    seed: 51773,
    region_prefix: "region",
  }
}

///|
pub struct LisaImageWindow {
  image_id : String
  vertices : Array[(Double, Double)]
  area : Double
  max_radius : Double
} derive(Debug)

///|
pub struct LisaImageRadii {
  image_id : String
  radii : Array[Double]
} derive(Debug)

///|
pub struct LisaCurveSet {
  cell_ids : Array[String]
  image_ids : Array[String]
  cell_types : Array[String]
  target_cell_types : Array[String]
  radii : Array[Double]
  feature_names : Array[String]
  values : Array[Array[Double]]
  effective_radii : Array[LisaImageRadii]
  windows : Array[LisaImageWindow]
  curve_kind : LisaCurveKind
} derive(Debug)

///|
pub fn LisaCurveSet::n_cells(self : LisaCurveSet) -> Int {
  self.cell_ids.length()
}

///|
pub fn LisaCurveSet::n_features(self : LisaCurveSet) -> Int {
  self.feature_names.length()
}

///|
pub fn LisaCurveSet::feature_index(
  self : LisaCurveSet,
  target_cell_type : String,
  radius : Double,
) -> Int {
  let mut type_index = -1
  for index in 0..<self.target_cell_types.length() {
    if self.target_cell_types[index] == target_cell_type {
      type_index = index
      break
    }
  }
  if type_index < 0 {
    return -1
  }
  let mut radius_index = -1
  for index in 0..<self.radii.length() {
    if self.radii[index] == radius {
      radius_index = index
      break
    }
  }
  if radius_index < 0 {
    -1
  } else {
    type_index * self.radii.length() + radius_index
  }
}

///|
pub fn LisaCurveSet::curve_for_cell(
  self : LisaCurveSet,
  cell_id : String,
) -> Array[Double]? {
  for index in 0..<self.cell_ids.length() {
    if self.cell_ids[index] == cell_id {
      return Some(self.values[index].copy())
    }
  }
  None
}

///|
pub fn LisaCurveSet::effective_radii_for_image(
  self : LisaCurveSet,
  image_id : String,
) -> Array[Double]? {
  for entry in self.effective_radii {
    if entry.image_id == image_id {
      return Some(entry.radii.copy())
    }
  }
  None
}

///|
pub struct LisaRegionEnrichment {
  cell_type : String
  region : String
  observed : Int
  expected : Double
  relative_frequency : Double
} derive(Debug)

///|
pub struct LisaRegionSummary {
  region : String
  size : Int
  dominant_cell_type : String
  maximum_enrichment : Double
} derive(Debug)

///|
pub struct LisaResult {
  curves : LisaCurveSet
  labels : Array[Int]
  regions : Array[String]
  centroids : Array[Array[Double]]
  cluster_sizes : Array[Int]
  inertia : Double
  silhouette : Double
  iterations : Int
  converged : Bool
  enrichment : Array[LisaRegionEnrichment]
  region_summaries : Array[LisaRegionSummary]
  config : LisaConfig
} derive(Debug)

///|
pub fn LisaResult::n_cells(self : LisaResult) -> Int {
  self.curves.n_cells()
}

///|
pub fn LisaResult::n_regions(self : LisaResult) -> Int {
  self.cluster_sizes.length()
}

///|
pub fn LisaResult::region_for_cell(
  self : LisaResult,
  cell_id : String,
) -> String? {
  for index in 0..<self.curves.cell_ids.length() {
    if self.curves.cell_ids[index] == cell_id {
      return Some(self.regions[index])
    }
  }
  None
}

///|
pub fn LisaResult::region_enrichment(
  self : LisaResult,
  cell_type : String,
  region : String,
) -> LisaRegionEnrichment? {
  for entry in self.enrichment {
    if entry.cell_type == cell_type && entry.region == region {
      return Some(entry)
    }
  }
  None
}

///|
pub fn LisaResult::top_enrichments(
  self : LisaResult,
  limit? : Int = 10,
  minimum_relative_frequency? : Double = 0.0,
) -> Array[LisaRegionEnrichment] raise LisaError {
  if limit < 0 {
    raise LisaError("lisaClust enrichment limit must be non-negative")
  }
  if !lisa_is_finite(minimum_relative_frequency) ||
    minimum_relative_frequency < 0.0 {
    raise LisaError("lisaClust enrichment threshold must be non-negative")
  }
  let selected : Array[LisaRegionEnrichment] = []
  for entry in self.enrichment {
    if entry.relative_frequency >= minimum_relative_frequency {
      selected.push(entry)
    }
  }
  selected.sort_by(fn(
    left : LisaRegionEnrichment,
    right : LisaRegionEnrichment,
  ) -> Int {
    if left.relative_frequency > right.relative_frequency {
      -1
    } else if left.relative_frequency < right.relative_frequency {
      1
    } else if left.region < right.region {
      -1
    } else if left.region > right.region {
      1
    } else if left.cell_type < right.cell_type {
      -1
    } else if left.cell_type > right.cell_type {
      1
    } else {
      0
    }
  })
  let output : Array[LisaRegionEnrichment] = []
  for index in 0..<limit.min(selected.length()) {
    output.push(selected[index])
  }
  output
}

///|
pub fn LisaResult::summary(self : LisaResult) -> String {
  "lisaClust(cells=" +
  self.n_cells().to_string() +
  ", features=" +
  self.curves.n_features().to_string() +
  ", regions=" +
  self.n_regions().to_string() +
  ", inertia=" +
  self.inertia.to_string() +
  ", silhouette=" +
  self.silhouette.to_string() +
  ")"
}

///|
pub struct LisaSpatialExperimentOutput {
  experiment : SpatialExperiment
  result : LisaResult
}

///|
fn lisa_unique_strings(values : Array[String]) -> Array[String] {
  let output : Array[String] = []
  let seen : Map[String, Bool] = Map([])
  for value in values {
    if !seen.contains(value) {
      seen[value] = true
      output.push(value)
    }
  }
  output
}

///|
fn lisa_group_images(
  cells : Array[LisaCell],
) -> (Array[String], Map[String, Array[Int]]) {
  let image_ids : Array[String] = []
  let groups : Map[String, Array[Int]] = Map([])
  for index in 0..<cells.length() {
    let image_id = cells[index].image_id
    match groups.get(image_id) {
      Some(indices) => indices.push(index)
      None => {
        image_ids.push(image_id)
        groups[image_id] = [index]
      }
    }
  }
  (image_ids, groups)
}

///|
fn lisa_validate_cells(cells : Array[LisaCell]) -> Unit raise LisaError {
  if cells.length() < 2 {
    raise LisaError("lisaClust requires at least two cells")
  }
  let identifiers : Map[String, Bool] = Map([])
  for cell in cells {
    if identifiers.contains(cell.cell_id) {
      raise LisaError("lisaClust cell IDs must be unique")
    }
    identifiers[cell.cell_id] = true
    if cell.image_id == "" || cell.cell_type == "" || cell.cell_id == "" {
      raise LisaError("lisaClust cell identifiers must not be empty")
    }
    if !lisa_is_finite(cell.x) || !lisa_is_finite(cell.y) {
      raise LisaError("lisaClust coordinates must be finite")
    }
  }
}

///|
fn lisa_cross(
  origin : (Double, Double),
  first : (Double, Double),
  second : (Double, Double),
) -> Double {
  (first.0 - origin.0) * (second.1 - origin.1) -
  (first.1 - origin.1) * (second.0 - origin.0)
}

///|
fn lisa_convex_hull(
  points : Array[(Double, Double)],
) -> Array[(Double, Double)] {
  let sorted = points.copy()
  sorted.sort_by(fn(left : (Double, Double), right : (Double, Double)) -> Int {
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
  let unique : Array[(Double, Double)] = []
  for point in sorted {
    if unique.length() == 0 || unique[unique.length() - 1] != point {
      unique.push(point)
    }
  }
  if unique.length() < 3 {
    return unique
  }
  let lower : Array[(Double, Double)] = []
  for point in unique {
    while lower.length() >= 2 &&
          lisa_cross(
            lower[lower.length() - 2],
            lower[lower.length() - 1],
            point,
          ) <=
          0.0 {
      ignore(lower.pop())
    }
    lower.push(point)
  }
  let upper : Array[(Double, Double)] = []
  let mut index = unique.length() - 1
  while index >= 0 {
    let point = unique[index]
    while upper.length() >= 2 &&
          lisa_cross(
            upper[upper.length() - 2],
            upper[upper.length() - 1],
            point,
          ) <=
          0.0 {
      ignore(upper.pop())
    }
    upper.push(point)
    index = index - 1
  }
  ignore(lower.pop())
  ignore(upper.pop())
  for point in upper {
    lower.push(point)
  }
  lower
}

///|
fn lisa_polygon_area(vertices : Array[(Double, Double)]) -> Double {
  let mut twice_area = 0.0
  for index in 0..<vertices.length() {
    let next = (index + 1) % vertices.length()
    twice_area = twice_area +
      vertices[index].0 * vertices[next].1 -
      vertices[next].0 * vertices[index].1
  }
  twice_area.abs() * 0.5
}

///|
fn lisa_expand_polygon(
  vertices : Array[(Double, Double)],
  padding : Double,
) -> Array[(Double, Double)] {
  if padding <= 0.0 {
    return vertices.copy()
  }
  let mut center_x = 0.0
  let mut center_y = 0.0
  for vertex in vertices {
    center_x = center_x + vertex.0
    center_y = center_y + vertex.1
  }
  center_x = center_x / vertices.length().to_double()
  center_y = center_y / vertices.length().to_double()
  let output : Array[(Double, Double)] = []
  for vertex in vertices {
    let dx = vertex.0 - center_x
    let dy = vertex.1 - center_y
    let distance = (dx * dx + dy * dy).sqrt()
    if distance > 1.0e-14 {
      let scale = (distance + padding) / distance
      output.push((center_x + dx * scale, center_y + dy * scale))
    } else {
      output.push(vertex)
    }
  }
  output
}

///|
fn lisa_image_window(
  image_id : String,
  cells : Array[LisaCell],
  indices : Array[Int],
  config : LisaConfig,
) -> LisaImageWindow raise LisaError {
  if indices.length() < 2 {
    raise LisaError("lisaClust images must contain at least two cells")
  }
  let mut min_x = cells[indices[0]].x
  let mut max_x = min_x
  let mut min_y = cells[indices[0]].y
  let mut max_y = min_y
  let points : Array[(Double, Double)] = []
  for index in indices {
    let cell = cells[index]
    min_x = min_x.min(cell.x)
    max_x = max_x.max(cell.x)
    min_y = min_y.min(cell.y)
    max_y = max_y.max(cell.y)
    points.push((cell.x, cell.y))
  }
  let width = max_x - min_x
  let height = max_y - min_y
  if width <= 1.0e-14 || height <= 1.0e-14 {
    raise LisaError(
      "lisaClust image windows must have positive width and height",
    )
  }
  let vertices = match config.window_kind {
    LisaRectangle =>
      [
        (min_x - config.window_padding, min_y - config.window_padding),
        (max_x + config.window_padding, min_y - config.window_padding),
        (max_x + config.window_padding, max_y + config.window_padding),
        (min_x - config.window_padding, max_y + config.window_padding),
      ]
    LisaConvexHull => {
      let hull = lisa_convex_hull(points)
      if hull.length() < 3 {
        raise LisaError("lisaClust convex windows require non-collinear cells")
      }
      lisa_expand_polygon(hull, config.window_padding)
    }
  }
  let area = lisa_polygon_area(vertices)
  if !lisa_is_finite(area) || area <= 1.0e-14 {
    raise LisaError("lisaClust image window area must be positive")
  }
  let mut vertex_min_x = vertices[0].0
  let mut vertex_max_x = vertices[0].0
  let mut vertex_min_y = vertices[0].1
  let mut vertex_max_y = vertices[0].1
  for vertex in vertices {
    vertex_min_x = vertex_min_x.min(vertex.0)
    vertex_max_x = vertex_max_x.max(vertex.0)
    vertex_min_y = vertex_min_y.min(vertex.1)
    vertex_max_y = vertex_max_y.max(vertex.1)
  }
  let max_radius = (vertex_max_x - vertex_min_x).min(
      vertex_max_y - vertex_min_y,
    ) /
    2.01
  LisaImageWindow::{ image_id, vertices, area, max_radius, }
}

///|
fn lisa_ray_boundary_distance(
  x : Double,
  y : Double,
  dx : Double,
  dy : Double,
  vertices : Array[(Double, Double)],
) -> Double {
  let mut closest = 1.0e300
  for index in 0..<vertices.length() {
    let next = (index + 1) % vertices.length()
    let ax = vertices[index].0
    let ay = vertices[index].1
    let ex = vertices[next].0 - ax
    let ey = vertices[next].1 - ay
    let denominator = dx * ey - dy * ex
    if denominator.abs() <= 1.0e-14 {
      continue
    }
    let qx = ax - x
    let qy = ay - y
    let distance = (qx * ey - qy * ex) / denominator
    let segment = (qx * dy - qy * dx) / denominator
    if distance >= 0.0 && segment >= -1.0e-10 && segment <= 1.0 + 1.0e-10 {
      closest = closest.min(distance)
    }
  }
  closest
}

///|
fn lisa_visible_disc_fraction(
  cell : LisaCell,
  radius : Double,
  window : LisaImageWindow,
  samples : Int,
) -> Double {
  let angle_step = 2.0 * @math.PI / samples.to_double()
  let mut area = 0.0
  for sample in 0..<samples {
    let angle = (sample.to_double() + 0.5) * angle_step
    let boundary = lisa_ray_boundary_distance(
      cell.x,
      cell.y,
      @math.cos(angle),
      @math.sin(angle),
      window.vertices,
    )
    let visible_radius = radius.min(boundary)
    area = area + 0.5 * visible_radius * visible_radius * angle_step
  }
  (area / (@math.PI * radius * radius)).max(1.0e-8).min(1.0)
}

///|
fn lisa_density_weights(
  cells : Array[LisaCell],
  indices : Array[Int],
  bandwidth : Double,
  min_density : Double,
) -> Array[Double] {
  let density = Array::make(indices.length(), 0.0)
  let denominator = 2.0 * bandwidth * bandwidth
  for target in 0..<indices.length() {
    let target_cell = cells[indices[target]]
    let mut value = 0.0
    for source in 0..<indices.length() {
      let source_cell = cells[indices[source]]
      let dx = target_cell.x - source_cell.x
      let dy = target_cell.y - source_cell.y
      value = value + @math.exp(-(dx * dx + dy * dy) / denominator)
    }
    density[target] = value
  }
  let mut mean = 0.0
  for value in density {
    mean = mean + value
  }
  mean = mean / density.length().to_double()
  for index in 0..<density.length() {
    density[index] = (density[index] / mean).max(min_density)
  }
  let mut clipped_mean = 0.0
  for value in density {
    clipped_mean = clipped_mean + value
  }
  clipped_mean = clipped_mean / density.length().to_double()
  let weights = Array::make(density.length(), 1.0)
  for index in 0..<density.length() {
    weights[index] = clipped_mean / density[index]
  }
  weights
}

///|
fn lisa_feature_names(
  target_cell_types : Array[String],
  radii : Array[Double],
) -> Array[String] {
  let output : Array[String] = []
  for cell_type in target_cell_types {
    for radius in radii {
      output.push(radius.to_string() + "_" + cell_type)
    }
  }
  output
}

///|
pub fn lisa_curves(
  cells : Array[LisaCell],
  config? : LisaConfig = LisaConfig::default(),
) -> LisaCurveSet raise LisaError {
  lisa_validate_cells(cells)
  let all_types : Array[String] = []
  for cell in cells {
    all_types.push(cell.cell_type)
  }
  let target_cell_types = lisa_unique_strings(all_types)
  let (image_order, groups) = lisa_group_images(cells)
  let values : Array[Array[Double]] = []
  for _ in cells {
    values.push(
      Array::make(target_cell_types.length() * config.radii.length(), 0.0),
    )
  }
  let windows : Array[LisaImageWindow] = []
  let effective_radii : Array[LisaImageRadii] = []
  for image_id in image_order {
    let indices = groups[image_id]
    let window = lisa_image_window(image_id, cells, indices, config)
    windows.push(window)
    let image_radii : Array[Double] = []
    for radius in config.radii {
      image_radii.push(radius.min(window.max_radius))
    }
    effective_radii.push(LisaImageRadii::{
      image_id,
      radii: image_radii.copy(),
    })
    let density_weights = lisa_density_weights(
      cells,
      indices,
      config.bandwidth,
      config.min_density,
    )
    let local_by_global : Map[Int, Int] = Map([])
    for local in 0..<indices.length() {
      local_by_global[indices[local]] = local
    }
    for source_local in 0..<indices.length() {
      let source_global = indices[source_local]
      for type_index in 0..<target_cell_types.length() {
        let target_type = target_cell_types[type_index]
        let mut target_count = 0
        for target_global in indices {
          if cells[target_global].cell_type == target_type {
            target_count = target_count + 1
          }
        }
        let intensity = target_count.to_double() / window.area
        for radius_index in 0..<image_radii.length() {
          let radius = image_radii[radius_index]
          let mut observed = 0.0
          for target_global in indices {
            if target_global == source_global ||
              cells[target_global].cell_type != target_type {
              continue
            }
            let dx = cells[target_global].x - cells[source_global].x
            let dy = cells[target_global].y - cells[source_global].y
            if dx * dx + dy * dy <= radius * radius + 1.0e-12 {
              observed = observed +
                density_weights[local_by_global[target_global]]
            }
          }
          let edge_fraction = if config.edge_correction {
            lisa_visible_disc_fraction(
              cells[source_global],
              radius,
              window,
              config.edge_samples,
            )
          } else {
            1.0
          }
          let expected = @math.PI * radius * radius * edge_fraction * intensity
          let feature = type_index * config.radii.length() + radius_index
          values[source_global][feature] = if expected <= 1.0e-14 {
            0.0
          } else {
            match config.curve_kind {
              LisaStandardizedK => (observed - expected) / expected.sqrt()
              LisaCenteredL => observed.sqrt() - expected.sqrt()
            }
          }
        }
      }
    }
  }
  let cell_ids : Array[String] = []
  let image_ids : Array[String] = []
  let cell_types : Array[String] = []
  for cell in cells {
    cell_ids.push(cell.cell_id)
    image_ids.push(cell.image_id)
    cell_types.push(cell.cell_type)
  }
  LisaCurveSet::{
    cell_ids,
    image_ids,
    cell_types,
    target_cell_types,
    radii: config.radii.copy(),
    feature_names: lisa_feature_names(target_cell_types, config.radii),
    values,
    effective_radii,
    windows,
    curve_kind: config.curve_kind,
  }
}

///|
fn lisa_squared_distance(left : Array[Double], right : Array[Double]) -> Double {
  let mut output = 0.0
  for index in 0..<left.length() {
    let difference = left[index] - right[index]
    output = output + difference * difference
  }
  output
}

///|
fn lisa_initialize_centroids(
  data : Array[Array[Double]],
  clusters : Int,
  first : Int,
) -> Array[Array[Double]] {
  let centroids : Array[Array[Double]] = [data[first].copy()]
  let selected = Array::make(data.length(), false)
  selected[first] = true
  while centroids.length() < clusters {
    let mut best_index = -1
    let mut best_distance = -1.0
    for point in 0..<data.length() {
      if selected[point] {
        continue
      }
      let mut nearest = lisa_squared_distance(data[point], centroids[0])
      for center in 1..<centroids.length() {
        nearest = nearest.min(
          lisa_squared_distance(data[point], centroids[center]),
        )
      }
      if nearest > best_distance {
        best_distance = nearest
        best_index = point
      }
    }
    if best_index < 0 {
      break
    }
    selected[best_index] = true
    centroids.push(data[best_index].copy())
  }
  centroids
}

///|
fn lisa_assign_clusters(
  data : Array[Array[Double]],
  centroids : Array[Array[Double]],
) -> (Array[Int], Double) {
  let labels = Array::make(data.length(), 0)
  let mut inertia = 0.0
  for point in 0..<data.length() {
    let mut best = 0
    let mut distance = lisa_squared_distance(data[point], centroids[0])
    for cluster in 1..<centroids.length() {
      let candidate = lisa_squared_distance(data[point], centroids[cluster])
      if candidate < distance {
        best = cluster
        distance = candidate
      }
    }
    labels[point] = best
    inertia = inertia + distance
  }
  (labels, inertia)
}

///|
fn lisa_update_centroids(
  data : Array[Array[Double]],
  labels : Array[Int],
  old_centroids : Array[Array[Double]],
) -> Array[Array[Double]] {
  let clusters = old_centroids.length()
  let dimensions = data[0].length()
  let output : Array[Array[Double]] = []
  for _ in 0..<clusters {
    output.push(Array::make(dimensions, 0.0))
  }
  let counts = Array::make(clusters, 0)
  for point in 0..<data.length() {
    let cluster = labels[point]
    counts[cluster] = counts[cluster] + 1
    for dimension in 0..<dimensions {
      output[cluster][dimension] = output[cluster][dimension] +
        data[point][dimension]
    }
  }
  for cluster in 0..<clusters {
    if counts[cluster] > 0 {
      for dimension in 0..<dimensions {
        output[cluster][dimension] = output[cluster][dimension] /
          counts[cluster].to_double()
      }
    } else {
      let mut farthest = 0
      let mut farthest_distance = -1.0
      for point in 0..<data.length() {
        let distance = lisa_squared_distance(
          data[point],
          old_centroids[labels[point]],
        )
        if distance > farthest_distance {
          farthest = point
          farthest_distance = distance
        }
      }
      output[cluster] = data[farthest].copy()
    }
  }
  output
}

///|
struct LisaKmeansFit {
  labels : Array[Int]
  centroids : Array[Array[Double]]
  inertia : Double
  iterations : Int
  converged : Bool
}

///|
fn lisa_kmeans_once(
  data : Array[Array[Double]],
  config : LisaConfig,
  start : Int,
) -> LisaKmeansFit {
  let first = (config.seed + start * 9973) % data.length()
  let mut centroids = lisa_initialize_centroids(data, config.n_clusters, first)
  let mut labels = Array::make(data.length(), -1)
  let mut iterations = 0
  let mut converged = false
  for iteration in 0..<config.max_iterations {
    let (next_labels, _) = lisa_assign_clusters(data, centroids)
    let next_centroids = lisa_update_centroids(data, next_labels, centroids)
    let mut maximum_shift = 0.0
    for cluster in 0..<centroids.length() {
      maximum_shift = maximum_shift.max(
        lisa_squared_distance(centroids[cluster], next_centroids[cluster]).sqrt(),
      )
    }
    let mut changed = false
    for point in 0..<labels.length() {
      if labels[point] != next_labels[point] {
        changed = true
      }
    }
    labels = next_labels
    centroids = next_centroids
    iterations = iteration + 1
    if !changed || maximum_shift < config.tolerance {
      converged = true
      break
    }
  }
  let (final_labels, inertia) = lisa_assign_clusters(data, centroids)
  LisaKmeansFit::{
    labels: final_labels,
    centroids,
    inertia,
    iterations,
    converged,
  }
}

///|
fn lisa_canonicalize_fit(fit : LisaKmeansFit) -> LisaKmeansFit {
  let mapping = Array::make(fit.centroids.length(), -1)
  let order : Array[Int] = []
  for label in fit.labels {
    if mapping[label] < 0 {
      mapping[label] = order.length()
      order.push(label)
    }
  }
  for label in 0..<fit.centroids.length() {
    if mapping[label] < 0 {
      mapping[label] = order.length()
      order.push(label)
    }
  }
  let labels = fit.labels.copy()
  for index in 0..<labels.length() {
    labels[index] = mapping[labels[index]]
  }
  let centroids : Array[Array[Double]] = []
  for old_label in order {
    centroids.push(fit.centroids[old_label].copy())
  }
  LisaKmeansFit::{ ..fit, labels, centroids, }
}

///|
fn lisa_silhouette(
  data : Array[Array[Double]],
  labels : Array[Int],
  clusters : Int,
) -> Double {
  if clusters <= 1 || data.length() <= 1 {
    return 0.0
  }
  let mut total = 0.0
  for point in 0..<data.length() {
    let own = labels[point]
    let sums = Array::make(clusters, 0.0)
    let counts = Array::make(clusters, 0)
    for other in 0..<data.length() {
      if other == point {
        continue
      }
      let cluster = labels[other]
      sums[cluster] = sums[cluster] +
        lisa_squared_distance(data[point], data[other]).sqrt()
      counts[cluster] = counts[cluster] + 1
    }
    if counts[own] == 0 {
      continue
    }
    let within = sums[own] / counts[own].to_double()
    let mut nearest = 1.0e300
    for cluster in 0..<clusters {
      if cluster != own && counts[cluster] > 0 {
        nearest = nearest.min(sums[cluster] / counts[cluster].to_double())
      }
    }
    if nearest < 1.0e299 {
      let denominator = within.max(nearest)
      if denominator > 0.0 {
        total = total + (nearest - within) / denominator
      }
    }
  }
  total / data.length().to_double()
}

///|
fn lisa_enrichment(
  curves : LisaCurveSet,
  labels : Array[Int],
  config : LisaConfig,
) -> (Array[LisaRegionEnrichment], Array[LisaRegionSummary]) {
  let cell_types = lisa_unique_strings(curves.cell_types)
  let counts : Array[Array[Int]] = []
  for _ in cell_types {
    counts.push(Array::make(config.n_clusters, 0))
  }
  let type_totals = Array::make(cell_types.length(), 0)
  let region_totals = Array::make(config.n_clusters, 0)
  for cell in 0..<curves.n_cells() {
    let mut type_index = 0
    for index in 0..<cell_types.length() {
      if cell_types[index] == curves.cell_types[cell] {
        type_index = index
        break
      }
    }
    let region = labels[cell]
    counts[type_index][region] = counts[type_index][region] + 1
    type_totals[type_index] = type_totals[type_index] + 1
    region_totals[region] = region_totals[region] + 1
  }
  let enrichment : Array[LisaRegionEnrichment] = []
  for type_index in 0..<cell_types.length() {
    for region in 0..<config.n_clusters {
      let expected = type_totals[type_index].to_double() *
        region_totals[region].to_double() /
        curves.n_cells().to_double()
      let relative_frequency = if expected > 0.0 {
        counts[type_index][region].to_double() / expected
      } else {
        0.0
      }
      enrichment.push(LisaRegionEnrichment::{
        cell_type: cell_types[type_index],
        region: config.region_prefix + "_" + (region + 1).to_string(),
        observed: counts[type_index][region],
        expected,
        relative_frequency,
      })
    }
  }
  let summaries : Array[LisaRegionSummary] = []
  for region in 0..<config.n_clusters {
    let mut dominant = ""
    let mut maximum = -1.0
    for entry in enrichment {
      if entry.region == config.region_prefix + "_" + (region + 1).to_string() &&
        entry.relative_frequency > maximum {
        dominant = entry.cell_type
        maximum = entry.relative_frequency
      }
    }
    summaries.push(LisaRegionSummary::{
      region: config.region_prefix + "_" + (region + 1).to_string(),
      size: region_totals[region],
      dominant_cell_type: dominant,
      maximum_enrichment: maximum.max(0.0),
    })
  }
  (enrichment, summaries)
}

///|
pub fn lisaclust_from_curves(
  curves : LisaCurveSet,
  config? : LisaConfig = LisaConfig::default(),
) -> LisaResult raise LisaError {
  if curves.n_cells() == 0 || curves.n_features() == 0 {
    raise LisaError("lisaClust clustering requires non-empty curves")
  }
  if curves.values.length() != curves.n_cells() {
    raise LisaError("lisaClust curve rows must match cells")
  }
  if config.n_clusters > curves.n_cells() {
    raise LisaError("lisaClust cluster count must not exceed the cell count")
  }
  for row in curves.values {
    if row.length() != curves.n_features() {
      raise LisaError("lisaClust curve matrix must be rectangular")
    }
    for value in row {
      if !lisa_is_finite(value) {
        raise LisaError("lisaClust curves must be finite")
      }
    }
  }
  let mut best : LisaKmeansFit? = None
  for start in 0..<config.n_starts {
    let fit = lisa_kmeans_once(curves.values, config, start)
    best = match best {
      None => Some(fit)
      Some(current) =>
        if fit.inertia < current.inertia {
          Some(fit)
        } else {
          Some(current)
        }
    }
  }
  let fit = lisa_canonicalize_fit(
    match best {
      Some(value) => value
      None => raise LisaError("lisaClust k-means did not produce a fit")
    },
  )
  let cluster_sizes = Array::make(config.n_clusters, 0)
  let regions : Array[String] = []
  for label in fit.labels {
    cluster_sizes[label] = cluster_sizes[label] + 1
    regions.push(config.region_prefix + "_" + (label + 1).to_string())
  }
  let silhouette = lisa_silhouette(curves.values, fit.labels, config.n_clusters)
  let (enrichment, region_summaries) = lisa_enrichment(
    curves,
    fit.labels,
    config,
  )
  LisaResult::{
    curves,
    labels: fit.labels,
    regions,
    centroids: fit.centroids,
    cluster_sizes,
    inertia: fit.inertia,
    silhouette,
    iterations: fit.iterations,
    converged: fit.converged,
    enrichment,
    region_summaries,
    config,
  }
}

///|
pub fn lisaclust(
  cells : Array[LisaCell],
  config? : LisaConfig = LisaConfig::default(),
) -> LisaResult raise LisaError {
  let curves = lisa_curves(cells, config~)
  lisaclust_from_curves(curves, config~)
}

///|
fn lisa_copy_string_map(input : Map[String, String]) -> Map[String, String] {
  let output : Map[String, String] = Map([], capacity=input.length())
  for key in input.keys() {
    output[key] = input[key]
  }
  output
}

///|
fn lisa_copy_spatial_experiment(
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
    output.row_data.push(lisa_copy_string_map(row))
  }
  for column in experiment.col_data {
    output.col_data.push(lisa_copy_string_map(column))
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
fn lisa_required_column(
  column : Map[String, String],
  key : String,
  row : Int,
) -> String raise LisaError {
  match column.get(key) {
    Some(value) =>
      if value == "" {
        raise LisaError(
          "lisaClust SpatialExperiment column '" +
          key +
          "' is empty at cell " +
          row.to_string(),
        )
      } else {
        value
      }
    None =>
      raise LisaError(
        "lisaClust SpatialExperiment is missing column '" + key + "'",
      )
  }
}

///|
pub fn lisaclust_spatial_experiment(
  experiment : SpatialExperiment,
  image_id_key? : String = "imageID",
  cell_type_key? : String = "cellType",
  cell_id_key? : String = "",
  region_key? : String = "region",
  config? : LisaConfig = LisaConfig::default(),
) -> LisaSpatialExperimentOutput raise LisaError {
  if experiment.col_data.length() == 0 {
    raise LisaError("lisaClust SpatialExperiment col_data must not be empty")
  }
  if experiment.col_data.length() != experiment.spatial_coords.length() {
    raise LisaError(
      "lisaClust SpatialExperiment coordinates must match col_data rows",
    )
  }
  if image_id_key == "" || cell_type_key == "" || region_key == "" {
    raise LisaError("lisaClust SpatialExperiment key names must not be empty")
  }
  let cells : Array[LisaCell] = []
  for row in 0..<experiment.col_data.length() {
    let column = experiment.col_data[row]
    let image_id = lisa_required_column(column, image_id_key, row)
    let cell_type = lisa_required_column(column, cell_type_key, row)
    let cell_id = if cell_id_key == "" {
      "cell_" + (row + 1).to_string()
    } else {
      lisa_required_column(column, cell_id_key, row)
    }
    cells.push(
      LisaCell::create(
        cell_id,
        image_id,
        cell_type,
        experiment.spatial_coords[row].x,
        experiment.spatial_coords[row].y,
      ),
    )
  }
  let result = lisaclust(cells, config~)
  let enriched = lisa_copy_spatial_experiment(experiment)
  for row in 0..<enriched.col_data.length() {
    enriched.col_data[row][region_key] = result.regions[row]
  }
  enriched.metadata["lisaclust_cells"] = result.n_cells().to_string()
  enriched.metadata["lisaclust_features"] = result.curves
    .n_features()
    .to_string()
  enriched.metadata["lisaclust_regions"] = result.n_regions().to_string()
  enriched.metadata["lisaclust_silhouette"] = result.silhouette.to_string()
  LisaSpatialExperimentOutput::{ experiment: enriched, result, }
}

///|
pub fn lisaclust_example_data() -> Array[LisaCell] {
  let cells : Array[LisaCell] = []
  for image in 0..<2 {
    for column in 0..<8 {
      for row in 0..<4 {
        let left = column < 4
        let cell_type = if left {
          if (column + row + image) % 3 == 0 {
            "B_cell"
          } else {
            "T_cell"
          }
        } else if (column + row + image) % 3 == 0 {
          "Stroma"
        } else {
          "Tumour"
        }
        let cell_id = "img" +
          (image + 1).to_string() +
          "_c" +
          column.to_string() +
          "_r" +
          row.to_string()
        let cell = LisaCell::create(
          cell_id,
          "image_" + (image + 1).to_string(),
          cell_type,
          column.to_double() * 2.0 + image.to_double() * 0.1,
          row.to_double() * 2.0 + (column % 2).to_double() * 0.15,
        ) catch {
          _ => abort("lisaClust example cell should be valid")
        }
        cells.push(cell)
      }
    }
  }
  cells
}
