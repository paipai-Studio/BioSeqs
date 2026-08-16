// Tests for the Bioconductor spicyR-inspired spatial association model.

///|
fn spicyr_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "spicyR value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn spicyr_test_config(
  edge_correction : Bool,
  include_zero_cells : Bool,
) -> @src.SpicyConfig {
  @src.SpicyConfig::create(
    radii=[1.0],
    edge_correction~,
    include_zero_cells~,
    window_padding=0.0,
    tolerance=0.01,
  ) catch {
    _ => abort("spicyR test configuration should be valid")
  }
}

///|
fn spicyr_test_cell(
  image_id : String,
  cell_type : String,
  x : Double,
  y : Double,
) -> @src.SpicyCell {
  @src.SpicyCell::create(image_id, cell_type, x, y) catch {
    _ => abort("spicyR test cell should be valid")
  }
}

///|
fn spicyr_test_image_metadata(
  image_id : String,
  condition : String,
  subject : String,
  covariates : Map[String, Double],
) -> @src.SpicyImageMetadata {
  @src.SpicyImageMetadata::create(image_id, condition, subject~, covariates~) catch {
    _ => abort("spicyR test metadata should be valid")
  }
}

///|
fn spicyr_test_pair(
  from_cell_type : String,
  to_cell_type : String,
) -> @src.SpicyPair {
  @src.SpicyPair::create(from_cell_type, to_cell_type) catch {
    _ => abort("spicyR test pair should be valid")
  }
}

///|
fn spicyr_test_metadata(
  count : Int,
  split : Int,
) -> Array[@src.SpicyImageMetadata] {
  let metadata : Array[@src.SpicyImageMetadata] = []
  for index in 0..<count {
    metadata.push(
      spicyr_test_image_metadata(
        "image_" + (index + 1).to_string(),
        if index < split {
          "control"
        } else {
          "treated"
        },
        "",
        Map([("age", 30.0 + (index % split).to_double() * 5.0)]),
      ),
    )
  }
  metadata
}

///|
fn spicyr_test_geometry_cells(
  metadata : Array[@src.SpicyImageMetadata],
  source_x : Double,
  source_y : Double,
  target_x : Double,
  target_y : Double,
) -> Array[@src.SpicyCell] {
  let cells : Array[@src.SpicyCell] = []
  for image in metadata {
    cells.push(spicyr_test_cell(image.image_id, "A", source_x, source_y))
    cells.push(spicyr_test_cell(image.image_id, "B", target_x, target_y))
    cells.push(spicyr_test_cell(image.image_id, "Frame", 0.0, 0.0))
    cells.push(spicyr_test_cell(image.image_id, "Frame", 10.0, 10.0))
  }
  cells
}

///|
fn spicyr_test_association(
  image_id : String,
  from_type : String,
  to_type : String,
  statistic : Double,
  count : Int,
) -> @src.SpicyAssociation {
  @src.SpicyAssociation::create(
    image_id,
    from_type,
    to_type,
    statistic,
    from_count=count,
    to_count=count,
  ) catch {
    _ => abort("spicyR test association should be valid")
  }
}

///|
fn spicyr_test_linear_input() -> (
  Array[@src.SpicyImageMetadata],
  Array[@src.SpicyAssociation],
) {
  let metadata = spicyr_test_metadata(8, 4)
  let associations : Array[@src.SpicyAssociation] = []
  let noise = [-0.2, 0.1, 0.0, 0.15, -0.1, 0.2, -0.05, 0.1]
  for index in 0..<metadata.length() {
    let treatment = if index < 4 { 0.0 } else { 2.0 }
    associations.push(
      spicyr_test_association(
        metadata[index].image_id,
        "A",
        "B",
        treatment + noise[index],
        index + 2,
      ),
    )
  }
  (metadata, associations)
}

///|
fn spicyr_test_random_input() -> (
  Array[@src.SpicyImageMetadata],
  Array[@src.SpicyAssociation],
) {
  let metadata : Array[@src.SpicyImageMetadata] = []
  let associations : Array[@src.SpicyAssociation] = []
  let subject_effects = [-1.0, -0.25, 0.5, 1.25, 2.0]
  for subject_index in 0..<subject_effects.length() {
    for condition_index in 0..<2 {
      let condition = if condition_index == 0 { "control" } else { "treated" }
      let image_id = "subject_" +
        (subject_index + 1).to_string() +
        "_" +
        condition
      metadata.push(
        spicyr_test_image_metadata(
          image_id,
          condition,
          "subject_" + (subject_index + 1).to_string(),
          Map([]),
        ),
      )
      let noise = if subject_index % 2 == 0 {
        condition_index.to_double() * 0.05
      } else {
        -condition_index.to_double() * 0.04
      }
      associations.push(
        spicyr_test_association(
          image_id,
          "A",
          "B",
          subject_effects[subject_index] +
          condition_index.to_double() * 1.5 +
          noise,
          8,
        ),
      )
    }
  }
  (metadata, associations)
}

///|
fn spicyr_test_two_pair_result() -> @src.SpicyResult {
  let (metadata, first) = spicyr_test_linear_input()
  let associations = first.copy()
  let weak = [0.1, -0.2, 0.05, 0.0, -0.1, 0.2, -0.05, 0.1]
  for index in 0..<metadata.length() {
    associations.push(
      spicyr_test_association(
        metadata[index].image_id,
        "C",
        "D",
        weak[index],
        5,
      ),
    )
  }
  @src.spicy_fit_associations(associations, metadata) catch {
    _ => abort("spicyR two-pair fit should succeed")
  }
}

///|
fn spicyr_test_spatial_experiment() -> @src.SpatialExperiment {
  let experiment = @src.SpatialExperiment::new()
  for image_index in 0..<6 {
    let image_id = "image_" + (image_index + 1).to_string()
    let condition = if image_index < 3 { "control" } else { "treated" }
    let target_x = if condition == "control" { 9.0 } else { 5.5 }
    let cell_types = ["A", "B", "Frame", "Frame"]
    let xs = [5.0, target_x, 0.0, 10.0]
    let ys = [5.0, if condition == "control" { 9.0 } else { 5.0 }, 0.0, 10.0]
    for cell_index in 0..<cell_types.length() {
      ignore(
        @src.se_add_col(
          experiment,
          Map([
            ("imageID", image_id),
            ("cellType", cell_types[cell_index]),
            ("condition", condition),
            ("subject", "subject_" + (image_index / 2 + 1).to_string()),
            ("age", (40 + image_index).to_string()),
          ]),
        ),
      )
      ignore(
        @src.se_add_spatial_coord(
          experiment,
          @src.SpatialCoord::new_2d(xs[cell_index], ys[cell_index]),
        ),
      )
    }
  }
  experiment
}

///|
test "spicyR: default configuration exposes upstream controls" {
  let config = @src.SpicyConfig::default()
  assert_eq(config.radii, [20.0, 50.0, 100.0])
  assert_true(config.edge_correction)
  assert_true(!config.include_zero_cells)
  assert_true(config.use_weights)
  assert_eq(config.weight_factor, 1.0)
  assert_eq(config.fdr_threshold, 0.05)
  assert_eq(config.max_iterations, 100)
}

///|
test "spicyR: custom configuration copies radii" {
  let radii = [1.0, 2.0, 3.0]
  let config = @src.SpicyConfig::create(
    radii~,
    edge_correction=false,
    include_zero_cells=true,
    use_weights=false,
    weight_factor=0.5,
    window_padding=2.0,
    fdr_threshold=0.2,
    reference_condition="control",
  ) catch {
    _ => abort("custom spicyR configuration should be valid")
  }
  radii[0] = 9.0
  assert_eq(config.radii[0], 1.0)
  assert_true(!config.edge_correction)
  assert_true(config.include_zero_cells)
  assert_true(!config.use_weights)
  assert_eq(config.reference_condition, "control")
}

///|
test "spicyR: constructors preserve typed values" {
  let cell = @src.SpicyCell::create("image", "T_cell", 1.5, 2.5) catch {
    _ => abort("spicyR cell should be valid")
  }
  let pair = @src.SpicyPair::create("T_cell", "Tumour") catch {
    _ => abort("spicyR pair should be valid")
  }
  let association = @src.SpicyAssociation::create(
    "image",
    "T_cell",
    "Tumour",
    1.25,
    from_count=4,
    to_count=6,
  ) catch {
    _ => abort("spicyR association should be valid")
  }
  assert_eq(cell.image_id, "image")
  assert_eq(cell.x, 1.5)
  assert_eq(pair.label(), "T_cell__Tumour")
  assert_eq(association.statistic, 1.25)
  assert_eq(association.from_count, 4)
}

///|
test "spicyR: metadata constructor copies covariates" {
  let covariates : Map[String, Double] = Map([("age", 52.0)])
  let metadata = @src.SpicyImageMetadata::create(
    "image",
    "control",
    subject="subject_1",
    covariates~,
  ) catch {
    _ => abort("spicyR metadata should be valid")
  }
  covariates["age"] = 90.0
  assert_eq(metadata.subject, "subject_1")
  assert_eq(metadata.covariates["age"], 52.0)
}

///|
test "spicyR: constructors reject empty and non-finite values" {
  let bad_cell = try {
    ignore(@src.SpicyCell::create("", "A", 0.0, 0.0))
    false
  } catch {
    SpicyError(_) => true
  }
  let bad_coordinate = try {
    ignore(@src.SpicyCell::create("image", "A", @double.not_a_number, 0.0))
    false
  } catch {
    SpicyError(_) => true
  }
  let bad_pair = try {
    ignore(@src.SpicyPair::create("A", ""))
    false
  } catch {
    SpicyError(_) => true
  }
  let bad_association = try {
    ignore(@src.SpicyAssociation::create("image", "A", "B", 0.0, from_count=-1))
    false
  } catch {
    SpicyError(_) => true
  }
  assert_true(bad_cell)
  assert_true(bad_coordinate)
  assert_true(bad_pair)
  assert_true(bad_association)
}

///|
test "spicyR: configuration rejects invalid radii" {
  let empty = try {
    ignore(@src.SpicyConfig::create(radii=[]))
    false
  } catch {
    SpicyError(_) => true
  }
  let non_positive = try {
    ignore(@src.SpicyConfig::create(radii=[0.0, 1.0]))
    false
  } catch {
    SpicyError(_) => true
  }
  let unsorted = try {
    ignore(@src.SpicyConfig::create(radii=[2.0, 1.0]))
    false
  } catch {
    SpicyError(_) => true
  }
  assert_true(empty)
  assert_true(non_positive)
  assert_true(unsorted)
}

///|
test "spicyR: configuration rejects invalid fitting controls" {
  let bad_fdr = try {
    ignore(@src.SpicyConfig::create(fdr_threshold=1.1))
    false
  } catch {
    SpicyError(_) => true
  }
  let bad_weight = try {
    ignore(@src.SpicyConfig::create(weight_factor=-1.0))
    false
  } catch {
    SpicyError(_) => true
  }
  let bad_iterations = try {
    ignore(@src.SpicyConfig::create(max_iterations=0))
    false
  } catch {
    SpicyError(_) => true
  }
  let bad_tolerance = try {
    ignore(@src.SpicyConfig::create(tolerance=1.0))
    false
  } catch {
    SpicyError(_) => true
  }
  assert_true(bad_fdr)
  assert_true(bad_weight)
  assert_true(bad_iterations)
  assert_true(bad_tolerance)
}

///|
test "spicyR: cross-L matches a deterministic rectangular window" {
  let metadata = spicyr_test_metadata(3, 1)
  let cells = spicyr_test_geometry_cells(metadata, 5.0, 5.0, 5.5, 5.0)
  let pair = spicyr_test_pair("A", "B")
  let associations = @src.spicy_pairwise(
    cells,
    metadata,
    pair,
    config=spicyr_test_config(false, false),
  ) catch {
    _ => abort("spicyR pairwise geometry should succeed")
  }
  let expected_l = (100.0 / 3.141592653589793).sqrt()
  assert_eq(associations.length(), 3)
  spicyr_test_close(associations[0].l_values[0], expected_l, 1.0e-10)
  spicyr_test_close(associations[0].statistic, expected_l - 1.0, 1.0e-10)
  assert_true(associations[0].observed)
  assert_true(associations[0].available)
}

///|
test "spicyR: edge correction increases boundary-source cross-L" {
  let metadata = spicyr_test_metadata(3, 1)
  let cells = spicyr_test_geometry_cells(metadata, 0.0, 0.0, 0.5, 0.0)
  let pair = spicyr_test_pair("A", "B")
  let plain = @src.spicy_pairwise(
    cells,
    metadata,
    pair,
    config=spicyr_test_config(false, false),
  ) catch {
    _ => abort("uncorrected spicyR pairwise fit should succeed")
  }
  let corrected = @src.spicy_pairwise(
    cells,
    metadata,
    pair,
    config=spicyr_test_config(true, false),
  ) catch {
    _ => abort("edge-corrected spicyR pairwise fit should succeed")
  }
  assert_true(corrected[0].l_values[0] > plain[0].l_values[0] * 1.9)
}

///|
test "spicyR: ordered pairs retain directional edge correction" {
  let metadata = spicyr_test_metadata(3, 1)
  let cells = spicyr_test_geometry_cells(metadata, 0.0, 0.0, 0.5, 0.0)
  let forward = @src.spicy_pairwise(
    cells,
    metadata,
    spicyr_test_pair("A", "B"),
    config=spicyr_test_config(true, false),
  ) catch {
    _ => abort("forward spicyR pair should succeed")
  }
  let reverse = @src.spicy_pairwise(
    cells,
    metadata,
    spicyr_test_pair("B", "A"),
    config=spicyr_test_config(true, false),
  ) catch {
    _ => abort("reverse spicyR pair should succeed")
  }
  assert_true(forward[0].l_values[0] > reverse[0].l_values[0])
}

///|
test "spicyR: same-type cross-L excludes each cell from itself" {
  let metadata = spicyr_test_metadata(3, 1)
  let cells : Array[@src.SpicyCell] = []
  for image in metadata {
    cells.push(spicyr_test_cell(image.image_id, "A", 5.0, 5.0))
    cells.push(spicyr_test_cell(image.image_id, "A", 5.5, 5.0))
    cells.push(spicyr_test_cell(image.image_id, "Frame", 0.0, 0.0))
    cells.push(spicyr_test_cell(image.image_id, "Frame", 10.0, 10.0))
  }
  let associations = @src.spicy_pairwise(
    cells,
    metadata,
    spicyr_test_pair("A", "A"),
    config=spicyr_test_config(false, false),
  ) catch {
    _ => abort("same-type spicyR pair should succeed")
  }
  let expected_l = (50.0 / 3.141592653589793).sqrt()
  spicyr_test_close(associations[0].l_values[0], expected_l, 1.0e-10)
}

///|
test "spicyR: a single same-type cell retains the Poisson baseline" {
  let metadata = spicyr_test_metadata(3, 1)
  let cells : Array[@src.SpicyCell] = []
  for image in metadata {
    cells.push(spicyr_test_cell(image.image_id, "A", 5.0, 5.0))
    cells.push(spicyr_test_cell(image.image_id, "Frame", 0.0, 0.0))
    cells.push(spicyr_test_cell(image.image_id, "Frame", 10.0, 10.0))
  }
  let associations = @src.spicy_pairwise(
    cells,
    metadata,
    spicyr_test_pair("A", "A"),
    config=spicyr_test_config(false, false),
  ) catch {
    _ => abort("single-cell same-type pair should succeed")
  }
  assert_true(!associations[0].observed)
  assert_true(associations[0].available)
  assert_eq(associations[0].l_values[0], 0.0)
  assert_eq(associations[0].statistic, -1.0)
}

///|
test "spicyR: absent cell types follow include-zero-cells policy" {
  let metadata = spicyr_test_metadata(3, 1)
  let cells : Array[@src.SpicyCell] = []
  for index in 0..<metadata.length() {
    let image_id = metadata[index].image_id
    cells.push(spicyr_test_cell(image_id, "A", 5.0, 5.0))
    if index > 0 {
      cells.push(spicyr_test_cell(image_id, "B", 5.5, 5.0))
    }
    cells.push(spicyr_test_cell(image_id, "Frame", 0.0, 0.0))
    cells.push(spicyr_test_cell(image_id, "Frame", 10.0, 10.0))
  }
  let pair = spicyr_test_pair("A", "B")
  let omitted = @src.spicy_pairwise(
    cells,
    metadata,
    pair,
    config=spicyr_test_config(false, false),
  ) catch {
    _ => abort("spicyR missing-pair calculation should succeed")
  }
  let retained = @src.spicy_pairwise(
    cells,
    metadata,
    pair,
    config=spicyr_test_config(false, true),
  ) catch {
    _ => abort("spicyR zero-cell calculation should succeed")
  }
  assert_true(!omitted[0].available)
  assert_true(retained[0].available)
  assert_eq(retained[0].statistic, -1.0)
}

///|
test "spicyR: radii are capped and duplicate caps collapse" {
  let metadata = spicyr_test_metadata(3, 1)
  let cells : Array[@src.SpicyCell] = []
  for image in metadata {
    cells.push(spicyr_test_cell(image.image_id, "A", 0.5, 2.0))
    cells.push(spicyr_test_cell(image.image_id, "B", 1.0, 2.0))
    cells.push(spicyr_test_cell(image.image_id, "Frame", 0.0, 0.0))
    cells.push(spicyr_test_cell(image.image_id, "Frame", 2.0, 4.0))
  }
  let config = @src.SpicyConfig::create(
    radii=[0.5, 5.0, 10.0],
    edge_correction=false,
  ) catch {
    _ => abort("radius-cap configuration should be valid")
  }
  let associations = @src.spicy_pairwise(
    cells,
    metadata,
    spicyr_test_pair("A", "B"),
    config~,
  ) catch {
    _ => abort("radius-cap calculation should succeed")
  }
  assert_eq(associations[0].radii.length(), 2)
  spicyr_test_close(associations[0].radii[1], 2.0 / 2.01, 1.0e-12)
}

///|
test "spicyR: pairwise validation rejects unknown cell types" {
  let metadata = spicyr_test_metadata(3, 1)
  let cells = spicyr_test_geometry_cells(metadata, 5.0, 5.0, 5.5, 5.0)
  let failed = try {
    ignore(
      @src.spicy_pairwise(cells, metadata, spicyr_test_pair("A", "missing")),
    )
    false
  } catch {
    SpicyError(_) => true
  }
  assert_true(failed)
}

///|
test "spicyR: pairwise validation rejects duplicate image metadata" {
  let metadata = spicyr_test_metadata(3, 1)
  metadata[2] = metadata[1]
  let cells = spicyr_test_geometry_cells(metadata, 5.0, 5.0, 5.5, 5.0)
  let failed = try {
    ignore(@src.spicy_pairwise(cells, metadata, spicyr_test_pair("A", "B")))
    false
  } catch {
    SpicyError(_) => true
  }
  assert_true(failed)
}

///|
test "spicyR: weighted linear model estimates a condition contrast" {
  let (metadata, associations) = spicyr_test_linear_input()
  let result = @src.spicy_fit_associations(associations, metadata) catch {
    _ => abort("spicyR weighted linear fit should succeed")
  }
  let pair = match result.pair("A", "B") {
    Some(value) => value
    None => abort("spicyR pair should be present")
  }
  let condition = match pair.condition("treated") {
    Some(value) => value
    None => abort("spicyR treated contrast should be present")
  }
  assert_true(pair.model is @src.SpicyModelKind::WeightedLinear)
  spicyr_test_close(condition.estimate, 2.0, 0.25)
  assert_true(condition.standard_error > 0.0)
  assert_true(condition.p_value >= 0.0 && condition.p_value <= 1.0)
}

///|
test "spicyR: count-derived precision weights normalize to mean one" {
  let (metadata, associations) = spicyr_test_linear_input()
  let result = @src.spicy_fit_associations(associations, metadata) catch {
    _ => abort("spicyR weighted fit should succeed")
  }
  let pair = result.pairs[0]
  let mut total = 0.0
  for association in pair.associations {
    total = total + association.weight
  }
  spicyr_test_close(
    total / pair.associations.length().to_double(),
    1.0,
    1.0e-12,
  )
  assert_true(pair.associations[0].weight < pair.associations[7].weight)
}

///|
test "spicyR: disabling precision weights assigns unit weights" {
  let (metadata, associations) = spicyr_test_linear_input()
  let config = @src.SpicyConfig::create(use_weights=false) catch {
    _ => abort("unweighted spicyR configuration should be valid")
  }
  let result = @src.spicy_fit_associations(associations, metadata, config~) catch {
    _ => abort("unweighted spicyR fit should succeed")
  }
  for association in result.pairs[0].associations {
    assert_eq(association.weight, 1.0)
  }
}

///|
test "spicyR: numeric covariates enter the fixed-effect design" {
  let (metadata, associations) = spicyr_test_linear_input()
  let result = @src.spicy_fit_associations(associations, metadata, covariate_names=[
    "age",
  ]) catch {
    _ => abort("spicyR covariate fit should succeed")
  }
  let condition = result.pairs[0].conditions[0]
  assert_true(condition.estimate == condition.estimate)
  assert_true(condition.standard_error > 0.0)
}

///|
test "spicyR: repeated subjects select the random-intercept model" {
  let (metadata, associations) = spicyr_test_random_input()
  let result = @src.spicy_fit_associations(associations, metadata) catch {
    _ => abort("spicyR random-intercept fit should succeed")
  }
  let pair = result.pairs[0]
  let condition = pair.conditions[0]
  assert_true(pair.model is @src.SpicyModelKind::WeightedRandomIntercept)
  spicyr_test_close(condition.estimate, 1.5, 0.15)
  assert_true(pair.random_intercept_variance > 0.1)
  assert_true(pair.residual_variance >= 0.0)
}

///|
test "spicyR: multiple conditions produce reference contrasts" {
  let metadata : Array[@src.SpicyImageMetadata] = []
  let associations : Array[@src.SpicyAssociation] = []
  let levels = ["control", "drug_a", "drug_b"]
  for level_index in 0..<levels.length() {
    for replicate in 0..<3 {
      let image_id = levels[level_index] + "_" + (replicate + 1).to_string()
      metadata.push(
        spicyr_test_image_metadata(image_id, levels[level_index], "", Map([])),
      )
      associations.push(
        spicyr_test_association(
          image_id,
          "A",
          "B",
          level_index.to_double() + (replicate - 1).to_double() * 0.1,
          6,
        ),
      )
    }
  }
  let config = @src.SpicyConfig::create(reference_condition="control") catch {
    _ => abort("multi-condition spicyR configuration should be valid")
  }
  let result = @src.spicy_fit_associations(associations, metadata, config~) catch {
    _ => abort("multi-condition spicyR fit should succeed")
  }
  assert_eq(result.pairs[0].conditions.length(), 2)
  let drug_a = match result.pairs[0].condition("drug_a") {
    Some(value) => value
    None => abort("drug_a contrast should be present")
  }
  let drug_b = match result.pairs[0].condition("drug_b") {
    Some(value) => value
    None => abort("drug_b contrast should be present")
  }
  spicyr_test_close(drug_a.estimate, 1.0, 0.05)
  spicyr_test_close(drug_b.estimate, 2.0, 0.05)
}

///|
test "spicyR: BH correction is bounded and rank monotone" {
  let result = spicyr_test_two_pair_result()
  let first = result.pairs[0].conditions[0]
  let second = result.pairs[1].conditions[0]
  assert_true(first.adjusted_p_value >= first.p_value)
  assert_true(second.adjusted_p_value >= second.p_value)
  assert_true(first.adjusted_p_value <= 1.0)
  assert_true(second.adjusted_p_value <= 1.0)
  if first.p_value <= second.p_value {
    assert_true(first.adjusted_p_value <= second.adjusted_p_value)
  } else {
    assert_true(second.adjusted_p_value <= first.adjusted_p_value)
  }
}

///|
test "spicyR: pair and association lookup preserve identities" {
  let result = spicyr_test_two_pair_result()
  let pair = match result.pair("A", "B") {
    Some(value) => value
    None => abort("spicyR A__B pair should exist")
  }
  let association = match pair.association("image_3") {
    Some(value) => value
    None => abort("spicyR image association should exist")
  }
  assert_eq(association.image_id, "image_3")
  assert_true(result.pair("missing", "B") is None)
  assert_true(pair.association("missing") is None)
}

///|
test "spicyR: top ranks pairs and handles zero limits" {
  let result = spicyr_test_two_pair_result()
  let top = result.top("treated", limit=2) catch {
    _ => abort("spicyR top query should succeed")
  }
  let first_p = top[0].conditions[0].p_value
  let second_p = top[1].conditions[0].p_value
  assert_true(first_p <= second_p)
  let empty = result.top("treated", limit=0) catch {
    _ => abort("spicyR empty top query should succeed")
  }
  assert_eq(empty.length(), 0)
}

///|
test "spicyR: significant filtering and summary use configured FDR" {
  let result = spicyr_test_two_pair_result()
  let all = result.significant("treated", threshold=1.0) catch {
    _ => abort("spicyR significant query should succeed")
  }
  assert_eq(all.length(), result.n_pairs())
  assert_eq(result.n_images(), 8)
  assert_true(result.summary().contains("spicyR(images=8, pairs=2"))
}

///|
test "spicyR: result queries reject absent contrasts and invalid limits" {
  let result = spicyr_test_two_pair_result()
  let bad_level = try {
    ignore(result.top("control"))
    false
  } catch {
    SpicyError(_) => true
  }
  let bad_limit = try {
    ignore(result.top("treated", limit=-1))
    false
  } catch {
    SpicyError(_) => true
  }
  let bad_threshold = try {
    ignore(result.significant("treated", threshold=1.1))
    false
  } catch {
    SpicyError(_) => true
  }
  assert_true(bad_level)
  assert_true(bad_limit)
  assert_true(bad_threshold)
}

///|
test "spicyR: fitting rejects a single condition" {
  let (metadata, associations) = spicyr_test_linear_input()
  for index in 0..<metadata.length() {
    metadata[index] = spicyr_test_image_metadata(
      metadata[index].image_id,
      "control",
      "",
      metadata[index].covariates,
    )
  }
  let failed = try {
    ignore(@src.spicy_fit_associations(associations, metadata))
    false
  } catch {
    SpicyError(_) => true
  }
  assert_true(failed)
}

///|
test "spicyR: fitting rejects duplicate pair-image associations" {
  let (metadata, associations) = spicyr_test_linear_input()
  associations.push(associations[0])
  let failed = try {
    ignore(@src.spicy_fit_associations(associations, metadata))
    false
  } catch {
    SpicyError(_) => true
  }
  assert_true(failed)
}

///|
test "spicyR: fitting rejects missing covariates" {
  let (metadata, associations) = spicyr_test_linear_input()
  let failed = try {
    ignore(
      @src.spicy_fit_associations(associations, metadata, covariate_names=[
        "missing",
      ]),
    )
    false
  } catch {
    SpicyError(_) => true
  }
  assert_true(failed)
}

///|
test "spicyR: direct cell analysis returns requested pairs" {
  let (cells, metadata, pairs) = @src.spicyr_example_data()
  let config = @src.SpicyConfig::create(
    radii=[0.5, 1.0, 2.0],
    edge_correction=true,
    reference_condition="control",
    tolerance=0.01,
  ) catch {
    _ => abort("spicyR example configuration should be valid")
  }
  let result = @src.spicyr(
    cells,
    metadata,
    pairs~,
    covariate_names=["age"],
    config~,
  ) catch {
    _ => abort("spicyR example analysis should succeed")
  }
  assert_eq(result.n_images(), 12)
  assert_eq(result.n_pairs(), 1)
  assert_true(
    result.pairs[0].model is @src.SpicyModelKind::WeightedRandomIntercept,
  )
  assert_true(result.pairs[0].conditions[0].estimate > 0.0)
}

///|
test "spicyR: automatic pair generation includes ordered pairs" {
  let (cells, metadata, _) = @src.spicyr_example_data()
  let config = @src.SpicyConfig::create(
    radii=[1.0],
    edge_correction=false,
    include_zero_cells=true,
    use_weights=false,
    tolerance=0.01,
  ) catch {
    _ => abort("automatic-pair configuration should be valid")
  }
  let result = @src.spicyr(cells, metadata, config~) catch {
    _ => abort("automatic spicyR analysis should succeed")
  }
  assert_eq(result.n_pairs(), 4)
  assert_true(result.pair("T_cell", "Tumour") is Some(_))
  assert_true(result.pair("Tumour", "T_cell") is Some(_))
}

///|
test "spicyR: direct analysis rejects duplicate requested pairs" {
  let (cells, metadata, pairs) = @src.spicyr_example_data()
  let duplicate = [pairs[0], pairs[0]]
  let failed = try {
    ignore(@src.spicyr(cells, metadata, pairs=duplicate))
    false
  } catch {
    SpicyError(_) => true
  }
  assert_true(failed)
}

///|
test "spicyR: SpatialExperiment integration preserves and enriches input" {
  let experiment = spicyr_test_spatial_experiment()
  experiment.metadata["source"] = "fixture"
  let pair = spicyr_test_pair("A", "B")
  let output = @src.spicyr_spatial_experiment(
    experiment,
    pairs=[pair],
    config=spicyr_test_config(false, false),
  ) catch {
    _ => abort("spicyR SpatialExperiment integration should succeed")
  }
  assert_eq(output.result.n_images(), 6)
  assert_eq(output.result.n_pairs(), 1)
  assert_eq(output.experiment.metadata["source"], "fixture")
  assert_eq(output.experiment.metadata["spicyr_images"], "6")
  assert_eq(output.experiment.metadata["spicyr_pairs"], "1")
  assert_true(!experiment.metadata.contains("spicyr_images"))
}

///|
test "spicyR: SpatialExperiment supports subject and numeric covariate keys" {
  let experiment = spicyr_test_spatial_experiment()
  let pair = spicyr_test_pair("A", "B")
  let output = @src.spicyr_spatial_experiment(
    experiment,
    subject_key="subject",
    covariate_keys=["age"],
    pairs=[pair],
    config=spicyr_test_config(false, false),
  ) catch {
    _ => abort("spicyR keyed SpatialExperiment fit should succeed")
  }
  assert_true(
    output.result.pairs[0].model is @src.SpicyModelKind::WeightedRandomIntercept,
  )
}

///|
test "spicyR: SpatialExperiment rejects coordinate mismatch" {
  let experiment = spicyr_test_spatial_experiment()
  ignore(experiment.spatial_coords.pop())
  let failed = try {
    ignore(@src.spicyr_spatial_experiment(experiment))
    false
  } catch {
    SpicyError(_) => true
  }
  assert_true(failed)
}

///|
test "spicyR: SpatialExperiment rejects missing required columns" {
  let experiment = spicyr_test_spatial_experiment()
  experiment.col_data[0].remove("cellType")
  let failed = try {
    ignore(@src.spicyr_spatial_experiment(experiment))
    false
  } catch {
    SpicyError(_) => true
  }
  assert_true(failed)
}

///|
test "spicyR: SpatialExperiment rejects image-varying metadata" {
  let experiment = spicyr_test_spatial_experiment()
  experiment.col_data[1]["condition"] = "other"
  let failed = try {
    ignore(@src.spicyr_spatial_experiment(experiment))
    false
  } catch {
    SpicyError(_) => true
  }
  assert_true(failed)
}

///|
test "spicyR: SpatialExperiment rejects nonnumeric covariates" {
  let experiment = spicyr_test_spatial_experiment()
  experiment.col_data[0]["age"] = "unknown"
  let failed = try {
    ignore(@src.spicyr_spatial_experiment(experiment, covariate_keys=["age"]))
    false
  } catch {
    SpicyError(_) => true
  }
  assert_true(failed)
}
