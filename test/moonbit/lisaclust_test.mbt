// Black-box tests for the Bioconductor lisaClust-inspired workflow.

///|
fn lisa_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "lisaClust value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn lisa_test_cell(
  cell_id : String,
  image_id : String,
  cell_type : String,
  x : Double,
  y : Double,
) -> @src.LisaCell {
  @src.LisaCell::create(cell_id, image_id, cell_type, x, y) catch {
    _ => abort("lisaClust test cell should be valid")
  }
}

///|
fn lisa_test_config(
  curve_kind? : @src.LisaCurveKind = @src.lisa_k_curve(),
  edge_correction? : Bool = false,
  window_kind? : @src.LisaWindowKind = @src.lisa_rectangle_window(),
  n_clusters? : Int = 2,
  region_prefix? : String = "region",
) -> @src.LisaConfig {
  @src.LisaConfig::create(
    radii=[1.5, 4.0],
    curve_kind~,
    window_kind~,
    bandwidth=1000.0,
    min_density=0.05,
    edge_correction~,
    window_padding=0.0,
    edge_samples=128,
    n_clusters~,
    n_starts=4,
    max_iterations=100,
    tolerance=1.0e-9,
    seed=7,
    region_prefix~,
  ) catch {
    _ => abort("lisaClust test configuration should be valid")
  }
}

///|
fn lisa_test_cells() -> Array[@src.LisaCell] {
  [
    lisa_test_cell("a1", "image_1", "A", 0.0, 0.0),
    lisa_test_cell("a2", "image_1", "A", 1.0, 1.0),
    lisa_test_cell("a3", "image_1", "A", 2.0, 0.5),
    lisa_test_cell("b1", "image_1", "B", 8.0, 8.0),
    lisa_test_cell("b2", "image_1", "B", 9.0, 9.0),
    lisa_test_cell("b3", "image_1", "B", 10.0, 8.5),
    lisa_test_cell("c1", "image_1", "C", 0.0, 10.0),
    lisa_test_cell("c2", "image_1", "C", 10.0, 0.0),
  ]
}

///|
fn lisa_test_exact_cells() -> Array[@src.LisaCell] {
  [
    lisa_test_cell("source", "image_1", "A", 2.0, 2.0),
    lisa_test_cell("near", "image_1", "B", 3.0, 2.0),
    lisa_test_cell("far", "image_1", "B", 8.0, 8.0),
    lisa_test_cell("corner_1", "image_1", "C", 0.0, 0.0),
    lisa_test_cell("corner_2", "image_1", "C", 10.0, 10.0),
  ]
}

///|
fn lisa_test_result() -> @src.LisaResult {
  @src.lisaclust(lisa_test_cells(), config=lisa_test_config()) catch {
    _ => abort("lisaClust test workflow should succeed")
  }
}

///|
fn lisa_test_experiment() -> @src.SpatialExperiment {
  let experiment = @src.SpatialExperiment::new()
  for cell in lisa_test_cells() {
    ignore(
      @src.se_add_col(
        experiment,
        Map(
          [
            ("cellID", cell.cell_id),
            ("imageID", cell.image_id),
            ("cellType", cell.cell_type),
          ],
        ),
      ),
    )
    ignore(
      @src.se_add_spatial_coord(
        experiment,
        @src.SpatialCoord::new_2d(cell.x, cell.y),
      ),
    )
  }
  experiment
}

///|
test "lisaClust: curve and window helpers expose both modes" {
  assert_true(
    @src.lisa_k_curve() is @src.LisaCurveKind::LisaStandardizedK,
  )
  assert_true(@src.lisa_l_curve() is @src.LisaCurveKind::LisaCenteredL)
  assert_true(
    @src.lisa_rectangle_window() is @src.LisaWindowKind::LisaRectangle,
  )
  assert_true(
    @src.lisa_convex_window() is @src.LisaWindowKind::LisaConvexHull,
  )
}

///|
test "lisaClust: default configuration follows upstream defaults" {
  let config = @src.LisaConfig::default()
  assert_eq(config.radii, [20.0, 50.0, 100.0])
  assert_true(config.curve_kind is @src.LisaCurveKind::LisaStandardizedK)
  assert_true(config.window_kind is @src.LisaWindowKind::LisaConvexHull)
  assert_eq(config.bandwidth, 100000.0)
  assert_eq(config.min_density, 0.05)
  assert_eq(config.edge_correction, true)
  assert_eq(config.n_clusters, 2)
  assert_eq(config.region_prefix, "region")
}

///|
test "lisaClust: custom configuration preserves controls" {
  let config = @src.LisaConfig::create(
    radii=[1.0, 3.0],
    curve_kind=@src.lisa_l_curve(),
    window_kind=@src.lisa_rectangle_window(),
    bandwidth=2.5,
    min_density=0.2,
    edge_correction=false,
    window_padding=1.0,
    edge_samples=64,
    n_clusters=3,
    n_starts=5,
    max_iterations=40,
    tolerance=1.0e-6,
    seed=9,
    region_prefix="domain",
  ) catch {
    _ => abort("custom lisaClust configuration should be valid")
  }
  assert_eq(config.radii, [1.0, 3.0])
  assert_true(config.curve_kind is @src.LisaCurveKind::LisaCenteredL)
  assert_eq(config.bandwidth, 2.5)
  assert_eq(config.edge_correction, false)
  assert_eq(config.n_clusters, 3)
  assert_eq(config.region_prefix, "domain")
}

///|
test "lisaClust: configuration defensively copies radii" {
  let radii = [1.0, 2.0]
  let config = @src.LisaConfig::create(radii~) catch {
    _ => abort("lisaClust radii should be valid")
  }
  radii[0] = 99.0
  assert_eq(config.radii, [1.0, 2.0])
}

///|
test "lisaClust: configuration rejects empty radii" {
  let failed = try {
    ignore(@src.LisaConfig::create(radii=[]))
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}

///|
test "lisaClust: configuration rejects non-increasing radii" {
  let failed = try {
    ignore(@src.LisaConfig::create(radii=[2.0, 2.0]))
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}

///|
test "lisaClust: configuration rejects non-finite radii" {
  let failed = try {
    ignore(@src.LisaConfig::create(radii=[1.0, @double.not_a_number]))
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}

///|
test "lisaClust: configuration rejects invalid density controls" {
  let bad_bandwidth = try {
    ignore(@src.LisaConfig::create(bandwidth=0.0))
    false
  } catch {
    LisaError(_) => true
  }
  let bad_minimum = try {
    ignore(@src.LisaConfig::create(min_density=0.0))
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(bad_bandwidth)
  assert_true(bad_minimum)
}

///|
test "lisaClust: configuration rejects invalid edge controls" {
  let bad_padding = try {
    ignore(@src.LisaConfig::create(window_padding=-1.0))
    false
  } catch {
    LisaError(_) => true
  }
  let bad_samples = try {
    ignore(@src.LisaConfig::create(edge_samples=16))
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(bad_padding)
  assert_true(bad_samples)
}

///|
test "lisaClust: configuration rejects invalid clustering controls" {
  let bad_clusters = try {
    ignore(@src.LisaConfig::create(n_clusters=0))
    false
  } catch {
    LisaError(_) => true
  }
  let bad_starts = try {
    ignore(@src.LisaConfig::create(n_starts=0))
    false
  } catch {
    LisaError(_) => true
  }
  let bad_iterations = try {
    ignore(@src.LisaConfig::create(max_iterations=0))
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(bad_clusters)
  assert_true(bad_starts)
  assert_true(bad_iterations)
}

///|
test "lisaClust: configuration rejects invalid tolerance seed and prefix" {
  let bad_tolerance = try {
    ignore(@src.LisaConfig::create(tolerance=1.0))
    false
  } catch {
    LisaError(_) => true
  }
  let bad_seed = try {
    ignore(@src.LisaConfig::create(seed=-1))
    false
  } catch {
    LisaError(_) => true
  }
  let bad_prefix = try {
    ignore(@src.LisaConfig::create(region_prefix=""))
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(bad_tolerance)
  assert_true(bad_seed)
  assert_true(bad_prefix)
}

///|
test "lisaClust: cell constructor validates identifiers" {
  let bad_cell = try {
    ignore(@src.LisaCell::create("", "image", "A", 0.0, 0.0))
    false
  } catch {
    LisaError(_) => true
  }
  let bad_image = try {
    ignore(@src.LisaCell::create("cell", "", "A", 0.0, 0.0))
    false
  } catch {
    LisaError(_) => true
  }
  let bad_type = try {
    ignore(@src.LisaCell::create("cell", "image", "", 0.0, 0.0))
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(bad_cell)
  assert_true(bad_image)
  assert_true(bad_type)
}

///|
test "lisaClust: cell constructor rejects non-finite coordinates" {
  let failed = try {
    ignore(
      @src.LisaCell::create(
        "cell",
        "image",
        "A",
        @double.not_a_number,
        0.0,
      ),
    )
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}

///|
test "lisaClust: curve matrix exposes stable dimensions and names" {
  let curves = @src.lisa_curves(
    lisa_test_cells(),
    config=lisa_test_config(),
  ) catch {
    _ => abort("lisaClust curves should compute")
  }
  assert_eq(curves.n_cells(), 8)
  assert_eq(curves.target_cell_types, ["A", "B", "C"])
  assert_eq(curves.radii, [1.5, 4.0])
  assert_eq(curves.n_features(), 6)
  assert_eq(curves.feature_names[0], "1.5_A")
  assert_eq(curves.feature_names[5], "4_C")
}

///|
test "lisaClust: feature lookup follows target-major radius order" {
  let curves = @src.lisa_curves(
    lisa_test_cells(),
    config=lisa_test_config(),
  ) catch {
    _ => abort("lisaClust curves should compute")
  }
  assert_eq(curves.feature_index("A", 1.5), 0)
  assert_eq(curves.feature_index("B", 4.0), 3)
  assert_eq(curves.feature_index("missing", 1.5), -1)
  assert_eq(curves.feature_index("A", 99.0), -1)
}

///|
test "lisaClust: exact standardized K uses observed-minus-expected scaling" {
  let config = @src.LisaConfig::create(
    radii=[2.0],
    window_kind=@src.lisa_rectangle_window(),
    bandwidth=1000000.0,
    edge_correction=false,
    window_padding=0.0,
    edge_samples=64,
  ) catch {
    _ => abort("exact K configuration should be valid")
  }
  let curves = @src.lisa_curves(lisa_test_exact_cells(), config~) catch {
    _ => abort("exact K curves should compute")
  }
  let feature = curves.feature_index("B", 2.0)
  let expected = @math.PI * 4.0 * 2.0 / 100.0
  let target = (1.0 - expected) / expected.sqrt()
  lisa_test_close(curves.values[0][feature], target, 1.0e-5)
}

///|
test "lisaClust: centered L uses square-root observed and expected counts" {
  let config = @src.LisaConfig::create(
    radii=[2.0],
    curve_kind=@src.lisa_l_curve(),
    window_kind=@src.lisa_rectangle_window(),
    bandwidth=1000000.0,
    edge_correction=false,
    window_padding=0.0,
    edge_samples=64,
  ) catch {
    _ => abort("exact L configuration should be valid")
  }
  let curves = @src.lisa_curves(lisa_test_exact_cells(), config~) catch {
    _ => abort("exact L curves should compute")
  }
  let feature = curves.feature_index("B", 2.0)
  let expected = @math.PI * 4.0 * 2.0 / 100.0
  lisa_test_close(curves.values[0][feature], 1.0 - expected.sqrt(), 1.0e-5)
}

///|
test "lisaClust: local curves exclude self matches" {
  let config = @src.LisaConfig::create(
    radii=[0.5],
    window_kind=@src.lisa_rectangle_window(),
    bandwidth=1000000.0,
    edge_correction=false,
    window_padding=0.0,
  ) catch {
    _ => abort("self exclusion configuration should be valid")
  }
  let curves = @src.lisa_curves(lisa_test_exact_cells(), config~) catch {
    _ => abort("self exclusion curves should compute")
  }
  let feature = curves.feature_index("A", 0.5)
  let expected = @math.PI * 0.25 / 100.0
  lisa_test_close(curves.values[0][feature], -expected.sqrt(), 1.0e-6)
}

///|
test "lisaClust: rectangle window area follows coordinate extent" {
  let curves = @src.lisa_curves(
    lisa_test_exact_cells(),
    config=lisa_test_config(),
  ) catch {
    _ => abort("rectangle lisaClust curves should compute")
  }
  assert_eq(curves.windows.length(), 1)
  lisa_test_close(curves.windows[0].area, 100.0, 1.0e-12)
}

///|
test "lisaClust: convex window uses hull area" {
  let cells = [
    lisa_test_cell("p1", "image", "A", 0.0, 0.0),
    lisa_test_cell("p2", "image", "A", 10.0, 0.0),
    lisa_test_cell("p3", "image", "B", 0.0, 10.0),
    lisa_test_cell("inside", "image", "B", 2.0, 2.0),
  ]
  let config = @src.LisaConfig::create(
    radii=[1.0],
    window_kind=@src.lisa_convex_window(),
    window_padding=0.0,
    edge_correction=false,
  ) catch {
    _ => abort("convex lisaClust configuration should be valid")
  }
  let curves = @src.lisa_curves(cells, config~) catch {
    _ => abort("convex lisaClust curves should compute")
  }
  lisa_test_close(curves.windows[0].area, 50.0, 1.0e-12)
  assert_eq(curves.windows[0].vertices.length(), 3)
}

///|
test "lisaClust: radii are capped at half the shortest window span" {
  let config = @src.LisaConfig::create(
    radii=[1.0, 100.0],
    window_kind=@src.lisa_rectangle_window(),
    window_padding=0.0,
    edge_correction=false,
  ) catch {
    _ => abort("radius cap configuration should be valid")
  }
  let curves = @src.lisa_curves(lisa_test_exact_cells(), config~) catch {
    _ => abort("radius capped curves should compute")
  }
  let radii = match curves.effective_radii_for_image("image_1") {
    Some(values) => values
    None => abort("effective radii should exist")
  }
  assert_eq(radii[0], 1.0)
  lisa_test_close(radii[1], 10.0 / 2.01, 1.0e-12)
}

///|
test "lisaClust: edge correction changes corner-cell expectation" {
  let cells = [
    lisa_test_cell("source", "image", "A", 0.0, 0.0),
    lisa_test_cell("near", "image", "B", 1.0, 1.0),
    lisa_test_cell("far", "image", "B", 8.0, 8.0),
    lisa_test_cell("corner", "image", "C", 10.0, 10.0),
  ]
  let without = @src.lisa_curves(
    cells,
    config=@src.LisaConfig::create(
      radii=[2.0],
      window_kind=@src.lisa_rectangle_window(),
      bandwidth=1000000.0,
      edge_correction=false,
      window_padding=0.0,
      edge_samples=720,
    ) catch {
      _ => abort("uncorrected edge configuration should be valid")
    },
  ) catch {
    _ => abort("uncorrected curves should compute")
  }
  let with_edge = @src.lisa_curves(
    cells,
    config=@src.LisaConfig::create(
      radii=[2.0],
      window_kind=@src.lisa_rectangle_window(),
      bandwidth=1000000.0,
      edge_correction=true,
      window_padding=0.0,
      edge_samples=720,
    ) catch {
      _ => abort("corrected edge configuration should be valid")
    },
  ) catch {
    _ => abort("corrected curves should compute")
  }
  let feature = without.feature_index("B", 2.0)
  assert_true(with_edge.values[0][feature] > without.values[0][feature])
}

///|
test "lisaClust: KDE weighting changes curves under inhomogeneous density" {
  let narrow = @src.LisaConfig::create(
    radii=[2.0],
    window_kind=@src.lisa_rectangle_window(),
    bandwidth=0.25,
    min_density=0.05,
    edge_correction=false,
    window_padding=0.0,
  ) catch {
    _ => abort("narrow KDE configuration should be valid")
  }
  let broad = @src.LisaConfig::create(
    radii=[2.0],
    window_kind=@src.lisa_rectangle_window(),
    bandwidth=1000000.0,
    min_density=0.05,
    edge_correction=false,
    window_padding=0.0,
  ) catch {
    _ => abort("broad KDE configuration should be valid")
  }
  let narrow_curves = @src.lisa_curves(
    lisa_test_exact_cells(),
    config=narrow,
  ) catch {
    _ => abort("narrow KDE curves should compute")
  }
  let broad_curves = @src.lisa_curves(
    lisa_test_exact_cells(),
    config=broad,
  ) catch {
    _ => abort("broad KDE curves should compute")
  }
  let feature = narrow_curves.feature_index("B", 2.0)
  assert_true(
    (narrow_curves.values[0][feature] - broad_curves.values[0][feature]).abs() >
    1.0e-6,
  )
}

///|
test "lisaClust: generated curves are finite" {
  let curves = @src.lisa_curves(
    lisa_test_cells(),
    config=lisa_test_config(edge_correction=true),
  ) catch {
    _ => abort("edge corrected curves should compute")
  }
  for row in curves.values {
    for value in row {
      assert_true(value == value)
      assert_true(value.abs() <= 1.0e300)
    }
  }
}

///|
test "lisaClust: curve lookup returns a defensive copy" {
  let curves = @src.lisa_curves(
    lisa_test_cells(),
    config=lisa_test_config(),
  ) catch {
    _ => abort("lisaClust curves should compute")
  }
  let curve = match curves.curve_for_cell("a1") {
    Some(value) => value
    None => abort("cell curve should exist")
  }
  curve[0] = 999.0
  assert_true(curves.values[0][0] != 999.0)
  assert_true(curves.curve_for_cell("missing") is None)
}

///|
test "lisaClust: duplicate cell IDs are rejected" {
  let cells = lisa_test_cells()
  cells.push(lisa_test_cell("a1", "image_1", "A", 4.0, 4.0))
  let failed = try {
    ignore(@src.lisa_curves(cells, config=lisa_test_config()))
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}

///|
test "lisaClust: fewer than two cells are rejected" {
  let failed = try {
    ignore(
      @src.lisa_curves(
        [lisa_test_cell("only", "image", "A", 0.0, 0.0)],
        config=lisa_test_config(),
      ),
    )
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}

///|
test "lisaClust: one-cell images are rejected" {
  let cells = lisa_test_cells()
  cells.push(lisa_test_cell("single", "image_2", "A", 2.0, 2.0))
  let failed = try {
    ignore(@src.lisa_curves(cells, config=lisa_test_config()))
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}

///|
test "lisaClust: degenerate rectangular windows are rejected" {
  let cells = [
    lisa_test_cell("one", "image", "A", 0.0, 0.0),
    lisa_test_cell("two", "image", "B", 1.0, 0.0),
  ]
  let failed = try {
    ignore(@src.lisa_curves(cells, config=lisa_test_config()))
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}

///|
test "lisaClust: collinear convex windows are rejected" {
  let cells = [
    lisa_test_cell("one", "image", "A", 0.0, 0.0),
    lisa_test_cell("two", "image", "B", 1.0, 1.0),
    lisa_test_cell("three", "image", "B", 2.0, 2.0),
  ]
  let config = @src.LisaConfig::create(
    radii=[1.0],
    window_kind=@src.lisa_convex_window(),
  ) catch {
    _ => abort("convex configuration should be valid")
  }
  let failed = try {
    ignore(@src.lisa_curves(cells, config~))
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}

///|
test "lisaClust: complete workflow returns deterministic regions" {
  let first = lisa_test_result()
  let second = lisa_test_result()
  assert_eq(first.labels, second.labels)
  assert_eq(first.regions, second.regions)
  assert_eq(first.cluster_sizes, second.cluster_sizes)
  assert_eq(first.centroids, second.centroids)
}

///|
test "lisaClust: complete workflow reports fit diagnostics" {
  let result = lisa_test_result()
  assert_eq(result.n_cells(), 8)
  assert_eq(result.n_regions(), 2)
  assert_eq(result.labels.length(), 8)
  assert_eq(result.centroids.length(), 2)
  assert_eq(result.cluster_sizes[0] + result.cluster_sizes[1], 8)
  assert_true(result.inertia >= 0.0)
  assert_true(result.silhouette >= -1.0)
  assert_true(result.silhouette <= 1.0)
  assert_true(result.iterations > 0)
  assert_true(result.converged)
}

///|
test "lisaClust: custom region prefix is propagated" {
  let result = @src.lisaclust(
    lisa_test_cells(),
    config=lisa_test_config(region_prefix="domain"),
  ) catch {
    _ => abort("custom prefix workflow should succeed")
  }
  for region in result.regions {
    assert_true(region.has_prefix("domain_"))
  }
}

///|
test "lisaClust: cell-to-region lookup handles known and unknown IDs" {
  let result = lisa_test_result()
  assert_true(result.region_for_cell("a1") is Some(_))
  assert_true(result.region_for_cell("missing") is None)
}

///|
test "lisaClust: region enrichment covers every type-region combination" {
  let result = lisa_test_result()
  assert_eq(result.enrichment.length(), 6)
  assert_eq(result.region_summaries.length(), 2)
  for summary in result.region_summaries {
    assert_true(summary.size > 0)
    assert_true(summary.dominant_cell_type != "")
    assert_true(summary.maximum_enrichment >= 0.0)
  }
}

///|
test "lisaClust: region enrichment uses observed over independence expectation" {
  let result = lisa_test_result()
  let entry = result.enrichment[0]
  if entry.expected > 0.0 {
    lisa_test_close(
      entry.relative_frequency,
      entry.observed.to_double() / entry.expected,
      1.0e-12,
    )
  }
  assert_true(
    result.region_enrichment(entry.cell_type, entry.region) is Some(_),
  )
  assert_true(result.region_enrichment("missing", entry.region) is None)
}

///|
test "lisaClust: top enrichments are sorted and filtered" {
  let result = lisa_test_result()
  let top = result.top_enrichments(
    limit=3,
    minimum_relative_frequency=0.0,
  ) catch {
    _ => abort("top enrichments should succeed")
  }
  assert_eq(top.length(), 3)
  assert_true(top[0].relative_frequency >= top[1].relative_frequency)
  assert_true(top[1].relative_frequency >= top[2].relative_frequency)
}

///|
test "lisaClust: top enrichment validates query controls" {
  let result = lisa_test_result()
  let bad_limit = try {
    ignore(result.top_enrichments(limit=-1))
    false
  } catch {
    LisaError(_) => true
  }
  let bad_threshold = try {
    ignore(result.top_enrichments(minimum_relative_frequency=-1.0))
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(bad_limit)
  assert_true(bad_threshold)
}

///|
test "lisaClust: summary reports cells features and regions" {
  let summary = lisa_test_result().summary()
  assert_true(summary.contains("lisaClust"))
  assert_true(summary.contains("cells=8"))
  assert_true(summary.contains("features=6"))
  assert_true(summary.contains("regions=2"))
}

///|
test "lisaClust: cluster count cannot exceed cells" {
  let curves = @src.lisa_curves(
    lisa_test_cells(),
    config=lisa_test_config(),
  ) catch {
    _ => abort("lisaClust curves should compute")
  }
  let config = lisa_test_config(n_clusters=9)
  let failed = try {
    ignore(@src.lisaclust_from_curves(curves, config~))
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}

///|
test "lisaClust: one-region clustering is supported" {
  let result = @src.lisaclust(
    lisa_test_cells(),
    config=lisa_test_config(n_clusters=1),
  ) catch {
    _ => abort("single-region lisaClust should succeed")
  }
  assert_eq(result.n_regions(), 1)
  assert_eq(result.cluster_sizes, [8])
  assert_eq(result.silhouette, 0.0)
}

///|
test "lisaClust: example data has two balanced images and four types" {
  let cells = @src.lisaclust_example_data()
  assert_eq(cells.length(), 64)
  let images : Map[String, Int] = Map([])
  let types : Map[String, Bool] = Map([])
  for cell in cells {
    images[cell.image_id] = match images.get(cell.image_id) {
      Some(count) => count + 1
      None => 1
    }
    types[cell.cell_type] = true
  }
  assert_eq(images["image_1"], 32)
  assert_eq(images["image_2"], 32)
  assert_eq(types.length(), 4)
}

///|
test "lisaClust: example data separates spatial compartments" {
  let config = @src.LisaConfig::create(
    radii=[2.5, 5.0],
    window_kind=@src.lisa_rectangle_window(),
    bandwidth=4.0,
    edge_correction=true,
    window_padding=0.1,
    edge_samples=64,
    n_clusters=2,
    n_starts=5,
    seed=11,
  ) catch {
    _ => abort("example lisaClust configuration should be valid")
  }
  let result = @src.lisaclust(@src.lisaclust_example_data(), config~) catch {
    _ => abort("example lisaClust workflow should succeed")
  }
  assert_eq(result.n_cells(), 64)
  assert_eq(result.n_regions(), 2)
  assert_true(result.silhouette > 0.0)
}

///|
test "lisaClust: SpatialExperiment adapter writes regions to a copy" {
  let experiment = lisa_test_experiment()
  let output = @src.lisaclust_spatial_experiment(
    experiment,
    cell_id_key="cellID",
    config=lisa_test_config(),
  ) catch {
    _ => abort("lisaClust SpatialExperiment integration should succeed")
  }
  assert_eq(output.result.n_cells(), 8)
  assert_eq(output.experiment.col_data.length(), 8)
  assert_true(output.experiment.col_data[0].contains("region"))
  assert_true(!experiment.col_data[0].contains("region"))
  assert_eq(output.experiment.metadata["lisaclust_regions"], "2")
  assert_eq(output.experiment.metadata["lisaclust_features"], "6")
}

///|
test "lisaClust: SpatialExperiment adapter supports generated cell IDs" {
  let output = @src.lisaclust_spatial_experiment(
    lisa_test_experiment(),
    config=lisa_test_config(),
  ) catch {
    _ => abort("generated lisaClust IDs should succeed")
  }
  assert_eq(output.result.curves.cell_ids[0], "cell_1")
  assert_eq(output.result.curves.cell_ids[7], "cell_8")
}

///|
test "lisaClust: SpatialExperiment adapter supports custom region key" {
  let output = @src.lisaclust_spatial_experiment(
    lisa_test_experiment(),
    cell_id_key="cellID",
    region_key="microenvironment",
    config=lisa_test_config(),
  ) catch {
    _ => abort("custom lisaClust region key should succeed")
  }
  assert_true(output.experiment.col_data[0].contains("microenvironment"))
  assert_true(!output.experiment.col_data[0].contains("region"))
}

///|
test "lisaClust: SpatialExperiment rejects empty col_data" {
  let failed = try {
    ignore(
      @src.lisaclust_spatial_experiment(
        @src.SpatialExperiment::new(),
        config=lisa_test_config(),
      ),
    )
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}

///|
test "lisaClust: SpatialExperiment rejects coordinate mismatch" {
  let experiment = lisa_test_experiment()
  ignore(experiment.spatial_coords.pop())
  let failed = try {
    ignore(
      @src.lisaclust_spatial_experiment(
        experiment,
        config=lisa_test_config(),
      ),
    )
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}

///|
test "lisaClust: SpatialExperiment rejects missing required columns" {
  let experiment = lisa_test_experiment()
  ignore(experiment.col_data[0].remove("cellType"))
  let failed = try {
    ignore(
      @src.lisaclust_spatial_experiment(
        experiment,
        config=lisa_test_config(),
      ),
    )
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}

///|
test "lisaClust: SpatialExperiment rejects duplicate supplied cell IDs" {
  let experiment = lisa_test_experiment()
  experiment.col_data[1]["cellID"] = experiment.col_data[0]["cellID"]
  let failed = try {
    ignore(
      @src.lisaclust_spatial_experiment(
        experiment,
        cell_id_key="cellID",
        config=lisa_test_config(),
      ),
    )
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}

///|
test "lisaClust: SpatialExperiment validates key names" {
  let failed = try {
    ignore(
      @src.lisaclust_spatial_experiment(
        lisa_test_experiment(),
        region_key="",
        config=lisa_test_config(),
      ),
    )
    false
  } catch {
    LisaError(_) => true
  }
  assert_true(failed)
}
