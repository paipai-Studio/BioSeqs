///|
fn dua_test_small_data() -> @src.DropletCountMatrix {
  @src.DropletCountMatrix::create(
    [
      [3, 2, 4, 3, 2, 3, 20, 18],
      [2, 3, 1, 2, 3, 4, 0, 0],
      [0, 0, 0, 0, 0, 0, 15, 0],
      [0, 0, 0, 0, 0, 0, 0, 15],
    ],
    feature_names=["ambient_a", "ambient_b", "marker_a", "marker_b"],
    barcodes=[
      "empty_1", "empty_2", "empty_3", "empty_4", "empty_5", "empty_6", "cell_a",
      "cell_b",
    ],
  ) catch {
    _ => abort("valid DropletUtils test data should build")
  }
}

///|
fn dua_test_config() -> @src.DropletUtilsConfig {
  @src.DropletUtilsConfig::create(
    lower=5,
    iterations=199,
    fdr_threshold=0.1,
    retain=30.0,
    rank_exclude=0,
    rank_window=0.25,
    seed=23,
  ) catch {
    _ => abort("valid DropletUtils configuration should build")
  }
}

///|
fn dua_test_result() -> @src.DropletEmptyDropsResult {
  @src.droplet_utils_empty_drops(
    dua_test_small_data(),
    config=dua_test_config(),
  ) catch {
    _ => abort("valid DropletUtils workflow should run")
  }
}

///|
fn dua_test_sce() -> @src.SingleCellExperiment {
  let data = dua_test_small_data()
  let assay : Array[Array[Double]] = []
  for feature in 0..<data.n_features() {
    let row : Array[Double] = []
    for barcode in 0..<data.n_barcodes() {
      row.push(data.counts[feature][barcode].to_double())
    }
    assay.push(row)
  }
  @src.SingleCellExperiment::new(assay, data.feature_names, data.barcodes)
}

///|
test "DropletUtils advanced constructor creates default names" {
  let data = @src.DropletCountMatrix::create([[1, 2], [3, 4]]) catch {
    _ => abort("valid matrix should build")
  }
  assert_eq(data.feature_names, ["feature1", "feature2"])
  assert_eq(data.barcodes, ["barcode1", "barcode2"])
  assert_eq(data.n_features(), 2)
  assert_eq(data.n_barcodes(), 2)
}

///|
test "DropletUtils advanced constructor copies matrix and names" {
  let counts = [[1, 2], [3, 4]]
  let features = ["g1", "g2"]
  let barcodes = ["b1", "b2"]
  let data = @src.DropletCountMatrix::create(
    counts,
    feature_names=features,
    barcodes~,
  ) catch {
    _ => abort("valid matrix should build")
  }
  counts[0][0] = 99
  features[0] = "changed"
  barcodes[0] = "changed"
  assert_eq(data.counts[0][0], 1)
  assert_eq(data.feature_names[0], "g1")
  assert_eq(data.barcodes[0], "b1")
}

///|
test "DropletUtils advanced copy_counts is defensive" {
  let data = dua_test_small_data()
  let copied = data.copy_counts()
  copied[0][0] = 999
  assert_eq(data.counts[0][0], 3)
}

///|
test "DropletUtils advanced computes feature and barcode totals" {
  let data = dua_test_small_data()
  assert_eq(data.barcode_totals(), [5, 5, 5, 5, 5, 7, 35, 33])
  assert_eq(data.feature_totals(), [55, 15, 15, 15])
}

///|
test "DropletUtils advanced rejects an empty matrix" {
  let raised = try {
    ignore(@src.DropletCountMatrix::create([]))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "DropletUtils advanced rejects an empty barcode dimension" {
  let raised = try {
    ignore(@src.DropletCountMatrix::create([[]]))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "DropletUtils advanced rejects ragged matrices" {
  let raised = try {
    ignore(@src.DropletCountMatrix::create([[1, 2], [3]]))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "DropletUtils advanced rejects negative counts" {
  let raised = try {
    ignore(@src.DropletCountMatrix::create([[1, -1], [3, 4]]))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "DropletUtils advanced validates identifier lengths" {
  let raised = try {
    ignore(
      @src.DropletCountMatrix::create([[1, 2], [3, 4]], feature_names=["g1"], barcodes=[
        "b1", "b2",
      ]),
    )
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "DropletUtils advanced rejects empty and duplicate identifiers" {
  let empty = try {
    ignore(
      @src.DropletCountMatrix::create([[1, 2], [3, 4]], feature_names=["g1", ""]),
    )
    false
  } catch {
    _ => true
  }
  let duplicate = try {
    ignore(
      @src.DropletCountMatrix::create([[1, 2], [3, 4]], barcodes=["b1", "b1"]),
    )
    false
  } catch {
    _ => true
  }
  assert_true(empty && duplicate)
}

///|
test "DropletUtils advanced default configuration matches portable workflow" {
  let config = @src.DropletUtilsConfig::default()
  assert_eq(config.lower, 100)
  assert_eq(config.iterations, 10000)
  assert_eq(config.fdr_threshold, 0.001)
  assert_eq(config.retain, -1.0)
  assert_eq(config.by_rank, 0)
  assert_eq(config.ignore, -1)
  assert_eq(config.alpha, -1.0)
  assert_eq(config.estimate_alpha, false)
  assert_eq(config.seed, 1)
}

///|
test "DropletUtils advanced configuration preserves explicit options" {
  let config = @src.DropletUtilsConfig::create(
    lower=12,
    iterations=75,
    fdr_threshold=0.2,
    retain=50.0,
    by_rank=3,
    ignore=2,
    test_ambient=true,
    alpha=4.0,
    estimate_alpha=true,
    rank_exclude=2,
    rank_window=0.5,
    gradient_threshold=-0.5,
    round_counts=false,
    seed=19,
  ) catch {
    _ => abort("valid explicit configuration should build")
  }
  assert_eq(config.lower, 12)
  assert_eq(config.alpha, 4.0)
  assert_true(config.test_ambient && config.estimate_alpha)
  assert_eq(config.round_counts, false)
  assert_eq(config.seed, 19)
}

///|
test "DropletUtils advanced rejects invalid count thresholds" {
  let lower = try {
    ignore(@src.DropletUtilsConfig::create(lower=-1))
    false
  } catch {
    _ => true
  }
  let iterations = try {
    ignore(@src.DropletUtilsConfig::create(iterations=0))
    false
  } catch {
    _ => true
  }
  let by_rank = try {
    ignore(@src.DropletUtilsConfig::create(by_rank=-1))
    false
  } catch {
    _ => true
  }
  let ignore_value = try {
    ignore(@src.DropletUtilsConfig::create(ignore=-2))
    false
  } catch {
    _ => true
  }
  assert_true(lower && iterations && by_rank && ignore_value)
}

///|
test "DropletUtils advanced rejects invalid probability configuration" {
  let fdr = try {
    ignore(@src.DropletUtilsConfig::create(fdr_threshold=1.1))
    false
  } catch {
    _ => true
  }
  let retain = try {
    ignore(@src.DropletUtilsConfig::create(retain=0.0))
    false
  } catch {
    _ => true
  }
  let alpha = try {
    ignore(@src.DropletUtilsConfig::create(alpha=0.0))
    false
  } catch {
    _ => true
  }
  let seed = try {
    ignore(@src.DropletUtilsConfig::create(seed=0))
    false
  } catch {
    _ => true
  }
  assert_true(fdr && retain && alpha && seed)
}

///|
test "DropletUtils advanced Good-Turing profile is normalized and positive" {
  let profile = @src.droplet_utils_good_turing([1, 1, 2, 3, 0, 0]) catch {
    _ => abort("Good-Turing profile should fit")
  }
  let mut total = 0.0
  for value in profile {
    assert_true(value > 0.0)
    total = total + value
  }
  assert_true((total - 1.0).abs() < 1.0e-10)
  assert_true((profile[0] - profile[1]).abs() < 1.0e-12)
  assert_true((profile[4] - profile[5]).abs() < 1.0e-12)
}

///|
test "DropletUtils advanced Good-Turing protects unobserved active features" {
  let profile = @src.droplet_utils_good_turing([17, 15, 0, 0]) catch {
    _ => abort("safe Good-Turing profile should fit")
  }
  assert_true(profile[2] > 0.0)
  assert_true(profile[3] > 0.0)
  assert_true((profile[2] - profile[3]).abs() < 1.0e-12)
}

///|
test "DropletUtils advanced Good-Turing rejects invalid inputs" {
  let empty = try {
    ignore(@src.droplet_utils_good_turing([]))
    false
  } catch {
    _ => true
  }
  let all_zero = try {
    ignore(@src.droplet_utils_good_turing([0, 0]))
    false
  } catch {
    _ => true
  }
  let negative = try {
    ignore(@src.droplet_utils_good_turing([1, -1]))
    false
  } catch {
    _ => true
  }
  assert_true(empty && all_zero && negative)
}

///|
test "DropletUtils advanced ambient profile uses lower threshold" {
  let ambient = @src.droplet_utils_ambient_profile(
    dua_test_small_data(),
    lower=5,
  ) catch {
    _ => abort("ambient profile should fit")
  }
  assert_eq(ambient.assumed_empty, [
    true, true, true, true, true, false, false, false,
  ])
  assert_eq(ambient.n_empty(), 5)
  assert_eq(ambient.total, 25)
  assert_eq(ambient.raw_counts, [14, 11, 0, 0])
}

///|
test "DropletUtils advanced ambient profile honors known-empty mask" {
  let mask = [true, false, true, false, true, false, false, false]
  let ambient = @src.droplet_utils_ambient_profile(
    dua_test_small_data(),
    lower=0,
    known_empty=mask,
  ) catch {
    _ => abort("known-empty ambient profile should fit")
  }
  assert_eq(ambient.assumed_empty, mask)
  assert_eq(ambient.n_empty(), 3)
  assert_eq(ambient.total, 15)
}

///|
test "DropletUtils advanced ambient profile supports by-rank selection" {
  let data = @src.DropletCountMatrix::create([[50, 40, 30, 20, 10]], barcodes=[
    "b1", "b2", "b3", "b4", "b5",
  ]) catch {
    _ => abort("valid rank data should build")
  }
  let ambient = @src.droplet_utils_ambient_profile(data, by_rank=2) catch {
    _ => abort("by-rank ambient profile should fit")
  }
  assert_eq(ambient.assumed_empty, [false, false, true, true, true])
  assert_eq(ambient.lower, 30)
  assert_eq(ambient.total, 60)
}

///|
test "DropletUtils advanced ambient profile rejects invalid masks" {
  let wrong_length = try {
    ignore(
      @src.droplet_utils_ambient_profile(dua_test_small_data(), known_empty=[
        true,
      ]),
    )
    false
  } catch {
    _ => true
  }
  let none_selected = try {
    ignore(
      @src.droplet_utils_ambient_profile(
        dua_test_small_data(),
        known_empty=Array::make(8, false),
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(wrong_length && none_selected)
}

///|
test "DropletUtils advanced ambient profile rejects empty ambient counts" {
  let raised = try {
    ignore(@src.droplet_utils_ambient_profile(dua_test_small_data(), lower=0))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "DropletUtils advanced barcode ranks average ties" {
  let data = @src.DropletCountMatrix::create([[100, 100, 60, 60, 30, 15, 7, 3]]) catch {
    _ => abort("valid rank data should build")
  }
  let ranks = @src.droplet_utils_barcode_ranks(
    data,
    lower=0,
    exclude_from=0,
    window=0.2,
  ) catch {
    _ => abort("barcode ranks should fit")
  }
  assert_eq(ranks.ranks[0], 1.5)
  assert_eq(ranks.ranks[1], 1.5)
  assert_eq(ranks.ranks[2], 3.5)
  assert_eq(ranks.ranks[3], 3.5)
}

///|
test "DropletUtils advanced barcode ranks return finite knees" {
  let data = @src.DropletCountMatrix::create([
    [1000, 800, 600, 350, 180, 90, 45, 20, 8, 3],
  ]) catch {
    _ => abort("valid rank curve should build")
  }
  let ranks = @src.droplet_utils_barcode_ranks(
    data,
    lower=0,
    exclude_from=0,
    window=0.25,
  ) catch {
    _ => abort("barcode-rank curve should fit")
  }
  assert_true(ranks.knee > 0.0 && ranks.knee <= 1000.0)
  assert_true(ranks.inflection > 0.0 && ranks.inflection <= 1000.0)
  assert_eq(ranks.totals[0], 1000)
}

///|
test "DropletUtils advanced barcode ranks reject insufficient curves" {
  let data = @src.DropletCountMatrix::create([[10, 10, 10]]) catch {
    _ => abort("valid tied rank data should build")
  }
  let raised = try {
    ignore(
      @src.droplet_utils_barcode_ranks(
        data,
        lower=0,
        exclude_from=0,
        window=0.2,
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "DropletUtils advanced multinomial log probability is exact" {
  let value = @src.droplet_utils_log_probability([2, 1], [0.25, 0.75]) catch {
    _ => abort("multinomial probability should compute")
  }
  assert_true((@math.exp(value) - 0.140625).abs() < 1.0e-10)
}

///|
test "DropletUtils advanced Dirichlet-multinomial probability is exact" {
  let value = @src.droplet_utils_log_probability([1, 1], [0.5, 0.5], alpha=2.0) catch {
    _ => abort("Dirichlet-multinomial probability should compute")
  }
  assert_true((@math.exp(value) - 1.0 / 3.0).abs() < 1.0e-10)
}

///|
test "DropletUtils advanced multinomial and overdispersed models differ" {
  let multinomial = @src.droplet_utils_log_probability([5, 0], [0.5, 0.5]) catch {
    _ => abort("multinomial probability should compute")
  }
  let overdispersed = @src.droplet_utils_log_probability(
    [5, 0],
    [0.5, 0.5],
    alpha=2.0,
  ) catch {
    _ => abort("Dirichlet-multinomial probability should compute")
  }
  assert_true(overdispersed > multinomial)
}

///|
test "DropletUtils advanced probability validates dimensions and counts" {
  let dimensions = try {
    ignore(@src.droplet_utils_log_probability([1], [0.5, 0.5]))
    false
  } catch {
    _ => true
  }
  let negative = try {
    ignore(@src.droplet_utils_log_probability([-1, 2], [0.5, 0.5]))
    false
  } catch {
    _ => true
  }
  assert_true(dimensions && negative)
}

///|
test "DropletUtils advanced probability validates ambient proportions" {
  let sum = try {
    ignore(@src.droplet_utils_log_probability([1, 1], [0.4, 0.4]))
    false
  } catch {
    _ => true
  }
  let negative = try {
    ignore(@src.droplet_utils_log_probability([1, 1], [-0.1, 1.1]))
    false
  } catch {
    _ => true
  }
  let impossible = try {
    ignore(@src.droplet_utils_log_probability([1, 0], [0.0, 1.0]))
    false
  } catch {
    _ => true
  }
  assert_true(sum && negative && impossible)
}

///|
test "DropletUtils advanced Monte Carlo p-value is reproducible" {
  let first = @src.droplet_utils_monte_carlo_p_value(
    [8, 2],
    [0.5, 0.5],
    iterations=199,
    seed=17,
  ) catch {
    _ => abort("Monte Carlo test should run")
  }
  let second = @src.droplet_utils_monte_carlo_p_value(
    [8, 2],
    [0.5, 0.5],
    iterations=199,
    seed=17,
  ) catch {
    _ => abort("Monte Carlo test should repeat")
  }
  assert_eq(first, second)
  assert_true(first.0 >= 1.0 / 200.0 && first.0 <= 1.0)
}

///|
test "DropletUtils advanced Monte Carlo uses Phipson-Smyth correction" {
  let result = @src.droplet_utils_monte_carlo_p_value(
    [30, 0],
    [0.5, 0.5],
    iterations=99,
    seed=5,
  ) catch {
    _ => abort("extreme Monte Carlo test should run")
  }
  assert_eq(result.1, true)
  assert_true((result.0 - 0.01).abs() < 1.0e-12)
}

///|
test "DropletUtils advanced Monte Carlo supports Dirichlet-multinomial null" {
  let result = @src.droplet_utils_monte_carlo_p_value(
    [7, 3],
    [0.5, 0.5],
    iterations=99,
    alpha=3.0,
    seed=7,
  ) catch {
    _ => abort("Dirichlet-multinomial Monte Carlo test should run")
  }
  assert_true(result.0 >= 0.01 && result.0 <= 1.0)
}

///|
test "DropletUtils advanced Monte Carlo rejects invalid controls" {
  let iterations = try {
    ignore(
      @src.droplet_utils_monte_carlo_p_value([1, 1], [0.5, 0.5], iterations=0),
    )
    false
  } catch {
    _ => true
  }
  let zero_total = try {
    ignore(@src.droplet_utils_monte_carlo_p_value([0, 0], [0.5, 0.5]))
    false
  } catch {
    _ => true
  }
  assert_true(iterations && zero_total)
}

///|
test "DropletUtils advanced alpha estimation is finite and positive" {
  let data = dua_test_small_data()
  let ambient = @src.droplet_utils_ambient_profile(data, lower=5) catch {
    _ => abort("ambient profile should fit")
  }
  let alpha = @src.droplet_utils_estimate_alpha(data, ambient) catch {
    _ => abort("alpha should fit")
  }
  assert_true(alpha >= 0.01 && alpha <= 10000.0)
  assert_true(alpha == alpha && alpha.abs() < 1.0e300)
}

///|
test "DropletUtils advanced alpha estimation validates dimensions" {
  let ambient = @src.droplet_utils_ambient_profile(
    dua_test_small_data(),
    lower=5,
  ) catch {
    _ => abort("ambient profile should fit")
  }
  let mismatched = @src.DropletCountMatrix::create([[1, 2]]) catch {
    _ => abort("mismatched matrix should build")
  }
  let raised = try {
    ignore(@src.droplet_utils_estimate_alpha(mismatched, ambient))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "DropletUtils advanced emptyDrops preserves upstream NA semantics" {
  let result = dua_test_result()
  for barcode in 0..<5 {
    assert_true(result.p_values[barcode] is None)
    assert_true(result.log_probabilities[barcode] is None)
    assert_true(result.limited[barcode] is None)
    assert_true(result.fdr[barcode] is None)
  }
  assert_true(result.p_values[5] is Some(_))
  assert_true(result.p_values[6] is Some(_))
  assert_true(result.p_values[7] is Some(_))
}

///|
test "DropletUtils advanced emptyDrops retains high-count barcodes" {
  let result = dua_test_result()
  assert_eq(result.always_retained, [
    false, false, false, false, false, false, true, true,
  ])
  assert_eq(result.calls[6], true)
  assert_eq(result.calls[7], true)
  match result.fdr[6] {
    Some(value) => assert_eq(value, 0.0)
    None => abort("retained barcode should have FDR")
  }
}

///|
test "DropletUtils advanced emptyDrops records Limited flags" {
  let result = dua_test_result()
  match result.limited[6] {
    Some(value) => assert_eq(value, true)
    None => abort("tested barcode should have Limited flag")
  }
  match result.p_values[6] {
    Some(value) => assert_true(value >= 1.0 / 200.0 && value <= 1.0)
    None => abort("tested barcode should have p-value")
  }
}

///|
test "DropletUtils advanced emptyDrops reports model diagnostics" {
  let result = dua_test_result()
  assert_eq(result.lower, 5)
  assert_eq(result.retain, 30.0)
  assert_eq(result.alpha, -1.0)
  assert_eq(result.alpha_estimated, false)
  assert_eq(result.iterations, 199)
  assert_eq(result.fdr_threshold, 0.1)
  assert_eq(result.ambient.n_empty(), 5)
}

///|
test "DropletUtils advanced emptyDrops supports ambient diagnostics" {
  let config = @src.DropletUtilsConfig::create(
    lower=5,
    iterations=49,
    fdr_threshold=0.2,
    retain=30.0,
    test_ambient=true,
    seed=11,
  ) catch {
    _ => abort("diagnostic configuration should build")
  }
  let result = @src.droplet_utils_empty_drops(dua_test_small_data(), config~) catch {
    _ => abort("diagnostic workflow should run")
  }
  assert_eq(result.n_tested(), 8)
  assert_true(result.p_values[0] is Some(_))
  assert_true(result.fdr[0] is None)
  assert_eq(result.calls[0], false)
}

///|
test "DropletUtils advanced emptyDrops supports ignore threshold" {
  let config = @src.DropletUtilsConfig::create(
    lower=5,
    iterations=49,
    fdr_threshold=0.2,
    retain=100.0,
    ignore=33,
    seed=11,
  ) catch {
    _ => abort("ignore configuration should build")
  }
  let result = @src.droplet_utils_empty_drops(dua_test_small_data(), config~) catch {
    _ => abort("ignore workflow should run")
  }
  assert_eq(result.n_tested(), 1)
  assert_true(result.p_values[6] is Some(_))
  assert_true(result.p_values[7] is None)
}

///|
test "DropletUtils advanced emptyDrops supports alpha estimation" {
  let config = @src.DropletUtilsConfig::create(
    lower=5,
    iterations=29,
    fdr_threshold=0.2,
    retain=30.0,
    estimate_alpha=true,
    seed=13,
  ) catch {
    _ => abort("estimated-alpha configuration should build")
  }
  let result = @src.droplet_utils_empty_drops(dua_test_small_data(), config~) catch {
    _ => abort("estimated-alpha workflow should run")
  }
  assert_eq(result.alpha_estimated, true)
  assert_true(result.alpha > 0.0)
}

///|
test "DropletUtils advanced emptyDrops supports known-empty droplets" {
  let mask = [true, true, true, true, true, false, false, false]
  let result = @src.droplet_utils_empty_drops(
    dua_test_small_data(),
    known_empty=mask,
    config=dua_test_config(),
  ) catch {
    _ => abort("known-empty workflow should run")
  }
  assert_eq(result.ambient.assumed_empty, mask)
  assert_eq(result.ambient.n_empty(), 5)
}

///|
test "DropletUtils advanced emptyDrops computes automatic knee retention" {
  let data = @src.droplet_utils_advanced_example() catch {
    _ => abort("advanced example should build")
  }
  let config = @src.DropletUtilsConfig::create(
    lower=15,
    iterations=19,
    fdr_threshold=0.2,
    rank_exclude=0,
    rank_window=0.2,
    seed=3,
  ) catch {
    _ => abort("automatic-retain configuration should build")
  }
  let result = @src.droplet_utils_empty_drops(data, config~) catch {
    _ => abort("automatic-retain workflow should run")
  }
  assert_true(result.retain > 0.0)
  let mut retained = 0
  for value in result.always_retained {
    if value {
      retained = retained + 1
    }
  }
  assert_true(retained > 0 && retained < data.n_barcodes())
}

///|
test "DropletUtils advanced result query helpers are consistent" {
  let result = dua_test_result()
  assert_eq(result.n_called(), 2)
  assert_eq(result.n_tested(), 3)
  assert_eq(result.index_of("cell_a"), Some(6))
  assert_eq(result.index_of("missing"), None)
  assert_eq(result.called_barcodes(), ["cell_a", "cell_b"])
  assert_true(result.summary().has_prefix("DropletUtils emptyDrops:"))
}

///|
test "DropletUtils advanced result filters called barcodes" {
  let filtered = dua_test_result().filter(dua_test_small_data()) catch {
    _ => abort("called barcodes should filter")
  }
  assert_eq(filtered.n_features(), 4)
  assert_eq(filtered.n_barcodes(), 2)
  assert_eq(filtered.barcodes, ["cell_a", "cell_b"])
  assert_eq(filtered.barcode_totals(), [35, 33])
}

///|
test "DropletUtils advanced result rejects mismatched filtering data" {
  let other = @src.DropletCountMatrix::create(
    [[1, 2], [3, 4]],
    feature_names=["g1", "g2"],
    barcodes=["x", "y"],
  ) catch {
    _ => abort("other matrix should build")
  }
  let raised = try {
    ignore(dua_test_result().filter(other))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "DropletUtils advanced example has ambient and cell populations" {
  let data = @src.droplet_utils_advanced_example() catch {
    _ => abort("advanced example should build")
  }
  assert_eq(data.n_features(), 12)
  assert_eq(data.n_barcodes(), 64)
  let totals = data.barcode_totals()
  assert_true(totals[0] < 20)
  assert_true(totals[40] > 80)
  assert_true(totals[52] > 80)
}

///|
test "DropletUtils advanced SingleCellExperiment writes result columns" {
  let output = @src.droplet_utils_empty_drops_sce(
    dua_test_sce(),
    output_prefix="droplet",
    config=dua_test_config(),
  ) catch {
    _ => abort("SingleCellExperiment workflow should run")
  }
  assert_true(output.experiment.col_data.contains("droplet.total"))
  assert_true(output.experiment.col_data.contains("droplet.pValue"))
  assert_true(output.experiment.col_data.contains("droplet.limited"))
  assert_true(output.experiment.col_data.contains("droplet.fdr"))
  assert_true(output.experiment.col_data.contains("droplet.class"))
  assert_true(output.experiment.col_data.contains("droplet.retained"))
  assert_eq(output.experiment.col_data["droplet.class"][6], "cell")
}

///|
test "DropletUtils advanced SingleCellExperiment writes ambient rows and metadata" {
  let output = @src.droplet_utils_empty_drops_sce(
    dua_test_sce(),
    output_prefix="droplet",
    config=dua_test_config(),
  ) catch {
    _ => abort("SingleCellExperiment workflow should run")
  }
  assert_eq(output.experiment.row_data["droplet.ambient"].length(), 4)
  assert_eq(output.experiment.row_data["droplet.ambientCount"].length(), 4)
  assert_eq(output.experiment.metadata["droplet.lower"], "5")
  assert_eq(output.experiment.metadata["droplet.called"], "2")
}

///|
test "DropletUtils advanced SingleCellExperiment does not mutate input" {
  let input = dua_test_sce()
  let output = @src.droplet_utils_empty_drops_sce(
    input,
    output_prefix="droplet",
    config=dua_test_config(),
  ) catch {
    _ => abort("SingleCellExperiment workflow should run")
  }
  assert_eq(input.col_data.contains("droplet.class"), false)
  output.experiment.assays["counts"][0][0] = 999.0
  assert_eq(input.assays["counts"][0][0], 3.0)
}

///|
test "DropletUtils advanced SingleCellExperiment parses known-empty column" {
  let input = dua_test_sce()
  input.col_data["known"] = [
    "empty", "empty", "empty", "empty", "empty", "cell", "cell", "cell",
  ]
  let output = @src.droplet_utils_empty_drops_sce(
    input,
    known_empty_column="known",
    config=dua_test_config(),
  ) catch {
    _ => abort("known-empty SingleCellExperiment workflow should run")
  }
  assert_eq(output.result.ambient.n_empty(), 5)
  assert_eq(output.result.ambient.assumed_empty[5], false)
}

///|
test "DropletUtils advanced SingleCellExperiment rounds counts by default" {
  let input = dua_test_sce()
  input.assays["counts"][0][0] = 3.4
  let output = @src.droplet_utils_empty_drops_sce(
    input,
    config=dua_test_config(),
  ) catch {
    _ => abort("rounded SingleCellExperiment workflow should run")
  }
  assert_eq(output.result.totals[0], 5)
}

///|
test "DropletUtils advanced SingleCellExperiment can require integer counts" {
  let input = dua_test_sce()
  input.assays["counts"][0][0] = 3.4
  let config = @src.DropletUtilsConfig::create(
    lower=5,
    iterations=19,
    fdr_threshold=0.2,
    retain=30.0,
    round_counts=false,
  ) catch {
    _ => abort("strict integer configuration should build")
  }
  let raised = try {
    ignore(@src.droplet_utils_empty_drops_sce(input, config~))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "DropletUtils advanced SingleCellExperiment rejects missing assay" {
  let raised = try {
    ignore(
      @src.droplet_utils_empty_drops_sce(
        dua_test_sce(),
        assay_name="missing",
        config=dua_test_config(),
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "DropletUtils advanced SingleCellExperiment validates known-empty values" {
  let input = dua_test_sce()
  input.col_data["known"] = [
    "empty", "empty", "empty", "empty", "empty", "unknown", "cell", "cell",
  ]
  let raised = try {
    ignore(
      @src.droplet_utils_empty_drops_sce(
        input,
        known_empty_column="known",
        config=dua_test_config(),
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "DropletUtils advanced SingleCellExperiment validates assay values" {
  let negative = dua_test_sce()
  negative.assays["counts"][0][0] = -1.0
  let negative_raised = try {
    ignore(
      @src.droplet_utils_empty_drops_sce(negative, config=dua_test_config()),
    )
    false
  } catch {
    _ => true
  }
  let non_finite = dua_test_sce()
  non_finite.assays["counts"][0][0] = @double.not_a_number
  let finite_raised = try {
    ignore(
      @src.droplet_utils_empty_drops_sce(non_finite, config=dua_test_config()),
    )
    false
  } catch {
    _ => true
  }
  assert_true(negative_raised && finite_raised)
}
