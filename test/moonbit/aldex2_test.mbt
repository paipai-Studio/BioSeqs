// Tests for the Bioconductor ALDEx2-inspired compositional inference module.

///|
fn aldex2_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "ALDEx2 value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn aldex2_test_config() -> @src.Aldex2Config {
  @src.Aldex2Config::create(
    mc_samples=32,
    denominator=@src.Aldex2Iqlr,
    seed=2026,
  ) catch {
    _ => abort("ALDEx2 test configuration should be valid")
  }
}

///|
fn aldex2_test_clr() -> @src.Aldex2Clr {
  let (counts, conditions, features, samples) = @src.aldex2_example_data()
  let config = aldex2_test_config()
  @src.aldex2_clr(
    counts,
    conditions,
    config~,
    feature_names=features,
    sample_names=samples,
  ) catch {
    Aldex2Error(message) => abort("ALDEx2 CLR should succeed: " + message)
  }
}

///|
fn aldex2_test_result() -> @src.Aldex2Result {
  let (counts, conditions, features, samples) = @src.aldex2_example_data()
  let config = aldex2_test_config()
  @src.aldex2(
    counts,
    conditions,
    config~,
    feature_names=features,
    sample_names=samples,
  ) catch {
    Aldex2Error(message) => abort("ALDEx2 example should succeed: " + message)
  }
}

///|
fn aldex2_test_effect(
  result : @src.Aldex2Result,
  name : String,
) -> @src.Aldex2EffectFeature {
  match result.effects.find_feature(name) {
    Some(value) => value
    None => abort("expected ALDEx2 effect feature: " + name)
  }
}

///|
fn aldex2_test_statistics(
  result : @src.Aldex2Result,
  name : String,
) -> @src.Aldex2TestFeature {
  match result.tests.find_feature(name) {
    Some(value) => value
    None => abort("expected ALDEx2 test feature: " + name)
  }
}

///|
fn aldex2_test_unwrap(values : Array[Double]?) -> Array[Double] {
  match values {
    Some(value) => value
    None => abort("expected an ALDEx2 posterior slice")
  }
}

///|
fn aldex2_test_median(values : Array[Double]) -> Double {
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
  if sorted.length() % 2 == 1 {
    sorted[sorted.length() / 2]
  } else {
    (sorted[sorted.length() / 2 - 1] + sorted[sorted.length() / 2]) / 2.0
  }
}

///|
fn aldex2_test_scale(
  sample_count : Int,
  mc_samples : Int,
  value : Double,
) -> Array[Array[Double]] {
  let scales : Array[Array[Double]] = []
  for _ in 0..<sample_count {
    scales.push(Array::make(mc_samples, value))
  }
  scales
}

///|
fn aldex2_test_experiment() -> @src.SummarizedExperiment {
  let (counts, _, _, _) = @src.aldex2_example_data()
  let assays : Map[String, Array[Array[Double]]] = Map([("counts", counts)])
  let metadata : Map[String, String] = Map([("source", "aldex2_test")])
  @src.summarized_experiment(assays, [], [], metadata)
}

///|
test "aldex2: default configuration matches upstream controls" {
  let config = @src.Aldex2Config::default()
  assert_eq(config.mc_samples, 128)
  assert_eq(config.prior, 0.5)
  assert_eq(config.denominator, @src.Aldex2All)
  assert_eq(config.denominator_indices, [])
  assert_eq(config.seed, 42)
  assert_false(config.paired)
  assert_eq(config.interval_level, 0.95)
}

///|
test "aldex2: custom configuration preserves controls" {
  let config = @src.Aldex2Config::create(
    mc_samples=16,
    prior=1.0,
    denominator=@src.Aldex2User,
    denominator_indices=[1, 3],
    seed=99,
    paired=true,
    interval_level=0.9,
  ) catch {
    _ => abort("custom ALDEx2 configuration should be valid")
  }
  assert_eq(config.mc_samples, 16)
  assert_eq(config.prior, 1.0)
  assert_eq(config.denominator_indices, [1, 3])
  assert_eq(config.seed, 99)
  assert_true(config.paired)
  assert_eq(config.interval_level, 0.9)
}

///|
test "aldex2: configuration rejects nonpositive Monte Carlo count" {
  let raised = try {
    ignore(@src.Aldex2Config::create(mc_samples=0))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: configuration rejects nonpositive prior" {
  let zero = try {
    ignore(@src.Aldex2Config::create(prior=0.0))
    false
  } catch {
    Aldex2Error(_) => true
  }
  let negative = try {
    ignore(@src.Aldex2Config::create(prior=-0.5))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(zero)
  assert_true(negative)
}

///|
test "aldex2: configuration rejects nonpositive seed" {
  let raised = try {
    ignore(@src.Aldex2Config::create(seed=0))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
  let config = @src.Aldex2Config::create(mc_samples=4, seed=2147483647) catch {
    _ => abort("maximum Park-Miller seed should be normalized")
  }
  let clr = @src.aldex2_clr(
    [[1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 4.0, 5.0]],
    ["A", "A", "B", "B"],
    config~,
  ) catch {
    _ => abort("normalized Park-Miller seed should support sampling")
  }
  assert_eq(clr.mc_sample_count(), 4)
}

///|
test "aldex2: configuration rejects invalid intervals" {
  let zero = try {
    ignore(@src.Aldex2Config::create(interval_level=0.0))
    false
  } catch {
    Aldex2Error(_) => true
  }
  let one = try {
    ignore(@src.Aldex2Config::create(interval_level=1.0))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(zero)
  assert_true(one)
}

///|
test "aldex2: user denominator requires feature indices" {
  let raised = try {
    ignore(@src.Aldex2Config::create(denominator=@src.Aldex2User))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: denominator indices require user mode" {
  let raised = try {
    ignore(
      @src.Aldex2Config::create(denominator=@src.Aldex2All, denominator_indices=[
        0,
      ]),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: example data uses feature by sample orientation" {
  let (counts, conditions, features, samples) = @src.aldex2_example_data()
  assert_eq(counts.length(), 12)
  assert_eq(counts[0].length(), 8)
  assert_eq(conditions.length(), 8)
  assert_eq(features.length(), 12)
  assert_eq(samples.length(), 8)
}

///|
test "aldex2: example labels identify intended signals" {
  let (_, conditions, features, samples) = @src.aldex2_example_data()
  assert_eq(features[0], "increased")
  assert_eq(features[1], "decreased")
  assert_eq(features[11], "all_zero")
  assert_eq(conditions[0], "control")
  assert_eq(conditions[4], "treated")
  assert_eq(samples[0], "C1")
  assert_eq(samples[7], "T4")
}

///|
test "aldex2: CLR dimensions match retained data" {
  let clr = aldex2_test_clr()
  assert_eq(clr.feature_count(), 11)
  assert_eq(clr.sample_count(), 8)
  assert_eq(clr.mc_sample_count(), 32)
}

///|
test "aldex2: all-zero features are removed" {
  let clr = aldex2_test_clr()
  assert_eq(clr.removed_feature_count(), 1)
  assert_eq(clr.original_feature_count, 12)
  assert_eq(clr.original_feature_indices.length(), 11)
  assert_eq(clr.original_feature_indices[10], 10)
}

///|
test "aldex2: supplied feature and sample names are retained" {
  let clr = aldex2_test_clr()
  assert_eq(clr.feature_names[0], "increased")
  assert_eq(clr.feature_names[10], "stable_6")
  assert_eq(clr.sample_names[0], "C1")
  assert_eq(clr.sample_names[7], "T4")
}

///|
test "aldex2: omitted names are generated deterministically" {
  let clr = @src.aldex2_clr(
    [[10.0, 11.0, 20.0, 21.0], [5.0, 6.0, 7.0, 8.0]],
    ["A", "A", "B", "B"],
    config=@src.Aldex2Config::create(mc_samples=4) catch {
      _ => abort("small ALDEx2 configuration should succeed")
    },
  ) catch {
    _ => abort("small ALDEx2 CLR should succeed")
  }
  assert_eq(clr.feature_names, ["feature_1", "feature_2"])
  assert_eq(clr.sample_names, ["sample_1", "sample_2", "sample_3", "sample_4"])
}

///|
test "aldex2: empty feature matrix is rejected" {
  let raised = try {
    ignore(@src.aldex2_clr([], ["A", "A", "B", "B"]))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: empty sample matrix is rejected" {
  let raised = try {
    ignore(@src.aldex2_clr([[]], []))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: condition count must match samples" {
  let raised = try {
    ignore(@src.aldex2_clr([[1.0, 2.0]], ["A"]))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: ragged count matrices are rejected" {
  let raised = try {
    ignore(
      @src.aldex2_clr([[1.0, 2.0, 3.0, 4.0], [1.0, 2.0]], ["A", "A", "B", "B"]),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: negative counts are rejected" {
  let raised = try {
    ignore(@src.aldex2_clr([[1.0, -1.0, 3.0, 4.0]], ["A", "A", "B", "B"]))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: fractional counts are rejected" {
  let raised = try {
    ignore(@src.aldex2_clr([[1.0, 2.5, 3.0, 4.0]], ["A", "A", "B", "B"]))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: nonfinite counts are rejected" {
  let nonfinite = 1.0 / 0.0
  let raised = try {
    ignore(@src.aldex2_clr([[1.0, nonfinite, 3.0, 4.0]], ["A", "A", "B", "B"]))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: all-zero data are rejected" {
  let raised = try {
    ignore(@src.aldex2_clr([[0.0, 0.0, 0.0, 0.0]], ["A", "A", "B", "B"]))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: empty sample libraries are rejected" {
  let raised = try {
    ignore(
      @src.aldex2_clr([[1.0, 0.0, 2.0, 3.0], [1.0, 0.0, 2.0, 3.0]], [
        "A", "A", "B", "B",
      ]),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: feature name count is validated" {
  let raised = try {
    ignore(
      @src.aldex2_clr([[1.0, 2.0, 3.0, 4.0]], ["A", "A", "B", "B"], feature_names=[
        "one", "two",
      ]),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: sample name count is validated" {
  let raised = try {
    ignore(
      @src.aldex2_clr([[1.0, 2.0, 3.0, 4.0]], ["A", "A", "B", "B"], sample_names=[
        "one",
      ]),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: empty names are rejected" {
  let feature = try {
    ignore(
      @src.aldex2_clr([[1.0, 2.0, 3.0, 4.0]], ["A", "A", "B", "B"], feature_names=[
        "",
      ]),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  let sample = try {
    ignore(
      @src.aldex2_clr([[1.0, 2.0, 3.0, 4.0]], ["A", "A", "B", "B"], sample_names=[
        "s1", "s2", " ", "s4",
      ]),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(feature)
  assert_true(sample)
}

///|
test "aldex2: duplicate names are rejected" {
  let feature = try {
    ignore(
      @src.aldex2_clr(
        [[1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 4.0, 5.0]],
        ["A", "A", "B", "B"],
        feature_names=["same", "same"],
      ),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  let sample = try {
    ignore(
      @src.aldex2_clr([[1.0, 2.0, 3.0, 4.0]], ["A", "A", "B", "B"], sample_names=[
        "s", "s", "t1", "t2",
      ]),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(feature)
  assert_true(sample)
}

///|
test "aldex2: empty condition labels are rejected" {
  let raised = try {
    ignore(@src.aldex2_clr([[1.0, 2.0, 3.0, 4.0]], ["A", " ", "B", "B"]))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: at least two conditions are required" {
  let raised = try {
    ignore(@src.aldex2_clr([[1.0, 2.0, 3.0, 4.0]], ["A", "A", "A", "A"]))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: paired analysis requires two conditions" {
  let config = @src.Aldex2Config::create(mc_samples=4, paired=true) catch {
    _ => abort("paired ALDEx2 configuration should succeed")
  }
  let raised = try {
    ignore(
      @src.aldex2_clr(
        [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]],
        ["A", "A", "B", "B", "C", "C"],
        config~,
      ),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: paired groups require equal sample counts" {
  let config = @src.Aldex2Config::create(mc_samples=4, paired=true) catch {
    _ => abort("paired ALDEx2 configuration should succeed")
  }
  let raised = try {
    ignore(
      @src.aldex2_clr(
        [[1.0, 2.0, 3.0, 4.0, 5.0]],
        ["A", "A", "B", "B", "B"],
        config~,
      ),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: scale sample row count is validated" {
  let config = @src.Aldex2Config::create(mc_samples=4) catch {
    _ => abort("scale ALDEx2 configuration should succeed")
  }
  let raised = try {
    ignore(
      @src.aldex2_clr([[1.0, 2.0, 3.0, 4.0]], ["A", "A", "B", "B"], config~, scale_samples=[
        [1.0, 1.0, 1.0, 1.0],
      ]),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: scale sample column count is validated" {
  let config = @src.Aldex2Config::create(mc_samples=4) catch {
    _ => abort("scale ALDEx2 configuration should succeed")
  }
  let raised = try {
    ignore(
      @src.aldex2_clr([[1.0, 2.0, 3.0, 4.0]], ["A", "A", "B", "B"], config~, scale_samples=[
        [1.0],
        [1.0],
        [1.0],
        [1.0],
      ]),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: scale samples must be finite and positive" {
  let config = @src.Aldex2Config::create(mc_samples=2) catch {
    _ => abort("scale ALDEx2 configuration should succeed")
  }
  let raised = try {
    ignore(
      @src.aldex2_clr([[1.0, 2.0, 3.0, 4.0]], ["A", "A", "B", "B"], config~, scale_samples=[
        [1.0, 1.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [1.0, 1.0],
      ]),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: user denominator validates bounds" {
  let config = @src.Aldex2Config::create(
    mc_samples=4,
    denominator=@src.Aldex2User,
    denominator_indices=[2],
  ) catch {
    _ => abort("user ALDEx2 configuration should succeed")
  }
  let raised = try {
    ignore(
      @src.aldex2_clr(
        [[1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 4.0, 5.0]],
        ["A", "A", "B", "B"],
        config~,
      ),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
  let zero_config = @src.Aldex2Config::create(
    mc_samples=4,
    denominator=@src.Aldex2User,
    denominator_indices=[0],
  ) catch {
    _ => abort("all-zero denominator configuration should succeed")
  }
  let zero_raised = try {
    ignore(
      @src.aldex2_clr(
        [[0.0, 0.0, 0.0, 0.0], [2.0, 3.0, 4.0, 5.0]],
        ["A", "A", "B", "B"],
        config=zero_config,
      ),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(zero_raised)
}

///|
test "aldex2: user denominator rejects duplicates" {
  let config = @src.Aldex2Config::create(
    mc_samples=4,
    denominator=@src.Aldex2User,
    denominator_indices=[0, 0],
  ) catch {
    _ => abort("user ALDEx2 configuration should succeed")
  }
  let raised = try {
    ignore(
      @src.aldex2_clr(
        [[1.0, 2.0, 3.0, 4.0], [2.0, 3.0, 4.0, 5.0]],
        ["A", "A", "B", "B"],
        config~,
      ),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: IQLR reports insufficient invariant features" {
  let config = @src.Aldex2Config::create(
    mc_samples=4,
    denominator=@src.Aldex2Iqlr,
  ) catch {
    _ => abort("IQLR ALDEx2 configuration should succeed")
  }
  let raised = try {
    ignore(
      @src.aldex2_clr(
        [[10.0, 11.0, 20.0, 21.0], [5.0, 6.0, 7.0, 8.0]],
        ["A", "A", "B", "B"],
        config~,
      ),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: LVHA reports insufficient shared features" {
  let config = @src.Aldex2Config::create(
    mc_samples=4,
    denominator=@src.Aldex2Lvha,
  ) catch {
    _ => abort("LVHA ALDEx2 configuration should succeed")
  }
  let raised = try {
    ignore(
      @src.aldex2_clr([[10.0, 11.0, 20.0, 21.0]], ["A", "A", "B", "B"], config~),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: posterior cube dimensions follow samples features and instances" {
  let clr = aldex2_test_clr()
  assert_eq(clr.dirichlet.length(), 8)
  assert_eq(clr.dirichlet[0].length(), 11)
  assert_eq(clr.dirichlet[0][0].length(), 32)
  assert_eq(clr.analysis.length(), 8)
  assert_eq(clr.analysis[0].length(), 11)
  assert_eq(clr.analysis[0][0].length(), 32)
}

///|
test "aldex2: every Dirichlet instance sums to one" {
  let clr = aldex2_test_clr()
  for sample in 0..<clr.sample_count() {
    for instance in 0..<clr.mc_sample_count() {
      let mut total = 0.0
      for feature in 0..<clr.feature_count() {
        total = total + clr.dirichlet[sample][feature][instance]
      }
      aldex2_test_close(total, 1.0, 1.0e-12)
    }
  }
}

///|
test "aldex2: Dirichlet draws are finite and positive" {
  let clr = aldex2_test_clr()
  for sample in 0..<clr.sample_count() {
    for feature in 0..<clr.feature_count() {
      for value in clr.dirichlet[sample][feature] {
        assert_true(value > 0.0)
        assert_true(value < 1.0)
        assert_true(value == value)
      }
    }
  }
}

///|
test "aldex2: all denominator centers each CLR instance at zero" {
  let config = @src.Aldex2Config::create(mc_samples=8) catch {
    _ => abort("all-denominator ALDEx2 configuration should succeed")
  }
  let clr = @src.aldex2_clr(
    [[10.0, 12.0, 20.0, 22.0], [5.0, 6.0, 8.0, 9.0], [3.0, 4.0, 5.0, 6.0]],
    ["A", "A", "B", "B"],
    config~,
  ) catch {
    _ => abort("all-denominator ALDEx2 CLR should succeed")
  }
  for sample in 0..<clr.sample_count() {
    for instance in 0..<clr.mc_sample_count() {
      let values = aldex2_test_unwrap(clr.sample_instance(sample, instance))
      let mut total = 0.0
      for value in values {
        total = total + value
      }
      aldex2_test_close(total, 0.0, 1.0e-11)
    }
  }
}

///|
test "aldex2: median denominator centers each instance at zero median" {
  let config = @src.Aldex2Config::create(
    mc_samples=8,
    denominator=@src.Aldex2Median,
  ) catch {
    _ => abort("median ALDEx2 configuration should succeed")
  }
  let clr = @src.aldex2_clr(
    [
      [10.0, 12.0, 20.0, 22.0],
      [5.0, 6.0, 8.0, 9.0],
      [3.0, 4.0, 5.0, 6.0],
      [9.0, 8.0, 7.0, 6.0],
    ],
    ["A", "A", "B", "B"],
    config~,
  ) catch {
    _ => abort("median ALDEx2 CLR should succeed")
  }
  for sample in 0..<clr.sample_count() {
    for instance in 0..<clr.mc_sample_count() {
      let values = aldex2_test_unwrap(clr.sample_instance(sample, instance))
      aldex2_test_close(aldex2_test_median(values), 0.0, 1.0e-12)
    }
  }
}

///|
test "aldex2: equal seeds produce identical posterior draws" {
  let first = aldex2_test_clr()
  let second = aldex2_test_clr()
  assert_eq(first.dirichlet, second.dirichlet)
  assert_eq(first.analysis, second.analysis)
}

///|
test "aldex2: different seeds produce different posterior draws" {
  let (counts, conditions, features, samples) = @src.aldex2_example_data()
  let first = @src.aldex2_clr(
    counts,
    conditions,
    config=@src.Aldex2Config::create(mc_samples=4, seed=1) catch {
      _ => abort("first ALDEx2 configuration should succeed")
    },
    feature_names=features,
    sample_names=samples,
  ) catch {
    _ => abort("first ALDEx2 CLR should succeed")
  }
  let second = @src.aldex2_clr(
    counts,
    conditions,
    config=@src.Aldex2Config::create(mc_samples=4, seed=2) catch {
      _ => abort("second ALDEx2 configuration should succeed")
    },
    feature_names=features,
    sample_names=samples,
  ) catch {
    _ => abort("second ALDEx2 CLR should succeed")
  }
  assert_true(first.dirichlet != second.dirichlet)
}

///|
test "aldex2: low-shape gamma path remains finite" {
  let config = @src.Aldex2Config::create(mc_samples=16, prior=0.1) catch {
    _ => abort("low-prior ALDEx2 configuration should succeed")
  }
  let clr = @src.aldex2_clr(
    [[0.0, 0.0, 1.0, 1.0], [1.0, 1.0, 0.0, 0.0]],
    ["A", "A", "B", "B"],
    config~,
  ) catch {
    _ => abort("low-shape ALDEx2 CLR should succeed")
  }
  for feature in 0..<clr.feature_count() {
    for value in clr.dirichlet[0][feature] {
      assert_true(value > 0.0)
      assert_true(value == value)
    }
  }
}

///|
test "aldex2: posterior slice queries validate bounds" {
  let clr = aldex2_test_clr()
  assert_true(clr.sample_instance(0, 0) is Some(_))
  assert_true(clr.feature_instance(0, 0) is Some(_))
  assert_true(clr.sample_instance(-1, 0) is None)
  assert_true(clr.sample_instance(0, 32) is None)
  assert_true(clr.feature_instance(11, 0) is None)
  assert_true(clr.feature_instance(0, -1) is None)
}

///|
test "aldex2: posterior slice queries use documented orientation" {
  let clr = aldex2_test_clr()
  let sample = aldex2_test_unwrap(clr.sample_instance(0, 0))
  let feature = aldex2_test_unwrap(clr.feature_instance(0, 0))
  assert_eq(sample.length(), 11)
  assert_eq(feature.length(), 8)
  aldex2_test_close(sample[0], feature[0], 0.0)
}

///|
test "aldex2: expected CLR query validates indices" {
  let clr = aldex2_test_clr()
  assert_true(clr.expected_clr(0, 0) is Some(_))
  assert_true(clr.expected_clr(-1, 0) is None)
  assert_true(clr.expected_clr(0, 8) is None)
}

///|
test "aldex2: expected proportions are bounded" {
  let clr = aldex2_test_clr()
  for feature in 0..<clr.feature_count() {
    let value = match clr.expected_proportion(feature, 0) {
      Some(value) => value
      None => abort("expected ALDEx2 proportion")
    }
    assert_true(value > 0.0)
    assert_true(value < 1.0)
  }
  assert_true(clr.expected_proportion(11, 0) is None)
}

///|
test "aldex2: IQLR finds shared invariant features in example data" {
  let clr = aldex2_test_clr()
  assert_true(clr.denominator_by_sample[0].length() >= 2)
  for sample in 1..<clr.sample_count() {
    assert_eq(clr.denominator_by_sample[sample], clr.denominator_by_sample[0])
  }
}

///|
test "aldex2: LVHA selects shared low-variance high-abundance features" {
  let config = @src.Aldex2Config::create(
    mc_samples=4,
    denominator=@src.Aldex2Lvha,
  ) catch {
    _ => abort("LVHA ALDEx2 configuration should succeed")
  }
  let clr = @src.aldex2_clr(
    [
      [100.0, 100.0, 100.0, 100.0],
      [90.0, 90.0, 90.0, 90.0],
      [80.0, 80.0, 80.0, 80.0],
      [70.0, 70.0, 70.0, 70.0],
      [60.0, 60.0, 60.0, 60.0],
      [50.0, 50.0, 50.0, 50.0],
      [40.0, 40.0, 40.0, 40.0],
      [30.0, 30.0, 30.0, 30.0],
    ],
    ["A", "A", "B", "B"],
    config~,
  ) catch {
    Aldex2Error(message) => abort("LVHA ALDEx2 CLR should succeed: " + message)
  }
  assert_true(clr.denominator_by_sample[0].length() >= 2)
  assert_eq(clr.denominator_by_sample[0], clr.denominator_by_sample[3])
}

///|
test "aldex2: zero denominator is condition specific" {
  let config = @src.Aldex2Config::create(
    mc_samples=4,
    denominator=@src.Aldex2Zero,
  ) catch {
    _ => abort("zero ALDEx2 configuration should succeed")
  }
  let clr = @src.aldex2_clr(
    [[10.0, 11.0, 0.0, 0.0], [0.0, 0.0, 12.0, 13.0], [5.0, 6.0, 7.0, 8.0]],
    ["A", "A", "B", "B"],
    config~,
  ) catch {
    _ => abort("zero-denominator ALDEx2 CLR should succeed")
  }
  assert_eq(clr.denominator_by_sample[0], [0, 2])
  assert_eq(clr.denominator_by_sample[1], [0, 2])
  assert_eq(clr.denominator_by_sample[2], [1, 2])
  assert_eq(clr.denominator_by_sample[3], [1, 2])
}

///|
test "aldex2: user denominator centers selected features" {
  let config = @src.Aldex2Config::create(
    mc_samples=4,
    denominator=@src.Aldex2User,
    denominator_indices=[1, 3],
  ) catch {
    _ => abort("user ALDEx2 configuration should succeed")
  }
  let clr = @src.aldex2_clr(
    [
      [0.0, 0.0, 0.0, 0.0],
      [10.0, 11.0, 20.0, 21.0],
      [5.0, 6.0, 8.0, 9.0],
      [3.0, 4.0, 5.0, 6.0],
    ],
    ["A", "A", "B", "B"],
    config~,
  ) catch {
    _ => abort("user-denominator ALDEx2 CLR should succeed")
  }
  assert_eq(clr.denominator_by_sample[0], [0, 2])
  assert_eq(clr.removed_feature_count(), 1)
  let values = aldex2_test_unwrap(clr.sample_instance(0, 0))
  aldex2_test_close((values[0] + values[2]) / 2.0, 0.0, 1.0e-12)
}

///|
test "aldex2: scale samples shift log abundance by log2 scale" {
  let config = @src.Aldex2Config::create(mc_samples=4, seed=7) catch {
    _ => abort("scale ALDEx2 configuration should succeed")
  }
  let counts = [
    [10.0, 11.0, 20.0, 21.0],
    [5.0, 6.0, 8.0, 9.0],
    [3.0, 4.0, 5.0, 6.0],
  ]
  let conditions = ["A", "A", "B", "B"]
  let unit = @src.aldex2_clr(
    counts,
    conditions,
    config~,
    scale_samples=aldex2_test_scale(4, 4, 1.0),
  ) catch {
    _ => abort("unit-scale ALDEx2 CLR should succeed")
  }
  let four = @src.aldex2_clr(
    counts,
    conditions,
    config~,
    scale_samples=aldex2_test_scale(4, 4, 4.0),
  ) catch {
    _ => abort("four-scale ALDEx2 CLR should succeed")
  }
  let unit_values = aldex2_test_unwrap(unit.sample_instance(0, 0))
  let four_values = aldex2_test_unwrap(four.sample_instance(0, 0))
  for feature in 0..<unit_values.length() {
    aldex2_test_close(four_values[feature] - unit_values[feature], 2.0, 1.0e-12)
  }
}

///|
test "aldex2: scale samples preserve within-sample log ratios" {
  let config = @src.Aldex2Config::create(mc_samples=4, seed=7) catch {
    _ => abort("scale ALDEx2 configuration should succeed")
  }
  let counts = [
    [10.0, 11.0, 20.0, 21.0],
    [5.0, 6.0, 8.0, 9.0],
    [3.0, 4.0, 5.0, 6.0],
  ]
  let conditions = ["A", "A", "B", "B"]
  let first = @src.aldex2_clr(
    counts,
    conditions,
    config~,
    scale_samples=aldex2_test_scale(4, 4, 1.0),
  ) catch {
    _ => abort("first scale ALDEx2 CLR should succeed")
  }
  let second = @src.aldex2_clr(
    counts,
    conditions,
    config~,
    scale_samples=aldex2_test_scale(4, 4, 10.0),
  ) catch {
    _ => abort("second scale ALDEx2 CLR should succeed")
  }
  let first_values = aldex2_test_unwrap(first.sample_instance(0, 0))
  let second_values = aldex2_test_unwrap(second.sample_instance(0, 0))
  aldex2_test_close(
    first_values[0] - first_values[1],
    second_values[0] - second_values[1],
    1.0e-12,
  )
}

///|
test "aldex2: CLR construction permits multi-condition designs" {
  let clr = @src.aldex2_clr(
    [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [6.0, 5.0, 4.0, 3.0, 2.0, 1.0]],
    ["A", "A", "B", "B", "C", "C"],
    config=@src.Aldex2Config::create(mc_samples=4) catch {
      _ => abort("multi-condition ALDEx2 configuration should succeed")
    },
  ) catch {
    _ => abort("multi-condition ALDEx2 CLR should succeed")
  }
  assert_eq(clr.group_names, ["A", "B", "C"])
  assert_eq(clr.group_indices, [[0, 1], [2, 3], [4, 5]])
}

///|
test "aldex2: t-test requires exactly two conditions" {
  let clr = @src.aldex2_clr(
    [[1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [6.0, 5.0, 4.0, 3.0, 2.0, 1.0]],
    ["A", "A", "B", "B", "C", "C"],
    config=@src.Aldex2Config::create(mc_samples=4) catch {
      _ => abort("multi-condition ALDEx2 configuration should succeed")
    },
  ) catch {
    _ => abort("multi-condition ALDEx2 CLR should succeed")
  }
  let raised = try {
    ignore(@src.aldex2_ttest(clr))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: tests require two replicates per condition" {
  let clr = @src.aldex2_clr(
    [[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]],
    ["A", "B", "B"],
    config=@src.Aldex2Config::create(mc_samples=4) catch {
      _ => abort("small ALDEx2 configuration should succeed")
    },
  ) catch {
    _ => abort("small ALDEx2 CLR should succeed")
  }
  let raised = try {
    ignore(@src.aldex2_ttest(clr))
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}

///|
test "aldex2: posterior tests preserve feature metadata" {
  let result = aldex2_test_result()
  assert_eq(result.tests.features.length(), 11)
  assert_eq(result.tests.features[0].feature_name, "increased")
  assert_eq(result.tests.group_a, "control")
  assert_eq(result.tests.group_b, "treated")
  assert_false(result.tests.paired)
  assert_eq(result.tests.mc_samples, 32)
}

///|
test "aldex2: posterior expected p-values are bounded" {
  let result = aldex2_test_result()
  for feature in result.tests.features {
    assert_true(feature.welch_p_value >= 0.0)
    assert_true(feature.welch_p_value <= 1.0)
    assert_true(feature.wilcoxon_p_value >= 0.0)
    assert_true(feature.wilcoxon_p_value <= 1.0)
  }
}

///|
test "aldex2: posterior BH values are bounded and conservative" {
  let result = aldex2_test_result()
  for feature in result.tests.features {
    assert_true(feature.welch_adjusted >= feature.welch_p_value - 1.0e-12)
    assert_true(feature.welch_adjusted <= 1.0)
    assert_true(feature.wilcoxon_adjusted >= feature.wilcoxon_p_value - 1.0e-12)
    assert_true(feature.wilcoxon_adjusted <= 1.0)
  }
}

///|
test "aldex2: strong example signals have expected effect directions" {
  let result = aldex2_test_result()
  assert_true(aldex2_test_effect(result, "increased").effect > 3.0)
  assert_true(aldex2_test_effect(result, "decreased").effect < -3.0)
  assert_true(aldex2_test_effect(result, "rare_increased").effect > 1.0)
}

///|
test "aldex2: strong example signals have small Welch eBH" {
  let result = aldex2_test_result()
  assert_true(aldex2_test_statistics(result, "increased").welch_adjusted < 0.01)
  assert_true(aldex2_test_statistics(result, "decreased").welch_adjusted < 0.01)
  assert_true(
    aldex2_test_statistics(result, "rare_increased").welch_adjusted < 0.1,
  )
}

///|
test "aldex2: stable feature remains near zero" {
  let result = aldex2_test_result()
  let stable = aldex2_test_effect(result, "stable_high")
  assert_true(stable.effect.abs() < 1.0)
  assert_true(stable.overlap > 0.2)
  assert_true(
    aldex2_test_statistics(result, "stable_high").welch_adjusted > 0.2,
  )
}

///|
test "aldex2: effect intervals contain posterior medians" {
  let result = aldex2_test_result()
  for feature in result.effects.features {
    assert_true(feature.effect_low <= feature.effect)
    assert_true(feature.effect_high >= feature.effect)
  }
}

///|
test "aldex2: effect overlap is a bounded sign probability" {
  let result = aldex2_test_result()
  for feature in result.effects.features {
    assert_true(feature.overlap > 0.0)
    assert_true(feature.overlap <= 0.5)
  }
  assert_true(
    aldex2_test_effect(result, "increased").overlap <
    aldex2_test_effect(result, "stable_high").overlap,
  )
}

///|
test "aldex2: paired analysis uses paired test and effect paths" {
  let (counts, conditions, features, samples) = @src.aldex2_example_data()
  let config = @src.Aldex2Config::create(
    mc_samples=32,
    denominator=@src.Aldex2Iqlr,
    seed=2026,
    paired=true,
  ) catch {
    _ => abort("paired ALDEx2 configuration should succeed")
  }
  let result = @src.aldex2(
    counts,
    conditions,
    config~,
    feature_names=features,
    sample_names=samples,
  ) catch {
    _ => abort("paired ALDEx2 analysis should succeed")
  }
  assert_true(result.tests.paired)
  assert_true(result.effects.paired)
  assert_true(aldex2_test_effect(result, "increased").effect > 1.0)
  assert_true(aldex2_test_effect(result, "decreased").effect < -1.0)
}

///|
test "aldex2: result selection combines BH and effect thresholds" {
  let result = aldex2_test_result()
  let selected = result.select(maximum_adjusted=0.1, minimum_effect=1.0)
  assert_true(selected.length() >= 3)
  for feature in selected {
    assert_true(feature.effect.abs() >= 1.0)
    assert_true(
      aldex2_test_statistics(result, feature.feature_name).welch_adjusted <= 0.1,
    )
  }
}

///|
test "aldex2: result ranking orders absolute effects descending" {
  let ranked = aldex2_test_result().ranked()
  assert_eq(ranked.length(), 11)
  for index in 1..<ranked.length() {
    assert_true(ranked[index - 1].effect.abs() >= ranked[index].effect.abs())
  }
  assert_true(
    ranked[0].feature_name == "increased" ||
    ranked[0].feature_name == "decreased",
  )
}

///|
test "aldex2: result summary reports analysis controls" {
  let summary = aldex2_test_result().summary(
    maximum_adjusted=0.1,
    minimum_effect=1.0,
  )
  assert_eq(summary.feature_count, 11)
  assert_eq(summary.sample_count, 8)
  assert_eq(summary.removed_feature_count, 1)
  assert_eq(summary.group_a, "control")
  assert_eq(summary.group_b, "treated")
  assert_eq(summary.mc_samples, 32)
  assert_eq(summary.denominator, "iqlr")
  assert_false(summary.scale_model)
  assert_true(summary.significant_welch_count >= 3)
  assert_true(summary.large_effect_count >= 3)
}

///|
test "aldex2: scaled analysis is recorded in summary" {
  let (counts, conditions, features, samples) = @src.aldex2_example_data()
  let config = @src.Aldex2Config::create(mc_samples=8) catch {
    _ => abort("scaled ALDEx2 configuration should succeed")
  }
  let result = @src.aldex2(
    counts,
    conditions,
    config~,
    feature_names=features,
    sample_names=samples,
    scale_samples=aldex2_test_scale(8, 8, 1.0),
  ) catch {
    _ => abort("scaled ALDEx2 analysis should succeed")
  }
  assert_true(result.summary().scale_model)
}

///|
test "aldex2: expected distance matrix matches sample count" {
  let distances = aldex2_test_clr().expected_distance()
  assert_eq(distances.length(), 8)
  for row in distances {
    assert_eq(row.length(), 8)
  }
}

///|
test "aldex2: expected distances are symmetric with zero diagonal" {
  let distances = aldex2_test_clr().expected_distance()
  for left in 0..<distances.length() {
    assert_eq(distances[left][left], 0.0)
    for right in 0..<distances.length() {
      assert_eq(distances[left][right], distances[right][left])
      assert_true(distances[left][right] >= 0.0)
    }
  }
  assert_true(distances[0][4] > 0.0)
}

///|
test "aldex2: TSV export contains upstream-style columns and rows" {
  let table = aldex2_test_result().to_tsv()
  assert_true(table.has_prefix("feature\trab.all\trab.group.a"))
  assert_true(table.contains("\twe.ep\twe.eBH\twi.ep\twi.eBH\n"))
  assert_true(table.contains("increased\t"))
  assert_true(table.contains("decreased\t"))
  assert_false(table.contains("all_zero\t"))
}

///|
test "aldex2: feature lookup returns none for unknown names" {
  let result = aldex2_test_result()
  assert_true(result.tests.find_feature("missing") is None)
  assert_true(result.effects.find_feature("missing") is None)
}

///|
test "aldex2: SummarizedExperiment output adds posterior assays" {
  let experiment = aldex2_test_experiment()
  let (_, conditions, features, samples) = @src.aldex2_example_data()
  let config = aldex2_test_config()
  let output = @src.aldex2_summarized_experiment(
    experiment,
    conditions,
    config~,
    feature_names=features,
    sample_names=samples,
  ) catch {
    _ => abort("ALDEx2 SummarizedExperiment integration should succeed")
  }
  for
    name in [
      "aldex2_effect", "aldex2_overlap", "aldex2_welch_bh", "aldex2_wilcoxon_bh",
    ] {
    let assay = match @src.se_assay(output.experiment, name) {
      Some(value) => value
      None => abort("expected ALDEx2 assay: " + name)
    }
    assert_eq(assay.length(), 12)
    assert_eq(assay[0].length(), 1)
  }
}

///|
test "aldex2: SummarizedExperiment reinserts filtered zero features" {
  let experiment = aldex2_test_experiment()
  let (_, conditions, features, samples) = @src.aldex2_example_data()
  let config = aldex2_test_config()
  let output = @src.aldex2_summarized_experiment(
    experiment,
    conditions,
    config~,
    feature_names=features,
    sample_names=samples,
  ) catch {
    _ => abort("ALDEx2 SummarizedExperiment integration should succeed")
  }
  let effect = match @src.se_assay(output.experiment, "aldex2_effect") {
    Some(value) => value
    None => abort("ALDEx2 effect assay missing")
  }
  let welch = match @src.se_assay(output.experiment, "aldex2_welch_bh") {
    Some(value) => value
    None => abort("ALDEx2 Welch assay missing")
  }
  assert_eq(effect[11], [0.0])
  assert_eq(welch[11], [1.0])
}

///|
test "aldex2: SummarizedExperiment integration is immutable" {
  let experiment = aldex2_test_experiment()
  let (_, conditions, features, samples) = @src.aldex2_example_data()
  let config = aldex2_test_config()
  let output = @src.aldex2_summarized_experiment(
    experiment,
    conditions,
    config~,
    feature_names=features,
    sample_names=samples,
  ) catch {
    _ => abort("ALDEx2 SummarizedExperiment integration should succeed")
  }
  assert_true(@src.se_assay(experiment, "aldex2_effect") is None)
  assert_true(experiment.metadata.get("aldex2_mc_samples") is None)
  assert_eq(output.experiment.metadata["source"], "aldex2_test")
  assert_eq(output.experiment.metadata["aldex2_mc_samples"], "32")
  assert_eq(output.experiment.metadata["aldex2_denominator"], "iqlr")
  assert_eq(output.experiment.metadata["aldex2_groups"], "control,treated")
}

///|
test "aldex2: SummarizedExperiment rejects missing assays" {
  let experiment = aldex2_test_experiment()
  let (_, conditions, features, samples) = @src.aldex2_example_data()
  let raised = try {
    ignore(
      @src.aldex2_summarized_experiment(
        experiment,
        conditions,
        assay_name="missing",
        feature_names=features,
        sample_names=samples,
      ),
    )
    false
  } catch {
    Aldex2Error(_) => true
  }
  assert_true(raised)
}
