///|
fn scuttle_test_close(
  left : Double,
  right : Double,
  tolerance : Double,
) -> Bool {
  (left - right).abs() <= tolerance
}

///|
fn scuttle_test_column_sums(matrix : Array[Array[Double]]) -> Array[Double] {
  let cells = if matrix.length() > 0 { matrix[0].length() } else { 0 }
  let output = Array::make(cells, 0.0)
  for row in matrix {
    for cell in 0..<cells {
      output[cell] = output[cell] + row[cell]
    }
  }
  output
}

///|
fn scuttle_test_total(matrix : Array[Array[Double]]) -> Double {
  let mut output = 0.0
  for row in matrix {
    for value in row {
      output = output + value
    }
  }
  output
}

///|
fn scuttle_test_outlier_value(value : Bool?) -> Bool {
  match value {
    Some(result) => result
    None => abort("expected an estimated outlier status")
  }
}

///|
fn scuttle_test_counts() -> Array[Array[Double]] {
  [
    [8.0, 4.0, 16.0, 8.0],
    [2.0, 6.0, 4.0, 12.0],
    [0.0, 2.0, 0.0, 4.0],
    [10.0, 8.0, 20.0, 16.0],
  ]
}

///|
fn scuttle_test_sce() -> @src.SingleCellExperiment {
  let counts = scuttle_test_counts()
  let sce = @src.SingleCellExperiment::new(counts, ["G1", "G2", "G3", "G4"], [
    "C1", "C2", "C3", "C4",
  ])
  sce.assays["logcounts"] = [
    [3.0, 2.0, 4.0, 3.0],
    [1.0, 2.5, 2.0, 3.5],
    [0.0, 1.0, 0.0, 2.0],
    [3.5, 3.0, 4.5, 4.0],
  ]
  sce.row_data["symbol"] = ["A", "B", "C", "D"]
  sce.col_data["batch"] = ["A", "A", "B", "B"]
  sce.reduced_dims["PCA"] = [[1.0, 0.0], [2.0, 0.0], [3.0, 1.0], [4.0, 1.0]]
  sce.metadata["project"] = "scuttle-test"
  sce.alternative_experiments["spike"] = @src.SingleCellExperiment::new(
    [[1.0, 2.0, 3.0, 4.0]],
    ["Spike1"],
    ["C1", "C2", "C3", "C4"],
  )
  sce
}

///|
test "scuttle: default outlier configuration matches upstream defaults" {
  let config = @src.ScuttleOutlierConfig::default()
  assert_eq(config.nmads, 3.0)
  assert_true(config.direction == @src.scuttle_outlier_both())
  assert_false(config.log_transform)
  assert_false(config.share_medians)
  assert_false(config.share_mads)
  assert_true(config.share_missing)
  assert_true(config.subset is None)
  assert_true(config.min_diff is None)
}

///|
test "scuttle: outlier configuration rejects negative nmads" {
  let mut raised = false
  ignore(@src.ScuttleOutlierConfig::create(nmads=-1.0)) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "scuttle: outlier configuration rejects negative min diff" {
  let mut raised = false
  ignore(@src.ScuttleOutlierConfig::create(min_diff=Some(-1.0))) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "scuttle: vanilla outlier thresholds use R MAD scaling" {
  let result = @src.scuttle_is_outlier([1.0, 2.0, 2.0, 3.0, 100.0]) catch {
    _ => abort("valid outlier detection should work")
  }
  assert_eq(result.batch_names, ["1"])
  assert_true(scuttle_test_close(result.medians[0], 2.0, 1.0e-12))
  assert_true(scuttle_test_close(result.mads[0], 1.4826, 1.0e-12))
  assert_true(scuttle_test_close(result.lower[0], -2.4478, 1.0e-10))
  assert_true(scuttle_test_close(result.higher[0], 6.4478, 1.0e-10))
  assert_false(scuttle_test_outlier_value(result.outliers[0]))
  assert_true(scuttle_test_outlier_value(result.outliers[4]))
  assert_eq(result.discarded_count(), 1)
}

///|
test "scuttle: lower-tail direction only calls small observations" {
  let config = @src.ScuttleOutlierConfig::create(
    nmads=1.0,
    direction=@src.scuttle_outlier_lower(),
  ) catch {
    _ => abort("valid lower-tail configuration should build")
  }
  let result = @src.scuttle_is_outlier([0.0, 10.0, 10.0, 11.0, 100.0], config~) catch {
    _ => abort("lower-tail outlier detection should work")
  }
  assert_true(scuttle_test_outlier_value(result.outliers[0]))
  assert_false(scuttle_test_outlier_value(result.outliers[4]))
  assert_true(result.higher[0] >= 1.0e299)
}

///|
test "scuttle: higher-tail direction only calls large observations" {
  let config = @src.ScuttleOutlierConfig::create(
    nmads=1.0,
    direction=@src.scuttle_outlier_higher(),
  ) catch {
    _ => abort("valid higher-tail configuration should build")
  }
  let result = @src.scuttle_is_outlier([0.0, 10.0, 10.0, 11.0, 100.0], config~) catch {
    _ => abort("higher-tail outlier detection should work")
  }
  assert_false(scuttle_test_outlier_value(result.outliers[0]))
  assert_true(scuttle_test_outlier_value(result.outliers[4]))
  assert_true(result.lower[0] <= -1.0e299)
}

///|
test "scuttle: nmads changes the outlier thresholds" {
  let config = @src.ScuttleOutlierConfig::create(nmads=5.0) catch {
    _ => abort("valid nmads configuration should build")
  }
  let result = @src.scuttle_is_outlier([1.0, 2.0, 2.0, 3.0, 100.0], config~) catch {
    _ => abort("outlier detection should work")
  }
  assert_true(scuttle_test_close(result.lower[0], -5.413, 1.0e-10))
  assert_true(scuttle_test_close(result.higher[0], 9.413, 1.0e-10))
}

///|
test "scuttle: min diff dominates a small MAD" {
  let config = @src.ScuttleOutlierConfig::create(nmads=1.0, min_diff=Some(10.0)) catch {
    _ => abort("valid min-diff configuration should build")
  }
  let result = @src.scuttle_is_outlier([1.0, 2.0, 2.0, 3.0, 20.0], config~) catch {
    _ => abort("outlier detection should work")
  }
  assert_eq(result.lower[0], -8.0)
  assert_eq(result.higher[0], 12.0)
  assert_true(scuttle_test_outlier_value(result.outliers[4]))
}

///|
test "scuttle: threshold subset is applied to all observations" {
  let config = @src.ScuttleOutlierConfig::create(
    nmads=1.0,
    subset=Some([0, 1, 2]),
  ) catch {
    _ => abort("valid subset configuration should build")
  }
  let result = @src.scuttle_is_outlier([1.0, 2.0, 3.0, 50.0], config~) catch {
    _ => abort("subset outlier detection should work")
  }
  assert_eq(result.medians[0], 2.0)
  assert_true(scuttle_test_outlier_value(result.outliers[3]))
}

///|
test "scuttle: explicit empty threshold subset yields missing statuses" {
  let config = @src.ScuttleOutlierConfig::create(subset=Some([])) catch {
    _ => abort("empty subset is a valid explicit subset")
  }
  let result = @src.scuttle_is_outlier([1.0, 2.0, 3.0], config~) catch {
    _ => abort("empty-subset detection should return missing statuses")
  }
  assert_false(result.threshold_valid[0])
  assert_true(result.outliers[0] is None)
  assert_true(result.outliers[2] is None)
}

///|
test "scuttle: threshold subset indices are validated" {
  let config = @src.ScuttleOutlierConfig::create(subset=Some([3])) catch {
    _ => abort("subset bounds are checked when data are available")
  }
  let mut raised = false
  ignore(@src.scuttle_is_outlier([1.0, 2.0, 3.0], config~)) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "scuttle: batch-specific thresholds are independent" {
  let config = @src.ScuttleOutlierConfig::create(nmads=1.0, batches=[
    "A", "A", "A", "B", "B", "B",
  ]) catch {
    _ => abort("valid batch configuration should build")
  }
  let result = @src.scuttle_is_outlier(
    [1.0, 2.0, 3.0, 100.0, 101.0, 102.0],
    config~,
  ) catch {
    _ => abort("batch-aware detection should work")
  }
  assert_eq(result.batch_names, ["A", "B"])
  assert_eq(result.medians, [2.0, 101.0])
  assert_false(scuttle_test_outlier_value(result.outliers[0]))
  assert_false(scuttle_test_outlier_value(result.outliers[5]))
}

///|
test "scuttle: batch names use sorted factor order" {
  let config = @src.ScuttleOutlierConfig::create(batches=["z", "a", "m", "z"]) catch {
    _ => abort("valid batches should build")
  }
  let result = @src.scuttle_is_outlier([1.0, 2.0, 3.0, 4.0], config~) catch {
    _ => abort("batch-aware detection should work")
  }
  assert_eq(result.batch_names, ["a", "m", "z"])
}

///|
test "scuttle: shared medians and MADs reproduce global thresholds" {
  let config = @src.ScuttleOutlierConfig::create(
    batches=["A", "A", "B", "B", "B"],
    share_medians=true,
    share_mads=true,
  ) catch {
    _ => abort("valid sharing configuration should build")
  }
  let result = @src.scuttle_is_outlier([1.0, 2.0, 3.0, 4.0, 100.0], config~) catch {
    _ => abort("shared outlier detection should work")
  }
  let global = @src.scuttle_is_outlier([1.0, 2.0, 3.0, 4.0, 100.0]) catch {
    _ => abort("global outlier detection should work")
  }
  assert_eq(result.medians, [global.medians[0], global.medians[0]])
  assert_eq(result.mads, [global.mads[0], global.mads[0]])
  assert_eq(result.lower, [global.lower[0], global.lower[0]])
}

///|
test "scuttle: shared median retains batch-specific MADs" {
  let config = @src.ScuttleOutlierConfig::create(
    nmads=1.0,
    batches=["A", "A", "A", "B", "B", "B"],
    share_medians=true,
  ) catch {
    _ => abort("valid median-sharing configuration should build")
  }
  let result = @src.scuttle_is_outlier(
    [1.0, 2.0, 3.0, 10.0, 20.0, 30.0],
    config~,
  ) catch {
    _ => abort("shared-median detection should work")
  }
  assert_eq(result.medians[0], result.medians[1])
  assert_true(result.mads[1] > result.mads[0])
}

///|
test "scuttle: shared MAD retains batch-specific medians" {
  let config = @src.ScuttleOutlierConfig::create(
    batches=["A", "A", "A", "B", "B", "B"],
    share_mads=true,
  ) catch {
    _ => abort("valid MAD-sharing configuration should build")
  }
  let result = @src.scuttle_is_outlier(
    [1.0, 2.0, 3.0, 100.0, 101.0, 102.0],
    config~,
  ) catch {
    _ => abort("shared-MAD detection should work")
  }
  assert_eq(result.medians, [2.0, 101.0])
  assert_eq(result.mads[0], result.mads[1])
}

///|
test "scuttle: missing batch thresholds are recovered by default" {
  let config = @src.ScuttleOutlierConfig::create(
    batches=["A", "A", "B", "B"],
    subset=Some([0, 1]),
  ) catch {
    _ => abort("valid missing-batch configuration should build")
  }
  let result = @src.scuttle_is_outlier([1.0, 3.0, 100.0, 200.0], config~) catch {
    _ => abort("missing batch should borrow thresholds")
  }
  assert_true(result.threshold_valid[0])
  assert_true(result.threshold_valid[1])
  assert_eq(result.lower[0], result.lower[1])
  assert_eq(result.higher[0], result.higher[1])
}

///|
test "scuttle: missing batch statuses can remain missing" {
  let config = @src.ScuttleOutlierConfig::create(
    batches=["A", "A", "B", "B"],
    subset=Some([0, 1]),
    share_missing=false,
  ) catch {
    _ => abort("valid no-sharing configuration should build")
  }
  let result = @src.scuttle_is_outlier([1.0, 3.0, 100.0, 200.0], config~) catch {
    _ => abort("missing batch should remain missing")
  }
  assert_true(result.threshold_valid[0])
  assert_false(result.threshold_valid[1])
  assert_true(result.outliers[2] is None)
}

///|
test "scuttle: batch vector length is validated" {
  let config = @src.ScuttleOutlierConfig::create(batches=["A"]) catch {
    _ => abort("configuration construction does not know metric length")
  }
  let mut raised = false
  ignore(@src.scuttle_is_outlier([1.0, 2.0], config~)) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "scuttle: log thresholds are returned on the original scale" {
  let config = @src.ScuttleOutlierConfig::create(nmads=1.0, log_transform=true) catch {
    _ => abort("valid log configuration should build")
  }
  let result = @src.scuttle_is_outlier([1.0, 2.0, 4.0, 8.0, 1024.0], config~) catch {
    _ => abort("log outlier detection should work")
  }
  assert_true(scuttle_test_close(result.medians[0], 4.0, 1.0e-10))
  assert_true(result.lower[0] > 1.0)
  assert_true(result.higher[0] < 16.0)
  assert_true(scuttle_test_outlier_value(result.outliers[4]))
}

///|
test "scuttle: zero can participate in log outlier detection" {
  let config = @src.ScuttleOutlierConfig::create(log_transform=true) catch {
    _ => abort("valid log configuration should build")
  }
  let result = @src.scuttle_is_outlier([0.0, 8.0, 8.0, 8.0, 8.0], config~) catch {
    _ => abort("zero should be represented on the log scale")
  }
  assert_true(scuttle_test_outlier_value(result.outliers[0]))
  assert_true(scuttle_test_close(result.medians[0], 8.0, 1.0e-10))
}

///|
test "scuttle: negative log metrics are rejected" {
  let config = @src.ScuttleOutlierConfig::create(log_transform=true) catch {
    _ => abort("valid log configuration should build")
  }
  let mut raised = false
  ignore(@src.scuttle_is_outlier([-1.0, 2.0], config~)) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "scuttle: empty metric returns an empty filter" {
  let result = @src.scuttle_is_outlier([]) catch {
    _ => abort("empty metrics should be supported")
  }
  assert_eq(result.outliers.length(), 0)
  assert_false(result.threshold_valid[0])
}

///|
test "scuttle: per-feature QC computes means and detection percentages" {
  let metrics = @src.scuttle_per_feature_qc([
    [0.0, 2.0, 4.0, 6.0],
    [1.0, 0.0, 0.0, 3.0],
  ]) catch {
    _ => abort("valid feature QC should work")
  }
  assert_eq(metrics.means, [3.0, 1.0])
  assert_eq(metrics.detected_percent, [75.0, 50.0])
}

///|
test "scuttle: feature detection limit is strict" {
  let metrics = @src.scuttle_per_feature_qc(
    [[0.0, 2.0, 4.0, 6.0]],
    detection_limit=2.0,
  ) catch {
    _ => abort("feature QC with a limit should work")
  }
  assert_eq(metrics.detected_percent, [50.0])
}

///|
test "scuttle: per-feature QC computes subset statistics" {
  let subset = @src.ScuttleNamedCellSubset::create("controls", [0, 1]) catch {
    _ => abort("valid subset should build")
  }
  let metrics = @src.scuttle_per_feature_qc(
    [[0.0, 2.0, 4.0, 6.0], [1.0, 0.0, 0.0, 3.0]],
    subsets=[subset],
  ) catch {
    _ => abort("subset feature QC should work")
  }
  assert_eq(metrics.subset_names, ["controls"])
  assert_eq(metrics.subset_means[0], [1.0, 0.5])
  assert_eq(metrics.subset_detected_percent[0], [50.0, 50.0])
  assert_true(
    scuttle_test_close(metrics.subset_ratios[0][0], 1.0 / 3.0, 1.0e-12),
  )
}

///|
test "scuttle: per-feature QC supports multiple cell subsets" {
  let first = @src.ScuttleNamedCellSubset::create("first", [0, 1]) catch {
    _ => abort("valid subset should build")
  }
  let last = @src.ScuttleNamedCellSubset::create("last", [2, 3]) catch {
    _ => abort("valid subset should build")
  }
  let metrics = @src.scuttle_per_feature_qc([[0.0, 2.0, 4.0, 6.0]], subsets=[
    first, last,
  ]) catch {
    _ => abort("multiple-subset feature QC should work")
  }
  assert_eq(metrics.subset_names, ["first", "last"])
  assert_eq(metrics.subset_means[0], [1.0])
  assert_eq(metrics.subset_means[1], [5.0])
}

///|
test "scuttle: zero feature mean has a finite zero subset ratio" {
  let subset = @src.ScuttleNamedCellSubset::create("empty", []) catch {
    _ => abort("empty named subset should build")
  }
  let metrics = @src.scuttle_per_feature_qc([[0.0, 0.0]], subsets=[subset]) catch {
    _ => abort("zero feature QC should work")
  }
  assert_eq(metrics.subset_ratios[0], [0.0])
}

///|
test "scuttle: repeated cell subset indices retain R indexing weights" {
  let subset = @src.ScuttleNamedCellSubset::create("weighted", [0, 0, 1]) catch {
    _ => abort("repeated indices are valid")
  }
  let metrics = @src.scuttle_per_feature_qc([[3.0, 9.0]], subsets=[subset]) catch {
    _ => abort("weighted subset QC should work")
  }
  assert_eq(metrics.subset_means[0], [5.0])
}

///|
test "scuttle: duplicate cell subset names are rejected" {
  let first = @src.ScuttleNamedCellSubset::create("same", [0]) catch {
    _ => abort("valid subset should build")
  }
  let second = @src.ScuttleNamedCellSubset::create("same", [1]) catch {
    _ => abort("valid subset should build")
  }
  let mut raised = false
  ignore(@src.scuttle_per_feature_qc([[1.0, 2.0]], subsets=[first, second])) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "scuttle: cell subset bounds are validated" {
  let subset = @src.ScuttleNamedCellSubset::create("bad", [2]) catch {
    _ => abort("bounds are checked with data")
  }
  let mut raised = false
  ignore(@src.scuttle_per_feature_qc([[1.0, 2.0]], subsets=[subset])) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "scuttle: per-feature QC rejects ragged matrices" {
  let mut raised = false
  ignore(@src.scuttle_per_feature_qc([[1.0, 2.0], [3.0]])) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "scuttle: per-feature QC rejects negative counts" {
  let mut raised = false
  ignore(@src.scuttle_per_feature_qc([[1.0, -1.0]])) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "scuttle: arbitrary feature sets are summed per cell" {
  let first = @src.ScuttleFeatureSet::create("first", [0, 2]) catch {
    _ => abort("valid feature set should build")
  }
  let second = @src.ScuttleFeatureSet::create("second", [1]) catch {
    _ => abort("valid feature set should build")
  }
  let result = @src.scuttle_aggregate_feature_sets(
    [[1.0, 2.0], [10.0, 20.0], [3.0, 4.0]],
    [first, second],
  ) catch {
    _ => abort("feature aggregation should work")
  }
  assert_eq(result.names, ["first", "second"])
  assert_eq(result.sizes, [2, 1])
  assert_eq(result.values, [[4.0, 6.0], [10.0, 20.0]])
}

///|
test "scuttle: overlapping feature sets are independent" {
  let left = @src.ScuttleFeatureSet::create("left", [0, 1]) catch {
    _ => abort("valid feature set should build")
  }
  let right = @src.ScuttleFeatureSet::create("right", [1, 2]) catch {
    _ => abort("valid feature set should build")
  }
  let result = @src.scuttle_aggregate_feature_sets([[1.0], [2.0], [4.0]], [
    left, right,
  ]) catch {
    _ => abort("overlapping aggregation should work")
  }
  assert_eq(result.values, [[3.0], [6.0]])
}

///|
test "scuttle: feature-set averages divide by set size" {
  let set = @src.ScuttleFeatureSet::create("mean", [0, 1]) catch {
    _ => abort("valid feature set should build")
  }
  let result = @src.scuttle_aggregate_feature_sets(
    [[2.0, 4.0], [4.0, 8.0]],
    [set],
    average=true,
  ) catch {
    _ => abort("average aggregation should work")
  }
  assert_eq(result.values, [[3.0, 6.0]])
}

///|
test "scuttle: detected features can be aggregated across sets" {
  let set = @src.ScuttleFeatureSet::create("detected", [0, 1, 2]) catch {
    _ => abort("valid feature set should build")
  }
  let result = @src.scuttle_aggregate_feature_sets(
    [[0.0, 2.0], [3.0, 1.0], [4.0, 0.0]],
    [set],
    detection_limit=Some(1.0),
  ) catch {
    _ => abort("detection aggregation should work")
  }
  assert_eq(result.values, [[2.0, 1.0]])
}

///|
test "scuttle: feature identifiers are aggregated in factor order" {
  let result = @src.scuttle_aggregate_features_by_ids(
    [[1.0], [2.0], [4.0], [8.0]],
    ["z", "a", "z", "m"],
  ) catch {
    _ => abort("identifier aggregation should work")
  }
  assert_eq(result.names, ["a", "m", "z"])
  assert_eq(result.values, [[2.0], [8.0], [5.0]])
}

///|
test "scuttle: empty feature identifiers are ignored" {
  let result = @src.scuttle_aggregate_features_by_ids([[1.0], [2.0], [4.0]], [
    "A", "", "A",
  ]) catch {
    _ => abort("empty identifiers should be ignored")
  }
  assert_eq(result.names, ["A"])
  assert_eq(result.values, [[5.0]])
}

///|
test "scuttle: feature identifier length is validated" {
  let mut raised = false
  ignore(@src.scuttle_aggregate_features_by_ids([[1.0], [2.0]], ["A"])) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "scuttle: duplicate feature-set names are rejected" {
  let first = @src.ScuttleFeatureSet::create("same", [0]) catch {
    _ => abort("valid feature set should build")
  }
  let second = @src.ScuttleFeatureSet::create("same", [1]) catch {
    _ => abort("valid feature set should build")
  }
  let mut raised = false
  ignore(@src.scuttle_aggregate_feature_sets([[1.0], [2.0]], [first, second])) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "scuttle: feature-set bounds are validated" {
  let set = @src.ScuttleFeatureSet::create("bad", [2]) catch {
    _ => abort("bounds are checked with data")
  }
  let mut raised = false
  ignore(@src.scuttle_aggregate_feature_sets([[1.0], [2.0]], [set])) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "scuttle: per-column downsampling has exact rounded totals" {
  let counts = scuttle_test_counts()
  let output = @src.scuttle_downsample_matrix(counts, 0.5, seed=42) catch {
    _ => abort("downsampling should work")
  }
  assert_eq(scuttle_test_column_sums(output), [10.0, 10.0, 20.0, 20.0])
}

///|
test "scuttle: downsampled entries never exceed sanitized inputs" {
  let counts = scuttle_test_counts()
  let output = @src.scuttle_downsample_matrix(counts, 0.75, seed=11) catch {
    _ => abort("downsampling should work")
  }
  for feature in 0..<counts.length() {
    for cell in 0..<counts[0].length() {
      assert_true(output[feature][cell] <= counts[feature][cell])
      assert_true(output[feature][cell] >= 0.0)
    }
  }
}

///|
test "scuttle: per-cell proportions can differ" {
  let counts = scuttle_test_counts()
  let output = @src.scuttle_downsample_columns(
    counts,
    [0.0, 0.25, 0.5, 1.0],
    seed=7,
  ) catch {
    _ => abort("vector downsampling should work")
  }
  assert_eq(scuttle_test_column_sums(output), [0.0, 5.0, 20.0, 40.0])
}

///|
test "scuttle: global downsampling has an exact matrix total" {
  let counts = scuttle_test_counts()
  let output = @src.scuttle_downsample_matrix(
    counts,
    0.3,
    by_column=false,
    seed=9,
  ) catch {
    _ => abort("global downsampling should work")
  }
  assert_eq(scuttle_test_total(output), 36.0)
}

///|
test "scuttle: downsampling rounds and clamps numeric inputs" {
  let output = @src.scuttle_downsample_matrix([[-2.0, 1.4], [3.6, 2.6]], 1.0) catch {
    _ => abort("input sanitization should work")
  }
  assert_eq(output, [[0.0, 1.0], [4.0, 3.0]])
}

///|
test "scuttle: equal seeds produce identical samples" {
  let counts = scuttle_test_counts()
  let first = @src.scuttle_downsample_matrix(counts, 0.45, seed=123) catch {
    _ => abort("first downsampling should work")
  }
  let second = @src.scuttle_downsample_matrix(counts, 0.45, seed=123) catch {
    _ => abort("second downsampling should work")
  }
  assert_eq(first, second)
}

///|
test "scuttle: different seeds change nontrivial samples" {
  let counts = scuttle_test_counts()
  let first = @src.scuttle_downsample_matrix(counts, 0.45, seed=1) catch {
    _ => abort("first downsampling should work")
  }
  let second = @src.scuttle_downsample_matrix(counts, 0.45, seed=2) catch {
    _ => abort("second downsampling should work")
  }
  assert_true(first != second)
}

///|
test "scuttle: zero and unit proportions are exact boundaries" {
  let counts = scuttle_test_counts()
  let zero = @src.scuttle_downsample_matrix(counts, 0.0) catch {
    _ => abort("zero downsampling should work")
  }
  let one = @src.scuttle_downsample_matrix(counts, 1.0) catch {
    _ => abort("unit downsampling should work")
  }
  assert_eq(scuttle_test_total(zero), 0.0)
  assert_eq(one, counts)
}

///|
test "scuttle: invalid downsampling proportions are rejected" {
  let mut raised = false
  ignore(@src.scuttle_downsample_matrix([[1.0]], 1.1)) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "scuttle: downsampling validates vector length and seed" {
  let mut bad_length = false
  ignore(@src.scuttle_downsample_columns([[1.0, 2.0]], [0.5])) catch {
    _ => bad_length = true
  }
  assert_true(bad_length)
  let mut bad_seed = false
  ignore(@src.scuttle_downsample_matrix([[1.0]], 0.5, seed=0)) catch {
    _ => bad_seed = true
  }
  assert_true(bad_seed)
}

///|
test "scuttle: batch downsampling targets the shallowest median" {
  let counts = [[5.0, 5.0, 10.0, 10.0], [5.0, 5.0, 10.0, 10.0]]
  let result = @src.scuttle_downsample_batches(
    counts,
    ["A", "A", "B", "B"],
    seed=4,
  ) catch {
    _ => abort("batch downsampling should work")
  }
  assert_eq(result.group_names, ["A", "B"])
  assert_eq(result.summaries, [10.0, 20.0])
  assert_eq(result.proportions, [1.0, 0.5])
  assert_eq(scuttle_test_column_sums(result.matrix), [10.0, 10.0, 10.0, 10.0])
}

///|
test "scuttle: global batch downsampling equalizes group totals" {
  let counts = [[5.0, 5.0, 10.0, 10.0], [5.0, 5.0, 10.0, 10.0]]
  let result = @src.scuttle_downsample_batches(
    counts,
    ["A", "A", "B", "B"],
    by_column=false,
    seed=4,
  ) catch {
    _ => abort("global batch downsampling should work")
  }
  let sums = scuttle_test_column_sums(result.matrix)
  assert_eq(sums[0] + sums[1], 20.0)
  assert_eq(sums[2] + sums[3], 20.0)
}

///|
test "scuttle: blocking equalizes batches only within each block" {
  let counts = [[10.0, 20.0, 100.0, 200.0]]
  let result = @src.scuttle_downsample_batches(
    counts,
    ["A", "B", "A", "B"],
    blocks=["sample1", "sample1", "sample2", "sample2"],
    seed=8,
  ) catch {
    _ => abort("blocked batch downsampling should work")
  }
  assert_eq(result.group_names, [
    "sample1|A", "sample1|B", "sample2|A", "sample2|B",
  ])
  assert_eq(result.proportions, [1.0, 0.5, 1.0, 0.5])
  assert_eq(scuttle_test_column_sums(result.matrix), [10.0, 10.0, 100.0, 100.0])
}

///|
test "scuttle: batch mean summary is selectable" {
  let counts = [[2.0, 10.0, 12.0, 24.0]]
  let result = @src.scuttle_downsample_batches(
    counts,
    ["A", "A", "B", "B"],
    summary_method=@src.scuttle_batch_mean(),
    seed=10,
  ) catch {
    _ => abort("mean batch downsampling should work")
  }
  assert_eq(result.summaries, [6.0, 18.0])
  assert_true(scuttle_test_close(result.proportions[1], 1.0 / 3.0, 1.0e-12))
}

///|
test "scuttle: geometric batch summary uses a pseudocount" {
  let counts = [[0.0, 3.0, 8.0, 8.0]]
  let result = @src.scuttle_downsample_batches(
    counts,
    ["A", "A", "B", "B"],
    summary_method=@src.scuttle_batch_geometric_mean(),
    seed=10,
  ) catch {
    _ => abort("geometric batch downsampling should work")
  }
  assert_eq(result.summaries[0], 2.0)
  assert_true(scuttle_test_close(result.summaries[1], 9.0, 1.0e-10))
}

///|
test "scuttle: batch downsampling validates annotations" {
  let mut bad_batch = false
  ignore(@src.scuttle_downsample_batches([[1.0, 2.0]], ["A"])) catch {
    _ => bad_batch = true
  }
  assert_true(bad_batch)
  let mut bad_block = false
  ignore(
    @src.scuttle_downsample_batches([[1.0, 2.0]], ["A", "B"], blocks=["one"]),
  ) catch {
    _ => bad_block = true
  }
  assert_true(bad_block)
}

///|
test "scuttle: feature QC SCE wrapper writes row metadata immutably" {
  let sce = scuttle_test_sce()
  let controls = @src.ScuttleNamedCellSubset::create("controls", [0, 1]) catch {
    _ => abort("valid subset should build")
  }
  let output = @src.scuttle_per_feature_qc_sce(sce, subsets=[controls]) catch {
    _ => abort("SCE feature QC should work")
  }
  assert_eq(@src.sce_get_row_data(sce, "scuttle.mean").length(), 0)
  assert_eq(
    @src.sce_get_row_data(output.experiment, "scuttle.mean").length(),
    4,
  )
  assert_eq(
    @src.sce_get_row_data(output.experiment, "scuttle.subsets.controls.ratio").length(),
    4,
  )
  output.experiment.assays["counts"][0][0] = 999.0
  assert_eq(sce.assays["counts"][0][0], 8.0)
}

///|
test "scuttle: feature QC SCE wrapper validates assay names" {
  let sce = scuttle_test_sce()
  let mut raised = false
  ignore(@src.scuttle_per_feature_qc_sce(sce, assay_name="missing")) catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "scuttle: SCE feature aggregation handles multiple assays" {
  let sce = scuttle_test_sce()
  let first = @src.ScuttleFeatureSet::create("PathwayA", [0, 1]) catch {
    _ => abort("valid feature set should build")
  }
  let second = @src.ScuttleFeatureSet::create("PathwayB", [2, 3]) catch {
    _ => abort("valid feature set should build")
  }
  let output = @src.scuttle_aggregate_feature_sets_sce(sce, [first, second], assay_names=[
    "counts", "logcounts",
  ]) catch {
    _ => abort("SCE feature aggregation should work")
  }
  assert_eq(output.row_names, ["PathwayA", "PathwayB"])
  assert_eq(output.assays["counts"], [
    [10.0, 10.0, 20.0, 20.0],
    [10.0, 10.0, 20.0, 20.0],
  ])
  assert_eq(output.assays["logcounts"].length(), 2)
  assert_eq(output.col_data["batch"], ["A", "A", "B", "B"])
  assert_eq(output.reduced_dims["PCA"].length(), 4)
  assert_eq(output.alternative_experiments["spike"].row_names, ["Spike1"])
  assert_eq(output.metadata["project"], "scuttle-test")
}

///|
test "scuttle: SCE feature aggregation discards gene metadata and deep copies" {
  let sce = scuttle_test_sce()
  let set = @src.ScuttleFeatureSet::create("all", [0, 1, 2, 3]) catch {
    _ => abort("valid feature set should build")
  }
  let output = @src.scuttle_aggregate_feature_sets_sce(sce, [set]) catch {
    _ => abort("SCE feature aggregation should work")
  }
  assert_eq(output.row_data.keys().collect().length(), 0)
  output.col_data["batch"][0] = "changed"
  output.reduced_dims["PCA"][0][0] = 999.0
  assert_eq(sce.col_data["batch"][0], "A")
  assert_eq(sce.reduced_dims["PCA"][0][0], 1.0)
}

///|
test "scuttle: SCE downsampling adds an assay without mutating input" {
  let sce = scuttle_test_sce()
  let output = @src.scuttle_downsample_sce(
    sce,
    0.5,
    output_assay="half",
    seed=42,
  ) catch {
    _ => abort("SCE downsampling should work")
  }
  assert_eq(@src.sce_get_assay(sce, "half").length(), 0)
  assert_eq(scuttle_test_column_sums(@src.sce_get_assay(output, "half")), [
    10.0, 10.0, 20.0, 20.0,
  ])
  assert_eq(output.metadata["scuttle.downsample.source"], "counts")
  output.assays["counts"][0][0] = 999.0
  assert_eq(sce.assays["counts"][0][0], 8.0)
  assert_eq(output.alternative_experiments["spike"].row_names, ["Spike1"])
}

///|
test "scuttle: SCE downsampling validates source assay" {
  let sce = scuttle_test_sce()
  let mut raised = false
  ignore(@src.scuttle_downsample_sce(sce, 0.5, assay_name="missing")) catch {
    _ => raised = true
  }
  assert_true(raised)
}
