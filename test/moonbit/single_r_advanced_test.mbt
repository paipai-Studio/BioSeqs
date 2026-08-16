///|
fn single_r_adv_test_example() -> (
  @src.SingleRAdvancedReference,
  @src.SingleRAdvancedData,
) {
  @src.single_r_advanced_example() catch {
    _ => abort("valid SingleR advanced example should build")
  }
}

///|
fn single_r_adv_test_result() -> @src.SingleRAdvancedResult {
  let (reference, data) = single_r_adv_test_example()
  @src.single_r_advanced(data, reference) catch {
    _ => abort("valid SingleR advanced example should classify")
  }
}

///|
fn single_r_adv_test_sce() -> @src.SingleCellExperiment {
  let (_, data) = single_r_adv_test_example()
  let experiment = @src.SingleCellExperiment::new(
    data.expression.map(fn(row) { row.copy() }),
    data.gene_names.copy(),
    data.cell_names.copy(),
  )
  experiment.assays["logcounts"] = data.expression.map(fn(row) { row.copy() })
  experiment.col_data["batch"] = ["x", "x", "y", "y", "y"]
  experiment.row_data["symbol"] = data.gene_names.copy()
  experiment.reduced_dims["PCA"] = [
    [0.0, 1.0],
    [1.0, 0.0],
    [0.5, 0.5],
    [0.2, 0.8],
    [0.8, 0.2],
  ]
  experiment.metadata["source"] = "test"
  experiment.alternative_experiments["spike"] = @src.SingleCellExperiment::new(
    [[1.0, 2.0, 3.0, 4.0, 5.0]],
    ["spike1"],
    data.cell_names.copy(),
  )
  experiment
}

///|
fn single_r_adv_test_second_reference() -> @src.SingleRAdvancedReference {
  let (reference, _) = single_r_adv_test_example()
  @src.SingleRAdvancedReference::create(
    reference.expression,
    [
      "T lymphocyte", "T lymphocyte", "B lymphocyte", "B lymphocyte", "Myeloid",
      "Myeloid", "Cytotoxic", "Cytotoxic",
    ],
    gene_names=reference.gene_names,
    sample_names=reference.sample_names,
    reference_name="immune_reference_two",
  ) catch {
    _ => abort("valid second SingleR reference should build")
  }
}

///|
test "SingleR advanced default configuration matches upstream controls" {
  let config = @src.SingleRAdvancedConfig::default()
  assert_eq(config.quantile, 0.8)
  assert_true(config.fine_tune)
  assert_eq(config.tune_threshold, 0.05)
  assert_true(config.prune)
  assert_eq(config.nmads, 3.0)
  assert_eq(config.minimum_delta_next, 0.0)
  assert_eq(config.marker_count, 0)
  assert_eq(config.combine_marker_count, 10)
}

///|
test "SingleR advanced configuration preserves explicit controls" {
  let config = @src.SingleRAdvancedConfig::create(
    quantile=0.6,
    fine_tune=false,
    tune_threshold=0.1,
    prune=false,
    nmads=2.0,
    minimum_delta_median=0.05,
    minimum_delta_next=0.02,
    marker_count=3,
    combine_marker_count=4,
  ) catch {
    _ => abort("valid SingleR advanced configuration should build")
  }
  assert_eq(config.quantile, 0.6)
  assert_eq(config.fine_tune, false)
  assert_eq(config.tune_threshold, 0.1)
  assert_eq(config.prune, false)
  assert_eq(config.nmads, 2.0)
  assert_eq(config.minimum_delta_median, 0.05)
  assert_eq(config.minimum_delta_next, 0.02)
  assert_eq(config.marker_count, 3)
  assert_eq(config.combine_marker_count, 4)
}

///|
test "SingleR advanced rejects invalid quantiles" {
  let lower = try {
    ignore(@src.SingleRAdvancedConfig::create(quantile=-0.1))
    false
  } catch {
    _ => true
  }
  let upper = try {
    ignore(@src.SingleRAdvancedConfig::create(quantile=1.1))
    false
  } catch {
    _ => true
  }
  assert_true(lower)
  assert_true(upper)
}

///|
test "SingleR advanced rejects invalid tuning and pruning controls" {
  let tuning = try {
    ignore(@src.SingleRAdvancedConfig::create(tune_threshold=-0.1))
    false
  } catch {
    _ => true
  }
  let nmads = try {
    ignore(@src.SingleRAdvancedConfig::create(nmads=-1.0))
    false
  } catch {
    _ => true
  }
  assert_true(tuning)
  assert_true(nmads)
}

///|
test "SingleR advanced rejects invalid marker controls" {
  let markers = try {
    ignore(@src.SingleRAdvancedConfig::create(marker_count=-1))
    false
  } catch {
    _ => true
  }
  let combine = try {
    ignore(@src.SingleRAdvancedConfig::create(combine_marker_count=0))
    false
  } catch {
    _ => true
  }
  assert_true(markers)
  assert_true(combine)
}

///|
test "SingleR advanced reference uses gene by sample orientation" {
  let (reference, _) = single_r_adv_test_example()
  assert_eq(reference.n_genes, 12)
  assert_eq(reference.n_samples, 8)
  assert_eq(reference.expression[0].length(), 8)
  assert_eq(reference.gene_names[0], "CD3D")
  assert_eq(reference.sample_names[7], "NK2")
}

///|
test "SingleR advanced preserves first occurrence label order" {
  let (reference, _) = single_r_adv_test_example()
  assert_eq(reference.label_names, ["T cell", "B cell", "Monocyte", "NK cell"])
}

///|
test "SingleR advanced reference defensively copies inputs" {
  let expression = [[1.0, 2.0], [3.0, 4.0]]
  let labels = ["A", "B"]
  let genes = ["g1", "g2"]
  let reference = @src.SingleRAdvancedReference::create(
    expression,
    labels,
    gene_names=genes,
  ) catch {
    _ => abort("valid reference should build")
  }
  expression[0][0] = 99.0
  labels[0] = "changed"
  genes[0] = "changed"
  assert_eq(reference.expression[0][0], 1.0)
  assert_eq(reference.labels[0], "A")
  assert_eq(reference.gene_names[0], "g1")
}

///|
test "SingleR advanced accepts finite negative normalized expression" {
  let reference = @src.SingleRAdvancedReference::create(
    [[-2.0, 1.0], [0.0, -1.0]],
    ["A", "B"],
  ) catch {
    _ => abort("finite normalized expression may be negative")
  }
  assert_eq(reference.expression[0][0], -2.0)
}

///|
test "SingleR advanced data uses gene by cell orientation" {
  let (_, data) = single_r_adv_test_example()
  assert_eq(data.n_genes, 12)
  assert_eq(data.n_cells, 5)
  assert_eq(data.expression[0].length(), 5)
  assert_eq(data.cell_names[4], "cell_ambiguous")
}

///|
test "SingleR advanced data defensively copies inputs" {
  let expression = [[1.0, 2.0], [3.0, 4.0]]
  let data = @src.SingleRAdvancedData::create(
    expression,
    gene_names=["g1", "g2"],
    cell_names=["c1", "c2"],
  ) catch {
    _ => abort("valid data should build")
  }
  expression[0][0] = 99.0
  let copied = data.copy_expression()
  copied[0][0] = 77.0
  assert_eq(data.expression[0][0], 1.0)
}

///|
test "SingleR advanced constructors reject empty matrices" {
  let reference = try {
    ignore(@src.SingleRAdvancedReference::create([], []))
    false
  } catch {
    _ => true
  }
  let data = try {
    ignore(@src.SingleRAdvancedData::create([[]]))
    false
  } catch {
    _ => true
  }
  assert_true(reference)
  assert_true(data)
}

///|
test "SingleR advanced constructors reject ragged matrices" {
  let reference = try {
    ignore(
      @src.SingleRAdvancedReference::create([[1.0, 2.0], [3.0]], ["A", "B"]),
    )
    false
  } catch {
    _ => true
  }
  let data = try {
    ignore(@src.SingleRAdvancedData::create([[1.0, 2.0], [3.0]]))
    false
  } catch {
    _ => true
  }
  assert_true(reference)
  assert_true(data)
}

///|
test "SingleR advanced constructors reject non-finite values" {
  let reference = try {
    ignore(
      @src.SingleRAdvancedReference::create([[1.0, 1.0e301], [2.0, 3.0]], [
        "A", "B",
      ]),
    )
    false
  } catch {
    _ => true
  }
  let data = try {
    ignore(@src.SingleRAdvancedData::create([[0.0 / 0.0, 1.0]]))
    false
  } catch {
    _ => true
  }
  assert_true(reference)
  assert_true(data)
}

///|
test "SingleR advanced reference validates label dimensions and content" {
  let dimensions = try {
    ignore(
      @src.SingleRAdvancedReference::create([[1.0, 2.0], [3.0, 4.0]], ["A"]),
    )
    false
  } catch {
    _ => true
  }
  let content = try {
    ignore(
      @src.SingleRAdvancedReference::create([[1.0, 2.0], [3.0, 4.0]], ["A", " "]),
    )
    false
  } catch {
    _ => true
  }
  assert_true(dimensions)
  assert_true(content)
}

///|
test "SingleR advanced constructors validate name dimensions" {
  let genes = try {
    ignore(
      @src.SingleRAdvancedData::create([[1.0, 2.0], [3.0, 4.0]], gene_names=[
        "g1",
      ]),
    )
    false
  } catch {
    _ => true
  }
  let samples = try {
    ignore(
      @src.SingleRAdvancedReference::create(
        [[1.0, 2.0], [3.0, 4.0]],
        ["A", "B"],
        sample_names=["s1"],
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(genes)
  assert_true(samples)
}

///|
test "SingleR advanced constructors reject duplicate names" {
  let genes = try {
    ignore(
      @src.SingleRAdvancedData::create([[1.0, 2.0], [3.0, 4.0]], gene_names=[
        "g", "g",
      ]),
    )
    false
  } catch {
    _ => true
  }
  let cells = try {
    ignore(
      @src.SingleRAdvancedData::create([[1.0, 2.0], [3.0, 4.0]], cell_names=[
        "c", "c",
      ]),
    )
    false
  } catch {
    _ => true
  }
  assert_true(genes)
  assert_true(cells)
}

///|
test "SingleR advanced Spearman handles positive negative and tied ranks" {
  let positive = @src.single_r_advanced_spearman([1.0, 2.0, 2.0, 4.0], [
    10.0, 20.0, 20.0, 40.0,
  ]) catch {
    _ => abort("valid tied Spearman vectors should work")
  }
  let negative = @src.single_r_advanced_spearman([1.0, 2.0, 3.0, 4.0], [
    4.0, 3.0, 2.0, 1.0,
  ]) catch {
    _ => abort("valid Spearman vectors should work")
  }
  assert_true(positive > 0.999)
  assert_true(negative < -0.999)
}

///|
test "SingleR advanced Spearman validates dimensions and finite values" {
  let dimensions = try {
    ignore(@src.single_r_advanced_spearman([1.0], [1.0, 2.0]))
    false
  } catch {
    _ => true
  }
  let finite = try {
    ignore(@src.single_r_advanced_spearman([1.0, 1.0e301], [1.0, 2.0]))
    false
  } catch {
    _ => true
  }
  assert_true(dimensions)
  assert_true(finite)
}

///|
test "SingleR advanced training aligns genes in test order" {
  let (reference, _) = single_r_adv_test_example()
  let training = @src.single_r_advanced_train(
    reference,
    ["GNLY", "missing", "CD3D", "MS4A1"],
    marker_count=1,
  ) catch {
    _ => abort("gene intersection should train")
  }
  assert_eq(training.common_gene_names, ["GNLY", "CD3D", "MS4A1"])
  assert_eq(training.test_gene_indices, [0, 2, 3])
  assert_eq(training.reference_gene_indices, [11, 0, 3])
}

///|
test "SingleR advanced training honors gene restrictions" {
  let (reference, data) = single_r_adv_test_example()
  let training = @src.single_r_advanced_train(
    reference,
    data.gene_names,
    restrict_genes=["CD3D", "MS4A1", "LYZ", "NKG7"],
    marker_count=1,
  ) catch {
    _ => abort("restricted training should work")
  }
  assert_eq(training.common_gene_names, ["CD3D", "MS4A1", "LYZ", "NKG7"])
}

///|
test "SingleR advanced training rejects insufficient shared genes" {
  let (reference, _) = single_r_adv_test_example()
  let raised = try {
    ignore(@src.single_r_advanced_train(reference, ["missing", "CD3D"]))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "SingleR advanced training computes automatic marker count" {
  let (reference, data) = single_r_adv_test_example()
  let training = @src.single_r_advanced_train(reference, data.gene_names) catch {
    _ => abort("automatic marker training should work")
  }
  assert_true(training.marker_count > 0)
  assert_true(training.common_markers.length() >= 2)
}

///|
test "SingleR advanced classic markers are directional" {
  let (reference, data) = single_r_adv_test_example()
  let training = @src.single_r_advanced_train(
    reference,
    data.gene_names,
    marker_count=2,
  ) catch {
    _ => abort("marker training should work")
  }
  let t_markers = training.markers_between("T cell", "B cell") catch {
    _ => abort("known labels should have markers")
  }
  let b_markers = training.markers_between("B cell", "T cell") catch {
    _ => abort("known labels should have markers")
  }
  assert_true(t_markers.contains("CD3D") || t_markers.contains("CD3E"))
  assert_true(b_markers.contains("MS4A1") || b_markers.contains("CD79A"))
}

///|
test "SingleR advanced self comparison has no markers" {
  let (reference, data) = single_r_adv_test_example()
  let training = @src.single_r_advanced_train(reference, data.gene_names) catch {
    _ => abort("marker training should work")
  }
  assert_eq(
    training.markers_between("T cell", "T cell") catch {
      _ => abort("known labels should query")
    },
    [],
  )
}

///|
test "SingleR advanced marker query rejects unknown labels" {
  let (reference, data) = single_r_adv_test_example()
  let training = @src.single_r_advanced_train(reference, data.gene_names) catch {
    _ => abort("marker training should work")
  }
  let raised = try {
    ignore(training.markers_between("missing", "T cell"))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "SingleR advanced classification locks training gene order" {
  let (reference, data) = single_r_adv_test_example()
  let training = @src.single_r_advanced_train(reference, data.gene_names) catch {
    _ => abort("training should work")
  }
  let reordered = @src.SingleRAdvancedData::create(
    data.expression,
    gene_names=[
      "CD3E", "CD3D", "TRAC", "MS4A1", "CD79A", "CD74", "LYZ", "S100A8", "FCGR3A",
      "NKG7", "KLRD1", "GNLY",
    ],
    cell_names=data.cell_names,
  ) catch {
    _ => abort("reordered data should build")
  }
  let raised = try {
    ignore(@src.single_r_advanced_classify(reordered, training))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "SingleR advanced classification reports cell by label scores" {
  let result = single_r_adv_test_result()
  assert_eq(result.n_cells(), 5)
  assert_eq(result.scores.length(), 5)
  assert_eq(result.scores[0].length(), 4)
  assert_eq(result.label_names.length(), 4)
}

///|
test "SingleR advanced classifies canonical immune profiles" {
  let result = single_r_adv_test_result()
  assert_eq(result.labels[0], "T cell")
  assert_eq(result.labels[1], "B cell")
  assert_eq(result.labels[2], "Monocyte")
  assert_eq(result.labels[3], "NK cell")
}

///|
test "SingleR advanced preserves pre-fine-tuning assignments" {
  let result = single_r_adv_test_result()
  assert_eq(result.first_labels[0], "T cell")
  assert_eq(result.first_labels[1], "B cell")
  assert_eq(result.first_labels[2], "Monocyte")
  assert_eq(result.first_labels[3], "NK cell")
}

///|
test "SingleR advanced quantile controls label aggregation" {
  let (reference, data) = single_r_adv_test_example()
  let low = @src.SingleRAdvancedConfig::create(
    quantile=0.0,
    fine_tune=false,
    prune=false,
  ) catch {
    _ => abort("low quantile configuration should build")
  }
  let high = @src.SingleRAdvancedConfig::create(
    quantile=1.0,
    fine_tune=false,
    prune=false,
  ) catch {
    _ => abort("high quantile configuration should build")
  }
  let low_result = @src.single_r_advanced(data, reference, config=low) catch {
    _ => abort("low quantile classification should work")
  }
  let high_result = @src.single_r_advanced(data, reference, config=high) catch {
    _ => abort("high quantile classification should work")
  }
  assert_true(high_result.scores[0][0] >= low_result.scores[0][0])
  assert_true(high_result.scores[1][1] >= low_result.scores[1][1])
}

///|
test "SingleR advanced fine tuning is active for ambiguous candidates" {
  let (reference, data) = single_r_adv_test_example()
  let config = @src.SingleRAdvancedConfig::create(
    tune_threshold=2.0,
    prune=false,
  ) catch {
    _ => abort("wide fine-tuning configuration should build")
  }
  let result = @src.single_r_advanced(data, reference, config~) catch {
    _ => abort("fine-tuned classification should work")
  }
  assert_true(result.fine_tuned[4])
  assert_true(result.delta_next[4] >= 0.0)
}

///|
test "SingleR advanced can disable fine tuning" {
  let (reference, data) = single_r_adv_test_example()
  let config = @src.SingleRAdvancedConfig::create(fine_tune=false, prune=false) catch {
    _ => abort("no-fine-tune configuration should build")
  }
  let result = @src.single_r_advanced(data, reference, config~) catch {
    _ => abort("classification without fine tuning should work")
  }
  for value in result.fine_tuned {
    assert_eq(value, false)
  }
  assert_eq(result.labels, result.first_labels)
}

///|
test "SingleR advanced delta from median matches score matrix" {
  let result = single_r_adv_test_result()
  let label = result.label_index(result.labels[0])
  let sorted = result.scores[0]
  let median = if sorted.length() == 4 {
    let values = sorted.copy()
    for index in 1..<values.length() {
      let current = values[index]
      let mut position = index
      while position > 0 && values[position - 1] > current {
        values[position] = values[position - 1]
        position = position - 1
      }
      values[position] = current
    }
    (values[1] + values[2]) / 2.0
  } else {
    0.0
  }
  assert_true(
    (result.delta_median[0] - (result.scores[0][label] - median)).abs() <
    1.0e-12,
  )
}

///|
test "SingleR advanced hard delta pruning uses NA semantics" {
  let (reference, data) = single_r_adv_test_example()
  let config = @src.SingleRAdvancedConfig::create(minimum_delta_next=3.0) catch {
    _ => abort("strict pruning configuration should build")
  }
  let result = @src.single_r_advanced(data, reference, config~) catch {
    _ => abort("strict pruning classification should work")
  }
  assert_eq(result.n_pruned(), data.n_cells)
  for value in result.pruned_labels {
    assert_true(value is None)
  }
}

///|
test "SingleR advanced can disable pruning" {
  let (reference, data) = single_r_adv_test_example()
  let config = @src.SingleRAdvancedConfig::create(prune=false) catch {
    _ => abort("unpruned configuration should build")
  }
  let result = @src.single_r_advanced(data, reference, config~) catch {
    _ => abort("unpruned classification should work")
  }
  assert_eq(result.n_pruned(), 0)
  for value in result.pruned_labels {
    assert_true(value is Some(_))
  }
}

///|
test "SingleR advanced computes per-label MAD thresholds" {
  let result = single_r_adv_test_result()
  let mut observed = 0
  for threshold in result.pruning_thresholds {
    if threshold is Some(_) {
      observed = observed + 1
    }
  }
  assert_true(observed > 0)
}

///|
test "SingleR advanced result query helpers are consistent" {
  let result = single_r_adv_test_result()
  assert_eq(result.label_index("T cell"), 0)
  assert_eq(result.label_index("missing"), -1)
  assert_true(result.assigned_score(0) > 0.0)
  let summary = result.summary()
  assert_eq(summary["total"], 5)
  assert_eq(
    summary["T cell"] +
    summary["B cell"] +
    summary["Monocyte"] +
    summary["NK cell"],
    5,
  )
}

///|
test "SingleR advanced assigned score validates cell indices" {
  let result = single_r_adv_test_result()
  let raised = try {
    ignore(result.assigned_score(99))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "SingleR advanced cluster annotation sums member profiles" {
  let (reference, data) = single_r_adv_test_example()
  let config = @src.SingleRAdvancedConfig::create(prune=false) catch {
    _ => abort("cluster configuration should build")
  }
  let result = @src.single_r_advanced_clusters(
    data,
    ["lymphoid", "lymphoid", "myeloid", "cytotoxic", "mixed"],
    reference,
    config~,
  ) catch {
    _ => abort("cluster annotation should work")
  }
  assert_eq(result.cell_names, ["lymphoid", "myeloid", "cytotoxic", "mixed"])
  assert_eq(result.n_cells(), 4)
  assert_eq(result.labels[1], "Monocyte")
  assert_eq(result.labels[2], "NK cell")
}

///|
test "SingleR advanced cluster annotation preserves first occurrence order" {
  let (reference, data) = single_r_adv_test_example()
  let result = @src.single_r_advanced_clusters(
    data,
    ["z", "a", "z", "b", "a"],
    reference,
  ) catch {
    _ => abort("cluster annotation should work")
  }
  assert_eq(result.cell_names, ["z", "a", "b"])
}

///|
test "SingleR advanced cluster annotation validates assignments" {
  let (reference, data) = single_r_adv_test_example()
  let dimensions = try {
    ignore(@src.single_r_advanced_clusters(data, ["a"], reference))
    false
  } catch {
    _ => true
  }
  let names = try {
    ignore(
      @src.single_r_advanced_clusters(data, ["a", "a", "b", "b", ""], reference),
    )
    false
  } catch {
    _ => true
  }
  assert_true(dimensions)
  assert_true(names)
}

///|
test "SingleR advanced multi-reference recomputation returns comparable scores" {
  let (reference, data) = single_r_adv_test_example()
  let second = single_r_adv_test_second_reference()
  let result = @src.single_r_advanced_combine(data, [reference, second]) catch {
    _ => abort("multi-reference classification should work")
  }
  assert_eq(result.n_cells(), 5)
  assert_eq(result.scores.length(), 5)
  assert_eq(result.scores[0].length(), 2)
  assert_eq(result.per_reference.length(), 2)
  assert_true(result.marker_gene_names.length() >= 2)
}

///|
test "SingleR advanced multi-reference records selected provenance" {
  let (reference, data) = single_r_adv_test_example()
  let second = single_r_adv_test_second_reference()
  let result = @src.single_r_advanced_combine(data, [reference, second]) catch {
    _ => abort("multi-reference classification should work")
  }
  assert_eq(result.references, ["immune_reference", "immune_reference_two"])
  for selected in result.selected_references {
    assert_true(result.references.contains(selected))
  }
  let summary = result.reference_summary()
  assert_eq(summary["immune_reference"] + summary["immune_reference_two"], 5)
}

///|
test "SingleR advanced multi-reference rejects too few references" {
  let (reference, data) = single_r_adv_test_example()
  let raised = try {
    ignore(@src.single_r_advanced_combine(data, [reference]))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "SingleR advanced multi-reference rejects duplicate names" {
  let (reference, data) = single_r_adv_test_example()
  let duplicate = @src.SingleRAdvancedReference::create(
    reference.expression,
    reference.labels,
    gene_names=reference.gene_names,
    sample_names=reference.sample_names,
    reference_name=reference.reference_name,
  ) catch {
    _ => abort("duplicate-name reference should build independently")
  }
  let raised = try {
    ignore(@src.single_r_advanced_combine(data, [reference, duplicate]))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "SingleR advanced multi-reference requires global gene intersection" {
  let (reference, data) = single_r_adv_test_example()
  let other = @src.SingleRAdvancedReference::create(
    [[1.0, 2.0], [2.0, 1.0]],
    ["X", "Y"],
    gene_names=["other1", "other2"],
    reference_name="other",
  ) catch {
    _ => abort("disjoint reference should build")
  }
  let raised = try {
    ignore(@src.single_r_advanced_combine(data, [reference, other]))
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "SingleR advanced SingleCellExperiment writes annotation columns" {
  let (reference, _) = single_r_adv_test_example()
  let output = @src.single_r_advanced_sce(
    single_r_adv_test_sce(),
    reference,
    output_prefix="singler215",
  ) catch {
    _ => abort("SingleR SCE annotation should work")
  }
  assert_true(output.experiment.col_data.contains("singler215.labels"))
  assert_true(output.experiment.col_data.contains("singler215.pruned.labels"))
  assert_true(output.experiment.col_data.contains("singler215.score"))
  assert_true(output.experiment.col_data.contains("singler215.delta.next"))
  assert_true(output.experiment.col_data.contains("singler215.delta.median"))
  assert_eq(output.experiment.col_data["singler215.labels"][0], "T cell")
}

///|
test "SingleR advanced SingleCellExperiment records provenance metadata" {
  let (reference, _) = single_r_adv_test_example()
  let output = @src.single_r_advanced_sce(
    single_r_adv_test_sce(),
    reference,
    output_prefix="singler215",
  ) catch {
    _ => abort("SingleR SCE annotation should work")
  }
  assert_eq(output.experiment.metadata["singler215.assay"], "logcounts")
  assert_eq(
    output.experiment.metadata["singler215.reference"],
    "immune_reference",
  )
  assert_eq(output.experiment.metadata["singler215.quantile"], "0.8")
  assert_eq(output.experiment.metadata["singler215.mode"], "cell")
  assert_eq(output.experiment.metadata["source"], "test")
}

///|
test "SingleR advanced SingleCellExperiment does not mutate input" {
  let (reference, _) = single_r_adv_test_example()
  let input = single_r_adv_test_sce()
  let output = @src.single_r_advanced_sce(
    input,
    reference,
    output_prefix="singler215",
  ) catch {
    _ => abort("SingleR SCE annotation should work")
  }
  assert_eq(input.col_data.contains("singler215.labels"), false)
  output.experiment.assays["logcounts"][0][0] = 999.0
  output.experiment.row_data["symbol"][0] = "changed"
  output.experiment.reduced_dims["PCA"][0][0] = 999.0
  output.experiment.alternative_experiments["spike"].assays["counts"][0][0] = 999.0
  assert_true(input.assays["logcounts"][0][0] != 999.0)
  assert_eq(input.row_data["symbol"][0], "CD3D")
  assert_eq(input.reduced_dims["PCA"][0][0], 0.0)
  assert_eq(input.alternative_experiments["spike"].assays["counts"][0][0], 1.0)
}

///|
test "SingleR advanced SingleCellExperiment expands cluster labels to cells" {
  let (reference, _) = single_r_adv_test_example()
  let input = single_r_adv_test_sce()
  input.col_data["cluster"] = ["a", "b", "c", "d", "a"]
  let output = @src.single_r_advanced_sce(
    input,
    reference,
    cluster_column="cluster",
    output_prefix="singler215",
  ) catch {
    _ => abort("cluster SCE annotation should work")
  }
  assert_eq(output.result.n_cells(), 4)
  assert_eq(output.experiment.col_data["singler215.labels"].length(), 5)
  assert_eq(
    output.experiment.col_data["singler215.labels"][0],
    output.experiment.col_data["singler215.labels"][4],
  )
  assert_eq(output.experiment.metadata["singler215.mode"], "cluster")
}

///|
test "SingleR advanced SingleCellExperiment supports custom assay" {
  let (reference, _) = single_r_adv_test_example()
  let output = @src.single_r_advanced_sce(
    single_r_adv_test_sce(),
    reference,
    assay_name="counts",
    output_prefix="custom",
  ) catch {
    _ => abort("custom assay SCE annotation should work")
  }
  assert_true(output.experiment.col_data.contains("custom.labels"))
  assert_eq(output.experiment.metadata["custom.assay"], "counts")
}

///|
test "SingleR advanced SingleCellExperiment rejects missing inputs" {
  let (reference, _) = single_r_adv_test_example()
  let assay = try {
    ignore(
      @src.single_r_advanced_sce(
        single_r_adv_test_sce(),
        reference,
        assay_name="missing",
      ),
    )
    false
  } catch {
    _ => true
  }
  let cluster = try {
    ignore(
      @src.single_r_advanced_sce(
        single_r_adv_test_sce(),
        reference,
        cluster_column="missing",
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(assay)
  assert_true(cluster)
}

///|
test "SingleR advanced SingleCellExperiment validates prefix" {
  let (reference, _) = single_r_adv_test_example()
  let raised = try {
    ignore(
      @src.single_r_advanced_sce(
        single_r_adv_test_sce(),
        reference,
        output_prefix=" ",
      ),
    )
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}
