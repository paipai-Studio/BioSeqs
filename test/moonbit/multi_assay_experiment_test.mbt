///|
test "experiment_basic" {
  let exp = @src.Experiment::new(
    name="RNA-seq",
    assay_type="counts",
    data=[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
    row_names=["gene1", "gene2"],
    col_names=["s1", "s2", "s3"],
  )
  assert_eq(@src.Experiment::nrow(exp), 2)
  assert_eq(@src.Experiment::ncol(exp), 3)
}

///|
test "experiment_empty" {
  let exp = @src.Experiment::new(
    name="empty",
    assay_type="empty",
    data=[],
    row_names=[],
    col_names=[],
  )
  assert_eq(@src.Experiment::nrow(exp), 0)
  assert_eq(@src.Experiment::ncol(exp), 0)
}

///|
test "sample_map_primary_ids" {
  let sm = @src.SampleMap::new(
    assays=["RNA-seq", "RNA-seq", "ChIP-seq"],
    primary_ids=["patient1", "patient2", "patient1"],
    colnames=["rna_s1", "rna_s2", "chip_s1"],
  )
  let ids = @src.SampleMap::primary_ids(sm)
  assert_eq(ids.length(), 2)
  assert_true(ids.contains("patient1"))
  assert_true(ids.contains("patient2"))
}

///|
test "multi_assay_experiment_basic" {
  let rna_exp = @src.Experiment::new(
    name="RNA-seq",
    assay_type="counts",
    data=[[10.0, 20.0], [30.0, 40.0]],
    row_names=["gene1", "gene2"],
    col_names=["rna_s1", "rna_s2"],
  )
  let chip_exp = @src.Experiment::new(
    name="ChIP-seq",
    assay_type="scores",
    data=[[1.5, 2.5]],
    row_names=["peak1"],
    col_names=["chip_s1", "chip_s2"],
  )
  let sm = @src.SampleMap::new(
    assays=["RNA-seq", "RNA-seq", "ChIP-seq", "ChIP-seq"],
    primary_ids=["p1", "p2", "p1", "p2"],
    colnames=["rna_s1", "rna_s2", "chip_s1", "chip_s2"],
  )
  let mae = @src.MultiAssayExperiment::new(
    experiments=[rna_exp, chip_exp],
    sample_map=sm,
    col_data=[],
  )
  assert_eq(@src.MultiAssayExperiment::num_experiments(mae), 2)
  assert_eq(@src.MultiAssayExperiment::num_samples(mae), 2)
}

///|
test "multi_assay_experiment_get_experiment" {
  let rna_exp = @src.Experiment::new(
    name="RNA-seq",
    assay_type="counts",
    data=[[1.0]],
    row_names=["g1"],
    col_names=["s1"],
  )
  let chip_exp = @src.Experiment::new(
    name="ChIP-seq",
    assay_type="scores",
    data=[[2.0]],
    row_names=["p1"],
    col_names=["s1"],
  )
  let sm = @src.SampleMap::new(assays=[], primary_ids=[], colnames=[])
  let mae = @src.MultiAssayExperiment::new(
    experiments=[rna_exp, chip_exp],
    sample_map=sm,
    col_data=[],
  )
  let names = @src.MultiAssayExperiment::experiment_names(mae)
  assert_eq(names.length(), 2)
  assert_true(names.contains("RNA-seq"))
  assert_true(names.contains("ChIP-seq"))
  match @src.MultiAssayExperiment::get_experiment(mae, "RNA-seq") {
    Some(_) => assert_true(true)
    None => assert_true(false)
  }
  match @src.MultiAssayExperiment::get_experiment(mae, "missing") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "multi_assay_experiment_common_samples" {
  let rna_exp = @src.Experiment::new(
    name="RNA-seq",
    assay_type="counts",
    data=[[1.0]],
    row_names=["g1"],
    col_names=["s1", "s2", "s3"],
  )
  let chip_exp = @src.Experiment::new(
    name="ChIP-seq",
    assay_type="scores",
    data=[[2.0]],
    row_names=["p1"],
    col_names=["s2", "s3", "s4"],
  )
  let sm = @src.SampleMap::new(assays=[], primary_ids=[], colnames=[])
  let mae = @src.MultiAssayExperiment::new(
    experiments=[rna_exp, chip_exp],
    sample_map=sm,
    col_data=[],
  )
  let common = @src.MultiAssayExperiment::common_samples(mae)
  // s2 and s3 are common
  assert_eq(common.length(), 2)
  assert_true(common.contains("s2"))
  assert_true(common.contains("s3"))
}

///|
test "multi_assay_experiment_subset_by_experiment" {
  let rna_exp = @src.Experiment::new(
    name="RNA-seq",
    assay_type="counts",
    data=[[1.0]],
    row_names=["g1"],
    col_names=["s1"],
  )
  let chip_exp = @src.Experiment::new(
    name="ChIP-seq",
    assay_type="scores",
    data=[[2.0]],
    row_names=["p1"],
    col_names=["s1"],
  )
  let snv_exp = @src.Experiment::new(
    name="SNV",
    assay_type="variants",
    data=[[3.0]],
    row_names=["v1"],
    col_names=["s1"],
  )
  let sm = @src.SampleMap::new(assays=[], primary_ids=[], colnames=[])
  let mae = @src.MultiAssayExperiment::new(
    experiments=[rna_exp, chip_exp, snv_exp],
    sample_map=sm,
    col_data=[],
  )
  let subset = @src.MultiAssayExperiment::subset_by_experiment(mae, [
    "RNA-seq", "SNV",
  ])
  assert_eq(@src.MultiAssayExperiment::num_experiments(subset), 2)
}

///|
test "multi_assay_experiment_subset_by_sample" {
  let rna_exp = @src.Experiment::new(
    name="RNA-seq",
    assay_type="counts",
    data=[[1.0, 2.0, 3.0]],
    row_names=["g1"],
    col_names=["s1", "s2", "s3"],
  )
  let sm = @src.SampleMap::new(
    assays=["RNA-seq", "RNA-seq", "RNA-seq"],
    primary_ids=["p1", "p2", "p3"],
    colnames=["s1", "s2", "s3"],
  )
  let mae = @src.MultiAssayExperiment::new(
    experiments=[rna_exp],
    sample_map=sm,
    col_data=[],
  )
  let subset = @src.MultiAssayExperiment::subset_by_sample(mae, ["p1", "p3"])
  match @src.MultiAssayExperiment::get_experiment(subset, "RNA-seq") {
    Some(exp) => {
      assert_eq(@src.Experiment::ncol(exp), 2)
      assert_eq(@src.Experiment::nrow(exp), 1)
    }
    None => assert_true(false)
  }
}

///|
test "multi_assay_experiment_assay_rownames_colnames" {
  let rna_exp = @src.Experiment::new(
    name="RNA-seq",
    assay_type="counts",
    data=[[1.0, 2.0], [3.0, 4.0]],
    row_names=["geneA", "geneB"],
    col_names=["sample1", "sample2"],
  )
  let sm = @src.SampleMap::new(assays=[], primary_ids=[], colnames=[])
  let mae = @src.MultiAssayExperiment::new(
    experiments=[rna_exp],
    sample_map=sm,
    col_data=[],
  )
  match @src.MultiAssayExperiment::assay(mae, "RNA-seq") {
    Some(d) => assert_eq(d.length(), 2)
    None => assert_true(false)
  }
  match @src.MultiAssayExperiment::rownames(mae, "RNA-seq") {
    Some(r) => {
      assert_eq(r.length(), 2)
      assert_true(r.contains("geneA"))
    }
    None => assert_true(false)
  }
  match @src.MultiAssayExperiment::colnames(mae, "RNA-seq") {
    Some(c) => {
      assert_eq(c.length(), 2)
      assert_true(c.contains("sample1"))
    }
    None => assert_true(false)
  }
}

///|
test "multi_assay_experiment_dims" {
  let rna_exp = @src.Experiment::new(
    name="RNA-seq",
    assay_type="counts",
    data=[[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
    row_names=["g1", "g2", "g3"],
    col_names=["s1", "s2"],
  )
  let chip_exp = @src.Experiment::new(
    name="ChIP-seq",
    assay_type="scores",
    data=[[2.0]],
    row_names=["p1"],
    col_names=["s1"],
  )
  let sm = @src.SampleMap::new(assays=[], primary_ids=[], colnames=[])
  let mae = @src.MultiAssayExperiment::new(
    experiments=[rna_exp, chip_exp],
    sample_map=sm,
    col_data=[],
  )
  let dims = @src.MultiAssayExperiment::dims(mae)
  assert_eq(dims.length(), 2)
  // Find RNA-seq dims
  let mut rna_nrow = 0
  let mut rna_ncol = 0
  for d in dims {
    if d.0 == "RNA-seq" {
      rna_nrow = d.1
      rna_ncol = d.2
    }
  }
  assert_eq(rna_nrow, 3) // nrow
  assert_eq(rna_ncol, 2) // ncol
}

///|
test "multi_assay_experiment_metadata" {
  let rna_exp = @src.Experiment::new(
    name="RNA-seq",
    assay_type="counts",
    data=[[1.0]],
    row_names=["g1"],
    col_names=["s1"],
  )
  let sm = @src.SampleMap::new(assays=[], primary_ids=[], colnames=[])
  let mae = @src.MultiAssayExperiment::new(
    experiments=[rna_exp],
    sample_map=sm,
    col_data=[],
  )
  let mae2 = @src.MultiAssayExperiment::add_metadata(mae, "study", "TCGA-BRCA")
  match @src.MultiAssayExperiment::get_metadata(mae2, "study") {
    Some(v) => assert_eq(v, "TCGA-BRCA")
    None => assert_true(false)
  }
  match @src.MultiAssayExperiment::get_metadata(mae2, "missing") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "multi_assay_experiment_sample_col_data" {
  let rna_exp = @src.Experiment::new(
    name="RNA-seq",
    assay_type="counts",
    data=[[1.0]],
    row_names=["g1"],
    col_names=["s1"],
  )
  let sm = @src.SampleMap::new(assays=[], primary_ids=[], colnames=[])
  let col_data : Array[(String, Map[String, String])] = [
    ("p1", Map([("tissue", "brain"), ("stage", "IV")])),
    ("p2", Map([("tissue", "liver"), ("stage", "II")])),
  ]
  let mae = @src.MultiAssayExperiment::new(
    experiments=[rna_exp],
    sample_map=sm,
    col_data~,
  )
  match @src.MultiAssayExperiment::sample_col_data(mae, "p1") {
    Some(d) =>
      match d.get("tissue") {
        Some(v) => assert_eq(v, "brain")
        None => assert_true(false)
      }
    None => assert_true(false)
  }
  match @src.MultiAssayExperiment::sample_col_data(mae, "p3") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "multi_assay_experiment_summary" {
  let rna_exp = @src.Experiment::new(
    name="RNA-seq",
    assay_type="counts",
    data=[[1.0]],
    row_names=["g1"],
    col_names=["s1"],
  )
  let sm = @src.SampleMap::new(assays=["RNA-seq"], primary_ids=["p1"], colnames=[
    "s1",
  ])
  let mae = @src.MultiAssayExperiment::new(
    experiments=[rna_exp],
    sample_map=sm,
    col_data=[],
  )
  let s = @src.MultiAssayExperiment::summary(mae)
  assert_true(s.contains("MultiAssayExperiment"))
  assert_true(s.contains("1 experiments"))
}

///|
test "create_multi_omics" {
  let rna_counts = [[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]]
  let rna_genes = ["gene1", "gene2"]
  let rna_samples = ["rna_s1", "rna_s2", "rna_s3"]
  let chip_scores = [[1.0, 2.0]]
  let chip_peaks = ["peak1"]
  let chip_samples = ["chip_s1", "chip_s2"]
  let sample_ids = ["p1", "p2", "p3"]
  let mae = @src.create_multi_omics(
    rna_counts, rna_genes, rna_samples, chip_scores, chip_peaks, chip_samples, sample_ids,
  )
  assert_eq(@src.MultiAssayExperiment::num_experiments(mae), 2)
  let names = @src.MultiAssayExperiment::experiment_names(mae)
  assert_true(names.contains("RNA-seq"))
  assert_true(names.contains("ChIP-seq"))
}
