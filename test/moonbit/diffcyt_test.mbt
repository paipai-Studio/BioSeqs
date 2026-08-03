///|
/// Tests for the diffcyt module (Bioconductor diffcyt).
///
/// Covers cytometry cell construction, FlowSOM clustering, DA/DS testing,
// design matrix, BH-FDR correction, and sample data validation.

// ===========================================================================
// 1. CytometryCell construction
// ===========================================================================

///|
test "dc_cell_creation" {
  let cell = @src.CytometryCell::new(
    "cell_1", "S1", "ctrl", [1.0, 2.0, 3.0]
  )
  assert_eq(cell.cell_id(), "cell_1")
  assert_eq(cell.sample_id(), "S1")
  assert_eq(cell.condition(), "ctrl")
  assert_eq(cell.n_markers(), 3)
  assert_eq(cell.marker_values()[0], 1.0)
  assert_eq(cell.cluster_id(), -1)
}

///|
test "dc_cell_set_cluster" {
  let cell = @src.CytometryCell::new("c1", "S1", "ctrl", [1.0, 2.0])
  cell.set_cluster_id(5)
  assert_eq(cell.cluster_id(), 5)
}

// ===========================================================================
// 2. Sample data
// ===========================================================================

///|
test "dc_sample_data_dimensions" {
  let cells = @src.diffcyt_sample_data()
  // 4 samples * 20 cells = 80 cells
  assert_eq(cells.length(), 80)
  assert_eq(cells[0].n_markers(), 3)
  assert_eq(cells[0].sample_id(), "S1")
  assert_eq(cells[0].condition(), "ctrl")
}

///|
test "dc_unique_samples_and_conditions" {
  let cells = @src.diffcyt_sample_data()
  let samples = @src.diffcyt_unique_samples(cells)
  let conditions = @src.diffcyt_unique_conditions(cells)
  assert_eq(samples.length(), 4)
  assert_eq(conditions.length(), 2)
  assert_true(conditions.contains("ctrl"))
  assert_true(conditions.contains("stim"))
}

// ===========================================================================
// 3. FlowSOM clustering
// ===========================================================================

///|
test "dc_cluster_assigns_all_cells" {
  let cells = @src.diffcyt_sample_data()
  let n_clusters = 5
  let codebooks = @src.diffcyt_cluster_cells(cells, n_clusters, 10)
  assert_eq(codebooks.length(), n_clusters)
  // All cells should have a valid cluster id.
  for cell in cells {
    assert_true(cell.cluster_id() >= 0)
    assert_true(cell.cluster_id() < n_clusters)
  }
}

///|
test "dc_cluster_codebooks_have_correct_dims" {
  let cells = @src.diffcyt_sample_data()
  let n_clusters = 3
  let codebooks = @src.diffcyt_cluster_cells(cells, n_clusters, 5)
  for cb in codebooks {
    assert_eq(cb.length(), 3) // 3 markers
  }
}

///|
test "dc_assign_clusters_new" {
  let cells = @src.diffcyt_sample_data()
  let codebooks = @src.diffcyt_cluster_cells(cells, 4, 10)
  // Re-assign using existing codebooks.
  let new_cells = [@src.CytometryCell::new("nc1", "S1", "ctrl", cells[0].marker_values())]
  @src.diffcyt_assign_clusters(new_cells, codebooks)
  assert_true(new_cells[0].cluster_id() >= 0)
  assert_true(new_cells[0].cluster_id() < 4)
}

// ===========================================================================
// 4. Cluster counts and medians
// ===========================================================================

///|
test "dc_calc_counts" {
  let cells = @src.diffcyt_sample_data()
  let _ = @src.diffcyt_cluster_cells(cells, 4, 10)
  let samples = @src.diffcyt_unique_samples(cells)
  let counts = @src.diffcyt_calc_counts(cells, 4, samples)
  assert_eq(counts.length(), 4) // 4 samples
  // Each sample should have 4 cluster counts.
  for row in counts {
    assert_eq(row.length(), 4)
  }
  // Total counts should sum to 20 per sample.
  for row in counts {
    let mut total = 0
    for c in row {
      total = total + c
    }
    assert_eq(total, 20)
  }
}

///|
test "dc_calc_medians" {
  let cells = @src.diffcyt_sample_data()
  let _ = @src.diffcyt_cluster_cells(cells, 3, 10)
  let samples = @src.diffcyt_unique_samples(cells)
  let medians = @src.diffcyt_calc_medians_by_cluster_marker(cells, 3, 3, samples)
  assert_eq(medians.length(), 4) // 4 samples
  // Each row has 3 clusters * 3 markers = 9 entries.
  for row in medians {
    assert_eq(row.length(), 9)
  }
}

// ===========================================================================
// 5. Design matrix
// ===========================================================================

///|
test "dc_design_matrix" {
  let samples = ["S1", "S2", "S3", "S4"]
  let conditions = ["ctrl", "ctrl", "stim", "stim"]
  let design = @src.diffcyt_create_design_matrix(samples, conditions)
  assert_eq(design.length(), 4)
  // 2 unique conditions => 2 columns (intercept + dummy).
  assert_eq(design[0].length(), 2)
  // Intercept is always 1.
  assert_eq(design[0][0], 1.0)
  // ctrl samples: dummy = 0 (reference), stim samples: dummy = 1.
  assert_eq(design[0][1], 0.0) // ctrl
  assert_eq(design[2][1], 1.0) // stim
}

// ===========================================================================
// 6. DA testing
// ===========================================================================

///|
test "dc_testDA_returns_results" {
  let cells = @src.diffcyt_sample_data()
  let _ = @src.diffcyt_cluster_cells(cells, 5, 10)
  let samples = @src.diffcyt_unique_samples(cells)
  let counts = @src.diffcyt_calc_counts(cells, 5, samples)
  let conditions : Array[String] = []
  let cond_map : Map[String, String] = Map::new()
  for cell in cells {
    cond_map.set(cell.sample_id, cell.condition)
  }
  for sid in samples {
    conditions.push(cond_map.get(sid).unwrap_or(""))
  }
  let design = @src.diffcyt_create_design_matrix(samples, conditions)
  let contrast_col = design[0].length() - 1
  let da_results = @src.diffcyt_testDA(counts, design, contrast_col, 5)
  assert_eq(da_results.length(), 5)
  // Check that p-values are in [0, 1].
  for r in da_results {
    assert_true(r.p_value() >= 0.0)
    assert_true(r.p_value() <= 1.0)
    assert_true(r.fdr() >= 0.0)
    assert_true(r.fdr() <= 1.0)
  }
}

///|
test "dc_testDA_fdr_monotonic" {
  let cells = @src.diffcyt_sample_data()
  let _ = @src.diffcyt_cluster_cells(cells, 4, 10)
  let samples = @src.diffcyt_unique_samples(cells)
  let counts = @src.diffcyt_calc_counts(cells, 4, samples)
  let conditions : Array[String] = []
  let cond_map : Map[String, String] = Map::new()
  for cell in cells {
    cond_map.set(cell.sample_id, cell.condition)
  }
  for sid in samples {
    conditions.push(cond_map.get(sid).unwrap_or(""))
  }
  let design = @src.diffcyt_create_design_matrix(samples, conditions)
  let da_results = @src.diffcyt_testDA(counts, design, 1, 4)
  // FDR should be >= p_value (BH correction inflates).
  for r in da_results {
    assert_true(r.fdr() >= r.p_value() - 0.0001)
  }
}

// ===========================================================================
// 7. DS testing
// ===========================================================================

///|
test "dc_testDS_returns_results" {
  let cells = @src.diffcyt_sample_data()
  let _ = @src.diffcyt_cluster_cells(cells, 3, 10)
  let samples = @src.diffcyt_unique_samples(cells)
  let conditions : Array[String] = []
  let cond_map : Map[String, String] = Map::new()
  for cell in cells {
    cond_map.set(cell.sample_id, cell.condition)
  }
  for sid in samples {
    conditions.push(cond_map.get(sid).unwrap_or(""))
  }
  let marker_names = ["marker_0", "marker_1", "marker_2"]
  let ds_results = @src.diffcyt_testDS(cells, samples, conditions, 3, 3, marker_names)
  // 3 clusters * 3 markers = 9 results (if 2 conditions).
  assert_true(ds_results.length() > 0)
  for r in ds_results {
    assert_true(r.p_value() >= 0.0)
    assert_true(r.p_value() <= 1.0)
    assert_true(r.fdr() >= 0.0)
    assert_true(r.fdr() <= 1.0)
  }
}

///|
test "dc_testDS_marker_names" {
  let cells = @src.diffcyt_sample_data()
  let _ = @src.diffcyt_cluster_cells(cells, 2, 10)
  let samples = @src.diffcyt_unique_samples(cells)
  let conditions : Array[String] = []
  let cond_map : Map[String, String] = Map::new()
  for cell in cells {
    cond_map.set(cell.sample_id, cell.condition)
  }
  for sid in samples {
    conditions.push(cond_map.get(sid).unwrap_or(""))
  }
  let marker_names = ["CD4", "CD8", "CD3"]
  let ds_results = @src.diffcyt_testDS(cells, samples, conditions, 2, 3, marker_names)
  for r in ds_results {
    assert_true(marker_names.contains(r.marker_name()))
  }
}

// ===========================================================================
// 8. Full pipeline
// ===========================================================================

///|
test "dc_run_full_pipeline" {
  let cells = @src.diffcyt_sample_data()
  let summary = @src.diffcyt_run(cells, 5, 0.05)
  assert_eq(summary.n_cells, 80)
  assert_eq(summary.n_clusters, 5)
  assert_eq(summary.n_markers, 3)
  assert_eq(summary.n_samples, 4)
  assert_eq(summary.n_conditions, 2)
  assert_eq(summary.fdr_threshold, 0.05)
  assert_true(summary.da_results.length() == 5)
  assert_true(summary.ds_results.length() > 0)
  // n_sig_da and n_sig_ds should be non-negative.
  assert_true(summary.n_sig_da() >= 0)
  assert_true(summary.n_sig_ds() >= 0)
}

// ===========================================================================
// 9. Top results extraction
// ===========================================================================

///|
test "dc_top_clusters_da" {
  let cells = @src.diffcyt_sample_data()
  let _ = @src.diffcyt_cluster_cells(cells, 5, 10)
  let samples = @src.diffcyt_unique_samples(cells)
  let counts = @src.diffcyt_calc_counts(cells, 5, samples)
  let conditions : Array[String] = []
  let cond_map : Map[String, String] = Map::new()
  for cell in cells {
    cond_map.set(cell.sample_id, cell.condition)
  }
  for sid in samples {
    conditions.push(cond_map.get(sid).unwrap_or(""))
  }
  let design = @src.diffcyt_create_design_matrix(samples, conditions)
  let da_results = @src.diffcyt_testDA(counts, design, 1, 5)
  let top3 = @src.diffcyt_top_clusters_da(da_results, 3)
  assert_eq(top3.length(), 3)
  // Should be sorted by p-value ascending.
  assert_true(top3[0].p_value() <= top3[1].p_value())
  assert_true(top3[1].p_value() <= top3[2].p_value())
}

///|
test "dc_top_table_ds" {
  let cells = @src.diffcyt_sample_data()
  let _ = @src.diffcyt_cluster_cells(cells, 3, 10)
  let samples = @src.diffcyt_unique_samples(cells)
  let conditions : Array[String] = []
  let cond_map : Map[String, String] = Map::new()
  for cell in cells {
    cond_map.set(cell.sample_id, cell.condition)
  }
  for sid in samples {
    conditions.push(cond_map.get(sid).unwrap_or(""))
  }
  let ds_results = @src.diffcyt_testDS(
    cells, samples, conditions, 3, 3, ["m0", "m1", "m2"]
  )
  let top5 = @src.diffcyt_top_table_ds(ds_results, 5)
  assert_true(top5.length() <= 5)
  if top5.length() >= 2 {
    assert_true(top5[0].p_value() <= top5[1].p_value())
  }
}

// ===========================================================================
// 10. Mark significant
// ===========================================================================

///|
test "dc_mark_significant" {
  let cells = @src.diffcyt_sample_data()
  let summary = @src.diffcyt_run(cells, 4, 0.5)
  // With high threshold, some results should be significant.
  let mut any_sig = false
  for r in summary.da_results {
    if r.significant() {
      any_sig = true
    }
  }
  // With FDR 0.5, at least some should pass.
  assert_true(any_sig || summary.da_results.length() == 0)
}

// ===========================================================================
// 11. Edge cases
// ===========================================================================

///|
test "dc_empty_cells" {
  let cells : Array[@src.CytometryCell] = []
  let codebooks = @src.diffcyt_cluster_cells(cells, 3, 10)
  assert_eq(codebooks.length(), 0)
}

///|
test "dc_zero_clusters" {
  let cells = @src.diffcyt_sample_data()
  let codebooks = @src.diffcyt_cluster_cells(cells, 0, 10)
  assert_eq(codebooks.length(), 0)
}

///|
test "dc_single_sample" {
  let cells = [@src.CytometryCell::new("c1", "S1", "ctrl", [1.0, 2.0, 3.0])]
  let samples = @src.diffcyt_unique_samples(cells)
  assert_eq(samples.length(), 1)
}
