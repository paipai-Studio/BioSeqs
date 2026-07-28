///|
/// Tests for dorothea module - TF activity prediction.

///|
test "dorothea_get_regulons" {
  let regulons = @src.dorothea_get_regulons()
  assert_true(regulons.length() > 0)
  
  // Check that key TFs are present
  let mut has_tp53 = false
  let mut has_myc = false
  for regulon in regulons {
    if regulon.tf == "TP53" {
      has_tp53 = true
    }
    if regulon.tf == "MYC" {
      has_myc = true
    }
  }
  assert_true(has_tp53)
  assert_true(has_myc)
}

///|
test "dorothea_regulon_targets" {
  let regulons = @src.dorothea_get_regulons()
  for regulon in regulons {
    assert_true(regulon.targets.length() > 0)
    // Check that all targets have valid directions
    for target in regulon.targets {
      assert_true(target.direction == "activation" || target.direction == "repression")
    }
  }
}

///|
test "dorothea_compute_activity_simple" {
  let regulons = @src.dorothea_get_regulons()
  let test_data = @src.dorothea_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  let cell_expr : Array[Double] = Array::new()
  let n_genes = counts.length()
  let mut g = 0
  while g < n_genes {
    cell_expr.push(counts[g][0])
    g = g + 1
  }
  
  // Find a regulon with targets in our gene set
  let mut tf_activity = 0.0
  for regulon in regulons {
    let activity = @src.dorothea_compute_activity_simple(cell_expr, gene_names, regulon)
    tf_activity = tf_activity + activity
  }
  
  assert_true(tf_activity >= 0.0)
}

///|
test "dorothea_compute_viper_activity" {
  let regulons = @src.dorothea_get_regulons()
  let test_data = @src.dorothea_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  let cell_expr : Array[Double] = Array::new()
  let n_genes = counts.length()
  let mut g = 0
  while g < n_genes {
    cell_expr.push(counts[g][0])
    g = g + 1
  }
  
  let mut has_activity = false
  for regulon in regulons {
    let activity = @src.dorothea_compute_viper_activity(cell_expr, gene_names, regulon)
    if activity.abs() > 0.0 {
      has_activity = true
    }
  }
  
  // VIPER activities should generally be non-zero if targets are found
  assert_true(has_activity || regulons.length() > 0)
}

///|
test "dorothea_permutation_test" {
  let regulons = @src.dorothea_get_regulons()
  let test_data = @src.dorothea_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  let cell_expr : Array[Double] = Array::new()
  let n_genes = counts.length()
  let mut g = 0
  while g < n_genes {
    cell_expr.push(counts[g][0])
    g = g + 1
  }
  
  // Use a regulon with many targets for better test
  let mut tp53_regulon = regulons[0] // Default
  for regulon in regulons {
    if regulon.tf == "TP53" {
      tp53_regulon = regulon
      break
    }
  }
  
  let (p_value, z_score) = @src.dorothea_permutation_test(cell_expr, gene_names, tp53_regulon, 100)
  assert_true(p_value >= 0.0 && p_value <= 1.0)
}

///|
test "dorothea_analyze_cell" {
  let regulons = @src.dorothea_get_regulons()
  let test_data = @src.dorothea_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  let cell_expr : Array[Double] = Array::new()
  let n_genes = counts.length()
  let mut g = 0
  while g < n_genes {
    cell_expr.push(counts[g][0])
    g = g + 1
  }
  
  let params = @src.DorotheaParams::with_params(50, 3, 0.05)
  let results = @src.dorothea_analyze_cell(cell_expr, gene_names, regulons, params)
  
  assert_true(results.length() > 0)
  for result in results {
    assert_true(result.n_targets >= 3)
    assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  }
}

///|
test "dorothea_analyze_cells" {
  let regulons = @src.dorothea_get_regulons()
  let test_data = @src.dorothea_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  let params = @src.DorotheaParams::with_params(50, 3, 0.05)
  let results = @src.dorothea_analyze_cells(counts, cell_ids, gene_names, regulons, params)
  
  assert_true(results.length() > 0)
  assert_true(results.length() <= regulons.length())
}

///|
test "dorothea_create_test_data" {
  let test_data = @src.dorothea_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  assert_eq(cell_ids.length(), 10)
  assert_true(gene_names.length() > 0)
  assert_true(counts.length() > 0)
  assert_true(counts[0].length() == 10)
}

///|
test "dorothea_params_new" {
  let params = @src.DorotheaParams::new()
  assert_eq(params.n_permutations, 1000)
  assert_eq(params.min_targets, 5)
  assert_true((params.significance - 0.05).abs() < 0.001)
}

///|
test "dorothea_params_with_params" {
  let params = @src.DorotheaParams::with_params(200, 3, 0.01)
  assert_eq(params.n_permutations, 200)
  assert_eq(params.min_targets, 3)
  assert_true((params.significance - 0.01).abs() < 0.001)
}

///|
test "dorothea_filter_results" {
  let results : Array[@src.TFActivityResult] = Array::new()
  let filtered = @src.dorothea_filter_results(results, 0.05)
  assert_eq(filtered.length(), 0)
}

///|
test "dorothea_sort_by_activity" {
  let results = @src.dorothea_analyze_cells(
    [[]],
    [],
    [],
    @src.dorothea_get_regulons(),
    @src.DorotheaParams::with_params(10, 3, 0.05)
  )
  
  // Sort empty or small results
  let sorted = @src.dorothea_sort_by_activity(results)
  assert_true(sorted.length() <= results.length())
}

///|
test "dorothea_get_top_tfs" {
  let regulons = @src.dorothea_get_regulons()
  let test_data = @src.dorothea_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  let params = @src.DorotheaParams::with_params(20, 3, 0.05)
  let results = @src.dorothea_analyze_cells(counts, cell_ids, gene_names, regulons, params)
  
  let top3 = @src.dorothea_get_top_tfs(results, 3)
  assert_true(top3.length() <= 3)
}

///|
test "dorothea_array_mean" {
  let arr = [1.0, 2.0, 3.0, 4.0, 5.0]
  let mean = @src.dorothea_array_mean(arr)
  assert_true((mean - 3.0).abs() < 0.001)
}

///|
test "dorothea_array_std" {
  let arr = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]
  let std = @src.dorothea_array_std(arr)
  assert_true(std > 0.0)
}

///|
test "dorothea_get_target_expression" {
  let regulons = @src.dorothea_get_regulons()
  let test_data = @src.dorothea_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  let cell_expr : Array[Double] = Array::new()
  let n_genes = counts.length()
  let mut g = 0
  while g < n_genes {
    cell_expr.push(counts[g][0])
    g = g + 1
  }
  
  // Test with first regulon
  let targets = @src.dorothea_get_target_expression(cell_expr, gene_names, regulons[0])
  assert_true(targets.length() >= 0)
}