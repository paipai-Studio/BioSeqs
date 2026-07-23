///|
test "tradeseq_create_trajectory_point" {
  let point = @src.TrajectoryPoint::new(
    "cell_001",
    0.5,
    [1.0, 2.0, 3.0],
    "control"
  )
  
  assert_eq(point.cell_id, "cell_001")
  assert_true((point.pseudotime - 0.5).abs() < 1.0e-10)
  assert_eq(point.expression.length(), 3)
  assert_eq(point.condition, "control")
}

///|
test "tradeseq_create_gene_expression" {
  let data = @src.GeneExpressionData::new(
    "GeneA",
    [1.0, 2.0, 3.0, 4.0],
    [0.1, 0.3, 0.5, 0.7],
    ["control", "control", "treatment", "treatment"]
  )
  
  assert_eq(data.gene_id, "GeneA")
  assert_eq(data.expression.length(), 4)
  assert_eq(data.pseudotime.length(), 4)
  assert_eq(data.condition.length(), 4)
}

///|
test "tradeseq_trajectory_data" {
  let mut data = @src.TrajectoryData::new()
  
  assert_eq(data.n_points, 0)
  assert_eq(data.n_genes, 0)
  
  data = data.add_gene("GeneA")
  assert_eq(data.n_genes, 1)
  
  let point = @src.TrajectoryPoint::new(
    "cell_001",
    0.5,
    [1.0],
    "control"
  )
  
  data = data.add_point(point)
  assert_eq(data.n_points, 1)
}

///|
test "tradeseq_fit_gam" {
  let gene_data = @src.GeneExpressionData::new(
    "GeneA",
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
    ["control", "control", "control", "treatment", "treatment", "treatment"]
  )
  
  let gam = @src.fit_gam(gene_data, 4)
  
  assert_eq(gam.gene_id, "GeneA")
  assert_true(gam.fitted_values.length() > 0)
}

///|
test "tradeseq_trade_test" {
  let gene_data = @src.GeneExpressionData::new(
    "GeneA",
    [1.0, 1.2, 1.5, 2.0, 2.5, 3.0],
    [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
    ["control", "control", "control", "treatment", "treatment", "treatment"]
  )
  
  let result = @src.trade_test_condition_effect(gene_data, 0.5, "control", "treatment")
  
  assert_eq(result.gene_id, "GeneA")
  assert_true(result.p_value >= 0.0)
  // Just check p_value is non-negative
}

///|
test "tradeseq_run_analysis" {
  let trajectory_data = @src.create_example_trajectory_data()
  
  assert_true(trajectory_data.n_points > 0)
  assert_true(trajectory_data.n_genes > 0)
  
  let result = @src.run_tradeseq_analysis(
    trajectory_data,
    trajectory_data.genes,
    trajectory_data.conditions,
    4,
    0.05
  )
  
  assert_eq(result.n_genes, 3)
  assert_true(result.n_significant >= 0)
  assert_true(result.n_significant <= result.n_genes)
}

///|
test "tradeseq_calculate_smooth" {
  let trajectory_data = @src.create_example_trajectory_data()
  
  let smooth = @src.calculate_gene_smooth(trajectory_data, "GeneA", 4, 20)
  
  assert_eq(smooth.gene_id, "GeneA")
  assert_eq(smooth.pseudotime.length(), 20)
  assert_eq(smooth.fitted.length(), 20)
}

///|
test "tradeseq_summary" {
  let trajectory_data = @src.create_example_trajectory_data()
  let result = @src.run_tradeseq_analysis(
    trajectory_data,
    trajectory_data.genes,
    trajectory_data.conditions,
    4,
    0.05
  )
  
  let summary = @src.tradeseq_summary(result)
  
  assert_true(summary.contains("tradeSeq Analysis Summary"))
  assert_true(summary.contains("Total genes tested:"))
}

///|
test "tradeseq_create_example" {
  let data = @src.create_example_trajectory_data()
  
  assert_true(data.n_points == 50)
  assert_true(data.n_genes == 3)
  assert_true(data.conditions.length() >= 2)
}

///|
test "tradeseq_condition_effect_zero" {
  let gene_data = @src.GeneExpressionData::new(
    "GeneA",
    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
    ["control", "control", "control", "treatment", "treatment", "treatment"]
  )
  
  let result = @src.trade_test_condition_effect(gene_data, 0.5, "control", "treatment")
  
  // No difference should give high p-value
  assert_true(result.p_value > 0.05)
}

///|
test "tradeseq_trajectory_empty" {
  let mut data = @src.TrajectoryData::new()
  data = data.add_gene("GeneA")
  
  let result = @src.run_tradeseq_analysis(
    data,
    ["GeneA"],
    ["control"],
    4,
    0.05
  )
  
  assert_eq(result.n_genes, 1)
}

///|
test "tradeseq_multiple_smooth" {
  let trajectory_data = @src.create_example_trajectory_data()
  
  let smooth_a = @src.calculate_gene_smooth(trajectory_data, "GeneA", 4, 15)
  let smooth_b = @src.calculate_gene_smooth(trajectory_data, "GeneB", 4, 15)
  
  assert_eq(smooth_a.gene_id, "GeneA")
  assert_eq(smooth_b.gene_id, "GeneB")
  assert_eq(smooth_a.fitted.length(), 15)
  assert_eq(smooth_b.fitted.length(), 15)
}
