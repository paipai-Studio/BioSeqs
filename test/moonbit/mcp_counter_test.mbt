///|
/// Tests for Bioconductor MCPcounter module - marker-based cell population
/// counting via geometric mean of marker gene expression.

// ---------------------------------------------------------------------------
// Data structure creation and accessors
// ---------------------------------------------------------------------------

test "mcp_cell_population_struct" {
  let pop = @src.McpCellPopulation::new("T cells", ["CD3D", "CD3E", "CD2"])
  assert_eq(pop.name(), "T cells")
  assert_eq(pop.markers().length(), 3)
}

test "mcp_result_accessors" {
  // Use mcp_run to get a result with known values
  let result = @src.mcp_run(
    ["G1", "G2"],
    ["S1", "S2"],
    [[4.0, 9.0], [4.0, 9.0]],
    populations=[@src.McpCellPopulation::new("PopA", ["G1", "G2"])],
  )
  assert_eq(result.populations().length(), 1)
  assert_eq(result.sample_names().length(), 2)
  assert_eq(result.scores().length(), 1)
  // Geometric mean of (4, 4) = 4.0 for S1
  assert_true((result.scores()[0][0] - 4.0).abs() < 0.001)
}

test "mcp_result_get_score" {
  // Build a result with known values via mcp_run
  // Geometric mean of (4, 9) = 6.0 for S1, (4, 4) = 4.0 for S2
  let result = @src.mcp_run(
    ["G1", "G2"],
    ["S1", "S2"],
    [[4.0, 4.0], [9.0, 4.0]],
    populations=[@src.McpCellPopulation::new("PopA", ["G1", "G2"])],
  )
  let s1 = result.get_score("PopA", "S1")
  assert_true((s1 - 6.0).abs() < 0.001)
  let s2 = result.get_score("PopA", "S2")
  assert_true((s2 - 4.0).abs() < 0.001)
  // Unknown population / sample returns 0.0
  assert_true(result.get_score("Unknown", "S1").abs() < 1.0e-12)
  assert_true(result.get_score("PopA", "Unknown").abs() < 1.0e-12)
}

test "mcp_result_get_population_scores" {
  // Build a result with two populations
  let result = @src.mcp_run(
    ["G1", "G2"],
    ["S1", "S2"],
    [[4.0, 9.0], [4.0, 9.0]],
    populations=[
      @src.McpCellPopulation::new("PopA", ["G1"]),
      @src.McpCellPopulation::new("PopB", ["G2"]),
    ],
  )
  let scores_a = result.get_population_scores("PopA")
  assert_eq(scores_a.length(), 2)
  assert_true((scores_a[0] - 4.0).abs() < 1.0e-9)
  assert_true((scores_a[1] - 9.0).abs() < 1.0e-9)
  let scores_b = result.get_population_scores("PopB")
  assert_eq(scores_b.length(), 2)
  assert_true((scores_b[0] - 4.0).abs() < 1.0e-9)
  // Unknown population returns empty array
  let empty_scores = result.get_population_scores("Unknown")
  assert_eq(empty_scores.length(), 0)
}

// ---------------------------------------------------------------------------
// Built-in populations
// ---------------------------------------------------------------------------

test "mcp_default_populations_count" {
  let pops = @src.mcp_default_populations()
  assert_eq(pops.length(), 10)
  // Each population should have at least one marker
  for pop in pops {
    assert_true(pop.markers().length() > 0)
  }
}

test "mcp_population_names" {
  let names = @src.mcp_population_names()
  assert_eq(names.length(), 10)
  assert_eq(names[0], "T cells")
  assert_eq(names[9], "Fibroblasts")
}

test "mcp_default_populations_marker_genes_present" {
  let pops = @src.mcp_default_populations()
  // T cells should contain CD3D
  let t_cells = pops[0]
  let mut has_cd3d = false
  for m in t_cells.markers() {
    if m == "CD3D" { has_cd3d = true; break }
  }
  assert_true(has_cd3d)
  // B lineage should contain CD19
  let b_cells = pops[3]
  let mut has_cd19 = false
  for m in b_cells.markers() {
    if m == "CD19" { has_cd19 = true; break }
  }
  assert_true(has_cd19)
}

// ---------------------------------------------------------------------------
// Geometric mean computation via mcp_run
// ---------------------------------------------------------------------------

test "mcp_run_geometric_mean_correctness" {
  // Single population with two markers, one sample
  // Geometric mean of (4.0, 9.0) = sqrt(4 * 9) = 6.0
  let pop = @src.McpCellPopulation::new("Test", ["G1", "G2"])
  let result = @src.mcp_run(
    ["G1", "G2"],
    ["S1"],
    [[4.0], [9.0]],
    populations=[pop],
  )
  assert_eq(result.populations().length(), 1)
  assert_eq(result.sample_names().length(), 1)
  let s = result.get_score("Test", "S1")
  assert_true((s - 6.0).abs() < 0.001)
}

test "mcp_run_handles_zero_expression" {
  // If all marker expressions are zero, score should be 0
  let pop = @src.McpCellPopulation::new("Test", ["G1", "G2"])
  let result = @src.mcp_run(
    ["G1", "G2"],
    ["S1"],
    [[0.0], [0.0]],
    populations=[pop],
  )
  assert_true(result.get_score("Test", "S1").abs() < 1.0e-12)
}

test "mcp_run_handles_missing_markers" {
  // If a marker is not in gene_names, it is skipped
  let pop = @src.McpCellPopulation::new("Test", ["G1", "MISSING"])
  let result = @src.mcp_run(
    ["G1"],
    ["S1"],
    [[4.0]],
    populations=[pop],
  )
  // Geometric mean of just [4.0] = 4.0
  let s = result.get_score("Test", "S1")
  assert_true((s - 4.0).abs() < 1.0e-9)
}

test "mcp_run_all_markers_missing" {
  // If no markers are found, score is 0
  let pop = @src.McpCellPopulation::new("Test", ["X1", "X2"])
  let result = @src.mcp_run(
    ["G1", "G2"],
    ["S1"],
    [[4.0], [9.0]],
    populations=[pop],
  )
  assert_true(result.get_score("Test", "S1").abs() < 1.0e-12)
}

test "mcp_run_multiple_samples" {
  // Two samples, one population
  let pop = @src.McpCellPopulation::new("Test", ["G1"])
  let result = @src.mcp_run(
    ["G1"],
    ["S1", "S2"],
    [[2.0, 8.0]],
    populations=[pop],
  )
  assert_eq(result.sample_names().length(), 2)
  assert_true((result.get_score("Test", "S1") - 2.0).abs() < 1.0e-9)
  assert_true((result.get_score("Test", "S2") - 8.0).abs() < 1.0e-9)
}

test "mcp_run_default_populations" {
  // Run with default populations on sample data
  let (genes, samples, matrix) = @src.mcp_sample_data()
  let result = @src.mcp_run(genes, samples, matrix)
  assert_eq(result.populations().length(), 10)
  assert_eq(result.sample_names().length(), 3)
  // All scores should be non-negative
  for row in result.scores() {
    for v in row {
      assert_true(v >= 0.0)
    }
  }
}

// ---------------------------------------------------------------------------
// Sample data correctness
// ---------------------------------------------------------------------------

test "mcp_sample_data_t_cell_rich_tumor_a" {
  // Tumor_A is T-cell rich: T cells (pop 0), CD8+ T cells (pop 1),
  // Cytotoxic lymphocytes (pop 2) should all score higher than B lineage (pop 3)
  let (genes, samples, matrix) = @src.mcp_sample_data()
  let result = @src.mcp_run(genes, samples, matrix)
  let t_cell_a = result.get_score("T cells", "Tumor_A")
  let b_cell_a = result.get_score("B lineage", "Tumor_A")
  assert_true(t_cell_a > b_cell_a)
}

test "mcp_sample_data_b_cell_rich_tumor_b" {
  // Tumor_B is B-cell rich: B lineage should score higher than T cells
  let (genes, samples, matrix) = @src.mcp_sample_data()
  let result = @src.mcp_run(genes, samples, matrix)
  let t_cell_b = result.get_score("T cells", "Tumor_B")
  let b_cell_b = result.get_score("B lineage", "Tumor_B")
  assert_true(b_cell_b > t_cell_b)
}

test "mcp_sample_data_fibroblast_rich_tumor_c" {
  // Tumor_C is fibroblast-rich: Fibroblasts (pop 9) should score higher than T cells
  let (genes, samples, matrix) = @src.mcp_sample_data()
  let result = @src.mcp_run(genes, samples, matrix)
  let t_cell_c = result.get_score("T cells", "Tumor_C")
  let fibro_c = result.get_score("Fibroblasts", "Tumor_C")
  assert_true(fibro_c > t_cell_c)
}

// ---------------------------------------------------------------------------
// to_string formatting
// ---------------------------------------------------------------------------

test "mcp_result_to_string" {
  let (genes, samples, matrix) = @src.mcp_sample_data()
  let result = @src.mcp_run(genes, samples, matrix)
  let s = result.to_string()
  assert_true(s.contains("MCPcounter"))
  assert_true(s.contains("Population"))
  assert_true(s.contains("Tumor_A"))
  assert_true(s.contains("T cells"))
}

// ---------------------------------------------------------------------------
// Empty / edge cases
// ---------------------------------------------------------------------------

test "mcp_run_empty_population_list" {
  let result = @src.mcp_run(
    ["G1"],
    ["S1"],
    [[1.0]],
    populations=[],
  )
  assert_eq(result.populations().length(), 0)
  assert_eq(result.sample_names().length(), 1)
}

test "mcp_run_no_samples" {
  let pop = @src.McpCellPopulation::new("Test", ["G1"])
  let result = @src.mcp_run(
    ["G1"],
    [],
    [],
    populations=[pop],
  )
  assert_eq(result.sample_names().length(), 0)
  assert_eq(result.scores().length(), 1) // 1 population row, but empty
  assert_eq(result.scores()[0].length(), 0)
}

test "mcp_run_geometric_mean_with_negative_values" {
  // Negative expression values are skipped (treated as invalid)
  let pop = @src.McpCellPopulation::new("Test", ["G1", "G2"])
  let result = @src.mcp_run(
    ["G1", "G2"],
    ["S1"],
    [[-1.0], [9.0]],
    populations=[pop],
  )
  // Only G2 contributes -> geometric mean = 9.0
  let s = result.get_score("Test", "S1")
  assert_true((s - 9.0).abs() < 1.0e-9)
}

test "mcp_run_custom_populations" {
  // Custom populations list with non-default markers
  let custom_pop = @src.McpCellPopulation::new("Custom", ["G1", "G2", "G3"])
  let result = @src.mcp_run(
    ["G1", "G2", "G3"],
    ["S1"],
    [[1.0], [4.0], [16.0]],
    populations=[custom_pop],
  )
  // Geometric mean of (1, 4, 16) = (1*4*16)^(1/3) = 64^(1/3) = 4.0
  let s = result.get_score("Custom", "S1")
  assert_true((s - 4.0).abs() < 0.001)
}
