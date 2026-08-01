///|
/// Tests for scmap single-cell cell type annotation.

// ---------------------------------------------------------------------------
// Data structure creation and accessors
// ---------------------------------------------------------------------------

test "scmap_reference_creation" {
  let ref = @src.scmap_sample_reference()
  assert_eq(ref.gene_names().length(), 10)
  assert_eq(ref.cell_names().length(), 9)
  assert_eq(ref.cell_types().length(), 9)
  assert_eq(ref.cell_types()[0], "T_cell")
  assert_eq(ref.cell_types()[3], "B_cell")
  assert_eq(ref.cell_types()[6], "NK_cell")
}

test "scmap_reference_unique_types" {
  let ref = @src.scmap_sample_reference()
  let types = ref.unique_cell_types()
  assert_eq(types.length(), 3)
}

test "scmap_query_creation" {
  let q = @src.scmap_sample_query()
  assert_eq(q.gene_names().length(), 10)
  assert_eq(q.cell_names().length(), 3)
  assert_eq(q.cell_names()[0], "query_T")
  assert_eq(q.cell_names()[1], "query_B")
  assert_eq(q.cell_names()[2], "query_unknown")
}

test "scmap_assignment_creation" {
  let a = @src.ScmapAssignment::new(
    "cell1", "T_cell", 0.85, 0.3, "cluster",
  )
  assert_eq(a.cell_name(), "cell1")
  assert_eq(a.assigned_type(), "T_cell")
  assert_eq(a.best_correlation(), 0.85)
  assert_eq(a.second_best_correlation(), 0.3)
  assert_eq(a.method(), "cluster")
}

// ---------------------------------------------------------------------------
// scmap-cluster method
// ---------------------------------------------------------------------------

test "scmap_cluster_classifies_t_cell" {
  let ref = @src.scmap_sample_reference()
  let query = @src.scmap_sample_query()
  let assignments = @src.scmap_cluster(ref, query, 0.5)
  assert_eq(assignments.length(), 3)
  // query_T should be classified as T_cell
  assert_eq(assignments[0].cell_name(), "query_T")
  assert_eq(assignments[0].assigned_type(), "T_cell")
  assert_eq(assignments[0].method(), "cluster")
}

test "scmap_cluster_classifies_b_cell" {
  let ref = @src.scmap_sample_reference()
  let query = @src.scmap_sample_query()
  let assignments = @src.scmap_cluster(ref, query, 0.5)
  // query_B should be classified as B_cell
  assert_eq(assignments[1].cell_name(), "query_B")
  assert_eq(assignments[1].assigned_type(), "B_cell")
}

test "scmap_cluster_high_threshold_unassigned" {
  let ref = @src.scmap_sample_reference()
  let query = @src.scmap_sample_query()
  // Very high threshold should result in unassigned
  let assignments = @src.scmap_cluster(ref, query, 0.99)
  for a in assignments {
    assert_eq(a.assigned_type(), "unassigned")
  }
}

test "scmap_cluster_correlation_is_valid" {
  let ref = @src.scmap_sample_reference()
  let query = @src.scmap_sample_query()
  let assignments = @src.scmap_cluster(ref, query, 0.5)
  for a in assignments {
    assert_true(a.best_correlation() >= -1.0)
    assert_true(a.best_correlation() <= 1.0)
  }
}

test "scmap_cluster_best_gt_second_best" {
  let ref = @src.scmap_sample_reference()
  let query = @src.scmap_sample_query()
  let assignments = @src.scmap_cluster(ref, query, 0.5)
  for a in assignments {
    assert_true(a.best_correlation() >= a.second_best_correlation())
  }
}

// ---------------------------------------------------------------------------
// scmap-cell method
// ---------------------------------------------------------------------------

test "scmap_cell_classifies_t_cell" {
  let ref = @src.scmap_sample_reference()
  let query = @src.scmap_sample_query()
  let assignments = @src.scmap_cell(ref, query, 3, 0.5)
  assert_eq(assignments.length(), 3)
  assert_eq(assignments[0].cell_name(), "query_T")
  assert_eq(assignments[0].method(), "cell")
}

test "scmap_cell_classifies_b_cell" {
  let ref = @src.scmap_sample_reference()
  let query = @src.scmap_sample_query()
  let assignments = @src.scmap_cell(ref, query, 3, 0.5)
  assert_eq(assignments[1].cell_name(), "query_B")
  assert_eq(assignments[1].assigned_type(), "B_cell")
}

test "scmap_cell_k_neighbours_parameter" {
  let ref = @src.scmap_sample_reference()
  let query = @src.scmap_sample_query()
  // k=1 should still work
  let assignments = @src.scmap_cell(ref, query, 1, 0.5)
  assert_eq(assignments.length(), 3)
}

test "scmap_cell_high_threshold_unassigned" {
  let ref = @src.scmap_sample_reference()
  let query = @src.scmap_sample_query()
  // Very high threshold fraction should make it harder to assign
  // query_unknown has low expression everywhere, so with k=9 (all ref cells)
  // the vote is split across types
  let assignments = @src.scmap_cell(ref, query, 9, 0.99)
  // query_unknown should be unassigned (no single type gets >99% of 9 votes = 9)
  // With 9 neighbours, need 9 votes for 1.0 fraction. Since ref has 3 of each type,
  // max votes = 3, fraction = 3/9 = 0.33 < 0.99
  assert_eq(assignments[2].assigned_type(), "unassigned")
}

// ---------------------------------------------------------------------------
// Summary utilities
// ---------------------------------------------------------------------------

test "scmap_summary_counts" {
  let ref = @src.scmap_sample_reference()
  let query = @src.scmap_sample_query()
  let assignments = @src.scmap_cluster(ref, query, 0.5)
  let summary = @src.scmap_summary(assignments)
  // Should have at least one type in summary
  assert_true(summary.length() >= 1)
  // Total should equal number of assignments
  let mut total = 0
  for entry in summary {
    let (_, count) = entry
    total = total + count
  }
  assert_eq(total, assignments.length())
}

test "scmap_assignment_to_string" {
  let a = @src.ScmapAssignment::new(
    "cell1", "T_cell", 0.85, 0.3, "cluster",
  )
  let s = a.to_string()
  assert_true(s.contains("cell1"))
  assert_true(s.contains("T_cell"))
  assert_true(s.contains("cluster"))
}

test "scmap_assignments_to_string" {
  let ref = @src.scmap_sample_reference()
  let query = @src.scmap_sample_query()
  let assignments = @src.scmap_cluster(ref, query, 0.5)
  let s = @src.scmap_assignments_to_string(assignments)
  assert_true(s.contains("scmap Assignments"))
  assert_true(s.contains("3 cells"))
}

// ---------------------------------------------------------------------------
// Sample data verification
// ---------------------------------------------------------------------------

test "scmap_sample_reference_expression_structure" {
  let ref = @src.scmap_sample_reference()
  // 10 genes × 9 cells
  assert_eq(ref.expression().length(), 10)
  assert_eq(ref.expression()[0].length(), 9)
  // T cell marker genes (GENE1-3) should be high in T cells (cols 0-2)
  assert_true(ref.expression()[0][0] > ref.expression()[0][3])
  // B cell marker genes (GENE4-6) should be high in B cells (cols 3-5)
  assert_true(ref.expression()[3][3] > ref.expression()[3][0])
}

test "scmap_sample_query_expression_structure" {
  let q = @src.scmap_sample_query()
  // 10 genes × 3 cells
  assert_eq(q.expression().length(), 10)
  assert_eq(q.expression()[0].length(), 3)
  // GENE1 should be high in query_T (col 0), low in query_B (col 1)
  assert_true(q.expression()[0][0] > q.expression()[0][1])
}

// ---------------------------------------------------------------------------
// Edge cases
// ---------------------------------------------------------------------------

test "scmap_cluster_single_query_cell" {
  let ref = @src.scmap_sample_reference()
  // Create a single-cell query using the T cell reference profile
  let single_expr : Array[Array[Double]] = []
  for g in 0..<10 {
    single_expr.push([ref.expression()[g][0]])
  }
  let q = @src.ScmapQuery::new(ref.gene_names(), ["single"], single_expr)
  let assignments = @src.scmap_cluster(ref, q, 0.5)
  assert_eq(assignments.length(), 1)
}

test "scmap_cell_single_query_cell" {
  let ref = @src.scmap_sample_reference()
  let q = @src.scmap_sample_query()
  let assignments = @src.scmap_cell(ref, q, 5, 0.5)
  assert_eq(assignments.length(), 3)
}
