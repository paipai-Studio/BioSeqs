///|
test "aucell_rank_genes" {
  let expression = [5.0, 2.0, 8.0, 1.0, 4.0]
  let ranks = @src.aucell_rank_genes(expression)
  // Should be sorted by expression value (ascending)
  assert_eq(ranks.length(), 5)
  // The gene with value 1.0 (index 3) should be first
  assert_eq(ranks[0], 3)
  // The gene with value 8.0 (index 2) should be last
  assert_eq(ranks[4], 2)
}

///|
test "aucell_rank_genes_sorted" {
  let expression = [1.0, 2.0, 3.0, 4.0, 5.0]
  let ranks = @src.aucell_rank_genes(expression)
  assert_eq(ranks.length(), 5)
  assert_eq(ranks[0], 0)
  assert_eq(ranks[1], 1)
  assert_eq(ranks[2], 2)
  assert_eq(ranks[3], 3)
  assert_eq(ranks[4], 4)
}

///|
test "aucell_compute_auc_basic" {
  // 5 genes, gene set contains genes at indices 0, 2, 4
  // Ranks are assumed to be sequential
  let total_genes = 5
  let gene_set_indices = [0, 2, 4]
  let ranks = [0, 1, 2, 3, 4] // Identity mapping

  let auc = @src.aucell_compute_auc(total_genes, gene_set_indices, ranks)
  assert_eq(auc >= 0.0 && auc <= 1.0, true)
}

///|
test "aucell_compute_auc_all_genes" {
  let total_genes = 3
  let gene_set_indices = [0, 1, 2]
  let ranks = [0, 1, 2]

  let auc = @src.aucell_compute_auc(total_genes, gene_set_indices, ranks)
  // When all genes are in the set, AUC should be reasonable (> 0.0)
  assert_eq(auc > 0.0, true)
}

///|
test "aucell_compute_auc_empty_set" {
  let total_genes = 5
  let gene_set_indices : Array[Int] = []
  let ranks = [0, 1, 2, 3, 4]

  let auc = @src.aucell_compute_auc(total_genes, gene_set_indices, ranks)
  assert_eq(auc, 0.0)
}

///|
test "aucell_gene_set" {
  let genes = ["GeneA", "GeneB", "GeneC"]
  let gs = @src.GeneSet::new("Set1", genes)
  assert_eq(gs.name, "Set1")
  assert_eq(gs.genes.length(), 3)
}

///|
test "aucell_data" {
  let expression = [[5.0, 3.0], [2.0, 8.0]]
  let cell_names = ["Cell1", "Cell2"]
  let gene_names = ["GeneA", "GeneB"]
  let gs1 = @src.GeneSet::new("Set1", ["GeneA"])
  let data = @src.AUCellData::new(expression, cell_names, gene_names, [gs1])
  assert_eq(data.expression.length(), 2)
  assert_eq(data.cell_names.length(), 2)
  assert_eq(data.gene_sets.length(), 1)
}

///|
test "aucell_run_basic" {
  // 2 cells × 3 genes
  let expression = [
    [10.0, 2.0, 5.0],
    [3.0, 8.0, 1.0],
  ]
  let cell_names = ["Cell1", "Cell2"]
  let gene_names = ["GeneA", "GeneB", "GeneC"]

  // Gene set with GeneA and GeneB
  let gs = @src.GeneSet::new("Signature", ["GeneA", "GeneB"])
  let data = @src.AUCellData::new(expression, cell_names, gene_names, [gs])

  let result = @src.run_aucell(data)
  assert_eq(result.cell_names.length(), 2)
  assert_eq(result.gene_set_names.length(), 1)
  assert_eq(result.score_matrix.length(), 2)
  assert_eq(result.score_matrix[0].length(), 1)

  // Cell1 has higher GeneA (10.0) and moderate GeneB (2.0)
  // Cell2 has lower GeneA (3.0) but higher GeneB (8.0)
  let cell1_score = result.score_matrix[0][0]
  let cell2_score = result.score_matrix[1][0]
  // Both should have valid scores
  assert_eq(cell1_score >= 0.0, true)
  assert_eq(cell2_score >= 0.0, true)
}

///|
test "aucell_run_multiple_sets" {
  let expression = [
    [5.0, 10.0, 3.0],
    [8.0, 2.0, 7.0],
  ]
  let cell_names = ["Cell1", "Cell2"]
  let gene_names = ["GeneA", "GeneB", "GeneC"]

  let gs1 = @src.GeneSet::new("Set1", ["GeneA", "GeneB"])
  let gs2 = @src.GeneSet::new("Set2", ["GeneB", "GeneC"])

  let data = @src.AUCellData::new(expression, cell_names, gene_names, [gs1, gs2])
  let result = @src.run_aucell(data)

  assert_eq(result.gene_set_names.length(), 2)
  assert_eq(result.score_matrix[0].length(), 2)
}

///|
test "aucell_normalize_scores" {
  let expression = [
    [10.0, 2.0],
    [3.0, 8.0],
    [5.0, 5.0],
  ]
  let cell_names = ["C1", "C2", "C3"]
  let gene_names = ["G1", "G2"]
  let gs = @src.GeneSet::new("Set1", ["G1"])
  let data = @src.AUCellData::new(expression, cell_names, gene_names, [gs])
  let result = @src.run_aucell(data)
  let normalized = @src.normalize_aucell_scores(result)

  assert_eq(normalized.score_matrix.length(), 3)
  assert_eq(normalized.score_matrix[0].length(), 1)

  // Normalized scores should be in [0, 1]
  let mut i = 0
  while i < 3 {
    let val = normalized.score_matrix[i][0]
    assert_eq(val >= 0.0 && val <= 1.0, true)
    i = i + 1
  }
}

///|
test "aucell_get_cell_scores" {
  let expression = [[10.0, 5.0], [2.0, 8.0]]
  let cell_names = ["C1", "C2"]
  let gene_names = ["G1", "G2"]
  let gs = @src.GeneSet::new("S1", ["G1"])
  let data = @src.AUCellData::new(expression, cell_names, gene_names, [gs])
  let result = @src.run_aucell(data)

  let scores = @src.get_cell_scores(result, "C1")
  assert_eq(scores.length(), 1)
}

///|
test "aucell_get_set_scores" {
  let expression = [[10.0, 5.0], [2.0, 8.0]]
  let cell_names = ["C1", "C2"]
  let gene_names = ["G1", "G2"]
  let gs = @src.GeneSet::new("S1", ["G1"])
  let data = @src.AUCellData::new(expression, cell_names, gene_names, [gs])
  let result = @src.run_aucell(data)

  let scores = @src.get_set_scores(result, "S1")
  assert_eq(scores.length(), 2)
}

///|
test "aucell_gene_not_found" {
  let expression = [[10.0, 5.0], [2.0, 8.0]]
  let cell_names = ["C1", "C2"]
  let gene_names = ["G1", "G2"]
  // Gene set with a gene not in the expression matrix
  let gs = @src.GeneSet::new("S1", ["G1", "NonExistent"])
  let data = @src.AUCellData::new(expression, cell_names, gene_names, [gs])
  let result = @src.run_aucell(data)

  // Should still work with only the found genes
  assert_eq(result.score_matrix.length(), 2)
  assert_eq(result.score_matrix[0][0] >= 0.0, true)
}
