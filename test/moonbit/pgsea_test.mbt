///|
/// Tests for pgsea module (parametric gene set enrichment analysis).

test "pgsea_score_matrix_dims" {
  let expr = [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [5.0, 4.0, 3.0, 2.0],
  ]
  let genes = ["G1", "G2", "G3"]
  let sets = [("S1", ["G1", "G2"]), ("S2", ["G3"])]
  let (matrix, names) = @src.pgsea_score_matrix(expr, genes, sets)
  assert_eq(names.length(), 2)
  assert_eq(matrix.length(), 2)
  assert_eq(matrix[0].length(), 4)
}

test "pgsea_run_differential_sets" {
  // G1, G2 up in treatment; G3 down in treatment
  // Different variances to avoid zero within-group variance after z-scoring
  let expr = [
    [1.0, 2.0, 3.0, 4.0],
    [1.5, 1.5, 3.5, 3.5],
    [4.0, 3.0, 2.0, 1.0],
  ]
  let genes = ["G1", "G2", "G3"]
  let sets = [("UP", ["G1", "G2"]), ("DOWN", ["G3"])]
  let labels = ["control", "control", "treat", "treat"]
  let results = @src.pgsea_run(expr, genes, sets, labels, "control", "treat")
  assert_eq(results.length(), 2)
  // UP set should have positive t (treatment higher); DOWN negative
  assert_true(results[0].t_stat > 0.0)
  assert_true(results[1].t_stat < 0.0)
}

test "pgsea_t_pvalue_small_for_large_t" {
  let p_small = @src.pgsea_t_pvalue(5.0, 10.0)
  let p_large = @src.pgsea_t_pvalue(0.5, 10.0)
  assert_true(p_small < 0.05)
  assert_true(p_small < p_large)
}

test "pgsea_score_one_zero_sd" {
  // single gene, constant expression -> sd = 0, score should be 0
  let expr = [[2.0, 2.0, 2.0]]
  let genes = ["G1"]
  let means = [2.0]
  let sds = [0.0]
  let score = @src.pgsea_score_one(expr, genes, ["G1"], 0, means, sds)
  assert_true(score.abs() < 1.0e-9)
}

test "pgsea_top_ranking" {
  let r1 = @src.PgseaResult::{ gene_set: "A", scores: [1.0], mean_control: 0.0, mean_treatment: 1.0, t_stat: 2.0, p_value: 0.1 }
  let r2 = @src.PgseaResult::{ gene_set: "B", scores: [1.0], mean_control: 0.0, mean_treatment: 1.0, t_stat: 3.0, p_value: 0.01 }
  let top = @src.pgsea_top([r1, r2], 1)
  assert_eq(top.length(), 1)
  assert_eq(top[0].gene_set, "B")
}

test "pgsea_export" {
  let r = @src.PgseaResult::{ gene_set: "S1", scores: [1.0], mean_control: 0.0, mean_treatment: 1.0, t_stat: 2.0, p_value: 0.1 }
  let tsv = @src.pgsea_export([r])
  assert_true(tsv.contains("gene_set"))
  assert_true(tsv.contains("S1"))
}
