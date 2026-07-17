///|
test "edger_dge_list" {
  let counts = [[10, 20, 30, 40], [50, 60, 70, 80], [90, 100, 110, 120]]
  let group = ["A", "A", "B", "B"]
  let genes = ["gene1", "gene2", "gene3"]
  let dge = @src.dge_list(counts, group, genes)

  assert_eq(dge.counts.length(), 3)
  assert_eq(dge.group.length(), 4)
  assert_eq(dge.genes.length(), 3)
}

///|
test "edger_calc_norm_factors" {
  let counts = [[100, 200, 300], [400, 500, 600]]
  let group = ["A", "B", "C"]
  let genes = ["gene1", "gene2"]
  let dge = @src.dge_list(counts, group, genes)
  let dge_norm = @src.calc_norm_factors(dge)

  assert_eq(dge_norm.norm_factors.length(), 3)
}

///|
test "edger_estimate_dispersion" {
  let counts = [[10, 20, 30], [40, 50, 60], [70, 80, 90]]
  let group = ["A", "B", "C"]
  let genes = ["gene1", "gene2", "gene3"]
  let dge = @src.dge_list(counts, group, genes)
  let disp = @src.estimate_dispersion(dge)

  assert_true(disp >= 0.0)
}

///|
test "edger_exact_test" {
  let counts = [[10, 20], [30, 40], [50, 60]]
  let group = ["A", "B"]
  let genes = ["gene1", "gene2", "gene3"]
  let dge = @src.dge_list(counts, group, genes)
  let dge_norm = @src.calc_norm_factors(dge)
  let results = @src.exact_test(dge_norm, "A", "B")

  assert_eq(results.genes.length(), 3)
  assert_eq(results.log_fc.length(), 3)
}

///|
test "edger_top_tags" {
  let counts = [[10, 20], [30, 40], [50, 60]]
  let group = ["A", "B"]
  let genes = ["gene1", "gene2", "gene3"]
  let dge = @src.dge_list(counts, group, genes)
  let dge_norm = @src.calc_norm_factors(dge)
  let results = @src.exact_test(dge_norm, "A", "B")
  let top = @src.top_tags(results, 2)

  assert_eq(top.genes.length(), 2)
}

///|
test "edger_glm_fit" {
  let counts = [[10, 20], [30, 40], [50, 60]]
  let group = ["A", "B"]
  let genes = ["gene1", "gene2", "gene3"]
  let dge = @src.dge_list(counts, group, genes)
  let design = [[1.0, 0.0], [1.0, 1.0]]
  let fit = @src.glm_fit(dge, design)

  assert_eq(fit.length(), 3)
}

///|
test "edger_glm_lrt" {
  let counts = [[10, 20], [30, 40], [50, 60]]
  let group = ["A", "B"]
  let genes = ["gene1", "gene2", "gene3"]
  let dge = @src.dge_list(counts, group, genes)
  let design = [[1.0, 0.0], [1.0, 1.0]]
  let results = @src.glm_lrt(dge, design, 1)

  assert_eq(results.genes.length(), 3)
}
