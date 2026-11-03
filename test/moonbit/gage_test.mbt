///|
/// Tests for GAGE module.

///|
test "gage_gene_set_creation" {
  let gs = @src.GageGeneSet::new("GS1", "Test Pathway", ["Gene1", "Gene2"])
  assert_eq(gs.id, "GS1")
  assert_eq(gs.name, "Test Pathway")
  assert_eq(gs.gene_ids.length(), 2)
}

///|
test "gage_run_basic" {
  let fcs = @src.gage_sample_fold_changes()
  let sets = @src.gage_sample_gene_sets()
  let results = @src.run_gage(fcs, sets)
  assert_eq(results.n_gene_sets, 3)
  assert_true(results.results.length() > 0)
}

///|
test "gage_mean_fc" {
  let fcs = @src.gage_sample_fold_changes()
  let sets = @src.gage_sample_gene_sets()
  let results = @src.run_gage(fcs, sets)
  // GS1 contains Gene1(2.0), Gene2(1.5), Gene3(1.8) => mean = 1.767
  let mut gs1_result : @src.GageResult? = None
  let mut i = 0
  while i < results.results.length() {
    if results.results[i].gene_set_id == "GS1" {
      gs1_result = Some(results.results[i])
    }
    i = i + 1
  }
  match gs1_result {
    Some(r) => {
      assert_true(r.mean_fc > 1.5)
      assert_true(r.mean_fc < 2.0)
      assert_true(r.is_up)
    }
    None => assert_true(false)
  }
}

///|
test "gage_downregulated_set" {
  let fcs = @src.gage_sample_fold_changes()
  let sets = @src.gage_sample_gene_sets()
  let results = @src.run_gage(fcs, sets)
  // GS2 contains Gene4(-1.2), Gene5(-0.8), Gene6(-1.5) => mean negative
  let mut gs2_result : @src.GageResult? = None
  let mut i = 0
  while i < results.results.length() {
    if results.results[i].gene_set_id == "GS2" {
      gs2_result = Some(results.results[i])
    }
    i = i + 1
  }
  match gs2_result {
    Some(r) => {
      assert_true(r.mean_fc < 0.0)
      assert_true(r.is_down)
    }
    None => assert_true(false)
  }
}

///|
test "gage_pvalue_range" {
  let fcs = @src.gage_sample_fold_changes()
  let sets = @src.gage_sample_gene_sets()
  let results = @src.run_gage(fcs, sets)
  let mut i = 0
  while i < results.results.length() {
    assert_true(results.results[i].p_val >= 0.0)
    assert_true(results.results[i].p_val <= 1.0)
    assert_true(results.results[i].p_adj >= 0.0)
    assert_true(results.results[i].p_adj <= 1.0)
    i = i + 1
  }
}

///|
test "gage_get_significant" {
  let fcs = @src.gage_sample_fold_changes()
  let sets = @src.gage_sample_gene_sets()
  let results = @src.run_gage(fcs, sets, fdr_threshold=1.0) // All significant
  let sig = results.get_significant()
  assert_eq(sig.length(), results.n_significant)
}

///|
test "gage_get_upregulated" {
  let fcs = @src.gage_sample_fold_changes()
  let sets = @src.gage_sample_gene_sets()
  let results = @src.run_gage(fcs, sets)
  let up = results.get_upregulated()
  assert_true(up.length() > 0)
}

///|
test "gage_get_downregulated" {
  let fcs = @src.gage_sample_fold_changes()
  let sets = @src.gage_sample_gene_sets()
  let results = @src.run_gage(fcs, sets)
  let down = results.get_downregulated()
  assert_true(down.length() > 0)
}

///|
test "gage_top_sets" {
  let fcs = @src.gage_sample_fold_changes()
  let sets = @src.gage_sample_gene_sets()
  let results = @src.run_gage(fcs, sets)
  let top = results.get_top_sets(2)
  assert_eq(top.length(), 2)
  // Top set should have smaller p_adj than second
  assert_true(top[0].p_adj <= top[1].p_adj)
}

///|
test "gage_summary" {
  let fcs = @src.gage_sample_fold_changes()
  let sets = @src.gage_sample_gene_sets()
  let results = @src.run_gage(fcs, sets)
  let s = results.summary()
  assert_true(s.length() > 0)
}

///|
test "gage_paired_test" {
  let fcs = @src.gage_sample_fold_changes()
  let sets = @src.gage_sample_gene_sets()
  let results = @src.run_gage(fcs, sets, paired=true)
  assert_true(results.paired)
}
