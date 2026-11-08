// ===========================================================================
// cobra3_analysis_test.mbt — 通量耦合 / 相平面 / 摘要 / 敏感性 / 模型差异
// ===========================================================================

///|
test "cobra3 analysis flux coupling" {
  let m = cobra3_toy()
  let pairs = @src.cobra3_flux_coupling(m)
  let mut r1_r2 = ""
  let mut ex_r1 = ""
  for p in pairs {
    if p.reaction_a == "R1" && p.reaction_b == "R2" {
      r1_r2 = p.kind
    }
    if p.reaction_a == "EX_A" && p.reaction_b == "R1" {
      ex_r1 = p.kind
    }
  }
  assert_eq(r1_r2, "fully_coupled")
  assert_eq(ex_r1, "directionally_coupled")
}

///|
test "cobra3 analysis phase plane" {
  let m = cobra3_toy()
  let pts = @src.cobra3_phase_plane(m, "R1", "R2", points=3)
  assert_eq(pts.length(), 9)
  // 对角线可行，目标 = R2
  let diag = pts[4]
  assert_true((diag.x_flux - 5.0).abs() < 1.0e-6)
  assert_true((diag.objective - 5.0).abs() < 1.0e-6)
  // 非对角不可行 -> 0
  let off = pts[1]
  assert_true(off.objective.abs() < 1.0e-9)
}

///|
test "cobra3 analysis model summary" {
  let m = cobra3_toy()
  let s = @src.cobra3_model_summary(m)
  assert_eq(s.n_metabolites, 2)
  assert_eq(s.n_reactions, 3)
  assert_eq(s.n_genes, 3)
  assert_eq(s.n_exchange, 1)
  assert_eq(s.objective, "R2")
  assert_true((s.growth - 10.0).abs() < 1.0e-6)
  let text = @src.cobra3_model_summary_string(m)
  assert_true(text.contains("reactions: 3"))
}

///|
test "cobra3 analysis reaction summary" {
  let m = cobra3_toy()
  let s = @src.cobra3_reaction_summary(m, "R2").unwrap()
  assert_true((s.flux - 10.0).abs() < 1.0e-6)
  assert_eq(s.gpr, "g3")
  assert_eq(s.participants.length(), 1)
}

///|
test "cobra3 analysis metabolite summary" {
  let m = cobra3_toy()
  let b = @src.cobra3_metabolite_summary(m, "B").unwrap()
  assert_eq(b.produced_by.length(), 1)
  assert_eq(b.produced_by[0], "R1")
  assert_eq(b.consumed_by.length(), 1)
  assert_eq(b.consumed_by[0], "R2")
}

///|
test "cobra3 analysis sensitivity sorted" {
  let m = cobra3_toy()
  let entries = @src.cobra3_sensitivity(m)
  assert_eq(entries.length(), 2)
  for i in 1..<entries.length() {
    assert_true(
      entries[i - 1].shadow_price.abs() >= entries[i].shadow_price.abs(),
    )
  }
}

///|
test "cobra3 analysis model diff" {
  let a = cobra3_toy()
  let b = @src.cobra3_model_set_reaction_bounds(a, "EX_A", -5.0, 1000.0)
  let c = @src.cobra3_model_remove_reaction(b, "R1")
  let d1 = @src.cobra3_model_diff(a, b)
  assert_eq(d1.changed_bounds.length(), 1)
  assert_eq(d1.changed_bounds[0], "EX_A")
  let d2 = @src.cobra3_model_diff(a, c)
  assert_eq(d2.removed_reactions.length(), 1)
  assert_eq(d2.removed_reactions[0], "R1")
}
