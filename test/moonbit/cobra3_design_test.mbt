// ===========================================================================
// cobra3_design_test.mbt — gapfill / minimal_medium / 包络 / MCS / 采样 / FastCC
// ===========================================================================

///|
fn cobra3_gap_base() -> @src.Cobra3Model {
  let mut m = @src.cobra3_model_new("gapbase")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -10.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R2", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "R2", 1.0)
  m
}

///|
test "cobra3 design gapfill restores growth" {
  let base = cobra3_gap_base()
  assert_true(@src.cobra3_model_optimize(base).objective_value.abs() < 1.0e-9)
  let mut universal = @src.cobra3_model_new("universal")
  universal = @src.cobra3_model_add_reaction(
    universal,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], 0.0, 1000.0, ""),
  )
  let res = @src.cobra3_gapfill(base, universal, lower_bound=0.05)
  assert_true(res.status is @src.Cobra3Status::Optimal)
  assert_eq(res.added_reactions.length(), 1)
  assert_eq(res.added_reactions[0], "R1")
  assert_true(res.objective_value >= 0.05 - 1.0e-6)
}

///|
test "cobra3 design minimal medium" {
  let m = cobra3_toy()
  let res = @src.cobra3_minimal_medium(m, min_growth=5.0, max_import=10.0)
  assert_true(res.status is @src.Cobra3Status::Optimal)
  assert_eq(res.n_components, 1)
  assert_true(res.medium["EX_A"] > 1.0e-9)
  assert_true(res.medium["EX_A"] <= 10.0 + 1.0e-9)
  assert_true(res.objective_value >= 5.0 - 1.0e-6)
}

///|
test "cobra3 design production envelope" {
  let m = cobra3_toy()
  let env = @src.cobra3_production_envelope(m, "R2", points=3)
  assert_eq(env.length(), 3)
  let mid = env[1]
  assert_true((mid.target_flux - 5.0).abs() < 1.0e-6)
  assert_true((mid.objective_max - 5.0).abs() < 1.0e-6)
}

///|
test "cobra3 design minimal cut sets" {
  let m = cobra3_toy()
  let mcs = @src.cobra3_minimal_cut_sets(m, threshold=0.01)
  assert_true(mcs.length() > 0)
  assert_true(mcs[0].size == 1)
}

///|
test "cobra3 design minimal cut sets includes disjoint pairs" {
  let mut m = @src.cobra3_model_new("mcs")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -20.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R2", [("A", -1.0), ("B", 1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R3", [("A", -1.0), ("C", 1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO", [("B", -1.0), ("C", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO", 1.0)
  let mcs = @src.cobra3_minimal_cut_sets(m, threshold=0.01)
  let mut has_r3 = false
  let mut has_pair = false
  for s in mcs {
    if s.size == 1 && s.reaction_ids[0] == "R3" {
      has_r3 = true
    }
    if s.size == 2 &&
      s.reaction_ids.contains("R1") &&
      s.reaction_ids.contains("R2") {
      has_pair = true
    }
  }
  assert_true(has_r3)
  assert_true(has_pair)
  // 含单割集元素的双组合不是极小割集
  for s in mcs {
    if s.size == 2 {
      assert_true(!s.reaction_ids.contains("R3"))
    }
  }
}

///|
test "cobra3 design minimal cut sets supports triple knockouts" {
  let mut m = @src.cobra3_model_new("mcs3")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -20.0, 1000.0, ""),
  )
  for id in ["P1", "P2", "P3"] {
    m = @src.cobra3_model_add_reaction(
      m,
      cobra3_make_rxn(id, [("A", -1.0), ("P", 1.0)], 0.0, 1000.0, ""),
    )
  }
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO", [("P", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO", 1.0)
  // 任意两个产 P 反应被敲除后仍有通路，生长不降，故 size<=2 找不到该割集
  let two = @src.cobra3_minimal_cut_sets(m, threshold=0.01, max_size=2)
  for s in two {
    assert_true(s.size <= 2)
  }
  let three = @src.cobra3_minimal_cut_sets(m, threshold=0.01, max_size=3)
  let mut has_triple = false
  for s in three {
    if s.size == 3 &&
      s.reaction_ids.contains("P1") &&
      s.reaction_ids.contains("P2") &&
      s.reaction_ids.contains("P3") {
      has_triple = true
    }
  }
  assert_true(has_triple)
  // max_candidates 限流：只允许测试一个组合
  let capped = @src.cobra3_minimal_cut_sets(
    m,
    threshold=0.01,
    max_size=3,
    max_candidates=1,
  )
  assert_true(capped.length() <= 1)
}

///|
test "cobra3 design minimal cut sets honors exclude" {
  let mut m = @src.cobra3_model_new("mcs_exclude")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -20.0, 1000.0, ""),
  )
  for id in ["P1", "P2", "P3"] {
    m = @src.cobra3_model_add_reaction(
      m,
      cobra3_make_rxn(id, [("A", -1.0), ("P", 1.0)], 0.0, 1000.0, ""),
    )
  }
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO", [("P", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO", 1.0)
  let all = @src.cobra3_minimal_cut_sets(m, threshold=0.01, max_size=3)
  let mut all_has_ex_a = false
  for s in all {
    if s.reaction_ids.contains("EX_A") {
      all_has_ex_a = true
    }
  }
  assert_true(all_has_ex_a)
  let mcs = @src.cobra3_minimal_cut_sets(m, threshold=0.01, max_size=3, exclude=[
    "EX_A",
  ])
  let mut has_bio = false
  let mut has_triple = false
  for s in mcs {
    assert_true(!s.reaction_ids.contains("EX_A"))
    if s.size == 1 && s.reaction_ids[0] == "BIO" {
      has_bio = true
    }
    if s.size == 3 &&
      s.reaction_ids.contains("P1") &&
      s.reaction_ids.contains("P2") &&
      s.reaction_ids.contains("P3") {
      has_triple = true
    }
  }
  assert_true(has_bio)
  assert_true(has_triple)
}

///|
test "cobra3 design minimal cut sets milp finds triple cut" {
  let mut m = @src.cobra3_model_new("mcs_milp3")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -20.0, 1000.0, ""),
  )
  for id in ["P1", "P2", "P3"] {
    m = @src.cobra3_model_add_reaction(
      m,
      cobra3_make_rxn(id, [("A", -1.0), ("P", 1.0)], 0.0, 1000.0, ""),
    )
  }
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO", [("P", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO", 1.0)
  // 单层 MILP 若把"存在低目标通量"写错，会被零通量平凡满足而返回空；
  // 这里排除 BIO/EX_A 后唯一最小割集必须同时敲掉三条并行通路。
  let mcs = @src.cobra3_minimal_cut_sets_milp(m, threshold=0.01, exclude=[
    "BIO", "EX_A",
  ])
  assert_eq(mcs.length(), 1)
  assert_eq(mcs[0].size, 3)
  assert_true(mcs[0].reaction_ids.contains("P1"))
  assert_true(mcs[0].reaction_ids.contains("P2"))
  assert_true(mcs[0].reaction_ids.contains("P3"))
  assert_true(mcs[0].growth <= 0.01 + 1.0e-6)
}

///|
test "cobra3 design minimal cut sets milp returns empty when no cut needed" {
  let m = cobra3_toy()
  let mcs = @src.cobra3_minimal_cut_sets_milp(m, threshold=100.0)
  assert_eq(mcs.length(), 0)
}

///|
test "cobra3 design nullspace dimension" {
  let m = cobra3_toy()
  let basis = @src.cobra3_nullspace_basis(m)
  assert_eq(basis.length(), 1)
  let s = @src.cobra3_model_stoichiometric_matrix(m)
  // S * v = 0
  for i in 0..<s.length() {
    let mut sum = 0.0
    for j in 0..<s[i].length() {
      sum = sum + s[i][j] * basis[0][j]
    }
    assert_true(sum.abs() < 1.0e-9)
  }
}

///|
test "cobra3 design sampling respects steady state" {
  let m = cobra3_toy()
  let samples = @src.cobra3_sample(m, n_samples=12, seed=7, warmup=4)
  assert_eq(samples.length(), 12)
  let s = @src.cobra3_model_stoichiometric_matrix(m)
  for sample in samples {
    for i in 0..<s.length() {
      let mut sum = 0.0
      for j in 0..<s[i].length() {
        sum = sum + s[i][j] * sample[j]
      }
      assert_true(sum.abs() < 1.0e-6)
    }
    assert_true(sample[0] >= -10.0 - 1.0e-9)
    assert_true(sample[2] >= -1.0e-9)
  }
}

///|
test "cobra3 design sample flux statistics" {
  let m = cobra3_toy()
  let samples = @src.cobra3_sample(m, n_samples=20, seed=3, warmup=5)
  let dist = @src.cobra3_sample_flux(m, samples)
  assert_eq(dist.length(), 3)
  assert_true(dist[2].minimum >= -1.0e-9)
  assert_true(dist[2].maximum <= 10.0 + 1.0e-9)
  assert_true(dist[2].mean >= 0.0)
}

///|
test "cobra3 design sampling diagnostics" {
  let samples : Array[Double] = Array(capacity=50)
  for i in 0..<50 {
    samples.push((i * 7 % 11).to_double())
  }
  let diag = @src.cobra3_sampling_diagnostics(samples)
  assert_true(diag.r_hat > 0.0)
  assert_true(diag.effective_sample_size > 0.0)
}

///|
fn cobra3_test_rand_next(state : Int32) -> Int32 {
  let s = (state * 1103515245 + 12345) % 2147483647
  if s < 0 {
    -s
  } else {
    s
  }
}

///|
test "cobra3 design ess near n for iid samples" {
  let n = 2000
  let samples : Array[Double] = Array(capacity=n)
  let mut state = 7
  for _ in 0..<n {
    state = cobra3_test_rand_next(state)
    samples.push(state.to_double() / 2147483647.0)
  }
  let diag = @src.cobra3_sampling_diagnostics(samples)
  assert_true(diag.effective_sample_size > n.to_double() * 0.3)
  assert_true(diag.effective_sample_size <= n.to_double() + 1.0e-9)
  assert_true(diag.geweke_z.abs() < 5.0)
}

///|
test "cobra3 design ess small for autocorrelated samples" {
  let n = 2000
  let samples : Array[Double] = Array(capacity=n)
  let mut state = 11
  let mut x = 0.0
  for _ in 0..<n {
    state = cobra3_test_rand_next(state)
    let noise = state.to_double() / 2147483647.0 - 0.5
    x = 0.95 * x + noise
    samples.push(x)
  }
  let diag = @src.cobra3_sampling_diagnostics(samples)
  assert_true(diag.effective_sample_size < n.to_double() * 0.25)
  assert_true(diag.effective_sample_size >= 1.0)
}

///|
test "cobra3 design geweke detects mean drift" {
  let n = 1000
  let samples : Array[Double] = Array(capacity=n)
  let mut state = 3
  for i in 0..<n {
    state = cobra3_test_rand_next(state)
    let noise = (state.to_double() / 2147483647.0 - 0.5) * 0.2
    let base = if i < n / 10 { 0.0 } else { 5.0 }
    samples.push(base + noise)
  }
  let diag = @src.cobra3_sampling_diagnostics(samples)
  assert_true(diag.geweke_z.abs() > 5.0)
}

///|
test "cobra3 design sampling diagnostics chains" {
  let chains : Array[Array[Double]] = Array(capacity=4)
  for _ in 0..<4 {
    let xs : Array[Double] = Array(capacity=50)
    for i in 0..<50 {
      xs.push((i * 7 % 11).to_double())
    }
    chains.push(xs)
  }
  let same = @src.cobra3_sampling_diagnostics_chains(chains)
  assert_true(same.r_hat >= 1.0 - 1.0e-9)
  assert_true(same.r_hat < 1.0 + 1.0e-9)
  let shifted : Array[Array[Double]] = Array(capacity=2)
  let a : Array[Double] = Array(capacity=20)
  let b : Array[Double] = Array(capacity=20)
  for i in 0..<20 {
    a.push(i.to_double())
    b.push((i + 100).to_double())
  }
  shifted.push(a)
  shifted.push(b)
  let diff = @src.cobra3_sampling_diagnostics_chains(shifted)
  assert_true(diff.r_hat > 1.5)
}

///|
test "cobra3 design sample chains shape" {
  let m = cobra3_toy()
  let chains = @src.cobra3_sample_chains(m, n_chains=2, n_samples=10, warmup=5)
  assert_eq(chains.length(), 2)
  assert_eq(chains[0].length(), 10)
  assert_eq(chains[0][0].length(), 3)
}

///|
test "cobra3 design fastcc keeps viable reactions" {
  let m = cobra3_toy()
  let res = @src.cobra3_fastcc(m)
  assert_eq(res.removed_reactions.length(), 0)
  assert_eq(res.consistent_model.reactions.length(), 3)
}

///|
test "cobra3 design detect loops on linear chain" {
  let m = cobra3_toy()
  let loops = @src.cobra3_detect_loops(m)
  assert_eq(loops.length(), 0)
}

///|
test "cobra3 design loopless fba keeps growth" {
  let m = cobra3_toy()
  let res = @src.cobra3_loopless_fba(m)
  assert_true(res.status is @src.Cobra3Status::Optimal)
  assert_true((res.objective_value - 10.0).abs() < 1.0e-6)
}

///|
test "cobra3 design loopless fba removes internal cycle" {
  let mut m = @src.cobra3_model_new("cycle")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -10.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R2", [("B", -1.0), ("A", 1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO", 1.0)
  let res = @src.cobra3_loopless_fba(m)
  assert_true(res.status is @src.Cobra3Status::Optimal)
  assert_true((res.objective_value - 10.0).abs() < 1.0e-5)
  let r2 = @src.cobra3_model_reaction_index(m, "R2")
  let r1 = @src.cobra3_model_reaction_index(m, "R1")
  assert_true(res.fluxes[r2].abs() < 1.0e-5)
  assert_true((res.fluxes[r1] - 10.0).abs() < 1.0e-5)
}
