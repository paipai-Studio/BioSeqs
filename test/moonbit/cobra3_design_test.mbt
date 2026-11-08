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
