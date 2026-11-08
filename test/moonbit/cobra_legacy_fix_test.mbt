// ===========================================================================
// cobra_legacy_fix_test.mbt — 旧 cobra 模块语义修复验证
// ===========================================================================

///|
fn legacy_branch_model() -> @src.FbaModel raise {
  let mut m = @src.fba_model_new(["A", "B"])
  m = m.add_reaction(
    {
      id: "EX_A",
      name: "EX_A",
      lower_bound: -10.0,
      upper_bound: 1000.0,
      objective: 0.0,
      genes: [],
    },
    [("A", -1.0)],
  )
  m = m.add_reaction(
    {
      id: "R1",
      name: "R1",
      lower_bound: 0.0,
      upper_bound: 10.0,
      objective: 0.0,
      genes: [],
    },
    [("A", -1.0), ("B", 1.0)],
  )
  m = m.add_reaction(
    {
      id: "BIO",
      name: "BIO",
      lower_bound: 0.0,
      upper_bound: 1000.0,
      objective: 1.0,
      genes: [],
    },
    [("B", -1.0)],
  )
  m
}

///|
test "legacy fva fraction of optimum now works" {
  let m = legacy_branch_model()
  let full = @src.cobra_flux_variability_analysis(m, fraction_of_optimum=1.0)
  let zero = @src.cobra_flux_variability_analysis(m, fraction_of_optimum=0.0)
  let mut full_min = -1.0
  let mut full_max = -1.0
  let mut zero_min = -1.0
  let mut zero_max = -1.0
  for f in full {
    if f.reaction_id == "R1" {
      full_min = f.minimum
      full_max = f.maximum
    }
  }
  for f in zero {
    if f.reaction_id == "R1" {
      zero_min = f.minimum
      zero_max = f.maximum
    }
  }
  assert_true((full_min - 10.0).abs() < 1.0e-5)
  assert_true((full_max - 10.0).abs() < 1.0e-5)
  assert_true(zero_min.abs() < 1.0e-5)
  assert_true((zero_max - 10.0).abs() < 1.0e-5)
}

///|
test "legacy model from json no longer empty" {
  let m3 = cobra3_toy()
  let text = @src.cobra3_model_to_json(m3)
  let fba = @src.cobra2_model_from_json(text)
  assert_eq(fba.reactions.length(), 3)
  assert_eq(fba.metabolites.length(), 2)
}

///|
test "legacy sbml read delegates to fbc parser" {
  let m3 = cobra3_toy()
  let xml = @src.cobra3_model_to_sbml(m3)
  let fba = @src.cobra_read_sbml(xml).unwrap()
  assert_eq(fba.reactions.length(), 3)
  let idx = {
    let mut found = -1
    for i in 0..<fba.reactions.length() {
      if fba.reactions[i].id == "R2" {
        found = i
      }
    }
    found
  }
  assert_true(idx >= 0)
  assert_true((fba.reactions[idx].objective - 1.0).abs() < 1.0e-9)
}

///|
test "legacy minimal medium uses cobra3 milp" {
  let m3 = cobra3_toy()
  let fba = @src.cobra3_to_fba_model(m3)
  let res = @src.cobra_minimal_medium(fba, min_objective=5.0, _max_import=10.0)
  assert_eq(res.n_components, 1)
}

///|
test "legacy enhanced solution uses exact duals" {
  let m3 = cobra3_toy()
  let fba = @src.cobra3_to_fba_model(m3)
  let sol = @src.fba_solve(fba)
  let enhanced = @src.cobra_enhanced_solution(fba, sol)
  assert_eq(enhanced.shadow_prices.length(), 2)
  assert_eq(enhanced.reduced_costs.length(), 3)
  assert_eq(enhanced.fluxes.length(), 3)
}
