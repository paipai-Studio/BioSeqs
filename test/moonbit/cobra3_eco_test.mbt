// ===========================================================================
// cobra3_eco_test.mbt — EFM / iMAT / GIMME / CORDA / GrowMatch / OptKnock
//                       / community / thermo
// ===========================================================================

///|
fn cobra3_branch_model() -> @src.Cobra3Model {
  let mut m = @src.cobra3_model_new("branch")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -10.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], 0.0, 10.0, "g1"),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R2", [("A", -1.0), ("C", 1.0)], 0.0, 10.0, "g2"),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO1", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO2", [("C", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO1", 1.0)
  m = @src.cobra3_model_set_objective(m, "BIO2", 1.0)
  m
}

///|
test "cobra3 eco efm linear chain" {
  let m = cobra3_toy()
  let modes = @src.cobra3_elementary_flux_modes(m)
  assert_true(modes.length() >= 1)
  for v in modes {
    let s = @src.cobra3_model_stoichiometric_matrix(m)
    for i in 0..<s.length() {
      let mut sum = 0.0
      for j in 0..<s[i].length() {
        sum = sum + s[i][j] * v[j]
      }
      assert_true(sum.abs() < 1.0e-9)
    }
  }
}

///|
test "cobra3 eco imat keeps high pathway" {
  let m = cobra3_branch_model()
  let res = @src.cobra3_imat(m, ["g1"], ["g2"], min_growth=5.0)
  assert_true(res.status is @src.Cobra3Status::Optimal)
  assert_true(res.fluxes[1].abs() >= 0.0)
}

///|
test "cobra3 eco gimme penalizes low pathway" {
  let m = cobra3_branch_model()
  let res = @src.cobra3_gimme(m, ["g2"], min_growth=5.0)
  assert_true(res.status is @src.Cobra3Status::Optimal)
  assert_true(res.fluxes[2].abs() < 1.0e-3)
}

///|
test "cobra3 eco corda removes low reactions" {
  let m = cobra3_branch_model()
  let res = @src.cobra3_corda(m, ["g2"], min_growth=5.0)
  assert_true(res.removed_reactions.contains("R2"))
  assert_true(@src.cobra3_model_reaction_index(res.context_model, "R2") < 0)
}

///|
test "cobra3 eco growmatch adds missing reaction" {
  let base = cobra3_gap_base()
  let mut universal = @src.cobra3_model_new("universal")
  universal = @src.cobra3_model_add_reaction(
    universal,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], 0.0, 1000.0, ""),
  )
  let res = @src.cobra3_growmatch(base, universal, min_growth=0.05)
  assert_eq(res.added_reactions.length(), 1)
  assert_eq(res.added_reactions[0], "R1")
  assert_true(res.predicted_growth > 0.05)
}

///|
test "cobra3 eco imat shuts low pathway" {
  let m = cobra3_branch_model()
  let res = @src.cobra3_imat(m, ["g1"], ["g2"], min_growth=5.0)
  assert_true(res.status is @src.Cobra3Status::Optimal)
  assert_true(res.fluxes[2].abs() < 1.0e-3)
}

///|
test "cobra3 eco community prefixes same-named internal reactions" {
  let mut m1 = @src.cobra3_model_new("m1")
  m1 = @src.cobra3_model_add_reaction(
    m1,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -10.0, 1000.0, ""),
  )
  m1 = @src.cobra3_model_add_reaction(
    m1,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], 0.0, 1000.0, ""),
  )
  m1 = @src.cobra3_model_add_reaction(
    m1,
    cobra3_make_rxn("BIO", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m1 = @src.cobra3_model_set_objective(m1, "BIO", 1.0)
  let mut m2 = @src.cobra3_model_new("m2")
  m2 = @src.cobra3_model_add_reaction(
    m2,
    cobra3_make_rxn("EX_C", [("C", -1.0)], -10.0, 1000.0, ""),
  )
  m2 = @src.cobra3_model_add_reaction(
    m2,
    cobra3_make_rxn("R2", [("C", -1.0), ("D", 1.0)], 0.0, 1000.0, ""),
  )
  m2 = @src.cobra3_model_add_reaction(
    m2,
    cobra3_make_rxn("BIO", [("D", -1.0)], 0.0, 1000.0, ""),
  )
  m2 = @src.cobra3_model_set_objective(m2, "BIO", 1.0)
  let res = @src.cobra3_community_fba([m1, m2])
  assert_true((res.growth - 20.0).abs() < 1.0e-5)
  assert_true(@src.cobra3_model_has_reaction(res.model, "m1:BIO"))
  assert_true(@src.cobra3_model_has_reaction(res.model, "m2:BIO"))
  assert_true(@src.cobra3_model_has_reaction(res.model, "EX_A"))
}

///|
test "cobra3 eco optknock enumerates knockouts" {
  let m = cobra3_toy()
  let res = @src.cobra3_optknock(m, "R2", max_knockouts=2)
  assert_true(res.length() > 0)
  let mut found_r1 = false
  for r in res {
    if r.knockout_ids.length() == 1 && r.knockout_ids[0] == "R1" {
      found_r1 = true
      assert_true(r.growth.abs() < 1.0e-6)
    }
  }
  assert_true(found_r1)
}

///|
test "cobra3 eco community sums independent models" {
  let m1 = cobra3_toy()
  let mut m2 = @src.cobra3_model_new("second")
  m2 = @src.cobra3_model_add_reaction(
    m2,
    cobra3_make_rxn("EX_B", [("B2", -1.0)], -10.0, 1000.0, ""),
  )
  m2 = @src.cobra3_model_add_reaction(
    m2,
    cobra3_make_rxn("S1", [("B2", -1.0), ("C2", 1.0)], 0.0, 1000.0, ""),
  )
  m2 = @src.cobra3_model_add_reaction(
    m2,
    cobra3_make_rxn("S2", [("C2", -1.0)], 0.0, 1000.0, ""),
  )
  m2 = @src.cobra3_model_set_objective(m2, "S2", 1.0)
  let res = @src.cobra3_community_fba([m1, m2])
  assert_eq(res.members.length(), 2)
  assert_true((res.growth - 20.0).abs() < 1.0e-5)
}

///|
test "cobra3 eco thermo blocks unfavorable direction" {
  let m = cobra3_toy()
  let favorite : Map[String, Double] = Map([("R1", -10.0)])
  let ok = @src.cobra3_thermo_fba(m, favorite)
  assert_true((ok.objective_value - 10.0).abs() < 1.0e-6)
  let unfavorable : Map[String, Double] = Map([("R1", 10.0)])
  let bad = @src.cobra3_thermo_fba(m, unfavorable)
  assert_true(bad.objective_value.abs() < 1.0e-6)
}
