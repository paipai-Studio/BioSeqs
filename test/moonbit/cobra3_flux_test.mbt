// ===========================================================================
// cobra3_flux_test.mbt — FBA/pFBA/FVA/敲除/MOMA/ROOM 测试
// ===========================================================================

///|
fn cobra3_make_rxn(
  id : String,
  participants : Array[(String, Double)],
  lb : Double,
  ub : Double,
  gpr_text : String,
) -> @src.Cobra3Reaction {
  let mets : Map[String, Double] = Map(participants)
  {
    id,
    name: id,
    metabolites: mets,
    lower_bound: lb,
    upper_bound: ub,
    objective: 0.0,
    subsystem: "",
    gpr: if gpr_text.length() == 0 {
      None
    } else {
      @src.cobra3_gpr_parse(gpr_text)
    },
    annotation: Map([], capacity=1),
  }
}

///|
fn cobra3_toy() -> @src.Cobra3Model {
  let mut m = @src.cobra3_model_new("toy")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -10.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], 0.0, 1000.0, "g1 or g2"),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R2", [("B", -1.0)], 0.0, 1000.0, "g3"),
  )
  m = @src.cobra3_model_set_objective(m, "R2", 1.0)
  m
}

///|
test "cobra3 fba toy growth" {
  let m = cobra3_toy()
  let sol = @src.cobra3_model_optimize(m)
  assert_true(sol.status is @src.Cobra3Status::Optimal)
  assert_true((sol.objective_value - 10.0).abs() < 1.0e-6)
  let fluxes = @src.cobra3_flux_map(m, sol)
  assert_true((fluxes["EX_A"] + 10.0).abs() < 1.0e-6)
  assert_true((fluxes["R1"] - 10.0).abs() < 1.0e-6)
  assert_true((fluxes["R2"] - 10.0).abs() < 1.0e-6)
}

///|
test "cobra3 fba slim optimize" {
  let m = cobra3_toy()
  assert_true((@src.cobra3_model_slim_optimize(m) - 10.0).abs() < 1.0e-6)
}

///|
test "cobra3 pfba minimizes total flux" {
  let m = cobra3_toy()
  let sol = @src.cobra3_pfba(m)
  assert_true(sol.status is @src.Cobra3Status::Optimal)
  assert_true((sol.objective_value - 10.0).abs() < 1.0e-6)
  let mut total = 0.0
  for v in sol.fluxes {
    total = total + v.abs()
  }
  assert_true((total - 30.0).abs() < 1.0e-4)
}

///|
test "cobra3 fva with fraction of optimum" {
  let mut m = @src.cobra3_model_new("branch")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -10.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], 0.0, 10.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R2", [("A", -1.0), ("C", 1.0)], 0.0, 10.0, ""),
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
  m = @src.cobra3_model_set_objective(m, "BIO2", 0.0)
  let full = @src.cobra3_fva(m, fraction_of_optimum=1.0)
  let zero = @src.cobra3_fva(m, fraction_of_optimum=0.0)
  let mut full_r1_min = 0.0
  let mut full_r1_max = 0.0
  let mut zero_r1_min = -1.0
  let mut zero_r1_max = -1.0
  for f in full {
    if f.reaction_id == "R1" {
      full_r1_min = f.minimum
      full_r1_max = f.maximum
    }
  }
  for f in zero {
    if f.reaction_id == "R1" {
      zero_r1_min = f.minimum
      zero_r1_max = f.maximum
    }
  }
  assert_true((full_r1_min - 10.0).abs() < 1.0e-5)
  assert_true((full_r1_max - 10.0).abs() < 1.0e-5)
  assert_true(zero_r1_min.abs() < 1.0e-5)
  assert_true((zero_r1_max - 10.0).abs() < 1.0e-5)
}

///|
test "cobra3 gene deletion respects gpr or" {
  let m = cobra3_toy()
  let res = @src.cobra3_single_gene_deletion(m)
  let mut g1_growth = -1.0
  let mut g2_growth = -1.0
  let mut g3_growth = -1.0
  for r in res {
    if r.id == "g1" {
      g1_growth = r.growth
    }
    if r.id == "g2" {
      g2_growth = r.growth
    }
    if r.id == "g3" {
      g3_growth = r.growth
    }
  }
  assert_true((g1_growth - 10.0).abs() < 1.0e-6)
  assert_true((g2_growth - 10.0).abs() < 1.0e-6)
  assert_true(g3_growth.abs() < 1.0e-6)
}

///|
test "cobra3 double gene deletion kills backup" {
  let m = cobra3_toy()
  let res = @src.cobra3_double_gene_deletion(m)
  let mut found = false
  let mut growth = -1.0
  for r in res {
    if r.id == "g1|g2" {
      found = true
      growth = r.growth
    }
  }
  assert_true(found)
  assert_true(growth.abs() < 1.0e-6)
  assert_eq(res.length(), 3)
}

///|
test "cobra3 reaction deletion" {
  let m = cobra3_toy()
  let res = @src.cobra3_single_reaction_deletion(m)
  let mut r1_growth = -1.0
  for r in res {
    if r.id == "R1" {
      r1_growth = r.growth
    }
  }
  assert_true(r1_growth.abs() < 1.0e-6)
}

///|
test "cobra3 quadratic moma" {
  let m = cobra3_toy()
  let res = @src.cobra3_moma(m, ["g1", "g2"])
  assert_true(res.status is @src.Cobra3Status::Optimal)
  assert_true(res.growth.abs() < 1.0e-4)
}

///|
test "cobra3 moma multi-target growth" {
  let m = cobra3_branch_model()
  let res = @src.cobra3_moma(m, ["g1"])
  assert_true(res.status is @src.Cobra3Status::Optimal)
  assert_true((res.growth - 10.0).abs() < 1.0e-4)
}

///|
test "cobra3 room keeps backup route" {
  let m = cobra3_toy()
  let res = @src.cobra3_room(m, ["g1"])
  assert_true(res.status is @src.Cobra3Status::Optimal)
  assert_true((res.growth - 10.0).abs() < 1.0e-4)
}

///|
test "cobra3 room enforces min growth" {
  let m = cobra3_branch_model()
  let ok = @src.cobra3_room(m, ["g1"], min_growth=9.0)
  assert_true(ok.status is @src.Cobra3Status::Optimal)
  assert_true(ok.growth >= 9.0 - 1.0e-6)
  let bad = @src.cobra3_room(m, ["g1"], min_growth=100.0)
  assert_true(bad.status is @src.Cobra3Status::Infeasible)
}

///|
test "cobra3 geometric fba" {
  let m = cobra3_toy()
  let res = @src.cobra3_geometric_fba(m)
  assert_true(res.status is @src.Cobra3Status::Optimal)
  // 99% 最优生长：解的通量应为 9.9，objective_value 应与通量一致
  assert_true((res.objective_value - 9.9).abs() < 1.0e-4)
  let mut total = 0.0
  for v in res.fluxes {
    total = total + v * v
  }
  // 99% 最优生长 -> 通量 9.9，Σv² = 3 * 9.9² = 294.03
  assert_true((total - 294.03).abs() < 0.1)
}

///|
test "cobra3 fva zero reaction in blocked model" {
  let m = cobra3_toy()
  let fva = @src.cobra3_fva(m, fraction_of_optimum=1.0)
  assert_eq(fva.length(), 3)
  for f in fva {
    assert_true((f.minimum - f.maximum).abs() < 1.0e-6)
  }
}

///|
test "cobra3 fva reports infeasible status" {
  let m = cobra3_toy()
  let fva = @src.cobra3_fva(m, fraction_of_optimum=1.5)
  assert_true(fva.length() > 0)
  for f in fva {
    assert_true(f.status is @src.Cobra3Status::Infeasible)
  }
}
