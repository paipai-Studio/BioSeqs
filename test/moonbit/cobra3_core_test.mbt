// ===========================================================================
// cobra3_core_test.mbt — COBRApy 对象模型 v2 测试
// ===========================================================================

///|
fn cobra3_test_model() -> @src.Cobra3Model {
  let mut model = @src.cobra3_model_new("toy")
  model = @src.cobra3_model_add_reaction_simple(
    model,
    "EX_A",
    [("A", 1.0)],
    -10.0,
    1000.0,
  )
  model = @src.cobra3_model_add_reaction_simple(
    model,
    "R1",
    [("A", -1.0), ("B", 1.0)],
    0.0,
    1000.0,
  )
  model = @src.cobra3_model_add_reaction_simple(
    model,
    "R2",
    [("B", -1.0)],
    0.0,
    1000.0,
  )
  model = @src.cobra3_model_set_objective(model, "R2", 1.0)
  model
}

///|
test "cobra3 core model lookup and counts" {
  let model = cobra3_test_model()
  assert_eq(model.metabolites.length(), 2)
  assert_eq(model.reactions.length(), 3)
  assert_true(@src.cobra3_model_reaction_index(model, "R1") == 1)
  assert_true(@src.cobra3_model_reaction_index(model, "nope") < 0)
  match @src.cobra3_model_get_reaction(model, "R1") {
    Some(r) => assert_eq(r.id, "R1")
    None => fail("R1 should exist")
  }
  assert_eq(@src.cobra3_model_objective_reaction(model).unwrap(), "R2")
}

///|
test "cobra3 core gpr parse and eval" {
  let gpr = @src.cobra3_gpr_parse("(g1 and g2) or g3").unwrap()
  let only_g3 : Map[String, Bool] = Map([("g3", true)])
  let both : Map[String, Bool] = Map([("g1", true), ("g2", true)])
  let only_g1 : Map[String, Bool] = Map([("g1", true)])
  assert_true(@src.cobra3_gpr_eval(gpr, only_g3))
  assert_true(@src.cobra3_gpr_eval(gpr, both))
  assert_true(!@src.cobra3_gpr_eval(gpr, only_g1))
  let genes = @src.cobra3_gpr_genes(gpr)
  assert_eq(genes.length(), 3)
  assert_eq(genes[0], "g1")
  assert_eq(genes[1], "g2")
  assert_eq(genes[2], "g3")
}

///|
test "cobra3 core gpr simplify and string" {
  let gpr = @src.cobra3_gpr_parse("g1 and (g2 and g3)").unwrap()
  let simple = @src.cobra3_gpr_simplify(gpr)
  assert_eq(@src.cobra3_gpr_to_string(simple), "(g1 and g2 and g3)")
}

///|
test "cobra3 core gpr knockout semantics" {
  let gpr = @src.cobra3_gpr_parse("g1 or (g2 and g3)").unwrap()
  let ko = @src.cobra3_gpr_knockout(gpr, "g1")
  let active : Map[String, Bool] = Map([("g2", true), ("g3", true)])
  assert_true(@src.cobra3_gpr_eval(ko, active))
  let ko2 = @src.cobra3_gpr_knockout(ko, "g2")
  let active2 : Map[String, Bool] = Map([("g3", true)])
  assert_true(!@src.cobra3_gpr_eval(ko2, active2))
}

///|
test "cobra3 core knockout reaction activity" {
  let r : @src.Cobra3Reaction = {
    id: "R",
    name: "R",
    metabolites: Map([], capacity=1),
    lower_bound: 0.0,
    upper_bound: 1000.0,
    objective: 0.0,
    subsystem: "",
    gpr: @src.cobra3_gpr_parse("A or B"),
    annotation: Map([], capacity=1),
  }
  assert_true(@src.cobra3_reaction_active_after_knockout(r, ["A"]))
  assert_true(!@src.cobra3_reaction_active_after_knockout(r, ["A", "B"]))
}

///|
test "cobra3 core add remove undo" {
  let model = cobra3_test_model()
  let model2 = @src.cobra3_model_remove_reaction(model, "R1")
  assert_eq(model2.reactions.length(), 2)
  assert_true(@src.cobra3_model_has_history(model2))
  let model3 = @src.cobra3_model_undo(model2)
  assert_eq(model3.reactions.length(), 3)
  assert_true(@src.cobra3_model_has_reaction(model3, "R1"))
}

///|
test "cobra3 core medium get set" {
  let model = cobra3_test_model()
  let medium = @src.cobra3_model_medium(model)
  assert_true((medium["EX_A"] - 10.0).abs() < 1.0e-9)
  let mg : Map[String, Double] = Map([("EX_A", 4.0)])
  let model2 = @src.cobra3_model_set_medium(model, mg)
  let medium2 = @src.cobra3_model_medium(model2)
  assert_true((medium2["EX_A"] - 4.0).abs() < 1.0e-9)
}

///|
test "cobra3 core formula and mass balance" {
  let counts = @src.cobra3_parse_formula("C6H12O6")
  assert_true((counts["C"] - 6.0).abs() < 1.0e-9)
  assert_true((counts["H"] - 12.0).abs() < 1.0e-9)
  assert_true((counts["O"] - 6.0).abs() < 1.0e-9)
  let counts2 = @src.cobra3_parse_formula("Ca(OH)2")
  assert_true((counts2["Ca"] - 1.0).abs() < 1.0e-9)
  assert_true((counts2["O"] - 2.0).abs() < 1.0e-9)
  assert_true((counts2["H"] - 2.0).abs() < 1.0e-9)
}

///|
test "cobra3 core mass balance balanced reaction" {
  let mut model = @src.cobra3_model_new("mb")
  let h2o : @src.Cobra3Metabolite = {
    id: "h2o_c",
    name: "h2o",
    formula: "H2O",
    charge: 0.0,
    compartment: "c",
    boundary: false,
    annotation: Map([], capacity=1),
  }
  let h : @src.Cobra3Metabolite = {
    id: "h_c",
    name: "h",
    formula: "H2",
    charge: 0.0,
    compartment: "c",
    boundary: false,
    annotation: Map([], capacity=1),
  }
  let o : @src.Cobra3Metabolite = {
    id: "o_c",
    name: "o",
    formula: "O2",
    charge: 0.0,
    compartment: "c",
    boundary: false,
    annotation: Map([], capacity=1),
  }
  model = @src.cobra3_model_add_metabolite(model, h2o)
  model = @src.cobra3_model_add_metabolite(model, h)
  model = @src.cobra3_model_add_metabolite(model, o)
  model = @src.cobra3_model_add_reaction_simple(
    model,
    "R",
    [("h2o_c", -1.0), ("h_c", 1.0), ("o_c", 0.5)],
    -1000.0,
    1000.0,
  )
  let mb = @src.cobra3_check_reaction_mass_balance(model, "R").unwrap()
  assert_true(mb.balanced)
}

///|
test "cobra3 core stoich matrix and fba conversion" {
  let model = cobra3_test_model()
  let s = @src.cobra3_model_stoichiometric_matrix(model)
  assert_eq(s.length(), 2)
  assert_eq(s[0].length(), 3)
  // A 行: EX_A +1, R1 -1
  assert_true((s[0][0] - 1.0).abs() < 1.0e-9)
  assert_true((s[0][1] + 1.0).abs() < 1.0e-9)
  let fba = @src.cobra3_to_fba_model(model)
  assert_eq(fba.reactions.length(), 3)
  let back = @src.cobra3_from_fba_model(fba)
  assert_eq(back.reactions.length(), 3)
  assert_eq(back.metabolites.length(), 2)
}

///|
test "cobra3 core exchange classification" {
  let model = cobra3_test_model()
  let ex = @src.cobra3_model_exchange_reactions(model)
  assert_eq(ex.length(), 1)
  assert_eq(ex[0], "EX_A")
}

///|
test "cobra3 core metadata preserved" {
  let met = @src.cobra3_metabolite_new("glc_c")
  assert_eq(met.compartment, "c")
  let gene = @src.cobra3_gene_new("g1")
  assert_true(gene.functional)
  let r = @src.cobra3_reaction_new("R")
  assert_eq(r.lower_bound, 0.0)
  assert_eq(r.upper_bound, 1000.0)
}

///|
test "cobra3 core set medium keeps upper bound" {
  let mut model = @src.cobra3_model_new("medium")
  model = @src.cobra3_model_add_reaction_simple(
    model,
    "EX_A",
    [("A", -1.0)],
    -10.0,
    500.0,
  )
  model = @src.cobra3_model_add_reaction_simple(
    model,
    "EX_B",
    [("B", -1.0)],
    -5.0,
    300.0,
  )
  let mg : Map[String, Double] = Map([("EX_A", 4.0)])
  let next = @src.cobra3_model_set_medium(model, mg)
  let a = @src.cobra3_model_get_reaction(next, "EX_A").unwrap()
  assert_eq(a.upper_bound, 500.0)
  assert_true((a.lower_bound + 4.0).abs() < 1.0e-9)
  let b = @src.cobra3_model_get_reaction(next, "EX_B").unwrap()
  assert_eq(b.upper_bound, 300.0)
  assert_eq(b.lower_bound, 0.0)
}

///|
test "cobra3 core add reaction undo is atomic" {
  let mut model = @src.cobra3_model_new("undo")
  let participants : Array[(String, Double)] = [("A", -1.0), ("B", 1.0)]
  let r : @src.Cobra3Reaction = {
    id: "R",
    name: "R",
    metabolites: Map(participants),
    lower_bound: 0.0,
    upper_bound: 1000.0,
    objective: 0.0,
    subsystem: "",
    gpr: @src.cobra3_gpr_parse("g1 and g2"),
    annotation: Map([], capacity=1),
  }
  model = @src.cobra3_model_add_reaction(model, r)
  assert_eq(model.metabolites.length(), 2)
  assert_eq(model.genes.length(), 2)
  let back = @src.cobra3_model_undo(model)
  assert_eq(back.reactions.length(), 0)
  assert_eq(back.metabolites.length(), 0)
  assert_eq(back.genes.length(), 0)
}

///|
test "cobra3 core gpr constants roundtrip" {
  let empty : Map[String, Bool] = Map([])
  let f = @src.cobra3_gpr_parse("false").unwrap()
  let t = @src.cobra3_gpr_parse("true").unwrap()
  assert_true(!@src.cobra3_gpr_eval(f, empty))
  assert_true(@src.cobra3_gpr_eval(t, empty))
  assert_eq(@src.cobra3_gpr_to_string(f), "false")
  assert_eq(@src.cobra3_gpr_to_string(t), "true")
}

///|
test "cobra3 core nonfunctional gene disables reaction" {
  let mut model = @src.cobra3_model_new("functional")
  model = @src.cobra3_model_add_gene(model, {
    id: "g1",
    name: "g1",
    annotation: Map([], capacity=1),
    functional: false,
  })
  let participants : Array[(String, Double)] = [("A", -1.0)]
  let r : @src.Cobra3Reaction = {
    id: "R",
    name: "R",
    metabolites: Map(participants),
    lower_bound: 0.0,
    upper_bound: 1000.0,
    objective: 0.0,
    subsystem: "",
    gpr: @src.cobra3_gpr_parse("g1"),
    annotation: Map([], capacity=1),
  }
  assert_true(!@src.cobra3_reaction_active_in_model(model, r, []))
  assert_true(@src.cobra3_reaction_active_after_knockout(r, []))
}

///|
test "cobra3 core copy deep copies maps" {
  let model = cobra3_test_model()
  let copy = @src.cobra3_model_copy(model)
  let r = @src.cobra3_model_get_reaction(copy, "R1").unwrap()
  r.metabolites["A"] = -99.0
  let orig = @src.cobra3_model_get_reaction(model, "R1").unwrap()
  assert_true((orig.metabolites["A"] + 1.0).abs() < 1.0e-9)
}
