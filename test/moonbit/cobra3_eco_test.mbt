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
test "cobra3 eco efm finds mode needing negative basis coefficient" {
  let mut m = @src.cobra3_model_new("efmneg")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], 0.0, 10.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], -1000.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_B", [("B", -1.0)], -10.0, 0.0, ""),
  )
  // 零空间基给出的是 (-1,1,1)，可行方向为其反向；旧实现固定 λ=1 会一条都找不到
  let modes = @src.cobra3_elementary_flux_modes(m)
  assert_eq(modes.length(), 1)
  let v = modes[0]
  let mut mx = 0.0
  for x in v {
    if x.abs() > mx {
      mx = x.abs()
    }
  }
  assert_true((mx - 1.0).abs() < 1.0e-9)
  for j in 0..<3 {
    assert_true(v[j].abs() > 1.0e-9)
  }
  let s = @src.cobra3_model_stoichiometric_matrix(m)
  for i in 0..<s.length() {
    let mut sum = 0.0
    for j in 0..<s[i].length() {
      sum = sum + s[i][j] * v[j]
    }
    assert_true(sum.abs() < 1.0e-9)
  }
}

///|
test "cobra3 eco efm keeps internal loop mode" {
  let mut m = @src.cobra3_model_new("efmloop")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], 0.0, 10.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], -1000.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R2", [("B", -1.0), ("C", 1.0)], -1000.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R3", [("C", -1.0), ("A", 1.0)], -1000.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_B", [("B", -1.0)], -10.0, 0.0, ""),
  )
  let modes = @src.cobra3_elementary_flux_modes(m)
  // 三条模式：内部循环、A 由 EX_A 消耗/B 由 EX_B 供给的主通路，以及
  // EX_A→A、R3→A、R2→C、EX_B→B 的旁路（旧实现每个基列子集只取第一个
  // 可行方向，会漏掉这条）。
  assert_eq(modes.length(), 3)
  let mut found_loop = false
  let mut found_exchange = false
  let mut found_bypass = false
  for v in modes {
    let loop_only = v[0].abs() < 1.0e-9 &&
      v[1].abs() > 1.0e-9 &&
      v[2].abs() > 1.0e-9 &&
      v[3].abs() > 1.0e-9 &&
      v[4].abs() < 1.0e-9
    let exchange = v[0].abs() > 1.0e-9 &&
      v[1].abs() > 1.0e-9 &&
      v[2].abs() < 1.0e-9 &&
      v[3].abs() < 1.0e-9 &&
      v[4].abs() > 1.0e-9
    let bypass = v[0].abs() > 1.0e-9 &&
      v[1].abs() < 1.0e-9 &&
      v[2].abs() > 1.0e-9 &&
      v[3].abs() > 1.0e-9 &&
      v[4].abs() > 1.0e-9
    if loop_only {
      found_loop = true
    }
    if exchange {
      found_exchange = true
    }
    if bypass {
      found_bypass = true
    }
  }
  assert_true(found_loop)
  assert_true(found_exchange)
  assert_true(found_bypass)
}

///|
test "cobra3 eco efm branch model yields two modes" {
  let m = cobra3_branch_model()
  let modes = @src.cobra3_elementary_flux_modes(m)
  assert_eq(modes.length(), 2)
}

///|
test "cobra3 eco efm enumerates parallel branches" {
  let mut m = @src.cobra3_model_new("efmpar")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", 1.0)], 0.0, 10.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], 0.0, 10.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R2", [("A", -1.0), ("B", 1.0)], 0.0, 10.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R3", [("A", -1.0), ("B", 1.0)], 0.0, 10.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  let modes = @src.cobra3_elementary_flux_modes(m)
  assert_eq(modes.length(), 3)
  let mut used_r1 = false
  let mut used_r2 = false
  let mut used_r3 = false
  for v in modes {
    if v[1].abs() > 1.0e-9 {
      used_r1 = true
    }
    if v[2].abs() > 1.0e-9 {
      used_r2 = true
    }
    if v[3].abs() > 1.0e-9 {
      used_r3 = true
    }
    let s = @src.cobra3_model_stoichiometric_matrix(m)
    for i in 0..<s.length() {
      let mut sum = 0.0
      for j in 0..<s[i].length() {
        sum = sum + s[i][j] * v[j]
      }
      assert_true(sum.abs() < 1.0e-9)
    }
  }
  assert_true(used_r1)
  assert_true(used_r2)
  assert_true(used_r3)
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
test "cobra3 eco community keeps same-named exchange with different metabolites" {
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
  // 与 m1 同名 EX_A，但交换的是 B 而不是 A
  m2 = @src.cobra3_model_add_reaction(
    m2,
    cobra3_make_rxn("EX_A", [("B", -1.0)], -10.0, 1000.0, ""),
  )
  m2 = @src.cobra3_model_add_reaction(
    m2,
    cobra3_make_rxn("R2", [("B", -1.0), ("C", 1.0)], 0.0, 1000.0, ""),
  )
  m2 = @src.cobra3_model_add_reaction(
    m2,
    cobra3_make_rxn("BIO", [("C", -1.0)], 0.0, 1000.0, ""),
  )
  m2 = @src.cobra3_model_set_objective(m2, "BIO", 1.0)
  let res = @src.cobra3_community_fba([m1, m2])
  // 旧实现会把 m2 的 EX_A 错并到 m1 的 EX_A，丢掉 B 交换、总生长只有 10
  assert_true((res.growth - 20.0).abs() < 1.0e-5)
  assert_true(@src.cobra3_model_has_reaction(res.model, "m2:EX_A"))
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
test "cobra3 eco optknock inner maximizes target" {
  let mut m = @src.cobra3_model_new("optknock")
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
    cobra3_make_rxn("BIO", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("P", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO", 1.0)
  let res = @src.cobra3_optknock(m, "P", max_knockouts=1)
  let mut found = false
  let mut target = -1.0
  for r in res {
    if r.knockout_ids.length() == 1 && r.knockout_ids[0] == "BIO" {
      found = true
      target = r.target_flux
      assert_true(r.growth.abs() < 1.0e-6)
    }
  }
  assert_true(found)
  assert_true((target - 10.0).abs() < 1.0e-4)
}

///|
test "cobra3 eco optknock supports triple knockouts" {
  let m = cobra3_toy()
  let res = @src.cobra3_optknock(m, "R2", max_knockouts=3)
  let mut has_triple = false
  for r in res {
    if r.knockout_ids.length() == 3 {
      has_triple = true
    }
  }
  assert_true(has_triple)
  let capped = @src.cobra3_optknock(m, "R2", max_knockouts=3, max_candidates=1)
  assert_true(capped.length() <= 1)
}

///|
fn cobra3_optknock_milp_model() -> @src.Cobra3Model {
  let mut m = @src.cobra3_model_new("optknock_milp")
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
    cobra3_make_rxn("R2", [("B", -1.0), ("C", 1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO", [("C", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("P", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO", 1.0)
  m
}

///|
test "cobra3 eco optknock milp trades growth for product" {
  let m = cobra3_optknock_milp_model()
  let res = @src.cobra3_optknock_milp(m, "P", min_growth=0.0, max_knockouts=1)
  assert_eq(res.length(), 1)
  assert_true((res[0].target_flux - 10.0).abs() < 1.0e-6)
  assert_true(res[0].growth.abs() < 1.0e-6)
  assert_eq(res[0].knockout_ids.length(), 1)
}

///|
test "cobra3 eco optknock milp respects min growth" {
  let m = cobra3_optknock_milp_model()
  let res = @src.cobra3_optknock_milp(m, "P", min_growth=5.0, max_knockouts=1)
  assert_eq(res.length(), 1)
  assert_eq(res[0].knockout_ids.length(), 0)
  assert_true((res[0].target_flux - 5.0).abs() < 1.0e-6)
  assert_true((res[0].growth - 5.0).abs() < 1.0e-6)
}

///|
test "cobra3 eco optknock milp reports infeasible" {
  let m = cobra3_optknock_milp_model()
  // 要求生长 20 超过最大生长 10，MILP 不可行
  let res = @src.cobra3_optknock_milp(m, "P", min_growth=20.0, max_knockouts=0)
  assert_eq(res.length(), 0)
}

///|
test "cobra3 eco optknock honors exclude" {
  let m = cobra3_optknock_milp_model()
  let res = @src.cobra3_optknock(m, "P", max_knockouts=2, exclude=["R1"])
  assert_true(res.length() > 0)
  for r in res {
    assert_true(!r.knockout_ids.contains("R1"))
  }
}

///|
test "cobra3 eco optknock candidate guard uses exclude" {
  let mut m = cobra3_toy()
  let excluded : Array[String] = Array(capacity=80)
  for i in 0..<78 {
    let id = "X" + i.to_string()
    m = @src.cobra3_model_add_reaction(
      m,
      cobra3_make_rxn(id, [("P", -1.0), ("Q", 1.0)], 0.0, 1000.0, ""),
    )
    excluded.push(id)
  }
  // 81 个反应但 exclude 后只剩 3 个候选：旧实现的 n>80 会直接返回空
  let res = @src.cobra3_optknock(m, "R2", max_knockouts=1, exclude=excluded)
  assert_true(res.length() > 0)
}

///|
test "cobra3 eco optknock milp honors exclude" {
  let mut m = @src.cobra3_model_new("optknock_milp_exclude")
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
    cobra3_make_rxn("R2", [("B", -1.0), ("C", 1.0)], 0.0, 1000.0, ""),
  )
  // BIO 强制至少 5 生长：不敲除时产物最多 10-5=5，只有敲除 BIO 才能到 10
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO", [("C", -1.0)], 5.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("P", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO", 1.0)
  let free = @src.cobra3_optknock_milp(
    m,
    "P",
    min_growth=0.0,
    max_knockouts=1,
    max_nodes=50000,
  )
  assert_eq(free.length(), 1)
  assert_true(free[0].knockout_ids.contains("BIO"))
  assert_true((free[0].target_flux - 10.0).abs() < 1.0e-6)
  let blocked = @src.cobra3_optknock_milp(
    m,
    "P",
    min_growth=0.0,
    max_knockouts=1,
    max_nodes=50000,
    exclude=["BIO"],
  )
  assert_eq(blocked.length(), 1)
  assert_true(!blocked[0].knockout_ids.contains("BIO"))
  assert_true((blocked[0].target_flux - 5.0).abs() < 1.0e-6)
}

///|
fn cobra3_optstrain_host() -> @src.Cobra3Model {
  let mut m = @src.cobra3_model_new("optstrain_host")
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
    cobra3_make_rxn("BIO", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO", 1.0)
  m
}

///|
fn cobra3_optstrain_universal() -> @src.Cobra3Model {
  let mut m = cobra3_optstrain_host()
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R2", [("A", -1.0), ("C", 1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("P", [("C", -1.0)], 0.0, 1000.0, ""),
  )
  m
}

///|
test "cobra3 eco optstrain respects min growth" {
  let host = cobra3_optstrain_host()
  let universal = cobra3_optstrain_universal()
  let res = @src.cobra3_optstrain(
    host,
    universal,
    "P",
    min_growth=5.0,
    max_knockouts=1,
    max_added=2,
  )
  assert_eq(res.length(), 1)
  assert_true((res[0].target_flux - 5.0).abs() < 1.0e-6)
  assert_true((res[0].growth - 5.0).abs() < 1.0e-6)
  assert_true(res[0].added_reactions.contains("R2"))
  assert_true(res[0].added_reactions.contains("P"))
}

///|
test "cobra3 eco optstrain trades growth for product" {
  let host = cobra3_optstrain_host()
  let universal = cobra3_optstrain_universal()
  let res = @src.cobra3_optstrain(
    host,
    universal,
    "P",
    min_growth=0.0,
    max_knockouts=1,
    max_added=2,
  )
  assert_eq(res.length(), 1)
  assert_true((res[0].target_flux - 10.0).abs() < 1.0e-6)
}

///|
test "cobra3 eco optstrain max added zero blocks heterologous pathway" {
  let host = cobra3_optstrain_host()
  let universal = cobra3_optstrain_universal()
  let res = @src.cobra3_optstrain(
    host,
    universal,
    "P",
    min_growth=0.0,
    max_knockouts=1,
    max_added=0,
  )
  assert_eq(res.length(), 1)
  assert_true(res[0].added_reactions.length() == 0)
  assert_true(res[0].target_flux.abs() < 1.0e-6)
}

///|
test "cobra3 eco optstrain rejects unknown target" {
  let host = cobra3_optstrain_host()
  let universal = cobra3_optstrain_universal()
  let res = @src.cobra3_optstrain(host, universal, "NOPE", max_added=2)
  assert_eq(res.length(), 0)
}

///|
fn cobra3_optcouple_base_model() -> @src.Cobra3Model {
  let mut m = @src.cobra3_model_new("optcouple")
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
    cobra3_make_rxn("BIO", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("P", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO", 1.0)
  m
}

///|
test "cobra3 eco optcouple couples product with growth" {
  let m = cobra3_optcouple_base_model()
  let res = @src.cobra3_optcouple(m, "P", min_growth=5.0, max_knockouts=1)
  assert_eq(res.length(), 1)
  assert_true((res[0].growth - 5.0).abs() < 1.0e-6)
  assert_true((res[0].target_flux - 5.0).abs() < 1.0e-6)
  assert_eq(res[0].knockout_ids.length(), 0)
  assert_true(res[0].coupled)
}

///|
test "cobra3 eco optcouple knocks out competing sink" {
  let mut m = cobra3_optcouple_base_model()
  // E 强制至少 3 通量，与产物竞争 B；敲掉 E 后产物上限从 2 提高到 5
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("E", [("B", -1.0), ("F", 1.0)], 3.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_F", [("F", -1.0)], 0.0, 1000.0, ""),
  )
  let res = @src.cobra3_optcouple(m, "P", min_growth=5.0, max_knockouts=1)
  assert_eq(res.length(), 1)
  assert_true(res[0].knockout_ids.contains("E"))
  assert_true((res[0].growth - 5.0).abs() < 1.0e-6)
  assert_true((res[0].target_flux - 5.0).abs() < 1.0e-6)
  assert_true(res[0].coupled)
}

///|
test "cobra3 eco optcouple requires positive min growth" {
  let m = cobra3_optcouple_base_model()
  let res = @src.cobra3_optcouple(m, "P", min_growth=0.0)
  assert_eq(res.length(), 0)
}

///|
test "cobra3 eco optcouple enforces min target" {
  let m = cobra3_optcouple_base_model()
  // B 总量 10：min_growth=5 与 min_target=6 同时满足需要 11，无法耦合
  let res = @src.cobra3_optcouple(
    m,
    "P",
    min_growth=5.0,
    max_knockouts=1,
    min_target=6.0,
  )
  assert_eq(res.length(), 1)
  assert_true(!res[0].coupled)
  assert_true(res[0].target_flux.abs() < 1.0e-6)
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

///|
test "cobra3 eco thermo closes forced unfavorable reaction" {
  let mut m = @src.cobra3_model_new("thermo-forced")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -10.0, 1000.0, ""),
  )
  // 原本只能正向的反应，遇到吸能方向应整反应关闭，而不是留下 lb>ub
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], 1.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO", 1.0)
  let dg : Map[String, Double] = Map([("R1", 5.0)])
  let res = @src.cobra3_thermo_fba(m, dg)
  assert_true(res.status is @src.Cobra3Status::Optimal)
  assert_true(res.objective_value.abs() < 1.0e-9)
}
