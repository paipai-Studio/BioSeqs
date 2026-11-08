///|
/// COBRApy v3 演示：对象模型 / FBA / pFBA / FVA / 敲除 / SBML & JSON.
fn main {
  println("=== BioSeqs COBRApy v3 (cobra3_*) ===")

  // 构建 toy 模型：葡萄糖摄取 -> F6P -> biomass
  let mut m = @src.cobra3_model_new("demo")
  m = @src.cobra3_model_add_reaction_simple(
    m,
    "EX_glc",
    [("glc_c", -1.0)],
    -10.0,
    1000.0,
  )
  m = @src.cobra3_model_add_reaction_simple(
    m,
    "PGI",
    [("glc_c", -1.0), ("f6p_c", 1.0)],
    0.0,
    1000.0,
  )
  m = @src.cobra3_model_add_reaction_simple(
    m,
    "BIOMASS",
    [("f6p_c", -1.0)],
    0.0,
    1000.0,
  )
  m = @src.cobra3_model_set_objective(m, "BIOMASS", 1.0)

  let sol = @src.cobra3_model_optimize(m)
  println("growth = " + sol.objective_value.to_string())
  let fmap = @src.cobra3_flux_map(m, sol)
  for id in ["EX_glc", "PGI", "BIOMASS"] {
    println("  flux " + id + " = " + fmap[id].to_string())
  }

  let pf = @src.cobra3_pfba(m)
  let mut total = 0.0
  for v in pf.fluxes {
    total = total + v.abs()
  }
  println("pFBA total |flux| = " + total.to_string())

  println("FVA (fraction=1.0):")
  for f in @src.cobra3_fva(m, fraction_of_optimum=1.0) {
    println(
      "  " +
      f.reaction_id +
      " [" +
      f.minimum.to_string() +
      ", " +
      f.maximum.to_string() +
      "]",
    )
  }

  println("single reaction deletion:")
  for d in @src.cobra3_single_reaction_deletion(m) {
    println("  " + d.id + " -> growth " + d.growth.to_string())
  }

  let xml = @src.cobra3_model_to_sbml(m)
  println("SBML bytes = " + xml.length().to_string())
  match @src.cobra3_read_sbml(xml) {
    Some(back) =>
      println(
        "SBML roundtrip reactions = " + back.reactions.length().to_string(),
      )
    None => println("SBML parse failed")
  }

  let json = @src.cobra3_model_to_json(m)
  match @src.cobra3_model_from_json(json) {
    Some(back) =>
      println(
        "JSON roundtrip reactions = " + back.reactions.length().to_string(),
      )
    None => println("JSON parse failed")
  }

  let medium = @src.cobra3_minimal_medium(m, min_growth=5.0, max_import=10.0)
  println("minimal medium components = " + medium.n_components.to_string())
}
