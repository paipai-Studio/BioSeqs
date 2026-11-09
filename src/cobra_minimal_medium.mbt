/// cobra_minimal_medium.mbt — Minimal medium search (COBRApy medium)
/// Find the smallest set of import fluxes that still enables growth.

///|
/// Minimal medium result.
pub struct CobraMinimalMediumResult {
  medium : Map[String, Double]
  objective_value : Double
  n_components : Int32
}

///|
/// Find minimal medium: minimize total import flux while maintaining growth.
///
/// Parameters:
/// - min_objective: minimum growth rate required
/// - max_import: maximum import flux per component
pub fn cobra_minimal_medium(
  model : FbaModel,
  min_objective~ : Double = 0.1,
  _max_import~ : Double = 1000.0,
  extracellular~ : String = "e",
) -> CobraMinimalMediumResult {
  let _ = extracellular
  let m3 = cobra3_from_fba_model(model)
  let res = cobra3_minimal_medium(
    m3,
    min_growth=min_objective,
    max_import=_max_import,
  )
  {
    medium: res.medium,
    objective_value: res.objective_value,
    n_components: res.n_components,
  }
}
