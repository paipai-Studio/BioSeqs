// scvi_autotune.mbt - Hyper-parameter search mirroring scvi.autotune
//
// scvi-tools uses Ray Tune / Optuna through ModelTuner. This implementation
// keeps the same surface (grid and random search over a parameter space) while
// staying dependency-free: callers supply an evaluation callback.

///|
pub struct ScviTuneConfig {
  n_trials : Int32
  seed : Int32
  metric : String
  greater_is_better : Bool
}

///|
pub fn ScviTuneConfig::new(
  n_trials : Int32,
  seed : Int32,
  metric : String,
  greater_is_better : Bool,
) -> ScviTuneConfig {
  { n_trials, seed, metric, greater_is_better, }
}

///|
pub fn ScviTuneConfig::default() -> ScviTuneConfig {
  { n_trials: 10, seed: 0, metric: "loss", greater_is_better: false, }
}

///|
pub struct ScviTuneTrial {
  params : Array[(String, Double)]
  score : Double
}

///|
pub struct ScviTuneResult {
  best_params : Array[(String, Double)]
  best_score : Double
  trials : Array[ScviTuneTrial]
}

///|
fn scvi_autotune_is_better(
  score : Double,
  best : Double,
  greater_is_better : Bool,
) -> Bool {
  if greater_is_better {
    score > best
  } else {
    score < best
  }
}

///|
fn scvi_autotune_initial_score(greater_is_better : Bool) -> Double {
  if greater_is_better {
    -1.0e30
  } else {
    1.0e30
  }
}

///|
/// Random search over `(name, low, high)` ranges.
pub fn scvi_autotune_random(
  ranges : Array[(String, Double, Double)],
  config : ScviTuneConfig,
  evaluate : (Array[(String, Double)]) -> Double,
) -> ScviTuneResult {
  let trials = Array(capacity=config.n_trials)
  let mut best_score = scvi_autotune_initial_score(config.greater_is_better)
  let mut best_params = Array(capacity=ranges.length())
  for i in 0..<ranges.length() {
    let (name, low, _high) = ranges[i]
    best_params.push((name, low))
  }
  for trial in 0..<config.n_trials {
    let params = Array(capacity=ranges.length())
    for i in 0..<ranges.length() {
      let (name, low, high) = ranges[i]
      let u = scvi_dist_uniform(config.seed + trial * 104729 + i * 7919)
      params.push((name, low + (high - low) * u))
    }
    let score = evaluate(params)
    trials.push({ params, score, })
    if scvi_autotune_is_better(score, best_score, config.greater_is_better) {
      best_score = score
      best_params = params
    }
  }
  { best_params, best_score, trials, }
}

///|
/// Exhaustive grid search. Each entry is `(name, candidate_values)`.
pub fn scvi_autotune_grid(
  grid : Array[(String, Array[Double])],
  greater_is_better : Bool,
  evaluate : (Array[(String, Double)]) -> Double,
) -> ScviTuneResult {
  let dims = grid.length()
  let trials = Array(capacity=0)
  let mut best_score = scvi_autotune_initial_score(greater_is_better)
  let mut best_params = Array(capacity=0)
  if dims == 0 {
    return { best_params, best_score, trials, }
  }
  let sizes = Array(capacity=dims)
  let mut total = 1
  for d in 0..<dims {
    let (_name, values) = grid[d]
    let size = values.length()
    sizes.push(size)
    total = total * size
  }
  if total == 0 {
    return { best_params, best_score, trials, }
  }
  for combo in 0..<total {
    let params = Array(capacity=dims)
    let mut rem = combo
    for d in 0..<dims {
      let (name, values) = grid[d]
      let idx = rem % sizes[d]
      rem = rem / sizes[d]
      params.push((name, values[idx]))
    }
    let score = evaluate(params)
    trials.push({ params, score, })
    if scvi_autotune_is_better(score, best_score, greater_is_better) {
      best_score = score
      best_params = params
    }
  }
  { best_params, best_score, trials, }
}

///|
/// scvi-tools `ModelTuner.tune` equivalent (random search).
pub fn scvi_autotune_model_tuner(
  ranges : Array[(String, Double, Double)],
  config : ScviTuneConfig,
  evaluate : (Array[(String, Double)]) -> Double,
) -> ScviTuneResult {
  scvi_autotune_random(ranges, config, evaluate)
}

///|
/// Extract one parameter value from a trial parameter list.
pub fn scvi_autotune_param(
  params : Array[(String, Double)],
  name : String,
) -> Double? {
  for i in 0..<params.length() {
    let (k, v) = params[i]
    if k == name {
      return Some(v)
    }
  }
  None
}
