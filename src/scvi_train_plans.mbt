// scvi_train_plans.mbt - Training plans mirroring scvi.train
//
// scvi-tools exposes dedicated Lightning plans for unsupervised, semi-
// supervised, classifier and adversarial training. MoonBit has no Lightning
// runtime, so this module provides the plan configuration, composite losses
// and a reusable training loop with early stopping.

///|
pub struct ScviTrainingPlanConfig {
  plan_type : String
  n_epochs : Int32
  learning_rate : Double
  early_stopping : Bool
  patience : Int32
  min_delta : Double
}

///|
pub fn ScviTrainingPlanConfig::new(
  plan_type : String,
  n_epochs : Int32,
  learning_rate : Double,
  early_stopping : Bool,
  patience : Int32,
  min_delta : Double,
) -> ScviTrainingPlanConfig {
  { plan_type, n_epochs, learning_rate, early_stopping, patience, min_delta, }
}

///|
pub fn ScviTrainingPlanConfig::default() -> ScviTrainingPlanConfig {
  ScviTrainingPlanConfig::new("unsupervised", 100, 0.001, true, 10, 1.0e-4)
}

///|
pub fn ScviTrainingPlanConfig::unsupervised(
  n_epochs : Int32,
  learning_rate : Double,
) -> ScviTrainingPlanConfig {
  ScviTrainingPlanConfig::new(
    "unsupervised", n_epochs, learning_rate, true, 10, 1.0e-4,
  )
}

///|
pub fn ScviTrainingPlanConfig::semisupervised(
  n_epochs : Int32,
  learning_rate : Double,
) -> ScviTrainingPlanConfig {
  ScviTrainingPlanConfig::new(
    "semisupervised", n_epochs, learning_rate, true, 10, 1.0e-4,
  )
}

///|
pub fn ScviTrainingPlanConfig::classifier(
  n_epochs : Int32,
  learning_rate : Double,
) -> ScviTrainingPlanConfig {
  ScviTrainingPlanConfig::new(
    "classifier", n_epochs, learning_rate, true, 10, 1.0e-4,
  )
}

///|
pub fn ScviTrainingPlanConfig::adversarial(
  n_epochs : Int32,
  learning_rate : Double,
) -> ScviTrainingPlanConfig {
  ScviTrainingPlanConfig::new(
    "adversarial", n_epochs, learning_rate, true, 10, 1.0e-4,
  )
}

///|
pub struct ScviPlanResult {
  epochs_run : Int32
  best_loss : Double
  final_loss : Double
  stopped_early : Bool
  history : Array[Double]
}

///|
/// Run a generic training loop. `step(epoch)` must perform one optimisation
/// step and return the training loss for that epoch.
pub fn scvi_training_plan_run(
  config : ScviTrainingPlanConfig,
  step : (Int32) -> Double,
) -> ScviPlanResult {
  let history = Array(capacity=config.n_epochs)
  let mut best = 1.0e30
  let mut best_epoch = 0
  let mut bad_epochs = 0
  let mut last = 0.0
  let mut stopped = false
  for epoch in 0..<config.n_epochs {
    let loss = step(epoch)
    history.push(loss)
    last = loss
    if loss < best - config.min_delta {
      best = loss
      best_epoch = epoch
      bad_epochs = 0
    } else {
      bad_epochs = bad_epochs + 1
    }
    if config.early_stopping && bad_epochs >= config.patience {
      stopped = true
      break
    }
  }
  ignore(best_epoch)
  {
    epochs_run: history.length(),
    best_loss: best,
    final_loss: last,
    stopped_early: stopped,
    history,
  }
}

///|
/// Composite semi-supervised loss: labelled cross-entropy plus an entropy
/// penalty on unlabelled cells.
pub fn scvi_semisupervised_loss(
  cross_entropy : Double,
  entropy : Double,
  alpha : Double,
) -> Double {
  cross_entropy + alpha * entropy
}

///|
/// Cross-entropy for a probability vector and an integer label.
pub fn scvi_classifier_loss(probs : Array[Double], label : Int32) -> Double {
  if probs.length() == 0 {
    return 0.0
  }
  let idx = if label < 0 {
    0
  } else if label >= probs.length() {
    probs.length() - 1
  } else {
    label
  }
  -@math.ln(scvi_dist_clamp(probs[idx]))
}

///|
/// Binary cross-entropy for adversarial domain prediction.
pub fn scvi_adversarial_loss(
  domain_prob : Double,
  domain_label : Int32,
) -> Double {
  let p = scvi_dist_clamp(domain_prob)
  let y = domain_label.to_double()
  -(y * @math.ln(p) + (1.0 - y) * @math.ln(1.0 - p))
}

///|
/// Entropy of a categorical distribution.
pub fn scvi_entropy(probs : Array[Double]) -> Double {
  let mut h = 0.0
  for i in 0..<probs.length() {
    let p = probs[i]
    if p > 1.0e-8 {
      h += -p * @math.ln(p)
    }
  }
  h
}

///|
/// Reduce a result's history to its most recent `n` losses.
pub fn scvi_plan_recent_losses(
  result : ScviPlanResult,
  n : Int32,
) -> Array[Double] {
  let total = result.history.length()
  let start = if n >= total { 0 } else { total - n }
  let out = Array(capacity=total - start)
  for i in start..<total {
    out.push(result.history[i])
  }
  out
}
