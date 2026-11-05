///|
/// Training infrastructure: TrainRunner, EarlyStopping, LRScheduler,
/// ModelCheckpoint, get_max_epochs_heuristic, TrainingHistory,
/// AdversarialClassifier, PerLabelSubsampling.
pub struct EarlyStoppingConfig {
  monitor : String
  patience : Int32
  min_delta : Double
  mode : String
}

///|
pub fn EarlyStoppingConfig::default() -> EarlyStoppingConfig {
  { monitor: "elbo", patience: 10, min_delta: 0.0, mode: "max", }
}

///|
pub struct EarlyStoppingState {
  config : EarlyStoppingConfig
  mut best_score : Double
  mut counter : Int32
  mut should_stop : Bool
}

///|
pub fn early_stopping_init(config : EarlyStoppingConfig) -> EarlyStoppingState {
  let init_score = match config.mode {
    "max" => -999999999999.0
    _ => 999999999999.0
  }
  { config, best_score: init_score, counter: 0, should_stop: false, }
}

///|
pub fn early_stopping_step(
  state : EarlyStoppingState,
  score : Double,
) -> EarlyStoppingState {
  let s = state
  let improved = match s.config.mode {
    "max" => score > s.best_score + s.config.min_delta
    _ => score < s.best_score - s.config.min_delta
  }
  if improved {
    s.best_score = score
    s.counter = 0
  } else {
    s.counter = s.counter + 1
    if s.counter >= s.config.patience {
      s.should_stop = true
    }
  }
  s
}

///|
pub enum LRScheduleType {
  ReduceLROnPlateau
  CosineAnnealing
  StepLR
  ExponentialLR
}

///|
pub struct LRScheduler {
  schedule_type : LRScheduleType
  initial_lr : Double
  mut current_lr : Double
  factor : Double
  patience : Int32
  min_lr : Double
  mut step_count : Int32
  mut best_loss : Double
  mut bad_count : Int32
}

///|
pub fn LRScheduler::new(
  schedule_type : LRScheduleType,
  initial_lr : Double,
) -> LRScheduler {
  {
    schedule_type,
    initial_lr,
    current_lr: initial_lr,
    factor: 0.5,
    patience: 5,
    min_lr: 0.000001,
    step_count: 0,
    best_loss: 999999999999.0,
    bad_count: 0,
  }
}

///|
pub fn lr_scheduler_step(scheduler : LRScheduler, loss : Double) -> LRScheduler {
  let s = scheduler
  s.step_count = s.step_count + 1
  match s.schedule_type {
    ReduceLROnPlateau =>
      if loss < s.best_loss {
        s.best_loss = loss
        s.bad_count = 0
      } else {
        s.bad_count = s.bad_count + 1
        if s.bad_count >= s.patience {
          s.current_lr = s.current_lr * s.factor
          if s.current_lr < s.min_lr {
            s.current_lr = s.min_lr
          }
          s.bad_count = 0
        }
      }
    CosineAnnealing => {
      let total_steps = 100
      s.current_lr = s.min_lr +
        (s.initial_lr - s.min_lr) *
        (
          1.0 +
          @math.cos(
            3.141592653589793 *
            s.step_count.to_double() /
            total_steps.to_double(),
          )
        ) /
        2.0
    }
    StepLR =>
      if s.step_count % 10 == 0 {
        s.current_lr = s.current_lr * s.factor
        if s.current_lr < s.min_lr {
          s.current_lr = s.min_lr
        }
      }
    ExponentialLR => {
      s.current_lr = s.initial_lr *
        @math.pow(s.factor, s.step_count.to_double())
      if s.current_lr < s.min_lr {
        s.current_lr = s.min_lr
      }
    }
  }
  s
}

///|
pub struct ModelCheckpoint {
  save_best : Bool
  mut best_score : Double
  save_path : String
  monitor : String
  mode : String
}

///|
pub fn ModelCheckpoint::new(save_path : String) -> ModelCheckpoint {
  {
    save_best: true,
    best_score: -999999999999.0,
    save_path,
    monitor: "elbo",
    mode: "max",
  }
}

///|
pub fn model_checkpoint_should_save(
  ckpt : ModelCheckpoint,
  score : Double,
) -> (ModelCheckpoint, Bool) {
  let c = ckpt
  let improved = match c.mode {
    "max" => score > c.best_score
    _ => score < c.best_score
  }
  if improved {
    c.best_score = score
    (c, true)
  } else {
    (c, false)
  }
}

///|
pub fn get_max_epochs_heuristic(n_cells : Int32, batch_size : Int32) -> Int32 {
  let n_batches = n_cells / batch_size
  if n_batches < 100 {
    400
  } else if n_batches < 1000 {
    200
  } else if n_batches < 10000 {
    100
  } else {
    50
  }
}

///|
pub struct TrainingHistory {
  loss_values : Array[Double]
  val_values : Array[Double]
  metric_values : Map[String, Array[Double]]
  epoch_times : Array[Double]
}

///|
pub fn TrainingHistory::new() -> TrainingHistory {
  { loss_values: [], val_values: [], metric_values: Map([]), epoch_times: [], }
}

///|
pub fn training_history_record_loss(
  history : TrainingHistory,
  loss : Double,
) -> TrainingHistory {
  let h = history
  h.loss_values.push(loss)
  h
}

///|
pub fn training_history_record_val(
  history : TrainingHistory,
  val : Double,
) -> TrainingHistory {
  let h = history
  h.val_values.push(val)
  h
}

///|
pub fn training_history_record_metric(
  history : TrainingHistory,
  name : String,
  value : Double,
) -> TrainingHistory {
  let h = history
  let existing = h.metric_values.get(name).unwrap_or([])
  let arr = existing
  arr.push(value)
  h.metric_values[name] = arr
  h
}

///|
pub fn training_history_record_time(
  history : TrainingHistory,
  time : Double,
) -> TrainingHistory {
  let h = history
  h.epoch_times.push(time)
  h
}

///|
pub fn training_history_summary(history : TrainingHistory) -> String {
  let n = history.loss_values.length()
  let mut min_loss = 999999999999.0
  let mut final_loss = 0.0
  for i in 0..<n {
    if history.loss_values[i] < min_loss {
      min_loss = history.loss_values[i]
    }
    final_loss = history.loss_values[i]
  }
  let mut total_time = 0.0
  for i in 0..<history.epoch_times.length() {
    total_time = total_time + history.epoch_times[i]
  }
  "Epochs: " +
  n.to_string() +
  ", Min loss: " +
  min_loss.to_string() +
  ", Final loss: " +
  final_loss.to_string() +
  ", Total time: " +
  total_time.to_string() +
  "s"
}

///|
pub struct TrainRunner {
  max_epochs : Int32
  batch_size : Int32
  mut early_stopping : EarlyStoppingState?
  mut lr_scheduler : LRScheduler?
  mut checkpoint : ModelCheckpoint?
  mut history : TrainingHistory
}

///|
pub fn TrainRunner::new(max_epochs : Int32, batch_size : Int32) -> TrainRunner {
  {
    max_epochs,
    batch_size,
    early_stopping: None,
    lr_scheduler: None,
    checkpoint: None,
    history: TrainingHistory::new(),
  }
}

///|
pub fn train_runner_run(
  runner : TrainRunner,
  loss_fn : (Int32) -> Double,
) -> TrainRunner {
  let r = runner
  for epoch in 0..<r.max_epochs {
    let loss = loss_fn(epoch)
    r.history = training_history_record_loss(r.history, loss)
    match r.early_stopping {
      Some(es) => {
        let es2 = early_stopping_step(es, loss)
        r.early_stopping = Some(es2)
        if es2.should_stop {
          break
        }
      }
      None => ()
    }
    match r.lr_scheduler {
      Some(sched) => {
        let sched2 = lr_scheduler_step(sched, loss)
        r.lr_scheduler = Some(sched2)
      }
      None => ()
    }
    match r.checkpoint {
      Some(ckpt) => {
        let (ckpt2, _) = model_checkpoint_should_save(ckpt, loss)
        r.checkpoint = Some(ckpt2)
      }
      None => ()
    }
  }
  r
}

// Adversarial classifier for semi-supervised training

///|
pub struct AdversarialClassifier {
  n_latent : Int32
  n_labels : Int32
  weights : Array[Array[Double]]
  bias : Array[Double]
  lambda_ : Double
}

///|
pub fn adversarial_classifier_init(
  n_latent : Int32,
  n_labels : Int32,
  lambda_ : Double,
) -> AdversarialClassifier {
  let w = Array(capacity=n_labels)
  for _i in 0..<n_labels {
    let row = Array(capacity=n_latent)
    for _j in 0..<n_latent {
      row.push(0.01)
    }
    w.push(row)
  }
  let b = Array(capacity=n_labels)
  for _i in 0..<n_labels {
    b.push(0.0)
  }
  { n_latent, n_labels, weights: w, bias: b, lambda_, }
}

///|
pub fn adversarial_classifier_predict(
  cls : AdversarialClassifier,
  z : Array[Double],
) -> Array[Double] {
  let scores = Array(capacity=cls.n_labels)
  for l in 0..<cls.n_labels {
    let mut score = cls.bias[l]
    for j in 0..<z.length() {
      if j < cls.weights[l].length() {
        score = score + z[j] * cls.weights[l][j]
      }
    }
    scores.push(score)
  }
  scores
}

///|
pub fn adversarial_classifier_train_step(
  cls : AdversarialClassifier,
  z : Array[Double],
  true_label : Int32,
) -> AdversarialClassifier {
  let c = cls
  let preds = adversarial_classifier_predict(c, z)
  let mut max_p = preds[0]
  for i in 1..<preds.length() {
    if preds[i] > max_p {
      max_p = preds[i]
    }
  }
  let mut sum_exp = 0.0
  for i in 0..<preds.length() {
    sum_exp = sum_exp + @math.exp(preds[i] - max_p)
  }
  for l in 0..<c.n_labels {
    let prob = @math.exp(preds[l] - max_p) / sum_exp
    let target = if l == true_label { 1.0 } else { 0.0 }
    let grad = (target - prob) * c.lambda_
    for j in 0..<z.length() {
      if j < c.weights[l].length() {
        c.weights[l][j] = c.weights[l][j] + 0.01 * grad * z[j]
      }
    }
    c.bias[l] = c.bias[l] + 0.01 * grad
  }
  c
}

// Per-label subsampling for balanced training

///|
pub fn per_label_subsample(
  labels : Array[Int32],
  n_per_label : Int32,
  seed : Int32,
) -> Array[Int32] {
  let n = labels.length()
  let mut rng = seed
  let label_indices = Map([])
  for i in 0..<n {
    let l = labels[i]
    let existing = label_indices.get(l).unwrap_or([])
    let arr = existing
    arr.push(i)
    label_indices[l] = arr
  }
  let result = Array(capacity=n_per_label * label_indices.length())
  for key in label_indices.keys() {
    let indices = label_indices[key]
    let shuffled = indices
    for i in 0..<shuffled.length() {
      rng = (rng * 1103515245 + 12345) % 2147483647
      let j = rng % shuffled.length()
      if j >= 0 {
        let tmp = shuffled[i]
        shuffled[i] = shuffled[j]
        shuffled[j] = tmp
      }
    }
    for i in 0..<n_per_label {
      if i < shuffled.length() {
        result.push(shuffled[i])
      }
    }
  }
  result
}

///|
pub fn _scvi_train_dummy() -> Unit {
  let _ = EarlyStoppingConfig::{
    monitor: "",
    patience: 0,
    min_delta: 0.0,
    mode: "",
  }
  let _ = ModelCheckpoint::{
    save_best: false,
    best_score: 0.0,
    save_path: "",
    monitor: "",
    mode: "",
  }
  let _ = TrainRunner::{
    max_epochs: 0,
    batch_size: 0,
    early_stopping: None,
    lr_scheduler: None,
    checkpoint: None,
    history: TrainingHistory::new(),
  }
  let _ = TrainingHistory::{
    loss_values: [],
    val_values: [],
    metric_values: Map([]),
    epoch_times: [],
  }
  let _ = AdversarialClassifier::{
    n_latent: 0,
    n_labels: 0,
    weights: [],
    bias: [],
    lambda_: 0.0,
  }
}

///|
pub fn _scvi_train2_dummy() -> Unit {
  let _ = ReduceLROnPlateau
  let _ = CosineAnnealing
  let _ = StepLR
  let _ = ExponentialLR
}
