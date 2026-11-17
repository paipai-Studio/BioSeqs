// scvi_train_plans.mbt — 训练计划：早停调度 + KL 退火 + n_obs 数据量渐增
//
// 对应 scvi-tools `JaxTrainingPlan` / `TrainRunner` 的三件稳定性关键件：
// - 早停（early_stopping + patience + min_delta）：验证损失不再改善即停；
// - KL warmup：kl_weight 从 0 线性或分步退火到 1，避免训练初期强正则把
//   后验压成先验（posterior collapse）；
// - n_obs 调度：每 epoch 使用的细胞数从 n_obs_start 起按 n_obs_step 逐阶
//   递增（每阶段 n_obs_phase_epochs 个 epoch），小样本先热身再全量。
// 另附 SCANVI/TotalVI 常用损失项（半监督加权、分类交叉熵、熵、对抗）。

///|
/// KL 权重退火策略。
///
/// - `KlConstant`：全程 kl_weight = 1（与旧 `scvi_train` 行为一致）
/// - `KlLinearWarmup(n)`：前 n 个 epoch 线性 0→1，之后保持 1
/// - `KlStepWarmup(n, phases)`：前 n 个 epoch 分 `phases` 级阶梯上升
pub(all) enum ScviKlWarmup {
  KlConstant
  KlLinearWarmup(Int32)
  KlStepWarmup(Int32, Int32)
}

///|
pub struct ScviTrainingPlanConfig {
  /// 任务类型标签（如 "unsupervised" / "semi-supervised"），由上层训练器解释
  plan_type : String
  max_epochs : Int32
  learning_rate : Double
  early_stopping : Bool
  patience : Int32
  min_delta : Double
  /// KL 退火策略
  kl_warmup : ScviKlWarmup
  /// 每 epoch 起始细胞数；<= 0 表示始终全量
  n_obs_start : Int32
  /// 每阶段递增的细胞数
  n_obs_step : Int32
  /// 每阶段持续的 epoch 数
  n_obs_phase_epochs : Int32
}

///|
pub fn ScviTrainingPlanConfig::new(
  plan_type : String,
  max_epochs : Int32,
  learning_rate : Double,
  early_stopping : Bool,
  patience : Int32,
  min_delta : Double,
) -> ScviTrainingPlanConfig {
  ScviTrainingPlanConfig::{
    plan_type,
    max_epochs,
    learning_rate,
    early_stopping,
    patience,
    min_delta,
    kl_warmup: KlConstant,
    n_obs_start: 0,
    n_obs_step: 0,
    n_obs_phase_epochs: 10,
  }
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
/// 完整构造：额外指定 KL 退火策略与 n_obs 数据量渐增调度。
pub fn ScviTrainingPlanConfig::with_plan(
  plan_type : String,
  max_epochs : Int32,
  learning_rate : Double,
  early_stopping : Bool,
  patience : Int32,
  min_delta : Double,
  kl_warmup : ScviKlWarmup,
  n_obs_start : Int32,
  n_obs_step : Int32,
  n_obs_phase_epochs : Int32,
) -> ScviTrainingPlanConfig {
  ScviTrainingPlanConfig::{
    plan_type,
    max_epochs,
    learning_rate,
    early_stopping,
    patience,
    min_delta,
    kl_warmup,
    n_obs_start,
    n_obs_step,
    n_obs_phase_epochs,
  }
}

///|
/// 第 `epoch` 个 epoch（0 基）的 kl_weight。warmup 长度 <= 0 视为无退火。
pub fn scvi_plan_kl_weight(warmup : ScviKlWarmup, epoch : Int32) -> Double {
  match warmup {
    KlConstant => 1.0
    KlLinearWarmup(n) =>
      if n <= 0 {
        1.0
      } else {
        let w = epoch.to_double() / n.to_double()
        if w < 1.0 {
          w
        } else {
          1.0
        }
      }
    KlStepWarmup(n, phases) =>
      if n <= 0 {
        1.0
      } else {
        let p = if phases < 1 { 1 } else { phases }
        // 当前所处阶梯：epoch*p/n 向下取整，warmup 结束后到顶
        let mut step = (epoch.to_double() * p.to_double() / n.to_double())
          .floor()
          .to_int()
        if step < 0 {
          step = 0
        }
        if step > p {
          step = p
        }
        step.to_double() / p.to_double()
      }
  }
}

///|
/// 第 `epoch` 个 epoch 使用的细胞数：`n_obs_start + 阶段数 * n_obs_step`，
/// 夹在 [1, n_cells] 内。`n_obs_start <= 0` 或 `n_obs_step <= 0` 时始终
/// 返回全量 `n_cells`。
pub fn scvi_plan_n_obs(
  config : ScviTrainingPlanConfig,
  epoch : Int32,
  n_cells : Int32,
) -> Int32 {
  if config.n_obs_start <= 0 || config.n_obs_step <= 0 || n_cells <= 0 {
    return n_cells
  }
  let phase_len = if config.n_obs_phase_epochs < 1 {
    1
  } else {
    config.n_obs_phase_epochs
  }
  let phase = epoch / phase_len
  let mut n = config.n_obs_start + phase * config.n_obs_step
  if n < 1 {
    n = 1
  }
  if n > n_cells {
    n = n_cells
  }
  n
}

///|
pub struct ScviTrainingPlanResult {
  /// 是否因验证损失停滞而提前停止
  stopped_early : Bool
  /// 实际执行的 epoch 数
  epochs_run : Int32
  /// 历史最优损失
  best_loss : Double
  /// 最后一个 epoch 的损失
  final_loss : Double
  /// 每个 epoch 的损失轨迹
  losses : Array[Double]
}

///|
/// 以通用 epoch 步进函数 `epoch_fn(epoch) -> loss` 执行训练计划：逐 epoch
/// 记录损失，`early_stopping` 开启时连续 `patience` 个 epoch 未见
/// `min_delta` 以上改善即提前停止。第一个 epoch 总是执行。
pub fn scvi_training_plan_run(
  config : ScviTrainingPlanConfig,
  epoch_fn : (Int32) -> Double,
) -> ScviTrainingPlanResult {
  let losses : Array[Double] = Array(capacity=config.max_epochs)
  let mut best = 1.0e300
  let mut final_loss = 1.0e300
  let mut wait = 0
  let mut epochs_run = 0
  let mut stopped_early = false
  for epoch in 0..<config.max_epochs {
    let loss = epoch_fn(epoch)
    losses.push(loss)
    epochs_run = epoch + 1
    final_loss = loss
    if loss < best - config.min_delta {
      best = loss
      wait = 0
    } else {
      wait = wait + 1
      if config.early_stopping && wait >= config.patience {
        stopped_early = true
        break
      }
    }
  }
  { stopped_early, epochs_run, best_loss: best, final_loss, losses, }
}

///|
/// 取训练轨迹的最后 `n` 个 epoch 损失（n 超过轨迹长度时返回全部）。
pub fn scvi_plan_recent_losses(
  result : ScviTrainingPlanResult,
  n : Int32,
) -> Array[Double] {
  let m = result.losses.length()
  let take = if n < 0 { 0 } else if n > m { m } else { n }
  let out : Array[Double] = Array(capacity=take)
  for i in (m - take)..<m {
    out.push(result.losses[i])
  }
  out
}

///|
/// SCANVI 半监督损失：重建项 + `cls_weight` × 分类项。
pub fn scvi_semisupervised_loss(
  recon_loss : Double,
  cls_loss : Double,
  cls_weight : Double,
) -> Double {
  recon_loss + cls_weight * cls_loss
}

///|
/// 分类交叉熵：目标类别的概率的负对数（概率截断在 [1e-12, 1]）。
pub fn scvi_classifier_loss(probs : Array[Double], target : Int32) -> Double {
  if target < 0 || target >= probs.length() {
    return 0.0
  }
  let p = probs[target]
  let p_safe = if p < 1.0e-12 { 1.0e-12 } else if p > 1.0 { 1.0 } else { p }
  -@math.ln(p_safe)
}

///|
/// 概率向量的 Shannon 熵（自然对数）。
pub fn scvi_entropy(probs : Array[Double]) -> Double {
  let mut h = 0.0
  for p in probs {
    if p > 0.0 {
      h = h - p * @math.ln(p)
    }
  }
  h
}

///|
/// 判别器二元交叉熵：`label` 为 1 时惩罚判真失败（-ln p），为 0 时惩罚
/// 判假失败（-ln(1-p)）。
pub fn scvi_adversarial_loss(pred : Double, label : Int32) -> Double {
  let p = if pred < 1.0e-12 {
    1.0e-12
  } else if pred > 1.0 - 1.0e-12 {
    1.0 - 1.0e-12
  } else {
    pred
  }
  if label == 1 {
    -@math.ln(p)
  } else {
    -@math.ln(1.0 - p)
  }
}

///|
/// 带训练计划的 scVI 训练：逐 epoch 按 `scvi_plan_kl_weight` 给 KL 项加权、
/// 按 `scvi_plan_n_obs` 取数据前缀子集训练，其余（Adam、梯度裁剪、库大小
/// 缩放的 NB 重建）与 `scvi_train` 完全一致。
pub fn scvi_train_with_plan(
  model : ScVIModel,
  x : Array[Array[Double]],
  config : ScviTrainingPlanConfig,
) -> ScVIModel {
  let m = model
  if x.length() == 0 {
    m.trained = true
    return m
  }
  let n_cells = x.length()
  let theta = Array(capacity=m.n_genes)
  for g in 0..<m.n_genes {
    theta.push(@math.exp(m.log_theta[g]))
  }
  // Adam 动量按模型参数形状预分配（注册顺序须与每 epoch 一致）
  let template = ScviAgTape::new()
  ignore(template.matrix_param(m.encoder_w1))
  ignore(template.vector_param(m.encoder_b1))
  ignore(template.matrix_param(m.encoder_w_mu))
  ignore(template.vector_param(m.encoder_b_mu))
  ignore(template.matrix_param(m.encoder_w_lv))
  ignore(template.vector_param(m.encoder_b_lv))
  ignore(template.matrix_param(m.decoder_w1))
  ignore(template.vector_param(m.decoder_b1))
  ignore(template.matrix_param(m.decoder_w_out))
  ignore(template.vector_param(m.decoder_b_out))
  let adam = scvi_adam_init(template, config.learning_rate)
  for epoch in 0..<config.max_epochs {
    let kl_w = scvi_plan_kl_weight(config.kl_warmup, epoch)
    let n_used = scvi_plan_n_obs(config, epoch, n_cells)
    let t = ScviAgTape::new()
    let w1 = t.matrix_param(m.encoder_w1)
    let b1 = t.vector_param(m.encoder_b1)
    let w_mu = t.matrix_param(m.encoder_w_mu)
    let b_mu = t.vector_param(m.encoder_b_mu)
    let w_lv = t.matrix_param(m.encoder_w_lv)
    let b_lv = t.vector_param(m.encoder_b_lv)
    let wd1 = t.matrix_param(m.decoder_w1)
    let bd1 = t.vector_param(m.decoder_b1)
    let w_out = t.matrix_param(m.decoder_w_out)
    let b_out = t.vector_param(m.decoder_b_out)
    let mut total = t.leaf([0.0])
    for i in 0..<n_used {
      let xi = x[i]
      let lib = scvi_array_sum(xi)
      let lib_size = if lib <= 0.0 { 1.0 } else { lib }
      let xin = t.leaf(xi)
      let hidden = t.relu(t.dense(xin, w1, b1))
      let mu = t.dense(hidden, w_mu, b_mu)
      let log_var = t.dense(hidden, w_lv, b_lv)
      let decoded = t.relu(t.dense(mu, wd1, bd1))
      let raw = t.dense(decoded, w_out, b_out)
      let rates = t.softmax(raw)
      let mean_counts = t.scale(rates, lib_size)
      let nll = t.nb_nll(mean_counts, xi, theta)
      let kl = t.scale(t.kl_divergence(mu, log_var), kl_w)
      total = t.add(total, t.add(nll, kl))
    }
    let loss = t.scale(total, 1.0 / n_used.to_double())
    t.backward(loss)
    ignore(t.clip_grad_norm(5.0))
    scvi_adam_step(t, adam)
  }
  m.trained = true
  m
}
