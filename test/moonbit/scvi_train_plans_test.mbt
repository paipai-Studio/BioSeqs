// scvi_train_plans_test.mbt — KL 退火 + n_obs 数据量渐增训练计划测试
// ===========================================================================

///|
test "scvi_plan_kl_weight_schedules" {
  // 恒定权重（默认，等价旧 scvi_train）
  assert_eq(@src.scvi_plan_kl_weight(@src.KlConstant, 0), 1.0)
  assert_eq(@src.scvi_plan_kl_weight(@src.KlConstant, 50), 1.0)
  // 线性 0→1
  let lin = @src.KlLinearWarmup(10)
  assert_eq(@src.scvi_plan_kl_weight(lin, 0), 0.0)
  assert_true(
    @moonbit.assert_num_close(@src.scvi_plan_kl_weight(lin, 5), 0.5, 1.0e-12),
  )
  assert_eq(@src.scvi_plan_kl_weight(lin, 10), 1.0)
  assert_eq(@src.scvi_plan_kl_weight(lin, 100), 1.0)
  // warmup 长度 <= 0 视为无退火
  assert_eq(@src.scvi_plan_kl_weight(@src.KlLinearWarmup(0), 3), 1.0)
  // 步进：warmup=10 分 4 级（每级 2.5 epoch）
  let step = @src.KlStepWarmup(10, 4)
  assert_eq(@src.scvi_plan_kl_weight(step, 0), 0.0)
  assert_true(
    @moonbit.assert_num_close(@src.scvi_plan_kl_weight(step, 3), 0.25, 1.0e-12),
  )
  assert_eq(@src.scvi_plan_kl_weight(step, 5), 0.5)
  assert_eq(@src.scvi_plan_kl_weight(step, 10), 1.0)
  assert_eq(@src.scvi_plan_kl_weight(step, 99), 1.0)
}

///|
test "scvi_plan_n_obs_ramp" {
  let cfg = @src.ScviTrainingPlanConfig::with_plan(
    "unsupervised",
    100,
    0.01,
    false,
    3,
    0.0,
    @src.KlConstant,
    10,
    10,
    5,
  )
  // 阶段 0（epoch 0-4）：10 个细胞
  assert_eq(@src.scvi_plan_n_obs(cfg, 0, 100), 10)
  assert_eq(@src.scvi_plan_n_obs(cfg, 4, 100), 10)
  // 阶段 1（epoch 5-9）：20 个
  assert_eq(@src.scvi_plan_n_obs(cfg, 5, 100), 20)
  assert_eq(@src.scvi_plan_n_obs(cfg, 9, 100), 20)
  // 阶段 2：30 个
  assert_eq(@src.scvi_plan_n_obs(cfg, 10, 100), 30)
  // 上限夹紧到全量
  assert_eq(@src.scvi_plan_n_obs(cfg, 100, 100), 100)
  // 数据比起始量还少 → 直接全量
  assert_eq(@src.scvi_plan_n_obs(cfg, 0, 8), 8)
  // 默认构造（n_obs_start=0）→ 始终全量
  let full = @src.ScviTrainingPlanConfig::new(
    "unsupervised", 100, 0.01, false, 3, 0.0,
  )
  assert_eq(@src.scvi_plan_n_obs(full, 0, 100), 100)
  assert_eq(@src.scvi_plan_n_obs(full, 42, 100), 100)
}

///|
test "scvi_train_with_plan_improves_elbo" {
  let x = [
    [5.0, 0.0, 2.0, 0.0, 1.0, 3.0],
    [0.0, 4.0, 1.0, 2.0, 0.0, 1.0],
    [3.0, 1.0, 0.0, 5.0, 2.0, 0.0],
    [1.0, 2.0, 3.0, 0.0, 4.0, 2.0],
  ]
  let model = @src.scvi_init_model(6, 1, @src.ScVIConfig::default())
  let mut before = 0.0
  for i in 0..<x.length() {
    before += @src.scvi_compute_elbo(model, x[i])
  }
  // 线性 KL warmup（前 10 个 epoch 权重 0→1），数据全量
  let cfg = @src.ScviTrainingPlanConfig::with_plan(
    "unsupervised",
    40,
    0.01,
    false,
    3,
    0.0,
    @src.KlLinearWarmup(10),
    0,
    0,
    10,
  )
  let trained = @src.scvi_train_with_plan(model, x, cfg)
  assert_true(trained.trained)
  let mut after = 0.0
  for i in 0..<x.length() {
    after += @src.scvi_compute_elbo(trained, x[i])
  }
  assert_true(after > before)
  // latent 有限且形状正确
  let latent = @src.scvi_get_latent_representation(trained, x)
  assert_eq(latent.length(), 4)
  assert_true(latent[0][0] == latent[0][0])
}

///|
test "scvi_train_with_plan_data_ramp_runs" {
  let x = [
    [5.0, 0.0, 2.0, 0.0, 1.0, 3.0],
    [0.0, 4.0, 1.0, 2.0, 0.0, 1.0],
    [3.0, 1.0, 0.0, 5.0, 2.0, 0.0],
    [1.0, 2.0, 3.0, 0.0, 4.0, 2.0],
  ]
  let model = @src.scvi_init_model(6, 1, @src.ScVIConfig::default())
  // 数据量渐增：2 → 4 → 4（夹紧），叠加步进 KL warmup
  let cfg = @src.ScviTrainingPlanConfig::with_plan(
    "unsupervised",
    30,
    0.01,
    false,
    3,
    0.0,
    @src.KlStepWarmup(6, 2),
    2,
    2,
    2,
  )
  let trained = @src.scvi_train_with_plan(model, x, cfg)
  assert_true(trained.trained)
  let latent = @src.scvi_get_latent_representation(trained, x)
  assert_eq(latent.length(), 4)
  for v in latent[0] {
    assert_true(v == v)
  }
}

///|
test "scvi_training_plan_early_stopping_and_losses" {
  // 损失 epoch0=1.0、之后 0.5：patience=3 → 连续 3 个不改善后提前停止
  let config = @src.ScviTrainingPlanConfig::new(
    "unsupervised", 50, 0.1, true, 3, 0.0,
  )
  let result = @src.scvi_training_plan_run(config, fn(epoch) {
    if epoch == 0 {
      1.0
    } else {
      0.5
    }
  })
  assert_true(result.stopped_early)
  assert_true(result.epochs_run < 50)
  assert_true(result.best_loss <= result.final_loss + 1.0e-12)
  assert_true((result.best_loss - 0.5).abs() < 1.0e-12)
  let recent = @src.scvi_plan_recent_losses(result, 2)
  assert_eq(recent.length(), 2)
  assert_true((recent[1] - 0.5).abs() < 1.0e-12)
  // 持续改善 → 不早停
  let result2 = @src.scvi_training_plan_run(config, fn(epoch) {
    1.0 / (epoch + 1).to_double()
  })
  assert_true(!result2.stopped_early)
  assert_eq(result2.epochs_run, 50)
  // 关闭早停 → 跑满
  let no_stop = @src.ScviTrainingPlanConfig::new(
    "unsupervised", 10, 0.1, false, 3, 0.0,
  )
  let result3 = @src.scvi_training_plan_run(no_stop, fn(_epoch) { 0.7 })
  assert_true(!result3.stopped_early)
  assert_eq(result3.epochs_run, 10)
  // 损失辅助项
  assert_true(@src.scvi_semisupervised_loss(1.0, 0.5, 0.2) == 1.1)
  assert_true(@src.scvi_classifier_loss([1.0, 0.0], 0) < 1.0e-6)
  assert_true(@src.scvi_entropy([1.0, 0.0]) < 1.0e-12)
  let adv = @src.scvi_adversarial_loss(0.9, 1)
  assert_true(adv > 0.0 && adv < 0.2)
}
