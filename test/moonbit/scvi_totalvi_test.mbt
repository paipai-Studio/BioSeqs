// scvi_totalvi_test.mbt — TotalVI 蛋白 fore/background 去卷积测试

///|
fn totalvi_test_model() -> @src.TotalVIModel {
  @src.totalvi_init(2, 2, @src.TotalVIConfig::new(4, 2, 1, 0.1, 0.01, 5))
}

///|
test "totalvi mixture posterior math" {
  // y 与 λ_fg 匹配 → 前景后验高；与 λ_bg 匹配 → 低
  let fg = @src.totalvi_mixture_posterior(10.0, 10.0, 1.0, 0.5)
  assert_true(fg > 0.99)
  let bg = @src.totalvi_mixture_posterior(10.0, 1.0, 10.0, 0.5)
  assert_true(bg < 0.01)
  // 似然相等时后验 = 先验
  assert_true(
    @moonbit.assert_num_close(
      @src.totalvi_mixture_posterior(5.0, 5.0, 5.0, 0.3),
      0.3,
      1.0e-9,
    ),
  )
  // 低先验压低后验
  let weak = @src.totalvi_mixture_posterior(10.0, 10.0, 1.0, 0.01)
  assert_true(weak < fg)
  assert_true(weak > 0.0)
  // 边界：负计数按 0、速率/先验钳制后仍在 (0,1)
  let edge = @src.totalvi_mixture_posterior(-1.0, 0.0, 0.0, 2.0)
  assert_true(edge > 0.0 && edge < 1.0)
}

///|
test "totalvi protein mixture and foreground probability" {
  let m = totalvi_test_model()
  let rna = [[1.0, 0.5], [0.2, 0.1]]
  let prot = [[3.0, 0.0], [1.0, 2.0]]
  let probs = @src.totalvi_get_foreground_probability(m, rna, prot)
  assert_eq(probs.length(), 2)
  assert_eq(probs[0].length(), 2)
  for row in probs {
    for p in row {
      assert_true(p >= 0.0 && p <= 1.0)
    }
  }
  // 与逐细胞 mixture 分解一致
  let z = @src.totalvi_encode(m, rna[0], prot[0])
  let mix = @src.totalvi_protein_mixture(m, z, prot[0])
  assert_eq(mix.length(), 2)
  assert_true(mix[0].foreground_rate > 0.0)
  assert_true(mix[0].background_rate > 0.0)
  assert_true(mix[0].prior > 0.0 && mix[0].prior < 1.0)
  assert_true(
    @moonbit.assert_num_close(
      mix[0].foreground_probability,
      probs[0][0],
      1.0e-12,
    ),
  )
  // 后验可由纯函数用模型参数复算
  let recomputed = @src.totalvi_mixture_posterior(
    prot[0][0],
    mix[0].foreground_rate,
    mix[0].background_rate,
    mix[0].prior,
  )
  assert_true(@moonbit.assert_num_close(recomputed, probs[0][0], 1.0e-12))
  // 校正表达 = 计数 × 后验，不超过原始计数
  let corrected = @src.totalvi_corrected_protein_expression(m, rna, prot)
  assert_true(
    @moonbit.assert_num_close(
      corrected[0][0],
      prot[0][0] * probs[0][0],
      1.0e-12,
    ),
  )
  assert_true(corrected[0][0] <= prot[0][0] + 1.0e-12)
}

///|
test "totalvi fit protein background EM" {
  let m = totalvi_test_model()
  // 蛋白 0 高计数（远超随机初始化的前景率 ~0.7）→ EM 把背景均值拉向数据
  // 蛋白 1 计数低 → 背景均值维持低位
  let n = 6
  let rna : Array[Array[Double]] = Array(capacity=n)
  let prot : Array[Array[Double]] = Array(capacity=n)
  for i in 0..<n {
    rna.push([1.0 + i.to_double() * 0.1, 0.5])
    prot.push([8.0, 1.0])
  }
  assert_true(@moonbit.assert_num_close(m.prot_bg_mean[0], 1.0, 1.0e-12))
  let probs = @src.totalvi_fit_protein_background(m, rna, prot, n_iters=5)
  assert_eq(probs.length(), n)
  assert_eq(probs[0].length(), 2)
  // 背景率 0 被拉向高计数数据（加权均值 ≈ 8）
  assert_true(m.prot_bg_mean[0] > 5.0)
  // 高计数蛋白的先验被压低（后验均值远小于 0.5 → logit 变负）
  assert_true(m.prot_mix_logit[0] < 0.0)
  // 上调背景后，高计数蛋白的前景后验低于低计数蛋白
  for i in 0..<n {
    for p in probs[i] {
      assert_true(p >= 0.0 && p <= 1.0)
    }
    assert_true(probs[i][0] < probs[i][1])
  }
}

///|
test "totalvi background edge cases" {
  let m = totalvi_test_model()
  // rna / protein 行数不一致 → 空
  assert_eq(
    @src.totalvi_get_foreground_probability(m, [[1.0, 1.0]], []).length(),
    0,
  )
  assert_eq(
    @src.totalvi_fit_protein_background(m, [[1.0, 1.0]], [
      [1.0, 1.0],
      [2.0, 2.0],
    ]).length(),
    0,
  )
  // 无蛋白模型 → 空
  let m0 = @src.totalvi_init(
    2,
    0,
    @src.TotalVIConfig::new(4, 2, 1, 0.1, 0.01, 5),
  )
  assert_eq(
    @src.totalvi_get_foreground_probability(m0, [[1.0, 1.0]], [[]]).length(),
    0,
  )
  // 单细胞、零计数仍可计算（后验由似然决定）
  let zero = @src.totalvi_get_foreground_probability(
    totalvi_test_model(),
    [[0.0, 0.0]],
    [[0.0, 0.0]],
  )
  assert_eq(zero.length(), 1)
  assert_eq(zero[0].length(), 2)
  for p in zero[0] {
    assert_true(p >= 0.0 && p <= 1.0)
  }
}
