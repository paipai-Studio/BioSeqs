// ===========================================================================
// gwas_regression_test.mbt — GWAS 回归关联闭环测试.
// ===========================================================================

///|
fn gr_variants() -> Array[@src.GwasVariantMeta] {
  [
    @src.GwasVariantMeta::{ chrom: "chr1", pos: 100, id: "rs0", },
    @src.GwasVariantMeta::{ chrom: "chr1", pos: 200, id: "rs1", },
    @src.GwasVariantMeta::{ chrom: "chr2", pos: 300, id: "rs2", },
  ]
}

///|
/// 40 样本: rs0 与表型强关联, rs1/rs2 独立.
fn gr_dosages() -> Array[Array[Double]] {
  let n = 40
  let dosages : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let is_case = i % 2 == 0
    let rank = i / 2
    let d0 = if is_case {
      if rank < 15 {
        2.0
      } else if rank < 18 {
        1.0
      } else {
        0.0
      }
    } else if rank < 15 {
      0.0
    } else if rank < 18 {
      1.0
    } else {
      2.0
    }
    let d1 = (i % 3).to_double()
    let d2 = (i * 2 % 3).to_double()
    dosages.push([d0, d1, d2])
  }
  dosages
}

///|
fn gr_phenotype_case_control() -> Array[Double] {
  let y : Array[Double] = Array(capacity=40)
  for i = 0; i < 40; i = i + 1 {
    y.push(if i % 2 == 0 { 1.0 } else { 0.0 })
  }
  y
}

///|
fn gr_phenotype_quantitative() -> Array[Double] {
  let dosages = gr_dosages()
  let y : Array[Double] = Array(capacity=40)
  for i = 0; i < 40; i = i + 1 {
    y.push(dosages[i][0] * 0.8 + (i % 5).to_double() * 0.05)
  }
  y
}

///|
fn gr_find(
  hits : Array[@src.GwasRegressionHit],
  id : String,
) -> @src.GwasRegressionHit? {
  for hit in hits {
    if hit.variant.id == id {
      return Some(hit)
    }
  }
  None
}

// ===========================================================================
// 逻辑回归关联
// ===========================================================================

///|
test "gwas_regression_logistic_flags_signal" {
  let hits = @src.gwas_regression_logistic(
    gr_dosages(),
    gr_phenotype_case_control(),
    gr_variants(),
  )
  assert_eq(hits.length(), 3)
  match gr_find(hits, "rs0") {
    Some(hit) => {
      assert_true(hit.beta > 0.0)
      assert_true(hit.p_value < 0.01)
    }
    None => assert_true(false)
  }
  match gr_find(hits, "rs1") {
    Some(hit) => assert_true(hit.p_value > 0.01)
    None => assert_true(false)
  }
}

// ===========================================================================
// 线性回归关联
// ===========================================================================

///|
test "gwas_regression_linear_recovers_effect" {
  let hits = @src.gwas_regression_linear(
    gr_dosages(),
    gr_phenotype_quantitative(),
    gr_variants(),
  )
  match gr_find(hits, "rs0") {
    Some(hit) => {
      assert_true((hit.beta - 0.8).abs() < 0.2)
      assert_true(hit.p_value < 0.01)
    }
    None => assert_true(false)
  }
}

// ===========================================================================
// PCA 协变量
// ===========================================================================

///|
test "gwas_regression_pca_covariates_shape" {
  let cov = @src.gwas_regression_pca_covariates(gr_dosages(), 2)
  assert_eq(cov.length(), 40)
  assert_eq(cov[0].length(), 2)
}

///|
test "gwas_regression_with_covariates_runs" {
  // 独立于基因型的协变量, 校正后信号应保留.
  let cov : Array[Array[Double]] = Array(capacity=40)
  for i = 0; i < 40; i = i + 1 {
    cov.push([(i % 4).to_double()])
  }
  let hits = @src.gwas_regression_logistic(
    gr_dosages(),
    gr_phenotype_case_control(),
    gr_variants(),
    covariates=Some(cov),
  )
  assert_eq(hits.length(), 3)
  match gr_find(hits, "rs0") {
    Some(hit) => assert_true(hit.p_value < 0.05)
    None => assert_true(false)
  }
}

///|
test "gwas_regression_pca_adjustment_runs" {
  // PCA 协变量(群体分层校正)应可正常参与回归.
  let cov = @src.gwas_regression_pca_covariates(gr_dosages(), 2)
  let hits = @src.gwas_regression_logistic(
    gr_dosages(),
    gr_phenotype_case_control(),
    gr_variants(),
    covariates=Some(cov),
  )
  assert_eq(hits.length(), 3)
  for hit in hits {
    assert_true(hit.p_value >= 0.0 && hit.p_value <= 1.0)
  }
}

// ===========================================================================
// 汇总 / 导出 / 可视化衔接
// ===========================================================================

///|
test "gwas_regression_summary_and_exports" {
  let hits = @src.gwas_regression_logistic(
    gr_dosages(),
    gr_phenotype_case_control(),
    gr_variants(),
  )
  let summary = @src.gwas_regression_summary(hits, genome_wide=0.05)
  assert_eq(summary.n_tested, 3)
  assert_true(summary.lambda_gc > 0.0)
  assert_true(summary.n_genome_wide >= 1)
  let tsv = @src.gwas_regression_tsv(hits)
  assert_true(tsv.contains("chrom\tpos\tid"))
  assert_true(tsv.contains("rs0"))
  let json = @src.gwas_regression_json(summary)
  assert_true(json.contains("\"lambda_gc\""))
  assert_true(json.contains("\"p_value\""))
  assert_true(@src.parse_json(json) is Some(_))
  assert_true(@src.gwas_regression_summary_line(summary).contains("tested=3"))
}

///|
test "gwas_regression_manhattan_and_qq_delegate" {
  let hits = @src.gwas_regression_logistic(
    gr_dosages(),
    gr_phenotype_case_control(),
    gr_variants(),
  )
  let summary = @src.gwas_regression_summary(hits)
  let manhattan = @src.gwas_regression_manhattan(summary)
  assert_eq(manhattan.points.length(), 3)
  let qq = @src.gwas_regression_qq(summary)
  assert_eq(qq.points.length(), 3)
  assert_true(qq.lambda_gc > 0.0)
}
