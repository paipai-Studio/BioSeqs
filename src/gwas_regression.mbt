///|
/// gwas_regression.mbt — GWAS 回归关联分析闭环(本地确定性实现).
///
/// 整合现有基础设施:
/// - 逐变异位点逻辑/线性回归(GLM IRLS + Wald 检验)估算效应量与 p 值;
/// - 可传入协变量(PCA 群体分层校正, 复用 skallel_pca);
/// - 汇总 λGC(复用 gwas_plot)与结果 TSV/JSON 导出;
/// - 输出可直接喂给 gwas_plot 的 Manhattan/QQ 数据层.
///
/// 数据不出本机.
pub(all) struct GwasVariantMeta {
  chrom : String
  pos : Int32
  id : String
}

///|
/// 单个变异位点的回归关联结果.
pub(all) struct GwasRegressionHit {
  variant : GwasVariantMeta
  /// 效应量(每增加一个 alt 等位基因).
  beta : Double
  /// 效应量标准误.
  se : Double
  /// Wald z 统计量.
  z : Double
  /// 双侧 p 值.
  p_value : Double
  /// 参与回归的已分型样本数.
  n_called : Int32
  /// alt 等位基因频率.
  alt_freq : Double
}

///|
/// 关联分析汇总.
pub(all) struct GwasRegressionSummary {
  hits : Array[GwasRegressionHit]
  /// 基因组膨胀因子 λGC.
  lambda_gc : Double
  /// 通过质控的检验位点数.
  n_tested : Int32
  /// 达全基因组显著阈值的位点数.
  n_genome_wide : Int32
  genome_wide_threshold : Double
}

///|
/// 从 PCA 主成分构造协变量矩阵(样本 × 主成分).
///
/// dosages 形状为 (n_samples × n_variants); 内部转置为 (n_variants × n_samples)
/// 交给 skallel_pca(Patterson 标度).
pub fn gwas_regression_pca_covariates(
  dosages : Array[Array[Double]],
  n_components : Int32,
) -> Array[Array[Double]] {
  let n_samples = dosages.length()
  let n_variants = if n_samples > 0 { dosages[0].length() } else { 0 }
  let gn : Array[Array[Int32]] = Array(capacity=n_variants)
  for v = 0; v < n_variants; v = v + 1 {
    let row : Array[Int32] = Array(capacity=n_samples)
    for s = 0; s < n_samples; s = s + 1 {
      let d = if v < dosages[s].length() { dosages[s][v] } else { 0.0 }
      row.push(if d < 0.0 { 0 } else { d.round().to_int() })
    }
    gn.push(row)
  }
  skallel_pca(gn, n_components).coords
}

///|
fn gwas_regression_core(
  dosages : Array[Array[Double]],
  phenotype : Array[Double],
  variants : Array[GwasVariantMeta],
  covariates : Array[Array[Double]],
  min_af : Double,
  min_call_rate : Double,
  logistic : Bool,
) -> Array[GwasRegressionHit] {
  let n = phenotype.length()
  let k = if covariates.length() > 0 { covariates[0].length() } else { 0 }
  let hits : Array[GwasRegressionHit] = Array(capacity=variants.length())
  for j = 0; j < variants.length(); j = j + 1 {
    let col : Array[Double] = Array::make(n, 0.0)
    let mut called = 0
    let mut sum = 0.0
    for s = 0; s < n; s = s + 1 {
      let d = if j < dosages[s].length() { dosages[s][j] } else { -1.0 }
      if d >= 0.0 {
        called = called + 1
        sum = sum + d
      }
      col[s] = d
    }
    let call_rate = if n > 0 { called.to_double() / n.to_double() } else { 0.0 }
    let mean_d = if called > 0 { sum / called.to_double() } else { 0.0 }
    let alt_freq = mean_d / 2.0
    let maf = if alt_freq < 1.0 - alt_freq { alt_freq } else { 1.0 - alt_freq }
    if call_rate < min_call_rate || maf <= 0.0 || maf < min_af {
      continue
    }
    for s = 0; s < n; s = s + 1 {
      if col[s] < 0.0 {
        col[s] = mean_d
      }
    }
    let x : Array[Array[Double]] = Array(capacity=n)
    for s = 0; s < n; s = s + 1 {
      let row : Array[Double] = [1.0, col[s]]
      for c = 0; c < k; c = c + 1 {
        row.push(covariates[s][c])
      }
      x.push(row)
    }
    let family = if logistic {
      glm_family_binomial()
    } else {
      glm_family_gaussian()
    }
    let fit = glm_fit_irls(x, phenotype.copy(), family)
    let wald = glm_wald_test(fit, 1)
    hits.push(GwasRegressionHit::{
      variant: variants[j],
      beta: fit.coefficients[1],
      se: fit.std_errors[1],
      z: wald.statistic(),
      p_value: wald.p_value(),
      n_called: called,
      alt_freq,
    })
  }
  hits
}

///|
/// 逻辑回归关联(二值表型).
pub fn gwas_regression_logistic(
  dosages : Array[Array[Double]],
  phenotype : Array[Double],
  variants : Array[GwasVariantMeta],
  covariates~ : Array[Array[Double]]? = None,
  min_af~ : Double = 0.0,
  min_call_rate~ : Double = 0.0,
) -> Array[GwasRegressionHit] {
  let cov = match covariates {
    Some(c) => c
    None => []
  }
  gwas_regression_core(
    dosages, phenotype, variants, cov, min_af, min_call_rate, true,
  )
}

///|
/// 线性回归关联(数量表型).
pub fn gwas_regression_linear(
  dosages : Array[Array[Double]],
  phenotype : Array[Double],
  variants : Array[GwasVariantMeta],
  covariates~ : Array[Array[Double]]? = None,
  min_af~ : Double = 0.0,
  min_call_rate~ : Double = 0.0,
) -> Array[GwasRegressionHit] {
  let cov = match covariates {
    Some(c) => c
    None => []
  }
  gwas_regression_core(
    dosages, phenotype, variants, cov, min_af, min_call_rate, false,
  )
}

///|
/// 关联结果 -> gwas_plot 标记数组.
pub fn gwas_regression_markers(
  hits : Array[GwasRegressionHit],
) -> Array[GwasMarker] {
  let out : Array[GwasMarker] = Array(capacity=hits.length())
  for hit in hits {
    out.push(GwasMarker::{
      chrom: hit.variant.chrom,
      pos: hit.variant.pos,
      p_value: hit.p_value,
    })
  }
  out
}

///|
/// 汇总关联结果: 计算 λGC 与显著位点数.
pub fn gwas_regression_summary(
  hits : Array[GwasRegressionHit],
  genome_wide~ : Double = 5.0e-8,
) -> GwasRegressionSummary {
  let markers = gwas_regression_markers(hits)
  let lambda = gwas_lambda_gc(markers)
  let mut sig = 0
  for hit in hits {
    if hit.p_value < genome_wide {
      sig = sig + 1
    }
  }
  {
    hits: hits.copy(),
    lambda_gc: lambda,
    n_tested: hits.length(),
    n_genome_wide: sig,
    genome_wide_threshold: genome_wide,
  }
}

///|
/// 汇总 -> Manhattan 图数据(衔接 gwas_plot).
pub fn gwas_regression_manhattan(
  summary : GwasRegressionSummary,
) -> ManhattanData {
  gwas_manhattan(gwas_regression_markers(summary.hits))
}

///|
/// 汇总 -> QQ 图数据(衔接 gwas_plot).
pub fn gwas_regression_qq(summary : GwasRegressionSummary) -> QqData {
  gwas_qq(gwas_regression_markers(summary.hits))
}

// ============================================================================
// 导出
// ============================================================================

///|
fn gr_fmt(x : Double) -> String {
  x.to_string()
}

///|
/// 关联结果 -> TSV 文本.
pub fn gwas_regression_tsv(hits : Array[GwasRegressionHit]) -> String {
  let sb = StringBuilder()
  sb.write_string("chrom\tpos\tid\talt_freq\tbeta\tse\tz\tp_value\tn\n")
  for hit in hits {
    sb.write_string(hit.variant.chrom)
    sb.write_char('\t')
    sb.write_string(hit.variant.pos.to_string())
    sb.write_char('\t')
    sb.write_string(hit.variant.id)
    sb.write_char('\t')
    sb.write_string(gr_fmt(hit.alt_freq))
    sb.write_char('\t')
    sb.write_string(gr_fmt(hit.beta))
    sb.write_char('\t')
    sb.write_string(gr_fmt(hit.se))
    sb.write_char('\t')
    sb.write_string(gr_fmt(hit.z))
    sb.write_char('\t')
    sb.write_string(gr_fmt(hit.p_value))
    sb.write_char('\t')
    sb.write_string(hit.n_called.to_string())
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// 关联结果 -> JSON 值.
pub fn gwas_regression_value(summary : GwasRegressionSummary) -> JsonValue {
  let arr : Array[JsonValue] = Array(capacity=summary.hits.length())
  for hit in summary.hits {
    arr.push(
      JsonObj(
        Map([
          ("chrom", JsonStr(hit.variant.chrom)),
          ("pos", JsonNum(hit.variant.pos.to_double())),
          ("id", JsonStr(hit.variant.id)),
          ("alt_freq", JsonNum(hit.alt_freq)),
          ("beta", JsonNum(hit.beta)),
          ("se", JsonNum(hit.se)),
          ("z", JsonNum(hit.z)),
          ("p_value", JsonNum(hit.p_value)),
          ("n", JsonNum(hit.n_called.to_double())),
        ]),
      ),
    )
  }
  JsonObj(
    Map([
      ("lambda_gc", JsonNum(summary.lambda_gc)),
      ("n_tested", JsonNum(summary.n_tested.to_double())),
      ("n_genome_wide", JsonNum(summary.n_genome_wide.to_double())),
      ("genome_wide_threshold", JsonNum(summary.genome_wide_threshold)),
      ("hits", JsonArr(arr)),
    ]),
  )
}

///|
/// 关联结果 -> JSON 文本.
pub fn gwas_regression_json(summary : GwasRegressionSummary) -> String {
  json_to_string(gwas_regression_value(summary))
}

///|
/// 关联分析摘要行.
pub fn gwas_regression_summary_line(summary : GwasRegressionSummary) -> String {
  "gwas-regression: tested=\{summary.n_tested}, genome-wide=\{summary.n_genome_wide}, lambda_gc=\{summary.lambda_gc}"
}
