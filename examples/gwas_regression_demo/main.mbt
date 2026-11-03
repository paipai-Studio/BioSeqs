///|
/// GWAS 回归关联分析闭环演示:
/// 合成队列(基因型 + 二值表型) -> 逐位点逻辑回归 -> λGC -> TSV/JSON ->
/// Manhattan/QQ 数据与 SVG.

///|
fn gw_demo_variants(n : Int32) -> Array[@src.GwasVariantMeta] {
  let out : Array[@src.GwasVariantMeta] = Array(capacity=n)
  for v = 0; v < n; v = v + 1 {
    let chrom = if v < n / 2 { "chr1" } else { "chr2" }
    out.push(@src.GwasVariantMeta::{
      chrom,
      pos: 1000 * (v + 1),
      id: "rs" + v.to_string(),
    })
  }
  out
}

///|
fn gw_demo_dosages(
  rng : @src.Rng,
  n : Int32,
  n_var : Int32,
) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=n_var)
    for _v = 0; _v < n_var; _v = _v + 1 {
      let a1 = if rng.randint(100) < 30 { 1 } else { 0 }
      let a2 = if rng.randint(100) < 30 { 1 } else { 0 }
      row.push((a1 + a2).to_double())
    }
    out.push(row)
  }
  out
}

///|
fn gw_demo_phenotype(
  rng : @src.Rng,
  dosages : Array[Array[Double]],
  causal : Int32,
) -> Array[Double] {
  let y : Array[Double] = Array(capacity=dosages.length())
  for i = 0; i < dosages.length(); i = i + 1 {
    let d = dosages[i][causal]
    let eta = -0.4 + 0.9 * d
    let p = 1.0 / (1.0 + @math.exp(-eta))
    y.push(if rng.next_double() < p { 1.0 } else { 0.0 })
  }
  y
}

///|
fn main {
  println("=== GWAS 回归关联分析闭环 ===")
  let n = 200
  let n_var = 20
  let causal = 5
  let rng = @src.Rng::new(20261030UL)
  let variants = gw_demo_variants(n_var)
  let dosages = gw_demo_dosages(rng, n, n_var)
  let phenotype = gw_demo_phenotype(rng, dosages, causal)
  println("样本数: \{n}, 位点数: \{n_var}, 因果位点: rs\{causal}")

  println("\n1. 逐位点逻辑回归")
  let hits = @src.gwas_regression_logistic(
    dosages,
    phenotype,
    variants,
    min_af=0.05,
  )
  let sorted = hits.copy()
  sorted.sort_by(fn(a, b) {
    if a.p_value < b.p_value {
      -1
    } else if a.p_value > b.p_value {
      1
    } else {
      0
    }
  })
  let shown = if sorted.length() < 5 { sorted.length() } else { 5 }
  for i = 0; i < shown; i = i + 1 {
    let hit = sorted[i]
    println(
      "  \{hit.variant.id} (\{hit.variant.chrom}:\{hit.variant.pos}) beta=\{hit.beta} se=\{hit.se} z=\{hit.z} p=\{hit.p_value}",
    )
  }

  println("\n2. 汇总(λGC / 显著位点)")
  let summary = @src.gwas_regression_summary(hits)
  println("  \{@src.gwas_regression_summary_line(summary)}")

  println("\n3. TSV 导出(前 3 行)")
  let tsv = @src.gwas_regression_tsv(hits)
  let lines = tsv.split("\n").to_array()
  let tshown = if lines.length() < 4 { lines.length() } else { 4 }
  for i = 0; i < tshown; i = i + 1 {
    println("  " + lines[i].to_owned())
  }

  println("\n4. JSON 导出")
  let json = @src.gwas_regression_json(summary)
  println("  JSON 长度 = \{json.length()} 字节")

  println("\n5. 可视化数据与 SVG")
  let manhattan = @src.gwas_regression_manhattan(summary)
  let qq = @src.gwas_regression_qq(summary)
  println("  Manhattan 点数: \{manhattan.points.length()}")
  println("  QQ 点数: \{qq.points.length()}, lambda_gc=\{qq.lambda_gc}")
  let manhattan_svg = @src.gwas_manhattan_svg(manhattan)
  let qq_svg = @src.gwas_qq_svg(qq)
  println("  Manhattan SVG 长度 = \{manhattan_svg.length()} 字节")
  println("  QQ SVG 长度 = \{qq_svg.length()} 字节")

  println(
    "\n演示完成: 基因型 + 表型 -> 回归关联 -> 膨胀因子 -> 导出 -> 可视化闭环",
  )
}
