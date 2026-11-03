///|
/// GWAS 可视化演示: Manhattan/QQ/λGC + 内嵌 SVG.
fn main {
  println("=== GWAS 可视化数据层 (Manhattan / QQ / λGC) ===")

  // 模拟 GWAS 摘要统计: 3 染色体 x 300 null 位点 + 3 显著位点
  let markers : Array[@src.GwasMarker] = []
  for c = 1; c <= 3; c = c + 1 {
    for i = 0; i < 100; i = i + 1 {
      let u = ((i * 37 + c * 101) % 100).to_double() / 100.0
      let p = if u <= 0.0 { 0.01 } else { u }
      markers.push({ chrom: "chr" + c.to_string(), pos: i * 1000, p_value: p, })
    }
  }
  markers.push({ chrom: "chr2", pos: 50000, p_value: 5.0e-9, })
  markers.push({ chrom: "chr1", pos: 3000, p_value: 1.0e-10, })
  markers.push({ chrom: "chr3", pos: 80000, p_value: 2.0e-7, })
  println("位 点 数: \{markers.length()} (303 null + 3 hit)")

  // 1. Manhattan 数据
  println("\n1. Manhattan 图数据")
  let manhattan = @src.gwas_manhattan(markers)
  println("  数据点: \{manhattan.points.length()}")
  println("  染色体边界: \{manhattan.chrom_boundaries.length()} 条")
  for b in manhattan.chrom_boundaries {
    println("    chr\{b.0}: x ∈ [\{b.1}, \{b.2}]")
  }
  println("  y_max: \{manhattan.y_max}")
  println(
    "  全基因组显著线: -log10(5e-8) = \{manhattan.genome_wide_line}",
  )
  println("  suggestive 线: -log10(1e-6) = \{manhattan.suggestive_line}")
  println("  显著标签: \{manhattan.labels.length()} 条")
  for l in manhattan.labels {
    println("    \{l.0}:\{l.1.to_string()} (-log10 p = \{l.2})")
  }

  // 2. λGC
  println("\n2. λGC 基因组膨胀因子")
  let lambda = @src.gwas_lambda_gc(markers)
  println("  λGC = \{lambda}")
  println("  解读: \{@src.gwas_lambda_interpretation(lambda)}")

  // 3. QQ 数据
  println("\n3. QQ plot 数据")
  let qq = @src.gwas_qq(markers)
  println("  配对点: \{qq.points.length()}")
  println("  λGC(内嵌): \{qq.lambda_gc}")
  println("  前 3 点(expected -> observed):")
  for i = 0; i < 3; i = i + 1 {
    let p = qq.points[i]
    println("    \{p.expected} -> \{p.observed} (带 [\{p.lower}, \{p.upper}])")
  }

  // 4. SVG 输出
  println("\n4. SVG 渲染(独立可打开的 .svg 文本)")
  let svg_m = @src.gwas_manhattan_svg(manhattan)
  let svg_q = @src.gwas_qq_svg(qq)
  println("  Manhattan SVG: \{svg_m.length()} 字符")
  println("  QQ SVG: \{svg_q.length()} 字符")
  println("  Manhattan 含散点: \{svg_m.contains("<circle")}")
  println("  QQ 含置信带: \{svg_q.contains("<path")}")
  println("  QQ 含 λGC 标注: \{svg_q.contains("λGC")}")

  println(
    "\n完成: Manhattan/QQ/λGC 数据层+SVG 渲染全链路本地完成.",
  )
}
