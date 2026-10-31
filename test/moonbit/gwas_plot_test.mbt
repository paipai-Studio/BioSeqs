// ===========================================================================
// gwas_plot_test.mbt — GWAS 可视化数据层测试(Manhattan/QQ/λGC/SVG)
// ===========================================================================

///|
/// 构造模拟 GWAS 摘要统计(均匀 p + 少数显著位点).
fn gwas_test_markers() -> Array[@src.GwasMarker] {
  let markers : Array[@src.GwasMarker] = []
  // chr1-3 各 100 个位点, 均匀 p
  for c = 1; c <= 3; c = c + 1 {
    for i = 0; i < 100; i = i + 1 {
      // 确定性伪均匀 p
      let u = ((i * 37 + c * 101) % 100).to_double() / 100.0
      let p = if u <= 0.0 { 0.01 } else { u }
      markers.push({ chrom: "chr" + c.to_string(), pos: i * 1000, p_value: p, })
    }
  }
  // 显著位点
  markers.push({ chrom: "chr2", pos: 50000, p_value: 5.0e-9, })
  markers.push({ chrom: "chr1", pos: 3000, p_value: 1.0e-10, })
  markers.push({ chrom: "chr3", pos: 80000, p_value: 2.0e-7, })
  markers
}

///|
test "gwas chrom rank ordering" {
  // 数字染色体按数值序
  let markers : Array[@src.GwasMarker] = [
    { chrom: "chr10", pos: 1, p_value: 0.5, },
    { chrom: "chr2", pos: 1, p_value: 0.5, },
    { chrom: "chrX", pos: 1, p_value: 0.5, },
    { chrom: "chr1", pos: 1, p_value: 0.5, },
  ]
  let data = @src.gwas_manhattan(markers)
  let first = data.points[0].chrom
  assert_eq(first, "chr1")
  assert_eq(data.points[1].chrom, "chr2")
  assert_eq(data.points[2].chrom, "chr10")
  assert_eq(data.points[3].chrom, "chrX")
}

///|
test "gwas manhattan basic structure" {
  let markers = gwas_test_markers()
  let data = @src.gwas_manhattan(markers)
  assert_eq(data.points.length(), 303)
  // 三条染色体边界
  assert_eq(data.chrom_boundaries.length(), 3)
  assert_eq(data.chrom_boundaries[0].0, "1")
  // 显著线: -log10(5e-8) ≈ 7.3
  assert_true(data.genome_wide_line > 7.2 && data.genome_wide_line < 7.4)
  // suggestive: -log10(1e-6) = 6
  assert_true(data.suggestive_line > 5.99 && data.suggestive_line < 6.01)
  // 最强位点 p=1e-10 -> -log10 = 10
  assert_true(data.y_max > 10.0)
}

///|
test "gwas manhattan genomic x monotonic" {
  let markers = gwas_test_markers()
  let data = @src.gwas_manhattan(markers)
  // 同染色体内 x 随位置单调; 染色体间递增
  let mut prev = -1.0
  for p in data.points {
    assert_true(p.genomic_x > prev || p.genomic_x >= 0.0)
    prev = p.genomic_x
  }
  // 染色体间 gap > 0
  let b1_end = data.chrom_boundaries[0].2
  let b2_start = data.chrom_boundaries[1].1
  assert_true(b2_start > b1_end)
}

///|
test "gwas manhattan significant flags and labels" {
  let markers = gwas_test_markers()
  let data = @src.gwas_manhattan(markers)
  // 2 个全基因组显著(5e-9 与 1e-10); 2e-7 仅 suggestive
  let mut n_sig = 0
  for p in data.points {
    if p.significant {
      n_sig = n_sig + 1
    }
  }
  assert_eq(n_sig, 2)
  // 标签含两条(默认 top 10, 此处 2 条)
  assert_eq(data.labels.length(), 2)
  // 最强标签是 p=1e-10 的 chr1:3000
  assert_eq(data.labels[0].0, "chr1")
  assert_eq(data.labels[0].1, 3000)
  assert_true(data.labels[0].2 > 9.99)
}

///|
test "gwas manhattan color alternation" {
  let markers = gwas_test_markers()
  let data = @src.gwas_manhattan(markers)
  // 染色体 1/3 同色, 2 交替
  let c1 = data.points[0].color_index
  let mut c2 = -1
  let mut c3 = -1
  for p in data.points {
    if p.chrom == "chr2" && c2 < 0 {
      c2 = p.color_index
    }
    if p.chrom == "chr3" && c3 < 0 {
      c3 = p.color_index
    }
  }
  assert_eq(c1, c3)
  assert_true(c1 != c2)
}

///|
test "gwas lambda gc calibrated null" {
  // 均匀 null p -> λGC 接近 1
  let markers : Array[@src.GwasMarker] = []
  for i = 0; i < 1000; i = i + 1 {
    let u = ((i * 61 + 7) % 1000).to_double() / 1000.0 + 0.0005
    markers.push({ chrom: "chr1", pos: i, p_value: u, })
  }
  let lambda = @src.gwas_lambda_gc(markers)
  assert_true(lambda > 0.8 && lambda < 1.2)
  let interp = @src.gwas_lambda_interpretation(lambda)
  assert_true(interp.contains("well-calibrated"))
}

///|
test "gwas lambda gc inflated signal" {
  // 全部小 p -> 膨胀
  let markers : Array[@src.GwasMarker] = []
  for i = 0; i < 100; i = i + 1 {
    markers.push({
      chrom: "chr1",
      pos: i,
      p_value: 0.001 + i.to_double() * 0.0001,
    })
  }
  let lambda = @src.gwas_lambda_gc(markers)
  assert_true(lambda > 1.5)
  let interp = @src.gwas_lambda_interpretation(lambda)
  assert_true(interp.contains("inflated"))
}

///|
test "gwas lambda gc empty" {
  assert_eq(@src.gwas_lambda_gc([]), 1.0)
}

///|
test "gwas qq basic structure" {
  let markers = gwas_test_markers()
  let qq = @src.gwas_qq(markers)
  assert_eq(qq.points.length(), 303)
  // null 大部分点位贴近 y=x 线(距离容差, 3 个显著位点除外)
  let mut on_diag = 0
  for p in qq.points {
    let d = p.observed - p.expected
    let ad = if d >= 0.0 { d } else { -d }
    if ad < 0.5 {
      on_diag = on_diag + 1
    }
  }
  assert_true(on_diag > 250)
  // λGC 已计算
  assert_true(qq.lambda_gc > 0.0)
}

///|
test "gwas qq sorted observed" {
  let markers = gwas_test_markers()
  let qq = @src.gwas_qq(markers)
  // observed 单调不减
  let mut prev = -1.0
  for p in qq.points {
    assert_true(p.observed >= prev)
    prev = p.observed
  }
  // expected 单调递增
  let mut prev_e = -1.0
  for p in qq.points {
    assert_true(p.expected > prev_e)
    prev_e = p.expected
  }
}

///|
test "gwas manhattan svg output" {
  let markers = gwas_test_markers()
  let data = @src.gwas_manhattan(markers)
  let svg = @src.gwas_manhattan_svg(data)
  assert_true(svg.has_prefix("<svg"))
  assert_true(svg.contains("Manhattan Plot"))
  assert_true(svg.contains("<circle"))
  assert_true(svg.contains("stroke=\"red\""))
  assert_true(svg.contains("</svg>"))
  // 染色体标签
  assert_true(svg.contains(">1<"))
  assert_true(svg.contains(">2<"))
}

///|
test "gwas qq svg output" {
  let markers = gwas_test_markers()
  let qq = @src.gwas_qq(markers)
  let svg = @src.gwas_qq_svg(qq)
  assert_true(svg.has_prefix("<svg"))
  assert_true(svg.contains("QQ Plot"))
  assert_true(svg.contains("λGC"))
  assert_true(svg.contains("<circle"))
  assert_true(svg.contains("</svg>"))
}

///|
test "gwas manhattan empty markers" {
  let data = @src.gwas_manhattan([])
  assert_eq(data.points.length(), 0)
  assert_eq(data.chrom_boundaries.length(), 0)
  assert_eq(data.labels.length(), 0)
  let svg = @src.gwas_manhattan_svg(data)
  assert_true(svg.has_prefix("<svg"))
}

///|
test "gwas qq empty markers" {
  let qq = @src.gwas_qq([])
  assert_eq(qq.points.length(), 0)
  assert_eq(qq.lambda_gc, 1.0)
}

///|
test "gwas lambda interpretation boundaries" {
  assert_true(@src.gwas_lambda_interpretation(0.85).contains("deflated"))
  assert_true(@src.gwas_lambda_interpretation(0.95).contains("well-calibrated"))
  assert_true(@src.gwas_lambda_interpretation(1.05).contains("well-calibrated"))
  assert_true(@src.gwas_lambda_interpretation(1.2).contains("inflated"))
}

///|
test "gwas manhattan single chromosome" {
  let markers : Array[@src.GwasMarker] = [
    { chrom: "chr1", pos: 100, p_value: 0.5, },
    { chrom: "chr1", pos: 200, p_value: 0.01, },
  ]
  let data = @src.gwas_manhattan(markers)
  assert_eq(data.chrom_boundaries.length(), 1)
  assert_eq(data.points.length(), 2)
  assert_eq(data.points[0].pos, 100)
  assert_eq(data.points[1].pos, 200)
  // 无显著
  assert_eq(data.labels.length(), 0)
}
