///|
/// gwas_plot.mbt — GWAS 可视化数据层(Manhattan / QQ / λGC + 内嵌 SVG 渲染).
///
/// 对齐 qqman / CMplot 的核心输出语义(本地无依赖):
/// - Manhattan: 染色体交替着色散点数据(按基因组位置排序, -log10 p),
///   全基因组与suggestive 两条显著线
/// - QQ plot: 观测 vs 期望 -log10 p 分位数配对, 置信带近似
/// - λGC: 中位数卡方膨胀因子 χ²median/0.4549, >1 膨胀 <1 收缩
/// - 内嵌 SVG 渲染: 散点/折线/轴/标签, 输出独立可打开的 .svg 文本
///
/// 输入为 (chrom, pos, p) 数组, 衔接 gwas_chisq_assoc 等现有统计.
/// 数据不出本机.
pub(all) struct GwasMarker {
  chrom : String
  pos : Int32
  p_value : Double
}

///|
/// Manhattan 图数据点.
pub(all) struct ManhattanPoint {
  chrom : String
  pos : Int32
  /// 累积基因组 x 坐标(bp)
  genomic_x : Double
  neg_log10_p : Double
  color_index : Int32
  significant : Bool
  suggestive : Bool
}

///|
/// Manhattan 图数据.
pub(all) struct ManhattanData {
  points : Array[ManhattanPoint]
  chrom_boundaries : Array[(String, Double, Double)]
  x_max : Double
  y_max : Double
  genome_wide_line : Double
  suggestive_line : Double
  /// 显著位点标签(前 top_hits 条)
  labels : Array[(String, Int32, Double)]
}

///|
/// QQ plot 数据点.
pub(all) struct QqPoint {
  expected : Double
  observed : Double
  lower : Double
  upper : Double
}

///|
/// QQ plot 数据.
pub(all) struct QqData {
  points : Array[QqPoint]
  lambda_gc : Double
  x_max : Double
  y_max : Double
}

///|
/// 染色体排序权重(数字优先, X/Y/M 其后).
fn gwas_chrom_rank(chrom : String) -> Int32 {
  if chrom.length() == 0 {
    return 999
  }
  // 去掉 chr 前缀后按本体排序
  let base : String = if chrom.length() > 3 && chrom.has_prefix("chr") {
    chrom[3:].to_owned()
  } else {
    chrom
  }
  if base.length() == 0 {
    return 999
  }
  let c = base[0].unsafe_to_char()
  if c >= '0' && c <= '9' {
    let n = parse_int(base)
    if n >= 1 && n <= 22 {
      return n
    }
    return 100 + n % 100
  }
  match base {
    "X" | "x" => 23
    "Y" | "y" => 24
    "MT" | "M" => 25
    _ => 200
  }
}

///|
/// 染色体显示名(去掉 chr 前缀).
fn gwas_chrom_label(chrom : String) -> String {
  if chrom.length() > 3 && chrom.has_prefix("chr") {
    chrom[3:].to_owned()
  } else {
    chrom
  }
}

///|
fn gwas_neg_log10(p : Double) -> Double {
  if p <= 0.0 {
    300.0
  } else if p >= 1.0 {
    0.0
  } else {
    -@math.ln(p) / @math.ln(10.0)
  }
}

///|
/// 构建 Manhattan 图数据.
/// genome_wide 默认 5e-8, suggestive 默认 1e-6; top_hits 条显著位点带标签.
pub fn gwas_manhattan(
  markers : Array[GwasMarker],
  genome_wide~ : Double = 5.0e-8,
  suggestive~ : Double = 1.0e-6,
  top_hits~ : Int32 = 10,
) -> ManhattanData {
  // 按染色体序与位置稳定排序
  let sorted : Array[GwasMarker] = Array(capacity=0)
  for m in markers {
    sorted.push(m)
  }
  sorted.sort_by(fn(a, b) {
    let ra = gwas_chrom_rank(a.chrom)
    let rb = gwas_chrom_rank(b.chrom)
    if ra != rb {
      ra - rb
    } else if a.pos != b.pos {
      a.pos - b.pos
    } else {
      compare(a.chrom, b.chrom)
    }
  })
  // 累积 x 坐标(染色体间留 gap)
  let gap : Double = 2.0e7 * 3.0
  let points : Array[ManhattanPoint] = Array(capacity=0)
  let boundaries : Array[(String, Double, Double)] = Array(capacity=0)
  let mut x : Double = 0.0
  let mut i = 0
  let n = sorted.length()
  while i < n {
    let chrom = sorted[i].chrom
    let chrom_start = x
    let mut j = i
    while j < n && sorted[j].chrom == chrom {
      let m = sorted[j]
      let px = x + m.pos.to_double()
      let nlp = gwas_neg_log10(m.p_value)
      points.push({
        chrom,
        pos: m.pos,
        genomic_x: px,
        neg_log10_p: nlp,
        color_index: gwas_chrom_rank(chrom) % 2,
        significant: m.p_value < genome_wide,
        suggestive: m.p_value < suggestive,
      })
      j = j + 1
    }
    // 该染色体长度跨度
    let span = if j > i { sorted[j - 1].pos.to_double() } else { 0.0 }
    let chrom_end = chrom_start + span
    boundaries.push((gwas_chrom_label(chrom), chrom_start, chrom_end))
    x = chrom_end + gap
    i = j
  }
  let mut y_max = 0.0
  for p in points {
    if p.neg_log10_p > y_max {
      y_max = p.neg_log10_p
    }
  }
  if y_max <= 0.0 {
    y_max = 1.0
  }
  y_max = y_max * 1.1 + 1.0
  // 显著标签(取最强 top_hits)
  let sig_points : Array[ManhattanPoint] = Array(capacity=0)
  for p in points {
    if p.significant {
      sig_points.push(p)
    }
  }
  sig_points.sort_by(fn(a, b) {
    if a.neg_log10_p > b.neg_log10_p {
      -1
    } else if a.neg_log10_p < b.neg_log10_p {
      1
    } else {
      compare(a.chrom, b.chrom)
    }
  })
  let labels : Array[(String, Int32, Double)] = Array(capacity=0)
  let limit = if sig_points.length() < top_hits {
    sig_points.length()
  } else {
    top_hits
  }
  for k = 0; k < limit; k = k + 1 {
    labels.push(
      (sig_points[k].chrom, sig_points[k].pos, sig_points[k].neg_log10_p),
    )
  }
  {
    points,
    chrom_boundaries: boundaries,
    x_max: x,
    y_max,
    genome_wide_line: gwas_neg_log10(genome_wide),
    suggestive_line: gwas_neg_log10(suggestive),
    labels,
  }
}

///|
/// λGC 基因组膨胀因子: χ²median / 0.4549(χ²1 中位数的卡方).
/// 期望 1.0; >1 通货膨胀(群体分层等), <1 收缩.
pub fn gwas_lambda_gc(markers : Array[GwasMarker]) -> Double {
  let n = markers.length()
  if n == 0 {
    return 1.0
  }
  let pvals : Array[Double] = Array(capacity=n)
  for m in markers {
    pvals.push(m.p_value)
  }
  pvals.sort()
  let median_p = if n % 2 == 1 {
    pvals[n / 2]
  } else {
    (pvals[n / 2 - 1] + pvals[n / 2]) / 2.0
  }
  if median_p <= 0.0 || median_p >= 1.0 {
    return 1.0
  }
  // 上尾语义: P(χ²1 > x) = median_p -> Φ(√x) = 1 - median_p/2
  // -> x = (probit(1 - median_p/2))²
  let z = gwas_probit(1.0 - median_p / 2.0)
  let chi2_med = z * z
  chi2_med / 0.4549364
}

///|
/// 标准正态分位数函数 probit(Beasley-Springer-Moro 简化).
fn gwas_probit(p : Double) -> Double {
  if p <= 0.0 {
    return -5.0
  }
  if p >= 1.0 {
    return 5.0
  }
  // 有理近似(Acklam 算法简化系数)
  let a = [
    -3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02,
    -3.066479806614716e+01, 2.506628277459239e+00,
  ]
  let b = [
    -5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01,
    -1.328068155288572e+01,
  ]
  let c = [
    -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00,
    4.374664141464968e+00, 2.938163982698783e+00,
  ]
  let d = [
    7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00,
  ]
  let p_low = 0.02425
  if p < p_low {
    let q = (-2.0 * @math.ln(p)).sqrt()
    let mut r = (
        ((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]
      ) /
      ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    if r < -5.0 {
      r = -5.0
    }
    return r
  }
  if p > 1.0 - p_low {
    let q = (-2.0 * @math.ln(1.0 - p)).sqrt()
    let mut r = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q +
      c[5]) /
      ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    if r > 5.0 {
      r = 5.0
    }
    return r
  }
  let q = p - 0.5
  let r = q * q
  (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) *
  q /
  (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
}

///|
/// 构建 QQ plot 数据(观测 vs 期望分位配对, 95% 置信带).
pub fn gwas_qq(markers : Array[GwasMarker]) -> QqData {
  let n = markers.length()
  if n == 0 {
    return { points: [], lambda_gc: 1.0, x_max: 1.0, y_max: 1.0, }
  }
  let pvals : Array[Double] = Array(capacity=n)
  for m in markers {
    pvals.push(m.p_value)
  }
  pvals.sort()
  // 降序配对(最大 p 在前 -> observed 从小到大, 与 expected 同向)
  let points : Array[QqPoint] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    // 期望分位(与 observed 同向: -log10 递增 -> p 递减)
    let expected_p = (n.to_double() - i.to_double() - 0.5) / n.to_double()
    let observed_p = pvals[n - 1 - i]
    // 95% 置信带: Beta(i+1, n-i) 分位近似 -> 用正态近似
    let se = gwas_qq_se(expected_p, n)
    points.push({
      expected: gwas_neg_log10(expected_p),
      observed: gwas_neg_log10(observed_p),
      lower: gwas_neg_log10(expected_p + 1.96 * se),
      upper: gwas_neg_log10(expected_p - 1.96 * se),
    })
  }
  let mut x_max = 0.0
  let mut y_max = 0.0
  for p in points {
    if p.expected > x_max {
      x_max = p.expected
    }
    if p.observed > y_max {
      y_max = y_max + 0.0
      if p.observed > y_max {
        y_max = p.observed
      }
    }
  }
  let lim = if x_max > y_max { x_max } else { y_max }
  if lim <= 0.0 {
    return {
      points,
      lambda_gc: gwas_lambda_gc(markers),
      x_max: 1.0,
      y_max: 1.0,
    }
  }
  {
    points,
    lambda_gc: gwas_lambda_gc(markers),
    x_max: lim * 1.05,
    y_max: lim * 1.05,
  }
}

///|
fn gwas_qq_se(p : Double, n : Int32) -> Double {
  // -log10 p 的正态近似标准误(dp/dp 分位方差)
  if p <= 0.0 || p >= 1.0 || n <= 0 {
    return 0.0
  }
  (p * (1.0 - p) / n.to_double()).sqrt()
}

// ---------------------------------------------------------------------------
// SVG 渲染
// ---------------------------------------------------------------------------

///|
fn gwas_svg_header(w : Int32, h : Int32, title : String) -> String {
  "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"\{w}\" height=\"\{h}\" viewBox=\"0 0 \{w} \{h}\">\n<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n<text x=\"\{w / 2}\" y=\"24\" text-anchor=\"middle\" font-size=\"14\" font-weight=\"bold\">\{title}</text>\n"
}

///|
fn gwas_fmt(x : Double) -> String {
  if x == x.trunc() && x > -1.0e12 && x < 1.0e12 {
    x.to_int().to_string()
  } else {
    let r = (x * 100.0).round() / 100.0
    r.to_string()
  }
}

///|
/// Manhattan 图 SVG 渲染.
pub fn gwas_manhattan_svg(
  data : ManhattanData,
  width~ : Int32 = 900,
  height~ : Int32 = 400,
) -> String {
  let sb = StringBuilder()
  sb.write_string(gwas_svg_header(width, height, "Manhattan Plot"))
  let left = 60.0
  let right = width.to_double() - 20.0
  let top = 40.0
  let bottom = height.to_double() - 40.0
  let plot_w = right - left
  let plot_h = bottom - top
  let x_scale = if data.x_max > 0.0 { plot_w / data.x_max } else { 1.0 }
  let y_scale = if data.y_max > 0.0 { plot_h / data.y_max } else { 1.0 }
  // 轴
  sb.write_string(
    "<line x1=\"\{gwas_fmt(left)}\" y1=\"\{gwas_fmt(bottom)}\" x2=\"\{gwas_fmt(right)}\" y2=\"\{gwas_fmt(bottom)}\" stroke=\"black\"/>\n",
  )
  sb.write_string(
    "<line x1=\"\{gwas_fmt(left)}\" y1=\"\{gwas_fmt(top)}\" x2=\"\{gwas_fmt(left)}\" y2=\"\{gwas_fmt(bottom)}\" stroke=\"black\"/>\n",
  )
  // y 轴刻度
  let y_ticks : Int32 = 5
  for i = 1; i <= y_ticks; i = i + 1 {
    let v = data.y_max * i.to_double() / y_ticks.to_double()
    let y = bottom - v * y_scale
    sb.write_string(
      "<line x1=\"\{gwas_fmt(left - 4)}\" y1=\"\{gwas_fmt(y)}\" x2=\"\{gwas_fmt(left)}\" y2=\"\{gwas_fmt(y)}\" stroke=\"black\"/>\n",
    )
    sb.write_string(
      "<text x=\"\{gwas_fmt(left - 8)}\" y=\"\{gwas_fmt(y + 4)}\" text-anchor=\"end\" font-size=\"10\">\{gwas_fmt(v)}</text>\n",
    )
  }
  sb.write_string(
    "<text x=\"16\" y=\"\{gwas_fmt((top + bottom) / 2.0)}\" font-size=\"11\" transform=\"rotate(-90 16 \{gwas_fmt((top + bottom) / 2.0)})\" text-anchor=\"middle\">-log10(p)</text>\n",
  )
  // 显著线
  let gw_y = bottom - data.genome_wide_line * y_scale
  if data.genome_wide_line < data.y_max {
    sb.write_string(
      "<line x1=\"\{gwas_fmt(left)}\" y1=\"\{gwas_fmt(gw_y)}\" x2=\"\{gwas_fmt(right)}\" y2=\"\{gwas_fmt(gw_y)}\" stroke=\"red\" stroke-dasharray=\"4,2\"/>\n",
    )
  }
  let sg_y = bottom - data.suggestive_line * y_scale
  if data.suggestive_line < data.y_max &&
    data.suggestive_line != data.genome_wide_line {
    sb.write_string(
      "<line x1=\"\{gwas_fmt(left)}\" y1=\"\{gwas_fmt(sg_y)}\" x2=\"\{gwas_fmt(right)}\" y2=\"\{gwas_fmt(sg_y)}\" stroke=\"blue\" stroke-dasharray=\"2,3\"/>\n",
    )
  }
  // 散点(染色体交替蓝/绿)
  for p in data.points {
    let x = left + p.genomic_x * x_scale
    let y = bottom - p.neg_log10_p * y_scale
    if x < left || x > right || y < top {
      continue
    }
    let color = if p.color_index == 0 { "#2c7fb8" } else { "#41ab5d" }
    let r = if p.significant { 3.0 } else { 2.0 }
    sb.write_string(
      "<circle cx=\"\{gwas_fmt(x)}\" cy=\"\{gwas_fmt(y)}\" r=\"\{gwas_fmt(r)}\" fill=\"\{color}\"/>\n",
    )
  }
  // 染色体标签(区域中心)
  for b in data.chrom_boundaries {
    let cx = left + (b.1 + b.2) / 2.0 * x_scale
    if cx >= left && cx <= right {
      sb.write_string(
        "<text x=\"\{gwas_fmt(cx)}\" y=\"\{gwas_fmt(bottom + 16)}\" text-anchor=\"middle\" font-size=\"9\">\{b.0}</text>\n",
      )
    }
  }
  // 显著标签
  for l in data.labels {
    // 找到点位置
    for p in data.points {
      if p.chrom == l.0 && p.pos == l.1 && p.neg_log10_p == l.2 {
        let x = left + p.genomic_x * x_scale
        let y = bottom - p.neg_log10_p * y_scale
        sb.write_string(
          "<text x=\"\{gwas_fmt(x)}\" y=\"\{gwas_fmt(y - 6)}\" text-anchor=\"middle\" font-size=\"9\" font-weight=\"bold\">\{gwas_chrom_label(l.0)}:\{l.1.to_string()}</text>\n",
        )
        break
      }
    }
  }
  sb.write_string("</svg>\n")
  sb.to_string()
}

///|
/// QQ plot SVG 渲染(含 y=x 参考线与置信带).
pub fn gwas_qq_svg(
  data : QqData,
  width~ : Int32 = 500,
  height~ : Int32 = 500,
) -> String {
  let sb = StringBuilder()
  sb.write_string(gwas_svg_header(width, height, "QQ Plot"))
  let left = 60.0
  let right = width.to_double() - 20.0
  let top = 40.0
  let bottom = height.to_double() - 50.0
  let plot_w = right - left
  let plot_h = bottom - top
  let x_scale = if data.x_max > 0.0 { plot_w / data.x_max } else { 1.0 }
  let y_scale = if data.y_max > 0.0 { plot_h / data.y_max } else { 1.0 }
  // 轴
  sb.write_string(
    "<line x1=\"\{gwas_fmt(left)}\" y1=\"\{gwas_fmt(bottom)}\" x2=\"\{gwas_fmt(right)}\" y2=\"\{gwas_fmt(bottom)}\" stroke=\"black\"/>\n",
  )
  sb.write_string(
    "<line x1=\"\{gwas_fmt(left)}\" y1=\"\{gwas_fmt(top)}\" x2=\"\{gwas_fmt(left)}\" y2=\"\{gwas_fmt(bottom)}\" stroke=\"black\"/>\n",
  )
  // y=x 参考线
  let diag_end_x = left + data.x_max * x_scale
  let diag_end_y = bottom - data.x_max * y_scale
  sb.write_string(
    "<line x1=\"\{gwas_fmt(left)}\" y1=\"\{gwas_fmt(bottom)}\" x2=\"\{gwas_fmt(diag_end_x)}\" y2=\"\{gwas_fmt(diag_end_y)}\" stroke=\"gray\" stroke-dasharray=\"4,3\"/>\n",
  )
  // 置信带(上界折线)
  let mut band_path = "M"
  let mut first = true
  for p in data.points {
    let x = left + p.expected * x_scale
    let yu = bottom - p.upper * y_scale
    if yu < top || yu > bottom {
      continue
    }
    if first {
      band_path = band_path + " \{gwas_fmt(x)} \{gwas_fmt(yu)}"
      first = false
    } else {
      band_path = band_path + " L \{gwas_fmt(x)} \{gwas_fmt(yu)}"
    }
  }
  // 下界反向
  let mut k = data.points.length() - 1
  while k >= 0 {
    let p = data.points[k]
    let x = left + p.expected * x_scale
    let yl = bottom - p.lower * y_scale
    if yl >= top && yl <= bottom {
      band_path = band_path + " L \{gwas_fmt(x)} \{gwas_fmt(yl)}"
    }
    k = k - 1
  }
  band_path = band_path + " Z"
  sb.write_string(
    "<path d=\"\{band_path}\" fill=\"#41ab5d\" fill-opacity=\"0.2\" stroke=\"none\"/>\n",
  )
  // 数据点
  for p in data.points {
    let x = left + p.expected * x_scale
    let y = bottom - p.observed * y_scale
    if x < left || x > right || y < top || y > bottom {
      continue
    }
    sb.write_string(
      "<circle cx=\"\{gwas_fmt(x)}\" cy=\"\{gwas_fmt(y)}\" r=\"2\" fill=\"#2c7fb8\"/>\n",
    )
  }
  // λGC 标注
  sb.write_string(
    "<text x=\"\{gwas_fmt(right - 8)}\" y=\"\{gwas_fmt(top + 14)}\" text-anchor=\"end\" font-size=\"11\">λGC = \{gwas_fmt(data.lambda_gc)}</text>\n",
  )
  // 轴标签
  sb.write_string(
    "<text x=\"\{gwas_fmt((left + right) / 2.0)}\" y=\"\{gwas_fmt(bottom + 30)}\" text-anchor=\"middle\" font-size=\"11\">Expected -log10(p)</text>\n",
  )
  sb.write_string(
    "<text x=\"16\" y=\"\{gwas_fmt((top + bottom) / 2.0)}\" font-size=\"11\" transform=\"rotate(-90 16 \{gwas_fmt((top + bottom) / 2.0)})\" text-anchor=\"middle\">Observed -log10(p)</text>\n",
  )
  sb.write_string("</svg>\n")
  sb.to_string()
}

///|
/// 基因组膨胀解读.
pub fn gwas_lambda_interpretation(lambda : Double) -> String {
  if lambda > 1.1 {
    "inflated (λGC > 1.1): 考虑群体分层校正(PCA/LMM)"
  } else if lambda < 0.9 {
    "deflated (λGC < 0.9): 检验保守或过校正"
  } else {
    "well-calibrated (0.9 <= λGC <= 1.1)"
  }
}
