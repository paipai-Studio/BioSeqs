// ============================================================================
// MDAnalysis 风格：Contact Map 概率（多帧统计）
// ============================================================================
//
// 给定 AtomGroup 的全部轨迹帧：
//   1. 每帧计算 cutoff 内的残基-残基原子对
//   2. 按残基汇总（同一残基内接触忽略）
//   3. 跨帧统计每个残基对接触的频率 → 接触概率矩阵
//
// 返回 R × R 的 Double 矩阵 + resid / resname 标签数组。
// ============================================================================

///|
/// 残基接触概率结果。
pub struct ContactMapResult {
  /// 残基总数 R
  n_res : Int32
  /// 残基 ID 数组（R 个）
  resids : Array[Int32]
  /// 残基名数组（R 个）
  resnames : Array[String]
  /// 残基在 AtomGroup 内的起始原子索引（R 个）
  res_start_idx : Array[Int32]
  /// 接触概率矩阵，shape `(R, R)`；对角元 = 0。
  prob : Array[Array[Double]]
  /// 每帧平均接触数（R²/2）随时间演化（可选统计）。
  frame_contact_counts : Array[Int32]
  /// 总帧数
  n_frames : Int32
}

///|
/// 从 AtomGroup 的轨迹统计残基-残基接触概率。
///
/// 参数：
///   - `cutoff~` 接触判据（默认 4.5 Å，与 MDAnalysis Contacts 一致）
///   - `res_per_atom~` false → 先按残基名/ID 自动识别残基边界；
///                     true  → 每个原子视为独立"残基"（原子级接触图）
///   - `pbc~` 是否用 PBC 最小镜像
pub fn AtomGroup::contact_map_probability(
  self : AtomGroup,
  cutoff~ : Double = 4.5,
  res_per_atom~ : Bool = false,
  pbc~ : Bool = false,
) -> ContactMapResult {
  let n_atoms = self.indices.length()
  let n_frames = self.universe.n_frames()
  let atoms = self.universe.atoms
  let zero = ContactMapResult::{
    n_res: 0,
    resids: Array(capacity=0),
    resnames: Array(capacity=0),
    res_start_idx: Array(capacity=0),
    prob: Array(capacity=0),
    frame_contact_counts: Array(capacity=0),
    n_frames,
  }
  if n_atoms < 2 || n_frames < 1 { return zero }

  // ===== 1. 建立残基索引 =====
  let resids : Array[Int32] = Array(capacity=n_atoms)
  let resnames : Array[String] = Array(capacity=n_atoms)
  let res_atom_idx : Array[Int32] = Array(capacity=n_atoms)
  let mut n_res = 0
  if res_per_atom {
    let mut i = 0
    while i < n_atoms {
      let idx = self.indices[i]
      resids.push(i)
      resnames.push(
        if idx < atoms.length() { atoms[idx].name } else { "" }
      )
      res_atom_idx.push(i)
      i = i + 1
    }
    n_res = n_atoms
  } else {
    // 按 resid + resname 边界切分
    let mut cur_resid = -1
    let mut cur_resname = ""
    let mut i = 0
    while i < n_atoms {
      let idx = self.indices[i]
      let rid = if idx < atoms.length() { atoms[idx].resid } else { -1 }
      let rn = if idx < atoms.length() { atoms[idx].resname } else { "" }
      let is_new = (rid != cur_resid) || (rn != cur_resname)
      if is_new {
        cur_resid = rid
        cur_resname = rn
        resids.push(n_res)
        resnames.push(rn)
        n_res = n_res + 1
      }
      res_atom_idx.push(n_res - 1)
      i = i + 1
    }
  }
  if n_res < 2 { return zero }

  let cutoff2 = cutoff * cutoff
  let box = if pbc { Some(self.universe.box) } else { None }

  // ===== 2. 接触计数累加 =====
  let count : Array[Array[Int32]] = Array(capacity=n_res)
  let mut i = 0
  while i < n_res {
    let row : Array[Int32] = Array(capacity=n_res)
    let mut j = 0
    while j < n_res { row.push(0); j = j + 1 }
    count.push(row)
    i = i + 1
  }
  let frame_contact_counts : Array[Int32] = Array(capacity=n_frames)
  let mut f = 0
  while f < n_frames {
    let fr = self.universe.frame(f)
    let mut f_count = 0
    // 逐原子对判断；同残基跳过
    i = 0
    while i < n_atoms {
      let ri = res_atom_idx[i]
      let idx_i = self.indices[i]
      let ci = if idx_i < fr.length() { fr[idx_i] } else { Vector3::new(0.0, 0.0, 0.0) }
      let mut j = i + 1
      while j < n_atoms {
        let rj = res_atom_idx[j]
        if ri == rj { j = j + 1; continue }
        let idx_j = self.indices[j]
        let cj = if idx_j < fr.length() { fr[idx_j] } else { Vector3::new(0.0, 0.0, 0.0) }
        let d2 = if pbc {
          mda_analysis_min_dist(ci, cj, box)
        } else {
          let dx = ci.x - cj.x; let dy = ci.y - cj.y; let dz = ci.z - cj.z
          dx * dx + dy * dy + dz * dz
        }
        if d2 <= cutoff2 {
          // 残基对 (min, max) 加一，保证唯一键
          let a = if ri < rj { ri } else { rj }
          let b = if ri < rj { rj } else { ri }
          count[a][b] = count[a][b] + 1
          f_count = f_count + 1
        }
        j = j + 1
      }
      i = i + 1
    }
    frame_contact_counts.push(f_count)
    f = f + 1
  }

  // ===== 3. 归一化 → 概率 =====
  let prob : Array[Array[Double]] = Array(capacity=n_res)
  i = 0
  while i < n_res {
    let row : Array[Double] = Array(capacity=n_res)
    let mut j = 0
    while j < n_res {
      if i == j { row.push(0.0) }
      else {
        let a = if i < j { i } else { j }
        let b = if i < j { j } else { i }
        row.push(count[a][b].to_double() / n_frames.to_double())
      }
      j = j + 1
    }
    prob.push(row)
    i = i + 1
  }

  // 收集残基起始索引
  let res_start_idx : Array[Int32] = Array(capacity=n_res)
  if res_per_atom {
    let mut k = 0
    while k < n_res { res_start_idx.push(k); k = k + 1 }
  } else {
    let mut k = 0
    let mut seen = -1
    while k < n_atoms {
      let r = res_atom_idx[k]
      if r != seen { res_start_idx.push(k); seen = r }
      k = k + 1
    }
  }

  ContactMapResult::{
    n_res,
    resids,
    resnames,
    res_start_idx,
    prob,
    frame_contact_counts,
    n_frames,
  }
}

// ============================================================================
// 查询方法
// ============================================================================

///|
/// 残基数。
pub fn ContactMapResult::n_residues(self : ContactMapResult) -> Int32 {
  self.n_res
}

///|
/// 接触概率矩阵（上三角 = 下三角镜像）。
pub fn ContactMapResult::probability_matrix(self : ContactMapResult) -> Array[Array[Double]] {
  self.prob
}

///|
/// 残基 ID 数组。
pub fn ContactMapResult::resid_ids(self : ContactMapResult) -> Array[Int32] {
  self.resids
}

///|
/// 残基名数组。
pub fn ContactMapResult::residue_names(self : ContactMapResult) -> Array[String] {
  self.resnames
}

///|
/// 最常出现的 Top-N 残基接触对（按概率降序）。
///
/// 返回 `Array[(Int32 i, Int32 j, Double prob)]`。
pub fn ContactMapResult::top_contacts(
  self : ContactMapResult,
  n_top~ : Int32 = 10,
) -> Array[(Int32, Int32, Double)] {
  let pairs : Array[(Int32, Int32, Double)] = Array(capacity=self.n_res * self.n_res / 2)
  let mut i = 0
  while i < self.n_res {
    let mut j = i + 1
    while j < self.n_res {
      let p = self.prob[i][j]
      if p > 0.0 { pairs.push((i, j, p)) }
      j = j + 1
    }
    i = i + 1
  }
  // 降序排序（简单插入排序，次数不大）
  let mut k = 1
  while k < pairs.length() {
    let key = pairs[k]
    let mut m = k - 1
    while m >= 0 {
      if pairs[m].2 < key.2 { pairs[m + 1] = pairs[m] } else { break }
      m = m - 1
    }
    pairs[m + 1] = key
    k = k + 1
  }
  let out : Array[(Int32, Int32, Double)] = Array(capacity=n_top)
  k = 0
  while k < pairs.length() && k < n_top {
    out.push(pairs[k])
    k = k + 1
  }
  out
}

///|
/// 每一帧的总接触数（按残基对统计，同残基内不计）。
pub fn ContactMapResult::frame_counts(self : ContactMapResult) -> Array[Int32] {
  self.frame_contact_counts
}

///|
/// 平均每帧接触数（残基对，同残基内不计）。
pub fn ContactMapResult::mean_frame_contacts(self : ContactMapResult) -> Double {
  let n = self.frame_contact_counts.length()
  if n == 0 { return 0.0 }
  let mut total = 0
  let mut i = 0
  while i < n { total = total + self.frame_contact_counts[i]; i = i + 1 }
  total.to_double() / n.to_double()
}

///|
/// 保存为 GROMACS-style xtc/ndx 之外的简单接触图 CSV（行列 resid）。
///
/// 格式：第一行 `,1,ALA,2,ALA,...`，后续行 `1,ALA,0,0.87,0.12,...`。
/// 返回 String，不自动写文件。
pub fn ContactMapResult::to_csv(self : ContactMapResult) -> String {
  let sb = StringBuilder()
  // Header
  sb.write_string(",")
  let mut j = 0
  while j < self.n_res {
    sb.write_string(self.resnames[j])
    sb.write_string(self.resids[j].to_string())
    if j + 1 < self.n_res { sb.write_string(",") }
    j = j + 1
  }
  sb.write_string("\n")
  let mut i = 0
  while i < self.n_res {
    sb.write_string(self.resnames[i])
    sb.write_string(self.resids[i].to_string())
    sb.write_string(",")
    j = 0
    while j < self.n_res {
      sb.write_string(self.prob[i][j].to_string())
      if j + 1 < self.n_res { sb.write_string(",") }
      j = j + 1
    }
    sb.write_string("\n")
    i = i + 1
  }
  sb.to_string()
}

// ============================================================================
// SVG / HTML 可视化
// ============================================================================

///|
/// 完整接触概率热力图（SVG only）。
///
/// 自包含 SVG，蓝→红渐变，含刻度、图例、tooltip。直接可浏览器打开。
/// 参数：
///   - `cell_size~` 每格像素（默认 6）
///   - `max_label~` 最多显示多少个残基名在轴上（默认 50，超过则抽样）
///   - `title~` SVG 标题（默认 "Contact Probability Heatmap"）
pub fn ContactMapResult::to_svg_heatmap(
  self : ContactMapResult,
  cell_size~ : Int32 = 6,
  max_label~ : Int32 = 50,
  title~ : String = "Contact Probability Heatmap",
) -> String {
  let n = self.n_res
  if n == 0 { return "<svg></svg>" }
  let cs = cell_size
  let pad_left = 80
  let pad_top = 80
  let pad_right = 180   // 留给 legend
  let pad_bottom = 80
  let width = pad_left + n * cs + pad_right
  let height = pad_top + n * cs + pad_bottom

  let sb = StringBuilder()
  sb.write_string("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"")
  sb.write_string(width.to_string())
  sb.write_string("\" height=\"")
  sb.write_string(height.to_string())
  sb.write_string("\" viewBox=\"0 0 ")
  sb.write_string(width.to_string())
  sb.write_string(" ")
  sb.write_string(height.to_string())
  sb.write_string("\">\n")
  // Styles
  sb.write_string("<style>\n")
  sb.write_string(".cm-cell:hover { stroke: #fff; stroke-width: 1.2; }\n")
  sb.write_string(".cm-label { font: 10px monospace; fill: #333; }\n")
  sb.write_string(".cm-title { font: bold 16px sans-serif; fill: #222; }\n")
  sb.write_string(".cm-legend-label { font: 10px monospace; fill: #333; }\n")
  sb.write_string(".cm-tooltip { font: 12px monospace; fill: #fff; }\n")
  sb.write_string("</style>\n")
  // Title
  sb.write_string("<text class=\"cm-title\" x=\"")
  sb.write_string((width / 2).to_string())
  sb.write_string("\" y=\"30\" text-anchor=\"middle\">")
  sb.write_string(title)
  sb.write_string("</text>\n")
  // Background
  sb.write_string("<rect x=\"")
  sb.write_string(pad_left.to_string())
  sb.write_string("\" y=\"")
  sb.write_string(pad_top.to_string())
  sb.write_string("\" width=\"")
  sb.write_string((n * cs).to_string())
  sb.write_string("\" height=\"")
  sb.write_string((n * cs).to_string())
  sb.write_string("\" fill=\"#fafafa\" stroke=\"#ccc\"/>\n")

  // Cell grid
  let max_label = if max_label > n { n } else { max_label }
  let step_x = if n <= max_label { 1 } else { (n / max_label) + 1 }
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      let p = self.prob[i][j]
      let color = mda_contact_prob_color(p)
      let x = pad_left + j * cs
      let y = pad_top + i * cs
      sb.write_string("<rect class=\"cm-cell\" x=\"")
      sb.write_string(x.to_string())
      sb.write_string("\" y=\"")
      sb.write_string(y.to_string())
      sb.write_string("\" width=\"")
      sb.write_string(cs.to_string())
      sb.write_string("\" height=\"")
      sb.write_string(cs.to_string())
      sb.write_string("\" fill=\"")
      sb.write_string(color)
      sb.write_string("\" data-i=\"")
      sb.write_string(i.to_string())
      sb.write_string("\" data-j=\"")
      sb.write_string(j.to_string())
      sb.write_string("\" data-p=\"")
      sb.write_string(p.to_string())
      sb.write_string("\"/>\n")
      j = j + 1
    }
    i = i + 1
  }

  // X labels
  sb.write_string("<g class=\"cm-xlabels\">\n")
  i = 0
  while i < n {
    if i % step_x == 0 {
      let x = pad_left + i * cs + cs / 2
      let y = pad_top - 8
      sb.write_string("<text class=\"cm-label\" x=\"")
      sb.write_string(x.to_string())
      sb.write_string("\" y=\"")
      sb.write_string(y.to_string())
      sb.write_string("\" text-anchor=\"end\" transform=\"rotate(-60 ")
      sb.write_string(x.to_string())
      sb.write_string(" ")
      sb.write_string(y.to_string())
      sb.write_string("\">")
      sb.write_string(self.resnames[i])
      sb.write_string(self.resids[i].to_string())
      sb.write_string("</text>\n")
    }
    i = i + 1
  }
  sb.write_string("</g>\n")
  // Y labels
  sb.write_string("<g class=\"cm-ylabels\">\n")
  i = 0
  while i < n {
    if i % step_x == 0 {
      let x = pad_left - 6
      let y = pad_top + i * cs + cs / 2
      sb.write_string("<text class=\"cm-label\" x=\"")
      sb.write_string(x.to_string())
      sb.write_string("\" y=\"")
      sb.write_string(y.to_string())
      sb.write_string("\" text-anchor=\"end\" dominant-baseline=\"middle\">")
      sb.write_string(self.resnames[i])
      sb.write_string(self.resids[i].to_string())
      sb.write_string("</text>\n")
    }
    i = i + 1
  }
  sb.write_string("</g>\n")

  // Color legend (右侧)
  let legend_x = pad_left + n * cs + 20
  let legend_y = pad_top
  let legend_h = n * cs
  let legend_w = 16
  sb.write_string("<g class=\"cm-legend\">\n")
  sb.write_string("<text class=\"cm-title\" x=\"")
  sb.write_string((legend_x + legend_w / 2).to_string())
  sb.write_string("\" y=\"")
  sb.write_string((legend_y - 16).to_string())
  sb.write_string("\" text-anchor=\"middle\" style=\"font-size:11px\">Probability</text>\n")
  let steps = 40
  let mut s = 0
  while s < steps {
    let t = s.to_double() / (steps.to_double() - 1.0)
    let val = t
    let color = mda_contact_prob_color(val)
    let rect_y = legend_y + (steps - 1 - s) * (legend_h.to_double() / steps.to_double()).to_int()
    let rect_h = (legend_h.to_double() / steps.to_double()).to_int() + 1
    sb.write_string("<rect x=\"")
    sb.write_string(legend_x.to_string())
    sb.write_string("\" y=\"")
    sb.write_string(rect_y.to_string())
    sb.write_string("\" width=\"")
    sb.write_string(legend_w.to_string())
    sb.write_string("\" height=\"")
    sb.write_string(rect_h.to_string())
    sb.write_string("\" fill=\"")
    sb.write_string(color)
    sb.write_string("\"/>\n")
    s = s + 1
  }
  sb.write_string("<rect x=\"")
  sb.write_string(legend_x.to_string())
  sb.write_string("\" y=\"")
  sb.write_string(legend_y.to_string())
  sb.write_string("\" width=\"")
  sb.write_string(legend_w.to_string())
  sb.write_string("\" height=\"")
  sb.write_string(legend_h.to_string())
  sb.write_string("\" fill=\"none\" stroke=\"#555\"/>\n")
  // Legend labels 1.0, 0.5, 0.0
  sb.write_string("<text class=\"cm-legend-label\" x=\"")
  sb.write_string((legend_x + legend_w + 4).to_string())
  sb.write_string("\" y=\"")
  sb.write_string((legend_y + 4).to_string())
  sb.write_string("\">1.0</text>\n")
  sb.write_string("<text class=\"cm-legend-label\" x=\"")
  sb.write_string((legend_x + legend_w + 4).to_string())
  sb.write_string("\" y=\"")
  sb.write_string((legend_y + legend_h / 2 + 4).to_string())
  sb.write_string("\">0.5</text>\n")
  sb.write_string("<text class=\"cm-legend-label\" x=\"")
  sb.write_string((legend_x + legend_w + 4).to_string())
  sb.write_string("\" y=\"")
  sb.write_string((legend_y + legend_h + 4).to_string())
  sb.write_string("\">0.0</text>\n")
  sb.write_string("</g>\n")
  sb.write_string("</svg>\n")
  sb.to_string()
}

///|
/// 完整 HTML 页面（自包含，内嵌 SVG + hover tooltip）。
///
/// 零外部依赖：直接保存为 .html 用浏览器打开。
pub fn ContactMapResult::to_html_heatmap(
  self : ContactMapResult,
  cell_size~ : Int32 = 6,
  title~ : String = "Contact Probability Heatmap",
) -> String {
  let n = self.n_res
  if n == 0 { return "<html><body>No data</body></html>" }
  let svg = self.to_svg_heatmap(cell_size=cell_size, title=title, max_label=60)
  let sb = StringBuilder()
  sb.write_string("<!DOCTYPE html>\n<html>\n<head>\n")
  sb.write_string("<meta charset=\"utf-8\">\n")
  sb.write_string("<title>")
  sb.write_string(title)
  sb.write_string("</title>\n")
  sb.write_string("<style>\n")
  sb.write_string("body { font-family: -apple-system, Segoe UI, sans-serif; margin: 20px; background: #fff; }\n")
  sb.write_string(".wrap { display: inline-block; position: relative; }\n")
  sb.write_string(".tooltip { position: absolute; background: rgba(0,0,0,0.85); color: #fff;\n")
  sb.write_string("  padding: 6px 10px; border-radius: 4px; font: 12px monospace;\n")
  sb.write_string("  pointer-events: none; opacity: 0; transition: opacity 0.15s;\n")
  sb.write_string("  z-index: 999; }\n")
  sb.write_string(".wrap:hover .tooltip { opacity: 1; }\n")
  sb.write_string(".info { margin: 20px 0; font-size: 13px; color: #555; }\n")
  sb.write_string(".info b { color: #222; }\n")
  sb.write_string("</style>\n</head>\n<body>\n")
  sb.write_string("<div class=\"info\">\n")
  sb.write_string("<b>")
  sb.write_string(self.n_frames.to_string())
  sb.write_string("</b> frames &nbsp;×&nbsp; <b>")
  sb.write_string(n.to_string())
  sb.write_string("</b> residues &nbsp;|&nbsp; <b>")
  sb.write_string(self.mean_frame_contacts().to_string())
  sb.write_string("</b> mean contacts / frame\n")
  sb.write_string("</div>\n")
  sb.write_string("<div class=\"wrap\">\n")
  sb.write_string(svg)
  sb.write_string("<div class=\"tooltip\" id=\"cm-tip\"></div>\n")
  sb.write_string("</div>\n")
  sb.write_string("<script>\n")
  sb.write_string("(function() {\n")
  sb.write_string("  var wrap = document.querySelector('.wrap');\n")
  sb.write_string("  var tip = document.getElementById('cm-tip');\n")
  sb.write_string("  var cells = document.querySelectorAll('.cm-cell');\n")
  sb.write_string("  cells.forEach(function(c) {\n")
  sb.write_string("    c.addEventListener('mouseenter', function(ev) {\n")
  sb.write_string("      var i = c.getAttribute('data-i');\n")
  sb.write_string("      var j = c.getAttribute('data-j');\n")
  sb.write_string("      var p = parseFloat(c.getAttribute('data-p'));\n")
  sb.write_string("      tip.innerHTML = '[' + i + ',' + j + ']  p=' + p.toFixed(3);\n")
  sb.write_string("      var r = wrap.getBoundingClientRect();\n")
  sb.write_string("      tip.style.left = (ev.clientX - r.left + 12) + 'px';\n")
  sb.write_string("      tip.style.top  = (ev.clientY - r.top  + 12) + 'px';\n")
  sb.write_string("    });\n")
  sb.write_string("    c.addEventListener('mousemove', function(ev) {\n")
  sb.write_string("      var r = wrap.getBoundingClientRect();\n")
  sb.write_string("      tip.style.left = (ev.clientX - r.left + 12) + 'px';\n")
  sb.write_string("      tip.style.top  = (ev.clientY - r.top  + 12) + 'px';\n")
  sb.write_string("    });\n")
  sb.write_string("    c.addEventListener('mouseleave', function() {\n")
  sb.write_string("      tip.style.opacity = 0;\n")
  sb.write_string("    });\n")
  sb.write_string("  });\n")
  sb.write_string("  wrap.addEventListener('mouseenter', function() { tip.style.opacity = 1; });\n")
  sb.write_string("})();\n")
  sb.write_string("</script>\n")
  sb.write_string("</body>\n</html>\n")
  sb.to_string()
}

// ============================================================================
// 内部配色：蓝→青→黄→红渐变（Viridis-style 四段插值）
// ============================================================================

///|
/// 把 [0,1] 概率映射成 #rrggbb 颜色字符串。
/// 0.0 = 深蓝, 0.5 = 青绿, 1.0 = 亮红。
fn mda_contact_prob_color(p : Double) -> String {
  let t = if p < 0.0 { 0.0 } else if p > 1.0 { 1.0 } else { p }
  // 四段 Viridis-like 插值
  let (r, g, b) : (Int32, Int32, Int32) = if t < 0.333 {
    let u = t / 0.333
    // 蓝(40,3,127) → 青(30,144,135)
    (
      (40.0 + (30.0 - 40.0) * u).to_int(),
      (3.0 + (144.0 - 3.0) * u).to_int(),
      (127.0 + (135.0 - 127.0) * u).to_int(),
    )
  } else if t < 0.666 {
    let u = (t - 0.333) / 0.333
    // 青(30,144,135) → 黄(253,231,37)
    (
      (30.0 + (253.0 - 30.0) * u).to_int(),
      (144.0 + (231.0 - 144.0) * u).to_int(),
      (135.0 + (37.0 - 135.0) * u).to_int(),
    )
  } else {
    let u = (t - 0.666) / 0.334
    // 黄(253,231,37) → 红(244,67,54)
    (
      (253.0 + (244.0 - 253.0) * u).to_int(),
      (231.0 + (67.0 - 231.0) * u).to_int(),
      (37.0 + (54.0 - 37.0) * u).to_int(),
    )
  }
  let sb = StringBuilder()
  sb.write_string("#")
  sb.write_string(mda_hex_byte(r))
  sb.write_string(mda_hex_byte(g))
  sb.write_string(mda_hex_byte(b))
  sb.to_string()
}

// ============================================================================
// A4: Native Contact Fraction 时间序列
// ============================================================================

///|
/// Native contact fraction 时间序列。
///
/// 定义：先用 frame 0 的 Cα-Cα 距离 ≤ 8 Å 定义一组 native contacts；
/// 然后每一帧统计这些 native 接触中还有多少保持 ≤ cutoff（默认 8 Å）。
/// 返回 Array[Double]，每帧一个 fraction。
///
/// 这是判断蛋白质折叠稳定性/路径的经典指标。
pub fn AtomGroup::native_contact_fraction(
  self : AtomGroup,
  cutoff~ : Double = 8.0,
  pbc~ : Bool = false,
) -> Array[Double] {
  let atoms = self.universe.atoms
  let n_frames = self.universe.n_frames()
  let box = if pbc { Some(self.universe.box) } else { None }
  let cutoff2 = cutoff * cutoff

  // 1. 识别残基边界 + 找每个残基的 Cα（或首原子兜底）
  let (starts, _resids) = mda_analysis_residue_boundaries(atoms, self.indices)
  let n_res = starts.length()
  let ca_idx : Array[Int32] = Array(capacity=n_res)
  let mut ri = 0
  while ri < n_res {
    let cur_s = starts[ri]
    let cur_e = if ri + 1 < n_res { starts[ri + 1] } else { self.indices.length() }
    // 优先找 CA，找不到就用第一个原子
    let mut found = -1
    let mut k = cur_s
    while k < cur_e && k < self.indices.length() {
      let idx = self.indices[k]
      if idx < atoms.length() && atoms[idx].name.to_upper() == "CA" {
        found = idx; break
      }
      k = k + 1
    }
    if found < 0 && cur_s < self.indices.length() { found = self.indices[cur_s] }
    ca_idx.push(found)
    ri = ri + 1
  }

  // 2. 在 frame 0 上找 native contacts（Cα-Cα ≤ cutoff）
  let fr0 = self.universe.frame(0)
  // native_contacts 用上三角 (i,j) 唯一对
  let native_i : Array[Int32] = Array(capacity=n_res * n_res / 2)
  let native_j : Array[Int32] = Array(capacity=n_res * n_res / 2)
  ri = 0
  while ri < n_res {
    if ca_idx[ri] < 0 || ca_idx[ri] >= fr0.length() { ri = ri + 1; continue }
    let ci = fr0[ca_idx[ri]]
    let mut rj = ri + 1
    while rj < n_res {
      if ca_idx[rj] < 0 || ca_idx[rj] >= fr0.length() { rj = rj + 1; continue }
      let cj = fr0[ca_idx[rj]]
      let d2 = if pbc {
        mda_analysis_min_dist(ci, cj, box)
      } else {
        let dx = ci.x - cj.x; let dy = ci.y - cj.y; let dz = ci.z - cj.z
        dx * dx + dy * dy + dz * dz
      }
      if d2 <= cutoff2 {
        native_i.push(ri)
        native_j.push(rj)
      }
      rj = rj + 1
    }
    ri = ri + 1
  }

  let n_native = native_i.length()
  let out : Array[Double] = Array(capacity=n_frames)
  if n_native == 0 {
    let mut f = 0
    while f < n_frames { out.push(0.0); f = f + 1 }
    return out
  }

  // 3. 逐帧统计保留比例
  let mut f = 0
  while f < n_frames {
    let fr = self.universe.frame(f)
    let mut kept = 0
    let mut k = 0
    while k < n_native {
      let a = native_i[k]
      let b = native_j[k]
      let ia = ca_idx[a]
      let ib = ca_idx[b]
      if ia >= 0 && ib >= 0 && ia < fr.length() && ib < fr.length() {
        let ca_pos = fr[ia]
        let cb_pos = fr[ib]
        let d2 = if pbc {
          mda_analysis_min_dist(ca_pos, cb_pos, box)
        } else {
          let dx = ca_pos.x - cb_pos.x; let dy = ca_pos.y - cb_pos.y; let dz = ca_pos.z - cb_pos.z
          dx * dx + dy * dy + dz * dz
        }
        if d2 <= cutoff2 { kept = kept + 1 }
      }
      k = k + 1
    }
    out.push(kept.to_double() / n_native.to_double())
    f = f + 1
  }
  out
}

fn mda_hex_byte(v : Int32) -> String {
  let x = if v < 0 { 0 } else if v > 255 { 255 } else { v }
  let hi = x / 16
  let lo = x % 16
  let sb = StringBuilder()
  sb.write_char(mda_hex_digit(hi))
  sb.write_char(mda_hex_digit(lo))
  sb.to_string()
}

fn mda_hex_digit(n : Int32) -> Char {
  if n < 10 { (n + 48).unsafe_to_char() } else { (n - 10 + 97).unsafe_to_char() }
}
