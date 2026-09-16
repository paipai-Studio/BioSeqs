///|
/// Bio.Graphics.ColorSpiral：HSV 颜色螺旋。
///
/// 在 HSV 颜色空间中沿螺旋线（色相旋转 + 明度递减）生成一组等高
/// 对比度、视觉上均匀分布的颜色，供图表系列、通路节点与基因组
/// 可视化配色。对应 Biopython Bio.Graphics.ColorSpiral。

/// 一个 RGB 颜色（0-255）。
pub(all) struct SpiralColor {
  r : Int
  g : Int
  b : Int
}

/// 颜色螺旋状态与参数。
pub(all) struct ColorSpiral {
  h0 : Double         // 起始色相 [0,1)
  hue_step : Double   // 每取一色色相旋转步长（圈的分数）
  s0 : Double         // 饱和度 [0,1]
  v0 : Double         // 起始明度 [0,1]
  v_decay : Double    // 每取一色明度衰减量
  vmin : Double       // 明度下限
  vmax : Double       // 明度上限
}

/// 构造一个颜色值。
pub fn SpiralColor::new(r : Int, g : Int, b : Int) -> SpiralColor {
  SpiralColor::{ r: r, g: g, b: b }
}

fn clamp255(v : Int) -> Int {
  if v < 0 {
    return 0
  }
  if v > 255 {
    return 255
  }
  v
}

/// 默认绿色起点的颜色螺旋（每 60 色完成 1 圈，明度缓慢下降）。
pub fn spiral_new() -> ColorSpiral {
  ColorSpiral::{
    h0: 0.3333333333333333, // 绿
    hue_step: 1.0 / 60.0,
    s0: 0.85,
    v0: 0.85,
    v_decay: 0.05,
    vmin: 0.30,
    vmax: 0.90,
  }
}

/// 自定义参数的颜色螺旋。
pub fn spiral_new_with(
  h0 : Double,
  hue_step : Double,
  s0 : Double,
  v0 : Double,
  v_decay : Double,
  vmin? : Double = 0.30,
  vmax? : Double = 0.90,
) -> ColorSpiral {
  ColorSpiral::{
    h0: h0, hue_step: hue_step, s0: s0, v0: v0, v_decay: v_decay,
    vmin: vmin, vmax: vmax,
  }
}

/// HSV 转 RGB（h/s/v 均 [0,1]）。
pub fn spiral_hsv_to_rgb(h : Double, s : Double, v : Double) -> SpiralColor {
  let mut hh = h - h.floor()
  if hh < 0.0 {
    hh = hh + 1.0
  }
  let i = (hh * 6.0).floor().to_int()
  let f = hh * 6.0 - i.to_double()
  let p = v * (1.0 - s)
  let q = v * (1.0 - f * s)
  let t = v * (1.0 - (1.0 - f) * s)
  let mut r = 0.0
  let mut g = 0.0
  let mut b = 0.0
  match i % 6 {
    0 => {
      r = v
      g = t
      b = p
    }
    1 => {
      r = q
      g = v
      b = p
    }
    2 => {
      r = p
      g = v
      b = t
    }
    3 => {
      r = p
      g = q
      b = v
    }
    4 => {
      r = t
      g = p
      b = v
    }
    _ => {
      r = v
      g = p
      b = q
    }
  }
  SpiralColor::{
    r: clamp255((r * 255.0).round().to_int()),
    g: clamp255((g * 255.0).round().to_int()),
    b: clamp255((b * 255.0).round().to_int()),
  }
}

/// 取螺旋中第 idx 个颜色（从 0 开始）。
pub fn spiral_color_at(sp : ColorSpiral, idx : Int) -> SpiralColor {
  let i = if idx < 0 { 0 } else { idx }
  let h = sp.h0 + i.to_double() * sp.hue_step
  let mut v = sp.v0 - i.to_double() * sp.v_decay
  if v < sp.vmin {
    v = sp.vmin
  }
  if v > sp.vmax {
    v = sp.vmax
  }
  spiral_hsv_to_rgb(h, sp.s0, v)
}

/// 取一串（n 个）颜色。
pub fn spiral_colors(sp : ColorSpiral, n : Int) -> Array[SpiralColor] {
  let out : Array[SpiralColor] = Array::new()
  for i in 0..<n {
    out.push(spiral_color_at(sp, i))
  }
  out
}

/// 当前取色位置的 HSV 表示（测试/调试用）。
pub fn spiral_hsv_at(sp : ColorSpiral, idx : Int) -> (Double, Double, Double) {
  let i = if idx < 0 { 0 } else { idx }
  let h = sp.h0 + i.to_double() * sp.hue_step
  let mut v = sp.v0 - i.to_double() * sp.v_decay
  if v < sp.vmin {
    v = sp.vmin
  }
  if v > sp.vmax {
    v = sp.vmax
  }
  (h, sp.s0, v)
}

fn hex_digit(v : Int) -> Char {
  if v < 10 {
    return Char::from_int('0'.to_int() + v)
  }
  Char::from_int('A'.to_int() + v - 10)
}

/// RGB 转十六进制字符串（如 "#FF0000"）。
pub fn spiral_hex(c : SpiralColor) -> String {
  let sb = StringBuilder()
  sb.write_string("#")
  sb.write_char(hex_digit(c.r / 16))
  sb.write_char(hex_digit(c.r % 16))
  sb.write_char(hex_digit(c.g / 16))
  sb.write_char(hex_digit(c.g % 16))
  sb.write_char(hex_digit(c.b / 16))
  sb.write_char(hex_digit(c.b % 16))
  sb.to_string()
}

/// RGB 转 CSS rgb() 字符串。
pub fn spiral_rgb(c : SpiralColor) -> String {
  "rgb(" + c.r.to_string() + "," + c.g.to_string() + "," + c.b.to_string() + ")"
}

/// 亮度（YIQ 感知亮度，0-255）。
pub fn spiral_luminance(c : SpiralColor) -> Double {
  (0.299 * c.r.to_double() + 0.587 * c.g.to_double() + 0.114 * c.b.to_double())
}

/// 对比度颜色（高亮底色上返回黑、深色底色上返回白），用于文字。
pub fn spiral_contrast(c : SpiralColor) -> SpiralColor {
  if spiral_luminance(c) > 140.0 {
    SpiralColor::{ r: 0, g: 0, b: 0 }
  } else {
    SpiralColor::{ r: 255, g: 255, b: 255 }
  }
}

/// SVG 色板（一行色块，可选标签）。
pub fn spiral_svg_palette(
  colors : Array[SpiralColor],
  labels? : Array[String],
) -> String {
  let sb = StringBuilder()
  sb.write_string(
    "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" +
    (colors.length() * 60).to_string() +
    "\" height=\"50\" viewBox=\"0 0 " +
    (colors.length() * 60).to_string() +
    " 50\">\n",
  )
  for i in 0..<colors.length() {
    let x = i * 60
    sb.write_string(
      "<rect x=\"" + x.to_string() + "\" y=\"5\" width=\"55\" height=\"30\" fill=\"" +
      spiral_hex(colors[i]) + "\" stroke=\"#000\" stroke-width=\"1\"/>\n",
    )
    if labels.is_some() && labels.unwrap().length() > i {
      sb.write_string(
        "<text x=\"" + (x + 27).to_string() + "\" y=\"47\" font-size=\"9\" text-anchor=\"middle\">" +
        labels.unwrap()[i] + "</text>\n",
      )
    }
  }
  sb.write_string("</svg>")
  sb.to_string()
}

/// 将颜色螺旋用于数值着色：值域 [lo, hi] → 色相渐变（蓝→红）。
pub fn spiral_scale_color(
  value : Double,
  lo : Double,
  hi : Double,
  h_lo? : Double = 0.666,  // 蓝
  h_hi? : Double = 0.0,    // 红
) -> SpiralColor {
  let mut t = if hi <= lo {
    0.5
  } else {
    (value - lo) / (hi - lo)
  }
  if t < 0.0 {
    t = 0.0
  }
  if t > 1.0 {
    t = 1.0
  }
  let h = h_lo + (h_hi - h_lo) * t
  spiral_hsv_to_rgb(h, 0.9, 0.9)
}