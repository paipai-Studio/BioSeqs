///|
/// Distribution plots rendered to SVG, ported from Biopython's
/// `Bio.Graphics.DistributionPlot` module.
///
/// The reference module exposes `Histogram`, `Scatter`, `LinePlot` and a
/// `DistributionPage` that lays several panels out on a grid and writes the
/// result to a file. This port keeps the same vocabulary while emitting an
/// SVG string (no external rendering backend is required), making it easy
/// to embed the figures in reports or test snapshots.

///|
/// A 2D point used by scatter and line plots.
pub struct PlotPoint {
  x : Double
  y : Double
} derive(Eq)

///|
pub fn PlotPoint::new(x~ : Double, y~ : Double) -> PlotPoint {
  PlotPoint::{ x, y, }
}

///|
/// A labeled series with an SVG stroke/fill color.
pub struct PlotSeries {
  points : Array[PlotPoint]
  color : String
  label : String
} derive(Eq)

///|
pub fn PlotSeries::new(
  points~ : Array[PlotPoint],
  color~ : String = "#1f77b4",
  label~ : String = "",
) -> PlotSeries {
  PlotSeries::{ points, color, label, }
}

///|
/// One histogram bin: half-open `[lo, hi)` with `count` observations.
pub struct HistBin {
  lo : Double
  hi : Double
  count : Int
} derive(Eq)

///|
/// Histogram data: bin edges plus counts, ready to render.
pub struct Histogram {
  bins : Array[HistBin]
  data_min : Double
  data_max : Double
  total : Int
} derive(Eq)

///|
/// Build a histogram from raw values.
///
/// `num_bins` defaults to 20 (the Biopython default). When every value is
/// equal the histogram collapses to a single bin spanning `[v, v+1)`.
pub fn Histogram::from_data(
  data : Array[Double],
  num_bins~ : Int = 20,
) -> Histogram {
  if data.is_empty() {
    return Histogram::{ bins: [], data_min: 0.0, data_max: 0.0, total: 0, }
  }
  let mut lo = data[0]
  let mut hi = data[0]
  for v in data {
    if v < lo {
      lo = v
    }
    if v > hi {
      hi = v
    }
  }
  let span = hi - lo
  let n = if num_bins < 1 { 1 } else { num_bins }
  if span == 0.0 {
    let bins = [HistBin::{ lo, hi: lo + 1.0, count: data.length(), }]
    return Histogram::{ bins, data_min: lo, data_max: hi, total: data.length(), }
  }
  let width = span / n.to_double()
  let counts : Array[Int] = Array::new()
  for _i in 0..<n {
    counts.push(0)
  }
  for v in data {
    let mut idx = ((v - lo) / width).to_int()
    if idx < 0 {
      idx = 0
    }
    if idx >= n {
      idx = n - 1
    }
    counts[idx] = counts[idx] + 1
  }
  let bins : Array[HistBin] = Array::new()
  for i in 0..<n {
    let bin_lo = lo + i.to_double() * width
    let bin_hi = if i == n - 1 { hi } else { bin_lo + width }
    bins.push(HistBin::{ lo: bin_lo, hi: bin_hi, count: counts[i], })
  }
  Histogram::{ bins, data_min: lo, data_max: hi, total: data.length(), }
}

///|
/// Maximum bin count (used to scale bar heights).
pub fn Histogram::max_count(self : Histogram) -> Int {
  let mut m = 0
  for bin in self.bins {
    if bin.count > m {
      m = bin.count
    }
  }
  m
}

///|
/// A scatter plot of one or more series.
pub struct ScatterPlot {
  series : Array[PlotSeries]
} derive(Eq)

///|
pub fn ScatterPlot::new(series~ : Array[PlotSeries]) -> ScatterPlot {
  ScatterPlot::{ series, }
}

///|
/// A line plot of one or more series (each rendered as a polyline).
pub struct LinePlot {
  series : Array[PlotSeries]
} derive(Eq)

///|
pub fn LinePlot::new(series~ : Array[PlotSeries]) -> LinePlot {
  LinePlot::{ series, }
}

///|
/// One panel of a distribution page.
pub(all) enum PlotPanel {
  HistPanel(String, Histogram)
  ScatterPanel(String, ScatterPlot)
  LinePanel(String, LinePlot)
} derive(Eq)

///|
/// A grid of panels rendered together, mirroring `DistributionPage`.
pub struct DistributionPage {
  title : String
  rows : Int
  cols : Int
  panels : Array[PlotPanel]
} derive(Eq)

///|
pub fn DistributionPage::new(
  title~ : String,
  rows~ : Int = 1,
  cols~ : Int = 1,
) -> DistributionPage {
  DistributionPage::{ title, rows, cols, panels: [], }
}

///|
/// Append a panel; panels fill the grid in row-major order.
pub fn DistributionPage::add_panel(
  self : DistributionPage,
  panel : PlotPanel,
) -> DistributionPage {
  let panels = Array::new()
  for existing in self.panels {
    panels.push(existing)
  }
  panels.push(panel)
  DistributionPage::{ ..self, panels, }
}

///|
fn plot_format_double(v : Double) -> String {
  let rounded = (v * 1000.0).to_int().to_double() / 1000.0
  rounded.to_string()
}

///|
/// "Nice" tick step: round `value` up to 1, 2 or 5 x 10^k.
fn plot_nice_num(value : Double, round~ : Bool = false) -> Double {
  if value == 0.0 {
    return 0.0
  }
  let mut v = if value < 0.0 { 0.0 - value } else { value }
  let mut exponent = 0
  while v >= 10.0 {
    v = v / 10.0
    exponent += 1
  }
  while v < 1.0 {
    v = v * 10.0
    exponent -= 1
  }
  let nice_fraction = if round {
    if v < 1.5 {
      1.0
    } else if v < 3.0 {
      2.0
    } else if v < 7.0 {
      5.0
    } else {
      10.0
    }
  } else {
    if v <= 1.0 {
      1.0
    } else if v <= 2.0 {
      2.0
    } else if v <= 5.0 {
      5.0
    } else {
      10.0
    }
  }
  let mut result = nice_fraction
  if exponent >= 0 {
    for _i in 0..<exponent {
      result = result * 10.0
    }
  } else {
    for _i in 0..<(0 - exponent) {
      result = result / 10.0
    }
  }
  result
}

///|
/// Generate roughly `target` tick values spanning `[lo, hi]`.
fn plot_ticks(lo : Double, hi : Double, target~ : Int = 5) -> Array[Double] {
  if lo == hi {
    return [lo]
  }
  let span = hi - lo
  let step = plot_nice_num(span / target.to_double(), round=true)
  let mut start = (lo / step).floor() * step
  let ticks : Array[Double] = Array::new()
  let mut iter_count = 0
  while start <= hi + step * 0.5 && iter_count < 64 {
    if start >= lo - step * 0.5 {
      ticks.push(start)
    }
    start += step
    iter_count += 1
  }
  ticks
}

///|
fn plot_append_axis(
  svg : StringBuilder,
  x0 : Double,
  y0 : Double,
  x1 : Double,
  y1 : Double,
) -> Unit {
  svg.write_string(
    "<line x1=\"" +
    plot_format_double(x0) +
    "\" y1=\"" +
    plot_format_double(y0) +
    "\" x2=\"" +
    plot_format_double(x1) +
    "\" y2=\"" +
    plot_format_double(y1) +
    "\" stroke=\"#333\" stroke-width=\"1\"/>\n",
  )
}

///|
fn plot_append_text(
  svg : StringBuilder,
  x : Double,
  y : Double,
  text : String,
  anchor~ : String = "middle",
  size~ : Int = 10,
) -> Unit {
  svg.write_string(
    "<text x=\"" +
    plot_format_double(x) +
    "\" y=\"" +
    plot_format_double(y) +
    "\" font-size=\"" +
    size.to_string() +
    "\" text-anchor=\"" +
    anchor +
    "\">" +
    text +
    "</text>\n",
  )
}

///|
/// Render a histogram panel to SVG inside the given rectangle.
pub fn Histogram::svg_rect(
  self : Histogram,
  x : Double,
  y : Double,
  w : Double,
  h : Double,
) -> String {
  let svg = StringBuilder::new()
  let margin_l = 8.0
  let margin_r = 4.0
  let margin_t = 4.0
  let margin_b = 14.0
  let plot_x = x + margin_l
  let plot_y = y + margin_t
  let plot_w = w - margin_l - margin_r
  let plot_h = h - margin_t - margin_b
  let max_count = self.max_count()
  let x_scale = if self.data_max == self.data_min {
    1.0
  } else {
    plot_w / (self.data_max - self.data_min)
  }
  let y_scale = if max_count == 0 { 0.0 } else { plot_h / max_count.to_double() }
  for bin in self.bins {
    let bx = plot_x + (bin.lo - self.data_min) * x_scale
    let bw = (bin.hi - bin.lo) * x_scale
    let bh = bin.count.to_double() * y_scale
    let by = plot_y + plot_h - bh
    svg.write_string(
      "<rect x=\"" +
      plot_format_double(bx) +
      "\" y=\"" +
      plot_format_double(by) +
      "\" width=\"" +
      plot_format_double(bw) +
      "\" height=\"" +
      plot_format_double(bh) +
      "\" fill=\"#1f77b4\" stroke=\"none\"/>\n",
    )
  }
  plot_append_axis(svg, plot_x, plot_y + plot_h, plot_x + plot_w, plot_y + plot_h)
  plot_append_axis(svg, plot_x, plot_y, plot_x, plot_y + plot_h)
  for tick in plot_ticks(self.data_min, self.data_max, target=4) {
    let tx = plot_x + (tick - self.data_min) * x_scale
    plot_append_text(svg, tx, plot_y + plot_h + 11.0, plot_format_double(tick))
  }
  svg.to_string()
}

///|
/// Render a scatter panel inside the given rectangle.
pub fn ScatterPlot::svg_rect(
  self : ScatterPlot,
  x : Double,
  y : Double,
  w : Double,
  h : Double,
) -> String {
  let svg = StringBuilder::new()
  let margin_l = 8.0
  let margin_r = 4.0
  let margin_t = 4.0
  let margin_b = 14.0
  let plot_x = x + margin_l
  let plot_y = y + margin_t
  let plot_w = w - margin_l - margin_r
  let plot_h = h - margin_t - margin_b
  let (x_lo, x_hi, y_lo, y_hi) = plot_bounds(self.series)
  let x_scale = if x_hi == x_lo { 1.0 } else { plot_w / (x_hi - x_lo) }
  let y_scale = if y_hi == y_lo { 1.0 } else { plot_h / (y_hi - y_lo) }
  for series in self.series {
    for point in series.points {
      let px = plot_x + (point.x - x_lo) * x_scale
      let py = plot_y + plot_h - (point.y - y_lo) * y_scale
      svg.write_string(
        "<circle cx=\"" +
        plot_format_double(px) +
        "\" cy=\"" +
        plot_format_double(py) +
        "\" r=\"2\" fill=\"" +
        series.color +
        "\"/>\n",
      )
    }
  }
  plot_append_axis(svg, plot_x, plot_y + plot_h, plot_x + plot_w, plot_y + plot_h)
  plot_append_axis(svg, plot_x, plot_y, plot_x, plot_y + plot_h)
  svg.to_string()
}

///|
/// Render a line panel inside the given rectangle.
pub fn LinePlot::svg_rect(
  self : LinePlot,
  x : Double,
  y : Double,
  w : Double,
  h : Double,
) -> String {
  let svg = StringBuilder::new()
  let margin_l = 8.0
  let margin_r = 4.0
  let margin_t = 4.0
  let margin_b = 14.0
  let plot_x = x + margin_l
  let plot_y = y + margin_t
  let plot_w = w - margin_l - margin_r
  let plot_h = h - margin_t - margin_b
  let (x_lo, x_hi, y_lo, y_hi) = plot_bounds(self.series)
  let x_scale = if x_hi == x_lo { 1.0 } else { plot_w / (x_hi - x_lo) }
  let y_scale = if y_hi == y_lo { 1.0 } else { plot_h / (y_hi - y_lo) }
  for series in self.series {
    if series.points.length() < 2 {
      continue
    }
    svg.write_string("<polyline points=\"")
    for i, point in series.points {
      let px = plot_x + (point.x - x_lo) * x_scale
      let py = plot_y + plot_h - (point.y - y_lo) * y_scale
      if i > 0 {
        svg.write_string(" ")
      }
      svg.write_string(plot_format_double(px))
      svg.write_string(",")
      svg.write_string(plot_format_double(py))
    }
    svg.write_string(
      "\" fill=\"none\" stroke=\"" + series.color + "\" stroke-width=\"1.5\"/>\n",
    )
  }
  plot_append_axis(svg, plot_x, plot_y + plot_h, plot_x + plot_w, plot_y + plot_h)
  plot_append_axis(svg, plot_x, plot_y, plot_x, plot_y + plot_h)
  svg.to_string()
}

///|
fn plot_bounds(series : Array[PlotSeries]) -> (Double, Double, Double, Double) {
  let mut x_lo = 0.0
  let mut x_hi = 0.0
  let mut y_lo = 0.0
  let mut y_hi = 0.0
  let mut seen = false
  for s in series {
    for p in s.points {
      if !seen {
        x_lo = p.x
        x_hi = p.x
        y_lo = p.y
        y_hi = p.y
        seen = true
      } else {
        if p.x < x_lo {
          x_lo = p.x
        }
        if p.x > x_hi {
          x_hi = p.x
        }
        if p.y < y_lo {
          y_lo = p.y
        }
        if p.y > y_hi {
          y_hi = p.y
        }
      }
    }
  }
  if !seen {
    return (0.0, 1.0, 0.0, 1.0)
  }
  if x_lo == x_hi {
    x_hi = x_lo + 1.0
  }
  if y_lo == y_hi {
    y_hi = y_lo + 1.0
  }
  (x_lo, x_hi, y_lo, y_hi)
}

///|
/// Render the full page to an SVG document of the given pixel size.
pub fn DistributionPage::to_svg(
  self : DistributionPage,
  width~ : Int = 600,
  height~ : Int = 400,
) -> String {
  let svg = StringBuilder::new()
  svg.write_string(
    "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" +
    width.to_string() +
    "\" height=\"" +
    height.to_string() +
    "\" viewBox=\"0 0 " +
    width.to_string() +
    " " +
    height.to_string() +
    "\">\n",
  )
  svg.write_string("<rect x=\"0\" y=\"0\" width=\"" + width.to_string() + "\" height=\"" + height.to_string() + "\" fill=\"white\"/>\n")
  let title_h = if self.title.length() > 0 { 18.0 } else { 0.0 }
  if title_h > 0.0 {
    plot_append_text(
      svg,
      width.to_double() / 2.0,
      14.0,
      self.title,
      anchor="middle",
      size=14,
    )
  }
  let grid_w = width.to_double()
  let grid_h = height.to_double() - title_h
  let cell_w = grid_w / self.cols.to_double()
  let cell_h = grid_h / self.rows.to_double()
  for index, panel in self.panels {
    let row = index / self.cols
    let col = index % self.cols
    if row >= self.rows {
      break
    }
    let cx = col.to_double() * cell_w
    let cy = title_h + row.to_double() * cell_h
    let body = match panel {
      HistPanel(title, hist) =>
        if title.length() > 0 {
          plot_append_text(
            svg,
            cx + cell_w / 2.0,
            cy + 12.0,
            title,
            size=11,
          )
          hist.svg_rect(cx, cy + 16.0, cell_w, cell_h - 16.0)
        } else {
          hist.svg_rect(cx, cy, cell_w, cell_h)
        }
      ScatterPanel(title, scatter) =>
        if title.length() > 0 {
          plot_append_text(
            svg,
            cx + cell_w / 2.0,
            cy + 12.0,
            title,
            size=11,
          )
          scatter.svg_rect(cx, cy + 16.0, cell_w, cell_h - 16.0)
        } else {
          scatter.svg_rect(cx, cy, cell_w, cell_h)
        }
      LinePanel(title, line) =>
        if title.length() > 0 {
          plot_append_text(
            svg,
            cx + cell_w / 2.0,
            cy + 12.0,
            title,
            size=11,
          )
          line.svg_rect(cx, cy + 16.0, cell_w, cell_h - 16.0)
        } else {
          line.svg_rect(cx, cy, cell_w, cell_h)
        }
    }
    svg.write_string(body)
  }
  svg.write_string("</svg>\n")
  svg.to_string()
}

///|
/// Convenience: render a single histogram to a standalone SVG document.
pub fn histogram_to_svg(
  hist : Histogram,
  width~ : Int = 400,
  height~ : Int = 240,
  title~ : String = "",
) -> String {
  let page = DistributionPage::new(title~, rows=1, cols=1).add_panel(
    HistPanel("", hist),
  )
  page.to_svg(width~, height~)
}

///|
/// Convenience: render a single scatter plot to a standalone SVG document.
pub fn scatter_to_svg(
  scatter : ScatterPlot,
  width~ : Int = 400,
  height~ : Int = 240,
  title~ : String = "",
) -> String {
  let page = DistributionPage::new(title~, rows=1, cols=1).add_panel(
    ScatterPanel("", scatter),
  )
  page.to_svg(width~, height~)
}

///|
/// Convenience: render a single line plot to a standalone SVG document.
pub fn line_to_svg(
  line : LinePlot,
  width~ : Int = 400,
  height~ : Int = 240,
  title~ : String = "",
) -> String {
  let page = DistributionPage::new(title~, rows=1, cols=1).add_panel(
    LinePanel("", line),
  )
  page.to_svg(width~, height~)
}