///|
fn main {
  println("=== Bio.Graphics.DistributionPlot Demo ===")

  println("\n1. Histogram of gene expression values")
  let data = [
    1.2, 1.5, 2.1, 2.3, 2.8, 3.1, 3.4, 3.9, 4.2, 4.8, 5.1, 5.5, 6.2, 6.8, 7.1,
    7.5, 8.2, 8.8, 9.1, 9.5,
  ]
  let hist = @src.Histogram::from_data(data, num_bins=8)
  println("  Bins: " + hist.bins.length().to_string() + ", total: " + hist.total.to_string())
  for bin in hist.bins {
    let bar = "█".repeat(bin.count)
    println(
      "  [" +
      bin.lo.to_string() +
      ", " +
      bin.hi.to_string() +
      "): " +
      bin.count.to_string() +
      " " +
      bar,
    )
  }
  let hist_svg = @src.histogram_to_svg(hist, width=400, height=240, title="Expression histogram")
  println("  SVG length: " + hist_svg.length().to_string() + " chars")

  println("\n2. Scatter plot of two genes")
  let scatter = @src.ScatterPlot::new(series=[
    @src.PlotSeries::new(
      points=[
        @src.PlotPoint::new(x=1.0, y=2.0),
        @src.PlotPoint::new(x=2.0, y=3.5),
        @src.PlotPoint::new(x=3.0, y=3.0),
        @src.PlotPoint::new(x=4.0, y=5.0),
        @src.PlotPoint::new(x=5.0, y=4.5),
      ],
      color="#1f77b4",
      label="Gene A",
    ),
    @src.PlotSeries::new(
      points=[
        @src.PlotPoint::new(x=1.0, y=1.0),
        @src.PlotPoint::new(x=2.0, y=1.5),
        @src.PlotPoint::new(x=3.0, y=2.5),
        @src.PlotPoint::new(x=4.0, y=2.0),
        @src.PlotPoint::new(x=5.0, y=3.5),
      ],
      color="#ff7f0e",
      label="Gene B",
    ),
  ])
  let scatter_svg = @src.scatter_to_svg(scatter, width=400, height=240, title="Gene scatter")
  println("  SVG length: " + scatter_svg.length().to_string() + " chars")

  println("\n3. Line plot of a time course")
  let line = @src.LinePlot::new(series=[
    @src.PlotSeries::new(
      points=[
        @src.PlotPoint::new(x=0.0, y=1.0),
        @src.PlotPoint::new(x=1.0, y=2.0),
        @src.PlotPoint::new(x=2.0, y=1.8),
        @src.PlotPoint::new(x=3.0, y=3.0),
        @src.PlotPoint::new(x=4.0, y=2.5),
        @src.PlotPoint::new(x=5.0, y=4.0),
      ],
      color="#2ca02c",
      label="Treatment",
    ),
  ])
  let line_svg = @src.line_to_svg(line, width=400, height=240, title="Time course")
  println("  SVG length: " + line_svg.length().to_string() + " chars")

  println("\n4. Distribution page with multiple panels")
  let page = @src.DistributionPage::new(
    title="QC report",
    rows=1,
    cols=2,
  ).add_panel(@src.PlotPanel::HistPanel("Histogram", hist)).add_panel(
    @src.PlotPanel::ScatterPanel("Scatter", scatter),
  )
  let page_svg = page.to_svg(width=600, height=300)
  println("  Page SVG length: " + page_svg.length().to_string() + " chars")
  println("  Contains <svg: " + page_svg.contains("<svg").to_string())
  println("  Contains </svg>: " + page_svg.contains("</svg>").to_string())
}