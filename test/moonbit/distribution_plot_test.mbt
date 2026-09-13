///|
/// Tests for the Biopython Bio.Graphics.DistributionPlot-compatible renderer.

///|
test "histogram bins count observations" {
  let hist = @src.Histogram::from_data(
    [1.0, 2.0, 3.0, 4.0, 5.0],
    num_bins=5,
  )
  assert_eq(hist.bins.length(), 5)
  for bin in hist.bins {
    assert_eq(bin.count, 1)
  }
  assert_eq(hist.total, 5)
}

///|
test "histogram collapses constant data" {
  let hist = @src.Histogram::from_data([7.0, 7.0, 7.0], num_bins=10)
  assert_eq(hist.bins.length(), 1)
  assert_eq(hist.bins[0].count, 3)
}

///|
test "histogram max count" {
  let hist = @src.Histogram::from_data(
    [0.0, 0.0, 0.0, 1.0, 2.0],
    num_bins=3,
  )
  assert_eq(hist.max_count(), 3)
}

///|
test "histogram empty data" {
  let hist = @src.Histogram::from_data([], num_bins=5)
  assert_eq(hist.bins.length(), 0)
  assert_eq(hist.total, 0)
}

///|
test "scatter plot svg contains circles" {
  let scatter = @src.ScatterPlot::new(series=[
    @src.PlotSeries::new(
      points=[@src.PlotPoint::new(x=0.0, y=0.0), @src.PlotPoint::new(
        x=1.0,
        y=1.0,
      )],
    ),
  ])
  let svg = @src.scatter_to_svg(scatter, width=200, height=200)
  assert_true(svg.contains("<svg"))
  assert_true(svg.contains("</svg>"))
  assert_true(svg.contains("<circle"))
}

///|
test "line plot svg contains polyline" {
  let line = @src.LinePlot::new(series=[
    @src.PlotSeries::new(
      points=[
        @src.PlotPoint::new(x=0.0, y=0.0),
        @src.PlotPoint::new(x=1.0, y=1.0),
        @src.PlotPoint::new(x=2.0, y=0.5),
      ],
      color="#ff0000",
    ),
  ])
  let svg = @src.line_to_svg(line, width=300, height=200)
  assert_true(svg.contains("<polyline"))
  assert_true(svg.contains("#ff0000"))
}

///|
test "histogram svg contains rects" {
  let hist = @src.Histogram::from_data(
    [1.0, 2.0, 3.0, 4.0, 5.0],
    num_bins=5,
  )
  let svg = @src.histogram_to_svg(hist, width=300, height=200)
  assert_true(svg.contains("<rect"))
  assert_true(svg.contains("<line"))
}

///|
test "distribution page lays out multiple panels" {
  let hist = @src.Histogram::from_data([1.0, 2.0, 3.0], num_bins=3)
  let scatter = @src.ScatterPlot::new(series=[
    @src.PlotSeries::new(points=[@src.PlotPoint::new(x=0.0, y=0.0)]),
  ])
  let page = @src.DistributionPage::new(
    title="Demo",
    rows=1,
    cols=2,
  ).add_panel(@src.PlotPanel::HistPanel("Hist", hist)).add_panel(
    @src.PlotPanel::ScatterPanel("Scatter", scatter),
  )
  let svg = page.to_svg(width=600, height=300)
  assert_true(svg.contains("Demo"))
  assert_true(svg.contains("Hist"))
  assert_true(svg.contains("Scatter"))
  assert_true(svg.contains("<rect"))
  assert_true(svg.contains("<circle"))
}

///|
test "plot series default color" {
  let series = @src.PlotSeries::new(points=[])
  assert_eq(series.color, "#1f77b4")
  assert_eq(series.label, "")
}