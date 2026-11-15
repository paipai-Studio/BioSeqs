///|
/// Tests for bioframe_ebi.mbt — HTML/SVG interval visualization.

// ---- bf_render_intervals_html ---------------------------------------------

test "bf_render_intervals_html: empty intervals → valid HTML with no rects" {
  let chromsizes = [("chr1", 1000)]
  let html = @src.bf_render_intervals_html([], chromsizes)
  assert_true(html.has_prefix("<!doctype html>"))
  assert_true(html.contains("<svg"))
  assert_true(!html.contains("class=\"interval\""))
}

test "bf_render_intervals_html: one chrom one interval → rect present" {
  let chromsizes = [("chr1", 1000)]
  let intervals = [@src.BedInterval::new("chr1", 100, 500)]
  let html = @src.bf_render_intervals_html(intervals, chromsizes)
  assert_true(html.contains("class=\"interval\""))
  assert_true(html.contains("chr1"))
}

test "bf_render_intervals_html: title appears in h2" {
  let chromsizes = [("chr1", 1000)]
  let html = @src.bf_render_intervals_html([], chromsizes, title="My Figure Title")
  assert_true(html.contains("<h2>"))
  assert_true(html.contains("My Figure Title"))
}

test "bf_render_intervals_html: multi-chromosome tracks" {
  let chromsizes = [("chr1", 1000), ("chr2", 2000)]
  let intervals = [
    @src.BedInterval::new("chr1", 100, 500),
    @src.BedInterval::new("chr2", 800, 1200),
  ]
  let html = @src.bf_render_intervals_html(intervals, chromsizes)
  assert_true(html.contains("chr1"))
  assert_true(html.contains("chr2"))
  assert_true(html.contains("class=\"interval\""))
}

test "bf_render_intervals_html: intervals on unknown chrom → skipped" {
  let chromsizes = [("chr1", 1000)]
  let intervals = [@src.BedInterval::new("chr99", 100, 500)]
  let html = @src.bf_render_intervals_html(intervals, chromsizes)
  assert_true(!html.contains("class=\"interval\""))
}

test "bf_render_intervals_html: produces self-closing HTML" {
  let html = @src.bf_render_intervals_html([], [("chr1", 100)])
  assert_true(html.has_suffix("</html>\n"))
}
