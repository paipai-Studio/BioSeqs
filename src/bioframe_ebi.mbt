/// bioframe_ebi.mbt — End-user Bioframe convenience utilities.
///
/// - bf_to_bigwig:  encode a bigWig file into bytes
/// - bf_to_bigbed:  encode a bigBed file into bytes
/// - bf_render_intervals_html:  produce a self-contained SVG track image
///
/// The heavy lifting (bigWig/bigBed encoding, zlib compression, R-tree
/// and B+ tree construction) lives in bigwig.mbt / bigbed.mbt; this module
/// is a thin ergonomic wrapper so callers can do:
///
///     let bytes = bf_to_bigwig(targets, intervals)
///     // ... write `bytes` to disk using whatever I/O layer you have ...
///     let html = bf_render_intervals_html(intervals, chromsizes)

// ================================================================
// 1. to_bigwig / to_bigbed  encoders
// ================================================================

///|
/// Encode a bigWig file into bytes.
///
/// Returns `Array[Int32]` (one byte per element).  Caller is responsible
/// for writing the bytes to disk or memory — this avoids coupling to any
/// particular I/O layer.  Raises BigWigError on invalid inputs.
pub fn bf_to_bigwig(
  targets : Array[BigBedTarget],
  intervals : Array[BigWigInterval],
) -> Array[Int32] raise BigWigError {
  bigwig_write(targets, intervals)
}

///|
/// Encode a bigBed file into bytes from minimal interval rows.
///
/// Internally delegates to BigBedRecord::create + bigbed_write.
/// Raises BigBedError on invalid inputs.
pub fn bf_to_bigbed(
  targets : Array[BigBedTarget],
  chroms : Array[String],
  starts : Array[Int32],
  ends : Array[Int32],
  names~ : Array[String] = [],
  scores~ : Array[Int32] = [],
  strands~ : Array[String] = [],
) -> Array[Int32] raise BigBedError {
  let n = chroms.length()
  let records : Array[BigBedRecord] = Array(capacity=n)
  for i in 0..<n {
    let name = if i < names.length() { names[i] } else { "." }
    let score = if i < scores.length() { scores[i] } else { 0 }
    let strand = if i < strands.length() { strands[i] } else { "." }
    records.push(
      BigBedRecord::create(
        chroms[i],
        starts[i],
        ends[i],
        name~,
        score~,
        strand~,
      ),
    )
  }
  bigbed_write(targets, records)
}

// ================================================================
// 2. Lightweight SVG / HTML visualization
// ================================================================

///|
/// Render BedInterval into a self-contained HTML page with inline SVG.
///
/// Returns a String (not bytes) so it can be inspected directly or saved
/// via whatever text I/O layer your platform exposes.
///
/// Parameters:
///   - intervals:            BedInterval array
///   - chromsizes:           [(chrom, length)] — used to set the x-axis
///                           for each track.  Intervals whose chrom is
///                           not present here are silently skipped.
///   - title:                optional figure caption
///   - width_px:             canvas width (default 900)
///   - height_per_chrom_px:  track height (default 24)
pub fn bf_render_intervals_html(
  intervals : Array[BedInterval],
  chromsizes : Array[(String, Int32)],
  title~ : String = "",
  width_px~ : Int32 = 900,
  height_per_chrom_px~ : Int32 = 24,
  padding_px~ : Int32 = 40,
) -> String {
  let n_chroms = chromsizes.length()
  let plot_h = n_chroms * height_per_chrom_px
  let total_h = padding_px * 2 + plot_h
  let total_w = width_px + padding_px * 2
  let y_base = padding_px + 10

  let sb = StringBuilder()
  sb.write("<!doctype html>\n<html><head>\n")
  sb.write("  <meta charset=\"utf-8\">\n")
  sb.write("  <title>Bioframe interval track</title>\n")
  sb.write("  <style>\n")
  sb.write("    body { font-family: sans-serif; padding: 16px; }\n")
  sb.write(
    "    .chrom-label { font-size: 11px; fill: #444; dominant-baseline: middle; }\n",
  )
  sb.write("    .chrom-axis  { stroke: #bbb; stroke-width: 1; }\n")
  sb.write("    .interval    { fill: #4c78a8; }\n")
  sb.write("  </style>\n")
  sb.write("</head><body>\n")
  if title.length() > 0 {
    sb.write("  <h2>")
    sb.write(title)
    sb.write("</h2>\n")
  }
  sb.write("  <svg width=\"")
  sb.write(total_w.to_string())
  sb.write("\" height=\"")
  sb.write(total_h.to_string())
  sb.write("\" xmlns=\"http://www.w3.org/2000/svg\">\n")

  // Chromosome labels + axes
  for i in 0..<n_chroms {
    let cs = chromsizes[i]
    let y = y_base + i * height_per_chrom_px
    sb.write("    <text x=\"8\" y=\"")
    sb.write((y + height_per_chrom_px / 2).to_string())
    sb.write("\" class=\"chrom-label\">")
    sb.write(cs.0)
    sb.write("</text>\n")
    sb.write("    <line x1=\"")
    sb.write(padding_px.to_string())
    sb.write("\" y1=\"")
    sb.write(y.to_string())
    sb.write("\" x2=\"")
    sb.write((padding_px + width_px).to_string())
    sb.write("\" y2=\"")
    sb.write(y.to_string())
    sb.write("\" class=\"chrom-axis\" />\n")
  }

  // Interval rectangles
  for iv in intervals {
    let mut row_idx = -1
    for j in 0..<n_chroms {
      if chromsizes[j].0 == iv.chrom {
        row_idx = j
        break
      }
    }
    if row_idx < 0 {
      continue
    }
    let chrom_len = chromsizes[row_idx].1
    if chrom_len <= 0 {
      continue
    }
    let y = y_base + row_idx * height_per_chrom_px
    let w64 = width_px.to_double()
    let cl64 = chrom_len.to_double()
    let sx64 = iv.start.to_double()
    let ex64 = iv.end.to_double()
    let p64 = padding_px.to_double()
    let x_start_d = p64 + sx64 / cl64 * w64
    let x_end_d = p64 + ex64 / cl64 * w64
    let w_d = x_end_d - x_start_d
    let w_clamped_d = if w_d < 1.0 { 1.0 } else { w_d }
    let rect_y = y + 4
    let rect_h = height_per_chrom_px - 8
    let xs = x_start_d.to_int().to_string()
    let xe = w_clamped_d.to_int().to_string()
    let ry = rect_y.to_string()
    let rh = rect_h.to_string()
    sb.write("    <rect x=\"")
    sb.write(xs)
    sb.write("\" y=\"")
    sb.write(ry)
    sb.write("\" width=\"")
    sb.write(xe)
    sb.write("\" height=\"")
    sb.write(rh)
    sb.write("\" class=\"interval\" />\n")
  }

  sb.write("  </svg>\n")
  sb.write("</body></html>\n")
  sb.to_string()
}
