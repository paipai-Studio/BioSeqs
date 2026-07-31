///|
/// Gviz module - Genome visualization track system.
/// Inspired by Bioconductor Gviz package.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/Gviz.html
///   - Hahne F, Durinck S (2024) Gviz: Plotting data and annotation information
///     along genomic coordinates. R package version 1.52.0.

///|
/// Track type for Gviz display.
pub enum GvizTrackType {
  AnnotationTrack
  GeneRegionTrack
  DataTrack
  IdeogramTrack
  GenomeAxisTrack
  SequenceTrack
}

///|
/// Create an AnnotationTrack type.
pub fn track_type_annotation() -> GvizTrackType {
  GvizTrackType::AnnotationTrack
}

///|
/// Create a GeneRegionTrack type.
pub fn track_type_gene_region() -> GvizTrackType {
  GvizTrackType::GeneRegionTrack
}

///|
/// Create a DataTrack type.
pub fn track_type_data() -> GvizTrackType {
  GvizTrackType::DataTrack
}

///|
/// Create an IdeogramTrack type.
pub fn gviz_track_type_ideogram() -> GvizTrackType {
  GvizTrackType::IdeogramTrack
}

///|
/// Create a GenomeAxisTrack type.
pub fn track_type_genome_axis() -> GvizTrackType {
  GvizTrackType::GenomeAxisTrack
}

///|
/// Create a SequenceTrack type.
pub fn track_type_sequence() -> GvizTrackType {
  GvizTrackType::SequenceTrack
}

///|
/// Display strand.
pub enum GvizStrand {
  Forward
  Reverse
  Unstranded
}

///|
/// Create Forward strand.
pub fn track_strand_forward() -> GvizStrand {
  GvizStrand::Forward
}

///|
/// Create Reverse strand.
pub fn track_strand_reverse() -> GvizStrand {
  GvizStrand::Reverse
}

///|
/// Create Unstranded.
pub fn track_strand_unstranded() -> GvizStrand {
  GvizStrand::Unstranded
}

///|
/// A genomic feature for annotation tracks.
pub struct GvizFeature {
  feature_id : String
  chromosome : String
  start : Int
  end_ : Int
  strand : GvizStrand
  feature_type : String
  label : String
}

///|
/// A genome track for visualization.
pub struct GvizTrack {
  track_name : String
  track_type : GvizTrackType
  chromosome : String
  start : Int
  end_ : Int
  features : Array[GvizFeature]
  data_values : Array[(Int, Double)]
  mut color : String
  mut display_label : String
}

///|
/// A genome region to display.
pub struct GvizRegion {
  chromosome : String
  start : Int
  end_ : Int
}

///|
/// The complete genome plot.
pub struct GvizPlot {
  tracks : Array[GvizTrack]
  region : GvizRegion
  title : String
  width : Int
  height : Int
}

///|
/// TrackType to string.
pub fn GvizTrackType::to_string(self : GvizTrackType) -> String {
  match self {
    GvizTrackType::AnnotationTrack => "AnnotationTrack"
    GvizTrackType::GeneRegionTrack => "GeneRegionTrack"
    GvizTrackType::DataTrack => "DataTrack"
    GvizTrackType::IdeogramTrack => "IdeogramTrack"
    GvizTrackType::GenomeAxisTrack => "GenomeAxisTrack"
    GvizTrackType::SequenceTrack => "SequenceTrack"
  }
}

///|
/// TrackStrand to string.
pub fn GvizStrand::to_string(self : GvizStrand) -> String {
  match self {
    GvizStrand::Forward => "+"
    GvizStrand::Reverse => "-"
    GvizStrand::Unstranded => "*"
  }
}

///|
/// Repeat a string n times.
fn gviz_repeat(s : String, n : Int) -> String {
  let mut result = ""
  let mut i = 0
  while i < n {
    result = result + s
    i = i + 1
  }
  result
}

///|
/// Create a GvizFeature.
pub fn gviz_feature(
  feature_id : String,
  chromosome : String,
  start : Int,
  end_ : Int,
  strand : GvizStrand,
  feature_type : String,
  label : String,
) -> GvizFeature {
  GvizFeature::{ feature_id, chromosome, start, end_, strand, feature_type, label }
}

///|
/// Create a GvizTrack.
pub fn gviz_track(
  track_name : String,
  track_type : GvizTrackType,
  chromosome : String,
  start : Int,
  end_ : Int,
) -> GvizTrack {
  GvizTrack::{
    track_name,
    track_type,
    chromosome,
    start,
    end_,
    features: Array::new(),
    data_values: Array::new(),
    color: "black",
    display_label: track_name,
  }
}

///|
/// Set track color.
pub fn GvizTrack::set_color(self : GvizTrack, color : String) -> Unit {
  self.color = color
}

///|
/// Set display label.
pub fn GvizTrack::set_label(self : GvizTrack, label : String) -> Unit {
  self.display_label = label
}

///|
/// Add a feature to a track.
pub fn GvizTrack::add_feature(self : GvizTrack, feature : GvizFeature) -> Unit {
  self.features.push(feature)
}

///|
/// Add a data point to a DataTrack.
pub fn GvizTrack::add_data_point(self : GvizTrack, position : Int, value : Double) -> Unit {
  self.data_values.push((position, value))
}

///|
/// Get number of features.
pub fn GvizTrack::get_n_features(self : GvizTrack) -> Int {
  self.features.length()
}

///|
/// Get number of data points.
pub fn GvizTrack::get_n_data_points(self : GvizTrack) -> Int {
  self.data_values.length()
}

///|
/// Get features within a specific region.
pub fn GvizTrack::get_features_in_region(self : GvizTrack, start : Int, end_ : Int) -> Array[GvizFeature] {
  let result : Array[GvizFeature] = Array::new()
  let mut i = 0
  while i < self.features.length() {
    let f = self.features[i]
    if f.start <= end_ && f.end_ >= start {
      result.push(f)
    }
    i = i + 1
  }
  result
}

///|
/// Get data values within a specific region.
pub fn GvizTrack::get_data_in_region(self : GvizTrack, start : Int, end_ : Int) -> Array[(Int, Double)] {
  let result : Array[(Int, Double)] = Array::new()
  let mut i = 0
  while i < self.data_values.length() {
    let (pos, val) = self.data_values[i]
    if pos >= start && pos <= end_ {
      result.push((pos, val))
    }
    i = i + 1
  }
  result
}

///|
/// Create a genome region.
pub fn gviz_region(chromosome : String, start : Int, end_ : Int) -> GvizRegion {
  GvizRegion::{ chromosome, start, end_ }
}

///|
/// Create a GvizPlot.
pub fn gviz_plot(
  region : GvizRegion,
  title? : String = "",
  width? : Int = 80,
  height? : Int = 25,
) -> GvizPlot {
  GvizPlot::{
    tracks: Array::new(),
    region,
    title,
    width,
    height,
  }
}

///|
/// Add a track to the plot.
pub fn GvizPlot::add_track(self : GvizPlot, track : GvizTrack) -> Unit {
  self.tracks.push(track)
}

///|
/// Get number of tracks.
pub fn GvizPlot::get_n_tracks(self : GvizPlot) -> Int {
  self.tracks.length()
}

///|
/// Get the plot region.
pub fn GvizPlot::get_region(self : GvizPlot) -> GvizRegion {
  self.region
}

///|
/// Render the plot as ASCII text.
pub fn GvizPlot::to_ascii(self : GvizPlot) -> String {
  let w = self.width
  let region_start = self.region.start
  let region_end = self.region.end_
  let region_len = region_end - region_start
  if region_len <= 0 {
    return "Invalid region: start >= end"
  }

  let mut result = ""
  result = result + self.title + "\n"
  result = result + "Chromosome: " + self.region.chromosome + " [" + region_start.to_string() + " - " + region_end.to_string() + "]\n"
  result = result + gviz_repeat("=", w) + "\n"

  // Process each track
  let mut ti = 0
  while ti < self.tracks.length() {
    let track = self.tracks[ti]
    result = result + "[" + track.track_type.to_string() + "] " + track.display_label + "\n"

    // Render track content
    match track.track_type {
      GvizTrackType::GenomeAxisTrack => {
        // Draw axis with position labels
        result = result + gviz_render_axis(region_start, region_end, w) + "\n"
      }
      GvizTrackType::AnnotationTrack | GvizTrackType::GeneRegionTrack => {
        result = result + gviz_render_features(track, region_start, region_end, w) + "\n"
      }
      GvizTrackType::DataTrack => {
        result = result + gviz_render_data(track, region_start, region_end, w) + "\n"
      }
      GvizTrackType::IdeogramTrack => {
        result = result + gviz_render_ideogram(track, region_start, region_end, w) + "\n"
      }
      GvizTrackType::SequenceTrack => {
        result = result + gviz_render_sequence(track, region_start, region_end, w) + "\n"
      }
    }
    result = result + "\n"
    ti = ti + 1
  }

  result
}

///|
/// Render genome axis.
fn gviz_render_axis(start : Int, end_ : Int, width : Int) -> String {
  let mut s = ""
  let len = end_ - start
  let step = len / (width - 1)
  if step < 1 {
    return "|" + gviz_repeat("-", width - 2) + "|"
  }

  // Draw the axis line
  s = s + "|"
  let mut i = 1
  while i < width - 1 {
    let pos = start + i * step
    if pos % (step * 5) < step {
      s = s + "+"
    } else {
      s = s + "-"
    }
    i = i + 1
  }
  s = s + "|"
  s = s + "\n"

  // Add position labels
  s = s + start.to_string()
  let label_end = end_.to_string()
  let padding = width - start.to_string().length() - label_end.length()
  let mut pi = 0
  while pi < padding {
    s = s + " "
    pi = pi + 1
  }
  s = s + label_end

  s
}

///|
/// Render features as boxes.
fn gviz_render_features(track : GvizTrack, region_start : Int, region_end : Int, width : Int) -> String {
  let region_len = (region_end - region_start).to_double()
  if region_len <= 0.0 {
    return ""
  }

  // Build the rendering grid
  let grid : Array[String] = Array::new()
  let mut i = 0
  while i < width {
    grid.push(" ")
    i = i + 1
  }

  // Place features
  let mut fi = 0
  while fi < track.features.length() {
    let f = track.features[fi]
    if f.end_ < region_start || f.start > region_end {
      fi = fi + 1
    } else {
      let rel_start = ((f.start - region_start).to_double() / region_len * (width - 1).to_double()).to_int()
      let rel_end = ((f.end_ - region_start).to_double() / region_len * (width - 1).to_double()).to_int()

      let s = if rel_start < 0 { 0 } else { rel_start }
      let e = if rel_end >= width { width - 1 } else { rel_end }

      let ch = match f.feature_type {
        "exon" => "█"
        "intron" => "-"
        "cds" => "█"
        "utr" => "░"
        "promoter" => ">"
        "enhancer" => "*"
        "gene" => "█"
        _ => "="
      }

      let mut ci = s
      while ci <= e {
        if ci >= 0 && ci < width {
          grid.set(ci, ch)
        }
        ci = ci + 1
      }

      // Add strand marker
      let mid = (s + e) / 2
      if mid >= 0 && mid < width {
        if f.strand is GvizStrand::Forward {
          // Don't overwrite if already set
        } else if f.strand is GvizStrand::Reverse {
          // Don't overwrite
        }
      }

      fi = fi + 1
    }
  }

  let row = grid.join("")

  // Build label line
  let label_row = gviz_render_labels(track, region_start, region_end, width)

  row + "\n" + label_row
}

///|
/// Render feature labels.
fn gviz_render_labels(track : GvizTrack, region_start : Int, region_end : Int, width : Int) -> String {
  let region_len = (region_end - region_start).to_double()
  if region_len <= 0.0 {
    return ""
  }

  let grid : Array[String] = Array::new()
  let mut i = 0
  while i < width {
    grid.push(" ")
    i = i + 1
  }

  let mut fi = 0
  while fi < track.features.length() {
    let f = track.features[fi]
    if f.end_ < region_start || f.start > region_end {
      fi = fi + 1
    } else {
      let mid_pos = (f.start + f.end_) / 2
      let rel_mid = ((mid_pos - region_start).to_double() / region_len * (width - 1).to_double()).to_int()
      if rel_mid >= 0 && rel_mid < width && f.label.length() > 0 {
        let label = f.label
        let mut li = 0
        while li < label.length() && rel_mid + li < width {
          grid.set(rel_mid + li, label.substring(start=li, end=li + 1))
          li = li + 1
        }
      }
      fi = fi + 1
    }
  }

  grid.join("")
}

///|
/// Render data track as a simple bar plot.
fn gviz_render_data(track : GvizTrack, region_start : Int, region_end : Int, width : Int) -> String {
  let data = track.get_data_in_region(region_start, region_end)
  if data.length() == 0 {
    return "(no data)"
  }

  // Find min and max
  let mut min_val = data[0].1
  let mut max_val = data[0].1
  let mut i = 1
  while i < data.length() {
    if data[i].1 < min_val { min_val = data[i].1 }
    if data[i].1 > max_val { max_val = data[i].1 }
    i = i + 1
  }

  let val_range = max_val - min_val
  if val_range < 1.0e-10 {
    return "(constant data: " + min_val.to_string() + ")"
  }

  let region_len = (region_end - region_start).to_double()

  // Build bar chart
  let bar_height = 5
  let grid : Array[Array[String]] = Array::new()
  let mut ri = 0
  while ri < bar_height {
    let row : Array[String] = Array::new()
    let mut ci = 0
    while ci < width {
      row.push(" ")
      ci = ci + 1
    }
    grid.push(row)
    ri = ri + 1
  }

  // Place data points as bars
  let mut di = 0
  while di < data.length() {
    let (pos, val) = data[di]
    let rel_pos = ((pos - region_start).to_double() / region_len * (width - 1).to_double()).to_int()
    let normalized = (val - min_val) / val_range
    let bar_h = (normalized * bar_height.to_double()).to_int()

    if rel_pos >= 0 && rel_pos < width {
      let mut hi = 0
      while hi < bar_h && hi < bar_height {
        grid[bar_height - 1 - hi].set(rel_pos, "█")
        hi = hi + 1
      }
    }
    di = di + 1
  }

  let mut result = ""
  let mut rri = 0
  while rri < bar_height {
    result = result + grid[rri].join("") + "\n"
    rri = rri + 1
  }

  result = result + "min: " + min_val.to_string() + ", max: " + max_val.to_string()
  result
}

///|
/// Render ideogram track.
fn gviz_render_ideogram(track : GvizTrack, region_start : Int, region_end : Int, width : Int) -> String {
  // Simplified ideogram: draw a chromosome bar with centromere marker
  let mut s = ""
  let mut i = 0
  while i < width {
    if i == width / 2 {
      s = s + "."
    } else if i < width / 2 {
      s = s + "="
    } else {
      s = s + "="
    }
    i = i + 1
  }
  s = s + "\n"
  s = s + "p" + gviz_repeat(" ", width / 2 - 2) + "cen" + gviz_repeat(" ", width / 2 - 2) + "q"
  s
}

///|
/// Render sequence track.
fn gviz_render_sequence(track : GvizTrack, region_start : Int, region_end : Int, width : Int) -> String {
  let region_len = region_end - region_start
  if region_len <= 0 {
    return ""
  }

  let mut s = ""
  let chars = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
  let display_len = if region_len < width { region_len } else { width }

  let mut i = 0
  while i < display_len {
    let idx = (region_start + i) % chars.length()
    s = s + chars.substring(start=idx, end=idx + 1)
    i = i + 1
  }

  s
}

///|
/// Generate summary string.
pub fn GvizPlot::summary(self : GvizPlot) -> String {
  let mut s = "=== Gviz Plot Summary ===\n".to_string()
  s = s + "Title: " + self.title + "\n"
  s = s + "Region: " + self.region.chromosome + ":" + self.region.start.to_string() + "-" + self.region.end_.to_string() + "\n"
  s = s + "Tracks: " + self.tracks.length().to_string() + "\n"
  let mut i = 0
  while i < self.tracks.length() {
    s = s + "  [" + (i + 1).to_string() + "] " + self.tracks[i].track_type.to_string() + ": " + self.tracks[i].display_label
    s = s + " (" + self.tracks[i].get_n_features().to_string() + " features, " + self.tracks[i].get_n_data_points().to_string() + " data points)\n"
    i = i + 1
  }
  s
}

///|
/// Create sample data for testing.
pub fn gviz_sample_plot() -> GvizPlot {
  let region = gviz_region("chr1", 1000, 5000)
  let plot = gviz_plot(region, title="Sample Genome Region", width=80, height=25)

  // Axis track
  let axis_track = gviz_track("axis", track_type_genome_axis(), "chr1", 1000, 5000)
  axis_track.set_label("Genome Axis")
  plot.add_track(axis_track)

  // Gene region track
  let gene_track = gviz_track("genes", track_type_gene_region(), "chr1", 1000, 5000)
  gene_track.set_label("Gene Annotations")
  gene_track.set_color("blue")
  gene_track.add_feature(gviz_feature("g1", "chr1", 1200, 1800, track_strand_forward(), "exon", "BRCA1"))
  gene_track.add_feature(gviz_feature("g2", "chr1", 2000, 3500, track_strand_reverse(), "exon", "TP53"))
  gene_track.add_feature(gviz_feature("g3", "chr1", 4000, 4500, track_strand_forward(), "exon", "EGFR"))
  plot.add_track(gene_track)

  // Data track
  let data_track = gviz_track("coverage", track_type_data(), "chr1", 1000, 5000)
  data_track.set_label("Coverage")
  data_track.set_color("red")
  let mut i = 0
  while i < 40 {
    let pos = 1000 + i * 100
    let val = 10.0 + @math.sin(i.to_double() * 0.5) * 5.0 + i.to_double() * 0.3
    data_track.add_data_point(pos, val)
    i = i + 1
  }
  plot.add_track(data_track)

  // Ideogram track
  let ideo_track = gviz_track("ideogram", gviz_track_type_ideogram(), "chr1", 1000, 5000)
  ideo_track.set_label("Ideogram")
  plot.add_track(ideo_track)

  plot
}
