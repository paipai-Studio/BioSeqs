///|
/// karyoploteR module - Karyotype visualization for genomic data.
/// Inspired by Bioconductor karyoploteR package.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/karyoploteR.html
///   - Gel B, Fernández E. (2024) karyoploteR: Plotting of karyotype regions with
///     additional data. R package version 1.22.0.

///|
/// Data track types for karyotype plots.
pub enum TrackType {
  Points
  Lines
  Bars
  Heatmap
  Ideogram
}

///|
/// Create a Points track type.
pub fn track_type_points() -> TrackType {
  TrackType::Points
}

///|
/// Create a Lines track type.
pub fn track_type_lines() -> TrackType {
  TrackType::Lines
}

///|
/// Create a Bars track type.
pub fn track_type_bars() -> TrackType {
  TrackType::Bars
}

///|
/// Create a Heatmap track type.
pub fn track_type_heatmap() -> TrackType {
  TrackType::Heatmap
}

///|
/// Create an Ideogram track type.
pub fn track_type_ideogram() -> TrackType {
  TrackType::Ideogram
}

///|
/// A genomic region for plotting.
pub struct KaryotypeRegion {
  chromosome : String
  start : Int
  end : Int
  label : String
  color : String
}

///|
/// A data point on a genomic track.
pub struct TrackPoint {
  chromosome : String
  position : Int
  value : Double
  label : String
}

///|
/// A karyotype data track.
pub struct KaryotypeTrack {
  track_id : String
  track_type : TrackType
  chromosome : String
  start : Int
  end : Int
  data : Array[TrackPoint]
  color : String
  y_min : Double
  y_max : Double
  label : String
}

///|
/// A chromosome ideogram representation.
pub struct IdeogramBand {
  chromosome : String
  start : Int
  end : Int
  band_name : String
  stain : String
}

///|
/// A karyotype plot with chromosomes and tracks.
pub struct KaryotypePlot {
  genome : String
  chromosomes : Array[String]
  chromosome_sizes : Map[String, Int]
  tracks : Array[KaryotypeTrack]
  regions : Array[KaryotypeRegion]
  ideograms : Array[IdeogramBand]
  title : String
}

///|
/// Track type to string.
pub fn TrackType::to_string(self : TrackType) -> String {
  match self {
    TrackType::Points => "points"
    TrackType::Lines => "lines"
    TrackType::Bars => "bars"
    TrackType::Heatmap => "heatmap"
    TrackType::Ideogram => "ideogram"
  }
}

///|
/// Create a new karyotype plot for a genome.
pub fn KaryotypePlot::new(genome : String) -> KaryotypePlot {
  let chromosomes : Array[String] = Array::new()
  let chromosome_sizes : Map[String, Int] = Map([], capacity=25)
  let tracks : Array[KaryotypeTrack] = Array::new()
  let regions : Array[KaryotypeRegion] = Array::new()
  let ideograms : Array[IdeogramBand] = Array::new()

  // Set up default human chromosomes
  if genome == "hg38" || genome == "hg19" || genome == "human" {
    let chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
      "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
      "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    let sizes = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
      159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
      114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616,
      64444167, 46709983, 50818468, 156040895, 57227415]

    let mut i = 0
    while i < chroms.length() {
      chromosomes.push(chroms[i])
      chromosome_sizes[chroms[i]] = sizes[i]
      i = i + 1
    }
  }

  KaryotypePlot::{
    genome,
    chromosomes,
    chromosome_sizes,
    tracks,
    regions,
    ideograms,
    title: "Karyotype Plot".to_string(),
  }
}

///|
/// Get the genome name.
pub fn KaryotypePlot::get_genome(plot : KaryotypePlot) -> String {
  plot.genome
}

///|
/// Get the list of chromosomes.
pub fn KaryotypePlot::get_chromosomes(plot : KaryotypePlot) -> Array[String] {
  plot.chromosomes
}

///|
/// Get chromosome size.
pub fn KaryotypePlot::get_chromosome_size(plot : KaryotypePlot, chr : String) -> Int {
  match plot.chromosome_sizes.get(chr) {
    Some(s) => s
    None => 0
  }
}

///|
/// Get the number of tracks.
pub fn KaryotypePlot::get_n_tracks(plot : KaryotypePlot) -> Int {
  plot.tracks.length()
}

///|
/// Get the number of regions.
pub fn KaryotypePlot::get_n_regions(plot : KaryotypePlot) -> Int {
  plot.regions.length()
}

///|
/// Add a track to the karyotype plot.
pub fn KaryotypePlot::add_track(plot : KaryotypePlot, track : KaryotypeTrack) -> KaryotypePlot {
  let tracks = plot.tracks.copy()
  tracks.push(track)
  KaryotypePlot::{
    genome: plot.genome,
    chromosomes: plot.chromosomes,
    chromosome_sizes: plot.chromosome_sizes,
    tracks,
    regions: plot.regions,
    ideograms: plot.ideograms,
    title: plot.title,
  }
}

///|
/// Add a region to the karyotype plot.
pub fn KaryotypePlot::add_region(plot : KaryotypePlot, region : KaryotypeRegion) -> KaryotypePlot {
  let regions = plot.regions.copy()
  regions.push(region)
  KaryotypePlot::{
    genome: plot.genome,
    chromosomes: plot.chromosomes,
    chromosome_sizes: plot.chromosome_sizes,
    tracks: plot.tracks,
    regions,
    ideograms: plot.ideograms,
    title: plot.title,
  }
}

///|
/// Create a new track.
pub fn KaryotypeTrack::new(track_id : String, track_type : TrackType, chromosome : String, start : Int, end : Int, label : String) -> KaryotypeTrack {
  KaryotypeTrack::{
    track_id,
    track_type,
    chromosome,
    start,
    end,
    data: Array::new(),
    color: "#3366CC".to_string(),
    y_min: 0.0,
    y_max: 1.0,
    label,
  }
}

///|
/// Add a data point to a track.
pub fn KaryotypeTrack::add_point(track : KaryotypeTrack, point : TrackPoint) -> KaryotypeTrack {
  let data = track.data.copy()
  data.push(point)
  KaryotypeTrack::{
    track_id: track.track_id,
    track_type: track.track_type,
    chromosome: track.chromosome,
    start: track.start,
    end: track.end,
    data,
    color: track.color,
    y_min: track.y_min,
    y_max: track.y_max,
    label: track.label,
  }
}

///|
/// Create a new track point.
pub fn TrackPoint::new(chromosome : String, position : Int, value : Double, label : String) -> TrackPoint {
  TrackPoint::{ chromosome, position, value, label }
}

///|
/// Create a new karyotype region.
pub fn KaryotypeRegion::new(chromosome : String, start : Int, end : Int, label : String, color : String) -> KaryotypeRegion {
  KaryotypeRegion::{ chromosome, start, end, label, color }
}

///|
/// Create a new ideogram band.
pub fn IdeogramBand::new(chromosome : String, start : Int, end : Int, band_name : String, stain : String) -> IdeogramBand {
  IdeogramBand::{ chromosome, start, end, band_name, stain }
}

///|
/// Set track color.
pub fn KaryotypeTrack::set_color(track : KaryotypeTrack, color : String) -> KaryotypeTrack {
  KaryotypeTrack::{
    track_id: track.track_id,
    track_type: track.track_type,
    chromosome: track.chromosome,
    start: track.start,
    end: track.end,
    data: track.data,
    color,
    y_min: track.y_min,
    y_max: track.y_max,
    label: track.label,
  }
}

///|
/// Set track y-range.
pub fn KaryotypeTrack::set_y_range(track : KaryotypeTrack, y_min : Double, y_max : Double) -> KaryotypeTrack {
  KaryotypeTrack::{
    track_id: track.track_id,
    track_type: track.track_type,
    chromosome: track.chromosome,
    start: track.start,
    end: track.end,
    data: track.data,
    color: track.color,
    y_min,
    y_max,
    label: track.label,
  }
}

///|
/// Get track data points.
pub fn KaryotypeTrack::get_data(track : KaryotypeTrack) -> Array[TrackPoint] {
  track.data
}

///|
/// Get number of data points.
pub fn KaryotypeTrack::get_n_points(track : KaryotypeTrack) -> Int {
  track.data.length()
}

///|
/// Filter a track to a specific chromosome.
pub fn KaryotypeTrack::filter_chromosome(track : KaryotypeTrack, chr : String) -> KaryotypeTrack {
  let filtered : Array[TrackPoint] = Array::new()
  let mut i = 0
  while i < track.data.length() {
    if track.data[i].chromosome == chr {
      filtered.push(track.data[i])
    }
    i = i + 1
  }
  KaryotypeTrack::{
    track_id: track.track_id,
    track_type: track.track_type,
    chromosome: chr,
    start: track.start,
    end: track.end,
    data: filtered,
    color: track.color,
    y_min: track.y_min,
    y_max: track.y_max,
    label: track.label,
  }
}

///|
/// Generate an ASCII representation of the karyotype plot.
pub fn KaryotypePlot::to_ascii(plot : KaryotypePlot, chr : String, width : Int) -> String {
  let w = if width < 40 { 60 } else { width }
  let chr_size = plot.get_chromosome_size(chr)
  if chr_size == 0 {
    return "Chromosome " + chr + " not found in genome " + plot.genome + "\n"
  }

  let mut result = "=== Karyotype Plot: " + plot.genome + " ===\n"
  result = result + "Chromosome: " + chr + " (size: " + chr_size.to_string() + " bp)\n"

  // Draw ideogram bar
  let mut bar = ""
  let mut i = 0
  while i < w {
    bar = bar + "="
    i = i + 1
  }
  result = result + bar + "\n"

  // Draw centromere marker
  let centromere_pos = ((chr_size.to_double() / 2.0) / chr_size.to_double() * (w.to_double() - 1.0)).to_int()
  let mut cent_bar = ""
  let mut j = 0
  while j < w {
    if j == centromere_pos {
      cent_bar = cent_bar + "*"
    } else {
      cent_bar = cent_bar + " "
    }
    j = j + 1
  }
  result = result + cent_bar + "\n"

  // Draw tick marks
  let mut ticks = ""
  let mut k = 0
  while k < w {
    if k == 0 || k == (w / 4) || k == (w / 2) || k == (3 * w / 4) || k == w - 1 {
      ticks = ticks + "|"
    } else {
      ticks = ticks + "-"
    }
    k = k + 1
  }
  result = result + ticks + "\n"

  // Position labels
  result = result + "0"
  let mut label_spaces = ""
  let mut ls = 0
  while ls < w - 10 {
    label_spaces = label_spaces + " "
    ls = ls + 1
  }
  let end_label = chr_size.to_string()
  result = result + label_spaces + end_label + "\n\n"

  // Draw track data
  let mut track_idx = 0
  while track_idx < plot.tracks.length() {
    let track = plot.tracks[track_idx]
    if track.chromosome == chr {
      result = result + "Track: " + track.label + " (" + track.track_type.to_string() + ")\n"
      let mut data_bar = ""
      let mut m = 0
      while m < w {
        data_bar = data_bar + "."
        m = m + 1
      }

      // Plot data points
      let mut p = 0
      while p < track.data.length() {
        let pt = track.data[p]
        let rel_pos = (pt.position.to_double() / chr_size.to_double() * (w.to_double() - 1.0)).to_int()
        let value = (pt.value - track.y_min) / (track.y_max - track.y_min)
        if rel_pos >= 0 && rel_pos < w {
          let ch = if value > 0.75 { "#" } else if value > 0.5 { "+" } else if value > 0.25 { "o" } else { "." }
          data_bar = data_bar.substring(start=0, end=rel_pos).to_string() + ch + data_bar.substring(start=rel_pos + 1, end=w).to_string()
        }
        p = p + 1
      }
      result = result + data_bar + "\n"
    }
    track_idx = track_idx + 1
  }

  result = result + "\n---\n"
  result
}

///|
/// Generate a summary string for the karyotype plot.
pub fn KaryotypePlot::summary(plot : KaryotypePlot) -> String {
  let mut s = "=== Karyotype Plot Summary ===\n".to_string()
  s = s + "Genome: " + plot.genome + "\n"
  s = s + "Chromosomes: " + plot.chromosomes.length().to_string() + "\n"
  s = s + "Tracks: " + plot.tracks.length().to_string() + "\n"
  s = s + "Regions: " + plot.regions.length().to_string() + "\n"

  let mut i = 0
  let total_bp = 0
  while i < plot.chromosomes.length() {
    let chr = plot.chromosomes[i]
    s = s + "  " + chr + ": " + plot.get_chromosome_size(chr).to_string() + " bp\n"
    i = i + 1
  }
  s
}

///|
/// Create a sample karyotype plot for human.
pub fn karyotype_sample() -> KaryotypePlot {
  let plot = KaryotypePlot::new("hg38")

  // Add some sample data tracks
  let mut gc_track = KaryotypeTrack::new("gc_content", TrackType::Lines, "chr1", 0, 248956422, "GC Content")
  gc_track = gc_track.set_color("#FF6600")

  let mut i = 0
  while i < 50 {
    let pos = i * 497913
    let value = 0.35 + 0.15 * ((i.to_double() / 50.0) * 2.0 - 1.0) * ((i.to_double() / 50.0) * 2.0 - 1.0)
    gc_track = gc_track.add_point(TrackPoint::new("chr1", pos, value, ""))
    i = i + 1
  }

  // Add a second track
  let mut snp_track = KaryotypeTrack::new("snp_density", TrackType::Bars, "chr1", 0, 248956422, "SNP Density")
  snp_track = snp_track.set_color("#003399")

  let mut j = 0
  while j < 30 {
    let pos = j * 829854
    let value = 0.2 + 0.6 * ((j.to_double() / 30.0))
    snp_track = snp_track.add_point(TrackPoint::new("chr1", pos, value, ""))
    j = j + 1
  }

  // Add regions of interest
  let region1 = KaryotypeRegion::new("chr1", 11000000, 12000000, "BRCA2 region", "#FF0000")
  let region2 = KaryotypeRegion::new("chr1", 28000000, 30000000, "TP53 region", "#00FF00")

  let mut result = plot
  result = result.add_track(gc_track)
  result = result.add_track(snp_track)
  result = result.add_region(region1)
  result = result.add_region(region2)

  result
}