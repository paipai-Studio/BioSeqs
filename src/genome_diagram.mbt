///|
/// GenomeDiagram - Genome diagram data model for visualization.
/// Inspired by Biopython's Bio.Graphics.GenomeDiagram module.
/// Provides data structures for representing genome diagrams
/// with tracks, features, and styling for SVG output.

///|
/// FeatureShape: the shape used to draw a feature in the diagram.
pub enum FeatureShape {
  Rectangle
  Arrow
  Diamond
  CrossedArrow
  Terminators
}

///|
pub fn FeatureShape::to_string(self : FeatureShape) -> String {
  match self {
    Rectangle => "rectangle"
    Arrow => "arrow"
    Diamond => "diamond"
    CrossedArrow => "crossed_arrow"
    Terminators => "terminators"
  }
}

///|
pub fn FeatureShape::from_string(s : String) -> FeatureShape {
  match s {
    "arrow" => Arrow
    "diamond" => Diamond
    "crossed_arrow" => CrossedArrow
    "terminators" => Terminators
    _ => Rectangle
  }
}

///|
/// Helper to create a Rectangle shape.
pub fn gd_rectangle_shape() -> FeatureShape { Rectangle }

///|
/// Helper to create an Arrow shape.
pub fn gd_arrow_shape() -> FeatureShape { Arrow }

///|
/// Helper to create a Diamond shape.
pub fn gd_diamond_shape() -> FeatureShape { Diamond }

///|
/// Helper to create a CrossedArrow shape.
pub fn gd_crossed_arrow_shape() -> FeatureShape { CrossedArrow }

///|
/// Helper to create a Terminators shape.
pub fn gd_terminators_shape() -> FeatureShape { Terminators }

///|
/// DiagramStyle: styling options for the diagram.
pub struct DiagramStyle {
  circular : Bool
  linear : Bool
  scale : Double
  border : Bool
  color_scheme : String
}

///|
pub fn DiagramStyle::new(
  circular? : Bool = false,
  linear? : Bool = true,
  scale? : Double = 1.0,
  border? : Bool = true,
  color_scheme? : String = "default",
) -> DiagramStyle {
  DiagramStyle::{ circular, linear, scale, border, color_scheme }
}

///|
pub fn DiagramStyle::circular(self : DiagramStyle) -> Bool { self.circular }

///|
pub fn DiagramStyle::linear(self : DiagramStyle) -> Bool { self.linear }

///|
pub fn DiagramStyle::scale(self : DiagramStyle) -> Double { self.scale }

///|
pub fn DiagramStyle::border(self : DiagramStyle) -> Bool { self.border }

///|
pub fn DiagramStyle::color_scheme(self : DiagramStyle) -> String { self.color_scheme }

///|
pub fn DiagramStyle::set_circular(self : DiagramStyle, v : Bool) -> DiagramStyle {
  DiagramStyle::{ ..self, circular: v, linear: !v }
}

///|
pub fn DiagramStyle::set_linear(self : DiagramStyle, v : Bool) -> DiagramStyle {
  DiagramStyle::{ ..self, linear: v, circular: !v }
}

///|
pub fn DiagramStyle::set_scale(self : DiagramStyle, v : Double) -> DiagramStyle {
  DiagramStyle::{ ..self, scale: v }
}

///|
pub fn DiagramStyle::set_border(self : DiagramStyle, v : Bool) -> DiagramStyle {
  DiagramStyle::{ ..self, border: v }
}

///|
pub fn DiagramStyle::set_color_scheme(self : DiagramStyle, v : String) -> DiagramStyle {
  DiagramStyle::{ ..self, color_scheme: v }
}

///|
/// DiagramFeature: a feature on the diagram (nucleotide level).
pub struct DiagramFeature {
  start : Int
  end : Int
  strand : String
  label : String
  color : String
  shape : FeatureShape
  depth : Int
}

///|
pub fn DiagramFeature::new(
  start : Int,
  end : Int,
  strand? : String = "+",
  label? : String = "",
  color? : String = "#45B7D1",
  shape? : FeatureShape = Rectangle,
  depth? : Int = 0,
) -> DiagramFeature {
  DiagramFeature::{ start, end, strand, label, color, shape, depth }
}

///|
pub fn DiagramFeature::start(self : DiagramFeature) -> Int { self.start }

///|
pub fn DiagramFeature::end(self : DiagramFeature) -> Int { self.end }

///|
pub fn DiagramFeature::strand(self : DiagramFeature) -> String { self.strand }

///|
pub fn DiagramFeature::label(self : DiagramFeature) -> String { self.label }

///|
pub fn DiagramFeature::color(self : DiagramFeature) -> String { self.color }

///|
pub fn DiagramFeature::shape(self : DiagramFeature) -> FeatureShape { self.shape }

///|
pub fn DiagramFeature::depth(self : DiagramFeature) -> Int { self.depth }

///|
pub fn DiagramFeature::with_color(self : DiagramFeature, color : String) -> DiagramFeature {
  DiagramFeature::{ ..self, color }
}

///|
pub fn DiagramFeature::with_label(self : DiagramFeature, label : String) -> DiagramFeature {
  DiagramFeature::{ ..self, label }
}

///|
pub fn DiagramFeature::with_depth(self : DiagramFeature, depth : Int) -> DiagramFeature {
  DiagramFeature::{ ..self, depth }
}

///|
pub fn DiagramFeature::overlaps(self : DiagramFeature, start : Int, end : Int) -> Bool {
  self.start < end && self.end > start
}

///|
pub fn DiagramFeature::length(self : DiagramFeature) -> Int {
  let diff = self.end - self.start
  if diff > 0 { diff } else { -diff }
}

///|
/// TrackFeature: a feature belonging to a track (amino acid level).
pub struct TrackFeature {
  start : Int
  end : Int
  strand : String
  label : String
  color : String
  shape : FeatureShape
}

///|
pub fn TrackFeature::new(
  start : Int,
  end : Int,
  strand? : String = "+",
  label? : String = "",
  color? : String = "#FF6B6B",
  shape? : FeatureShape = Rectangle,
) -> TrackFeature {
  TrackFeature::{ start, end, strand, label, color, shape }
}

///|
pub fn TrackFeature::start(self : TrackFeature) -> Int { self.start }

///|
pub fn TrackFeature::end(self : TrackFeature) -> Int { self.end }

///|
pub fn TrackFeature::strand(self : TrackFeature) -> String { self.strand }

///|
pub fn TrackFeature::label(self : TrackFeature) -> String { self.label }

///|
pub fn TrackFeature::color(self : TrackFeature) -> String { self.color }

///|
pub fn TrackFeature::shape(self : TrackFeature) -> FeatureShape { self.shape }

///|
pub fn TrackFeature::with_color(self : TrackFeature, color : String) -> TrackFeature {
  TrackFeature::{ ..self, color }
}

///|
pub fn TrackFeature::with_label(self : TrackFeature, label : String) -> TrackFeature {
  TrackFeature::{ ..self, label }
}

///|
pub fn TrackFeature::length(self : TrackFeature) -> Int {
  let diff = self.end - self.start
  if diff > 0 { diff } else { -diff }
}

///|
/// Track: a track on the diagram containing track features.
pub struct Track {
  name : String
  features : Array[TrackFeature]
  is_labeled : Bool
  color_scheme : String
}

///|
pub fn Track::new(
  name : String,
  is_labeled? : Bool = true,
  color_scheme? : String = "default",
) -> Track {
  Track::{ name, features: [], is_labeled, color_scheme }
}

///|
pub fn Track::name(self : Track) -> String { self.name }

///|
pub fn Track::features(self : Track) -> Array[TrackFeature] { self.features }

///|
pub fn Track::is_labeled(self : Track) -> Bool { self.is_labeled }

///|
pub fn Track::color_scheme(self : Track) -> String { self.color_scheme }

///|
pub fn Track::feature_count(self : Track) -> Int {
  self.features.length()
}

///|
pub fn Track::add_feature(self : Track, feature : TrackFeature) -> Track {
  let new_features : Array[TrackFeature] = []
  for f in self.features {
    new_features.push(f)
  }
  new_features.push(feature)
  Track::{ ..self, features: new_features }
}

///|
pub fn Track::set_is_labeled(self : Track, v : Bool) -> Track {
  Track::{ ..self, is_labeled: v }
}

///|
pub fn Track::set_color_scheme(self : Track, v : String) -> Track {
  Track::{ ..self, color_scheme: v }
}

///|
/// Diagram: the main genome diagram structure.
pub struct Diagram {
  tracks : Array[Track]
  features : Array[DiagramFeature]
  name : String
  start : Int
  end : Int
  style : DiagramStyle
}

///|
pub fn Diagram::new(
  name : String,
  start : Int,
  end : Int,
) -> Diagram {
  Diagram::{
    tracks: [],
    features: [],
    name,
    start,
    end,
    style: DiagramStyle::new(),
  }
}

///|
pub fn Diagram::name(self : Diagram) -> String { self.name }

///|
pub fn Diagram::start(self : Diagram) -> Int { self.start }

///|
pub fn Diagram::end(self : Diagram) -> Int { self.end }

///|
pub fn Diagram::style(self : Diagram) -> DiagramStyle { self.style }

///|
pub fn Diagram::tracks(self : Diagram) -> Array[Track] { self.tracks }

///|
pub fn Diagram::features(self : Diagram) -> Array[DiagramFeature] { self.features }

///|
pub fn Diagram::length(self : Diagram) -> Int {
  let diff = self.end - self.start
  if diff > 0 { diff } else { -diff }
}

///|
pub fn Diagram::track_count(self : Diagram) -> Int {
  self.tracks.length()
}

///|
pub fn Diagram::feature_count(self : Diagram) -> Int {
  let mut count = self.features.length()
  for track in self.tracks {
    count = count + track.feature_count()
  }
  count
}

///|
pub fn Diagram::with_style(self : Diagram, style : DiagramStyle) -> Diagram {
  Diagram::{ ..self, style }
}

///|
pub fn Diagram::add_track(self : Diagram, track_name : String) -> Diagram {
  let track = Track::new(track_name)
  let new_tracks : Array[Track] = []
  for t in self.tracks {
    new_tracks.push(t)
  }
  new_tracks.push(track)
  Diagram::{ ..self, tracks: new_tracks }
}

///|
pub fn Diagram::add_feature_to_diagram(
  self : Diagram,
  start : Int,
  end : Int,
  strand : String,
  label : String,
  shape : FeatureShape,
) -> Diagram {
  let feature = DiagramFeature::new(
    start, end, strand=strand, label=label, shape=shape,
  )
  let new_features : Array[DiagramFeature] = []
  for f in self.features {
    new_features.push(f)
  }
  new_features.push(feature)
  Diagram::{ ..self, features: new_features }
}

///|
pub fn Diagram::add_feature_to_track(
  self : Diagram,
  track_name : String,
  start : Int,
  end : Int,
  strand : String,
  label : String,
) -> Diagram {
  let feature = TrackFeature::new(
    start, end, strand=strand, label=label,
  )
  let new_tracks : Array[Track] = []
  for track in self.tracks {
    if track.name == track_name {
      new_tracks.push(track.add_feature(feature))
    } else {
      new_tracks.push(track)
    }
  }
  Diagram::{ ..self, tracks: new_tracks }
}

///|
pub fn Diagram::get_track_features(self : Diagram, track_name : String) -> Array[TrackFeature] {
  for track in self.tracks {
    if track.name == track_name {
      return track.features
    }
  }
  []
}

///|
pub fn Diagram::find_overlapping(
  self : Diagram,
  start : Int,
  end : Int,
) -> Array[DiagramFeature] {
  let result : Array[DiagramFeature] = []
  for feature in self.features {
    if feature.overlaps(start, end) {
      result.push(feature)
    }
  }
  result
}

///|
pub fn Diagram::label_features_auto(self : Diagram, min_size : Int) -> Diagram {
  let new_features : Array[DiagramFeature] = []
  for feature in self.features {
    let f = if feature.length() > min_size && feature.label() == "" {
      feature.with_label("feature_" + feature.start.to_string())
    } else {
      feature
    }
    new_features.push(f)
  }
  Diagram::{ ..self, features: new_features }
}

///|
pub fn Diagram::to_svg(self : Diagram, width : Int, height : Int) -> String {
  let len = self.length()
  if len == 0 {
    return ""
  }
  let w = width.to_double()
  let pad = 30.0
  let usable_w = w - 2.0 * pad
  let scale = usable_w / len.to_double()

  let svg = StringBuilder::new()
  svg.write_string(
    "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"",
  )
  svg.write_string(width.to_string())
  svg.write_string("\" height=\"")
  svg.write_string(height.to_string())
  svg.write_string("\" viewBox=\"0 0 ")
  svg.write_string(width.to_string())
  svg.write_string(" ")
  svg.write_string(height.to_string())
  svg.write_string("\">\n")

  if self.style.border {
    svg.write_string(
      "<rect x=\"1\" y=\"1\" width=\"",
    )
    svg.write_string((width - 2).to_string())
    svg.write_string("\" height=\"")
    svg.write_string((height - 2).to_string())
    svg.write_string(
      "\" fill=\"none\" stroke=\"#333\" stroke-width=\"2\"/>\n",
    )
  }

  let name_y = 18.0
  svg.write_string(
    "<text x=\"",
  )
  svg.write_string(pad.to_string())
  svg.write_string("\" y=\"")
  svg.write_string(name_y.to_string())
  svg.write_string(
    "\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#333\">",
  )
  svg.write_string(self.name)
  svg.write_string("</text>\n")

  let backbone_y = height.to_double() - pad
  svg.write_string(
    "<line x1=\"",
  )
  svg.write_string(pad.to_string())
  svg.write_string("\" y1=\"")
  svg.write_string(backbone_y.to_string())
  svg.write_string("\" x2=\"")
  svg.write_string((w - pad).to_string())
  svg.write_string("\" y2=\"")
  svg.write_string(backbone_y.to_string())
  svg.write_string(
    "\" stroke=\"#666\" stroke-width=\"2\"/>\n",
  )

  let tick_count = 5
  let tick_width = usable_w / tick_count.to_double()
  let mut ti = 0
  while ti <= tick_count {
    let x = pad + ti.to_double() * tick_width
    let pos = self.start + (ti.to_double() * len.to_double() / tick_count.to_double()).to_int()
    svg.write_string(
      "<line x1=\"",
    )
    svg.write_string(x.to_string())
    svg.write_string("\" y1=\"")
    svg.write_string((backbone_y - 5.0).to_string())
    svg.write_string("\" x2=\"")
    svg.write_string(x.to_string())
    svg.write_string("\" y2=\"")
    svg.write_string((backbone_y + 5.0).to_string())
    svg.write_string(
      "\" stroke=\"#666\" stroke-width=\"1\"/>\n",
    )
    svg.write_string(
      "<text x=\"",
    )
    svg.write_string(x.to_string())
    svg.write_string("\" y=\"")
    svg.write_string((backbone_y + 18.0).to_string())
    svg.write_string(
      "\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#666\" text-anchor=\"middle\">",
    )
    svg.write_string(pos.to_string())
    svg.write_string("</text>\n")
    ti = ti + 1
  }

  let track_height = 20.0
  let track_gap = 5.0
  let top_tracks_y = name_y + 30.0

  let mut track_idx = 0
  for track in self.tracks {
    let ty = top_tracks_y + track_idx.to_double() * (track_height + track_gap)
    svg.write_string(
      "<text x=\"",
    )
    svg.write_string(pad.to_string())
    svg.write_string("\" y=\"")
    svg.write_string((ty + 12.0).to_string())
    svg.write_string(
      "\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#555\">",
    )
    svg.write_string(track.name)
    svg.write_string("</text>\n")

    for feat in track.features {
      let fx = pad + (feat.start - self.start).to_double() * scale
      let mut fw = (feat.end - feat.start).to_double() * scale
      if fw < 1.0 { fw = 1.0 }
      let fy = ty + 2.0
      let fh = track_height - 4.0
      let color = feat.color
      match feat.shape {
        Arrow => {
          let tip_x = fx + fw
          svg.write_string(
            "<polygon points=\"",
          )
          svg.write_string(fx.to_string())
          svg.write_string(",")
          svg.write_string((fy).to_string())
          svg.write_string(" ")
          svg.write_string((fx + fw * 0.7).to_string())
          svg.write_string(",")
          svg.write_string((fy).to_string())
          svg.write_string(" ")
          svg.write_string(tip_x.to_string())
          svg.write_string(",")
          svg.write_string((fy + fh / 2.0).to_string())
          svg.write_string(" ")
          svg.write_string((fx + fw * 0.7).to_string())
          svg.write_string(",")
          svg.write_string((fy + fh).to_string())
          svg.write_string(" ")
          svg.write_string(fx.to_string())
          svg.write_string(",")
          svg.write_string((fy + fh).to_string())
          svg.write_string("\" fill=\"")
          svg.write_string(color)
          svg.write_string("\" stroke=\"#333\" stroke-width=\"0.5\"/>\n")
        }
        _ => {
          svg.write_string(
            "<rect x=\"",
          )
          svg.write_string(fx.to_string())
          svg.write_string("\" y=\"")
          svg.write_string(fy.to_string())
          svg.write_string("\" width=\"")
          svg.write_string(fw.to_string())
          svg.write_string("\" height=\"")
          svg.write_string(fh.to_string())
          svg.write_string("\" fill=\"")
          svg.write_string(color)
          svg.write_string("\" stroke=\"#333\" stroke-width=\"0.5\"/>\n")
        }
      }
      if feat.label.length() > 0 && track.is_labeled {
        let label_x = fx + fw / 2.0
        svg.write_string(
          "<text x=\"",
        )
        svg.write_string(label_x.to_string())
        svg.write_string("\" y=\"")
        svg.write_string((fy + fh / 2.0 + 4.0).to_string())
        svg.write_string(
          "\" font-family=\"sans-serif\" font-size=\"9\" fill=\"#333\" text-anchor=\"middle\">",
        )
        svg.write_string(feat.label)
        svg.write_string("</text>\n")
      }
    }
    track_idx = track_idx + 1
  }

  let feat_start_y = top_tracks_y + track_idx.to_double() * (track_height + track_gap) + 10.0
  let mut fi = 0
  for feat in self.features {
    let fx = pad + (feat.start - self.start).to_double() * scale
    let mut fw = (feat.end - feat.start).to_double() * scale
    if fw < 1.0 { fw = 1.0 }
    let fy = feat_start_y + (fi.to_double() * 18.0)
    let fh = 14.0
    let color = feat.color
    let depth_offset = feat.depth.to_double() * 3.0
    match feat.shape {
      Arrow => {
        let tip_x = fx + fw
        svg.write_string(
          "<polygon points=\"",
        )
        svg.write_string(fx.to_string())
        svg.write_string(",")
        svg.write_string((fy + depth_offset).to_string())
        svg.write_string(" ")
        svg.write_string((fx + fw * 0.7).to_string())
        svg.write_string(",")
        svg.write_string((fy + depth_offset).to_string())
        svg.write_string(" ")
        svg.write_string(tip_x.to_string())
        svg.write_string(",")
        svg.write_string((fy + fh / 2.0 + depth_offset).to_string())
        svg.write_string(" ")
        svg.write_string((fx + fw * 0.7).to_string())
        svg.write_string(",")
        svg.write_string((fy + fh + depth_offset).to_string())
        svg.write_string(" ")
        svg.write_string(fx.to_string())
        svg.write_string(",")
        svg.write_string((fy + fh + depth_offset).to_string())
        svg.write_string("\" fill=\"")
        svg.write_string(color)
        svg.write_string("\" stroke=\"#333\" stroke-width=\"0.5\"/>\n")
      }
      Diamond => {
        let cx = fx + fw / 2.0
        let cy = fy + fh / 2.0 + depth_offset
        svg.write_string(
          "<polygon points=\"",
        )
        svg.write_string(cx.to_string())
        svg.write_string(",")
        svg.write_string((cy - fh / 2.0).to_string())
        svg.write_string(" ")
        svg.write_string((fx + fw).to_string())
        svg.write_string(",")
        svg.write_string(cy.to_string())
        svg.write_string(" ")
        svg.write_string(cx.to_string())
        svg.write_string(",")
        svg.write_string((cy + fh / 2.0).to_string())
        svg.write_string(" ")
        svg.write_string(fx.to_string())
        svg.write_string(",")
        svg.write_string(cy.to_string())
        svg.write_string("\" fill=\"")
        svg.write_string(color)
        svg.write_string("\" stroke=\"#333\" stroke-width=\"0.5\"/>\n")
      }
      CrossedArrow => {
        let tip_x = fx + fw
        svg.write_string(
          "<line x1=\"",
        )
        svg.write_string(fx.to_string())
        svg.write_string("\" y1=\"")
        svg.write_string((fy + depth_offset).to_string())
        svg.write_string("\" x2=\"")
        svg.write_string(tip_x.to_string())
        svg.write_string("\" y2=\"")
        svg.write_string((fy + fh / 2.0 + depth_offset).to_string())
        svg.write_string(
          "\" stroke=\"",
        )
        svg.write_string(color)
        svg.write_string("\" stroke-width=\"2\"/>\n")
        svg.write_string(
          "<line x1=\"",
        )
        svg.write_string(fx.to_string())
        svg.write_string("\" y1=\"")
        svg.write_string((fy + fh + depth_offset).to_string())
        svg.write_string("\" x2=\"")
        svg.write_string(tip_x.to_string())
        svg.write_string("\" y2=\"")
        svg.write_string((fy + fh / 2.0 + depth_offset).to_string())
        svg.write_string(
          "\" stroke=\"",
        )
        svg.write_string(color)
        svg.write_string("\" stroke-width=\"2\"/>\n")
      }
      Terminators => {
        svg.write_string(
          "<line x1=\"",
        )
        svg.write_string(fx.to_string())
        svg.write_string("\" y1=\"")
        svg.write_string((fy + depth_offset).to_string())
        svg.write_string("\" x2=\"")
        svg.write_string(fx.to_string())
        svg.write_string("\" y2=\"")
        svg.write_string((fy + fh + depth_offset).to_string())
        svg.write_string(
          "\" stroke=\"",
        )
        svg.write_string(color)
        svg.write_string("\" stroke-width=\"3\"/>\n")
        svg.write_string(
          "<line x1=\"",
        )
        svg.write_string((fx + fw).to_string())
        svg.write_string("\" y1=\"")
        svg.write_string((fy + depth_offset).to_string())
        svg.write_string("\" x2=\"")
        svg.write_string((fx + fw).to_string())
        svg.write_string("\" y2=\"")
        svg.write_string((fy + fh + depth_offset).to_string())
        svg.write_string(
          "\" stroke=\"",
        )
        svg.write_string(color)
        svg.write_string("\" stroke-width=\"3\"/>\n")
      }
      Rectangle => {
        svg.write_string(
          "<rect x=\"",
        )
        svg.write_string(fx.to_string())
        svg.write_string("\" y=\"")
        svg.write_string((fy + depth_offset).to_string())
        svg.write_string("\" width=\"")
        svg.write_string(fw.to_string())
        svg.write_string("\" height=\"")
        svg.write_string(fh.to_string())
        svg.write_string("\" fill=\"")
        svg.write_string(color)
        svg.write_string("\" stroke=\"#333\" stroke-width=\"0.5\"/>\n")
      }
    }
    if feat.label.length() > 0 {
      let label_x = fx + fw / 2.0
      svg.write_string(
        "<text x=\"",
      )
      svg.write_string(label_x.to_string())
      svg.write_string("\" y=\"")
      svg.write_string((fy + fh + depth_offset + 12.0).to_string())
      svg.write_string(
        "\" font-family=\"sans-serif\" font-size=\"9\" fill=\"#333\" text-anchor=\"middle\">",
      )
      svg.write_string(feat.label)
      svg.write_string("</text>\n")
    }
    fi = fi + 1
  }

  svg.write_string("</svg>\n")
  svg.to_string()
}

///|
/// Create a new diagram with the given name and coordinate range.
pub fn gd_create_diagram(name : String, start : Int, end : Int) -> Diagram {
  Diagram::new(name, start, end)
}

///|
/// Add a new track to the diagram.
pub fn gd_add_track(diagram : Diagram, track_name : String) -> Diagram {
  diagram.add_track(track_name)
}

///|
/// Add a diagram-level feature to the diagram.
pub fn gd_add_feature(
  diagram : Diagram,
  start : Int,
  end : Int,
  strand : String,
  label : String,
  shape : FeatureShape,
) -> Diagram {
  diagram.add_feature_to_diagram(start, end, strand, label, shape)
}

///|
/// Add a feature to a named track in the diagram.
pub fn gd_add_track_feature(
  diagram : Diagram,
  track_name : String,
  start : Int,
  end : Int,
  strand : String,
  label : String,
) -> Diagram {
  diagram.add_feature_to_track(track_name, start, end, strand, label)
}

///|
/// Set the diagram's style.
pub fn gd_set_style(diagram : Diagram, style : DiagramStyle) -> Diagram {
  diagram.with_style(style)
}

///|
/// Get all track features for a named track.
pub fn gd_get_features(diagram : Diagram, track_name : String) -> Array[TrackFeature] {
  diagram.get_track_features(track_name)
}

///|
/// Find all diagram features overlapping the given range.
pub fn gd_find_overlapping_features(
  diagram : Diagram,
  start : Int,
  end : Int,
) -> Array[DiagramFeature] {
  diagram.find_overlapping(start, end)
}

///|
/// Get the length of the diagram (end - start).
pub fn gd_diagram_length(diagram : Diagram) -> Int {
  diagram.length()
}

///|
/// Get the number of tracks in the diagram.
pub fn gd_track_count(diagram : Diagram) -> Int {
  diagram.track_count()
}

///|
/// Get the total number of features in the diagram.
pub fn gd_feature_count(diagram : Diagram) -> Int {
  diagram.feature_count()
}

///|
/// Generate an SVG string representation of the diagram.
pub fn gd_to_svg_string(diagram : Diagram, width : Int, height : Int) -> String {
  diagram.to_svg(width, height)
}

///|
/// Set the color of a diagram feature.
pub fn gd_set_feature_color(feature : DiagramFeature, color : String) -> DiagramFeature {
  feature.with_color(color)
}

///|
/// Auto-label features larger than min_size that have no label.
pub fn gd_label_features(diagram : Diagram, min_size : Int) -> Diagram {
  diagram.label_features_auto(min_size)
}