///|
/// Bio.Graphics.Chromosome - Chromosome-level visualization with SVG output.
/// Inspired by Biopython's Bio.Graphics.Chromosome module.
///
/// Provides data structures and SVG generation for visualizing
/// chromosomes with gene features, markers, and regions.
/// Supports linear and circular chromosome layouts.

// ============================================================================
// Data Structures
// ============================================================================

///|
/// Orientation of a feature on the chromosome.
pub enum ChrOrientation {
  Forward
  Reverse
  None
} derive(Eq)

///|
/// Create a Forward orientation.
pub fn chr_orientation_forward() -> ChrOrientation { Forward }

///|
/// Create a Reverse orientation.
pub fn chr_orientation_reverse() -> ChrOrientation { Reverse }

///|
/// Create a None orientation.
pub fn chr_orientation_none() -> ChrOrientation { None }

///|
/// Type of a genomic feature.
pub enum ChrFeatureType {
  Gene
  Exon
  Intron
  CDS
  UTR
  Promoter
  Enhancer
  Repetitive
  Marker
  Breakpoint
} derive(Eq)

///|
/// Create a Gene feature type.
pub fn chr_feature_type_gene() -> ChrFeatureType { Gene }

///|
/// Create an Exon feature type.
pub fn chr_feature_type_exon() -> ChrFeatureType { Exon }

///|
/// Create an Intron feature type.
pub fn chr_feature_type_intron() -> ChrFeatureType { Intron }

///|
/// Create a CDS feature type.
pub fn chr_feature_type_cds() -> ChrFeatureType { CDS }

///|
/// Create a UTR feature type.
pub fn chr_feature_type_utr() -> ChrFeatureType { UTR }

///|
/// Create a Promoter feature type.
pub fn chr_feature_type_promoter() -> ChrFeatureType { Promoter }

///|
/// A genomic feature on a chromosome.
pub struct ChrFeature {
  name : String
  start : Double
  end : Double
  orientation : ChrOrientation
  feature_type : ChrFeatureType
  color : String
  label : String
  height : Double
}

///|
/// A region of interest or highlight.
pub struct ChrRegion {
  name : String
  start : Double
  end : Double
  color : String
  opacity : Double
  label : String
}

///|
/// A chromosome arm (p or q).
pub struct ChrArm {
  arm_type : String
  start : Double
  end : Double
  centromere : Double
}

///|
/// A chromosome with features.
pub struct Chromosome {
  name : String
  length : Double
  centromere_pos : Double
  is_circular : Bool
  features : Array[ChrFeature]
  regions : Array[ChrRegion]
  bands : Array[(String, Double, Double, String)]
}

///|
/// Layout configuration for chromosome drawing.
pub struct ChrLayout {
  width : Double
  height : Double
  margin_left : Double
  margin_right : Double
  margin_top : Double
  margin_bottom : Double
  chromosome_thickness : Double
  label_font_size : Double
  show_labels : Bool
  show_centromere : Bool
  show_bands : Bool
  color_scheme : String
}

///|
/// SVG element for chromosome visualization.
pub struct ChrSVGElement {
  element_type : String
  attributes : Map[String, String]
  content : String
  children : Array[ChrSVGElement]
}

///|
/// Complete chromosome diagram.
pub struct ChrDiagram {
  chromosomes : Array[Chromosome]
  layout : ChrLayout
  title : String
  svg_elements : Array[ChrSVGElement]
}

// ============================================================================
// Constructors
// ============================================================================

///|
/// Create a chromosome feature.
pub fn ChrFeature::new(
  name~ : String,
  start~ : Double,
  end~ : Double,
  orientation? : ChrOrientation = None,
  feature_type? : ChrFeatureType = ChrFeatureType::Gene,
  color? : String = "#45B7D1",
  label? : String = "",
  height? : Double = 0.3,
) -> ChrFeature {
  ChrFeature::{
    name,
    start,
    end,
    orientation,
    feature_type,
    color,
    label: if label.length() > 0 { label } else { name },
    height,
  }
}

///|
/// Create a chromosome region.
pub fn ChrRegion::new(
  name~ : String,
  start~ : Double,
  end~ : Double,
  color? : String = "#FFD700",
  opacity? : Double = 0.3,
  label? : String = "",
) -> ChrRegion {
  ChrRegion::{
    name,
    start,
    end,
    color,
    opacity,
    label: if label.length() > 0 { label } else { name },
  }
}

///|
/// Create a chromosome.
pub fn Chromosome::new(
  name~ : String,
  length~ : Double,
  centromere_pos? : Double = 0.0,
  is_circular? : Bool = false,
) -> Chromosome {
  Chromosome::{
    name,
    length,
    centromere_pos,
    is_circular,
    features: [],
    regions: [],
    bands: [],
  }
}

///|
/// Create layout configuration.
pub fn ChrLayout::new(
  width? : Double = 800.0,
  height? : Double = 600.0,
  margin_left? : Double = 60.0,
  margin_right? : Double = 60.0,
  margin_top? : Double = 40.0,
  margin_bottom? : Double = 40.0,
  chromosome_thickness? : Double = 10.0,
  label_font_size? : Double = 12.0,
  show_labels? : Bool = true,
  show_centromere? : Bool = true,
  show_bands? : Bool = false,
  color_scheme? : String = "default",
) -> ChrLayout {
  ChrLayout::{
    width,
    height,
    margin_left,
    margin_right,
    margin_top,
    margin_bottom,
    chromosome_thickness,
    label_font_size,
    show_labels,
    show_centromere,
    show_bands,
    color_scheme,
  }
}

///|
/// Create a chromosome diagram.
pub fn ChrDiagram::new(
  chromosomes~ : Array[Chromosome],
  layout? : ChrLayout = ChrLayout::new(),
  title? : String = "",
) -> ChrDiagram {
  ChrDiagram::{
    chromosomes,
    layout,
    title,
    svg_elements: [],
  }
}

// ============================================================================
// Feature operations
// ============================================================================

///|
/// Add a feature to a chromosome.
pub fn Chromosome::add_feature(self : Chromosome, feature : ChrFeature) -> Chromosome {
  Chromosome::{
    name: self.name,
    length: self.length,
    centromere_pos: self.centromere_pos,
    is_circular: self.is_circular,
    features: self.features + [feature],
    regions: self.regions,
    bands: self.bands,
  }
}

///|
/// Add a region to a chromosome.
pub fn Chromosome::add_region(self : Chromosome, region : ChrRegion) -> Chromosome {
  Chromosome::{
    name: self.name,
    length: self.length,
    centromere_pos: self.centromere_pos,
    is_circular: self.is_circular,
    features: self.features,
    regions: self.regions + [region],
    bands: self.bands,
  }
}

///|
/// Add a band (stain pattern) to a chromosome.
pub fn Chromosome::add_band(
  self : Chromosome,
  name : String,
  start : Double,
  end : Double,
  stain : String,
) -> Chromosome {
  Chromosome::{
    name: self.name,
    length: self.length,
    centromere_pos: self.centromere_pos,
    is_circular: self.is_circular,
    features: self.features,
    regions: self.regions,
    bands: self.bands + [(name, start, end, stain)],
  }
}

///|
/// Get all features of a specific type.
pub fn Chromosome::features_by_type(
  self : Chromosome,
  feature_type : ChrFeatureType,
) -> Array[ChrFeature] {
  let result : Array[ChrFeature] = []
  for f in self.features {
    if f.feature_type == feature_type {
      result.push(f)
    }
  }
  result
}

///|
/// Get features overlapping a region.
pub fn Chromosome::features_in_region(
  self : Chromosome,
  start : Double,
  end : Double,
) -> Array[ChrFeature] {
  let result : Array[ChrFeature] = []
  for f in self.features {
    if f.start <= end && f.end >= start {
      result.push(f)
    }
  }
  result
}

///|
/// Get the GC band color based on stain intensity.
pub fn chr_band_color(stain : String) -> String {
  match stain {
    "gneg" => "#FFFFFF"
    "gpos25" => "#F0F0F0"
    "gpos50" => "#D0D0D0"
    "gpos75" => "#A0A0A0"
    "gpos100" => "#606060"
    "acen" => "#FF4444"
    "gvar" => "#44AAFF"
    "stalk" => "#888888"
    _ => "#E8E8E8"
  }
}

// ============================================================================
// SVG generation
// ============================================================================

///|
/// Helper: max of two doubles.
fn chr_max(a : Double, b : Double) -> Double {
  if a > b { a } else { b }
}

///|
/// Helper: min of two doubles.
fn chr_min(a : Double, b : Double) -> Double {
  if a < b { a } else { b }
}

///|
/// Helper: abs of double.
fn chr_abs(x : Double) -> Double {
  if x < 0.0 { -x } else { x }
}

///|
/// Helper: PI constant.
fn chr_pi() -> Double { 3.14159265358979323846 }

///|
/// Helper: create SVG attributes map.
fn chr_svg_attrs(attrs : Array[(String, String)]) -> Map[String, String] {
  let map : Map[String, String] = Map([], capacity=attrs.length())
  let mut i = 0
  while i < attrs.length() {
    let entry = attrs[i]
    let k = entry.0
    let v = entry.1
    map[k] = v
    i = i + 1
  }
  map
}

///|
/// Helper: create SVG element with positional arguments for convenience.
fn chr_svg_elem(
  element_type : String,
  attributes : Map[String, String],
  content : String,
  children : Array[ChrSVGElement],
) -> ChrSVGElement {
  ChrSVGElement::{
    element_type,
    attributes,
    content,
    children,
  }
}

///|
/// Generate SVG for a single chromosome (linear).
pub fn chr_chromosome_to_svg(
  chr : Chromosome,
  layout : ChrLayout,
  x_offset : Double,
  y_offset : Double,
  scale : Double,
) -> Array[ChrSVGElement] {
  let elements : Array[ChrSVGElement] = []

  let thickness = layout.chromosome_thickness
  let chrom_length = chr.length * scale

  // Chromosome body
  let body_rect = chr_svg_elem("rect", chr_svg_attrs([
    ("x", (x_offset).to_string()),
    ("y", (y_offset - thickness / 2.0).to_string()),
    ("width", chrom_length.to_string()),
    ("height", thickness.to_string()),
    ("fill", "#E0E0E0"),
    ("stroke", "#333333"),
    ("stroke-width", "1"),
    ("rx", "3"),
    ("ry", "3"),
  ]), "", [])
  elements.push(body_rect)

  // Centromere
  if layout.show_centromere && chr.centromere_pos > 0.0 {
    let cen_x = x_offset + chr.centromere_pos * scale
    let cen_width = thickness * 1.5
    let centromere = chr_svg_elem("ellipse", chr_svg_attrs([
      ("cx", cen_x.to_string()),
      ("cy", y_offset.to_string()),
      ("rx", (cen_width / 2.0).to_string()),
      ("ry", (thickness / 2.0 + 2.0).to_string()),
      ("fill", "#CC0000"),
      ("stroke", "#990000"),
      ("stroke-width", "1"),
    ]), "", [])
    elements.push(centromere)
  }

  // Bands
  if layout.show_bands {
    for band in chr.bands {
      let (_, bstart, bend, stain) = band
      let bx = x_offset + bstart * scale
      let bw = (bend - bstart) * scale
      let color = chr_band_color(stain)
      let band_elem = chr_svg_elem("rect", chr_svg_attrs([
        ("x", bx.to_string()),
        ("y", (y_offset - thickness / 2.0).to_string()),
        ("width", bw.to_string()),
        ("height", thickness.to_string()),
        ("fill", color),
        ("stroke", "none"),
      ]), "", [])
      elements.push(band_elem)
    }
  }

  // Features
  for feat in chr.features {
    let fx = x_offset + feat.start * scale
    let fw = chr_max((feat.end - feat.start) * scale, 2.0)
    let fh = thickness * feat.height

    // Feature shape based on orientation
    let feat_color = feat.color
    let arrow_size = chr_min(fw, 5.0)

    if feat.orientation == Forward {
      // Arrow pointing right
      let poly_points = [
        (fx).to_string() + "," + (y_offset - fh / 2.0).to_string(),
        (fx + fw - arrow_size).to_string() + "," + (y_offset - fh / 2.0).to_string(),
        (fx + fw).to_string() + "," + (y_offset).to_string(),
        (fx + fw - arrow_size).to_string() + "," + (y_offset + fh / 2.0).to_string(),
        (fx).to_string() + "," + (y_offset + fh / 2.0).to_string(),
      ]
      let points_str = poly_points.join(" ")
      let feat_elem = chr_svg_elem("polygon", chr_svg_attrs([
        ("points", points_str),
        ("fill", feat_color),
        ("stroke", "#333333"),
        ("stroke-width", "0.5"),
        ("opacity", "0.8"),
      ]), "", [])
      elements.push(feat_elem)
    } else if feat.orientation == Reverse {
      // Arrow pointing left
      let poly_points = [
        (fx + fw).to_string() + "," + (y_offset - fh / 2.0).to_string(),
        (fx + arrow_size).to_string() + "," + (y_offset - fh / 2.0).to_string(),
        (fx).to_string() + "," + (y_offset).to_string(),
        (fx + arrow_size).to_string() + "," + (y_offset + fh / 2.0).to_string(),
        (fx + fw).to_string() + "," + (y_offset + fh / 2.0).to_string(),
      ]
      let points_str = poly_points.join(" ")
      let feat_elem = chr_svg_elem("polygon", chr_svg_attrs([
        ("points", points_str),
        ("fill", feat_color),
        ("stroke", "#333333"),
        ("stroke-width", "0.5"),
        ("opacity", "0.8"),
      ]), "", [])
      elements.push(feat_elem)
    } else {
      // No orientation - rectangle
      let feat_elem = chr_svg_elem("rect", chr_svg_attrs([
        ("x", fx.to_string()),
        ("y", (y_offset - fh / 2.0).to_string()),
        ("width", fw.to_string()),
        ("height", fh.to_string()),
        ("fill", feat_color),
        ("stroke", "#333333"),
        ("stroke-width", "0.5"),
        ("opacity", "0.8"),
      ]), "", [])
      elements.push(feat_elem)
    }

    // Label
    if layout.show_labels && feat.label.length() > 0 {
      let label_x = fx + fw / 2.0
      let label_y = y_offset - thickness / 2.0 - 4.0
      let text_elem = chr_svg_elem("text", chr_svg_attrs([
        ("x", label_x.to_string()),
        ("y", label_y.to_string()),
        ("font-size", layout.label_font_size.to_string()),
        ("text-anchor", "middle"),
        ("fill", "#333333"),
      ]), feat.label, [])
      elements.push(text_elem)
    }
  }

  // Regions (highlighted areas)
  for region in chr.regions {
    let rx = x_offset + region.start * scale
    let rw = (region.end - region.start) * scale
    let region_elem = chr_svg_elem("rect", chr_svg_attrs([
      ("x", rx.to_string()),
      ("y", (y_offset - thickness / 2.0 - 15.0).to_string()),
      ("width", rw.to_string()),
      ("height", (thickness + 30.0 + 15.0).to_string()),
      ("fill", region.color),
      ("opacity", region.opacity.to_string()),
      ("stroke", region.color),
      ("stroke-width", "2"),
      ("rx", "2"),
      ("ry", "2"),
    ]), "", [])
    elements.push(region_elem)

    if layout.show_labels && region.label.length() > 0 {
      let label_x = rx + rw / 2.0
      let label_y = y_offset - thickness / 2.0 - 8.0
      let text_elem = chr_svg_elem("text", chr_svg_attrs([
        ("x", label_x.to_string()),
        ("y", label_y.to_string()),
        ("font-size", (layout.label_font_size - 2.0).to_string()),
        ("text-anchor", "middle"),
        ("fill", "#666666"),
      ]), region.label, [])
      elements.push(text_elem)
    }
  }

  // Chromosome name label
  let name_x = x_offset - 10.0
  let name_y = y_offset + 4.0
  let name_elem = chr_svg_elem("text", chr_svg_attrs([
    ("x", name_x.to_string()),
    ("y", name_y.to_string()),
    ("font-size", (layout.label_font_size + 2.0).to_string()),
    ("text-anchor", "end"),
    ("font-weight", "bold"),
    ("fill", "#333333"),
  ]), chr.name, [])
  elements.push(name_elem)

  elements
}

///|
/// Generate SVG for a circular chromosome.
pub fn chr_circular_chromosome_to_svg(
  chr : Chromosome,
  layout : ChrLayout,
  cx : Double,
  cy : Double,
  radius : Double,
  scale : Double,
) -> Array[ChrSVGElement] {
  let elements : Array[ChrSVGElement] = []

  // Chromosome circle
  let _circumference = chr.length * scale
  let r = radius
  let chr_circle = chr_svg_elem("circle", chr_svg_attrs([
    ("cx", cx.to_string()),
    ("cy", cy.to_string()),
    ("r", r.to_string()),
    ("fill", "none"),
    ("stroke", "#E0E0E0"),
    ("stroke-width", (layout.chromosome_thickness).to_string()),
  ]), "", [])
  elements.push(chr_circle)

  // Features as arcs
  for feat in chr.features {
    let start_angle = (feat.start / chr.length) * 2.0 * chr_pi()
    let end_angle = (feat.end / chr.length) * 2.0 * chr_pi()

    let x1 = cx + r * @math.cos(start_angle)
    let y1 = cy + r * @math.sin(start_angle)
    let x2 = cx + r * @math.cos(end_angle)
    let y2 = cy + r * @math.sin(end_angle)

    let large_arc = if end_angle - start_angle > chr_pi() { 1 } else { 0 }
    let sweep = 1

    let d = "M " + x1.to_string() + " " + y1.to_string() +
      " A " + r.to_string() + " " + r.to_string() +
      " 0 " + large_arc.to_string() + " " + sweep.to_string() +
      " " + x2.to_string() + " " + y2.to_string()

    let arc_elem = chr_svg_elem("path", chr_svg_attrs([
      ("d", d),
      ("fill", "none"),
      ("stroke", feat.color),
      ("stroke-width", (layout.chromosome_thickness * feat.height).to_string()),
      ("opacity", "0.8"),
    ]), "", [])
    elements.push(arc_elem)
  }

  // Centromere
  if chr.centromere_pos > 0.0 {
    let cen_angle = (chr.centromere_pos / chr.length) * 2.0 * chr_pi()
    let cen_x = cx + r * @math.cos(cen_angle)
    let cen_y = cy + r * @math.sin(cen_angle)
    let cen_elem = chr_svg_elem("circle", chr_svg_attrs([
      ("cx", cen_x.to_string()),
      ("cy", cen_y.to_string()),
      ("r", (layout.chromosome_thickness / 2.0 + 3.0).to_string()),
      ("fill", "#CC0000"),
    ]), "", [])
    elements.push(cen_elem)
  }

  // Chromosome name at center
  if layout.show_labels {
    let name_elem = chr_svg_elem("text", chr_svg_attrs([
      ("x", cx.to_string()),
      ("y", cy.to_string()),
      ("font-size", (layout.label_font_size + 4.0).to_string()),
      ("text-anchor", "middle"),
      ("font-weight", "bold"),
      ("fill", "#333333"),
    ]), chr.name, [])
    elements.push(name_elem)
  }

  elements
}

// ============================================================================
// SVG assembly
// ============================================================================

///|
/// Convert SVG element to string.
pub fn ChrSVGElement::to_svg(self : ChrSVGElement) -> String {
  let mut attrs = ""
  for key in self.attributes.keys() {
    attrs = attrs + " " + key + "=\"" + self.attributes[key] + "\""
  }

  let mut result = "<" + self.element_type + attrs
  if self.element_type == "text" || self.element_type == "title" || self.element_type == "desc" {
    result = result + ">" + self.content + "</" + self.element_type + ">"
  } else if self.children.length() > 0 {
    result = result + ">"
    for child in self.children {
      result = result + child.to_svg()
    }
    result = result + "</" + self.element_type + ">"
  } else {
    result = result + "/>"
  }
  result
}

///|
/// Generate complete SVG string for a chromosome diagram.
pub fn ChrDiagram::to_svg(self : ChrDiagram) -> String {
  let layout = self.layout
  let n_chr = self.chromosomes.length()

  // Calculate scale
  let max_length = chr_max_length(self.chromosomes)
  let available_width = layout.width - layout.margin_left - layout.margin_right
  let scale = if max_length > 0.0 { available_width / max_length } else { 1.0 }

  // Calculate spacing
  let total_height_needed = n_chr.to_double() * (layout.chromosome_thickness + 80.0)
  let available_height = layout.height - layout.margin_top - layout.margin_bottom
  let y_spacing = if n_chr > 0 { available_height / (n_chr + 1).to_double() } else { 50.0 }

  // Build SVG elements
  let svg_children : Array[ChrSVGElement] = []

  // Title
  if self.title.length() > 0 {
    let title_elem = chr_svg_elem("text", chr_svg_attrs([
      ("x", (layout.width / 2.0).to_string()),
      ("y", (25.0).to_string()),
      ("font-size", "18"),
      ("text-anchor", "middle"),
      ("font-weight", "bold"),
      ("fill", "#333333"),
    ]), self.title, [])
    svg_children.push(title_elem)
  }

  // Chromosomes
  let mut i = 0
  while i < n_chr {
    let chr = self.chromosomes[i]
    let y_offset = layout.margin_top + (i.to_double() + 1.0) * y_spacing + 30.0
    let x_offset = layout.margin_left + 50.0

    if chr.is_circular {
      let radius = layout.chromosome_thickness * 8.0
      let cx = x_offset + radius + 20.0
      let cy = y_offset
      let chr_elements = chr_circular_chromosome_to_svg(chr, layout, cx, cy, radius, scale)
      for elem in chr_elements {
        svg_children.push(elem)
      }
    } else {
      let chr_elements = chr_chromosome_to_svg(chr, layout, x_offset, y_offset, scale)
      for elem in chr_elements {
        svg_children.push(elem)
      }
    }
    i = i + 1
  }

  // SVG root element
  let svg_open = "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" +
    layout.width.to_string() + "\" height=\"" +
    layout.height.to_string() + "\" viewBox=\"0 0 " +
    layout.width.to_string() + " " +
    layout.height.to_string() + "\">"

  let mut svg_content = ""
  for child in svg_children {
    svg_content = svg_content + child.to_svg()
  }

  svg_open + svg_content + "</svg>"
}

///|
/// Get the maximum chromosome length.
fn chr_max_length(chromosomes : Array[Chromosome]) -> Double {
  let mut max = 0.0
  for chr in chromosomes {
    if chr.length > max {
      max = chr.length
    }
  }
  max
}

// ============================================================================
// Helper: create common chromosome diagrams
// ============================================================================

///|
/// Create a human karyotype-like diagram with sample features.
pub fn chr_create_human_karyotype() -> ChrDiagram {
  let layout = ChrLayout::new(
    width=900.0,
    height=700.0,
    chromosome_thickness=12.0,
    show_labels=true,
    show_centromere=true,
  )

  let chromosomes : Array[Chromosome] = []

  // Add a sample chromosome 1
  let chr1 = Chromosome::new(name="chr1", length=248956422.0, centromere_pos=121535434.0)
  let chr1 = chr1.add_feature(ChrFeature::new(name="BRCA2", start=32300000.0, end=32400268.0, orientation=Forward, feature_type=ChrFeatureType::Gene, color="#FF6B6B", label="BRCA2"))
  let chr1 = chr1.add_feature(ChrFeature::new(name="DMD", start=31100000.0, end=33300000.0, orientation=Forward, feature_type=ChrFeatureType::Gene, color="#4ECDC4", label="DMD"))
  let chr1 = chr1.add_region(ChrRegion::new(name="BRCA1 region", start=43000000.0, end=43200000.0, color="#FFE66D", opacity=0.3, label="BRCA1 locus"))
  chromosomes.push(chr1)

  // Add chromosome 7
  let chr7 = Chromosome::new(name="chr7", length=159345973.0, centromere_pos=54495348.0)
  let chr7 = chr7.add_feature(ChrFeature::new(name="EGFR", start=55086714.0, end=55275019.0, orientation=Forward, feature_type=ChrFeatureType::Gene, color="#95E1D3", label="EGFR"))
  let chr7 = chr7.add_feature(ChrFeature::new(name="CFTR", start=117287160.0, end=117595905.0, orientation=Reverse, feature_type=ChrFeatureType::Gene, color="#F38181", label="CFTR"))
  chromosomes.push(chr7)

  // Add chromosome 11
  let chr11 = Chromosome::new(name="chr11", length=135086622.0, centromere_pos=51000000.0)
  let chr11 = chr11.add_feature(ChrFeature::new(name="INS", start=2159000.0, end=2164000.0, orientation=Forward, feature_type=ChrFeatureType::Gene, color="#AA96DA", label="INS"))
  let chr11 = chr11.add_feature(ChrFeature::new(name="MYCN", start=1590000.0, end=1595000.0, orientation=Forward, feature_type=ChrFeatureType::Gene, color="#FCBAD3", label="MYCN"))
  chromosomes.push(chr11)

  ChrDiagram::new(chromosomes=chromosomes, layout=layout, title="Human Chromosome Map")
}

///|
/// Create a bacterial circular chromosome diagram.
pub fn chr_create_bacterial_chromosome() -> ChrDiagram {
  let layout = ChrLayout::new(
    width=600.0,
    height=600.0,
    chromosome_thickness=8.0,
    show_labels=true,
  )

  let ecoli = Chromosome::new(name="E. coli K-12", length=4600000.0, centromere_pos=2300000.0, is_circular=true)
  let ecoli = ecoli.add_feature(ChrFeature::new(name="rrl", start=404000.0, end=405000.0, orientation=Forward, feature_type=ChrFeatureType::Gene, color="#FF6B6B", label="rrl"))
  let ecoli = ecoli.add_feature(ChrFeature::new(name="rpoB", start=430000.0, end=435000.0, orientation=Forward, feature_type=ChrFeatureType::Gene, color="#4ECDC4", label="rpoB"))
  let ecoli = ecoli.add_feature(ChrFeature::new(name="lacZ", start=3600000.0, end=3603000.0, orientation=Forward, feature_type=ChrFeatureType::Gene, color="#95E1D3", label="lacZ"))
  let ecoli = ecoli.add_feature(ChrFeature::new(name="trpA", start=4000000.0, end=4001000.0, orientation=Reverse, feature_type=ChrFeatureType::Gene, color="#F38181", label="trpA"))
  let ecoli = ecoli.add_region(ChrRegion::new(name="Pathogenicity Island", start=1000000.0, end=1200000.0, color="#FFE66D", opacity=0.2, label="PAI-1"))

  ChrDiagram::new(chromosomes=[ecoli], layout=layout, title="E. coli K-12 Circular Chromosome")
}

// ============================================================================
// Analysis functions
// ============================================================================

///|
/// Get total number of features across all chromosomes.
pub fn ChrDiagram::total_features(self : ChrDiagram) -> Int {
  let mut count = 0
  for chr in self.chromosomes {
    count = count + chr.features.length()
  }
  count
}

///|
/// Get features grouped by type across all chromosomes.
pub fn ChrDiagram::features_by_type(self : ChrDiagram) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=10)
  for chr in self.chromosomes {
    for feat in chr.features {
      let type_name = chr_feature_type_name(feat.feature_type)
      if counts.contains(type_name) {
        counts[type_name] = counts[type_name] + 1
      } else {
        counts[type_name] = 1
      }
    }
  }
  counts
}

///|
/// Get chromosome summary as string.
pub fn ChrDiagram::summary(self : ChrDiagram) -> String {
  let mut result = "Chromosome Diagram: " + self.title + "\n"
  result = result + "Chromosomes: " + self.chromosomes.length().to_string() + "\n"
  result = result + "Total features: " + self.total_features().to_string() + "\n"
  result = result + "Max length: " + chr_max_length(self.chromosomes).to_string() + " bp\n"
  result
}

///|
/// Convert feature type to string name.
fn chr_feature_type_name(ft : ChrFeatureType) -> String {
  match ft {
    Gene => "gene"
    Exon => "exon"
    Intron => "intron"
    CDS => "cds"
    UTR => "utr"
    Promoter => "promoter"
    Enhancer => "enhancer"
    Repetitive => "repetitive"
    Marker => "marker"
    Breakpoint => "breakpoint"
  }
}

///|
/// Count features of a specific type across all chromosomes.
pub fn ChrDiagram::count_by_type(self : ChrDiagram, ft : ChrFeatureType) -> Int {
  let mut count = 0
  for chr in self.chromosomes {
    for feat in chr.features {
      if feat.feature_type == ft {
        count = count + 1
      }
    }
  }
  count
}

///|
/// Get all feature names across chromosomes.
pub fn ChrDiagram::all_feature_names(self : ChrDiagram) -> Array[String] {
  let names : Array[String] = []
  for chr in self.chromosomes {
    for feat in chr.features {
      names.push(feat.name)
    }
  }
  names
}
