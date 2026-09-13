///|
/// KEGG pathway SVG renderer, ported from Biopython's
/// `Bio.Graphics.KGML_vis` module.
///
/// The reference module renders KGML pathways to PDF/PNG via ReportLab. This
/// port renders to an SVG string instead, keeping the same visual semantics:
/// gene/enzyme/ortholog entries are rectangles filled with their `bgcolor`,
/// compounds are circles, maps are dashed rectangles, relations are colored
/// connector lines (by relation type), and reactions are arrowed edges
/// between substrate and product compounds.

///|
/// Rendering options for a KGML canvas.
pub(all) struct KgmlCanvasOptions {
  show_labels : Bool
  show_relations : Bool
  show_reactions : Bool
  margin : Double
  relation_colors : Map[String, String]
}

///|
/// Default options: labels on, relations and reactions on, 10 px margin,
/// and the canonical KEGG relation-type color palette.
pub fn kgml_default_options() -> KgmlCanvasOptions {
  let colors : Map[String, String] = Map([])
  colors.set("ECrel", "#1f77b4")
  colors.set("PPrel", "#2ca02c")
  colors.set("GErel", "#9467bd")
  colors.set("PCrel", "#ff7f0e")
  colors.set("maplink", "#7f7f7f")
  KgmlCanvasOptions::{
    show_labels: true,
    show_relations: true,
    show_reactions: true,
    margin: 10.0,
    relation_colors: colors,
  }
}

///|
/// A KGML pathway canvas ready to render.
pub struct KgmlCanvas {
  pathway : KgmlPathway
  options : KgmlCanvasOptions
}

///|
pub fn KgmlCanvas::new(
  pathway : KgmlPathway,
  options~ : KgmlCanvasOptions = kgml_default_options(),
) -> KgmlCanvas {
  KgmlCanvas::{ pathway, options, }
}

///|
fn kgml_format_double(v : Double) -> String {
  let rounded = (v * 100.0).to_int().to_double() / 100.0
  rounded.to_string()
}

///|
/// Center point of an entry's graphics, in KGML coordinates.
fn kgml_entry_center(entry : KgmlEntry) -> (Double, Double)? {
  match entry.graphics() {
    Some(g) => Some((g.x() + g.width() / 2.0, g.y() + g.height() / 2.0))
    None => None
  }
}

///|
/// Bounding box of all entry graphics: `(min_x, min_y, max_x, max_y)`.
fn kgml_bbox(pathway : KgmlPathway) -> (Double, Double, Double, Double) {
  let mut min_x = 0.0
  let mut min_y = 0.0
  let mut max_x = 0.0
  let mut max_y = 0.0
  let mut seen = false
  for entry in pathway.entries() {
    match entry.graphics() {
      Some(g) => {
        let x1 = g.x()
        let y1 = g.y()
        let x2 = g.x() + g.width()
        let y2 = g.y() + g.height()
        if !seen {
          min_x = x1
          min_y = y1
          max_x = x2
          max_y = y2
          seen = true
        } else {
          if x1 < min_x {
            min_x = x1
          }
          if y1 < min_y {
            min_y = y1
          }
          if x2 > max_x {
            max_x = x2
          }
          if y2 > max_y {
            max_y = y2
          }
        }
      }
      None => ()
    }
  }
  if !seen {
    (0.0, 0.0, 1.0, 1.0)
  } else {
    (min_x, min_y, max_x, max_y)
  }
}

///|
/// Build a `compound-id -> (cx, cy)` map from all entries whose `name`
/// field carries space-separated KEGG identifiers.
fn kgml_name_centers(
  pathway : KgmlPathway,
) -> Map[String, (Double, Double)] {
  let centers : Map[String, (Double, Double)] = Map([])
  for entry in pathway.entries() {
    match kgml_entry_center(entry) {
      Some(center) =>
        for name_view in entry.name().split(" ") {
          let name = name_view.trim().to_owned()
          if name.length() > 0 && !(centers.contains(name)) {
            centers.set(name, center)
          }
        }
      None => ()
    }
  }
  centers
}

///|
fn kgml_write_text(
  svg : StringBuilder,
  x : Double,
  y : Double,
  text : String,
  size~ : Int = 9,
  color~ : String = "#000000",
) -> Unit {
  svg.write_string(
    "<text x=\"" +
    kgml_format_double(x) +
    "\" y=\"" +
    kgml_format_double(y) +
    "\" font-size=\"" +
    size.to_string() +
    "\" fill=\"" +
    color +
    "\" text-anchor=\"middle\">" +
    text +
    "</text>\n",
  )
}

///|
fn kgml_write_line(
  svg : StringBuilder,
  x1 : Double,
  y1 : Double,
  x2 : Double,
  y2 : Double,
  color~ : String = "#000000",
  width~ : Double = 1.0,
  dash~ : String = "",
) -> Unit {
  svg.write_string(
    "<line x1=\"" +
    kgml_format_double(x1) +
    "\" y1=\"" +
    kgml_format_double(y1) +
    "\" x2=\"" +
    kgml_format_double(x2) +
    "\" y2=\"" +
    kgml_format_double(y2) +
    "\" stroke=\"" +
    color +
    "\" stroke-width=\"" +
    kgml_format_double(width) +
    "\"" +
    (if dash.length() > 0 {
      " stroke-dasharray=\"" + dash + "\""
    } else {
      ""
    }) +
    "/>\n",
  )
}

///|
/// Arrowhead polygon at `(x, y)` pointing from `(from_x, from_y)`.
fn kgml_write_arrowhead(
  svg : StringBuilder,
  from_x : Double,
  from_y : Double,
  x : Double,
  y : Double,
  color~ : String = "#000000",
) -> Unit {
  let dx = x - from_x
  let dy = y - from_y
  let len = (dx * dx + dy * dy).sqrt()
  if len == 0.0 {
    return
  }
  let ux = dx / len
  let uy = dy / len
  let px = 0.0 - uy
  let py = ux
  let size = 6.0
  let tip_x = x
  let tip_y = y
  let base_x = x - ux * size
  let base_y = y - uy * size
  let left_x = base_x + px * size * 0.5
  let left_y = base_y + py * size * 0.5
  let right_x = base_x - px * size * 0.5
  let right_y = base_y - py * size * 0.5
  svg.write_string(
    "<polygon points=\"" +
    kgml_format_double(tip_x) +
    "," +
    kgml_format_double(tip_y) +
    " " +
    kgml_format_double(left_x) +
    "," +
    kgml_format_double(left_y) +
    " " +
    kgml_format_double(right_x) +
    "," +
    kgml_format_double(right_y) +
    "\" fill=\"" +
    color +
    "\"/>\n",
  )
}

///|
/// Render the pathway to an SVG document.
pub fn KgmlCanvas::to_svg(
  self : KgmlCanvas,
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
  svg.write_string(
    "<rect x=\"0\" y=\"0\" width=\"" +
    width.to_string() +
    "\" height=\"" +
    height.to_string() +
    "\" fill=\"white\"/>\n",
  )
  let (min_x, min_y, max_x, max_y) = kgml_bbox(self.pathway)
  let bbox_w = max_x - min_x
  let bbox_h = max_y - min_y
  let margin = self.options.margin
  let avail_w = width.to_double() - 2.0 * margin
  let avail_h = height.to_double() - 2.0 * margin
  let scale = if bbox_w == 0.0 || bbox_h == 0.0 {
    1.0
  } else if avail_w / bbox_w < avail_h / bbox_h {
    avail_w / bbox_w
  } else {
    avail_h / bbox_h
  }
  let tx = margin - min_x * scale
  let ty = margin - min_y * scale
  fn sx(gx : Double) -> Double {
    gx * scale + tx
  }
  fn sy(gy : Double) -> Double {
    gy * scale + ty
  }
  // Title
  if self.pathway.title().length() > 0 {
    kgml_write_text(
      svg,
      width.to_double() / 2.0,
      14.0,
      self.pathway.title(),
      size=13,
      color="#222",
    )
  }
  // Entries
  for entry in self.pathway.entries() {
    match entry.graphics() {
      Some(g) => {
        let x = sx(g.x())
        let y = sy(g.y())
        let w = g.width() * scale
        let h = g.height() * scale
        let etype = entry.etype()
        let gtype = g.gtype()
        let fill = g.bgcolor()
        let stroke = g.fgcolor()
        if gtype == "circle" || etype == "compound" {
          let r = if w < h { w / 2.0 } else { h / 2.0 }
          svg.write_string(
            "<circle cx=\"" +
            kgml_format_double(x + w / 2.0) +
            "\" cy=\"" +
            kgml_format_double(y + h / 2.0) +
            "\" r=\"" +
            kgml_format_double(r) +
            "\" fill=\"" +
            fill +
            "\" stroke=\"" +
            stroke +
            "\" stroke-width=\"1\"/>\n",
          )
        } else if etype == "map" {
          svg.write_string(
            "<rect x=\"" +
            kgml_format_double(x) +
            "\" y=\"" +
            kgml_format_double(y) +
            "\" width=\"" +
            kgml_format_double(w) +
            "\" height=\"" +
            kgml_format_double(h) +
            "\" fill=\"" +
            fill +
            "\" stroke=\"" +
            stroke +
            "\" stroke-width=\"1\" stroke-dasharray=\"4,2\"/>\n",
          )
        } else {
          svg.write_string(
            "<rect x=\"" +
            kgml_format_double(x) +
            "\" y=\"" +
            kgml_format_double(y) +
            "\" width=\"" +
            kgml_format_double(w) +
            "\" height=\"" +
            kgml_format_double(h) +
            "\" fill=\"" +
            fill +
            "\" stroke=\"" +
            stroke +
            "\" stroke-width=\"1\"/>\n",
          )
        }
        if self.options.show_labels && g.name().length() > 0 {
          kgml_write_text(svg, x + w / 2.0, y + h / 2.0 + 3.0, g.name())
        }
      }
      None => ()
    }
  }
  // Relations
  if self.options.show_relations {
    for relation in self.pathway.relations() {
      let e1 = self.pathway.get_entry(relation.entry1())
      let e2 = self.pathway.get_entry(relation.entry2())
      match (e1, e2) {
        (Some(entry1), Some(entry2)) =>
          match (kgml_entry_center(entry1), kgml_entry_center(entry2)) {
            (Some((x1, y1)), Some((x2, y2))) => {
              let color = self.options.relation_colors.get(
                relation.rtype(),
              ).unwrap_or("#888888")
              kgml_write_line(
                svg,
                sx(x1),
                sy(y1),
                sx(x2),
                sy(y2),
                color~,
                width=1.0,
              )
            }
            _ => ()
          }
        _ => ()
      }
    }
  }
  // Reactions
  if self.options.show_reactions {
    let centers = kgml_name_centers(self.pathway)
    for reaction in self.pathway.reactions() {
      let arrow_color = "#444444"
      let reversible = reaction.rtype() == "reversible"
      for substrate in reaction.substrates() {
        match centers.get(substrate) {
          Some((sx1, sy1)) =>
            for product in reaction.products() {
              match centers.get(product) {
                Some((px, py)) => {
                  kgml_write_line(
                    svg,
                    sx(sx1),
                    sy(sy1),
                    sx(px),
                    sy(py),
                    color=arrow_color,
                    width=1.2,
                  )
                  kgml_write_arrowhead(
                    svg,
                    sx(sx1),
                    sy(sy1),
                    sx(px),
                    sy(py),
                    color=arrow_color,
                  )
                  if reversible {
                    kgml_write_arrowhead(
                      svg,
                      sx(px),
                      sy(py),
                      sx(sx1),
                      sy(sy1),
                      color=arrow_color,
                    )
                  }
                }
                None => ()
              }
            }
          None => ()
        }
      }
    }
  }
  svg.write_string("</svg>\n")
  svg.to_string()
}

///|
/// Convenience: render a KGML pathway string directly to SVG.
///
/// Returns `None` when the KGML fails to parse.
pub fn kgml_to_svg(
  content : String,
  width~ : Int = 600,
  height~ : Int = 400,
) -> String? {
  match parse_kgml(content) {
    Some(pathway) =>
      Some(KgmlCanvas::new(pathway).to_svg(width~, height~))
    None => None
  }
}

///|
/// Count entries by type for a quick visual summary.
pub fn kgml_entry_summary(pathway : KgmlPathway) -> Array[(String, Int)] {
  let counts : Map[String, Int] = Map([])
  let order : Array[String] = Array::new()
  for entry in pathway.entries() {
    let t = entry.etype()
    if !(counts.contains(t)) {
      counts.set(t, 0)
      order.push(t)
    }
    counts.set(t, counts.get(t).unwrap_or(0) + 1)
  }
  let result : Array[(String, Int)] = Array::new()
  for t in order {
    result.push((t, counts.get(t).unwrap_or(0)))
  }
  result
}