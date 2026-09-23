///|
/// KEGG pathway visualization (pathview-style).
///
/// Inspired by the R/Bioconductor package pathview (Luo & Brouwer 2013), this
/// module maps gene-level log2 fold changes onto KEGG pathway nodes and
/// renders a color-coded SVG pathway diagram (green = down, red = up).
///
/// The pathway topology (nodes, edges, gene-to-node mapping) is provided by
/// the caller for self-containment; no network access is required.
///
/// Reference: Luo W, Brouwer C (2013) Bioinformatics 29(14):1830-1.

///|
/// An RGB color.
pub(all) struct PathviewColor {
  red : Int
  green : Int
  blue : Int
}

///|
/// A pathway node (gene box, compound, or group).
pub(all) struct PathviewNode {
  node_id : String
  label : String
  node_type : String   // "gene", "compound", "map", "group"
  x : Double
  y : Double
  width : Double
  height : Double
}

///|
/// A pathway edge.
pub(all) struct PathviewEdge {
  source : String
  target : String
  edge_type : String
}

///|
/// A pathway (nodes + edges + gene-to-node mapping).
pub(all) struct PathviewPathway {
  pathway_id : String
  nodes : Array[PathviewNode]
  edges : Array[PathviewEdge]
  gene_to_node : Map[String, String]
}

///|
/// Per-node summary of mapped gene values.
pub(all) struct PathviewNodeStat {
  node_id : String
  mean_value : Double
  num_up : Int
  num_down : Int
  num_genes : Int
}

///|
/// Create an empty pathway.
pub fn pathview_new(pathway_id : String) -> PathviewPathway {
  PathviewPathway::{
    pathway_id: pathway_id,
    nodes: [],
    edges: [],
    gene_to_node: Map([], capacity=32),
  }
}

///|
/// Add a node to the pathway.
pub fn pathview_add_node(pw : PathviewPathway, node : PathviewNode) -> PathviewPathway {
  let nodes : Array[PathviewNode] = Array(capacity=0)
  for n in pw.nodes {
    nodes.push(n)
  }
  nodes.push(node)
  PathviewPathway::{ pathway_id: pw.pathway_id, nodes, edges: pw.edges, gene_to_node: pw.gene_to_node, }
}

///|
/// Add an edge to the pathway.
pub fn pathview_add_edge(pw : PathviewPathway, edge : PathviewEdge) -> PathviewPathway {
  let edges : Array[PathviewEdge] = Array(capacity=0)
  for e in pw.edges {
    edges.push(e)
  }
  edges.push(edge)
  PathviewPathway::{ pathway_id: pw.pathway_id, nodes: pw.nodes, edges, gene_to_node: pw.gene_to_node, }
}

///|
/// Map a gene to a pathway node.
pub fn pathview_map_gene(pw : PathviewPathway, gene : String, node_id : String) -> PathviewPathway {
  let g2n : Map[String, String] = pw.gene_to_node
  g2n.set(gene, node_id)
  PathviewPathway::{ pathway_id: pw.pathway_id, nodes: pw.nodes, edges: pw.edges, gene_to_node: g2n, }
}

///|
/// Find a node by ID.
pub fn pathview_find_node(pw : PathviewPathway, node_id : String) -> PathviewNode? {
  for n in pw.nodes {
    if n.node_id == node_id {
      return Some(n)
    }
  }
  None
}

///|
/// Map a continuous log2 fold change to an RGB color (green -> gray -> red).
/// Values are clamped to +-limit (default 3.0).
pub fn pathview_logfc_color(logfc : Double, limit? : Double = 3.0) -> PathviewColor {
  let lim = if limit > 0.0 { limit } else { 3.0 }
  let mut t = logfc / lim
  if t > 1.0 { t = 1.0 }
  if t < -1.0 { t = -1.0 }
  if t < 0.0 {
    // gray (128) -> green (0,255,0)
    let f = 1.0 + t  // ranges 0..1 as t goes -1..0
    let r = (128.0 * f).to_int()
    let g = (128.0 + 127.0 * (1.0 - f)).to_int()
    let b = (128.0 * f).to_int()
    PathviewColor::{ red: r, green: g, blue: b }
  } else {
    // gray (128) -> red (255,0,0)
    let f = t  // 0..1
    let r = (128.0 + 127.0 * f).to_int()
    let g = (128.0 * (1.0 - f)).to_int()
    let b = (128.0 * (1.0 - f)).to_int()
    PathviewColor::{ red: r, green: g, blue: b }
  }
}

///|
/// Convert a color to an SVG `rgb(r,g,b)` string.
pub fn pathview_color_rgb(c : PathviewColor) -> String {
  "rgb(" + c.red.to_string() + "," + c.green.to_string() + "," + c.blue.to_string() + ")"
}

///|
/// Convert a color to a hex `#rrggbb` string.
pub fn pathview_color_hex(c : PathviewColor) -> String {
  "#" + pathview_hex2(c.red) + pathview_hex2(c.green) + pathview_hex2(c.blue)
}

///|
/// Two-digit hex for an integer 0..255.
fn pathview_hex2(v : Int) -> String {
  let digits = "0123456789ABCDEF"
  let hi = v / 16
  let lo = v % 16
  substring(digits, hi, 1) + substring(digits, lo, 1)
}

///|
/// Compute per-node statistics from gene-level values.
///
/// Parameters:
///   pw: the pathway (provides gene -> node mapping)
///   gene_vals: array of (gene, value) pairs (value = log2 fold change)
///
/// Returns one PathviewNodeStat per node that has at least one mapped gene.
pub fn pathview_node_stats(
  pw : PathviewPathway,
  gene_vals : Array[(String, Double)],
) -> Array[PathviewNodeStat] {
  // aggregate values per node
  let node_sums : Map[String, Double] = Map([], capacity=16)
  let node_counts : Map[String, Int] = Map([], capacity=16)
  let node_up : Map[String, Int] = Map([], capacity=16)
  let node_down : Map[String, Int] = Map([], capacity=16)
  for entry in gene_vals {
    let g = entry.0
    let v = entry.1
    match pw.gene_to_node.get(g) {
      Some(node_id) => {
        node_sums.set(node_id, node_sums.get(node_id).unwrap_or(0.0) + v)
        node_counts.set(node_id, node_counts.get(node_id).unwrap_or(0) + 1)
        if v > 0.0 {
          node_up.set(node_id, node_up.get(node_id).unwrap_or(0) + 1)
        } else if v < 0.0 {
          node_down.set(node_id, node_down.get(node_id).unwrap_or(0) + 1)
        }
      }
      None => ()
    }
  }
  let results : Array[PathviewNodeStat] = Array(capacity=0)
  for node_id in node_sums.keys().collect() {
    let cnt = node_counts.get(node_id).unwrap_or(0)
    let mean = if cnt > 0 { node_sums.get(node_id).unwrap_or(0.0) / cnt.to_double() } else { 0.0 }
    results.push(PathviewNodeStat::{
      node_id: node_id,
      mean_value: mean,
      num_up: node_up.get(node_id).unwrap_or(0),
      num_down: node_down.get(node_id).unwrap_or(0),
      num_genes: cnt,
    })
  }
  results
}

///|
/// Render the pathway as an SVG string with gene values mapped to node colors.
///
/// Parameters:
///   pw: the pathway
///   gene_vals: array of (gene, value) pairs (log2 fold change)
///   limit: absolute value at which color saturates (default 3.0)
pub fn pathview_render_svg(
  pw : PathviewPathway,
  gene_vals : Array[(String, Double)],
  limit? : Double = 3.0,
) -> String {
  let stats = pathview_node_stats(pw, gene_vals)
  // map node_id -> mean value
  let node_mean : Map[String, Double] = Map([], capacity=16)
  for s in stats {
    node_mean.set(s.node_id, s.mean_value)
  }
  let mut minx = 1000.0
  let mut miny = 1000.0
  let mut maxx = -1000.0
  let mut maxy = -1000.0
  for n in pw.nodes {
    if n.x < minx { minx = n.x }
    if n.y < miny { miny = n.y }
    if n.x + n.width > maxx { maxx = n.x + n.width }
    if n.y + n.height > maxy { maxy = n.y + n.height }
  }
  if minx > maxx { minx = 0.0; maxx = 100.0 }
  if miny > maxy { miny = 0.0; maxy = 100.0 }
  let w = maxx - minx + 20.0
  let h = maxy - miny + 20.0
  let sb = StringBuilder()
  sb.write_string("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" + w.to_string() + "\" height=\"" + h.to_string() + "\" viewBox=\"" + minx.to_string() + " " + miny.to_string() + " " + w.to_string() + " " + h.to_string() + "\">\n")
  sb.write_string("  <rect x=\"" + minx.to_string() + "\" y=\"" + miny.to_string() + "\" width=\"" + w.to_string() + "\" height=\"" + h.to_string() + "\" fill=\"white\"/>\n")
  // Draw edges
  for e in pw.edges {
    let a = pathview_find_node(pw, e.source)
    let b = pathview_find_node(pw, e.target)
    match (a, b) {
      (Some(na), Some(nb)) => {
        sb.write_string("  <line x1=\"" + (na.x + na.width / 2.0).to_string() + "\" y1=\"" + (na.y + na.height / 2.0).to_string() + "\" x2=\"" + (nb.x + nb.width / 2.0).to_string() + "\" y2=\"" + (nb.y + nb.height / 2.0).to_string() + "\" stroke=\"gray\" stroke-width=\"1\"/>\n")
      }
      _ => ()
    }
  }
  // Draw nodes
  for n in pw.nodes {
    let fill = match node_mean.get(n.node_id) {
      Some(v) => pathview_color_rgb(pathview_logfc_color(v, limit=limit))
      None => "lightgray"
    }
    sb.write_string("  <rect x=\"" + n.x.to_string() + "\" y=\"" + n.y.to_string() + "\" width=\"" + n.width.to_string() + "\" height=\"" + n.height.to_string() + "\" fill=\"" + fill + "\" stroke=\"black\" stroke-width=\"1\"/>\n")
    sb.write_string("  <text x=\"" + (n.x + 2.0).to_string() + "\" y=\"" + (n.y + n.height / 2.0 + 4.0).to_string() + "\" font-size=\"10\">" + n.label + "</text>\n")
  }
  sb.write_string("</svg>\n")
  sb.to_string()
}

///|
/// Render a color scale legend (green -> gray -> red) as an SVG string.
pub fn pathview_legend(limit? : Double = 3.0) -> String {
  let sb = StringBuilder()
  let n = 9
  let step = 20.0
  sb.write_string("<svg width=\"" + ((n + 2).to_double() * step).to_string() + "\" height=\"40\">\n")
  let mut i = 0
  while i <= n {
    let t = -1.0 + 2.0 * i.to_double() / n.to_double()
    let val = t * (if limit > 0.0 { limit } else { 3.0 })
    let c = pathview_logfc_color(val, limit=limit)
    sb.write_string("  <rect x=\"" + (i.to_double() * step).to_string() + "\" y=\"5\" width=\"" + step.to_string() + "\" height=\"15\" fill=\"" + pathview_color_rgb(c) + "\"/>\n")
    i = i + 1
  }
  sb.write_string("</svg>\n")
  sb.to_string()
}

///|
/// Count up-regulated, down-regulated, and unchanged genes.
pub fn pathview_gene_tally(gene_vals : Array[(String, Double)]) -> (Int, Int) {
  let mut up = 0
  let mut down = 0
  for entry in gene_vals {
    let v = entry.1
    if v > 0.0 {
      up = up + 1
    } else if v < 0.0 {
      down = down + 1
    }
  }
  (up, down)
}
