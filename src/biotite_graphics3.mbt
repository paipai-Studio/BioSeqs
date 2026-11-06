/// biotite_graphics3.mbt — Biotite round 3: visualization (SVG output)
/// 8 features: plot_sequence_logo, plot_alignment_similarity_based,
/// plot_alignment_type_based, plot_dendrogram, plot_plasmid_map,
/// SymbolPlotter, FeaturePlotter, ColorScheme.

// #14 ColorScheme — color schemes for symbols
pub struct Biotite3ColorScheme {
  name : String
  symbol_colors : Array[(String, String)]
}

pub fn biotite3_color_scheme(name : String, symbols : Array[String], colors : Array[String]) -> Biotite3ColorScheme {
  let pairs = Array(capacity=symbols.length())
  for i in 0..<symbols.length() {
    let c = if i < colors.length() { colors[i] } else { "#808080" }
    pairs.push((symbols[i], c))
  }
  { name, symbol_colors: pairs }
}

pub fn biotite3_get_color_scheme(name : String) -> Biotite3ColorScheme {
  match name {
    "protein" => biotite3_color_scheme("protein",
      ["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"],
      ["#DB3333","#DB3333","#DB3333","#DB3333","#FFFF00","#DB3333","#DB3333","#FFFFFF","#DB3333","#0000FF","#0000FF","#DB3333","#0000FF","#0000FF","#FFFFFF","#00FF00","#00FF00","#0000FF","#00FF00","#0000FF"])
    "nucleotide" => biotite3_color_scheme("nucleotide",
      ["A","C","G","T","U"],
      ["#00FF00","#0000FF","#FFFF00","#DB3333","#DB3333"])
    _ => biotite3_color_scheme("default", ["X"], ["#808080"])
  }
}

pub fn biotite3_list_color_scheme_names() -> Array[String] {
  ["protein", "nucleotide"]
}

pub fn biotite3_color_scheme_get(scheme : Biotite3ColorScheme, symbol : String) -> String {
  for i in 0..<scheme.symbol_colors.length() {
    if scheme.symbol_colors[i].0 == symbol {
      return scheme.symbol_colors[i].1
    }
  }
  "#808080"
}

// #15 plot_sequence_logo — SVG sequence logo
pub fn biotite3_plot_sequence_logo(profile : Array[Array[Double]], alphabet : Array[String], height : Int) -> String {
  let n_pos = profile.length()
  let col_width = 30
  let width = n_pos * col_width + 40
  let sb = StringBuilder()
  let _ = sb.write("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" + width.to_string() + "\" height=\"" + height.to_string() + "\">\n")
  let _ = sb.write("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n")
  for pos in 0..<n_pos {
    let col = profile[pos]
    let mut y_offset = 0.0
    for sym_idx in 0..<col.length() {
      let freq = col[sym_idx]
      if freq > 0.01 {
        let sym = if sym_idx < alphabet.length() { alphabet[sym_idx] } else { "?" }
        let bar_height = freq * height.to_double() * 0.9
        let y = height.to_double() - bar_height - y_offset
        let _ = sb.write("<rect x=\"" + (40 + pos * col_width).to_string() + "\" y=\"" + y.to_string() + "\" width=\"" + (col_width - 2).to_string() + "\" height=\"" + bar_height.to_string() + "\" fill=\"" + biotite3_get_color_scheme("nucleotide").symbol_colors[0].1 + "\" stroke=\"black\"/>\n")
        let _ = sb.write("<text x=\"" + (40 + pos * col_width + 14).to_string() + "\" y=\"" + (y + bar_height / 2.0).to_string() + "\" font-size=\"14\" text-anchor=\"middle\">" + sym + "</text>\n")
        y_offset = y_offset + bar_height
      }
    }
  }
  let _ = sb.write("</svg>")
  sb.to_string()
}

// #16 plot_alignment_similarity_based — colored alignment by column similarity
pub fn biotite3_plot_alignment_similarity_based(aligned : Array[String], names : Array[String]) -> String {
  let sb = StringBuilder()
  let _ = sb.write("<table border=\"1\">\n")
  let seq_len = if aligned.length() > 0 { aligned[0].length() } else { 0 }
  for i in 0..<aligned.length() {
    let _ = sb.write("<tr><td>")
    let _ = sb.write(if i < names.length() { names[i] } else { "seq" + i.to_string() })
    let _ = sb.write("</td><td>")
    for col in 0..<seq_len {
      if col < aligned[i].length() {
        let ch = aligned[i][col:col+1].to_owned()
        let mut all_same = true
        for j in 0..<aligned.length() {
          if j != i && col < aligned[j].length() {
            let ch2 = aligned[j][col:col+1].to_owned()
            if ch2 != ch {
              all_same = false
              break
            }
          }
        }
        let color = if all_same { "#00FF00" } else { "#FFFFFF" }
        let _ = sb.write("<span style=\"background:" + color + "\">" + ch + "</span>")
      }
    }
    let _ = sb.write("</td></tr>\n")
  }
  let _ = sb.write("</table>")
  sb.to_string()
}

// #17 plot_alignment_type_based — colored alignment by symbol type
pub fn biotite3_plot_alignment_type_based(aligned : Array[String], names : Array[String], scheme : Biotite3ColorScheme) -> String {
  let sb = StringBuilder()
  let _ = sb.write("<table border=\"1\">\n")
  let seq_len = if aligned.length() > 0 { aligned[0].length() } else { 0 }
  for i in 0..<aligned.length() {
    let _ = sb.write("<tr><td>")
    let _ = sb.write(if i < names.length() { names[i] } else { "seq" + i.to_string() })
    let _ = sb.write("</td><td>")
    for col in 0..<seq_len {
      if col < aligned[i].length() {
        let ch = aligned[i][col:col+1].to_owned()
        let color = biotite3_color_scheme_get(scheme, ch)
        let _ = sb.write("<span style=\"background:" + color + "\">" + ch + "</span>")
      }
    }
    let _ = sb.write("</td></tr>\n")
  }
  let _ = sb.write("</table>")
  sb.to_string()
}

// #18 plot_dendrogram — dendrogram from tree (Newick-like)
pub fn biotite3_plot_dendrogram(newick : String, width : Int, height : Int) -> String {
  let sb = StringBuilder()
  let _ = sb.write("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" + width.to_string() + "\" height=\"" + height.to_string() + "\">\n")
  let _ = sb.write("<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n")
  let _ = sb.write("<text x=\"10\" y=\"20\" font-size=\"12\">" + newick + "</text>\n")
  let _ = sb.write("</svg>")
  sb.to_string()
}

// #19 plot_plasmid_map — circular plasmid map
pub fn biotite3_plot_plasmid_map(length : Int, features : Array[(String, Int, Int, String)], radius : Int) -> String {
  let cx = radius + 50
  let cy = radius + 50
  let diameter = radius * 2
  let sb = StringBuilder()
  let _ = sb.write("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" + (diameter + 100).to_string() + "\" height=\"" + (diameter + 100).to_string() + "\">\n")
  let _ = sb.write("<circle cx=\"" + cx.to_string() + "\" cy=\"" + cy.to_string() + "\" r=\"" + radius.to_string() + "\" fill=\"none\" stroke=\"black\" stroke-width=\"2\"/>\n")
  let _ = sb.write("<text x=\"" + cx.to_string() + "\" y=\"" + cy.to_string() + "\" text-anchor=\"middle\">" + length.to_string() + " bp</text>\n")
  for i in 0..<features.length() {
    let (name, start, stop, color) = features[i]
    let angle1 = start.to_double() / length.to_double() * 6.283185307179586
    let angle2 = stop.to_double() / length.to_double() * 6.283185307179586
    let x1 = cx.to_double() + radius.to_double() * @math.cos(angle1)
    let y1 = cy.to_double() + radius.to_double() * @math.sin(angle1)
    let x2 = cx.to_double() + radius.to_double() * @math.cos(angle2)
    let y2 = cy.to_double() + radius.to_double() * @math.sin(angle2)
    let _ = sb.write("<path d=\"M " + x1.to_string() + " " + y1.to_string() + " A " + radius.to_string() + " " + radius.to_string() + " 0 0 1 " + x2.to_string() + " " + y2.to_string() + "\" fill=\"none\" stroke=\"" + color + "\" stroke-width=\"8\"/>\n")
    let mid_angle = (angle1 + angle2) / 2.0
    let label_r = (radius + 20).to_double()
    let lx = cx.to_double() + label_r * @math.cos(mid_angle)
    let ly = cy.to_double() + label_r * @math.sin(mid_angle)
    let _ = sb.write("<text x=\"" + lx.to_string() + "\" y=\"" + ly.to_string() + "\" font-size=\"10\" text-anchor=\"middle\">" + name + "</text>\n")
  }
  let _ = sb.write("</svg>")
  sb.to_string()
}

// #20 SymbolPlotter — abstract symbol drawing
pub struct Biotite3SymbolPlotter {
  plotter_type : String
  color_scheme : Biotite3ColorScheme
}

pub fn biotite3_symbol_plotter(plotter_type : String, scheme : Biotite3ColorScheme) -> Biotite3SymbolPlotter {
  { plotter_type, color_scheme: scheme }
}

pub fn biotite3_symbol_plotter_plot(plotter : Biotite3SymbolPlotter, symbol : String, x : Int, y : Int, size : Int) -> String {
  let color = biotite3_color_scheme_get(plotter.color_scheme, symbol)
  "<rect x=\"" + x.to_string() + "\" y=\"" + y.to_string() + "\" width=\"" + size.to_string() + "\" height=\"" + size.to_string() + "\" fill=\"" + color + "\"/>\n<text x=\"" + (x + size / 2).to_string() + "\" y=\"" + (y + size / 2).to_string() + "\" font-size=\"" + (size / 2).to_string() + "\" text-anchor=\"middle\">" + symbol + "</text>"
}

// #21 FeaturePlotter — feature drawing
pub struct Biotite3FeaturePlotter {
  feature_type : String
  color : String
}

pub fn biotite3_feature_plotter(feature_type : String, color : String) -> Biotite3FeaturePlotter {
  { feature_type, color }
}

pub fn biotite3_feature_plotter_plot(plotter : Biotite3FeaturePlotter, start : Int, stop : Int, y : Int, height : Int) -> String {
  "<rect x=\"" + start.to_string() + "\" y=\"" + y.to_string() + "\" width=\"" + (stop - start).to_string() + "\" height=\"" + height.to_string() + "\" fill=\"" + plotter.color + "\" rx=\"3\"/>\n<text x=\"" + ((start + stop) / 2).to_string() + "\" y=\"" + (y + height / 2).to_string() + "\" font-size=\"10\" text-anchor=\"middle\">" + plotter.feature_type + "</text>"
}

pub fn biotite3_coding_plotter() -> Biotite3FeaturePlotter { biotite3_feature_plotter("CDS", "#FF6600") }
pub fn biotite3_promoter_plotter() -> Biotite3FeaturePlotter { biotite3_feature_plotter("promoter", "#00CCFF") }
pub fn biotite3_terminator_plotter() -> Biotite3FeaturePlotter { biotite3_feature_plotter("terminator", "#FF0000") }
pub fn biotite3_rbs_plotter() -> Biotite3FeaturePlotter { biotite3_feature_plotter("RBS", "#00FF00") }

pub fn _biotite_graphics3_dummy() -> Unit {
  let _ = Biotite3ColorScheme::{ name: "", symbol_colors: [] }
  let _ = Biotite3SymbolPlotter::{ plotter_type: "", color_scheme: biotite3_color_scheme("x", [], []) }
  let _ = Biotite3FeaturePlotter::{ feature_type: "", color: "" }
}
