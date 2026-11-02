///|
/// Visualization utilities for structure and sequence data.
///
/// Biotite integrates with matplotlib for visualization. This module
/// provides data preparation helpers that generate plot-ready output
/// (text-based representations, ASCII plots, and data for plotting).

///|
/// Generate an ASCII art secondary structure diagram.
///
/// H = helix, E = strand, C = coil, space = gap
pub fn ss_ascii_diagram(
  sequence : String,
  ss : String,
  width~ : Int32 = 60,
) -> String {
  let n = sequence.length()
  let mut result = ""
  let mut pos = 0
  while pos < n {
    let end = if pos + width < n { pos + width } else { n }
    // Sequence line
    let mut seq_line = ""
    for i in pos..<end {
      seq_line = seq_line + sequence.unsafe_get(i).unsafe_to_char().to_string()
    }
    // SS line
    let mut ss_line = ""
    for i in pos..<end {
      if i < ss.length() {
        ss_line = ss_line + ss.unsafe_get(i).unsafe_to_char().to_string()
      } else {
        ss_line = ss_line + " "
      }
    }
    // Position ruler
    let mut ruler = ""
    let ruler_start = pos + 1
    let ruler_end = end
    let ruler_str = ruler_start.to_string() + "-" + ruler_end.to_string()
    ruler = ruler_str
    result = result + ruler + "\n"
    result = result + seq_line + "\n"
    result = result + ss_line + "\n"
    result = result + "\n"
    pos = pos + width
  }
  result
}

///|
/// Generate a simple ASCII bar chart for numeric data.
///
/// Each value gets a bar of proportional length using '#' characters.
pub fn ascii_bar_chart(
  labels : Array[String],
  values : Array[Double],
  max_width~ : Int32 = 50,
) -> String {
  if values.length() == 0 {
    return ""
  }
  let mut max_val = 0.0
  let mut min_val = 0.0
  for v in values {
    if v > max_val {
      max_val = v
    }
    if v < min_val {
      min_val = v
    }
  }
  let range = if max_val > 0.0 { max_val } else { 1.0 }
  let mut result = ""
  for i in 0..<values.length() {
    let label = if i < labels.length() { labels[i] } else { "" }
    let normalized = values[i] / range
    let bar_len = (normalized * max_width.to_double()).floor().to_int()
    let mut bar = ""
    for _ in 0..<bar_len {
      bar = bar + "#"
    }
    // Pad label to 10 chars
    let mut padded_label = label
    while padded_label.length() < 10 {
      padded_label = padded_label + " "
    }
    result = result +
      padded_label +
      " |" +
      bar +
      " " +
      values[i].to_string() +
      "\n"
  }
  result
}

///|
/// Generate a heatmap representation of a distance matrix.
///
/// Uses characters: ' ' (0-20%), '.' (20-40%), ':' (40-60%),
/// '!' (60-80%), '#' (80-100%).
pub fn distance_matrix_heatmap(
  matrix : Array[Array[Double]],
  max_distance~ : Double = 20.0,
) -> String {
  if matrix.length() == 0 {
    return ""
  }
  let mut result = ""
  for row in matrix {
    let mut line = ""
    for val in row {
      let normalized = if val > max_distance { 1.0 } else { val / max_distance }
      let c = if normalized < 0.2 {
        ' '
      } else if normalized < 0.4 {
        '.'
      } else if normalized < 0.6 {
        ':'
      } else if normalized < 0.8 {
        '!'
      } else {
        '#'
      }
      line = line + c.to_string()
    }
    result = result + line + "\n"
  }
  result
}

///|
/// Generate a Ramachandran plot in ASCII format.
///
/// Plots phi/psi angles as a 2D grid of characters.
pub fn ramachandran_ascii(
  phi_psi_pairs : Array[(Double, Double)],
  width~ : Int32 = 72,
  height~ : Int32 = 24,
) -> String {
  // Grid: phi from -180 to 180, psi from -180 to 180
  let grid : Array[Array[Char]] = Array::make(height, Array::make(width, ' '))
  for pair in phi_psi_pairs {
    let (phi, psi) = pair
    // Map phi [-180, 180] to [0, width)
    let x = ((phi + 180.0) / 360.0 * width.to_double()).floor().to_int()
    // Map psi [-180, 180] to [height-1, 0] (inverted y)
    let y = height -
      1 -
      ((psi + 180.0) / 360.0 * height.to_double()).floor().to_int()
    if x >= 0 && x < width && y >= 0 && y < height {
      grid[y][x] = '*'
    }
  }
  let mut result = ""
  // Top border
  result = result + "    -180      0     180 (phi)\n"
  for row in grid {
    let mut line = ""
    for c in row {
      line = line + c.to_string()
    }
    result = result + line + "\n"
  }
  result
}

///|
/// Generate a sequence conservation plot from alignment.
///
/// Uses Shannon entropy at each position to compute conservation.
/// Returns a text-based plot with conservation scores.
pub fn conservation_plot(
  alignment : Array[String],
  width~ : Int32 = 60,
) -> String {
  if alignment.length() == 0 {
    return ""
  }
  let n = alignment[0].length()
  let n_seqs = alignment.length().to_double()
  let mut result = ""
  let mut pos = 0
  while pos < n {
    let end = if pos + width < n { pos + width } else { n }
    // Conservation scores
    let mut scores = ""
    for i in pos..<end {
      let counts : Map[UInt16, Int32] = Map([], capacity=32)
      for seq in alignment {
        if i < seq.length() {
          let c = seq.unsafe_get(i)
          counts[c] = counts.get(c).unwrap_or(0) + 1
        }
      }
      // Max frequency / total = conservation
      let mut max_count = 0
      let keys = counts.keys().collect()
      for k in keys {
        if counts[k] > max_count {
          max_count = counts[k]
        }
      }
      let conservation = max_count.to_double() / n_seqs
      let c = if conservation > 0.9 {
        '*'
      } else if conservation > 0.7 {
        ':'
      } else if conservation > 0.5 {
        '.'
      } else {
        ' '
      }
      scores = scores + c.to_string()
    }
    result = result + scores + "\n"
    pos = pos + width
  }
  result
}
