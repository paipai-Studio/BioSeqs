///|
/// Graphics - Bioinformatics graphics and visualization.
/// Inspired by Biopython Bio.Graphics.
/// Provides sequence logos, alignment visualization, and sequence feature plotting.

///|
pub struct SeqLogo {
  sequences : Array[String]
  height : Int
  width : Int
  colors : Map[String, String]
}

///|
pub fn SeqLogo::new(sequences : Array[String]) -> SeqLogo {
  SeqLogo::{
    sequences,
    height: 50,
    width: 400,
    colors: get_default_colors()
  }
}

///|
pub struct LogoColumn {
  position : Int
  letter_frequencies : Map[String, Double]
  information_content : Double
  letters : Array[(String, Double, Double)]
}

///|
pub fn LogoColumn::new(position : Int) -> LogoColumn {
  LogoColumn::{
    position,
    letter_frequencies: Map([], capacity=4),
    information_content: 0.0,
    letters: Array::new()
  }
}

///|
pub struct AlignmentPlot {
  sequences : Array[(String, String)]
  highlight_regions : Array[(Int, Int, String)]
  consensus : String
  show_consensus : Bool
}

///|
pub fn AlignmentPlot::new(sequences : Array[(String, String)]) -> AlignmentPlot {
  AlignmentPlot::{
    sequences,
    highlight_regions: Array::new(),
    consensus: "",
    show_consensus: true
  }
}

///|
pub struct FeaturePlot {
  features : Array[(String, Int, Int, String, String)]
  sequence_length : Int
  colors : Map[String, String]
}

///|
pub fn FeaturePlot::new(features : Array[(String, Int, Int, String, String)], sequence_length : Int) -> FeaturePlot {
  FeaturePlot::{
    features,
    sequence_length,
    colors: get_feature_colors()
  }
}

///|
pub fn get_default_colors() -> Map[String, String] {
  let colors = Map([], capacity=20)
  
  colors.set("A", "#FF0000")
  colors.set("T", "#0000FF")
  colors.set("C", "#00FF00")
  colors.set("G", "#FFA500")
  
  colors.set("N", "#808080")
  colors.set("-", "#CCCCCC")
  
  colors.set("R", "#FF0000")
  colors.set("K", "#FF0000")
  colors.set("D", "#0000FF")
  colors.set("E", "#0000FF")
  colors.set("Q", "#00FF00")
  colors.set("N", "#00FF00")
  colors.set("H", "#FFA500")
  colors.set("Y", "#FFA500")
  colors.set("S", "#FFFF00")
  colors.set("T", "#FFFF00")
  colors.set("P", "#FF00FF")
  colors.set("G", "#FF00FF")
  colors.set("A", "#808080")
  colors.set("V", "#808080")
  colors.set("L", "#808080")
  colors.set("I", "#808080")
  colors.set("M", "#808080")
  colors.set("F", "#808080")
  colors.set("W", "#808080")
  colors.set("C", "#808080")
  
  colors
}

///|
pub fn get_feature_colors() -> Map[String, String] {
  let colors = Map([], capacity=10)
  
  colors.set("exon", "#FF6B6B")
  colors.set("intron", "#4ECDC4")
  colors.set("CDS", "#45B7D1")
  colors.set("promoter", "#96CEB4")
  colors.set("UTR", "#FFEAA7")
  colors.set("start_codon", "#DDA0DD")
  colors.set("stop_codon", "#FF7F50")
  colors.set("repeat", "#98D8C8")
  colors.set("gene", "#F7DC6F")
  colors.set("misc", "#BB8FCE")
  
  colors
}

///|
pub fn seqlogo_calculate_columns(sequences : Array[String]) -> Array[LogoColumn] {
  let columns : Array[LogoColumn] = Array::new()
  
  if sequences.length() == 0 {
    return columns
  }
  
  let seq_len = sequences[0].length()
  let n_seqs = sequences.length()
  
  let mut pos = 0
  while pos < seq_len {
    let col = LogoColumn::new(pos + 1)
    let freqs = col.letter_frequencies.copy()
    
    let counts = Map([
      ("A", 0), ("T", 0), ("C", 0), ("G", 0), ("N", 0)
    ], capacity=5)
    
    let mut seq_idx = 0
    while seq_idx < n_seqs {
      if pos < sequences[seq_idx].length() {
        let base = char_at(sequences[seq_idx], pos)
        if counts.contains(base) {
          counts.set(base, counts.get(base).unwrap() + 1)
        } else {
          counts.set(base, 1)
        }
      }
      seq_idx = seq_idx + 1
    }
    
    let mut total = 0
    let count_keys = counts.keys().collect()
    let mut i = 0
    while i < count_keys.length() {
      total = total + counts.get(count_keys[i]).unwrap()
      i = i + 1
    }
    
    if total > 0 {
      let mut ic = 0.0
      let ic_max = @math.log2(n_seqs.to_double())
      
      let mut j = 0
      while j < count_keys.length() {
        let key = count_keys[j]
        let count = counts.get(key).unwrap()
        let freq = count.to_double() / total.to_double()
        freqs.set(key, freq)
        
        if freq > 0.0 {
          ic = ic - freq * @math.log2(freq)
        }
        j = j + 1
      }
      
      let info_content = ic / ic_max * ic_max
      
      let sorted_letters : Array[(String, Double)] = Array::new()
      let mut k = 0
      while k < count_keys.length() {
        let key = count_keys[k]
        if counts.get(key).unwrap() > 0 {
          sorted_letters.push((key, freqs.get(key).unwrap()))
        }
        k = k + 1
      }
      
      let mut l = 0
      while l < sorted_letters.length() - 1 {
        let mut m = l + 1
        while m < sorted_letters.length() {
          if sorted_letters[l].1 < sorted_letters[m].1 {
            let tmp = sorted_letters[l]
            sorted_letters[l] = sorted_letters[m]
            sorted_letters[m] = tmp
          }
          m = m + 1
        }
        l = l + 1
      }
      
      let letters : Array[(String, Double, Double)] = Array::new()
      let mut cum_height = 0.0
      let mut n = 0
      while n < sorted_letters.length() {
        let height = sorted_letters[n].1 * info_content
        letters.push((sorted_letters[n].0, cum_height, height))
        cum_height = cum_height + height
        n = n + 1
      }
      
      columns.push(LogoColumn::{
        position: pos + 1,
        letter_frequencies: freqs,
        information_content: info_content,
        letters: letters
      })
    }
    
    pos = pos + 1
  }
  
  columns
}

///|
pub fn seqlogo_generate_ascii(logo : SeqLogo) -> String {
  let columns = seqlogo_calculate_columns(logo.sequences)
  
  let mut output = ""
  
  let max_ic = @math.log2(logo.sequences.length().to_double())
  let scale = logo.height.to_double() / max_ic
  
  let mut row = logo.height - 1
  while row >= 0 {
    let mut line = ""
    let mut col_idx = 0
    while col_idx < columns.length() {
      let col = columns[col_idx]
      let mut found = false
      let mut letter_idx = 0
      while letter_idx < col.letters.length() {
        let (letter, start, height) = col.letters[letter_idx]
        let scaled_start = start * scale
        let scaled_height = height * scale
        
        if row.to_double() >= scaled_start && row.to_double() < scaled_start + scaled_height {
          line = line + letter
          found = true
          break
        }
        letter_idx = letter_idx + 1
      }
      if !found {
        line = line + " "
      }
      col_idx = col_idx + 1
    }
    if line.length() > 0 {
      output = output + line + "\n"
    }
    row = row - 1
  }
  
  output = output + "Position: "
  let mut pos = 0
  while pos < columns.length() {
    output = output + (pos + 1).to_string() + " "
    pos = pos + 1
  }
  output = output + "\n"
  
  output
}

///|
pub fn alignment_plot_generate_ascii(plot : AlignmentPlot) -> String {
  let mut output = ""
  
  if plot.show_consensus {
    let consensus = graphics_calculate_consensus(plot.sequences)
    output = output + "Consensus: " + consensus + "\n\n"
  }
  
  let mut seq_idx = 0
  while seq_idx < plot.sequences.length() {
    let (name, seq) = plot.sequences[seq_idx]
    
    let mut line = ""
    if name.length() > 12 {
      line = line + substring(name, 0, 12)
    } else {
      line = line + name
      let mut pad = 0
      while pad < 12 - name.length() {
        line = line + " "
        pad = pad + 1
      }
    }
    line = line + " " + seq + "\n"
    output = output + line
    
    seq_idx = seq_idx + 1
  }
  
  output
}

///|
fn graphics_calculate_consensus(sequences : Array[(String, String)]) -> String {
  if sequences.length() == 0 {
    return ""
  }
  
  let seq_len = sequences[0].1.length()
  let n_seqs = sequences.length()
  
  let mut consensus = ""
  
  let mut pos = 0
  while pos < seq_len {
    let counts = Map([
      ("A", 0), ("T", 0), ("C", 0), ("G", 0), ("-", 0)
    ], capacity=5)
    
    let mut seq_idx = 0
    while seq_idx < n_seqs {
      if pos < sequences[seq_idx].1.length() {
        let base = char_at(sequences[seq_idx].1, pos)
        if counts.contains(base) {
          counts.set(base, counts.get(base).unwrap() + 1)
        }
      }
      seq_idx = seq_idx + 1
    }
    
    let mut max_count = 0
    let mut max_base = "N"
    let count_keys = counts.keys().collect()
    let mut i = 0
    while i < count_keys.length() {
      let key = count_keys[i]
      let count = counts.get(key).unwrap()
      if count > max_count {
        max_count = count
        max_base = key
      }
      i = i + 1
    }
    
    consensus = consensus + max_base
    pos = pos + 1
  }
  
  consensus
}

///|
pub fn feature_plot_generate_ascii(plot : FeaturePlot) -> String {
  let mut output = ""
  
  output = output + "Sequence length: " + plot.sequence_length.to_string() + "\n\n"
  
  output = output + "Legend:\n"
  let color_keys = plot.colors.keys().collect()
  let mut i = 0
  while i < color_keys.length() {
    let key = color_keys[i]
    output = output + "  [" + key + "] " + plot.colors[key] + "\n"
    i = i + 1
  }
  output = output + "\n"
  
  output = output + "Feature track:\n"
  output = output + "|"
  let mut pos = 0
  while pos < plot.sequence_length {
    if pos % 10 == 0 {
      output = output + "-"
    } else {
      output = output + "."
    }
    pos = pos + 1
  }
  output = output + "|\n"
  
  let mut feat_idx = 0
  while feat_idx < plot.features.length() {
    let (name, start, end, feat_type, strand) = plot.features[feat_idx]
    
    let mut line = ""
    let mut pos = 0
    while pos < plot.sequence_length {
      if pos >= start && pos < end {
        line = line + "="
      } else {
        line = line + " "
      }
      pos = pos + 1
    }
    
    let strand_sym = if strand == "-" { "<-" } else { "->" }
    output = output + strand_sym + " [" + name + "] " + line + "\n"
    
    feat_idx = feat_idx + 1
  }
  
  output = output + "\nPositions: 0"
  pos = 0
  while pos <= plot.sequence_length {
    if pos % 20 == 0 && pos > 0 {
      output = output + " ".repeat(pos.to_string().length() - 1) + pos.to_string()
    }
    pos = pos + 1
  }
  output = output + "\n"
  
  output
}

///|
pub fn create_example_seqlogo() -> SeqLogo {
  let seqs = [
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA"
  ]
  SeqLogo::new(seqs)
}

///|
pub fn create_example_alignment_plot() -> AlignmentPlot {
  let seqs = [
    ("Human", "ATGGCCGGTAAATGA"),
    ("Mouse", "ATGGCCGGTAAATGA"),
    ("Rat", "ATGGCCGGTAATAGA"),
    ("Chicken", "ATGGCCGGTAATGGA")
  ]
  AlignmentPlot::new(seqs)
}

///|
pub fn create_example_feature_plot() -> FeaturePlot {
  let features = [
    ("GeneA", 0, 100, "gene", "+"),
    ("Exon1", 0, 30, "exon", "+"),
    ("Intron1", 30, 70, "intron", "+"),
    ("Exon2", 70, 100, "exon", "+"),
    ("Promoter", -20, 0, "promoter", "+"),
    ("GeneB", 150, 250, "gene", "-"),
    ("CDS", 160, 240, "CDS", "-")
  ]
  FeaturePlot::new(features, 300)
}