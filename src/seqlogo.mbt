///|
/// seqLogo module - Sequence logo computation and rendering.
/// Inspired by the Bioconductor seqLogo package.
///
/// Computes position weight matrices (PWMs) and sequence logos where the
/// height of each letter represents its information content contribution.
/// Provides:
/// - PWM construction from aligned sequences
/// - Per-position information content (with custom backgrounds)
/// - Logo computation with stacked, colored letters
/// - ASCII-art and tabular rendering
/// - Consensus sequence extraction

// ============================================================================
// Private helpers
// ============================================================================

///|
/// Return the row index (0=A, 1=C, 2=G, 3=T) for a single-base string.
fn seqlogo_base_index(ch : String) -> Int {
  if ch == "A" || ch == "a" {
    0
  } else if ch == "C" || ch == "c" {
    1
  } else if ch == "G" || ch == "g" {
    2
  } else if ch == "T" || ch == "t" {
    3
  } else {
    -1
  }
}

///|
/// Return the canonical base character for a row index.
fn seqlogo_base_char(base : Int) -> Char {
  if base == 0 {
    'A'
  } else if base == 1 {
    'C'
  } else if base == 2 {
    'G'
  } else if base == 3 {
    'T'
  } else {
    'N'
  }
}

///|
/// Return the lowercase variant of an uppercase DNA base character.
fn seqlogo_to_lower(c : Char) -> Char {
  let code = c.to_int()
  if code >= 65 && code <= 90 {
    Int::unsafe_to_char(code + 32)
  } else {
    c
  }
}

///|
/// Sort (frequency, base) pairs by frequency descending (selection sort).
fn seqlogo_sort_desc(
  pairs : Array[(Double, Int)],
) -> Array[(Double, Int)] {
  let sorted = pairs.copy()
  let n = sorted.length()
  let mut i = 0
  while i < n - 1 {
    let mut max_idx = i
    let mut j = i + 1
    while j < n {
      if sorted[j].0 > sorted[max_idx].0 {
        max_idx = j
      }
      j = j + 1
    }
    if max_idx != i {
      let tmp = sorted[i]
      sorted[i] = sorted[max_idx]
      sorted[max_idx] = tmp
    }
    i = i + 1
  }
  sorted
}

// ============================================================================
// Data structures
// ============================================================================

///|
/// A position weight matrix for DNA (4 rows: A, C, G, T; L columns).
pub struct SeqLogoPwm {
  matrix : Array[Array[Double]]
  width : Int
  alphabet_size : Int
  alphabet : Array[String]
}

///|
/// Construct a PWM from a 4 x L matrix (row 0 = A, 1 = C, 2 = G, 3 = T).
pub fn SeqLogoPwm::new(matrix : Array[Array[Double]]) -> SeqLogoPwm {
  let alphabet_size = matrix.length()
  let width = if alphabet_size > 0 {
    matrix[0].length()
  } else {
    0
  }
  SeqLogoPwm::{
    matrix,
    width,
    alphabet_size,
    alphabet: ["A", "C", "G", "T"],
  }
}

///|
/// Build a PWM from a set of equal-length aligned DNA sequences.
pub fn SeqLogoPwm::from_sequences(
  sequences : Array[String],
) -> SeqLogoPwm {
  let n = sequences.length()
  if n == 0 {
    return SeqLogoPwm::new(Array::new())
  }
  let width = sequences[0].length()
  let counts : Array[Array[Double]] = Array::new()
  for _ in 0..<4 {
    counts.push(Array::make(width, 0.0))
  }
  for seq in sequences {
    for pos in 0..<width {
      if pos < seq.length() {
        let ch = seq[pos:pos + 1].to_owned()
        let base_idx = seqlogo_base_index(ch)
        if base_idx >= 0 {
          counts[base_idx][pos] = counts[base_idx][pos] + 1.0
        }
      }
    }
  }
  let matrix : Array[Array[Double]] = Array::new()
  for _ in 0..<4 {
    matrix.push(Array::make(width, 0.0))
  }
  let nd = n.to_double()
  for base in 0..<4 {
    for pos in 0..<width {
      matrix[base][pos] = counts[base][pos] / nd
    }
  }
  SeqLogoPwm::new(matrix)
}

///|
/// Number of positions in the PWM.
pub fn SeqLogoPwm::width(self : SeqLogoPwm) -> Int {
  self.width
}

///|
/// Alphabet size (number of rows).
pub fn SeqLogoPwm::alphabet_size(self : SeqLogoPwm) -> Int {
  self.alphabet_size
}

///|
/// Alphabet labels.
pub fn SeqLogoPwm::alphabet(self : SeqLogoPwm) -> Array[String] {
  self.alphabet
}

///|
/// Frequency of `base` (0-3) at `pos`.
pub fn SeqLogoPwm::get(self : SeqLogoPwm, base : Int, pos : Int) -> Double {
  if base >= 0 &&
    base < self.alphabet_size &&
    pos >= 0 &&
    pos < self.width {
    self.matrix[base][pos]
  } else {
    0.0
  }
}

///|
/// Human-readable representation of the PWM.
pub fn SeqLogoPwm::to_string(self : SeqLogoPwm) -> String {
  let buf = StringBuilder::new()
  buf.write_string("SeqLogoPwm(width=")
  buf.write_string(self.width.to_string())
  buf.write_string(")\n")
  for base in 0..<self.alphabet_size {
    buf.write_string(self.alphabet[base])
    buf.write_string(": ")
    for pos in 0..<self.width {
      if pos > 0 {
        buf.write_string("\t")
      }
      buf.write_string(self.matrix[base][pos].to_string())
    }
    buf.write_string("\n")
  }
  buf.to_string()
}

///|
/// A single rendered letter with its height and color.
pub struct SeqLogoLetter {
  letter : Char
  height : Double
  color : String
}

///|
/// Construct a SeqLogoLetter.
pub fn SeqLogoLetter::new(
  letter : Char,
  height : Double,
  color : String,
) -> SeqLogoLetter {
  SeqLogoLetter::{ letter, height, color }
}

///|
/// The base character.
pub fn SeqLogoLetter::letter(self : SeqLogoLetter) -> Char {
  self.letter
}

///|
/// Stacked height (in bits) of the letter.
pub fn SeqLogoLetter::height(self : SeqLogoLetter) -> Double {
  self.height
}

///|
/// Hex color string for the letter.
pub fn SeqLogoLetter::color(self : SeqLogoLetter) -> String {
  self.color
}

///|
/// A single column of the sequence logo.
pub struct SeqLogoColumn {
  letters : Array[SeqLogoLetter]
  info_content : Double
  position : Int
}

///|
/// Construct a SeqLogoColumn (letters start empty; use add_letter).
pub fn SeqLogoColumn::new(
  position : Int,
  info_content : Double,
) -> SeqLogoColumn {
  SeqLogoColumn::{
    letters: Array::new(),
    info_content,
    position,
  }
}

///|
/// Append a letter to this column.
pub fn SeqLogoColumn::add_letter(
  self : SeqLogoColumn,
  letter : SeqLogoLetter,
) -> Unit {
  self.letters.push(letter)
}

///|
/// Letters in the column, ordered from highest to lowest.
pub fn SeqLogoColumn::letters(self : SeqLogoColumn) -> Array[SeqLogoLetter] {
  self.letters
}

///|
/// Information content (bits) for this column.
pub fn SeqLogoColumn::info_content(self : SeqLogoColumn) -> Double {
  self.info_content
}

///|
/// Position index of this column.
pub fn SeqLogoColumn::position(self : SeqLogoColumn) -> Int {
  self.position
}

///|
/// A complete sequence logo.
pub struct SeqLogoResult {
  columns : Array[SeqLogoColumn]
  total_ic : Double
  width : Int
  alphabet : Array[String]
}

///|
/// Construct a SeqLogoResult, computing total IC and width from the columns.
pub fn SeqLogoResult::new(
  columns : Array[SeqLogoColumn],
) -> SeqLogoResult {
  let mut total_ic = 0.0
  for col in columns {
    total_ic = total_ic + col.info_content()
  }
  SeqLogoResult::{
    columns,
    total_ic,
    width: columns.length(),
    alphabet: ["A", "C", "G", "T"],
  }
}

///|
/// Columns of the logo.
pub fn SeqLogoResult::columns(self : SeqLogoResult) -> Array[SeqLogoColumn] {
  self.columns
}

///|
/// Total information content across all columns.
pub fn SeqLogoResult::total_ic(self : SeqLogoResult) -> Double {
  self.total_ic
}

///|
/// Number of columns in the logo.
pub fn SeqLogoResult::width(self : SeqLogoResult) -> Int {
  self.width
}

///|
/// Alphabet used by the logo.
pub fn SeqLogoResult::alphabet(self : SeqLogoResult) -> Array[String] {
  self.alphabet
}

// ============================================================================
// Information content
// ============================================================================

///|
/// Maximum information content for an alphabet of `alphabet_size` symbols:
/// log2(alphabet_size). For DNA this is 2.0 bits.
pub fn seqlogo_max_information_content(alphabet_size : Int) -> Double {
  if alphabet_size <= 1 {
    return 0.0
  }
  @math.log2(alphabet_size.to_double())
}

///|
/// Per-position information content with a uniform background.
pub fn seqlogo_information_content(pwm : SeqLogoPwm) -> Array[Double] {
  let bg = [0.25, 0.25, 0.25, 0.25]
  seqlogo_information_content_bg(pwm, bg)
}

///|
/// Per-position information content with a custom background distribution.
/// IC(j) = sum_i pwm[i][j] * log2(pwm[i][j] / bg[i]).
pub fn seqlogo_information_content_bg(
  pwm : SeqLogoPwm,
  background : Array[Double],
) -> Array[Double] {
  let width = pwm.width()
  let ic = Array::make(width, 0.0)
  for pos in 0..<width {
    let mut s = 0.0
    for base in 0..<pwm.alphabet_size() {
      let f = pwm.get(base, pos)
      if f > 0.0 && base < background.length() && background[base] > 0.0 {
        s = s + f * @math.log2(f / background[base])
      }
    }
    ic[pos] = s
  }
  ic
}

///|
/// Total information content summed over all positions.
pub fn seqlogo_total_information_content(pwm : SeqLogoPwm) -> Double {
  let ic_arr = seqlogo_information_content(pwm)
  let mut total = 0.0
  for v in ic_arr {
    if v > 0.0 {
      total = total + v
    }
  }
  total
}

// ============================================================================
// Color scheme
// ============================================================================

///|
/// Hex color for a DNA base character.
/// A -> green, C -> blue, G -> orange, T -> red.
pub fn seqlogo_base_color(base : Char) -> String {
  if base == 'A' {
    "#00CC00"
  } else if base == 'C' {
    "#0000CC"
  } else if base == 'G' {
    "#FFB300"
  } else if base == 'T' {
    "#CC0000"
  } else {
    "#000000"
  }
}

// ============================================================================
// Logo computation
// ============================================================================

///|
/// Compute a sequence logo with a uniform background.
pub fn seqlogo_compute_logo(pwm : SeqLogoPwm) -> SeqLogoResult {
  let bg = [0.25, 0.25, 0.25, 0.25]
  seqlogo_compute_logo_bg(pwm, bg)
}

///|
/// Compute a sequence logo with a custom background distribution.
/// For each position the information content is computed, bases are sorted by
/// frequency (highest first), and each letter's height is IC * frequency.
pub fn seqlogo_compute_logo_bg(
  pwm : SeqLogoPwm,
  background : Array[Double],
) -> SeqLogoResult {
  let ic_arr = seqlogo_information_content_bg(pwm, background)
  let max_ic = seqlogo_max_information_content(pwm.alphabet_size())
  let columns : Array[SeqLogoColumn] = Array::new()
  for pos in 0..<pwm.width() {
    let raw_ic = ic_arr[pos]
    let ic = if raw_ic < 0.0 {
      0.0
    } else if max_ic > 0.0 && raw_ic > max_ic {
      max_ic
    } else {
      raw_ic
    }
    let pairs : Array[(Double, Int)] = Array::new()
    for base in 0..<pwm.alphabet_size() {
      pairs.push((pwm.get(base, pos), base))
    }
    let sorted = seqlogo_sort_desc(pairs)
    let col = SeqLogoColumn::new(pos, ic)
    for i in 0..<sorted.length() {
      let freq = sorted[i].0
      let base = sorted[i].1
      let base_char = seqlogo_base_char(base)
      let height = ic * freq
      let color = seqlogo_base_color(base_char)
      col.add_letter(SeqLogoLetter::new(base_char, height, color))
    }
    columns.push(col)
  }
  SeqLogoResult::new(columns)
}

// ============================================================================
// Visualization
// ============================================================================

///|
/// Render a sequence logo as ASCII art.
/// Letters are stacked per column with heights proportional to their
/// information contribution. Uppercase letters mark large contributions and
/// lowercase letters mark small ones. `max_height` is the row count used for
/// the tallest possible column (2 bits for DNA).
pub fn seqlogo_to_ascii(
  logo : SeqLogoResult,
  max_height : Int,
) -> String {
  if max_height <= 0 {
    return ""
  }
  let max_ic = seqlogo_max_information_content(4)
  let scale = if max_ic > 0.0 {
    max_height.to_double() / max_ic
  } else {
    0.0
  }
  // Build a per-column stack of characters from bottom (index 0) to top.
  let stacks : Array[Array[Char]] = Array::new()
  for col in logo.columns() {
    let stack : Array[Char] = Array::new()
    let letters = col.letters()
    // Iterate from lowest frequency (bottom) to highest (top).
    let mut i = letters.length() - 1
    while i >= 0 {
      let letter = letters[i]
      let rows = (letter.height() * scale + 0.5).to_int()
      let r = if rows < 0 { 0 } else { rows }
      let ch = if letter.height() >= 0.5 {
        letter.letter()
      } else {
        seqlogo_to_lower(letter.letter())
      }
      for _ in 0..<r {
        stack.push(ch)
      }
      i = i - 1
    }
    stacks.push(stack)
  }
  let width = logo.width()
  let buf = StringBuilder::new()
  let mut row = max_height - 1
  while row >= 0 {
    let mut col_idx = 0
    while col_idx < width {
      let stack = stacks[col_idx]
      if row < stack.length() {
        buf.write_char(stack[row])
      } else {
        buf.write_char(' ')
      }
      col_idx = col_idx + 1
    }
    buf.write_char('\n')
    row = row - 1
  }
  buf.to_string()
}

///|
/// Render a sequence logo as a tab-separated table.
/// Columns: position, base, height, IC.
pub fn seqlogo_to_text_table(logo : SeqLogoResult) -> String {
  let buf = StringBuilder::new()
  buf.write_string("position\tbase\theight\tIC\n")
  for col in logo.columns() {
    let pos_str = col.position().to_string()
    let ic_str = col.info_content().to_string()
    for letter in col.letters() {
      buf.write_string(pos_str)
      buf.write_string("\t")
      buf.write_char(letter.letter())
      buf.write_string("\t")
      buf.write_string(letter.height().to_string())
      buf.write_string("\t")
      buf.write_string(ic_str)
      buf.write_string("\n")
    }
  }
  buf.to_string()
}

// ============================================================================
// Utility
// ============================================================================

///|
/// Return the consensus sequence (most frequent base at each position).
pub fn seqlogo_consensus_sequence(pwm : SeqLogoPwm) -> String {
  let buf = StringBuilder::new()
  for pos in 0..<pwm.width() {
    let mut max_freq = -1.0
    let mut max_base = 0
    for base in 0..<pwm.alphabet_size() {
      let f = pwm.get(base, pos)
      if f > max_freq {
        max_freq = f
        max_base = base
      }
    }
    buf.write_char(seqlogo_base_char(max_base))
  }
  buf.to_string()
}

///|
/// Build a textual summary of a sequence logo.
pub fn seqlogo_logo_summary(logo : SeqLogoResult) -> String {
  let buf = StringBuilder::new()
  buf.write_string("Sequence Logo Summary\n")
  buf.write_string("======================\n")
  buf.write_string("Width: ")
  buf.write_string(logo.width().to_string())
  buf.write_string("\nTotal IC: ")
  buf.write_string(logo.total_ic().to_string())
  buf.write_string(" bits\nMean IC per position: ")
  let mean_ic = if logo.width() > 0 {
    logo.total_ic() / logo.width().to_double()
  } else {
    0.0
  }
  buf.write_string(mean_ic.to_string())
  buf.write_string(" bits\n")
  buf.write_string("Per-position IC:\n")
  for col in logo.columns() {
    buf.write_string("  pos ")
    buf.write_string(col.position().to_string())
    buf.write_string(": IC=")
    buf.write_string(col.info_content().to_string())
    let mut top_char = ' '
    let mut top_height = -1.0
    for letter in col.letters() {
      if letter.height() > top_height {
        top_height = letter.height()
        top_char = letter.letter()
      }
    }
    buf.write_string(" top=")
    buf.write_char(top_char)
    buf.write_string("\n")
  }
  buf.to_string()
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// A sample PWM of 10 positions with varying degrees of conservation.
pub fn seqlogo_sample_pwm() -> SeqLogoPwm {
  // 4 rows (A, C, G, T) x 10 columns.
  let matrix = [
    // A
    [
      0.90, 0.10, 0.25, 0.05, 0.20, 0.01, 0.40, 0.85, 0.15, 0.30,
    ],
    // C
    [
      0.03, 0.10, 0.25, 0.90, 0.30, 0.02, 0.40, 0.05, 0.15, 0.30,
    ],
    // G
    [
      0.04, 0.70, 0.25, 0.03, 0.30, 0.03, 0.10, 0.05, 0.55, 0.20,
    ],
    // T
    [
      0.03, 0.10, 0.25, 0.02, 0.20, 0.94, 0.10, 0.05, 0.15, 0.20,
    ],
  ]
  SeqLogoPwm::new(matrix)
}

///|
/// Twenty aligned DNA sequences of length 10 used to build a sample PWM.
pub fn seqlogo_sample_sequences() -> Array[String] {
  [
    "AGACATGAGC",
    "AGACATGAGG",
    "CGACGTGAGT",
    "AGTCATGAGC",
    "GGACATGAGC",
    "AGGCATGAGC",
    "AGACGTGATC",
    "AGACATGTGC",
    "AGACATGCGC",
    "TGACATGAGC",
    "AGACATGAGA",
    "AGACATGAGC",
    "CGACATGAGG",
    "AGATATGAGC",
    "AGACACGAGC",
    "AGAGATGAGC",
    "AGACATGAGG",
    "AGACATGACC",
    "AGACATGTGC",
    "AAACATGAGC",
  ]
}
