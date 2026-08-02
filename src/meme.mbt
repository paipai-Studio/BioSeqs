// Bio.motifs.meme - MEME motif format parsing.
// Inspired by Biopython's Bio.motifs.meme module.
//
// Parses MEME Minimal Motif Format files containing position-specific
// probability matrices (PSPMs) for discovered motifs. Each motif includes
// a letter-probability matrix, motif name, E-value, and background frequencies.
//
// Reference:
//   Bailey TL et al. (2009) MEME Suite: tools for motif discovery and searching.
//   Nucleic Acids Research 37:W202-W208.
//
// Minimal MEME Motif Format:
//   MEME version 4
//   ALPHABET= ACGT
//   strands: + -
//   Background letter frequencies
//   A 0.25 C 0.25 G 0.25 T 0.25
//   MOTIF motif_id alt_id
//   letter-probability matrix: alength= 4 w= 8 nsites= 20 E= 1.1e-03
//     0.250000  0.250000  0.250000  0.250000
//     ...
//   URL http://example.com/motif1

// ===========================================================================
// Data structures
// ===========================================================================

///|
/// A MEME motif with position-specific probability matrix (PSPM).
pub struct MemeMotif {
  pub name : String // motif identifier (e.g., "motif1" or "1")
  pub mut alt_name : String // alternative name (may be empty)
  pub alphabet : String // e.g., "ACGT" or "ACDEFGHIKLMNPQRSTVWY"
  pub mut width : Int // motif width (number of positions)
  pub mut nsites : Int // number of sites used to build the motif
  pub mut evalue : Double // E-value of the motif
  pub pspm : Array[Array[Double]] // position-specific probability matrix (width x alphabet_len)
  pub mut url : String // optional motif URL
  pub mut strands : String // strand info (e.g., "+ -" or "+")
}

///|
/// Create a new MemeMotif.
pub fn MemeMotif::new(
  name : String,
  pspm : Array[Array[Double]],
  alphabet~ : String = "ACGT",
) -> MemeMotif {
  let width = if pspm.length() > 0 { pspm.length() } else { 0 }
  MemeMotif::{
    name,
    alt_name: "",
    alphabet,
    width,
    nsites: 0,
    evalue: 0.0,
    pspm,
    url: "",
    strands: "+ -",
  }
}

///|
/// Get the consensus sequence of the motif.
pub fn MemeMotif::consensus(self : MemeMotif) -> String {
  let result = StringBuilder::new()
  let alphabet = self.alphabet_to_array()
  let mut pos = 0
  while pos < self.pspm.length() {
    let row = self.pspm[pos]
    let mut max_prob = -1.0
    let mut max_idx = 0
    let mut i = 0
    while i < row.length() && i < alphabet.length() {
      if row[i] > max_prob {
        max_prob = row[i]
        max_idx = i
      }
      i = i + 1
    }
    result.write_string(alphabet[max_idx])
    pos = pos + 1
  }
  result.to_string()
}

///|
/// Convert alphabet string to array of single-character strings.
pub fn MemeMotif::alphabet_to_array(self : MemeMotif) -> Array[String] {
  let result = Array::new()
  let mut i = 0
  while i < self.alphabet.length() {
    result.push(self.alphabet[i:i + 1].to_owned())
    i = i + 1
  }
  result
}

///|
/// Get the probability of a specific letter at a position.
pub fn MemeMotif::probability(
  self : MemeMotif,
  letter : String,
  position : Int,
) -> Double {
  if position < 0 || position >= self.pspm.length() {
    return 0.0
  }
  let idx = meme_index_of(self.alphabet, letter)
  if idx < 0 || idx >= self.pspm[position].length() {
    return 0.0
  }
  self.pspm[position][idx]
}

///|
/// Get the information content (bits) of the motif.
pub fn MemeMotif::information_content(self : MemeMotif) -> Double {
  let alphabet_len = self.alphabet.length().to_double()
  let background = 1.0 / alphabet_len
  let mut total_ic = 0.0
  for row in self.pspm {
    for prob in row {
      if prob > 0.0 {
        total_ic = total_ic + prob * @math.ln(prob / background) / @math.ln(2.0)
      }
    }
  }
  total_ic
}

///|
/// A MEME record holding parsed motifs and metadata.
pub struct MemeRecord {
  pub mut version : String // MEME version string (e.g., "4")
  pub mut alphabet : String // alphabet (e.g., "ACGT")
  pub mut strands : String // strand info
  pub background_frequencies : Array[(String, Double)] // background letter frequencies
  pub motifs : Array[MemeMotif] // list of motifs
}

///|
/// Create a new empty MemeRecord.
pub fn MemeRecord::new() -> MemeRecord {
  MemeRecord::{
    version: "4",
    alphabet: "ACGT",
    strands: "+ -",
    background_frequencies: Array::new(),
    motifs: Array::new(),
  }
}

///|
/// Get the number of motifs.
pub fn MemeRecord::num_motifs(self : MemeRecord) -> Int {
  self.motifs.length()
}

///|
/// Find a motif by name.
pub fn MemeRecord::find(self : MemeRecord, name : String) -> MemeMotif? {
  for m in self.motifs {
    if m.name == name || m.alt_name == name {
      return Some(m)
    }
  }
  None
}

// ===========================================================================
// Parsing
// ===========================================================================

///|
/// Parse MEME minimal motif format text into a MemeRecord.
pub fn meme_parse(text : String) -> MemeRecord {
  let record = MemeRecord::new()
  let lines = text.split("\n").to_array()
  let mut i = 0
  let n = lines.length()
  let mut in_background = false
  let mut in_matrix = false
  let mut matrix_rows : Array[Array[Double]] = Array::new()
  let mut expected_width = 0
  // Motif building variables
  let mut has_motif = false
  let mut m_name = ""
  let mut m_alt_name = ""
  let mut m_alphabet = "ACGT"
  let mut m_nsites = 0
  let mut m_evalue = 0.0
  let mut m_url = ""
  let mut m_strands = "+ -"
  while i < n {
    let line = lines[i].trim().to_string()
    if line.length() == 0 {
      i = i + 1
      continue
    }
    // Version line
    if line.starts_with("MEME version") {
      let parts = line.split(" ").to_array()
      if parts.length() >= 3 {
        record.version = parts[2].trim().to_string()
      }
      i = i + 1
      continue
    }
    // Alphabet line
    if line.starts_with("ALPHABET=") {
      record.alphabet = line.substring(start=9).trim().to_string()
      i = i + 1
      continue
    }
    // Strands line
    if line.starts_with("strands:") {
      record.strands = line.substring(start=8).trim().to_string()
      i = i + 1
      continue
    }
    // Background letter frequencies header
    if line.starts_with("Background letter frequencies") {
      in_background = true
      i = i + 1
      continue
    }
    // Parse background frequencies
    if in_background {
      in_background = false
      let parts = line.split(" ").to_array()
      let mut j = 0
      while j + 1 < parts.length() {
        let letter = parts[j].trim().to_string()
        let freq_str = parts[j + 1].trim().to_string()
        let freq_opt = parse_double(freq_str)
        let freq = match freq_opt {
          Some(v) => v
          None => 0.0
        }
        if freq > 0.0 {
          record.background_frequencies.push((letter, freq))
        }
        j = j + 2
      }
      i = i + 1
      continue
    }
    // Helper to finalize current motif
    fn finalize_motif(
      record : MemeRecord,
      has_motif : Bool,
      m_name : String,
      m_alt_name : String,
      m_alphabet : String,
      m_nsites : Int,
      m_evalue : Double,
      m_url : String,
      m_strands : String,
      matrix_rows : Array[Array[Double]],
    ) -> Unit {
      if has_motif {
        let m = MemeMotif::new(m_name, matrix_rows, alphabet=m_alphabet)
        m.alt_name = m_alt_name
        m.nsites = m_nsites
        m.evalue = m_evalue
        m.url = m_url
        m.strands = m_strands
        m.width = matrix_rows.length()
        record.motifs.push(m)
      }
    }
    // MOTIF line
    if line.starts_with("MOTIF") {
      // If we were building a motif, save it first
      if has_motif {
        finalize_motif(
          record,
          has_motif,
          m_name,
          m_alt_name,
          m_alphabet,
          m_nsites,
          m_evalue,
          m_url,
          m_strands,
          matrix_rows,
        )
      }
      in_matrix = false
      matrix_rows = Array::new()
      let parts = line.split(" ").to_array()
      m_name = if parts.length() >= 2 { parts[1].trim().to_string() } else { "" }
      m_alt_name = if parts.length() >= 3 {
        parts[2].trim().to_string()
      } else {
        ""
      }
      m_alphabet = record.alphabet
      m_strands = record.strands
      m_nsites = 0
      m_evalue = 0.0
      m_url = ""
      has_motif = true
      i = i + 1
      continue
    }
    // letter-probability matrix header
    if line.starts_with("letter-probability matrix:") {
      in_matrix = true
      matrix_rows = Array::new()
      // Parse alength, w, nsites, E from the header
      let header = line.substring(start=26)
      expected_width = parse_meme_int_field(header, "w=")
      let _alength = parse_meme_int_field(header, "alength=")
      m_nsites = parse_meme_int_field(header, "nsites=")
      m_evalue = parse_meme_double_field(header, "E=")
      i = i + 1
      continue
    }
    // URL line (handle even after matrix is finalized, while still in motif)
    if has_motif && line.starts_with("URL") {
      m_url = line.substring(start=4).trim().to_string()
      in_matrix = false
      finalize_motif(
        record,
        has_motif,
        m_name,
        m_alt_name,
        m_alphabet,
        m_nsites,
        m_evalue,
        m_url,
        m_strands,
        matrix_rows,
      )
      has_motif = false
      i = i + 1
      continue
    }
    // Matrix data rows
    if in_matrix && has_motif {
      // Try to parse as matrix row (starts with spaces and contains numbers)
      let trimmed = line.trim().to_string()
      if trimmed.length() > 0 {
        let parts = trimmed.split(" ").to_array()
        let row : Array[Double] = Array::new()
        for p in parts {
          let s = p.trim().to_string()
          if s.length() > 0 {
            let val_opt = parse_double(s)
            let val = match val_opt {
              Some(v) => v
              None => 0.0
            }
            if val > 0.0 || s == "0" || s.starts_with("0.") {
              row.push(val)
            }
          }
        }
        if row.length() > 0 {
          matrix_rows.push(row)
        }
      }
      // Check if we've read enough rows - just stop matrix mode,
      // don't finalize yet (URL may follow)
      if expected_width > 0 && matrix_rows.length() >= expected_width {
        in_matrix = false
      }
      i = i + 1
      continue
    }
    i = i + 1
  }
  // Save last motif if still building
  if has_motif {
    let m = MemeMotif::new(m_name, matrix_rows, alphabet=m_alphabet)
    m.alt_name = m_alt_name
    m.nsites = m_nsites
    m.evalue = m_evalue
    m.url = m_url
    m.strands = m_strands
    m.width = matrix_rows.length()
    record.motifs.push(m)
  }
  record
}

///|
/// Find the index of a substring in a string, returns -1 if not found.
fn meme_index_of(s : String, sub : String) -> Int {
  let n = s.length()
  let kn = sub.length()
  if kn == 0 || kn > n {
    return -1
  }
  let mut i = 0
  while i <= n - kn {
    let mut matched = true
    let mut j = 0
    while j < kn {
      if s.unsafe_get(i + j) != sub.unsafe_get(j) {
        matched = false
        break
      }
      j = j + 1
    }
    if matched {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Parse an integer field from a MEME header line (e.g., "w= 8").
fn parse_meme_int_field(text : String, field : String) -> Int {
  let idx = meme_index_of(text, field)
  if idx >= 0 {
    let rest = text.substring(start=idx + field.length()).trim().to_string()
    let num_str = StringBuilder::new()
    let mut i = 0
    while i < rest.length() {
      let ch = rest[i:i + 1].to_string()
      if ch >= "0" && ch <= "9" {
        num_str.write_string(ch)
        i = i + 1
      } else {
        break
      }
    }
    let s = num_str.to_string()
    if s.length() > 0 {
      return parse_int(s)
    }
    return 0
  }
  0
}

///|
/// Parse a double field from a MEME header line (e.g., "E= 1.1e-03").
fn parse_meme_double_field(text : String, field : String) -> Double {
  let idx = meme_index_of(text, field)
  if idx >= 0 {
    let rest = text.substring(start=idx + field.length()).trim().to_string()
    let num_str = StringBuilder::new()
    let mut i = 0
    while i < rest.length() {
      let ch = rest[i:i + 1].to_string()
      if (ch >= "0" && ch <= "9") || ch == "." || ch == "-" || ch == "+" ||
        ch == "e" ||
        ch == "E" {
        num_str.write_string(ch)
        i = i + 1
      } else {
        break
      }
    }
    let s = num_str.to_string()
    if s.length() > 0 {
      let v_opt = parse_double(s)
      let v = match v_opt {
        Some(x) => x
        None => 0.0
      }
      if v > 0.0 || s == "0" || s.starts_with("0.") {
        return v
      }
      return 0.0
    }
    return 0.0
  }
  0.0
}

// ===========================================================================
// Serialization
// ===========================================================================

///|
/// Serialize a MemeMotif to MEME minimal motif format.
pub fn meme_motif_to_string(motif : MemeMotif) -> String {
  let sb = StringBuilder::new()
  sb.write_string("MOTIF " + motif.name)
  if motif.alt_name.length() > 0 {
    sb.write_string(" " + motif.alt_name)
  }
  sb.write_char('\n')
  // Letter-probability matrix header
  let alphabet_len = motif.alphabet.length()
  sb.write_string(
    "letter-probability matrix: alength= " +
    alphabet_len.to_string() +
    " w= " +
    motif.width.to_string() +
    " nsites= " +
    motif.nsites.to_string() +
    " E= " +
    motif.evalue.to_string() +
    "\n",
  )
  // Matrix rows
  for row in motif.pspm {
    let mut first = true
    for val in row {
      if !first {
        sb.write_string("  ")
      }
      sb.write_string(val.to_string())
      first = false
    }
    sb.write_char('\n')
  }
  // URL
  if motif.url.length() > 0 {
    sb.write_string("URL " + motif.url + "\n")
  }
  sb.to_string()
}

///|
/// Serialize a MemeRecord to MEME minimal motif format.
pub fn meme_record_to_string(record : MemeRecord) -> String {
  let sb = StringBuilder::new()
  sb.write_string("MEME version " + record.version + "\n\n")
  sb.write_string("ALPHABET= " + record.alphabet + "\n\n")
  if record.strands.length() > 0 {
    sb.write_string("strands: " + record.strands + "\n\n")
  }
  sb.write_string("Background letter frequencies\n")
  let mut bf_idx = 0
  while bf_idx < record.background_frequencies.length() {
    let entry = record.background_frequencies[bf_idx]
    let letter = entry.0
    let freq = entry.1
    sb.write_string(letter + " " + freq.to_string() + " ")
    bf_idx = bf_idx + 1
  }
  sb.write_char('\n')
  for motif in record.motifs {
    sb.write_char('\n')
    sb.write_string(meme_motif_to_string(motif))
  }
  sb.to_string()
}

// ===========================================================================
// Sample data
// ===========================================================================

///|
/// Create a sample MEME format text with 2 DNA motifs.
pub fn meme_sample_text() -> String {
  "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\nMOTIF motif1 FirstMotif\nletter-probability matrix: alength= 4 w= 6 nsites= 20 E= 1.1e-03\n  0.750000  0.000000  0.000000  0.250000\n  0.000000  0.000000  0.750000  0.250000\n  0.800000  0.000000  0.200000  0.000000\n  0.000000  0.900000  0.000000  0.100000\n  0.700000  0.100000  0.100000  0.100000\n  0.850000  0.050000  0.050000  0.050000\n\nMOTIF motif2 SecondMotif\nletter-probability matrix: alength= 4 w= 5 nsites= 15 E= 5.2e-02\n  0.100000  0.100000  0.100000  0.700000\n  0.000000  0.800000  0.000000  0.200000\n  0.600000  0.000000  0.300000  0.100000\n  0.000000  0.000000  0.900000  0.100000\n  0.200000  0.200000  0.200000  0.400000\n"
}

///|
/// Create a sample MemeMotif for testing.
pub fn meme_sample_motif() -> MemeMotif {
  let pspm = [
    [0.75, 0.0, 0.0, 0.25],
    [0.0, 0.0, 0.75, 0.25],
    [0.8, 0.0, 0.2, 0.0],
    [0.0, 0.9, 0.0, 0.1],
    [0.7, 0.1, 0.1, 0.1],
    [0.85, 0.05, 0.05, 0.05],
  ]
  let m = MemeMotif::new("motif1", pspm)
  m.alt_name = "FirstMotif"
  m.nsites = 20
  m.evalue = 0.0011
  m
}
