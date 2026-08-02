// Bio.SubsMat.FreqTable - Frequency table for sequence analysis.
//
// Inspired by Biopython's Bio.SubsMat.FreqTable module, which provides a
// class to handle frequency tables or letter count files. This pure-MoonBit
// module provides the same core functionality: build a frequency table from
// either a count dictionary or a frequency dictionary, parse count/freq text
// files, and query letters' frequencies/counts. Useful as input to
// substitution-matrix construction (e.g. expected frequency table for
// log-odds matrices) and for general sequence composition analysis.
//
// Reference:
//   https://biopython.org/docs/1.79/api/Bio.SubsMat.FreqTable.html

// ============================================================================
// FreqTableType - input dictionary semantics
// ============================================================================

///|
/// Identifies whether an input dictionary holds raw counts or
/// pre-computed frequencies.
pub enum FreqTableType {
  COUNT
  FREQ
}

///|
/// Construct a COUNT table-type tag.
pub fn freq_table_count_type() -> FreqTableType {
  FreqTableType::COUNT
}

///|
/// Construct a FREQ table-type tag.
pub fn freq_table_freq_type() -> FreqTableType {
  FreqTableType::FREQ
}

// ============================================================================
// FreqTable - the frequency table
// ============================================================================

///|
/// A frequency table for sequence alphabet letters.
///
/// `data` always holds letter frequencies (summing to ~1.0 when the input
/// was a COUNT table). `count` holds the raw counts when available, and is
/// empty when the table was constructed directly from frequencies.
pub struct FreqTable {
  alphabet : Array[String]
  data : Map[String, Double] // letter -> frequency
  count : Map[String, Int] // letter -> count (may be empty)
}

///|
/// Build a `FreqTable` from an input dictionary.
///
/// When `dict_type` is `COUNT`, the values in `in_dict` are treated as raw
/// integer counts and frequencies are derived as count / total. When
/// `dict_type` is `FREQ`, the values are taken as frequencies directly and
/// `count` is left empty.
///
/// If `alphabet` is empty, it is inferred from the dictionary keys (sorted).
pub fn FreqTable::new(
  in_dict : Map[String, Double],
  dict_type : FreqTableType,
  alphabet : Array[String]?
) -> FreqTable {
  let data : Map[String, Double] = Map([], capacity=16)
  let count : Map[String, Int] = Map([], capacity=16)
  let alph : Array[String] = Array::new()
  // Determine the alphabet.
  match alphabet {
    Some(a) =>
      for ch in a {
        alph.push(ch)
      }
    None => for k in in_dict.keys() {
      alph.push(k)
    }
  }
  // Sort alphabet for deterministic order.
  alph.sort_by(fn(a : String, b : String) -> Int {
    if a < b { -1 } else if a > b { 1 } else { 0 }
  })
  match dict_type {
    COUNT => {
      // Compute total count first.
      let mut total = 0.0
      for k in in_dict.keys() {
        total = total + in_dict[k]
      }
      if total <= 0.0 {
        // Degenerate: all zero counts; set frequencies to 0.
        for k in in_dict.keys() {
          data[k] = 0.0
          count[k] = in_dict[k].to_int()
        }
      } else {
        for k in in_dict.keys() {
          let c = in_dict[k]
          count[k] = c.to_int()
          data[k] = c / total
        }
      }
    }
    FREQ =>
      for k in in_dict.keys() {
        data[k] = in_dict[k]
      }
  }
  // Ensure all alphabet letters have a frequency entry (default 0.0).
  // Only populate `count` for COUNT tables, so FREQ tables report
  // `has_counts()` as false.
  for ch in alph {
    if !data.contains(ch) {
      data[ch] = 0.0
    }
  }
  match dict_type {
    COUNT =>
      for ch in alph {
        if !count.contains(ch) {
          count[ch] = 0
        }
      }
    FREQ => ()
  }
  FreqTable::{ alphabet: alph, data, count }
}

///|
/// The alphabet letters used as indices into the table.
pub fn FreqTable::alphabet(self : FreqTable) -> Array[String] {
  self.alphabet
}

///|
/// The frequency dictionary (letter -> frequency).
pub fn FreqTable::data(self : FreqTable) -> Map[String, Double] {
  self.data
}

///|
/// The count dictionary (letter -> count). May be empty when the table was
/// constructed from a frequency file.
pub fn FreqTable::count(self : FreqTable) -> Map[String, Int] {
  self.count
}

///|
/// Look up the frequency for a letter. Returns 0.0 for letters not in the
/// table.
pub fn FreqTable::frequency(self : FreqTable, letter : String) -> Double {
  self.data.get_or_default(letter, 0.0)
}

///|
/// Look up the count for a letter. Returns 0 for letters not in the table
/// or when no counts are available.
pub fn FreqTable::count_of(self : FreqTable, letter : String) -> Int {
  self.count.get_or_default(letter, 0)
}

///|
/// True if the table stores raw counts (i.e. `count` is non-empty).
pub fn FreqTable::has_counts(self : FreqTable) -> Bool {
  self.count.length() > 0
}

///|
/// Total number of letters represented by the underlying counts.
/// Returns 0 when no counts are stored.
pub fn FreqTable::total_count(self : FreqTable) -> Int {
  let mut total = 0
  for k in self.count.keys() {
    total = total + self.count[k]
  }
  total
}

///|
/// Number of letters in the alphabet.
pub fn FreqTable::size(self : FreqTable) -> Int {
  self.alphabet.length()
}

///|
/// True if the table contains the given letter.
pub fn FreqTable::contains(self : FreqTable, letter : String) -> Bool {
  self.data.contains(letter)
}

///|
/// All letters in the table (sorted alphabetically).
pub fn FreqTable::letters(self : FreqTable) -> Array[String] {
  let out : Array[String] = Array::new()
  for ch in self.alphabet {
    out.push(ch)
  }
  out
}

///|
/// Pretty-print the table as a multi-line string, one letter per line:
///   "<letter> <frequency>" (and "<count>" appended when counts exist).
pub fn FreqTable::to_string(self : FreqTable) -> String {
  let sb = StringBuilder::new()
  for ch in self.alphabet {
    let freq = self.data.get_or_default(ch, 0.0)
    sb.write_string(ch)
    sb.write_char('\t')
    sb.write_string(freq.to_string())
    if self.has_counts() {
      let c = self.count.get_or_default(ch, 0)
      sb.write_char('\t')
      sb.write_string(c.to_string())
    }
    sb.write_char('\n')
  }
  sb.to_string()
}

// ============================================================================
// Parsers: read_count / read_freq
// ============================================================================

///|
/// Split a text blob into non-empty trimmed tokens, handling whitespace
/// (spaces, tabs, newlines) uniformly.
fn freq_table_tokenize(text : String) -> Array[String] {
  let tokens : Array[String] = Array::new()
  let n = text.length()
  let mut i = 0
  let mut buf = StringBuilder::new()
  while i < n {
    let ch = substring(text, i, 1)
    if ch == " " || ch == "\t" || ch == "\n" || ch == "\r" {
      if buf.to_string().length() > 0 {
        tokens.push(buf.to_string())
        buf = StringBuilder::new()
      }
    } else {
      buf.write_string(ch)
    }
    i = i + 1
  }
  let tail = buf.to_string()
  if tail.length() > 0 {
    tokens.push(tail)
  }
  tokens
}

///|
/// Parse a count file's text contents into a `FreqTable`.
///
/// The expected format is whitespace-separated `<letter> <count>` pairs,
/// one or more per line. Counts are parsed as integers (decimals are
/// truncated). The alphabet is inferred from the letters seen.
///
/// Example input:
///   "A 50\nC 37\nG 23\nT 58"
pub fn freq_table_read_count(text : String) -> FreqTable {
  let tokens = freq_table_tokenize(text)
  let in_dict : Map[String, Double] = Map([], capacity=16)
  let mut i = 0
  while i + 1 < tokens.length() {
    let letter = tokens[i]
    let count_str = tokens[i + 1]
    let count_val = parse_double(count_str)
    match count_val {
      Some(v) => in_dict[letter] = v
      None => in_dict[letter] = 0.0
    }
    i = i + 2
  }
  FreqTable::new(in_dict, FreqTableType::COUNT, None)
}

///|
/// Parse a frequency file's text contents into a `FreqTable`.
///
/// The expected format is whitespace-separated `<letter> <frequency>` pairs.
/// Frequencies are parsed as doubles. The alphabet is inferred from the
/// letters seen. The `count` dictionary is left empty (the source file
/// contains only frequencies).
///
/// Example input:
///   "A 0.2976\nC 0.2202\nG 0.1369\nT 0.3452"
pub fn freq_table_read_freq(text : String) -> FreqTable {
  let tokens = freq_table_tokenize(text)
  let in_dict : Map[String, Double] = Map([], capacity=16)
  let mut i = 0
  while i + 1 < tokens.length() {
    let letter = tokens[i]
    let freq_str = tokens[i + 1]
    let freq_val = parse_double(freq_str)
    match freq_val {
      Some(v) => in_dict[letter] = v
      None => in_dict[letter] = 0.0
    }
    i = i + 2
  }
  FreqTable::new(in_dict, FreqTableType::FREQ, None)
}

// ============================================================================
// Convenience constructors
// ============================================================================

///|
/// Build a `FreqTable` from a count dictionary (Map[String, Int]).
pub fn freq_table_from_counts(
  counts : Map[String, Int]
) -> FreqTable {
  let in_dict : Map[String, Double] = Map([], capacity=counts.size())
  for k in counts.keys() {
    in_dict[k] = counts[k].to_double()
  }
  FreqTable::new(in_dict, FreqTableType::COUNT, None)
}

///|
/// Build a `FreqTable` from a frequency dictionary (Map[String, Double]).
pub fn freq_table_from_freqs(
  freqs : Map[String, Double]
) -> FreqTable {
  FreqTable::new(freqs, FreqTableType::FREQ, None)
}

///|
/// Build a `FreqTable` by counting letters in a sequence string.
/// Non-alphabet characters are ignored (treated as separate letters only if
/// they appear in the sequence; the user may pass an explicit alphabet to
/// restrict counting).
pub fn freq_table_from_sequence(
  sequence : String,
  alphabet : Array[String]?
) -> FreqTable {
  let counts : Map[String, Int] = Map([], capacity=16)
  let restrict : Map[String, Bool] = Map([], capacity=0)
  match alphabet {
    Some(a) =>
      for ch in a {
        restrict[ch] = true
      }
    None => ()
  }
  let has_restrict = restrict.size() > 0
  let n = sequence.length()
  let mut i = 0
  while i < n {
    let ch = substring(sequence, i, 1)
    if has_restrict && !restrict.contains(ch) {
      i = i + 1
      continue
    }
    let cur = counts.get_or_default(ch, 0)
    counts[ch] = cur + 1
    i = i + 1
  }
  freq_table_from_counts(counts)
}

// ============================================================================
// Statistical helpers
// ============================================================================

///|
/// Compute the Shannon entropy (in bits) of the frequency distribution.
/// H = -sum(p_i * log2(p_i)). Letters with zero frequency are skipped.
pub fn freq_table_entropy(table : FreqTable) -> Double {
  let mut h = 0.0
  for k in table.data.keys() {
    let p = table.data[k]
    if p > 0.0 {
      h = h - p * (@math.ln(p) / @math.ln(2.0))
    }
  }
  h
}

///|
/// Relative entropy (Kullback-Leibler divergence) D_KL(table || other),
/// summed over letters present in `table`. Letters missing in `other` are
/// treated as having frequency 0 (and contribute 0 to the sum, since
/// 0 * log(0/0) is taken as 0).
pub fn freq_table_relative_entropy(
  table : FreqTable,
  other : FreqTable
) -> Double {
  let mut d = 0.0
  for k in table.data.keys() {
    let p = table.data[k]
    if p <= 0.0 {
      continue
    }
    let q = other.data.get_or_default(k, 0.0)
    if q <= 0.0 {
      continue
    }
    d = d + p * (@math.ln(p / q) / @math.ln(2.0))
  }
  d
}

///|
/// Normalize the frequency dictionary so that values sum to 1.0. Returns a
/// new `FreqTable`. Counts are preserved when present.
pub fn freq_table_normalize(table : FreqTable) -> FreqTable {
  let mut total = 0.0
  for k in table.data.keys() {
    total = total + table.data[k]
  }
  let out : Map[String, Double] = Map([], capacity=table.data.size())
  for k in table.data.keys() {
    if total > 0.0 {
      out[k] = table.data[k] / total
    } else {
      out[k] = 0.0
    }
  }
  let alph = table.alphabet
  let result = FreqTable::new(out, FreqTableType::FREQ, Some(alph))
  // Preserve counts if the source had them.
  if table.has_counts() {
    for k in table.count.keys() {
      result.count[k] = table.count[k]
    }
  }
  result
}

///|
/// Jensen-Shannon distance between two frequency tables, using the average
/// distribution M = 0.5 * (P + Q). Returns JS = 0.5 * KL(P||M) + 0.5 * KL(Q||M).
/// Both tables should share the same alphabet; letters missing from one are
/// treated as 0 in the other.
pub fn freq_table_js_distance(
  table : FreqTable,
  other : FreqTable
) -> Double {
  // Build the union of letters.
  let letters : Map[String, Bool] = Map([], capacity=32)
  for k in table.data.keys() {
    letters[k] = true
  }
  for k in other.data.keys() {
    letters[k] = true
  }
  let mut js = 0.0
  for k in letters.keys() {
    let p = table.data.get_or_default(k, 0.0)
    let q = other.data.get_or_default(k, 0.0)
    let m = 0.5 * (p + q)
    if p > 0.0 && m > 0.0 {
      js = js + 0.5 * p * (@math.ln(p / m) / @math.ln(2.0))
    }
    if q > 0.0 && m > 0.0 {
      js = js + 0.5 * q * (@math.ln(q / m) / @math.ln(2.0))
    }
  }
  if js < 0.0 {
    js = 0.0
  }
  js.sqrt()
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// A sample DNA count file string (whitespace-separated), matching the
/// Biopython documentation example.
pub fn freq_table_sample_count_text() -> String {
  "A 50\nC 37\nG 23\nT 58"
}

///|
/// A sample DNA frequency file string (whitespace-separated), matching the
/// Biopython documentation example.
pub fn freq_table_sample_freq_text() -> String {
  "A 0.2976\nC 0.2202\nG 0.1369\nT 0.3452"
}

///|
/// A sample protein count file string for the 20 standard amino acids.
pub fn freq_table_sample_protein_count_text() -> String {
  "A 21\nC 5\nD 18\nE 22\nF 13\nG 25\nH 7\nI 16\nK 17\nL 26\nM 9\nN 14\nP 15\nQ 12\nR 16\nS 22\nT 19\nV 20\nW 4\nY 11"
}

///|
/// Build a sample DNA `FreqTable` from counts, useful for demos and tests.
pub fn freq_table_sample_dna() -> FreqTable {
  freq_table_read_count(freq_table_sample_count_text())
}

///|
/// Build a sample DNA `FreqTable` from frequencies, useful for demos and tests.
pub fn freq_table_sample_dna_freq() -> FreqTable {
  freq_table_read_freq(freq_table_sample_freq_text())
}

///|
/// Build a sample protein `FreqTable` from counts.
pub fn freq_table_sample_protein() -> FreqTable {
  freq_table_read_count(freq_table_sample_protein_count_text())
}
