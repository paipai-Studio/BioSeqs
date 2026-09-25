///|
/// pyfaidx sequence-level utilities.
///
/// Port of the sequence-oriented helpers from pyfaidx (0.9.0.4) that sit on
/// top of the `.fai` random-access layer already provided by `faidx.mbt`:
/// IUPAC complement / reverse / reverse-complement, GC content (strict,
/// IUPAC-weighted), FASTA line wrapping, and BED / UCSC region-string parsing.

///|
pub suberror PyfaidxError {
  PyfaidxError(String)
}

///|
fn[T] pyfaidx_fail(msg : String) -> T raise PyfaidxError {
  raise PyfaidxError(msg)
}

///|
fn pyfaidx_uppercase(s : String) -> String {
  let buf = StringBuilder()
  for i = 0; i < s.length(); i = i + 1 {
    let c = s.unsafe_get(i).to_int()
    let u = if c >= 97 && c <= 122 { c - 32 } else { c }
    buf.write_char((u).unsafe_to_char())
  }
  buf.to_string()
}

///|
fn pyfaidx_is_ws(c : Char) -> Bool {
  c == ' ' || c == '\t' || c == '\n' || c == '\r'
}

///|
fn pyfaidx_ws_tokens(s : String) -> Array[String] {
  let tokens : Array[String] = []
  let n = s.length()
  let mut i = 0
  while i < n {
    while i < n && pyfaidx_is_ws(s.unsafe_get(i).unsafe_to_char()) {
      i = i + 1
    }
    if i >= n {
      break
    }
    let start = i
    while i < n && !pyfaidx_is_ws(s.unsafe_get(i).unsafe_to_char()) {
      i = i + 1
    }
    tokens.push(s[start:i].to_owned())
  }
  tokens
}

///|
/// Find the first occurrence of a byte/char, or -1.
fn pyfaidx_index_of_char(s : String, target : Int) -> Int {
  for i = 0; i < s.length(); i = i + 1 {
    if s.unsafe_get(i).to_int() == target {
      return i
    }
  }
  -1
}

// ---------------------------------------------------------------------------
// IUPAC complement / reverse
// ---------------------------------------------------------------------------

///|
/// The IUPAC complement translation table used by pyfaidx.
fn pyfaidx_complement_source() -> String {
  "ACTGNactgnYRWSKMDVHBXyrwskmdvhbx"
}

///|
fn pyfaidx_complement_target() -> String {
  "TGACNtgacnRYWSMKHBDVXrywsmkhbdvx"
}

///|
/// Complement a DNA sequence using the full IUPAC ambiguity alphabet.
/// Raises on any character outside the alphabet (matching pyfaidx's `complement`).
pub fn pyfaidx_complement(seq : String) -> String raise PyfaidxError {
  let source = pyfaidx_complement_source()
  let target = pyfaidx_complement_target()
  let buf = StringBuilder()
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq.unsafe_get(i).unsafe_to_char()
    let mut found = false
    for j = 0; j < source.length(); j = j + 1 {
      if source.unsafe_get(j).unsafe_to_char() == c {
        buf.write_char(target.unsafe_get(j).unsafe_to_char())
        found = true
        break
      }
    }
    if !found {
      pyfaidx_fail("Sequence contains non-DNA character at position \{i + 1}")
    }
  }
  buf.to_string()
}

///|
/// Reverse a sequence string.
pub fn pyfaidx_reverse(seq : String) -> String {
  let buf = StringBuilder()
  let mut i = seq.length() - 1
  while i >= 0 {
    buf.write_char(seq.unsafe_get(i).unsafe_to_char())
    i = i - 1
  }
  buf.to_string()
}

///|
/// Reverse complement of a sequence (complement, then reverse).
pub fn pyfaidx_reverse_complement(seq : String) -> String raise PyfaidxError {
  pyfaidx_reverse(pyfaidx_complement(seq))
}

// ---------------------------------------------------------------------------
// GC content
// ---------------------------------------------------------------------------

///|
/// GC content (`(G + g + C + c) / length`), denominator is the full sequence.
pub fn pyfaidx_gc(seq : String) -> Double {
  let mut gc = 0.0
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq.unsafe_get(i).unsafe_to_char()
    if c == 'G' || c == 'g' || c == 'C' || c == 'c' {
      gc = gc + 1.0
    }
  }
  if seq.length() == 0 {
    0.0
  } else {
    gc / seq.length().to_double()
  }
}

///|
/// GC content ignoring non-ACGT characters (`(G + C) / len(filtered)`).
pub fn pyfaidx_gc_strict(seq : String) -> Double {
  let upper = pyfaidx_uppercase(seq)
  let mut gc = 0.0
  let mut n = 0
  for i = 0; i < upper.length(); i = i + 1 {
    let c = upper.unsafe_get(i).unsafe_to_char()
    if c == 'A' || c == 'C' || c == 'G' || c == 'T' {
      n = n + 1
      if c == 'G' || c == 'C' {
        gc = gc + 1.0
      }
    }
  }
  if n == 0 {
    0.0
  } else {
    gc / n.to_double()
  }
}

///|
/// GC content with IUPAC ambiguity codes weighted by their GC likelihood.
pub fn pyfaidx_gc_iupac(seq : String) -> Double {
  let upper = pyfaidx_uppercase(seq)
  let mut gc = 0.0
  let mut n = 0
  for i = 0; i < upper.length(); i = i + 1 {
    let c = upper.unsafe_get(i).unsafe_to_char()
    // filter: ACGTMRWSYKVHDBN
    if c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'M' || c == 'R' ||
    c == 'W' || c == 'S' || c == 'Y' || c == 'K' || c == 'V' || c == 'H' ||
    c == 'D' || c == 'B' || c == 'N' {
      n = n + 1
      if c == 'S' || c == 'C' || c == 'G' {
        gc = gc + 1.0
      } else if c == 'B' || c == 'V' {
        gc = gc + 0.67
      } else if c == 'M' || c == 'R' || c == 'Y' || c == 'K' {
        gc = gc + 0.5
      } else if c == 'H' || c == 'D' {
        gc = gc + 0.33
      } else if c == 'N' {
        gc = gc + 0.25
      }
    }
  }
  if n == 0 {
    0.0
  } else {
    gc / n.to_double()
  }
}

// ---------------------------------------------------------------------------
// Line wrapping + region parsing
// ---------------------------------------------------------------------------

///|
/// Wrap a sequence into lines of `n` characters, each followed by a newline.
pub fn pyfaidx_wrap_sequence(n : Int, seq : String) -> String {
  let buf = StringBuilder()
  let len = seq.length()
  let mut i = 0
  while i < len {
    let end = if i + n < len { i + n } else { len }
    buf.write_string(seq[i:end].to_owned())
    buf.write_char('\n')
    i = end
  }
  buf.to_string()
}

///|
/// Parse a BED entry into `(name, start, end)`; a `#`-comment yields
/// `(None, None, None)`.
pub fn pyfaidx_bed_split(entry : String) -> (String?, Int?, Int?) raise PyfaidxError {
  if entry.length() > 0 && entry.unsafe_get(0).unsafe_to_char() == '#' {
    return (None, None, None)
  }
  let tokens = pyfaidx_ws_tokens(entry)
  if tokens.length() < 3 {
    pyfaidx_fail("Malformed BED entry")
  }
  (Some(tokens[0]), Some(parse_int(tokens[1])), Some(parse_int(tokens[2])))
}

///|
/// Parse a UCSC region string `name:start-end` (1-based inclusive) into
/// `(name, start-1, end)` (0-based half-open). Missing interval → `None`.
pub fn pyfaidx_ucsc_split(region : String) -> (String, Int?, Int?) {
  let colon = pyfaidx_index_of_char(region, ':'.to_int())
  if colon < 0 {
    return (region, None, None)
  }
  let rname = region[0:colon].to_owned()
  let interval = region[colon + 1:region.length()].to_owned()
  let dash = pyfaidx_index_of_char(interval, '-'.to_int())
  if dash < 0 {
    return (rname, None, None)
  }
  let start = parse_int(interval[0:dash].to_owned()) - 1
  let end = parse_int(interval[dash + 1:interval.length()].to_owned())
  (rname, Some(start), Some(end))
}