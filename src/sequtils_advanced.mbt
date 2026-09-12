///|
/// Bio.SeqUtils advanced functions, mirroring Biopython's
/// Bio/SeqUtils/__init__.py (gc_fraction, GC123, nt_search,
/// six_frame_translations) as of Biopython 1.88.
///
/// xGC_skew is intentionally omitted: it only draws a Tkinter GUI and has no
/// computational result beyond GC_skew (already provided by seq_utils.mbt).

///|
/// GC "weight" of each IUPAC nucleotide, mirroring Bio.SeqUtils._gc_values.
/// G/C count as 1.0, A/T/U/W as 0.0, two-way ambiguity as 0.5 and three-way
/// ambiguity as 1/3 or 2/3.
let sequtils_gc_values : Map[String, Double] = Map([
  ("G", 1.0),
  ("C", 1.0),
  ("A", 0.0),
  ("T", 0.0),
  ("U", 0.0),
  ("S", 1.0), // Strong interaction (3 H bonds) (G or C)
  ("W", 0.0), // Weak interaction (2 H bonds) (A or T)
  ("M", 0.5), // Amino (A or C)
  ("R", 0.5), // Purine (A or G)
  ("Y", 0.5), // Pyrimidine (T or C)
  ("K", 0.5), // Keto (G or T)
  ("V", 2.0 / 3.0), // Not T or U (A or C or G)
  ("B", 2.0 / 3.0), // Not A (C or G or T)
  ("H", 1.0 / 3.0), // Not G (A or C or T)
  ("D", 1.0 / 3.0), // Not C (A or G or T)
  ("X", 0.5), // Any nucleotide (A or C or G or T)
  ("N", 0.5), // Any nucleotide (A or C or G or T)
])

///|
/// Count the (case-sensitive) occurrences of every character of `chars`
/// in `seq`. Mirrors the sum(seq.count(x) ...) idiom in Biopython.
fn sequtils_count_chars(seq : String, chars : String) -> Int {
  let counts : Map[Int, Int] = Map([], capacity=64)
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq.unsafe_get(i).to_int()
    counts[c] = counts.get(c).or(0) + 1
  }
  let mut total = 0
  for j = 0; j < chars.length(); j = j + 1 {
    total = total + counts.get(chars.unsafe_get(j).to_int()).or(0)
  }
  total
}

///|
/// Calculate G+C percentage in seq (float between 0 and 1).
///
/// Ambiguous nucleotides in this context are those different from ATCGSWU
/// (S is G or C, and W is A or T).
///
/// - `ambiguous="remove"` (default): only count GCS and only include ACTGSWU
///   when calculating the sequence length. Equivalent to removing all
///   characters in the set BDHKMNRVXY before calculating the GC content.
/// - `ambiguous="ignore"`: treat only unambiguous nucleotides (GCS) as
///   counting towards the GC percentage, but include all characters when
///   calculating the sequence length.
/// - `ambiguous="weighted"`: use a "mean" value when counting ambiguous
///   characters, e.g. G/C count as 1, N/X as 0.5, D as 0.33 etc. See
///   Bio.SeqUtils._gc_values for the full list.
///
/// Raises `SeqError` for any other value of the ambiguous parameter.
/// Returns zero for an empty sequence.
pub fn gc_fraction(
  seq : String,
  ambiguous? : String = "remove",
) -> Double raise SeqError {
  let gc = sequtils_count_chars(seq, "CGScgs").to_double()
  let length = if ambiguous == "remove" {
    (sequtils_count_chars(seq, "CGScgs") + sequtils_count_chars(seq, "ATWUatwu")).to_double()
  } else if ambiguous == "ignore" || ambiguous == "weighted" {
    seq.length().to_double()
  } else {
    raise SeqError("ambiguous value '" + ambiguous + "' not recognized")
  }
  let mut gc_total = gc
  if ambiguous == "weighted" {
    for x in ["B", "D", "H", "K", "M", "N", "R", "V", "X", "Y"] {
      let c = sequtils_count_chars(seq, x + x.to_lower()).to_double()
      gc_total = gc_total + c * sequtils_gc_values.get(x).or(0.0)
    }
  }
  if length == 0.0 {
    0.0
  } else {
    gc_total / length
  }
}

///|
/// Result of [gc123]: total GC content plus GC content per codon position,
/// all as percentages between 0 and 100.
pub struct GC123Result {
  gc_all : Double
  gc1 : Double
  gc2 : Double
  gc3 : Double
} derive(Eq, Show)

///|
/// Calculate G+C content: total, for first, second and third codon positions.
///
/// Copes with mixed case sequences, but does NOT deal with ambiguous
/// nucleotides (only A/T/G/C are counted). Mirrors Biopython's GC123.
///
/// Note: Biopython raises ZeroDivisionError on an empty sequence; this
/// implementation returns all zeros instead.
pub fn gc123(seq : String) -> GC123Result {
  // counts[nt][pos]: nt 0..3 = A, T, G, C; pos 0..2 = codon position.
  let counts : FixedArray[FixedArray[Int]] = FixedArray::make(
    4,
    FixedArray::make(3, 0),
  )
  // Each row needs its own array; allocate them individually to avoid aliasing.
  for k = 0; k < 4; k = k + 1 {
    counts[k] = FixedArray::make(3, 0)
  }
  let nts : Array[Int] = [
    'A'.to_int(),
    'T'.to_int(),
    'G'.to_int(),
    'C'.to_int(),
  ]
  let length = seq.length()
  let mut i = 0
  while i < length {
    let mut pos = 0
    while pos < 3 && i + pos < length {
      let ci = seq.unsafe_get(i + pos).to_int()
      let upper = if ci >= 97 && ci <= 122 { ci - 32 } else { ci }
      for k = 0; k < 4; k = k + 1 {
        if upper == nts[k] {
          counts[k][pos] = counts[k][pos] + 1
        }
      }
      pos = pos + 1
    }
    i = i + 3
  }
  // Position-wise GC percentage
  let gc_at = fn(pos : Int) -> Double {
    let n = counts[0][pos] + counts[1][pos] + counts[2][pos] + counts[3][pos]
    if n == 0 {
      0.0
    } else {
      (counts[2][pos] + counts[3][pos]).to_double() * 100.0 / n.to_double()
    }
  }
  let gcall_gc = counts[2][0] +
    counts[2][1] +
    counts[2][2] +
    counts[3][0] +
    counts[3][1] +
    counts[3][2]
  let nall = counts[0][0] +
    counts[0][1] +
    counts[0][2] +
    counts[1][0] +
    counts[1][1] +
    counts[1][2] +
    counts[2][0] +
    counts[2][1] +
    counts[2][2] +
    counts[3][0] +
    counts[3][1] +
    counts[3][2]
  let gc_all = if nall == 0 {
    0.0
  } else {
    gcall_gc.to_double() * 100.0 / nall.to_double()
  }
  GC123Result::{ gc_all, gc1: gc_at(0), gc2: gc_at(1), gc3: gc_at(2), }
}

///|
/// Result of [nt_search]: the regex-style IUPAC pattern (as built by
/// Biopython) and the 0-based start positions of the (overlapping) matches.
pub struct NtSearchResult {
  pattern : String
  positions : Array[Int]
} derive(Eq, Show)

///|
/// IUPAC ambiguous DNA values, mirroring Bio.Data.IUPACData.ambiguous_dna_values.
/// Single-letter values mean the base itself; multi-letter values are the
/// set of bases the code may stand for.
let sequtils_ambiguous_dna_values : Map[String, String] = Map([
  ("A", "A"),
  ("C", "C"),
  ("G", "G"),
  ("T", "T"),
  ("M", "AC"),
  ("R", "AG"),
  ("W", "AT"),
  ("S", "CG"),
  ("Y", "CT"),
  ("K", "GT"),
  ("V", "ACG"),
  ("H", "ACT"),
  ("D", "AGT"),
  ("B", "CGT"),
  ("X", "GATC"),
  ("N", "GATC"),
])

///|
/// Search for a DNA subseq in seq, using ambiguous values
/// (like N = A or T or C or G, R = A or G etc.).
///
/// Searches only on the forward strand and is case sensitive (Biopython
/// builds an uppercase pattern that is matched literally against seq).
/// Overlapping matches are reported, and matching restarts one base after
/// the start of the previous match, mirroring Biopython's nt_search loop.
///
/// The returned [NtSearchResult] holds the regex-style pattern (e.g.
/// "[GATC]GA" for subseq "NGA") and the 0-based positions of all matches.
/// Raises `SeqError` for subseq characters that are not IUPAC DNA nucleotide
/// codes (Biopython raises KeyError in that case).
pub fn nt_search(
  seq : String,
  subseq : String,
) -> NtSearchResult raise SeqError {
  let pattern_buf = StringBuilder()
  let basesets : Array[Array[Int]] = []
  for i = 0; i < subseq.length(); i = i + 1 {
    let ci = subseq.unsafe_get(i).to_int()
    let upper = if ci >= 97 && ci <= 122 { ci - 32 } else { ci }
    let key = Char::from_int(upper).to_string()
    match sequtils_ambiguous_dna_values.get(key) {
      Some(value) => {
        if value.length() == 1 {
          pattern_buf.write_string(value)
        } else {
          pattern_buf.write_string("[" + value + "]")
        }
        let set : Array[Int] = []
        for j = 0; j < value.length(); j = j + 1 {
          set.push(value.unsafe_get(j).to_int())
        }
        basesets.push(set)
      }
      None =>
        raise SeqError(
          "not an IUPAC DNA nucleotide: '" +
          Char::from_int(ci).to_string() +
          "'",
        )
    }
  }
  let pattern = pattern_buf.to_string()
  let positions : Array[Int] = []
  let plen = basesets.length()
  if plen > 0 {
    let limit = seq.length() - plen
    let mut pos = 0
    while pos <= limit {
      let mut matched = true
      for j = 0; j < plen; j = j + 1 {
        let sc = seq.unsafe_get(pos + j).to_int()
        let set = basesets[j]
        let mut found = false
        for k = 0; k < set.length(); k = k + 1 {
          if sc == set[k] {
            found = true
            break
          }
        }
        if !found {
          matched = false
          break
        }
      }
      if matched {
        positions.push(pos)
      }
      // Biopython continues searching from pos + 1, allowing overlaps.
      pos = pos + 1
    }
  }
  NtSearchResult::{ pattern, positions, }
}

///|
/// Format a double with two decimal places (Python "%0.2f" semantics for
/// the non-negative values used here).
fn sequtils_fmt2(x : Double) -> String {
  let scaled = (x * 100.0 + 0.5).floor().to_int()
  let int_part = scaled / 100
  let frac = scaled - int_part * 100
  int_part.to_string() +
  "." +
  (if frac < 10 { "0" + frac.to_string() } else { frac.to_string() })
}

///|
/// Right-align an integer in a field of width 5 (Python "%5d").
fn sequtils_pad5(n : Int) -> String {
  let s = n.to_string()
  let buf = StringBuilder()
  for _i = 0; _i < 5 - s.length(); _i = _i + 1 {
    buf.write_char(' ')
  }
  buf.write_string(s)
  buf.to_string()
}

///|
/// Reverse a string (UTF-16 code units, safe for ASCII sequences).
fn sequtils_reverse_str(s : String) -> String {
  let n = s.length()
  let buf : FixedArray[UInt16] = FixedArray::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    buf[n - 1 - i] = s.unsafe_get(i)
  }
  unsafe_from_uint16_fixedarray(buf)
}

///|
/// Return pretty string showing the 6 frame translations and GC content.
///
/// Nice looking 6 frame translation with GC content - code from xbbtools,
/// similar to DNA Strider's six-frame translation. Mirrors Biopython's
/// six_frame_translations with the standard genetic code (table 1).
///
/// Raises `SeqError` if `genetic_code` is not 1 (only the standard table is
/// bundled) or if the sequence contains characters that cannot be translated.
pub fn six_frame_translations(
  seq : String,
  genetic_code? : Int = 1,
) -> String raise SeqError {
  if genetic_code != 1 {
    raise SeqError(
      "genetic code " +
      genetic_code.to_string() +
      " is not supported; only the standard code (table 1) is available",
    )
  }
  let s = Seq::new(seq)
  // Biopython uses the RNA complement tables when the sequence contains u,
  // and the DNA tables otherwise.
  let mut has_u = false
  for i = 0; i < length_of(seq); i = i + 1 {
    let c = seq.unsafe_get(i)
    if c == 'u'.to_int().to_uint16() || c == 'U'.to_int().to_uint16() {
      has_u = true
      break
    }
  }
  let anti = (if has_u {
      s.reverse_complement_rna()
    } else {
      s.reverse_complement()
    }).data
  let comp = (if has_u { s.complement_rna() } else { s.complement() }).data
  let length = seq.length()
  // frames[0..2] hold the +1/+2/+3 translations of seq, frames[3..5] hold the
  // reversed translations (frames -1/-2/-3 of Biopython's dict).
  let frames : FixedArray[String] = FixedArray::make(6, "")
  for i = 0; i < 3; i = i + 1 {
    let fragment_length = 3 * ((length - i) / 3)
    let fwd = Seq::new(seq.substring(start=i, end=i + fragment_length)).translate().data
    let rev_src = Seq::new(anti.substring(start=i, end=i + fragment_length)).translate().data
    frames[i] = fwd
    frames[3 + i] = sequtils_reverse_str(rev_src)
  }

  // Create header
  let short = if length > 20 {
    seq.substring(start=0, end=10) + " ... " + seq.substring(start=length - 10)
  } else {
    seq
  }
  let mut header = "GC_Frame:"
  for nt in ['a', 't', 'g', 'c'] {
    let upper = nt.to_int() - 32
    let mut count = 0
    for i = 0; i < length; i = i + 1 {
      if seq.unsafe_get(i).to_int() == upper {
        count = count + 1
      }
    }
    header = header + " " + nt.to_string() + ":" + count.to_string()
  }
  let gc = 100.0 * gc_fraction(seq, ambiguous="ignore")
  header = header +
    "\nSequence: " +
    short.to_lower() +
    ", " +
    length.to_string() +
    " nt, " +
    sequtils_fmt2(gc) +
    " %GC\n\n\n"
  let res = StringBuilder()
  res.write_string(header)

  // A frame slice formatted with the double-space separator, optionally
  // indented by `indent` leading spaces. Mirrors (" "*ind + "  ".join(...)).
  let frame_line = fn(frame : String, start : Int, indent : Int) -> String {
    let stop = if start + 20 < frame.length() {
      start + 20
    } else {
      frame.length()
    }
    let buf = StringBuilder()
    for _k = 0; _k < indent; _k = _k + 1 {
      buf.write_char(' ')
    }
    for j = start; j < stop; j = j + 1 {
      if j > start {
        buf.write_char(' ')
        buf.write_char(' ')
      }
      buf.write_char(frame.unsafe_get(j).unsafe_to_char())
    }
    buf.to_string()
  }

  let mut i = 0
  while i < length {
    let stop = if i + 60 < length { i + 60 } else { length }
    let subseq = seq.substring(start=i, end=stop)
    let csubseq = comp.substring(start=i, end=stop)
    let p = i / 3
    res.write_string((i + 1).to_string() + "/" + (i / 3 + 1).to_string() + "\n")
    res.write_string(frame_line(frames[2], p, 2) + "\n")
    res.write_string(frame_line(frames[1], p, 1) + "\n")
    res.write_string(frame_line(frames[0], p, 0) + "\n")
    // seq
    res.write_string(subseq.to_lower() + sequtils_pad5(gc.to_int()) + " %\n")
    res.write_string(csubseq.to_lower() + "\n")
    // - frames
    res.write_string(frame_line(frames[4], p, 0) + "\n")
    res.write_string(frame_line(frames[3], p, 1) + "\n")
    res.write_string(frame_line(frames[5], p, 2) + "\n\n")
    i = i + 60
  }
  res.to_string()
}

///|
/// Length helper kept local to avoid depending on String internals.
fn length_of(s : String) -> Int {
  s.length()
}
