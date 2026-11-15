// ============================================================================
// CIGAR utility functions — helpers shared across BAM operations.
// Pure functions, zero struct / enum modification.
// ============================================================================

///|
/// Total query-consuming length of a single CIGAR op (M/I/S/=/X).
/// Mismatches (N / H / P) are excluded. Matches samtools `bam_cigar_opqlen`.
pub fn cigar_op_qlen(op : BamCigarOp) -> Int32 {
  match op {
    BamCigarOp::Match
    | BamCigarOp::Insert
    | BamCigarOp::SoftClip
    | BamCigarOp::Equal
    | BamCigarOp::Diff => 1
    _ => 0
  }
}

///|
/// Total reference-consuming length of a single CIGAR op (M/D/N/=/X).
/// Insertions and clips are excluded. Matches samtools `bam_cigar_oprlen`.
pub fn cigar_op_rlen(op : BamCigarOp) -> Int32 {
  match op {
    BamCigarOp::Match
    | BamCigarOp::Delete
    | BamCigarOp::Skip
    | BamCigarOp::Equal
    | BamCigarOp::Diff => 1
    _ => 0
  }
}

///|
/// Op is query-consuming (qlen > 0).
pub fn cigar_op_consumes_query(op : BamCigarOp) -> Bool {
  cigar_op_qlen(op) == 1
}

///|
/// Op is reference-consuming (rlen > 0).
pub fn cigar_op_consumes_reference(op : BamCigarOp) -> Bool {
  cigar_op_rlen(op) == 1
}

///|
/// Op consumes neither query nor reference (H / P).
pub fn cigar_op_is_clipping(op : BamCigarOp) -> Bool {
  op == BamCigarOp::SoftClip || op == BamCigarOp::HardClip
}

///|
/// Is the CIGAR operation an insertion relative to the reference (I)?
pub fn cigar_op_is_insertion(op : BamCigarOp) -> Bool {
  op == BamCigarOp::Insert
}

///|
/// Is the CIGAR operation a deletion relative to the reference (D / N)?
/// `is_deletion` is true for both deletion and skip (intron).
pub fn cigar_op_is_deletion(op : BamCigarOp) -> Bool {
  op == BamCigarOp::Delete || op == BamCigarOp::Skip
}

///|
/// Map char code → BamCigarOp. `=5` op code is '='; 'X' maps to Diff.
/// Returns None for unknown chars. Matches samtools cigar table.
pub fn cigar_from_char(c : Char) -> BamCigarOp? {
  match c {
    'M' => Some(BamCigarOp::Match)
    'I' => Some(BamCigarOp::Insert)
    'D' => Some(BamCigarOp::Delete)
    'N' => Some(BamCigarOp::Skip)
    'S' => Some(BamCigarOp::SoftClip)
    'H' => Some(BamCigarOp::HardClip)
    'P' => Some(BamCigarOp::Pad)
    '=' => Some(BamCigarOp::Equal)
    'X' => Some(BamCigarOp::Diff)
    _ => None
  }
}

///|
/// Map BamCigarOp → char code.
pub fn cigar_to_char(op : BamCigarOp) -> Char {
  match op {
    BamCigarOp::Match => 'M'
    BamCigarOp::Insert => 'I'
    BamCigarOp::Delete => 'D'
    BamCigarOp::Skip => 'N'
    BamCigarOp::SoftClip => 'S'
    BamCigarOp::HardClip => 'H'
    BamCigarOp::Pad => 'P'
    BamCigarOp::Equal => '='
    BamCigarOp::Diff => 'X'
  }
}

///|
/// Map BamCigarOp → samtools op code byte (0=M, 1=I, 2=D, ..., 8=X).
pub fn cigar_to_code(op : BamCigarOp) -> UInt16 {
  match op {
    BamCigarOp::Match => 0
    BamCigarOp::Insert => 1
    BamCigarOp::Delete => 2
    BamCigarOp::Skip => 3
    BamCigarOp::SoftClip => 4
    BamCigarOp::HardClip => 5
    BamCigarOp::Pad => 6
    BamCigarOp::Equal => 7
    BamCigarOp::Diff => 8
  }
}

///|
/// Map samtools op code byte → BamCigarOp. Returns None for unknown codes.
pub fn cigar_from_code(code : UInt16) -> BamCigarOp? {
  if code == 0 { Some(BamCigarOp::Match) }
  else if code == 1 { Some(BamCigarOp::Insert) }
  else if code == 2 { Some(BamCigarOp::Delete) }
  else if code == 3 { Some(BamCigarOp::Skip) }
  else if code == 4 { Some(BamCigarOp::SoftClip) }
  else if code == 5 { Some(BamCigarOp::HardClip) }
  else if code == 6 { Some(BamCigarOp::Pad) }
  else if code == 7 { Some(BamCigarOp::Equal) }
  else if code == 8 { Some(BamCigarOp::Diff) }
  else { None }
}

///|
/// Parse a standard text CIGAR string like `"50M10I30M"` into an element
/// array. Ignores trailing soft/hard clips gracefully. Returns an empty
/// array on empty input.
pub fn cigar_parse(text : String) -> Array[BamCigarElement] {
  let out : Array[BamCigarElement] = Array(capacity=0)
  if text.length() == 0 {
    return out
  }
  let mut num_start = 0
  let mut i = 0
  while i < text.length() {
    let c = text.unsafe_get(i).unsafe_to_char()
    // digits accumulate; non-digit closes one op
    if (c >= '0' && c <= '9') || (c >= 'A' && c <= 'Z') || c == '=' {
      if (c >= 'A' && c <= 'Z') || c == '=' {
        let op = cigar_from_char(c)
        match op {
          Some(o) => {
            let num_str = text[num_start:i].to_owned()
            let len = parse_int(num_str)
            out.push(BamCigarElement::{ op: o, length: len, })
          }
          None => ()
        }
        num_start = i + 1
      }
    }
    i = i + 1
  }
  out
}

///|
/// Serialize a CIGAR element array back into a text string like
/// `"50M10I30M"`. Hard/soft clips are preserved.
pub fn cigar_to_string(cigar : Array[BamCigarElement]) -> String {
  let sb = StringBuilder()
  for el in cigar {
    sb.write_string(el.length.to_string())
    sb.write_char(cigar_to_char(el.op))
  }
  sb.to_string()
}

///|
/// Sum query-consuming length over an element array. Matches Pysam
/// `AlignedSegment.infer_query_length` (does **not** include H clips).
pub fn cigar_query_length(cigar : Array[BamCigarElement]) -> Int32 {
  let mut len = 0
  for el in cigar {
    len = len + el.length * cigar_op_qlen(el.op)
  }
  len
}

///|
/// Sum reference-consuming length over an element array (M/D/N/=/X).
/// Matches samtools `bam_cigar2rlen` / Pysam
/// `AlignedSegment.reference_length`.
pub fn cigar_reference_length(cigar : Array[BamCigarElement]) -> Int32 {
  let mut len = 0
  for el in cigar {
    len = len + el.length * cigar_op_rlen(el.op)
  }
  len
}

///|
/// Sum read-consuming length including H clips (M/I/S/H/=/X). This is what
/// samtools `bam_cigar2qlen` uses; matches Pysam
/// `AlignedSegment.infer_read_length`.
pub fn cigar_read_length_with_hardclip(cigar : Array[BamCigarElement]) -> Int32 {
  let mut len = 0
  for el in cigar {
    match el.op {
      BamCigarOp::Match
      | BamCigarOp::Insert
      | BamCigarOp::SoftClip
      | BamCigarOp::HardClip
      | BamCigarOp::Equal
      | BamCigarOp::Diff => len = len + el.length
      _ => ()
    }
  }
  len
}

///|
/// Return a new CIGAR array with leading / trailing hard clips stripped.
/// Soft clips are preserved. Useful before writing CIGAR to text / BAM.
pub fn cigar_strip_hardclips(cigar : Array[BamCigarElement]) -> Array[BamCigarElement] {
  // Find first non-H element
  let mut start = 0
  while start < cigar.length() && cigar[start].op == BamCigarOp::HardClip {
    start = start + 1
  }
  // Find last non-H element
  let mut end = cigar.length()
  while end > start && cigar[end - 1].op == BamCigarOp::HardClip {
    end = end - 1
  }
  if start == 0 && end == cigar.length() {
    return cigar
  }
  let out : Array[BamCigarElement] = Array(capacity=end - start)
  let mut i = start
  while i < end {
    out.push(cigar[i])
    i = i + 1
  }
  out
}

///|
/// Return a new CIGAR array that collapses adjacent ops of the same type.
/// e.g. `[10M, 5M, 2I, 3D]` → `[15M, 2I, 3D]`.
pub fn cigar_merge_adjacent(cigar : Array[BamCigarElement]) -> Array[BamCigarElement] {
  let out : Array[BamCigarElement] = Array(capacity=cigar.length())
  if cigar.length() == 0 {
    return out
  }
  let mut cur = cigar[0]
  let mut i = 1
  while i < cigar.length() {
    if cigar[i].op == cur.op {
      cur = BamCigarElement::{ op: cur.op, length: cur.length + cigar[i].length, }
    } else {
      out.push(cur)
      cur = cigar[i]
    }
    i = i + 1
  }
  out.push(cur)
  out
}

///|
/// Produce an array of `(query_position, ref_position)` from a CIGAR, where
/// aligned positions (M/=/X) are emitted one pair per base, and insertion
/// positions are emitted at the next ref base. Deletions are skipped. Useful
/// for coordinate translation.
pub fn cigar_coordinate_pairs(
  cigar : Array[BamCigarElement],
  query_start : Int32,
  ref_start : Int32,
) -> Array[(Int32, Int32)] {
  let out : Array[(Int32, Int32)] = Array(capacity=0)
  let mut qpos = query_start
  let mut rpos = ref_start
  for el in cigar {
    match el.op {
      BamCigarOp::Match
      | BamCigarOp::Equal
      | BamCigarOp::Diff => {
        let mut k = 0
        while k < el.length {
          out.push((qpos, rpos))
          qpos = qpos + 1
          rpos = rpos + 1
          k = k + 1
        }
      }
      BamCigarOp::Insert => {
        let mut k = 0
        while k < el.length {
          out.push((qpos, rpos)) // insert at current ref base
          qpos = qpos + 1
          k = k + 1
        }
      }
      BamCigarOp::Delete | BamCigarOp::Skip => {
        rpos = rpos + el.length
      }
      BamCigarOp::SoftClip => {
        qpos = qpos + el.length
      }
      BamCigarOp::HardClip | BamCigarOp::Pad => ()
    }
  }
  out
}

///|
/// Advance `(qpos, rpos)` walking forward through the CIGAR until the
/// reference position `target_ref` is reached or passed. Returns the new
/// pair. Used by CIGAR-aware region lookups.
pub fn cigar_seek_ref(
  cigar : Array[BamCigarElement],
  qpos0 : Int32,
  rpos0 : Int32,
  target_ref : Int32,
) -> (Int32, Int32) {
  let mut qpos = qpos0
  let mut rpos = rpos0
  for el in cigar {
    let rdelta = el.length * cigar_op_rlen(el.op)
    if rpos + rdelta > target_ref {
      break
    }
    qpos = qpos + el.length * cigar_op_qlen(el.op)
    rpos = rpos + rdelta
  }
  (qpos, rpos)
}
