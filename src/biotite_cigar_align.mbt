/// biotite_cigar_align.mbt — CIGAR string <-> Alignment conversion
/// Inspired by biotite.sequence.align.CigarOp and cigar operations

///|
/// CIGAR operation types.
pub enum BiotiteCigarOp {
  BiotiteCigarMatch // M: match or mismatch (aligned)
  BiotiteCigarInsert // I: insertion in read (gap in reference)
  BiotiteCigarDelete // D: deletion in read (gap in reference)
  BiotiteCigarSoftClip // S: soft clip (present in read but not aligned)
  BiotiteCigarHardClip // H: hard clip (not present in read)
  BiotiteCigarN // N: skip (reference skip, gap in read)
  BiotiteCigarEqual // =: exact match
  BiotiteCigarDiff // X: mismatch
  BiotiteCigarPad // P: padding
}

///|
/// A CIGAR element: operation + length.
pub struct BiotiteCigarElement {
  op : BiotiteCigarOp
  length : Int32
}

///|
/// Convert CIGAR operation to character.
pub fn BiotiteCigarOp::to_char(self : BiotiteCigarOp) -> String {
  match self {
    BiotiteCigarMatch => "M"
    BiotiteCigarInsert => "I"
    BiotiteCigarDelete => "D"
    BiotiteCigarSoftClip => "S"
    BiotiteCigarHardClip => "H"
    BiotiteCigarN => "N"
    BiotiteCigarEqual => "="
    BiotiteCigarDiff => "X"
    BiotiteCigarPad => "P"
  }
}

///|
/// Parse a character to CIGAR operation.
pub fn BiotiteCigarOp::from_char(ch : String) -> BiotiteCigarOp? {
  match ch {
    "M" => Some(BiotiteCigarMatch)
    "I" => Some(BiotiteCigarInsert)
    "D" => Some(BiotiteCigarDelete)
    "S" => Some(BiotiteCigarSoftClip)
    "H" => Some(BiotiteCigarHardClip)
    "N" => Some(BiotiteCigarN)
    "=" => Some(BiotiteCigarEqual)
    "X" => Some(BiotiteCigarDiff)
    "P" => Some(BiotiteCigarPad)
    _ => None
  }
}

///|
/// CIGAR string representation.
pub struct BiotiteCigar {
  elements : Array[BiotiteCigarElement]
}

///|
/// Parse a CIGAR string (e.g. "3M1I2M2D4M") into elements.
pub fn BiotiteCigar::from_string(cigar_str : String) -> BiotiteCigar {
  let elements : Array[BiotiteCigarElement] = Array(capacity=0)
  let n = cigar_str.length()
  let mut i = 0
  while i < n {
    // Parse number
    let mut num = 0
    while i < n && cigar_str[i].to_int() >= 48 && cigar_str[i].to_int() <= 57 {
      num = num * 10 + (cigar_str[i].to_int() - 48)
      i = i + 1
    }
    // Parse operation
    if i < n {
      let ch = cigar_str[i:i + 1].to_owned()
      let op = BiotiteCigarOp::from_char(ch)
      if op is Some(o) {
        elements.push({ op: o, length: num, })
      }
      i = i + 1
    }
  }
  { elements, }
}

///|
/// Convert CIGAR to string representation.
pub fn BiotiteCigar::to_string(self : BiotiteCigar) -> String {
  let sb = StringBuilder()
  for i = 0; i < self.elements.length(); i = i + 1 {
    let elem = self.elements[i]
    let _ = sb.write(elem.length.to_string() + elem.op.to_char())
  }
  sb.to_string()
}

///|
/// Calculate the reference length consumed by the CIGAR.
/// M, D, N, =, X consume reference; I, S, H, P do not.
pub fn BiotiteCigar::reference_length(self : BiotiteCigar) -> Int32 {
  let mut total = 0
  for i = 0; i < self.elements.length(); i = i + 1 {
    let elem = self.elements[i]
    match elem.op {
      BiotiteCigarMatch
      | BiotiteCigarDelete
      | BiotiteCigarN
      | BiotiteCigarEqual
      | BiotiteCigarDiff => total = total + elem.length
      _ => ()
    }
  }
  total
}

///|
/// Calculate the read length consumed by the CIGAR.
/// M, I, S, =, X consume read; D, N, H, P do not.
pub fn BiotiteCigar::read_length(self : BiotiteCigar) -> Int32 {
  let mut total = 0
  for i = 0; i < self.elements.length(); i = i + 1 {
    let elem = self.elements[i]
    match elem.op {
      BiotiteCigarMatch
      | BiotiteCigarInsert
      | BiotiteCigarSoftClip
      | BiotiteCigarEqual
      | BiotiteCigarDiff => total = total + elem.length
      _ => ()
    }
  }
  total
}

///|
/// Calculate the alignment length (number of aligned columns).
/// Only M, =, X, I, D consume alignment columns.
pub fn BiotiteCigar::alignment_length(self : BiotiteCigar) -> Int32 {
  let mut total = 0
  for i = 0; i < self.elements.length(); i = i + 1 {
    let elem = self.elements[i]
    match elem.op {
      BiotiteCigarMatch
      | BiotiteCigarInsert
      | BiotiteCigarDelete
      | BiotiteCigarEqual
      | BiotiteCigarDiff => total = total + elem.length
      _ => ()
    }
  }
  total
}

///|
/// Convert an alignment (two aligned sequences with '-' for gaps) to CIGAR.
pub fn biotite_alignment_to_cigar(
  ref_seq : String,
  read_seq : String,
) -> BiotiteCigar {
  let elements : Array[BiotiteCigarElement] = Array(capacity=0)
  let n = ref_seq.length()
  let mut current_op : BiotiteCigarOp = BiotiteCigarMatch
  let mut current_len = 0
  for i = 0; i < n; i = i + 1 {
    let r = ref_seq[i:i + 1].to_owned()
    let q = read_seq[i:i + 1].to_owned()
    let op = if r == "-" {
      BiotiteCigarInsert
    } else if q == "-" {
      BiotiteCigarDelete
    } else if r == q {
      BiotiteCigarEqual
    } else {
      BiotiteCigarDiff
    }
    if i == 0 {
      current_op = op
      current_len = 1
    } else if op.to_char() == current_op.to_char() {
      current_len = current_len + 1
    } else {
      elements.push({ op: current_op, length: current_len, })
      current_op = op
      current_len = 1
    }
  }
  if current_len > 0 {
    elements.push({ op: current_op, length: current_len, })
  }
  { elements, }
}

///|
/// Apply CIGAR to a reference and read to produce aligned sequences.
/// Returns (aligned_ref, aligned_read) with '-' for gaps.
pub fn BiotiteCigar::to_alignment(
  self : BiotiteCigar,
  ref_seq : String,
  read_seq : String,
) -> (String, String) {
  let ref_sb = StringBuilder()
  let read_sb = StringBuilder()
  let mut ref_pos = 0
  let mut read_pos = 0
  for i = 0; i < self.elements.length(); i = i + 1 {
    let elem = self.elements[i]
    match elem.op {
      BiotiteCigarMatch | BiotiteCigarEqual | BiotiteCigarDiff =>
        for _j = 0; _j < elem.length; _j = _j + 1 {
          if ref_pos < ref_seq.length() {
            let _ = ref_sb.write(ref_seq[ref_pos:ref_pos + 1].to_owned())
          }
          if read_pos < read_seq.length() {
            let _ = read_sb.write(read_seq[read_pos:read_pos + 1].to_owned())
          }
          ref_pos = ref_pos + 1
          read_pos = read_pos + 1
        }
      BiotiteCigarInsert =>
        for _j = 0; _j < elem.length; _j = _j + 1 {
          let _ = ref_sb.write_char('-')
          if read_pos < read_seq.length() {
            let _ = read_sb.write(read_seq[read_pos:read_pos + 1].to_owned())
          }
          read_pos = read_pos + 1
        }
      BiotiteCigarDelete | BiotiteCigarN =>
        for _j = 0; _j < elem.length; _j = _j + 1 {
          if ref_pos < ref_seq.length() {
            let _ = ref_sb.write(ref_seq[ref_pos:ref_pos + 1].to_owned())
          }
          let _ = read_sb.write_char('-')
          ref_pos = ref_pos + 1
        }
      BiotiteCigarSoftClip =>
        // Soft clip: present in read but not in alignment
        read_pos = read_pos + elem.length
      BiotiteCigarHardClip =>
        // Hard clip: not present in read
        ()
      BiotiteCigarPad => ()
    }
  }
  (ref_sb.to_string(), read_sb.to_string())
}

///|
/// Use unused constructors.
fn _biotite_cigar_use_constructors() -> Unit {
  ignore(BiotiteCigarPad)
}
