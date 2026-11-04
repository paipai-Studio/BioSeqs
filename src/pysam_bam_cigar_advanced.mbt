/// pysam_bam_cigar_advanced.mbt — Advanced CIGAR operations
/// Inspired by pysam AlignedSegment.cigar manipulation

///|
/// CIGAR operation type.
pub enum PysamCigarOp {
  M // match/mismatch
  I // insertion to reference
  D // deletion from reference
  N // reference skip
  S // soft clip
  H // hard clip
  P // padding
  Eq // sequence match
  X // sequence mismatch
}

///|
/// A single CIGAR element.
pub struct PysamCigarElement {
  length : Int32
  op : PysamCigarOp
}

///|
/// Parse CIGAR string into elements.
pub fn pysam_parse_cigar(cigar : String) -> Array[PysamCigarElement] {
  let result : Array[PysamCigarElement] = Array(capacity=0)
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      num = num * 10 + (c - 48)
    } else {
      let op = match c {
        77 => PysamCigarOp::M
        73 => PysamCigarOp::I
        68 => PysamCigarOp::D
        78 => PysamCigarOp::N
        83 => PysamCigarOp::S
        72 => PysamCigarOp::H
        80 => PysamCigarOp::P
        61 => PysamCigarOp::Eq
        88 => PysamCigarOp::X
        _ => PysamCigarOp::M
      }
      result.push({ length: num, op, })
      num = 0
    }
  }
  result
}

///|
/// Format CIGAR elements to string.
pub fn pysam_cigar_to_string(elements : Array[PysamCigarElement]) -> String {
  let sb = StringBuilder()
  for i = 0; i < elements.length(); i = i + 1 {
    let e = elements[i]
    let _ = sb.write(e.length.to_string())
    let c = match e.op {
      PysamCigarOp::M => "M"
      PysamCigarOp::I => "I"
      PysamCigarOp::D => "D"
      PysamCigarOp::N => "N"
      PysamCigarOp::S => "S"
      PysamCigarOp::H => "H"
      PysamCigarOp::P => "P"
      PysamCigarOp::Eq => "="
      PysamCigarOp::X => "X"
    }
    let _ = sb.write(c)
  }
  sb.to_string()
}

///|
/// Get query (read) consuming length from CIGAR.
pub fn pysam_cigar_query_len(cigar : String) -> Int32 {
  let elements = pysam_parse_cigar(cigar)
  let mut total = 0
  for i = 0; i < elements.length(); i = i + 1 {
    let e = elements[i]
    match e.op {
      PysamCigarOp::M
      | PysamCigarOp::I
      | PysamCigarOp::S
      | PysamCigarOp::Eq
      | PysamCigarOp::X => total = total + e.length
      _ => ()
    }
  }
  total
}

///|
/// Get reference consuming length from CIGAR.
pub fn pysam_cigar_ref_len_adv(cigar : String) -> Int32 {
  let elements = pysam_parse_cigar(cigar)
  let mut total = 0
  for i = 0; i < elements.length(); i = i + 1 {
    let e = elements[i]
    match e.op {
      PysamCigarOp::M
      | PysamCigarOp::D
      | PysamCigarOp::N
      | PysamCigarOp::Eq
      | PysamCigarOp::X => total = total + e.length
      _ => ()
    }
  }
  total
}

///|
/// Get aligned length (M + = + X only).
pub fn pysam_cigar_aligned_len(cigar : String) -> Int32 {
  let elements = pysam_parse_cigar(cigar)
  let mut total = 0
  for i = 0; i < elements.length(); i = i + 1 {
    let e = elements[i]
    match e.op {
      PysamCigarOp::M | PysamCigarOp::Eq | PysamCigarOp::X =>
        total = total + e.length
      _ => ()
    }
  }
  total
}

///|
/// Count clip bases (S + H).
pub fn pysam_cigar_clip_len(cigar : String) -> Int32 {
  let elements = pysam_parse_cigar(cigar)
  let mut total = 0
  for i = 0; i < elements.length(); i = i + 1 {
    let e = elements[i]
    match e.op {
      PysamCigarOp::S | PysamCigarOp::H => total = total + e.length
      _ => ()
    }
  }
  total
}

///|
/// Invert CIGAR (for reverse complement alignment).
pub fn pysam_invert_cigar(cigar : String) -> String {
  let elements = pysam_parse_cigar(cigar)
  // Reverse order, keep operations same
  let reversed : Array[PysamCigarElement] = Array(capacity=0)
  for i = elements.length() - 1; i >= 0; i = i - 1 {
    reversed.push(elements[i])
  }
  pysam_cigar_to_string(reversed)
}

///|
/// Soft-clip left side of CIGAR.
pub fn pysam_softclip_left(cigar : String, clip_len : Int32) -> String {
  let elements = pysam_parse_cigar(cigar)
  let sb = StringBuilder()
  let _ = sb.write(clip_len.to_string() + "S")
  for i = 0; i < elements.length(); i = i + 1 {
    let e = elements[i]
    let c = match e.op {
      PysamCigarOp::M => "M"
      PysamCigarOp::I => "I"
      PysamCigarOp::D => "D"
      PysamCigarOp::N => "N"
      PysamCigarOp::S => "S"
      PysamCigarOp::H => "H"
      PysamCigarOp::P => "P"
      PysamCigarOp::Eq => "="
      PysamCigarOp::X => "X"
    }
    let _ = sb.write(e.length.to_string() + c)
  }
  sb.to_string()
}

///|
/// Count operations of a specific type.
pub fn pysam_cigar_count_op(cigar : String, op : PysamCigarOp) -> Int32 {
  let elements = pysam_parse_cigar(cigar)
  let mut total = 0
  for i = 0; i < elements.length(); i = i + 1 {
    // Compare by checking the string representation
    let matches = match (elements[i].op, op) {
      (PysamCigarOp::M, PysamCigarOp::M) => true
      (PysamCigarOp::I, PysamCigarOp::I) => true
      (PysamCigarOp::D, PysamCigarOp::D) => true
      (PysamCigarOp::N, PysamCigarOp::N) => true
      (PysamCigarOp::S, PysamCigarOp::S) => true
      (PysamCigarOp::H, PysamCigarOp::H) => true
      (PysamCigarOp::P, PysamCigarOp::P) => true
      (PysamCigarOp::Eq, PysamCigarOp::Eq) => true
      (PysamCigarOp::X, PysamCigarOp::X) => true
      _ => false
    }
    if matches {
      total = total + elements[i].length
    }
  }
  total
}
