///|
/// pyfaidx Sequence model: a named, possibly-coordinate-tagged subsequence of
/// a FASTA record, with Python-slice semantics (`Sequence.__getitem__`).
///
/// Ported from pyfaidx 0.9.0.4 `Sequence`. `start`/`end` are the record
/// coordinates (1-based or 0-based, tracked lazily across slices); when a
/// slice drops coordinate information (e.g. a step other than ±1) they become
/// `None`. The `-seq` operator is the reverse complement.

///|
pub struct PyfaidxSequence {
  name : String
  seq : String
  start : Int?
  end : Int?
  comp : Bool
}

///|
pub fn pyfaidx_seq_new(
  name : String,
  seq : String,
  start : Int?,
  end : Int?,
  comp : Bool,
) -> PyfaidxSequence {
  PyfaidxSequence::{ name, seq, start, end, comp }
}

///|
/// Python `slice.indices(length)` semantics. Returns normalized
/// `(start, stop, step)`; for a negative step the `stop` may be `-1` (the
/// "past the start" sentinel). `step` is never 0.
pub fn pyfaidx_slice_indices(
  start0 : Int?,
  stop0 : Int?,
  step0 : Int?,
  length : Int,
) -> (Int, Int, Int) raise PyfaidxError {
  let step = match step0 {
    Some(s) => s
    None => 1
  }
  if step == 0 {
    raise PyfaidxError("slice step cannot be zero")
  }
  let (defstart, defstop) = if step < 0 { (length - 1, -1) } else { (0, length) }
  // CPython only wraps/clamps an *explicit* index; the None defaults
  // (including the -1 sentinel for a negative step) are taken as-is.
  let start = match start0 {
    Some(s) => {
      let mut v = s
      if v < 0 {
        v = v + length
      }
      if v < 0 {
        v = if step < 0 { -1 } else { 0 }
      }
      if v >= length {
        v = if step < 0 { length - 1 } else { length }
      }
      v
    }
    None => defstart
  }
  let stop = match stop0 {
    Some(s) => {
      let mut v = s
      if v < 0 {
        v = v + length
      }
      if v < 0 {
        v = if step < 0 { -1 } else { 0 }
      }
      if v >= length {
        v = if step < 0 { length - 1 } else { length }
      }
      v
    }
    None => defstop
  }
  (start, stop, step)
}

///|
/// Extract a substring given normalized (start, stop, step) indices.
fn pyfaidx_string_slice(s : String, sstart : Int, sstop : Int, sstep : Int) -> String {
  let buf = StringBuilder()
  if sstep > 0 {
    let mut k = sstart
    while k < sstop {
      buf.write_char(s.unsafe_get(k).unsafe_to_char())
      k = k + sstep
    }
  } else {
    let mut k = sstart
    while k > sstop {
      buf.write_char(s.unsafe_get(k).unsafe_to_char())
      k = k + sstep
    }
  }
  buf.to_string()
}

///|
/// Coordinate correction factor (``-1`` for 1-based records, ``0`` for 0-based).
fn pyfaidx_correction_factor(seq : PyfaidxSequence) -> Int raise PyfaidxError {
  match (seq.start, seq.end) {
    (Some(s), Some(e)) => {
      let l = seq.seq.length()
      let span = (e - s).abs()
      if l == 0 {
        0
      } else if l == span + 1 {
        -1
      } else if l == span {
        0
      } else {
        raise PyfaidxError(
          "Coordinates imply a different length than Sequence.seq",
        )
      }
    }
    _ => 0
  }
}

///|
/// Integer indexing (`Sequence[i]`), with negative-index wrapping.
pub fn pyfaidx_seq_getitem(seq : PyfaidxSequence, index : Int) -> PyfaidxSequence {
  let mut n = index
  if n < 0 {
    n = n + seq.seq.length()
  }
  let ch = pyfaidx_string_slice(seq.seq, n, n + 1, 1)
  match seq.start {
    Some(s) =>
      if s != 0 {
        pyfaidx_seq_new(seq.name, ch, Some(s + n), Some(s + n), seq.comp)
      } else {
        pyfaidx_seq_new(seq.name, ch, None, None, seq.comp)
      }
    None => pyfaidx_seq_new(seq.name, ch, None, None, seq.comp)
  }
}

///|
/// Slice indexing (`Sequence[start:stop:step]`).
pub fn pyfaidx_seq_slice(
  seq : PyfaidxSequence,
  start0 : Int?,
  stop0 : Int?,
  step0 : Int?,
) -> PyfaidxSequence raise PyfaidxError {
  let (sstart, sstop, sstep) = pyfaidx_slice_indices(start0, stop0, step0, seq.seq.length())
  let sliced = pyfaidx_string_slice(seq.seq, sstart, sstop, sstep)
  match (seq.start, seq.end) {
    (Some(self_start), Some(self_end)) => {
      if sstep.abs() > 1 {
        pyfaidx_seq_new(seq.name, sliced, None, None, seq.comp)
      } else if sstep == -1 {
        let mut stop_adj = sstop
        if stop_adj == -1 {
          stop_adj = 0
        }
        pyfaidx_seq_new(
          seq.name,
          sliced,
          Some(self_end - stop_adj),
          Some(self_start + stop_adj),
          seq.comp,
        )
      } else {
        let cf = pyfaidx_correction_factor(seq)
        let (ns, ne) = if self_start > self_end {
          (self_start - sstart, self_start - sstop - cf)
        } else {
          (self_start + sstart, self_start + sstop + cf)
        }
        pyfaidx_seq_new(seq.name, sliced, Some(ns), Some(ne), seq.comp)
      }
    }
    _ => pyfaidx_seq_new(seq.name, sliced, None, None, seq.comp)
  }
}

///|
/// Complement of a sequence (toggles the complement flag).
pub fn pyfaidx_seq_complement(seq : PyfaidxSequence) -> PyfaidxSequence raise PyfaidxError {
  let comp = pyfaidx_complement(seq.seq)
  pyfaidx_seq_new(seq.name, comp, seq.start, seq.end, !seq.comp)
}

///|
/// Reverse of a sequence (``seq[::-1]``).
pub fn pyfaidx_seq_reverse(seq : PyfaidxSequence) -> PyfaidxSequence raise PyfaidxError {
  pyfaidx_seq_slice(seq, None, None, Some(-1))
}

///|
/// Reverse complement (``-seq``).
pub fn pyfaidx_seq_neg(seq : PyfaidxSequence) -> PyfaidxSequence raise PyfaidxError {
  pyfaidx_seq_complement(pyfaidx_seq_reverse(seq))
}

///|
/// `fancy_name`: `name[:start-end][ (complement)]`.
pub fn pyfaidx_seq_fancy_name(seq : PyfaidxSequence) -> String {
  let mut name = seq.name
  match (seq.start, seq.end) {
    (Some(s), Some(e)) => name = name + ":" + s.to_string() + "-" + e.to_string()
    _ => ()
  }
  if seq.comp {
    name = name + " (complement)"
  }
  name
}

///|
/// `orientation`: 1 forward, -1 reverse, None otherwise.
pub fn pyfaidx_seq_orientation(seq : PyfaidxSequence) -> Int? {
  match (seq.start, seq.end) {
    (Some(s), Some(e)) => {
      if s < e && not(seq.comp) {
        Some(1)
      } else if s > e && seq.comp {
        Some(-1)
      } else {
        None
      }
    }
    _ => None
  }
}