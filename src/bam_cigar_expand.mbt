/// Expand CIGAR operations to per-base operations.
///
/// pysam equivalent: iterating `pysam.AlignedSegment.cigar` and expanding
/// each (op, length) into a flat sequence of single-base operations.
///
/// For example, "3M1I2D" expands to [M, M, M, I, D, D].

///|
/// A single-base CIGAR operation (expanded form)
pub struct ExpandedCigarOp {
  op : BamCigarOp
  query_pos : Int32
  ref_pos : Int32
}

///|
pub fn ExpandedCigarOp::op(self : ExpandedCigarOp) -> BamCigarOp {
  self.op
}

///|
pub fn ExpandedCigarOp::query_pos(self : ExpandedCigarOp) -> Int32 {
  self.query_pos
}

///|
pub fn ExpandedCigarOp::ref_pos(self : ExpandedCigarOp) -> Int32 {
  self.ref_pos
}

///|
/// Expand a CIGAR array into per-base operations with position tracking.
///
/// query_pos starts at 0 (after soft-clip handling depends on include_clip).
/// ref_pos starts at rec.pos.
/// SoftClip and HardClip operations are included if include_clip is true.
pub fn bam_expand_cigar(
  rec : BamRecord,
  include_clip~ : Bool = false,
) -> Array[ExpandedCigarOp] {
  let result : Array[ExpandedCigarOp] = []
  let mut qpos : Int32 = 0
  let mut rpos = rec.pos
  for el in rec.cigar {
    let oplen = el.length
    match el.op {
      BamCigarOp::Match | BamCigarOp::Equal | BamCigarOp::Diff =>
        for k = 0; k < oplen; k = k + 1 {
          result.push({ op: el.op, query_pos: qpos, ref_pos: rpos, })
          qpos = qpos + 1
          rpos = rpos + 1
        }
      BamCigarOp::Insert =>
        for k = 0; k < oplen; k = k + 1 {
          result.push({ op: el.op, query_pos: qpos, ref_pos: rpos, })
          qpos = qpos + 1
        }
      BamCigarOp::Delete =>
        for k = 0; k < oplen; k = k + 1 {
          result.push({ op: el.op, query_pos: -1, ref_pos: rpos, })
          rpos = rpos + 1
        }
      BamCigarOp::Skip =>
        for k = 0; k < oplen; k = k + 1 {
          result.push({ op: el.op, query_pos: -1, ref_pos: rpos, })
          rpos = rpos + 1
        }
      BamCigarOp::SoftClip =>
        if include_clip {
          for k = 0; k < oplen; k = k + 1 {
            result.push({ op: el.op, query_pos: qpos, ref_pos: -1, })
            qpos = qpos + 1
          }
        } else {
          qpos = qpos + oplen
        }
      BamCigarOp::HardClip =>
        if include_clip {
          for k = 0; k < oplen; k = k + 1 {
            result.push({ op: el.op, query_pos: -1, ref_pos: -1, })
          }
        }
      BamCigarOp::Pad =>
        for k = 0; k < oplen; k = k + 1 {
          result.push({ op: el.op, query_pos: qpos, ref_pos: rpos, })
        }
    }
  }
  result
}

///|
/// Expand CIGAR to a string of single-char operations.
/// e.g., "3M1I2D" -> "MMMIID"
pub fn bam_expand_cigar_string(rec : BamRecord) -> String {
  let buf = StringBuilder()
  for el in rec.cigar {
    let oplen = el.length
    let c = match el.op {
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
    for k = 0; k < oplen; k = k + 1 {
      buf.write_char(c)
    }
  }
  buf.to_string()
}

///|
/// Count operations by type in a CIGAR array.
/// Returns a map from operation char to count.
pub fn bam_cigar_op_counts(rec : BamRecord) -> Map[Char, Int32] {
  let counts : Map[Char, Int32] = Map([])
  for el in rec.cigar {
    let c = match el.op {
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
    let cur = match counts.get(c) {
      Some(v) => v
      None => 0
    }
    counts[c] = cur + el.length
  }
  counts
}
