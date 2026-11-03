/// Reconstruct aligned sequence with gaps from CIGAR and query sequence.
///
/// pysam equivalent: `pysam.AlignedSegment.get_aligned_sequence()`
/// Returns the query sequence with insertion bases and deletion gaps ('-')
/// inserted at the positions dictated by the CIGAR string, so the result
/// aligns position-by-position with the reference.

///|
/// Result of aligned sequence reconstruction
pub struct AlignedSeqResult {
  /// Query sequence with gaps for deletions (no soft-clip/hard-clip bases)
  query_aligned : String
  /// Reference-aligned positions: (query_index, ref_index) for each aligned column
  aligned_pairs : Array[(Int32, Int32?)]
  /// Length of the aligned sequence (including gaps)
  aligned_length : Int32
}

///|
pub fn AlignedSeqResult::query_aligned(self : AlignedSeqResult) -> String {
  self.query_aligned
}

///|
pub fn AlignedSeqResult::aligned_pairs(
  self : AlignedSeqResult,
) -> Array[(Int32, Int32?)] {
  self.aligned_pairs
}

///|
pub fn AlignedSeqResult::aligned_length(self : AlignedSeqResult) -> Int32 {
  self.aligned_length
}

///|
/// Reconstruct the aligned query sequence with gaps from CIGAR.
///
/// This inserts '-' for deletions and includes insertion bases, producing
/// a string that aligns column-by-column with the reference (excluding skips).
/// Soft-clipped and hard-clipped bases are excluded.
pub fn bam_get_aligned_sequence(rec : BamRecord) -> AlignedSeqResult {
  let buf = StringBuilder()
  let pairs : Array[(Int32, Int32?)] = []
  let mut qpos : Int32 = 0
  let mut rpos = rec.pos
  for el in rec.cigar {
    let oplen = el.length
    match el.op {
      BamCigarOp::Match | BamCigarOp::Equal | BamCigarOp::Diff =>
        for k = 0; k < oplen; k = k + 1 {
          if qpos < rec.seq.length() {
            buf.write_char(
              Int32::unsafe_to_char(rec.seq.unsafe_get(qpos).to_int()),
            )
            pairs.push((qpos, Some(rpos)))
          }
          qpos = qpos + 1
          rpos = rpos + 1
        }
      BamCigarOp::Insert =>
        for k = 0; k < oplen; k = k + 1 {
          if qpos < rec.seq.length() {
            buf.write_char(
              Int32::unsafe_to_char(rec.seq.unsafe_get(qpos).to_int()),
            )
            pairs.push((qpos, None))
          }
          qpos = qpos + 1
        }
      BamCigarOp::Delete =>
        for k = 0; k < oplen; k = k + 1 {
          buf.write_string("-")
          pairs.push((-1, Some(rpos)))
          rpos = rpos + 1
        }
      BamCigarOp::SoftClip | BamCigarOp::HardClip =>
        // Skip soft/hard clipped bases
        if el.op == BamCigarOp::SoftClip {
          qpos = qpos + oplen
        }
      BamCigarOp::Skip =>
        // N (skip) - advance reference but no gap in query
        rpos = rpos + oplen
      BamCigarOp::Pad => ()
    }
  }
  {
    query_aligned: buf.to_string(),
    aligned_pairs: pairs,
    aligned_length: buf.to_string().length(),
  }
}

///|
/// Get only the reference-aligned positions (no insertions, no deletions).
/// Returns array of (query_pos, ref_pos) for M/=/X operations.
pub fn bam_get_reference_aligned(rec : BamRecord) -> Array[(Int32, Int32)] {
  let result : Array[(Int32, Int32)] = []
  let mut qpos : Int32 = 0
  let mut rpos = rec.pos
  for el in rec.cigar {
    let oplen = el.length
    match el.op {
      BamCigarOp::Match | BamCigarOp::Equal | BamCigarOp::Diff =>
        for k = 0; k < oplen; k = k + 1 {
          result.push((qpos, rpos))
          qpos = qpos + 1
          rpos = rpos + 1
        }
      BamCigarOp::Insert | BamCigarOp::SoftClip => qpos = qpos + oplen
      BamCigarOp::Delete | BamCigarOp::Skip => rpos = rpos + oplen
      BamCigarOp::HardClip | BamCigarOp::Pad => ()
    }
  }
  result
}

///|
/// Get the query sequence with only soft-clipped bases removed.
/// (Convenience wrapper — same as rec.query_sequence() but clips soft-clipped ends.)
pub fn bam_get_query_no_clip(rec : BamRecord) -> String {
  let mut start : Int32 = 0
  let mut end : Int32 = rec.seq.length()
  for el in rec.cigar {
    match el.op {
      BamCigarOp::SoftClip =>
        if start == 0 {
          start = el.length
        } else {
          end = rec.seq.length() - el.length
        }
      _ => ()
    }
  }
  if start >= end {
    return ""
  }
  let buf = StringBuilder()
  for i = start; i < end; i = i + 1 {
    if i >= 0 && i < rec.seq.length() {
      buf.write_char(Int32::unsafe_to_char(rec.seq.unsafe_get(i).to_int()))
    }
  }
  buf.to_string()
}
