/// Parse MD tag to extract mismatch positions and types.
///
/// pysam equivalent: `pysam.AlignedSegment.get_aligned_pairs()` with
/// mismatch detection, or parsing the MD tag manually.
///
/// The MD tag format: runs of matching bases as numbers, followed by
/// mismatching reference bases (uppercase letters), with ^prefix for deletions.
/// Example: "10A5^C3" means 10 matches, 1 mismatch (A), 5 matches,
/// 1 deletion (C), 3 matches.

///|
/// A mismatch event extracted from the MD tag
pub struct MdMismatch {
  /// 0-based reference position
  ref_pos : Int32
  /// Type of mismatch: 'S' for substitution, 'D' for deletion
  kind : Char
  /// The reference base(s) at this position (for substitution: 1 base; for deletion: the deleted sequence)
  ref_base : String
  /// The query base at this position (for deletions, this is '-')
  query_base : Char
}

///|
pub fn MdMismatch::ref_pos(self : MdMismatch) -> Int32 {
  self.ref_pos
}

///|
pub fn MdMismatch::kind(self : MdMismatch) -> Char {
  self.kind
}

///|
pub fn MdMismatch::ref_base(self : MdMismatch) -> String {
  self.ref_base
}

///|
pub fn MdMismatch::query_base(self : MdMismatch) -> Char {
  self.query_base
}

///|
/// Parse an MD tag string and extract mismatch information.
///
/// alignment_start: 0-based reference position where the alignment begins
/// query_seq: the query sequence (used to determine the query base at mismatches)
/// cigar: the CIGAR array (used to map reference positions to query positions)
///
/// Returns an array of mismatches sorted by reference position.
pub fn parse_md_tag(
  md_tag : String,
  alignment_start : Int32,
  query_seq : String,
  cigar : FixedArray[BamCigarElement],
) -> Array[MdMismatch] {
  let mismatches : Array[MdMismatch] = []
  let mut ref_pos = alignment_start
  let mut _qpos : Int32 = 0
  // Build a reference-to-query position map from CIGAR (for M/=/X ops)
  let ref_to_query : Map[Int32, Int32] = Map([])
  let mut rq_qpos : Int32 = 0
  let mut rq_rpos = alignment_start
  for el in cigar {
    let oplen = el.length
    match el.op {
      BamCigarOp::Match | BamCigarOp::Equal | BamCigarOp::Diff =>
        for k = 0; k < oplen; k = k + 1 {
          ref_to_query[rq_rpos] = rq_qpos
          rq_qpos = rq_qpos + 1
          rq_rpos = rq_rpos + 1
        }
      BamCigarOp::Insert | BamCigarOp::SoftClip => rq_qpos = rq_qpos + oplen
      BamCigarOp::Delete | BamCigarOp::Skip => rq_rpos = rq_rpos + oplen
      BamCigarOp::HardClip | BamCigarOp::Pad => ()
    }
  }
  // Parse MD tag
  let mut i = 0
  let md_len = md_tag.length()
  while i < md_len {
    let ch = md_tag.unsafe_get(i).to_int()
    if ch >= 48 && ch <= 57 {
      // Digit: run of matches
      let num_str = StringBuilder()
      while i < md_len {
        let d = md_tag.unsafe_get(i).to_int()
        if d >= 48 && d <= 57 {
          num_str.write_char(Int32::unsafe_to_char(d))
          i = i + 1
        } else {
          break
        }
      }
      let match_run = parse_int(num_str.to_string())
      ref_pos = ref_pos + match_run
    } else if ch == 94 {
      // '^' deletion: followed by one or more reference bases
      i = i + 1
      let del_buf = StringBuilder()
      while i < md_len {
        let d = md_tag.unsafe_get(i).to_int()
        if (d >= 65 && d <= 90) || (d >= 97 && d <= 122) {
          del_buf.write_char(Int32::unsafe_to_char(d))
          i = i + 1
        } else {
          break
        }
      }
      let del_seq = del_buf.to_string()
      mismatches.push({
        ref_pos,
        kind: 'D',
        ref_base: del_seq,
        query_base: '-',
      })
      ref_pos = ref_pos + del_seq.length()
    } else {
      // Single mismatching base (A/C/G/T/N)
      let ref_base = Int32::unsafe_to_char(ch)
      let qb = match ref_to_query.get(ref_pos) {
        Some(qp) =>
          if qp >= 0 && qp < query_seq.length() {
            Int32::unsafe_to_char(query_seq.unsafe_get(qp).to_int())
          } else {
            'N'
          }
        None => 'N'
      }
      mismatches.push({
        ref_pos,
        kind: 'S',
        ref_base: Char::to_string(ref_base),
        query_base: qb,
      })
      ref_pos = ref_pos + 1
      i = i + 1
    }
  }
  mismatches
}

///|
/// Count total mismatches from an MD tag (substitutions + deletions).
pub fn count_md_mismatches(md_tag : String) -> Int32 {
  let mut count : Int32 = 0
  let mut i = 0
  let len = md_tag.length()
  while i < len {
    let ch = md_tag.unsafe_get(i).to_int()
    if ch >= 48 && ch <= 57 {
      // Skip number
      while i < len {
        let d = md_tag.unsafe_get(i).to_int()
        if d >= 48 && d <= 57 {
          i = i + 1
        } else {
          break
        }
      }
    } else if ch == 94 {
      // Deletion: count each deleted base
      i = i + 1
      while i < len {
        let d = md_tag.unsafe_get(i).to_int()
        if (d >= 65 && d <= 90) || (d >= 97 && d <= 122) {
          count = count + 1
          i = i + 1
        } else {
          break
        }
      }
    } else {
      // Substitution
      count = count + 1
      i = i + 1
    }
  }
  count
}

///|
/// Extract only substitution positions from MD tag (ignoring deletions).
/// Returns array of (ref_pos, ref_base) pairs.
pub fn md_substitution_positions(
  md_tag : String,
  alignment_start : Int32,
) -> Array[(Int32, Char)] {
  let result : Array[(Int32, Char)] = []
  let mut ref_pos = alignment_start
  let mut i = 0
  let len = md_tag.length()
  while i < len {
    let ch = md_tag.unsafe_get(i).to_int()
    if ch >= 48 && ch <= 57 {
      let num_str = StringBuilder()
      while i < len {
        let d = md_tag.unsafe_get(i).to_int()
        if d >= 48 && d <= 57 {
          num_str.write_char(Int32::unsafe_to_char(d))
          i = i + 1
        } else {
          break
        }
      }
      ref_pos = ref_pos + parse_int(num_str.to_string())
    } else if ch == 94 {
      i = i + 1
      while i < len {
        let d = md_tag.unsafe_get(i).to_int()
        if (d >= 65 && d <= 90) || (d >= 97 && d <= 122) {
          i = i + 1
        } else {
          break
        }
      }
    } else {
      result.push((ref_pos, Int32::unsafe_to_char(ch)))
      ref_pos = ref_pos + 1
      i = i + 1
    }
  }
  result
}
