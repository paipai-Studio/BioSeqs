/// pysam_bam_mismatch.mbt — Mismatch analysis between read and reference
/// Inspired by pysam/samtools calmd + mismatch counting

///|
/// Mismatch information for a read.
pub struct PysamMismatchInfo {
  n_mismatch : Int32
  n_match : Int32
  n_insertion : Int32
  n_deletion : Int32
  n_soft_clip : Int32
  mismatch_positions : Array[Int32] // positions in read
  mismatch_bases : Array[(Int32, Int32, Int32)] // (read_pos, read_base, ref_base)
}

///|
/// Analyze mismatches between a read and reference.
pub fn pysam_analyze_mismatch(
  read_seq : String,
  ref_seq : String,
  cigar : String,
  read_pos : Int32, // 0-based position in reference
) -> PysamMismatchInfo {
  let mut n_mismatch = 0
  let mut n_match = 0
  let mut n_insertion = 0
  let mut n_deletion = 0
  let mut n_soft_clip = 0
  let mismatch_positions : Array[Int32] = Array(capacity=0)
  let mismatch_bases : Array[(Int32, Int32, Int32)] = Array(capacity=0)
  let mut ref_idx = read_pos
  let mut read_idx = 0
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      num = num * 10 + (c - 48)
    } else {
      match c {
        77 => // M
          for j = 0; j < num; j = j + 1 {
            let rb = if read_idx < read_seq.length() {
              read_seq[read_idx].to_int()
            } else {
              63
            }
            let fb = if ref_idx >= 0 && ref_idx < ref_seq.length() {
              ref_seq[ref_idx].to_int()
            } else {
              63
            }
            if rb != fb {
              n_mismatch = n_mismatch + 1
              mismatch_positions.push(read_idx)
              mismatch_bases.push((read_idx, rb, fb))
            } else {
              n_match = n_match + 1
            }
            read_idx = read_idx + 1
            ref_idx = ref_idx + 1
          }
        61 => { // =
          n_match = n_match + num
          read_idx = read_idx + num
          ref_idx = ref_idx + num
        }
        88 => // X
          for j = 0; j < num; j = j + 1 {
            let rb = if read_idx < read_seq.length() {
              read_seq[read_idx].to_int()
            } else {
              63
            }
            let fb = if ref_idx >= 0 && ref_idx < ref_seq.length() {
              ref_seq[ref_idx].to_int()
            } else {
              63
            }
            n_mismatch = n_mismatch + 1
            mismatch_positions.push(read_idx)
            mismatch_bases.push((read_idx, rb, fb))
            read_idx = read_idx + 1
            ref_idx = ref_idx + 1
          }
        73 => { // I
          n_insertion = n_insertion + num
          read_idx = read_idx + num
        }
        68 => { // D
          n_deletion = n_deletion + num
          ref_idx = ref_idx + num
        }
        83 => { // S
          n_soft_clip = n_soft_clip + num
          read_idx = read_idx + num
        }
        _ => ()
      }
      num = 0
    }
  }
  {
    n_mismatch,
    n_match,
    n_insertion,
    n_deletion,
    n_soft_clip,
    mismatch_positions,
    mismatch_bases,
  }
}

///|
/// Compute mismatch rate.
pub fn pysam_mismatch_rate(info : PysamMismatchInfo) -> Double {
  let aligned = info.n_match + info.n_mismatch
  if aligned == 0 {
    return 0.0
  }
  info.n_mismatch.to_double() / aligned.to_double()
}

///|
/// Compute edit distance (mismatches + insertions + deletions).
pub fn pysam_edit_distance(info : PysamMismatchInfo) -> Int32 {
  info.n_mismatch + info.n_insertion + info.n_deletion
}

///|
/// Generate MD tag from mismatch analysis.
pub fn pysam_generate_md_tag(
  read_seq : String,
  ref_seq : String,
  cigar : String,
  read_pos : Int32,
) -> String {
  let sb = StringBuilder()
  let mut ref_idx = read_pos
  let mut read_idx = 0
  let mut match_run = 0
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      num = num * 10 + (c - 48)
    } else {
      match c {
        77 | 61 | 88 =>
          for j = 0; j < num; j = j + 1 {
            let rb = if read_idx < read_seq.length() {
              read_seq[read_idx].to_int()
            } else {
              63
            }
            let fb = if ref_idx >= 0 && ref_idx < ref_seq.length() {
              ref_seq[ref_idx].to_int()
            } else {
              63
            }
            if rb == fb && (c == 77 || c == 61) {
              match_run = match_run + 1
            } else {
              if match_run > 0 {
                let _ = sb.write(match_run.to_string())
                match_run = 0
              }
              if c == 77 || c == 88 {
                let _ = sb.write_char(fb.unsafe_to_char())
              }
            }
            read_idx = read_idx + 1
            ref_idx = ref_idx + 1
          }
        68 => {
          if match_run > 0 {
            let _ = sb.write(match_run.to_string())
            match_run = 0
          }
          let _ = sb.write_char('^')
          for j = 0; j < num; j = j + 1 {
            if ref_idx >= 0 && ref_idx < ref_seq.length() {
              let _ = sb.write_char(ref_seq[ref_idx].to_int().unsafe_to_char())
            }
            ref_idx = ref_idx + 1
          }
        }
        73 => read_idx = read_idx + num
        83 => read_idx = read_idx + num
        _ => ()
      }
      num = 0
    }
  }
  if match_run > 0 {
    let _ = sb.write(match_run.to_string())
  }
  sb.to_string()
}
