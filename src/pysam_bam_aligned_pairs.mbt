/// pysam_bam_aligned_pairs.mbt — Aligned pairs, blocks, coordinate mapping, revcomp
/// Items 5-8: get_aligned_pairs, get_blocks, query_to_reference, reverse_complement

///|
/// Get aligned pairs (query_pos, ref_pos) from CIGAR.
/// skip_softclip=true excludes soft-clipped positions.
pub fn pysam_get_aligned_pairs(
  cigar : String,
  ref_start : Int32,
  skip_softclip~ : Bool = false,
) -> Array[(Int32, Int32)] {
  let pairs : Array[(Int32, Int32)] = Array(capacity=0)
  let mut ref_pos = ref_start
  let mut query_pos = 0
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      num = num * 10 + (c - 48)
    } else {
      match c {
        77 | 61 | 88 => { // M, =, X
          for j = 0; j < num; j = j + 1 {
            pairs.push((query_pos + j, ref_pos + j))
          }
          query_pos = query_pos + num
          ref_pos = ref_pos + num
        }
        73 => // I
          query_pos = query_pos + num
        68 => // D
          ref_pos = ref_pos + num
        78 => // N
          ref_pos = ref_pos + num
        83 => { // S
          if !skip_softclip {
            for j = 0; j < num; j = j + 1 {
              pairs.push((query_pos + j, -1))
            }
          }
          query_pos = query_pos + num
        }
        _ => ()
      }
      num = 0
    }
  }
  pairs
}

///|
/// Get contiguous alignment blocks [(ref_start, ref_end), ...].
/// Skips N (intron) and D (deletion) gaps.
pub fn pysam_get_blocks(
  cigar : String,
  ref_start : Int32,
) -> Array[(Int32, Int32)] {
  let blocks : Array[(Int32, Int32)] = Array(capacity=0)
  let mut ref_pos = ref_start
  let mut block_start = -1
  let mut in_block = false
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      num = num * 10 + (c - 48)
    } else {
      match c {
        77 | 61 | 88 => { // M, =, X
          if !in_block {
            block_start = ref_pos
            in_block = true
          }
          ref_pos = ref_pos + num
        }
        68 | 78 => { // D, N: gap
          if in_block {
            blocks.push((block_start, ref_pos))
            in_block = false
          }
          ref_pos = ref_pos + num
        }
        _ => ()
      }
      num = 0
    }
  }
  if in_block {
    blocks.push((block_start, ref_pos))
  }
  blocks
}

///|
/// Map a query position to reference position.
pub fn pysam_query_to_reference(
  cigar : String,
  ref_start : Int32,
  query_pos : Int32,
) -> Int32? {
  let mut ref_pos = ref_start
  let mut q_pos = 0
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      num = num * 10 + (c - 48)
    } else {
      match c {
        77 | 61 | 88 => { // M, =, X
          if query_pos >= q_pos && query_pos < q_pos + num {
            return Some(ref_pos + (query_pos - q_pos))
          }
          q_pos = q_pos + num
          ref_pos = ref_pos + num
        }
        73 => { // I
          if query_pos >= q_pos && query_pos < q_pos + num {
            return None // insertion has no ref position
          }
          q_pos = q_pos + num
        }
        68 => // D
          ref_pos = ref_pos + num
        78 => // N
          ref_pos = ref_pos + num
        83 => { // S
          if query_pos >= q_pos && query_pos < q_pos + num {
            return None // soft-clip has no ref position
          }
          q_pos = q_pos + num
        }
        _ => ()
      }
      num = 0
    }
  }
  None
}

///|
/// Map a reference position to query position.
pub fn pysam_reference_to_query(
  cigar : String,
  ref_start : Int32,
  ref_pos_target : Int32,
) -> Int32? {
  let mut ref_pos = ref_start
  let mut q_pos = 0
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      num = num * 10 + (c - 48)
    } else {
      match c {
        77 | 61 | 88 => { // M, =, X
          if ref_pos_target >= ref_pos && ref_pos_target < ref_pos + num {
            return Some(q_pos + (ref_pos_target - ref_pos))
          }
          q_pos = q_pos + num
          ref_pos = ref_pos + num
        }
        73 => // I
          q_pos = q_pos + num
        68 | 78 => // D, N
          ref_pos = ref_pos + num
        83 => // S
          q_pos = q_pos + num
        _ => ()
      }
      num = 0
    }
  }
  None
}

///|
/// Reverse complement a DNA sequence.
fn pysam_revcomp_seq(seq : String) -> String {
  let n = seq.length()
  let result : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    let comp = match seq[i].to_int() {
      65 => 84 // A->T
      84 => 65 // T->A
      71 => 67 // G->C
      67 => 71 // C->G
      97 => 116 // a->t
      116 => 97 // t->a
      103 => 99 // g->c
      99 => 103 // c->g
      x => x // N and others stay
    }
    result[n - 1 - i] = comp
  }
  let sb = StringBuilder()
  for i = 0; i < n; i = i + 1 {
    let _ = sb.write_char(result[i].unsafe_to_char())
  }
  sb.to_string()
}

///|
/// Reverse a quality string.
fn pysam_revqual(qual : String) -> String {
  let n = qual.length()
  let sb = StringBuilder()
  for i = n - 1; i >= 0; i = i - 1 {
    let _ = sb.write_char(qual[i].to_int().unsafe_to_char())
  }
  sb.to_string()
}

///|
/// Reverse CIGAR string (element order reversed).
fn pysam_rev_cigar(cigar : String) -> String {
  // Parse elements, reverse order
  let elements : Array[String] = Array(capacity=0)
  let mut num_str = ""
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      num_str = num_str + c.unsafe_to_char().to_string()
    } else {
      elements.push(num_str + c.unsafe_to_char().to_string())
      num_str = ""
    }
  }
  let sb = StringBuilder()
  for i = elements.length() - 1; i >= 0; i = i - 1 {
    let _ = sb.write(elements[i])
  }
  sb.to_string()
}

///|
/// Reverse complement a SAM read (for reverse-strand alignment).
pub fn pysam_revcomp_read(sam_line : String) -> String {
  let fields = sam_line.split("\t").to_array()
  if fields.length() < 11 {
    return sam_line
  }
  let flag = parse_int(fields[1].to_owned())
  let is_reverse = (flag & 16) != 0
  if !is_reverse {
    return sam_line // only revcomp reverse-strand reads
  }
  let seq = fields[9].to_owned()
  let qual = fields[10].to_owned()
  let cigar = fields[5].to_owned()
  // Reverse complement sequence
  fields[9] = pysam_revcomp_seq(seq)
  // Reverse quality
  if qual != "*" {
    fields[10] = pysam_revqual(qual)
  }
  // Reverse CIGAR
  fields[5] = pysam_rev_cigar(cigar)
  // Rebuild line
  let sb = StringBuilder()
  for i = 0; i < fields.length(); i = i + 1 {
    if i > 0 {
      let _ = sb.write_char('\t')
    }
    let _ = sb.write(fields[i])
  }
  sb.to_string()
}
