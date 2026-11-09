// Pysam third-round: Post-alignment processing (items 16-20)
// Overlap clipping, BQSR, adapter trimming, quality sliding window, mate fix

///|
pub struct Pysam3ReadPair {
  read1_seq : String
  read1_qual : String
  read1_start : Int32
  read1_cigar : Array[(Int32, Int32)] // (op, len)
  read2_seq : String
  read2_qual : String
  read2_start : Int32
  read2_cigar : Array[(Int32, Int32)]
  read2_is_reverse : Bool
}

///|
pub fn pysam3_clip_overlap(pair : Pysam3ReadPair) -> Pysam3ReadPair {
  // Calculate read1 end position
  let mut r1_end = pair.read1_start
  for i in 0..<pair.read1_cigar.length() {
    let (op, length) = pair.read1_cigar[i]
    if op == 0 || op == 2 || op == 3 || op == 7 || op == 8 {
      r1_end = r1_end + length
    }
  }
  // Check if reads overlap
  if pair.read2_start >= r1_end {
    return pair
  }
  let overlap = r1_end - pair.read2_start
  if overlap <= 0 {
    return pair
  }
  // Clip overlap from 3' end of read1
  let r1_keep = pair.read1_seq.length() - overlap
  let new_r1_seq = if r1_keep > 0 {
    pair.read1_seq[0:r1_keep].to_owned()
  } else {
    ""
  }
  let new_r1_qual = if r1_keep > 0 {
    pair.read1_qual[0:r1_keep].to_owned()
  } else {
    ""
  }
  // Clip overlap from 5' end of read2
  let new_r2_seq = if overlap < pair.read2_seq.length() {
    pair.read2_seq[overlap:pair.read2_seq.length()].to_owned()
  } else {
    ""
  }
  let new_r2_qual = if overlap < pair.read2_qual.length() {
    pair.read2_qual[overlap:pair.read2_qual.length()].to_owned()
  } else {
    ""
  }
  {
    ..pair,
    read1_seq: new_r1_seq,
    read1_qual: new_r1_qual,
    read2_seq: new_r2_seq,
    read2_qual: new_r2_qual,
  }
}

// --- BQSR (Base Quality Score Recalibration) ---

///|
pub struct Pysam3BQSRTable {
  // Recalibration table: (qual_bin, observed_errors, total_observations)
  entries : Array[(Int32, Int32, Int32)]
}

///|
pub fn pysam3_bqsr_empty_table(num_bins : Int32) -> Pysam3BQSRTable {
  let entries : Array[(Int32, Int32, Int32)] = []
  for i in 0..<num_bins {
    entries.push((i, 0, 0))
  }
  { entries, }
}

///|
pub fn pysam3_bqsr_update(
  table : Pysam3BQSRTable,
  qual_bin : Int32,
  is_error : Bool,
) -> Pysam3BQSRTable {
  if qual_bin < 0 || qual_bin >= table.entries.length() {
    return table
  }
  let (bin, errors, total) = table.entries[qual_bin]
  let new_entry = (bin, if is_error { errors + 1 } else { errors }, total + 1)
  let new_entries = []
  for i in 0..<table.entries.length() {
    if i == qual_bin {
      new_entries.push(new_entry)
    } else {
      new_entries.push(table.entries[i])
    }
  }
  { entries: new_entries, }
}

///|
pub fn pysam3_bqsr_recalibrate(table : Pysam3BQSRTable, qual : Int32) -> Int32 {
  if qual < 0 || qual >= table.entries.length() {
    return qual
  }
  let (_, errors, total) = table.entries[qual]
  if total == 0 {
    return qual
  }
  let empirical_err = errors.to_double() / total.to_double()
  if empirical_err <= 0.0 {
    return qual
  }
  let empirical_qual = (-10.0 * @math.log10(empirical_err)).to_int()
  // Blend original and empirical
  (qual + empirical_qual) / 2
}

// --- Adapter trimming ---

///|
pub fn pysam3_trim_adapter(
  seq : String,
  qual : String,
  adapter : String,
  min_overlap : Int32,
) -> (String, String) {
  let seq_len = seq.length()
  let adp_len = adapter.length()
  let mut best_trim = 0
  let _trim_end = seq_len - min_overlap + 1
  for i in 0..<_trim_end {
    let max_cmp = if adp_len < seq_len - i { adp_len } else { seq_len - i }
    if max_cmp < min_overlap {
      break
    }
    let mut matches = 0
    for j in 0..<max_cmp {
      if seq[i + j] == adapter[j] {
        matches = matches + 1
      }
    }
    let mismatch_rate = 1.0 - matches.to_double() / max_cmp.to_double()
    if mismatch_rate < 0.1 && i > 0 {
      best_trim = i
      break
    }
  }
  if best_trim == 0 {
    return (seq, qual)
  }
  (seq[0:best_trim].to_owned(), qual[0:best_trim].to_owned())
}

// --- Quality sliding window trim (3' end) ---

///|
pub fn pysam3_trim_qual_window_3prime(
  seq : String,
  qual : String,
  window_size : Int32,
  qual_threshold : Int32,
) -> (String, String) {
  let n = seq.length()
  if n < window_size {
    return (seq, qual)
  }
  let quals = pysam3_phred33_to_quality(qual)
  let mut cut_pos = n
  let _win_end = n - window_size + 1
  for i in 0..<_win_end {
    let mut sum = 0
    for j in 0..<window_size {
      sum = sum + quals[i + j]
    }
    let avg = sum / window_size
    if avg < qual_threshold {
      cut_pos = i
      break
    }
  }
  if cut_pos == 0 {
    return ("", "")
  }
  (seq[0:cut_pos].to_owned(), qual[0:cut_pos].to_owned())
}

///|
pub fn pysam3_trim_qual_window_5prime(
  seq : String,
  qual : String,
  window_size : Int32,
  qual_threshold : Int32,
) -> (String, String) {
  let n = seq.length()
  if n < window_size {
    return (seq, qual)
  }
  let quals = pysam3_phred33_to_quality(qual)
  let mut cut_pos = 0
  let _n_ws = n - window_size + 1
  for i in 0..<_n_ws {
    let mut sum = 0
    for j in 0..<window_size {
      sum = sum + quals[i + j]
    }
    let avg = sum / window_size
    if avg >= qual_threshold {
      cut_pos = i
      break
    }
  }
  if cut_pos >= n {
    return ("", "")
  }
  (seq[cut_pos:n].to_owned(), qual[cut_pos:n].to_owned())
}

///|
pub fn pysam3_trim_qual_both_ends(
  seq : String,
  qual : String,
  window_size : Int32,
  qual_threshold : Int32,
) -> (String, String) {
  let (s1, q1) = pysam3_trim_qual_window_5prime(
    seq, qual, window_size, qual_threshold,
  )
  pysam3_trim_qual_window_3prime(s1, q1, window_size, qual_threshold)
}

// --- Mate coordinate fix ---

///|
pub struct Pysam3MateInfo {
  read_name : String
  flag : Int32
  ref_id : Int32
  pos : Int32
  mate_ref_id : Int32
  mate_pos : Int32
  insert_size : Int32
}

///|
pub fn pysam3_fix_mate(
  read1 : Pysam3MateInfo,
  read2 : Pysam3MateInfo,
) -> (Pysam3MateInfo, Pysam3MateInfo) {
  // Fix mate coordinates: read1's mate = read2's position, and vice versa
  let new_r1 = {
    ..read1,
    mate_ref_id: read2.ref_id,
    mate_pos: read2.pos,
    insert_size: if read1.ref_id == read2.ref_id {
      read2.pos - read1.pos
    } else {
      0
    },
  }
  let new_r2 = {
    ..read2,
    mate_ref_id: read1.ref_id,
    mate_pos: read1.pos,
    insert_size: if read1.ref_id == read2.ref_id {
      read1.pos - read2.pos
    } else {
      0
    },
  }
  (new_r1, new_r2)
}

///|
pub fn pysam3_fill_mate_info(
  reads : Array[Pysam3MateInfo],
) -> Array[Pysam3MateInfo] {
  // Build name -> index map for first and second reads
  let result = []
  for i in 0..<reads.length() {
    result.push(reads[i])
  }
  for i in 0..<reads.length() {
    let _j_start = i + 1
    for j in _j_start..<reads.length() {
      if reads[i].read_name == reads[j].read_name {
        let (fixed_i, fixed_j) = pysam3_fix_mate(reads[i], reads[j])
        result[i] = fixed_i
        result[j] = fixed_j
      }
    }
  }
  result
}

///|
pub fn pysam3_read_pair_dummy() -> Pysam3ReadPair {
  {
    read1_seq: "",
    read1_qual: "",
    read1_start: 0,
    read1_cigar: [],
    read2_seq: "",
    read2_qual: "",
    read2_start: 0,
    read2_cigar: [],
    read2_is_reverse: false,
  }
}

///|
pub fn pysam3_mate_info_dummy() -> Pysam3MateInfo {
  {
    read_name: "",
    flag: 0,
    ref_id: 0,
    pos: 0,
    mate_ref_id: 0,
    mate_pos: 0,
    insert_size: 0,
  }
}
