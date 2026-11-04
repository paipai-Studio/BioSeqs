/// pysam_bam_baq.mbt — Base Alignment Quality (BAQ)
/// Inspired by pysam/samtools baq

///|
/// Compute BAQ scores for a read.
/// BAQ adjusts base qualities based on alignment to reference.
/// Returns adjusted quality scores.
pub fn pysam_baq_compute(
  read_seq : String,
  ref_seq : String,
  read_qual : Array[Int32],
  read_pos : Int32, // offset in reference
) -> Array[Int32] {
  let n = read_seq.length()
  let result : Array[Int32] = Array::make(n, 0)
  if n == 0 || read_qual.length() < n {
    return result
  }
  // Simple BAQ: for each base, compute a quality adjustment
  // based on whether the base matches the reference
  for i = 0; i < n; i = i + 1 {
    let read_base = read_seq[i].to_int()
    let ref_idx = read_pos + i
    let ref_base = if ref_idx >= 0 && ref_idx < ref_seq.length() {
      ref_seq[ref_idx].to_int()
    } else {
      -1
    }
    if ref_base == read_base {
      // Match: keep original quality
      result[i] = read_qual[i]
    } else {
      // Mismatch: reduce quality
      result[i] = if read_qual[i] > 10 { read_qual[i] - 10 } else { 0 }
    }
  }
  result
}

///|
/// Apply BAQ to SAM records.
/// Adds BQ tag with adjusted qualities.
pub fn pysam_baq_apply(sam_text : String, ref_seq : String) -> String {
  let lines = sam_text.split("\n").to_array()
  let sb = StringBuilder()
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 {
      continue
    }
    if line[0].to_int() == 64 {
      let _ = sb.write(line + "\n")
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 11 {
      let _ = sb.write(line + "\n")
      continue
    }
    let pos = parse_int(fields[3].to_owned()) - 1
    let seq = fields[9].to_owned()
    let qual_str = fields[10].to_owned()
    // Parse quality scores
    let quals : Array[Int32] = Array::make(qual_str.length(), 0)
    for q = 0; q < qual_str.length(); q = q + 1 {
      quals[q] = qual_str[q].to_int() - 33
    }
    // Compute BAQ
    let baq = pysam_baq_compute(seq, ref_seq, quals, pos)
    // Create BQ tag
    let bq_sb = StringBuilder()
    for q = 0; q < baq.length(); q = q + 1 {
      let _ = bq_sb.write_char((baq[q] + 33).unsafe_to_char())
    }
    // Add BQ tag
    let new_line = line + "\tBQ:Z:" + bq_sb.to_string()
    let _ = sb.write(new_line + "\n")
  }
  sb.to_string()
}

///|
/// Calculate mean base quality.
pub fn pysam_mean_base_quality(qual_str : String) -> Double {
  if qual_str.length() == 0 {
    return 0.0
  }
  let mut total = 0
  for i = 0; i < qual_str.length(); i = i + 1 {
    total = total + (qual_str[i].to_int() - 33)
  }
  total.to_double() / qual_str.length().to_double()
}

///|
/// Quality score distribution.
pub fn pysam_quality_distribution(qual_str : String) -> Array[(Int32, Int32)] {
  let max_q = 60
  let hist : Array[Int32] = Array::make(max_q + 1, 0)
  for i = 0; i < qual_str.length(); i = i + 1 {
    let q = qual_str[i].to_int() - 33
    if q >= 0 && q <= max_q {
      hist[q] = hist[q] + 1
    }
  }
  let result : Array[(Int32, Int32)] = Array(capacity=0)
  for q = 0; q <= max_q; q = q + 1 {
    if hist[q] > 0 {
      result.push((q, hist[q]))
    }
  }
  result
}
