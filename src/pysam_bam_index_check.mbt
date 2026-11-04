/// pysam_bam_index_check.mbt — BAM index validation and statistics
/// Inspired by pysam/samtools idxstats, index verification

///|
/// Index statistics per reference.
pub struct PysamIdxStats {
  contig : String
  length : Int32
  mapped : Int32
  unmapped : Int32
}

///|
/// Compute idxstats from SAM text.
pub fn pysam_idxstats(sam_text : String) -> Array[PysamIdxStats] {
  let lines = sam_text.split("\n").to_array()
  // Parse @SQ headers for contig lengths
  let contigs : Array[String] = Array(capacity=0)
  let lengths : Array[Int32] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() != 64 {
      break
    }
    // @SQ\tSN:chr1\tLN:248956422
    if line.has_prefix("@SQ") {
      let fields = line.split("\t").to_array()
      let mut sn = ""
      let mut ln = 0
      for f = 1; f < fields.length(); f = f + 1 {
        let tag = fields[f].to_owned()
        if tag.has_prefix("SN:") {
          sn = tag[3:].to_owned()
        } else if tag.has_prefix("LN:") {
          ln = parse_int(tag[3:].to_owned())
        }
      }
      if sn.length() > 0 {
        contigs.push(sn)
        lengths.push(ln)
      }
    }
  }
  // Count mapped/unmapped per contig
  let mapped : Array[Int32] = Array::make(contigs.length(), 0)
  let unmapped : Array[Int32] = Array::make(contigs.length(), 0)
  let mut total_unmapped = 0
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 3 {
      continue
    }
    let flag = parse_int(fields[1].to_owned())
    let rname = fields[2].to_owned()
    if (flag & 4) != 0 || rname == "*" {
      total_unmapped = total_unmapped + 1
      continue
    }
    // Find contig index
    for j = 0; j < contigs.length(); j = j + 1 {
      if contigs[j] == rname {
        mapped[j] = mapped[j] + 1
        break
      }
    }
  }
  // Build result
  let result : Array[PysamIdxStats] = Array(capacity=0)
  for j = 0; j < contigs.length(); j = j + 1 {
    result.push({
      contig: contigs[j],
      length: lengths[j],
      mapped: mapped[j],
      unmapped: unmapped[j],
    })
  }
  // Add unmapped entry
  result.push({ contig: "*", length: 0, mapped: 0, unmapped: total_unmapped, })
  result
}

///|
/// Check if index is consistent with SAM.
pub fn pysam_check_index(sam_text : String) -> Bool {
  // Simple check: verify @SQ headers exist
  let lines = sam_text.split("\n").to_array()
  let mut has_sq = false
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 {
      continue
    }
    if line[0].to_int() != 64 {
      break
    }
    if line.has_prefix("@SQ") {
      has_sq = true
      break
    }
  }
  has_sq
}

///|
/// Get reference sequence lengths from headers.
pub fn pysam_get_ref_lengths(sam_text : String) -> Array[(String, Int32)] {
  let lines = sam_text.split("\n").to_array()
  let result : Array[(String, Int32)] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() != 64 {
      break
    }
    if line.has_prefix("@SQ") {
      let fields = line.split("\t").to_array()
      let mut sn = ""
      let mut ln = 0
      for f = 1; f < fields.length(); f = f + 1 {
        let tag = fields[f].to_owned()
        if tag.has_prefix("SN:") {
          sn = tag[3:].to_owned()
        } else if tag.has_prefix("LN:") {
          ln = parse_int(tag[3:].to_owned())
        }
      }
      if sn.length() > 0 {
        result.push((sn, ln))
      }
    }
  }
  result
}

///|
/// Count total records in SAM.
pub fn pysam_count_records(sam_text : String) -> Int32 {
  let lines = sam_text.split("\n").to_array()
  let mut count = 0
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() > 0 && line[0].to_int() != 64 {
      count = count + 1
    }
  }
  count
}
