/// pysam_bam_coverage.mbt — Per-position coverage/depth calculation
/// Inspired by pysam pileup and samtools depth

///|
/// Coverage at a single position.
pub struct PysamCoveragePosition {
  contig : String
  pos : Int32 // 0-based
  depth : Int32 // number of reads covering this position
}

///|
/// Calculate per-position coverage from SAM records.
pub fn pysam_bam_coverage(
  sam_text : String,
  contig : String,
  start : Int32,
  end : Int32,
) -> Array[PysamCoveragePosition] {
  let lines = sam_text.split("\n").to_array()
  let depth : Array[Int32] = Array::make(end - start, 0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 6 {
      continue
    }
    let rname = fields[2].to_owned()
    if rname != contig {
      continue
    }
    let flag = parse_int(fields[1].to_owned())
    if (flag & 4) != 0 { // unmapped
      continue
    }
    let pos = parse_int(fields[3].to_owned()) - 1
    let cigar = fields[5].to_owned()
    // Get reference positions consumed by this read
    let ref_positions = pysam_cov_cigar_positions(pos, cigar)
    for j = 0; j < ref_positions.length(); j = j + 1 {
      let p = ref_positions[j]
      if p >= start && p < end {
        depth[p - start] = depth[p - start] + 1
      }
    }
  }
  let result : Array[PysamCoveragePosition] = Array(capacity=0)
  for p = start; p < end; p = p + 1 {
    result.push({ contig, pos: p, depth: depth[p - start], })
  }
  result
}

///|
/// Get reference positions from CIGAR.
fn pysam_cov_cigar_positions(start : Int32, cigar : String) -> Array[Int32] {
  let positions : Array[Int32] = Array(capacity=0)
  let mut ref_pos = start
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      num = num * 10 + (c - 48)
    } else {
      match c {
        77 | 61 | 88 => {
          // M, =, X: consume both query and reference
          for j = 0; j < num; j = j + 1 {
            positions.push(ref_pos + j)
          }
          ref_pos = ref_pos + num
        }
        68 =>
          // D: consume reference only
          ref_pos = ref_pos + num
        78 =>
          // N: reference skip
          ref_pos = ref_pos + num
        _ => () // I, S, H, P: don't consume reference
      }
      num = 0
    }
  }
  positions
}

///|
/// Calculate mean coverage for a region.
pub fn pysam_mean_coverage(
  sam_text : String,
  contig : String,
  start : Int32,
  end : Int32,
) -> Double {
  let cov = pysam_bam_coverage(sam_text, contig, start, end)
  if cov.length() == 0 {
    return 0.0
  }
  let mut total = 0
  for i = 0; i < cov.length(); i = i + 1 {
    total = total + cov[i].depth
  }
  total.to_double() / cov.length().to_double()
}

///|
/// Find positions with coverage above threshold.
pub fn pysam_coverage_above(
  sam_text : String,
  contig : String,
  start : Int32,
  end : Int32,
  threshold : Int32,
) -> Array[PysamCoveragePosition] {
  let cov = pysam_bam_coverage(sam_text, contig, start, end)
  let result : Array[PysamCoveragePosition] = Array(capacity=0)
  for i = 0; i < cov.length(); i = i + 1 {
    if cov[i].depth >= threshold {
      result.push(cov[i])
    }
  }
  result
}

///|
/// Calculate coverage histogram (depth -> count of positions).
pub fn pysam_coverage_histogram(
  sam_text : String,
  contig : String,
  start : Int32,
  end : Int32,
) -> Array[(Int32, Int32)] {
  let cov = pysam_bam_coverage(sam_text, contig, start, end)
  // Find max depth
  let mut max_depth = 0
  for i = 0; i < cov.length(); i = i + 1 {
    if cov[i].depth > max_depth {
      max_depth = cov[i].depth
    }
  }
  let hist : Array[Int32] = Array::make(max_depth + 1, 0)
  for i = 0; i < cov.length(); i = i + 1 {
    hist[cov[i].depth] = hist[cov[i].depth] + 1
  }
  let result : Array[(Int32, Int32)] = Array(capacity=0)
  for d = 0; d <= max_depth; d = d + 1 {
    result.push((d, hist[d]))
  }
  result
}
