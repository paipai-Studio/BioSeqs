/// pysam_bam_mate.mbt — Mate lookup, region parsing, downsample, multi-region
/// Items 1-4: mate查找, BAM降采样, 区域字符串解析, 多区域合并迭代

///|
/// Find mate of a read in SAM text by qname.
pub fn pysam_find_mate(
  sam_text : String,
  qname : String,
  flag : Int32,
) -> String? {
  let is_read1 = (flag & 64) != 0
  let lines = sam_text.split("\n").to_array()
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 2 {
      continue
    }
    let curr_qname = fields[0].to_owned()
    if curr_qname != qname {
      continue
    }
    let curr_flag = parse_int(fields[1].to_owned())
    let curr_is_read1 = (curr_flag & 64) != 0
    // Return the mate (opposite read number)
    if is_read1 && !curr_is_read1 {
      return Some(line)
    }
    if !is_read1 && curr_is_read1 {
      return Some(line)
    }
  }
  None
}

///|
/// Parsed region from string like "chr1:100-200" or "chr1:100-200(+)".
pub struct PysamParsedRegion {
  contig : String
  start : Int32 // 0-based, inclusive
  end : Int32 // exclusive
  strand : Int32 // 0=unknown, 1=forward, -1=reverse
}

///|
/// Parse a region string.
/// Formats: "chr1", "chr1:100-200", "chr1:100-200(+)", "chr1:100-200(-)"
pub fn pysam_parse_region_str(region : String) -> PysamParsedRegion? {
  // Find colon
  let colon_idx = region.find(":")
  match colon_idx {
    None => return Some({ contig: region, start: 0, end: -1, strand: 0, })
    Some(ci) => {
      let contig = region[0:ci].to_owned()
      let rest = region[ci + 1:].to_owned()
      // Find strand
      let mut strand = 0
      let mut range_str = rest
      if rest.has_suffix("(+)") {
        strand = 1
        range_str = rest[0:rest.length() - 3].to_owned()
      } else if rest.has_suffix("(-)") {
        strand = -1
        range_str = rest[0:rest.length() - 3].to_owned()
      }
      // Parse start-end
      let dash_idx = range_str.find("-")
      match dash_idx {
        None => return None
        Some(di) => {
          let start = parse_int(range_str[0:di].to_owned()) - 1
          let end = parse_int(range_str[di + 1:].to_owned())
          return Some({ contig, start, end, strand, })
        }
      }
    }
  }
}

///|
/// Merge overlapping regions and sort.
pub fn pysam_merge_regions(
  regions : Array[(String, Int32, Int32)],
) -> Array[(String, Int32, Int32)] {
  if regions.length() <= 1 {
    return regions
  }
  let sorted = regions.copy()
  sorted.sort_by(fn(a, b) {
    if a.0 < b.0 {
      -1
    } else if a.0 > b.0 {
      1
    } else if a.1 < b.1 {
      -1
    } else if a.1 > b.1 {
      1
    } else {
      0
    }
  })
  let merged : Array[(String, Int32, Int32)] = Array(capacity=0)
  merged.push(sorted[0])
  for i = 1; i < sorted.length(); i = i + 1 {
    let last = merged[merged.length() - 1]
    let curr = sorted[i]
    if curr.0 == last.0 && curr.1 <= last.2 {
      let new_end = if curr.2 > last.2 { curr.2 } else { last.2 }
      merged[merged.length() - 1] = (last.0, last.1, new_end)
    } else {
      merged.push(curr)
    }
  }
  merged
}

///|
/// Downsample SAM by fraction (0.0 to 1.0).
/// Uses simple hash-based sampling for reproducibility.
pub fn pysam_downsample(
  sam_text : String,
  fraction : Double,
  seed~ : Int32 = 42,
) -> String {
  if fraction >= 1.0 {
    return sam_text
  }
  let lines = sam_text.split("\n").to_array()
  let sb = StringBuilder()
  let threshold = (fraction * 2147483647.0).to_int()
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
    if fields.length() < 1 {
      let _ = sb.write(line + "\n")
      continue
    }
    let qname = fields[0].to_owned()
    // Simple hash
    let mut hash = seed
    for c = 0; c < qname.length(); c = c + 1 {
      hash = (hash * 31 + qname[c].to_int()) % 2147483647
    }
    if hash < 0 {
      hash = -hash
    }
    if hash < threshold {
      let _ = sb.write(line + "\n")
    }
  }
  sb.to_string()
}

///|
/// Downsample to exactly N reads.
pub fn pysam_downsample_n(
  sam_text : String,
  n : Int32,
  seed~ : Int32 = 42,
) -> String {
  ignore(seed)
  let lines = sam_text.split("\n").to_array()
  // Count records
  let mut total = 0
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() > 0 && line[0].to_int() != 64 {
      total = total + 1
    }
  }
  if total <= n {
    return sam_text
  }
  // Select every k-th record
  let k = total / n
  let sb = StringBuilder()
  let mut count = 0
  let mut kept = 0
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 {
      continue
    }
    if line[0].to_int() == 64 {
      let _ = sb.write(line + "\n")
      continue
    }
    if count % k == 0 && kept < n {
      let _ = sb.write(line + "\n")
      kept = kept + 1
    }
    count = count + 1
  }
  sb.to_string()
}

///|
/// Fetch from multiple merged regions.
pub fn pysam_fetch_multi_regions(
  sam_text : String,
  regions : Array[(String, Int32, Int32)],
) -> Array[String] {
  let merged = pysam_merge_regions(regions)
  let result : Array[String] = Array(capacity=0)
  let lines = sam_text.split("\n").to_array()
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 4 {
      continue
    }
    let rname = fields[2].to_owned()
    let pos = parse_int(fields[3].to_owned()) - 1
    for r = 0; r < merged.length(); r = r + 1 {
      let (contig, start, end) = merged[r]
      if rname == contig && pos >= start && pos < end {
        result.push(line)
        break
      }
    }
  }
  result
}
