///|
/// htslib region parsing and indexing (region.c + regidx.c port).
/// Provides region list parsing ("chr1:100-200,chr2:300-400") and
/// a region index with efficient overlap queries.

///|
/// Parse an integer from a string, returning None on failure.
fn hts_parse_int(s : String) -> Int? {
  if s.length() == 0 {
    return None
  }
  let mut result = 0
  let mut i = 0
  let mut neg = false
  if s[0] == '-' {
    neg = true
    i = 1
  }
  if i >= s.length() {
    return None
  }
  while i < s.length() {
    let c = s[i]
    if c < '0' || c > '9' {
      return None
    }
    result = result * 10 + (c.to_int() - '0'.to_int())
    i = i + 1
  }
  Some(if neg { -result } else { result })
}

///|
/// A single genomic interval (0-based, half-open [beg, end)).
pub struct HtsRegion {
  chrom : String
  beg : Int // 0-based start
  end : Int // 0-based end (exclusive)
}

///|
/// Parse a single region string like "chr1:100-200" or "chr1".
/// Returns the parsed region, or None if the string is malformed.
/// Coordinates in the input are 1-based inclusive; internally stored
/// as 0-based half-open [beg, end).
pub fn hts_parse_region_str(s : String) -> HtsRegion? {
  // Find colon
  let mut colon_idx = -1
  for i = 0; i < s.length(); i = i + 1 {
    if s[i] == ':' {
      colon_idx = i
      break
    }
  }
  if colon_idx == -1 {
    return Some(HtsRegion::{ chrom: s, beg: 0, end: 0x7FFFFFFF })
  }
  let chrom = s[:colon_idx].to_owned()
  let range_str = s[colon_idx + 1:].to_owned()
  // Find dash in range_str
  let mut dash_idx = -1
  for i = 0; i < range_str.length(); i = i + 1 {
    if range_str[i] == '-' {
      dash_idx = i
      break
    }
  }
  if dash_idx == -1 {
    // chr1:100 → single position
    match hts_parse_int(range_str.to_string()) {
      None => None
      Some(p) => Some(HtsRegion::{ chrom, beg: p - 1, end: p })
    }
  } else {
    let start_str = range_str[:dash_idx].to_owned()
    let stop_str = range_str[dash_idx + 1:].to_owned()
    match (hts_parse_int(start_str), hts_parse_int(stop_str)) {
      (Some(s_val), Some(e_val)) => Some(HtsRegion::{ chrom, beg: s_val - 1, end: e_val })
      _ => None
    }
  }
}

///|
/// Parse a comma-separated list of regions.
/// Example: "chr1:100-200,chr2:300-400,chr3" → array of HtsRegion.
pub fn hts_parse_region_list(s : String) -> Array[HtsRegion] {
  let result : Array[HtsRegion] = []
  // Split by comma manually
  let parts : Array[String] = []
  let mut start = 0
  for i = 0; i < s.length(); i = i + 1 {
    if s[i] == ',' {
      parts.push(s[start:i].to_string())
      start = i + 1
    }
  }
  parts.push(s[start:].to_string())
  for part in parts {
    match hts_parse_region_str(part) {
      Some(r) => result.push(r)
      None => ()
    }
  }
  result
}

///|
/// Region index for efficient overlap queries.
/// Stores regions grouped by chromosome, sorted and merged.
pub struct HtsRegIdx {
  // chromosome name → list of merged intervals
  chroms : Map[String, Array[(Int, Int)]]
}

///|
/// Create an empty region index.
pub fn HtsRegIdx::new() -> HtsRegIdx {
  HtsRegIdx::{ chroms: Map::new() }
}

///|
/// Insert a single region into the index.
pub fn HtsRegIdx::insert(self : HtsRegIdx, chrom : String, beg : Int, end : Int) -> Unit {
  let intervals = match self.chroms.get(chrom) {
    Some(arr) => arr
    None => {
      let arr : Array[(Int, Int)] = []
      self.chroms[chrom] = arr
      arr
    }
  }
  intervals.push((beg, end))
}

///|
/// Insert a region from an HtsRegion value.
pub fn HtsRegIdx::insert_region(self : HtsRegIdx, r : HtsRegion) -> Unit {
  self.insert(r.chrom, r.beg, r.end)
}

///|
/// Insert regions from a comma-separated string.
pub fn HtsRegIdx::insert_list(self : HtsRegIdx, s : String) -> Unit {
  for r in hts_parse_region_list(s) {
    self.insert_region(r)
  }
}

///|
/// Insert a region from a BED-style line (chrom\tbeg\tend, 0-based).
pub fn HtsRegIdx::insert_bed_line(self : HtsRegIdx, line : String) -> Unit {
  // Split by tab
  let fields : Array[String] = []
  let mut start = 0
  for i = 0; i < line.length(); i = i + 1 {
    if line[i] == '\t' {
      fields.push(line[start:i].to_string())
      start = i + 1
    }
  }
  fields.push(line[start:].to_owned())
  if fields.length() >= 3 {
    let chrom = fields[0]
    let beg = hts_parse_int(fields[1])
    let end = hts_parse_int(fields[2])
    match (beg, end) {
      (Some(b), Some(e)) => self.insert(chrom, b, e)
      _ => ()
    }
  } else if fields.length() == 1 && fields[0].length() > 0 {
    self.insert(fields[0], 0, 0x7FFFFFFF)
  }
}

///|
/// Sort and merge overlapping or adjacent intervals for a chromosome.
fn HtsRegIdx::compact(self : HtsRegIdx, chrom : String) -> Unit {
  let intervals = match self.chroms.get(chrom) {
    Some(arr) => arr
    None => return
  }
  if intervals.length() <= 1 {
    return
  }
  // Sort by beg, then end
  intervals.sort_by(fn(a, b) {
    let (a0, a1) = a
    let (b0, b1) = b
    if a0 != b0 {
      a0 - b0
    } else {
      a1 - b1
    }
  })
  // Merge overlapping/adjacent
  let (i0b, i0e) = intervals[0]
  let merged : Array[(Int, Int)] = [(i0b, i0e)]
  for i = 1; i < intervals.length(); i = i + 1 {
    let (lbeg, lend) = merged[merged.length() - 1]
    let (cbeg, cend) = intervals[i]
    if lend < cbeg {
      merged.push((cbeg, cend))
    } else if lend < cend {
      merged[merged.length() - 1] = (lbeg, cend)
    }
  }
  self.chroms[chrom] = merged
}

///|
/// Compact all chromosomes (sort + merge intervals).
pub fn HtsRegIdx::compact_all(self : HtsRegIdx) -> Unit {
  let keys : Array[String] = []
  for k in self.chroms.keys() {
    keys.push(k)
  }
  for k in keys {
    self.compact(k)
  }
}

///|
/// Query whether a region [beg, end) overlaps any indexed region on chrom.
pub fn HtsRegIdx::overlap(
  self : HtsRegIdx,
  chrom : String,
  beg : Int,
  end : Int,
) -> Bool {
  let intervals = match self.chroms.get(chrom) {
    Some(arr) => arr
    None => return false
  }
  for i = 0; i < intervals.length(); i = i + 1 {
    let (ibeg, iend) = intervals[i]
    if beg < iend && end > ibeg {
      return true
    }
  }
  false
}

///|
/// Query overlap and return all matching intervals.
pub fn HtsRegIdx::overlapping(
  self : HtsRegIdx,
  chrom : String,
  beg : Int,
  end : Int,
) -> Array[(Int, Int)] {
  let result : Array[(Int, Int)] = []
  let intervals = match self.chroms.get(chrom) {
    Some(arr) => arr
    None => return result
  }
  for i = 0; i < intervals.length(); i = i + 1 {
    let (ibeg, iend) = intervals[i]
    if beg < iend && end > ibeg {
      result.push((ibeg, iend))
    }
  }
  result
}

///|
/// Get the number of regions for a given chromosome.
pub fn HtsRegIdx::seq_nregs(self : HtsRegIdx, chrom : String) -> Int {
  match self.chroms.get(chrom) {
    Some(arr) => arr.length()
    None => 0
  }
}

///|
/// Get the total number of regions across all chromosomes.
pub fn HtsRegIdx::nregs(self : HtsRegIdx) -> Int {
  let mut total = 0
  for v in self.chroms.values() {
    total = total + v.length()
  }
  total
}

///|
/// Get all chromosome names in the index.
pub fn HtsRegIdx::seq_names(self : HtsRegIdx) -> Array[String] {
  let names : Array[String] = []
  for k in self.chroms.keys() {
    names.push(k)
  }
  names
}

///|
/// Build a region index from a comma-separated region string.
pub fn hts_regidx_from_string(s : String) -> HtsRegIdx {
  let idx = HtsRegIdx::new()
  idx.insert_list(s)
  idx.compact_all()
  idx
}

///|
/// Build a region index from BED content (multi-line string).
pub fn hts_regidx_from_bed(content : String) -> HtsRegIdx {
  let idx = HtsRegIdx::new()
  // Split by newline
  let lines : Array[String] = []
  let mut start = 0
  for i = 0; i < content.length(); i = i + 1 {
    if content[i] == '\n' {
      lines.push(content[start:i].to_string())
      start = i + 1
    }
  }
  lines.push(content[start:].to_owned())
  for line in lines {
    let trimmed = line.trim()
    if trimmed.length() > 0 && trimmed[0] != '#' {
      idx.insert_bed_line(trimmed.to_string())
    }
  }
  idx.compact_all()
  idx
}
