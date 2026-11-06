// Bioframe third-round: format conversions (items 26-29)
// bed_to_gff, gff_to_bed, bedgraph_to_bed, ucsc_track_parse

///|
pub fn bf3_bed_to_gff(
  bed : Array[(String, Int32, Int32, String, Double, String)],
) -> Array[String] {
  // Convert BED6 to GFF3 lines
  // BED: (chrom, start, end, name, score, strand)
  // GFF3: seqid, source, type, start(1-based), end, score, strand, phase, attributes
  let lines : Array[String] = []
  for i in 0..<bed.length() {
    let (chrom, start, end_, name, score, strand) = bed[i]
    let gff_start = start + 1 // BED is 0-based, GFF is 1-based
    let strand_ = if strand == "" { "." } else { strand }
    let attrs = "ID=" + name + ";Name=" + name
    lines.push(
      chrom +
      "\t" +
      "BED" +
      "\t" +
      "region" +
      "\t" +
      gff_start.to_string() +
      "\t" +
      end_.to_string() +
      "\t" +
      score.to_string() +
      "\t" +
      strand_ +
      "\t" +
      "." +
      "\t" +
      attrs,
    )
  }
  lines
}

///|
pub fn bf3_gff_to_bed(
  gff_lines : Array[String],
) -> Array[(String, Int32, Int32, String, Double, String)] {
  // Convert GFF3 lines to BED6
  // GFF3: seqid, source, type, start(1-based), end, score, strand, phase, attributes
  let result : Array[(String, Int32, Int32, String, Double, String)] = []
  for i in 0..<gff_lines.length() {
    let line = gff_lines[i]
    if line.length() == 0 || line[0].to_int() == 35 { // skip empty and comments (#)
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() >= 9 {
      let chrom = fields[0].to_owned()
      let start = parse_int(fields[3].to_owned()) - 1 // GFF 1-based → BED 0-based
      let end_ = parse_int(fields[4].to_owned())
      let score = parse_double(fields[5].to_owned()).unwrap_or(0.0)
      let strand = fields[6].to_owned()
      // Extract name from attributes (ID=xxx;Name=yyy)
      let attrs = fields[8].to_owned()
      let name = {
        let parts = attrs.split(";").to_array()
        let mut n = "region"
        for p in 0..<parts.length() {
          let part = parts[p].to_owned()
          if part.has_prefix("Name=") {
            n = part[5:].to_owned()
            break
          }
          if part.has_prefix("ID=") && n == "region" {
            n = part[3:].to_owned()
          }
        }
        n
      }
      result.push((chrom, start, end_, name, score, strand))
    }
  }
  result
}

///|
pub fn bf3_bedgraph_to_bed(
  bedgraph : Array[(String, Int32, Int32, Double)],
  threshold : Double,
) -> Array[(String, Int32, Int32, Double)] {
  // Convert bedGraph to BED by thresholding: keep regions where value >= threshold
  let result : Array[(String, Int32, Int32, Double)] = []
  for i in 0..<bedgraph.length() {
    let (chrom, start, end_, value) = bedgraph[i]
    if value >= threshold {
      result.push((chrom, start, end_, value))
    }
  }
  result
}

///|
pub fn bf3_bedgraph_to_bed_regions(
  bedgraph : Array[(String, Int32, Int32, Double)],
  threshold : Double,
) -> Array[(String, Int32, Int32)] {
  // Threshold bedGraph and merge adjacent regions above threshold into BED
  let above : Array[(String, Int32, Int32)] = []
  for i in 0..<bedgraph.length() {
    let (chrom, start, end_, value) = bedgraph[i]
    if value >= threshold {
      above.push((chrom, start, end_))
    }
  }
  // Merge adjacent/overlapping
  if above.is_empty() {
    return []
  }
  let sorted = above.copy()
  sorted.sort_by(fn(a, b) {
    let (ca, sa, _) = a
    let (cb, sb, _) = b
    if ca < cb {
      -1
    } else if ca > cb {
      1
    } else if sa < sb {
      -1
    } else if sa > sb {
      1
    } else {
      0
    }
  })
  let merged : Array[(String, Int32, Int32)] = []
  let (m_chrom, m_start, m_end) = sorted[0]
  let mut cur_chrom = m_chrom
  let mut cur_start = m_start
  let mut cur_end = m_end
  for i in 1..<sorted.length() {
    let (chrom, start, end_) = sorted[i]
    if chrom == cur_chrom && start <= cur_end {
      if end_ > cur_end {
        cur_end = end_
      }
    } else {
      merged.push((cur_chrom, cur_start, cur_end))
      cur_chrom = chrom
      cur_start = start
      cur_end = end_
    }
  }
  merged.push((cur_chrom, cur_start, cur_end))
  merged
}

///|
pub struct Bf3UcscTrack {
  track_type : String
  name : String
  description : String
  visibility : String
  color : String
  url : String
  extra_params : Map[String, String]
}

///|
pub fn bf3_ucsc_track_parse(track_line : String) -> Bf3UcscTrack? {
  // Parse a UCSC track declaration line
  // Example: track name="myTrack" description="My track" visibility=2 color=255,0,0
  if not(track_line.has_prefix("track")) {
    return None
  }
  let params : Map[String, String] = Map([])
  let mut track_type = ""
  let mut name = ""
  let mut description = ""
  let mut visibility = ""
  let mut color = ""
  let mut url = ""
  // Simple parser: split by space, handle key=value and key="quoted value"
  let rest = track_line[6:].to_owned() // skip "track "
  let mut i = 0
  let chars = rest
  while i < chars.length() {
    // Skip whitespace
    while i < chars.length() && chars[i].to_int() == 32 {
      i = i + 1
    }
    if i >= chars.length() {
      break
    }
    // Read key
    let key_start = i
    while i < chars.length() &&
          chars[i].to_int() != 61 &&
          chars[i].to_int() != 32 {
      i = i + 1
    }
    let key = chars[key_start:i].to_owned()
    if i >= chars.length() || chars[i].to_int() != 61 {
      // No value, skip
      continue
    }
    i = i + 1 // skip '='
    // Read value
    let mut value = ""
    if i < chars.length() && chars[i].to_int() == 34 {
      // Quoted value
      i = i + 1 // skip opening quote
      let val_start = i
      while i < chars.length() && chars[i].to_int() != 34 {
        i = i + 1
      }
      value = chars[val_start:i].to_owned()
      if i < chars.length() {
        i = i + 1 // skip closing quote
      }
    } else {
      let val_start = i
      while i < chars.length() && chars[i].to_int() != 32 {
        i = i + 1
      }
      value = chars[val_start:i].to_owned()
    }
    match key {
      "type" => track_type = value
      "name" => name = value
      "description" => description = value
      "visibility" => visibility = value
      "color" => color = value
      "url" => url = value
      _ => params[key] = value
    }
  }
  Some({
    track_type,
    name,
    description,
    visibility,
    color,
    url,
    extra_params: params,
  })
}
