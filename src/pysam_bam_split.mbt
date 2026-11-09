/// pysam_bam_split.mbt — Split SAM/BAM by chromosome or tag
/// Inspired by pysam/samtools split

///|
/// Split SAM by chromosome into separate SAM texts.
/// Returns (contig_name, sam_text) pairs.
pub fn pysam_split_by_chrom(sam_text : String) -> Array[(String, String)] {
  let lines = sam_text.split("\n").to_array()
  // Collect headers
  let headers : Array[String] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 {
      continue
    }
    if line[0].to_int() == 64 {
      headers.push(line)
    } else {
      break
    }
  }
  // Group records by chromosome
  let chroms : Array[String] = Array(capacity=0)
  let records : Array[Array[String]] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 3 {
      continue
    }
    let rname = fields[2].to_owned()
    let mut found = false
    for j = 0; j < chroms.length(); j = j + 1 {
      if chroms[j] == rname {
        records[j].push(line)
        found = true
        break
      }
    }
    if !found {
      chroms.push(rname)
      records.push([line])
    }
  }
  // Build output
  let result : Array[(String, String)] = Array(capacity=0)
  for j = 0; j < chroms.length(); j = j + 1 {
    let sb = StringBuilder()
    for h = 0; h < headers.length(); h = h + 1 {
      let _ = sb.write(headers[h] + "\n")
    }
    for r = 0; r < records[j].length(); r = r + 1 {
      let _ = sb.write(records[j][r] + "\n")
    }
    result.push((chroms[j], sb.to_string()))
  }
  result
}

///|
/// Split SAM by read group.
pub fn pysam_split_by_readgroup(sam_text : String) -> Array[(String, String)] {
  let lines = sam_text.split("\n").to_array()
  let headers : Array[String] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 {
      continue
    }
    if line[0].to_int() == 64 {
      headers.push(line)
    } else {
      break
    }
  }
  let groups : Array[String] = Array(capacity=0)
  let records : Array[Array[String]] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    // Find RG tag
    let fields = line.split("\t").to_array()
    let mut rg = "unknown"
    for f = 0; f < fields.length(); f = f + 1 {
      let tag = fields[f].to_owned()
      if tag.has_prefix("RG:Z:") {
        rg = tag[5:].to_owned()
        break
      }
    }
    let mut found = false
    for j = 0; j < groups.length(); j = j + 1 {
      if groups[j] == rg {
        records[j].push(line)
        found = true
        break
      }
    }
    if !found {
      groups.push(rg)
      records.push([line])
    }
  }
  let result : Array[(String, String)] = Array(capacity=0)
  for j = 0; j < groups.length(); j = j + 1 {
    let sb = StringBuilder()
    for h = 0; h < headers.length(); h = h + 1 {
      let _ = sb.write(headers[h] + "\n")
    }
    for r = 0; r < records[j].length(); r = r + 1 {
      let _ = sb.write(records[j][r] + "\n")
    }
    result.push((groups[j], sb.to_string()))
  }
  result
}

///|
/// Split SAM by mapped/unmapped.
pub fn pysam_split_mapped(sam_text : String) -> (String, String) {
  let lines = sam_text.split("\n").to_array()
  let headers_sb = StringBuilder()
  let mapped_sb = StringBuilder()
  let unmapped_sb = StringBuilder()
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 {
      continue
    }
    if line[0].to_int() == 64 {
      let _ = headers_sb.write(line + "\n")
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 2 {
      continue
    }
    let flag = parse_int(fields[1].to_owned())
    if (flag & 4) != 0 {
      let _ = unmapped_sb.write(line + "\n")
    } else {
      let _ = mapped_sb.write(line + "\n")
    }
  }
  let header = headers_sb.to_string()
  (header + mapped_sb.to_string(), header + unmapped_sb.to_string())
}
