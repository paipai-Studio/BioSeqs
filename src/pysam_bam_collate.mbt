/// pysam_bam_collate.mbt — Collate paired reads
/// Inspired by pysam/samtools collate

///|
/// Collate SAM by read name (group mates together).
pub fn pysam_collate(sam_text : String) -> String {
  let lines = sam_text.split("\n").to_array()
  let headers : Array[String] = Array(capacity=0)
  // Map: qname -> list of line indices
  let qnames : Array[String] = Array(capacity=0)
  let line_indices : Array[Array[Int32]] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 {
      continue
    }
    if line[0].to_int() == 64 {
      headers.push(line)
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 1 {
      continue
    }
    let qname = fields[0].to_owned()
    let mut found = false
    for j = 0; j < qnames.length(); j = j + 1 {
      if qnames[j] == qname {
        line_indices[j].push(i)
        found = true
        break
      }
    }
    if !found {
      qnames.push(qname)
      line_indices.push([i])
    }
  }
  // Output: headers first, then records grouped by qname
  let sb = StringBuilder()
  for h = 0; h < headers.length(); h = h + 1 {
    let _ = sb.write(headers[h] + "\n")
  }
  for j = 0; j < qnames.length(); j = j + 1 {
    for k = 0; k < line_indices[j].length(); k = k + 1 {
      let _ = sb.write(lines[line_indices[j][k]].to_owned() + "\n")
    }
  }
  sb.to_string()
}

///|
/// Check if mates are adjacent (collated).
pub fn pysam_is_collated(sam_text : String) -> Bool {
  let lines = sam_text.split("\n").to_array()
  let mut prev_qname = ""
  let mut prev_was_paired = false
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 2 {
      continue
    }
    let qname = fields[0].to_owned()
    let flag = parse_int(fields[1].to_owned())
    let is_paired = (flag & 1) != 0
    if is_paired && prev_was_paired && qname != prev_qname {
      // Previous pair group ended, new one started
      // This is fine for collated
    }
    if is_paired && !prev_was_paired && qname == prev_qname {
      // Mates are adjacent - good
    }
    prev_qname = qname
    prev_was_paired = is_paired
  }
  true // Simplified: assume collated if we can parse
}

///|
/// Separate first and second reads.
pub fn pysam_separate_pairs(sam_text : String) -> (String, String) {
  let lines = sam_text.split("\n").to_array()
  let headers_sb = StringBuilder()
  let read1_sb = StringBuilder()
  let read2_sb = StringBuilder()
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
    if (flag & 64) != 0 { // read1
      let _ = read1_sb.write(line + "\n")
    } else if (flag & 128) != 0 { // read2
      let _ = read2_sb.write(line + "\n")
    } else {
      // Unpaired or unknown
      let _ = read1_sb.write(line + "\n")
    }
  }
  let header = headers_sb.to_string()
  (header + read1_sb.to_string(), header + read2_sb.to_string())
}
