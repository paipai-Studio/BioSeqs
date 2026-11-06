/// pysam_bam_markdup.mbt — Mark duplicate reads
/// Inspired by pysam/samtools markdup

///|
/// Mark duplicates in SAM records.
/// A duplicate is identified by same position, same CIGAR, same mate position.
/// Returns SAM text with duplicate flag set.
pub fn pysam_markdup(sam_text : String) -> String {
  let lines = sam_text.split("\n").to_array()
  let sb = StringBuilder()
  // Track seen (rname, pos, cigar, mate_pos) keys
  let keys : Array[String] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 {
      let _ = sb.write("\n")
      continue
    }
    if line[0].to_int() == 64 { // header
      let _ = sb.write(line + "\n")
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 11 {
      let _ = sb.write(line + "\n")
      continue
    }
    let flag = parse_int(fields[1].to_owned())
    let rname = fields[2].to_owned()
    let pos = fields[3].to_owned()
    let cigar = fields[5].to_owned()
    let pnext = fields[7].to_owned()
    // Skip unmapped or secondary/supplementary
    if (flag & 4) != 0 || (flag & 256) != 0 || (flag & 2048) != 0 {
      let _ = sb.write(line + "\n")
      continue
    }
    // Create key
    let key = rname + ":" + pos + ":" + cigar + ":" + pnext
    let mut is_dup = false
    for j = 0; j < keys.length(); j = j + 1 {
      if keys[j] == key {
        is_dup = true
        break
      }
    }
    if is_dup {
      // Set duplicate flag (1024)
      let new_flag = flag | 1024
      fields[1] = new_flag.to_string()
      let new_line = pysam_join_fields(fields)
      let _ = sb.write(new_line + "\n")
    } else {
      keys.push(key)
      let _ = sb.write(line + "\n")
    }
  }
  sb.to_string()
}

///|
/// Remove duplicate reads (filter out marked duplicates).
pub fn pysam_rmdup(sam_text : String) -> String {
  let lines = sam_text.split("\n").to_array()
  let sb = StringBuilder()
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 {
      continue
    }
    if line[0].to_int() == 64 { // header
      let _ = sb.write(line + "\n")
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 2 {
      let _ = sb.write(line + "\n")
      continue
    }
    let flag = parse_int(fields[1].to_owned())
    if (flag & 1024) != 0 { // duplicate
      continue // skip
    }
    let _ = sb.write(line + "\n")
  }
  sb.to_string()
}

///|
/// Count duplicates in SAM.
pub fn pysam_count_duplicates(sam_text : String) -> Int32 {
  let lines = sam_text.split("\n").to_array()
  let mut count = 0
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 2 {
      continue
    }
    let flag = parse_int(fields[1].to_owned())
    if (flag & 1024) != 0 {
      count = count + 1
    }
  }
  count
}

///|
/// Join tab-separated fields.
fn pysam_join_fields(fields : Array[StringView]) -> String {
  let sb = StringBuilder()
  for i = 0; i < fields.length(); i = i + 1 {
    if i > 0 {
      let _ = sb.write_char('\t')
    }
    let _ = sb.write(fields[i])
  }
  sb.to_string()
}
