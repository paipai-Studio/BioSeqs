/// pysam_bam_supplementary.mbt — Supplementary alignment handling
/// Inspired by pysam split-read / supplementary alignment support

///|
/// Supplementary alignment info.
pub struct PysamSupplementaryAlign {
  qname : String
  rname : String
  pos : Int32
  cigar : String
  mapq : Int32
  is_reverse : Bool
}

///|
/// Extract supplementary alignments from SAM.
pub fn pysam_extract_supplementary(
  sam_text : String,
) -> Array[PysamSupplementaryAlign] {
  let lines = sam_text.split("\n").to_array()
  let result : Array[PysamSupplementaryAlign] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 6 {
      continue
    }
    let flag = parse_int(fields[1].to_owned())
    if (flag & 2048) == 0 { // not supplementary
      continue
    }
    let qname = fields[0].to_owned()
    let rname = fields[2].to_owned()
    let pos = parse_int(fields[3].to_owned()) - 1
    let cigar = fields[5].to_owned()
    let mapq = parse_int(fields[4].to_owned())
    let is_reverse = (flag & 16) != 0
    result.push({ qname, rname, pos, cigar, mapq, is_reverse, })
  }
  result
}

///|
/// Find primary alignment for a supplementary.
pub fn pysam_find_primary(
  sam_text : String,
  qname : String,
) -> PysamSupplementaryAlign? {
  let lines = sam_text.split("\n").to_array()
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 6 {
      continue
    }
    let curr_qname = fields[0].to_owned()
    if curr_qname != qname {
      continue
    }
    let flag = parse_int(fields[1].to_owned())
    if (flag & 2048) != 0 || (flag & 256) != 0 { // skip supplementary/secondary
      continue
    }
    let rname = fields[2].to_owned()
    let pos = parse_int(fields[3].to_owned()) - 1
    let cigar = fields[5].to_owned()
    let mapq = parse_int(fields[4].to_owned())
    let is_reverse = (flag & 16) != 0
    return Some({ qname, rname, pos, cigar, mapq, is_reverse, })
  }
  None
}

///|
/// Reconstruct full read from primary + supplementary alignments.
pub fn pysam_reconstruct_read(
  primary : PysamSupplementaryAlign,
  supplementary : Array[PysamSupplementaryAlign],
) -> Array[(String, Int32, Int32)] {
  // Returns (rname, start, end) segments
  let segments : Array[(String, Int32, Int32)] = Array(capacity=0)
  // Primary segment
  let primary_len = pysam_sa_ref_len(primary.cigar)
  segments.push((primary.rname, primary.pos, primary.pos + primary_len))
  // Supplementary segments
  for i = 0; i < supplementary.length(); i = i + 1 {
    let sa = supplementary[i]
    let sa_len = pysam_sa_ref_len(sa.cigar)
    segments.push((sa.rname, sa.pos, sa.pos + sa_len))
  }
  // Sort by position
  segments.sort_by(fn(a, b) {
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
  segments
}

///|
/// Reference consuming length from CIGAR.
fn pysam_sa_ref_len(cigar : String) -> Int32 {
  let mut total = 0
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      num = num * 10 + (c - 48)
    } else {
      if c == 77 || c == 68 || c == 78 || c == 61 || c == 88 {
        total = total + num
      }
      num = 0
    }
  }
  total
}

///|
/// Count supplementary alignments in SAM.
pub fn pysam_count_supplementary(sam_text : String) -> Int32 {
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
    if (flag & 2048) != 0 {
      count = count + 1
    }
  }
  count
}

///|
/// Remove supplementary alignments from SAM.
pub fn pysam_remove_supplementary(sam_text : String) -> String {
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
    if fields.length() < 2 {
      let _ = sb.write(line + "\n")
      continue
    }
    let flag = parse_int(fields[1].to_owned())
    if (flag & 2048) == 0 {
      let _ = sb.write(line + "\n")
    }
  }
  sb.to_string()
}

///|
/// Parse SA tag (supplementary alignments).
pub fn pysam_parse_sa_tag(sa_tag : String) -> Array[PysamSupplementaryAlign] {
  // SA tag format: rname,pos,strand,CIGAR,mapQ,NM;...
  let result : Array[PysamSupplementaryAlign] = Array(capacity=0)
  let entries = sa_tag.split(";").to_array()
  for i = 0; i < entries.length(); i = i + 1 {
    let entry = entries[i].to_owned()
    if entry.length() == 0 {
      continue
    }
    let parts = entry.split(",").to_array()
    if parts.length() < 6 {
      continue
    }
    let rname = parts[0].to_owned()
    let pos = parse_int(parts[1].to_owned()) - 1
    let strand = parts[2].to_owned()
    let cigar = parts[3].to_owned()
    let mapq = parse_int(parts[4].to_owned())
    let is_reverse = strand == "-"
    result.push({ qname: "", rname, pos, cigar, mapq, is_reverse, })
  }
  result
}
