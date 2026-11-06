/// pysam_bam_template_length.mbt — Template length computation
/// Inspired by pysam AlignedSegment.template_length

///|
/// Compute template length (insert size) for a pair.
/// Uses outer coordinate distance.
pub fn pysam_template_length(
  pos1 : Int32,
  cigar1 : String,
  pos2 : Int32,
  cigar2 : String,
) -> Int32 {
  let end1 = pos1 + pysam_tl_ref_len(cigar1)
  let end2 = pos2 + pysam_tl_ref_len(cigar2)
  if pos1 <= pos2 {
    end2 - pos1
  } else {
    -(end1 - pos2)
  }
}

///|
/// Reference consuming length from CIGAR.
fn pysam_tl_ref_len(cigar : String) -> Int32 {
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
/// Infer template length from a single read's TLEN field.
pub fn pysam_infer_template_length(tlen : Int32) -> Int32 {
  if tlen < 0 {
    -tlen
  } else {
    tlen
  }
}

///|
/// Check if template length is consistent within a pair.
pub fn pysam_check_template_length(tlen1 : Int32, tlen2 : Int32) -> Bool {
  tlen1 == -tlen2
}

///|
/// Compute template length distribution statistics.
pub fn pysam_template_length_distribution(
  sam_text : String,
  bin_size~ : Int32 = 10,
) -> Array[(Int32, Int32)] {
  let lines = sam_text.split("\n").to_array()
  // Collect positive template lengths from first-in-pair
  let lengths : Array[Int32] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 9 {
      continue
    }
    let flag = parse_int(fields[1].to_owned())
    if (flag & 1) == 0 || (flag & 64) == 0 || (flag & 4) != 0 {
      continue
    }
    let tlen = parse_int(fields[8].to_owned())
    if tlen > 0 {
      lengths.push(tlen)
    }
  }
  if lengths.length() == 0 {
    return []
  }
  lengths.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  // Build histogram
  let hist : Array[(Int32, Int32)] = Array(capacity=0)
  let mut i = 0
  while i < lengths.length() {
    let bin = lengths[i] / bin_size * bin_size
    let mut count = 0
    while i < lengths.length() && lengths[i] / bin_size * bin_size == bin {
      count = count + 1
      i = i + 1
    }
    hist.push((bin, count))
  }
  hist
}

///|
/// Validate template lengths: check mate consistency.
pub fn pysam_validate_template_lengths(
  sam_text : String,
) -> Array[(String, Int32)] {
  let lines = sam_text.split("\n").to_array()
  let errors : Array[(String, Int32)] = Array(capacity=0)
  // Map qname -> tlen
  let qnames : Array[String] = Array(capacity=0)
  let tlens : Array[Int32] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 9 {
      continue
    }
    let qname = fields[0].to_owned()
    let tlen = parse_int(fields[8].to_owned())
    let mut found = false
    for j = 0; j < qnames.length(); j = j + 1 {
      if qnames[j] == qname {
        if tlens[j] != -tlen && tlen != 0 && tlens[j] != 0 {
          errors.push((qname, tlen))
        }
        found = true
        break
      }
    }
    if !found {
      qnames.push(qname)
      tlens.push(tlen)
    }
  }
  errors
}
