/// pysam_bam_adjacent.mbt — Adjacent/overlapping read pair analysis
/// Inspired by pysam overlapping fragment detection

///|
/// Check if two reads overlap on the reference.
pub fn pysam_reads_overlap(
  pos1 : Int32,
  cigar1 : String,
  pos2 : Int32,
  cigar2 : String,
) -> Bool {
  let end1 = pos1 + pysam_adj_ref_len(cigar1)
  let end2 = pos2 + pysam_adj_ref_len(cigar2)
  pos1 < end2 && pos2 < end1
}

///|
/// Reference consuming length.
fn pysam_adj_ref_len(cigar : String) -> Int32 {
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
/// Calculate overlap length between two reads.
pub fn pysam_overlap_length(
  pos1 : Int32,
  cigar1 : String,
  pos2 : Int32,
  cigar2 : String,
) -> Int32 {
  let end1 = pos1 + pysam_adj_ref_len(cigar1)
  let end2 = pos2 + pysam_adj_ref_len(cigar2)
  let overlap_start = if pos1 > pos2 { pos1 } else { pos2 }
  let overlap_end = if end1 < end2 { end1 } else { end2 }
  if overlap_end > overlap_start {
    overlap_end - overlap_start
  } else {
    0
  }
}

///|
/// Find overlapping read pairs in SAM.
pub fn pysam_find_overlapping_pairs(
  sam_text : String,
) -> Array[(String, String)] {
  let lines = sam_text.split("\n").to_array()
  // Group by qname
  let qnames : Array[String] = Array(capacity=0)
  let positions : Array[Int32] = Array(capacity=0)
  let cigars : Array[String] = Array(capacity=0)
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
    if (flag & 1) == 0 || (flag & 4) != 0 {
      continue
    }
    let qname = fields[0].to_owned()
    let pos = parse_int(fields[3].to_owned()) - 1
    let cigar = fields[5].to_owned()
    qnames.push(qname)
    positions.push(pos)
    cigars.push(cigar)
  }
  // Find overlapping pairs
  let result : Array[(String, String)] = Array(capacity=0)
  for i = 0; i < qnames.length(); i = i + 1 {
    for j = i + 1; j < qnames.length(); j = j + 1 {
      if qnames[i] == qnames[j] {
        if pysam_reads_overlap(positions[i], cigars[i], positions[j], cigars[j]) {
          result.push((qnames[i], "overlap"))
        }
      }
    }
  }
  result
}

///|
/// Check if read pair is adjacent (no gap between mates).
pub fn pysam_is_adjacent(pos1 : Int32, cigar1 : String, pos2 : Int32) -> Bool {
  let end1 = pos1 + pysam_adj_ref_len(cigar1)
  end1 == pos2 || pos2 + 1 == end1
}

///|
/// Count adjacent read pairs.
pub fn pysam_count_adjacent_pairs(sam_text : String) -> Int32 {
  let lines = sam_text.split("\n").to_array()
  // Group by qname
  let qnames : Array[String] = Array(capacity=0)
  let positions : Array[Int32] = Array(capacity=0)
  let cigars : Array[String] = Array(capacity=0)
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
    if (flag & 1) == 0 || (flag & 4) != 0 {
      continue
    }
    let qname = fields[0].to_owned()
    let pos = parse_int(fields[3].to_owned()) - 1
    let cigar = fields[5].to_owned()
    qnames.push(qname)
    positions.push(pos)
    cigars.push(cigar)
  }
  let mut count = 0
  for i = 0; i < qnames.length(); i = i + 1 {
    for j = i + 1; j < qnames.length(); j = j + 1 {
      if qnames[i] == qnames[j] {
        if pysam_is_adjacent(positions[i], cigars[i], positions[j]) ||
          pysam_is_adjacent(positions[j], cigars[j], positions[i]) {
          count = count + 1
        }
      }
    }
  }
  count
}
