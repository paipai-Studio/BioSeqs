/// pysam_bam_fetch_region.mbt — Region-based BAM queries
/// Inspired by pysam AlignmentFile.fetch(contig, start, stop)

///|
/// A BAM record summary for region queries.
pub struct PysamBamRegionRecord {
  qname : String
  flag : Int32
  rname : String
  pos : Int32 // 0-based leftmost mapping position
  mapq : Int32
  cigar : String
  rnext : String
  pnext : Int32
  tlen : Int32
  seq : String
  qual : String
}

///|
/// Fetch records overlapping a genomic region.
/// Parses SAM text and returns records within [start, end) on contig.
pub fn pysam_fetch_region(
  sam_text : String,
  contig : String,
  start : Int32,
  end : Int32,
) -> Array[PysamBamRegionRecord] {
  let lines = sam_text.split("\n").to_array()
  let results : Array[PysamBamRegionRecord] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 { // '@'
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 11 {
      continue
    }
    let rname = fields[2].to_owned()
    if rname != contig && rname != "*" {
      continue
    }
    let pos = parse_int(fields[3].to_owned()) - 1 // SAM is 1-based
    if pos < 0 {
      continue
    }
    // Check overlap with [start, end)
    // Need to estimate read end from CIGAR
    let cigar = fields[5].to_owned()
    let read_len = pysam_cigar_ref_len(cigar)
    let read_end = pos + read_len
    // Overlap: pos < end && read_end > start
    if pos < end && read_end > start {
      let qname = fields[0].to_owned()
      let flag = parse_int(fields[1].to_owned())
      let mapq = parse_int(fields[4].to_owned())
      let rnext = fields[6].to_owned()
      let pnext = parse_int(fields[7].to_owned()) - 1
      let tlen = parse_int(fields[8].to_owned())
      let seq = fields[9].to_owned()
      let qual = fields[10].to_owned()
      results.push({
        qname,
        flag,
        rname,
        pos,
        mapq,
        cigar,
        rnext,
        pnext,
        tlen,
        seq,
        qual,
      })
    }
  }
  results
}

///|
/// Fetch records from multiple regions.
pub fn pysam_fetch_regions(
  sam_text : String,
  regions : Array[(String, Int32, Int32)],
) -> Array[PysamBamRegionRecord] {
  let all : Array[PysamBamRegionRecord] = Array(capacity=0)
  for i = 0; i < regions.length(); i = i + 1 {
    let (contig, start, end) = regions[i]
    let recs = pysam_fetch_region(sam_text, contig, start, end)
    for j = 0; j < recs.length(); j = j + 1 {
      all.push(recs[j])
    }
  }
  all
}

///|
/// Count records overlapping a region.
pub fn pysam_count_region(
  sam_text : String,
  contig : String,
  start : Int32,
  end : Int32,
) -> Int32 {
  pysam_fetch_region(sam_text, contig, start, end).length()
}

///|
/// Calculate reference-consuming length from CIGAR string.
fn pysam_cigar_ref_len(cigar : String) -> Int32 {
  let mut total = 0
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 { // digit
      num = num * 10 + (c - 48)
    } else {
      // M, D, N, =, X consume reference
      if c == 77 || c == 68 || c == 78 || c == 61 || c == 88 {
        total = total + num
      }
      num = 0
    }
  }
  total
}

///|
/// Fetch until N records found (head).
pub fn pysam_head(sam_text : String, n : Int32) -> Array[PysamBamRegionRecord] {
  let lines = sam_text.split("\n").to_array()
  let results : Array[PysamBamRegionRecord] = Array(capacity=0)
  let mut count = 0
  for i = 0; i < lines.length() && count < n; i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 11 {
      continue
    }
    let qname = fields[0].to_owned()
    let flag = parse_int(fields[1].to_owned())
    let rname = fields[2].to_owned()
    let pos = parse_int(fields[3].to_owned()) - 1
    let mapq = parse_int(fields[4].to_owned())
    let cigar = fields[5].to_owned()
    let rnext = fields[6].to_owned()
    let pnext = parse_int(fields[7].to_owned()) - 1
    let tlen = parse_int(fields[8].to_owned())
    let seq = fields[9].to_owned()
    let qual = fields[10].to_owned()
    results.push({
      qname,
      flag,
      rname,
      pos,
      mapq,
      cigar,
      rnext,
      pnext,
      tlen,
      seq,
      qual,
    })
    count = count + 1
  }
  results
}
