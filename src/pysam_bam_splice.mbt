/// pysam_bam_splice.mbt — Spliced alignment (CIGAR N operations) support
/// Inspired by pysam RNA-seq spliced alignment handling

///|
/// A splice junction.
pub struct PysamSpliceJunction {
  contig : String
  donor : Int32 // end of exon 1 (0-based, exclusive)
  acceptor : Int32 // start of exon 2 (0-based)
  strand : Bool // true = forward
  n_reads : Int32 // number of reads supporting this junction
}

///|
/// Extract splice junctions from SAM.
pub fn pysam_extract_junctions(sam_text : String) -> Array[PysamSpliceJunction] {
  let lines = sam_text.split("\n").to_array()
  let junctions : Array[PysamSpliceJunction] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 6 {
      continue
    }
    let rname = fields[2].to_owned()
    let flag = parse_int(fields[1].to_owned())
    let pos = parse_int(fields[3].to_owned()) - 1
    let cigar = fields[5].to_owned()
    let is_reverse = (flag & 16) != 0
    // Find N operations in CIGAR
    let mut ref_pos = pos
    let mut num = 0
    for c = 0; c < cigar.length(); c = c + 1 {
      let ch = cigar[c].to_int()
      if ch >= 48 && ch <= 57 {
        num = num * 10 + (ch - 48)
      } else {
        if ch == 78 { // N: splice
          let donor = ref_pos
          let acceptor = ref_pos + num
          // Check if junction already exists
          let mut found = false
          for j = 0; j < junctions.length(); j = j + 1 {
            if junctions[j].contig == rname &&
              junctions[j].donor == donor &&
              junctions[j].acceptor == acceptor {
              junctions[j] = {
                contig: junctions[j].contig,
                donor: junctions[j].donor,
                acceptor: junctions[j].acceptor,
                strand: junctions[j].strand,
                n_reads: junctions[j].n_reads + 1,
              }
              found = true
              break
            }
          }
          if !found {
            junctions.push({
              contig: rname,
              donor,
              acceptor,
              strand: !is_reverse,
              n_reads: 1,
            })
          }
          ref_pos = ref_pos + num
        }
        if ch == 77 || ch == 68 || ch == 78 || ch == 61 || ch == 88 {
          ref_pos = ref_pos + num
        }
        num = 0
      }
    }
  }
  junctions
}

///|
/// Get exons from a spliced CIGAR.
pub fn pysam_get_exons(pos : Int32, cigar : String) -> Array[(Int32, Int32)] {
  let exons : Array[(Int32, Int32)] = Array(capacity=0)
  let mut ref_pos = pos
  let mut exon_start = pos
  let mut in_exon = false
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      num = num * 10 + (c - 48)
    } else {
      match c {
        77 | 61 | 88 => {
          if !in_exon {
            exon_start = ref_pos
            in_exon = true
          }
          ref_pos = ref_pos + num
        }
        78 => {
          if in_exon {
            exons.push((exon_start, ref_pos))
            in_exon = false
          }
          ref_pos = ref_pos + num
        }
        68 => ref_pos = ref_pos + num
        _ => ()
      }
      num = 0
    }
  }
  if in_exon {
    exons.push((exon_start, ref_pos))
  }
  exons
}

///|
/// Count spliced reads in SAM.
pub fn pysam_count_spliced(sam_text : String) -> Int32 {
  let lines = sam_text.split("\n").to_array()
  let mut count = 0
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 6 {
      continue
    }
    let cigar = fields[5].to_owned()
    // Check for N in CIGAR
    for c = 0; c < cigar.length(); c = c + 1 {
      if cigar[c].to_int() == 78 { // 'N'
        count = count + 1
        break
      }
    }
  }
  count
}

///|
/// Filter junctions by minimum read support.
pub fn pysam_filter_junctions(
  junctions : Array[PysamSpliceJunction],
  min_reads~ : Int32 = 2,
) -> Array[PysamSpliceJunction] {
  let result : Array[PysamSpliceJunction] = Array(capacity=0)
  for i = 0; i < junctions.length(); i = i + 1 {
    if junctions[i].n_reads >= min_reads {
      result.push(junctions[i])
    }
  }
  result
}
