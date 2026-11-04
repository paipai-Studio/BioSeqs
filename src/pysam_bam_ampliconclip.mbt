/// pysam_bam_ampliconclip.mbt — Primer clipping from reads
/// Inspired by pysam/samtools ampliconclip

///|
/// Primer information.
pub struct PysamPrimer {
  name : String
  chrom : String
  start : Int32 // 0-based
  end : Int32 // exclusive
  strand : Bool // true = forward, false = reverse
}

///|
/// Clip primers from SAM reads.
/// Returns SAM text with primer sequences soft-clipped.
pub fn pysam_ampliconclip(
  sam_text : String,
  primers : Array[PysamPrimer],
) -> String {
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
    if fields.length() < 11 {
      let _ = sb.write(line + "\n")
      continue
    }
    let rname = fields[2].to_owned()
    let pos = parse_int(fields[3].to_owned()) - 1
    let flag = parse_int(fields[1].to_owned())
    let is_reverse = (flag & 16) != 0
    // Find matching primer
    let mut clip_len = 0
    for p = 0; p < primers.length(); p = p + 1 {
      let primer = primers[p]
      if primer.chrom != rname {
        continue
      }
      if !is_reverse && primer.strand {
        // Forward read, forward primer: clip from start
        if pos >= primer.start && pos < primer.end {
          clip_len = primer.end - pos
        }
      }
      if is_reverse && !primer.strand {
        // Reverse read, reverse primer: clip from end
        // Simplified: just mark clip length
        if pos >= primer.start && pos < primer.end {
          clip_len = primer.end - primer.start
        }
      }
    }
    if clip_len > 0 {
      // Soft-clip the primer region
      let cigar = fields[5].to_owned()
      let seq = fields[9].to_owned()
      let qual = fields[10].to_owned()
      let new_cigar = clip_len.to_string() + "S" + cigar
      fields[5] = new_cigar.to_string()
      // Adjust position
      fields[3] = (pos + clip_len + 1).to_string()
      // Truncate sequence
      if seq.length() > clip_len {
        fields[9] = seq[clip_len:].to_owned()
      }
      if qual.length() > clip_len && qual != "*" {
        fields[10] = qual[clip_len:].to_owned()
      }
      let new_line = pysam_ac_join(fields)
      let _ = sb.write(new_line + "\n")
    } else {
      let _ = sb.write(line + "\n")
    }
  }
  sb.to_string()
}

///|
/// Join tab-separated fields.
fn pysam_ac_join(fields : Array[StringView]) -> String {
  let sb = StringBuilder()
  for i = 0; i < fields.length(); i = i + 1 {
    if i > 0 {
      let _ = sb.write_char('\t')
    }
    let _ = sb.write(fields[i])
  }
  sb.to_string()
}

///|
/// Count reads with primers clipped.
pub fn pysam_count_clipped(sam_text : String) -> Int32 {
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
    // Check for soft-clip at start
    if cigar.length() > 0 && cigar[0].to_int() >= 48 && cigar[0].to_int() <= 57 {
      // Starts with a number, check if followed by S
      let mut j = 0
      while j < cigar.length() &&
            cigar[j].to_int() >= 48 &&
            cigar[j].to_int() <= 57 {
        j = j + 1
      }
      if j < cigar.length() && cigar[j].to_int() == 83 { // 'S'
        count = count + 1
      }
    }
  }
  count
}

///|
/// Suppress unused constructor warnings.
pub fn pysam_bam_unused_suppress_1() -> Unit {
  ignore(PysamPrimer::{ name: "", chrom: "", start: 0, end: 0, strand: true, })
}
