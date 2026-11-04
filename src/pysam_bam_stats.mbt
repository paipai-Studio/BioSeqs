/// pysam_bam_stats.mbt — Comprehensive BAM statistics
/// Inspired by pysam/samtools stats: SN, RN, FF sections

///|
/// BAM statistics summary.
pub struct PysamBamStats {
  n_reads : Int32
  n_mapped : Int32
  n_unmapped : Int32
  n_paired : Int32
  n_proper_pair : Int32
  n_duplicate : Int32
  n_secondary : Int32
  n_supplementary : Int32
  n_forward : Int32
  n_reverse : Int32
  n_read1 : Int32
  n_read2 : Int32
  n_qcfail : Int32
  avg_mapq : Double
  avg_seq_len : Double
  bases_mapped : Int32
  error_rate : Double
}

///|
/// Calculate comprehensive BAM statistics from SAM text.
pub fn pysam_bam_stats(sam_text : String) -> PysamBamStats {
  let lines = sam_text.split("\n").to_array()
  let mut n_reads = 0
  let mut n_mapped = 0
  let mut n_unmapped = 0
  let mut n_paired = 0
  let mut n_proper_pair = 0
  let mut n_duplicate = 0
  let mut n_secondary = 0
  let mut n_supplementary = 0
  let mut n_forward = 0
  let mut n_reverse = 0
  let mut n_read1 = 0
  let mut n_read2 = 0
  let mut n_qcfail = 0
  let mut total_mapq = 0
  let mut total_seq_len = 0
  let mut bases_mapped = 0
  let mut total_nm = 0
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 11 {
      continue
    }
    n_reads = n_reads + 1
    let flag = parse_int(fields[1].to_owned())
    let rname = fields[2].to_owned()
    let mapq = parse_int(fields[4].to_owned())
    let cigar = fields[5].to_owned()
    let seq = fields[9].to_owned()
    // Flag checks
    let is_paired = (flag & 1) != 0
    let is_proper = (flag & 2) != 0
    let is_unmapped = (flag & 4) != 0
    let is_reverse = (flag & 16) != 0
    let is_read1 = (flag & 64) != 0
    let is_read2 = (flag & 128) != 0
    let is_secondary = (flag & 256) != 0
    let is_qcfail = (flag & 512) != 0
    let is_duplicate = (flag & 1024) != 0
    let is_supplementary = (flag & 2048) != 0
    if is_paired {
      n_paired = n_paired + 1
    }
    if is_proper {
      n_proper_pair = n_proper_pair + 1
    }
    if is_unmapped || rname == "*" {
      n_unmapped = n_unmapped + 1
    } else {
      n_mapped = n_mapped + 1
    }
    if is_reverse {
      n_reverse = n_reverse + 1
    } else {
      n_forward = n_forward + 1
    }
    if is_read1 {
      n_read1 = n_read1 + 1
    }
    if is_read2 {
      n_read2 = n_read2 + 1
    }
    if is_secondary {
      n_secondary = n_secondary + 1
    }
    if is_qcfail {
      n_qcfail = n_qcfail + 1
    }
    if is_duplicate {
      n_duplicate = n_duplicate + 1
    }
    if is_supplementary {
      n_supplementary = n_supplementary + 1
    }
    total_mapq = total_mapq + mapq
    total_seq_len = total_seq_len + seq.length()
    // Bases mapped from CIGAR (M, =, X)
    bases_mapped = bases_mapped + pysam_stats_cigar_match_len(cigar)
    // NM tag
    for f = 11; f < fields.length(); f = f + 1 {
      let tag = fields[f].to_owned()
      if tag.has_prefix("NM:i:") {
        total_nm = total_nm + parse_int(tag[5:].to_owned())
      }
    }
  }
  let avg_mapq = if n_reads > 0 {
    total_mapq.to_double() / n_reads.to_double()
  } else {
    0.0
  }
  let avg_seq_len = if n_reads > 0 {
    total_seq_len.to_double() / n_reads.to_double()
  } else {
    0.0
  }
  let error_rate = if bases_mapped > 0 {
    total_nm.to_double() / bases_mapped.to_double()
  } else {
    0.0
  }
  {
    n_reads,
    n_mapped,
    n_unmapped,
    n_paired,
    n_proper_pair,
    n_duplicate,
    n_secondary,
    n_supplementary,
    n_forward,
    n_reverse,
    n_read1,
    n_read2,
    n_qcfail,
    avg_mapq,
    avg_seq_len,
    bases_mapped,
    error_rate,
  }
}

///|
/// Count M/=/X operations in CIGAR.
fn pysam_stats_cigar_match_len(cigar : String) -> Int32 {
  let mut total = 0
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      num = num * 10 + (c - 48)
    } else {
      // M, =, X consume query
      if c == 77 || c == 61 || c == 88 {
        total = total + num
      }
      num = 0
    }
  }
  total
}

///|
/// Format stats as text (like samtools stats output).
pub fn pysam_bam_stats_text(stats : PysamBamStats) -> String {
  let sb = StringBuilder()
  let _ = sb.write("SN\tNumber of reads:\t" + stats.n_reads.to_string() + "\n")
  let _ = sb.write(
    "SN\tNumber of mapped reads:\t" + stats.n_mapped.to_string() + "\n",
  )
  let _ = sb.write(
    "SN\tNumber of unmapped reads:\t" + stats.n_unmapped.to_string() + "\n",
  )
  let _ = sb.write(
    "SN\tNumber of paired reads:\t" + stats.n_paired.to_string() + "\n",
  )
  let _ = sb.write(
    "SN\tNumber of properly paired:\t" + stats.n_proper_pair.to_string() + "\n",
  )
  let _ = sb.write(
    "SN\tNumber of duplicates:\t" + stats.n_duplicate.to_string() + "\n",
  )
  let _ = sb.write(
    "SN\tNumber of secondary:\t" + stats.n_secondary.to_string() + "\n",
  )
  let _ = sb.write(
    "SN\tNumber of supplementary:\t" + stats.n_supplementary.to_string() + "\n",
  )
  let _ = sb.write("SN\tAverage MAPQ:\t" + stats.avg_mapq.to_string() + "\n")
  let _ = sb.write(
    "SN\tAverage sequence length:\t" + stats.avg_seq_len.to_string() + "\n",
  )
  let _ = sb.write(
    "SN\tBases mapped:\t" + stats.bases_mapped.to_string() + "\n",
  )
  let _ = sb.write("SN\tError rate:\t" + stats.error_rate.to_string() + "\n")
  sb.to_string()
}
