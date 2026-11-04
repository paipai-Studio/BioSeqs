/// pysam_bam_read_pair.mbt — Read pair validation and statistics
/// Inspired by pysam paired-end read handling

///|
/// Read pair statistics.
pub struct PysamReadPairStats {
  n_total : Int32
  n_paired : Int32
  n_proper_pair : Int32
  n_both_mapped : Int32
  n_only_first_mapped : Int32
  n_only_second_mapped : Int32
  n_both_unmapped : Int32
  n_same_strand : Int32
  n_opposite_strand : Int32
  n_discordant : Int32 // not proper pair but both mapped
}

///|
/// Analyze read pairs in SAM.
pub fn pysam_read_pair_stats(sam_text : String) -> PysamReadPairStats {
  let lines = sam_text.split("\n").to_array()
  let mut n_total = 0
  let mut n_paired = 0
  let mut n_proper_pair = 0
  let mut n_both_mapped = 0
  let mut n_only_first_mapped = 0
  let mut n_only_second_mapped = 0
  let mut n_both_unmapped = 0
  let mut n_same_strand = 0
  let mut n_opposite_strand = 0
  let mut n_discordant = 0
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 7 {
      continue
    }
    n_total = n_total + 1
    let flag = parse_int(fields[1].to_owned())
    let rnext = fields[6].to_owned()
    let is_paired = (flag & 1) != 0
    let is_proper = (flag & 2) != 0
    let is_unmapped = (flag & 4) != 0
    let is_mate_unmapped = (flag & 8) != 0
    let is_reverse = (flag & 16) != 0
    let is_mate_reverse = (flag & 32) != 0
    let is_read1 = (flag & 64) != 0
    let is_read2 = (flag & 128) != 0
    if !is_paired {
      continue
    }
    n_paired = n_paired + 1
    if is_proper {
      n_proper_pair = n_proper_pair + 1
    }
    if !is_unmapped && !is_mate_unmapped {
      n_both_mapped = n_both_mapped + 1
      if !is_proper {
        n_discordant = n_discordant + 1
      }
    } else if !is_unmapped && is_mate_unmapped {
      if is_read1 {
        n_only_first_mapped = n_only_first_mapped + 1
      } else if is_read2 {
        n_only_second_mapped = n_only_second_mapped + 1
      }
    } else if is_unmapped && !is_mate_unmapped {
      if is_read1 {
        n_only_first_mapped = n_only_first_mapped + 1
      } else if is_read2 {
        n_only_second_mapped = n_only_second_mapped + 1
      }
    } else {
      n_both_unmapped = n_both_unmapped + 1
    }
    if !is_unmapped && !is_mate_unmapped {
      if is_reverse == is_mate_reverse {
        n_same_strand = n_same_strand + 1
      } else {
        n_opposite_strand = n_opposite_strand + 1
      }
    }
    ignore(rnext)
  }
  {
    n_total,
    n_paired,
    n_proper_pair,
    n_both_mapped,
    n_only_first_mapped,
    n_only_second_mapped,
    n_both_unmapped,
    n_same_strand,
    n_opposite_strand,
    n_discordant,
  }
}

///|
/// Check if a read is properly paired.
pub fn pysam_is_proper_pair(flag : Int32) -> Bool {
  (flag & 1) != 0 && (flag & 2) != 0
}

///|
/// Check if a read is a singleton (one mate mapped, other unmapped).
pub fn pysam_is_singleton(flag : Int32) -> Bool {
  (flag & 1) != 0 &&
  (((flag & 4) != 0 && (flag & 8) == 0) || ((flag & 4) == 0 && (flag & 8) != 0))
}

///|
/// Extract discordant reads (both mapped but not proper pair).
pub fn pysam_extract_discordant(sam_text : String) -> String {
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
      continue
    }
    let flag = parse_int(fields[1].to_owned())
    let is_paired = (flag & 1) != 0
    let is_proper = (flag & 2) != 0
    let is_unmapped = (flag & 4) != 0
    let is_mate_unmapped = (flag & 8) != 0
    if is_paired && !is_proper && !is_unmapped && !is_mate_unmapped {
      let _ = sb.write(line + "\n")
    }
  }
  sb.to_string()
}
