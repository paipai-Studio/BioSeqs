// Pysam third-round: Coordinate and region operations (items 21-24)
// BAM→BED, intron detection, reference fetch, position-specific depth

///|
pub struct Pysam3BedRecord {
  chrom : String
  start : Int32 // 0-based
  end_ : Int32 // 1-based exclusive
  name : String
  score : Int32
  strand : String // "+" or "-"
}

///|
pub fn pysam3_bam_to_bed(
  chrom : String,
  ref_start : Int32,
  cigar : Array[(Int32, Int32)],
  read_name : String,
  mapq : Int32,
  is_reverse : Bool,
) -> Pysam3BedRecord {
  let mut ref_end = ref_start
  for i in 0..<cigar.length() {
    let (op, length) = cigar[i]
    if op == 0 || op == 2 || op == 3 || op == 7 || op == 8 {
      ref_end = ref_end + length
    }
  }
  {
    chrom,
    start: ref_start,
    end_: ref_end,
    name: read_name,
    score: mapq,
    strand: if is_reverse {
      "-"
    } else {
      "+"
    },
  }
}

///|
pub fn pysam3_bed_to_string(rec : Pysam3BedRecord) -> String {
  let sb = StringBuilder()
  let _ = sb.write(rec.chrom)
  let _ = sb.write("\t")
  let _ = sb.write(rec.start.to_string())
  let _ = sb.write("\t")
  let _ = sb.write(rec.end_.to_string())
  let _ = sb.write("\t")
  let _ = sb.write(rec.name)
  let _ = sb.write("\t")
  let _ = sb.write(rec.score.to_string())
  let _ = sb.write("\t")
  let _ = sb.write(rec.strand)
  sb.to_string()
}

///|
pub fn pysam3_bam_to_bed_multi(
  alignments : Array[
    (String, Int32, Array[(Int32, Int32)], String, Int32, Bool),
  ],
) -> Array[Pysam3BedRecord] {
  // (chrom, ref_start, cigar, read_name, mapq, is_reverse)
  let result : Array[Pysam3BedRecord] = []
  for i in 0..<alignments.length() {
    let (chrom, ref_start, cigar, read_name, mapq, is_reverse) = alignments[i]
    result.push(
      pysam3_bam_to_bed(chrom, ref_start, cigar, read_name, mapq, is_reverse),
    )
  }
  result
}

// --- Intron detection from RNA-seq alignments ---

///|
pub struct Pysam3Intron {
  chrom : String
  start : Int32 // 0-based
  end_ : Int32 // 1-based exclusive
  strand : String
  support : Int32 // number of reads supporting this intron
}

///|
pub fn pysam3_find_introns(
  chrom : String,
  ref_start : Int32,
  cigar : Array[(Int32, Int32)],
  is_reverse : Bool,
) -> Array[Pysam3Intron] {
  let introns : Array[Pysam3Intron] = []
  let mut pos = ref_start
  for i in 0..<cigar.length() {
    let (op, length) = cigar[i]
    if op == 3 {
      // N: reference skip = intron
      introns.push({
        chrom,
        start: pos,
        end_: pos + length,
        strand: if is_reverse {
          "-"
        } else {
          "+"
        },
        support: 1,
      })
    }
    if op == 0 || op == 2 || op == 3 || op == 7 || op == 8 {
      pos = pos + length
    }
  }
  introns
}

///|
pub fn pysam3_find_introns_multi(
  alignments : Array[(String, Int32, Array[(Int32, Int32)], Bool)],
) -> Array[Pysam3Intron] {
  // Merge introns from multiple alignments, summing support
  let all_introns : Array[Pysam3Intron] = []
  for i in 0..<alignments.length() {
    let (chrom, ref_start, cigar, is_reverse) = alignments[i]
    let introns = pysam3_find_introns(chrom, ref_start, cigar, is_reverse)
    for j in 0..<introns.length() {
      all_introns.push(introns[j])
    }
  }
  // Merge duplicates
  let merged : Array[Pysam3Intron] = []
  for i in 0..<all_introns.length() {
    let intron = all_introns[i]
    let mut found = false
    for j in 0..<merged.length() {
      let m = merged[j]
      if m.chrom == intron.chrom &&
        m.start == intron.start &&
        m.end_ == intron.end_ {
        merged[j] = { ..m, support: m.support + 1, }
        found = true
        break
      }
    }
    if not(found) {
      merged.push(intron)
    }
  }
  merged
}

// --- Reference sequence fetching from indexed FASTA ---

///|
pub fn pysam3_fetch_reference(
  ref_seq : String,
  start : Int32,
  end_ : Int32,
) -> String {
  let n = ref_seq.length()
  let s = if start < 0 { 0 } else { start }
  let e = if end_ > n { n } else { end_ }
  if s >= e {
    return ""
  }
  ref_seq[s:e].to_owned()
}

///|
pub fn pysam3_fetch_reference_complement(
  ref_seq : String,
  start : Int32,
  end_ : Int32,
) -> String {
  let sub = pysam3_fetch_reference(ref_seq, start, end_)
  let sb = StringBuilder()
  for i in 0..<sub.length() {
    let c = sub[i].to_int()
    let comp = match c {
      65 | 97 => 84 // A→T
      84 | 116 => 65 // T→A
      67 | 99 => 71 // C→G
      71 | 103 => 67 // G→C
      _ => 78 // N→N
    }
    let _ = sb.write_char(Int32::unsafe_to_char(comp))
  }
  sb.to_string()
}

///|
pub fn pysam3_fetch_reference_uppercase(
  ref_seq : String,
  start : Int32,
  end_ : Int32,
) -> String {
  pysam3_fetch_reference(ref_seq, start, end_).to_upper()
}

// --- Position-specific depth ---

///|
pub struct Pysam3DepthQuery {
  // Simple structure for position-specific depth queries
  positions : Array[Int32] // sorted reference positions of read starts
  lengths : Array[Int32] // read lengths (reference consuming)
}

///|
pub fn pysam3_depth_at_position(
  read_starts : Array[Int32],
  read_lengths : Array[Int32],
  target_pos : Int32,
) -> Int32 {
  let mut depth = 0
  for i in 0..<read_starts.length() {
    if i < read_lengths.length() {
      let start = read_starts[i]
      let end_ = start + read_lengths[i]
      if target_pos >= start && target_pos < end_ {
        depth = depth + 1
      }
    }
  }
  depth
}

///|
pub fn pysam3_depth_range(
  read_starts : Array[Int32],
  read_lengths : Array[Int32],
  start : Int32,
  end_ : Int32,
) -> Array[(Int32, Int32)] {
  // (position, depth) for each position in [start, end)
  let result : Array[(Int32, Int32)] = []
  for pos in start..<end_ {
    result.push((pos, pysam3_depth_at_position(read_starts, read_lengths, pos)))
  }
  result
}

///|
pub fn pysam3_depth_summary(
  read_starts : Array[Int32],
  read_lengths : Array[Int32],
  start : Int32,
  end_ : Int32,
) -> (Int32, Int32, Double, Int32) {
  // (min_depth, max_depth, mean_depth, total_positions)
  let depths = pysam3_depth_range(read_starts, read_lengths, start, end_)
  if depths.length() == 0 {
    return (0, 0, 0.0, 0)
  }
  let mut min_d = depths[0].1
  let mut max_d = depths[0].1
  let mut sum = 0
  for i in 0..<depths.length() {
    let d = depths[i].1
    if d < min_d {
      min_d = d
    }
    if d > max_d {
      max_d = d
    }
    sum = sum + d
  }
  let mean = sum.to_double() / depths.length().to_double()
  (min_d, max_d, mean, depths.length())
}

///|
pub fn pysam3_depth_query_dummy() -> Pysam3DepthQuery {
  { positions: [], lengths: [], }
}
