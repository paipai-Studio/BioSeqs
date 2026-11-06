// Pysam third-round: Flagstat and statistics (items 6-9)
// flagstat, full flagstat categories, phred utilities, mapq histogram

///|
pub struct Pysam3FlagStat {
  total : Int32
  qc_fail : Int32
  duplicate : Int32
  secondary : Int32
  supplementary : Int32
  mapped : Int32
  paired : Int32
  read1 : Int32
  read2 : Int32
  properly_paired : Int32
  mate_mapped : Int32
  singleton : Int32
  mate_diff_chr : Int32
  mate_diff_chr_mapq5 : Int32
}

///|
pub fn pysam3_flagstat_empty() -> Pysam3FlagStat {
  {
    total: 0,
    qc_fail: 0,
    duplicate: 0,
    secondary: 0,
    supplementary: 0,
    mapped: 0,
    paired: 0,
    read1: 0,
    read2: 0,
    properly_paired: 0,
    mate_mapped: 0,
    singleton: 0,
    mate_diff_chr: 0,
    mate_diff_chr_mapq5: 0,
  }
}

///|
fn pysam3_has_flag(flag : Int32, bit : Int32) -> Bool {
  flag / (1 << bit) % 2 == 1
}

///|
pub fn pysam3_flagstat_count(
  stats : Pysam3FlagStat,
  flag : Int32,
  mapq : Int32,
  mate_ref_id : Int32,
  ref_id : Int32,
) -> Pysam3FlagStat {
  let s = { ..stats, total: stats.total + 1, }
  let s = if pysam3_has_flag(flag, 9) {
    { ..s, qc_fail: s.qc_fail + 1, }
  } else {
    s
  }
  let s = if pysam3_has_flag(flag, 11) {
    { ..s, duplicate: s.duplicate + 1, }
  } else {
    s
  }
  let s = if pysam3_has_flag(flag, 8) {
    { ..s, secondary: s.secondary + 1, }
  } else {
    s
  }
  let s = if pysam3_has_flag(flag, 12) {
    { ..s, supplementary: s.supplementary + 1, }
  } else {
    s
  }
  let is_mapped = not(pysam3_has_flag(flag, 2))
  let s = if is_mapped { { ..s, mapped: s.mapped + 1, } } else { s }
  let is_paired = pysam3_has_flag(flag, 0)
  let s = if is_paired { { ..s, paired: s.paired + 1, } } else { s }
  let s = if pysam3_has_flag(flag, 6) {
    { ..s, read1: s.read1 + 1, }
  } else {
    s
  }
  let s = if pysam3_has_flag(flag, 7) {
    { ..s, read2: s.read2 + 1, }
  } else {
    s
  }
  let is_proper = pysam3_has_flag(flag, 1)
  let s = if is_proper {
    { ..s, properly_paired: s.properly_paired + 1, }
  } else {
    s
  }
  let mate_mapped_ok = is_paired && not(pysam3_has_flag(flag, 3))
  let s = if mate_mapped_ok {
    { ..s, mate_mapped: s.mate_mapped + 1, }
  } else {
    s
  }
  let is_singleton = is_paired && not(is_mapped) && mate_mapped_ok
  let s = if is_singleton { { ..s, singleton: s.singleton + 1, } } else { s }
  let mate_diff = is_mapped && mate_mapped_ok && mate_ref_id != ref_id
  let s = if mate_diff {
    { ..s, mate_diff_chr: s.mate_diff_chr + 1, }
  } else {
    s
  }
  let mate_diff_q5 = mate_diff && mapq >= 5
  let s = if mate_diff_q5 {
    { ..s, mate_diff_chr_mapq5: s.mate_diff_chr_mapq5 + 1, }
  } else {
    s
  }
  s
}

///|
pub fn pysam3_flagstat_to_text(stats : Pysam3FlagStat) -> String {
  let sb = StringBuilder()
  let lines = [
    (stats.total, "total"),
    (stats.secondary, "secondary"),
    (stats.supplementary, "supplementary"),
    (stats.duplicate, "duplicates"),
    (stats.mapped, "mapped"),
    (stats.paired, "paired in sequencing"),
    (stats.read1, "read1"),
    (stats.read2, "read2"),
    (stats.properly_paired, "properly paired"),
    (stats.mate_mapped, "with itself and mate mapped"),
    (stats.singleton, "singletons"),
    (stats.mate_diff_chr, "with mate mapped to a different chr"),
    (stats.mate_diff_chr_mapq5, "with mate mapped to a different chr (mapQ>=5)"),
  ]
  for i in 0..<lines.length() {
    let (count, label) = lines[i]
    let pct = if stats.total > 0 {
      count.to_double() / stats.total.to_double() * 100.0
    } else {
      0.0
    }
    let _ = sb.write(count.to_string())
    let _ = sb.write(" + ")
    let _ = sb.write("0 in total (")
    let _ = sb.write(pct.to_string())
    let _ = sb.write(") : ")
    let _ = sb.write(label)
    let _ = sb.write("\n")
  }
  sb.to_string()
}

// --- Phred score utilities ---

///|
pub fn pysam3_phred_to_prob(phred : Double) -> Double {
  @math.pow(10.0, -phred / 10.0)
}

///|
pub fn pysam3_prob_to_phred(prob : Double) -> Double {
  if prob <= 0.0 {
    60.0
  } else if prob >= 1.0 {
    0.0
  } else {
    -10.0 * @math.log10(prob)
  }
}

///|
pub fn pysam3_phred_to_error_rate(phred : Double) -> Double {
  @math.pow(10.0, -phred / 10.0)
}

///|
pub fn pysam3_error_rate_to_phred(error_rate : Double) -> Double {
  if error_rate <= 0.0 {
    60.0
  } else {
    -10.0 * @math.log10(error_rate)
  }
}

///|
pub fn pysam3_phred_to_sanger_char(phred : Int32) -> Int32 {
  phred + 33
}

///|
pub fn pysam3_sanger_char_to_phred(ch : Int32) -> Int32 {
  ch - 33
}

///|
pub fn pysam3_phred_to_solexa_char(phred : Int32) -> Int32 {
  phred + 64
}

///|
pub fn pysam3_solexa_char_to_phred(ch : Int32) -> Int32 {
  ch - 64
}

///|
pub fn pysam3_solexa_to_phred(solexa : Double) -> Double {
  10.0 * @math.log10(@math.pow(10.0, solexa / 10.0) + 1.0)
}

///|
pub fn pysam3_phred_to_solexa(phred : Double) -> Double {
  10.0 * @math.log10(@math.pow(10.0, phred / 10.0) - 1.0)
}

///|
pub fn pysam3_phred33_to_quality(qual_str : String) -> Array[Int32] {
  let quals : Array[Int32] = []
  for i in 0..<qual_str.length() {
    quals.push(qual_str[i].to_int() - 33)
  }
  quals
}

///|
pub fn pysam3_quality_to_phred33(quals : Array[Int32]) -> String {
  let sb = StringBuilder()
  for i in 0..<quals.length() {
    let _ = sb.write_char(Int32::unsafe_to_char(quals[i] + 33))
  }
  sb.to_string()
}

///|
pub fn pysam3_mean_phred(qual_str : String) -> Double {
  if qual_str.length() == 0 {
    return 0.0
  }
  let mut sum = 0
  for i in 0..<qual_str.length() {
    sum = sum + qual_str[i].to_int() - 33
  }
  sum.to_double() / qual_str.length().to_double()
}

// --- Mapping quality histogram ---

///|
pub fn pysam3_mapq_histogram(
  mapqs : Array[Int32],
  bin_size : Int32,
) -> Array[(Int32, Int32, Int32)] {
  // returns (bin_start, bin_end, count)
  if mapqs.length() == 0 {
    return []
  }
  let mut max_val = 0
  for i in 0..<mapqs.length() {
    if mapqs[i] > max_val {
      max_val = mapqs[i]
    }
  }
  let num_bins = max_val / bin_size + 1
  let bins : Array[Int32] = []
  for _ in 0..<num_bins {
    bins.push(0)
  }
  for i in 0..<mapqs.length() {
    let bin_idx = mapqs[i] / bin_size
    bins[bin_idx] = bins[bin_idx] + 1
  }
  let result : Array[(Int32, Int32, Int32)] = []
  for i in 0..<num_bins {
    result.push((i * bin_size, (i + 1) * bin_size - 1, bins[i]))
  }
  result
}

///|
pub fn pysam3_base_qual_histogram(
  qual_strs : Array[String],
  bin_size : Int32,
) -> Array[(Int32, Int32, Int32)] {
  let mut max_qual = 0
  let all_quals : Array[Int32] = []
  for s in 0..<qual_strs.length() {
    let qs = qual_strs[s]
    for i in 0..<qs.length() {
      let q = qs[i].to_int() - 33
      all_quals.push(q)
      if q > max_qual {
        max_qual = q
      }
    }
  }
  pysam3_mapq_histogram(all_quals, bin_size)
}
