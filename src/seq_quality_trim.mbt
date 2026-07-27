// NGS quality trimming, adapter removal, and sequence quality control.
//
// Provides functions for preprocessing high-throughput sequencing reads,
// inspired by tools like Trimmomatic and Bioconductor's ShortRead package.
// Supports quality-based trimming, adapter removal, poly-A/poly-T trimming,
// length filtering, GC content filtering, and batch processing of FASTQ reads.

///|
/// Represents a single FASTQ read with its identifier, sequence, and quality string.
pub struct FastqRead {
  id : String
  sequence : String
  quality : String
}

///|
/// Construct a FastqRead with the given id, sequence, and quality string.
pub fn FastqRead::new(
  id : String,
  sequence : String,
  quality : String,
) -> FastqRead {
  FastqRead::{ id, sequence, quality }
}

///|
/// Represents the result of trimming a single read, including the original
/// and trimmed sequences/qualities, and a description of the trim type.
pub struct TrimResult {
  read_id : String
  original_seq : String
  trimmed_seq : String
  original_qual : String
  trimmed_qual : String
  trim_type : String
}

///|
/// Construct a TrimResult with the given fields.
pub fn TrimResult::new(
  read_id : String,
  original_seq : String,
  trimmed_seq : String,
  original_qual : String,
  trimmed_qual : String,
  trim_type : String,
) -> TrimResult {
  TrimResult::{
    read_id,
    original_seq,
    trimmed_seq,
    original_qual,
    trimmed_qual,
    trim_type,
  }
}

///|
/// Configuration for adapter trimming.
///
/// - adapters: Array of adapter sequences to search for
/// - min_overlap: Minimum overlap length for adapter detection
/// - allowed_mismatches: Number of allowed mismatches in adapter alignment
pub struct AdapterConfig {
  adapters : Array[String]
  min_overlap : Int
  allowed_mismatches : Int
}

///|
/// Construct an AdapterConfig with the given parameters.
pub fn AdapterConfig::new(
  adapters : Array[String],
  min_overlap : Int,
  allowed_mismatches : Int,
) -> AdapterConfig {
  AdapterConfig::{ adapters, min_overlap, allowed_mismatches }
}

///|
/// Configuration for quality trimming.
///
/// - min_quality: Minimum base quality score (Phred) to keep a base
/// - min_length: Minimum read length after trimming to keep
/// - window_size: Sliding window size for quality filtering
/// - avg_quality: Minimum average quality for the sliding window
pub struct QualityConfig {
  min_quality : Int
  min_length : Int
  window_size : Int
  avg_quality : Int
}

///|
/// Construct a QualityConfig with the given parameters.
pub fn QualityConfig::new(
  min_quality : Int,
  min_length : Int,
  window_size : Int,
  avg_quality : Int,
) -> QualityConfig {
  QualityConfig::{ min_quality, min_length, window_size, avg_quality }
}

///|
/// Statistics from a batch trimming operation.
///
/// Provides counts of total, kept, trimmed, and discarded reads,
/// along with total base pairs and GC content before and after trimming.
pub struct TrimStatistics {
  total_reads : Int
  kept_reads : Int
  trimmed_reads : Int
  discarded_reads : Int
  total_bp_before : Int
  total_bp_after : Int
  gc_content_before : Double
  gc_content_after : Double
}

///|
/// Construct a TrimStatistics with the given fields.
pub fn TrimStatistics::new(
  total_reads : Int,
  kept_reads : Int,
  trimmed_reads : Int,
  discarded_reads : Int,
  total_bp_before : Int,
  total_bp_after : Int,
  gc_content_before : Double,
  gc_content_after : Double,
) -> TrimStatistics {
  TrimStatistics::{
    total_reads,
    kept_reads,
    trimmed_reads,
    discarded_reads,
    total_bp_before,
    total_bp_after,
    gc_content_before,
    gc_content_after,
  }
}

///|
/// Convert a quality character (Phred+33) to its numeric quality score.
fn sqt_phred(c : UInt16) -> Int {
  c.to_int() - 33
}

///|
/// Check if a base character is G or C (guanosine or cytosine).
fn sqt_is_gc(base : UInt16) -> Bool {
  base == 71 || base == 67 || base == 103 || base == 99
}

///|
/// Compute the GC content of a sequence string as a Double between 0.0 and 1.0.
fn sqt_gc_content(seq : String) -> Double {
  let len = seq.length()
  if len == 0 {
    return 0.0
  }
  let mut gc = 0
  let mut i = 0
  while i < len {
    if sqt_is_gc(seq.unsafe_get(i)) {
      gc = gc + 1
    }
    i = i + 1
  }
  gc.to_double() / len.to_double()
}

///|
/// Quality-based trimming from both 5' and 3' ends.
///
/// Trims low-quality bases from both ends of the sequence. Bases with
/// quality scores below `min_qual` are removed. Returns a tuple of
/// (trimmed_sequence, trimmed_quality).
pub fn sqt_quality_trim(
  seq : String,
  qual : String,
  min_qual? : Int = 20,
) -> (String, String) {
  let len = seq.length()
  if len == 0 {
    return ("", "")
  }

  let mut start = 0
  while start < len {
    let q = sqt_phred(qual.unsafe_get(start))
    if q >= min_qual {
      break
    }
    start = start + 1
  }

  let mut end = len
  while end > start {
    let q = sqt_phred(qual.unsafe_get(end - 1))
    if q >= min_qual {
      break
    }
    end = end - 1
  }

  if start >= end {
    return ("", "")
  }

  (seq[start:end].to_owned(), qual[start:end].to_owned())
}

///|
/// Remove adapter sequence from a read.
///
/// Searches for the adapter in the sequence. When found at a position,
/// the sequence is truncated up to that position (adapter and everything
/// after is removed). Supports a configurable number of allowed mismatches.
/// Returns (trimmed_sequence, trimmed_quality).
pub fn sqt_trim_adapter(
  seq : String,
  qual : String,
  adapter : String,
  allowed_mismatches? : Int = 0,
) -> (String, String) {
  let seq_len = seq.length()
  let adapter_len = adapter.length()

  if adapter_len == 0 || seq_len == 0 {
    return (seq, qual)
  }

  if seq_len < adapter_len {
    return (seq, qual)
  }

  let mut best_pos = -1
  let mut best_matches = 0

  let min_overlap = adapter_len / 2
  let mut i = 0
  while i <= seq_len - min_overlap {
    let mut mismatches = 0
    let mut matches = 0
    let mut j = 0
    while j < adapter_len && i + j < seq_len {
      if seq.unsafe_get(i + j) == adapter.unsafe_get(j) {
        matches = matches + 1
      } else {
        mismatches = mismatches + 1
        if mismatches > allowed_mismatches {
          break
        }
      }
      j = j + 1
    }

    if matches >= min_overlap && mismatches <= allowed_mismatches {
      if matches > best_matches {
        best_matches = matches
        best_pos = i
      }
    }
    i = i + 1
  }

  if best_pos >= 0 && best_matches >= min_overlap {
    (seq[0:best_pos].to_owned(), qual[0:best_pos].to_owned())
  } else {
    (seq, qual)
  }
}

///|
/// Trim poly-A and poly-T tails from both ends of a sequence.
///
/// Removes consecutive A's and T's from the 5' and 3' ends of the
/// sequence. Returns (trimmed_sequence, trimmed_quality).
pub fn sqt_trim_poly_a(
  seq : String,
  qual : String,
  min_tail? : Int = 3,
) -> (String, String) {
  let len = seq.length()
  if len == 0 {
    return ("", "")
  }

  let mut start = 0
  let mut leading_at = 0
  while start < len {
    let base = seq.unsafe_get(start)
    if base == 65 || base == 84 {
      start = start + 1
      leading_at = leading_at + 1
    } else {
      break
    }
  }

  let mut end = len
  let mut trailing_at = 0
  while end > start {
    let base = seq.unsafe_get(end - 1)
    if base == 65 || base == 84 {
      end = end - 1
      trailing_at = trailing_at + 1
    } else {
      break
    }
  }

  if leading_at < min_tail {
    start = 0
  }
  if trailing_at < min_tail {
    end = len
  }

  if start >= end {
    return ("", "")
  }

  (seq[start:end].to_owned(), qual[start:end].to_owned())
}

///|
/// Discard reads shorter than `min_length`.
///
/// Returns `Some((sequence, quality))` if the read meets the minimum length,
/// or `None` if it should be discarded.
pub fn sqt_trim_by_length(
  seq : String,
  qual : String,
  min_length? : Int = 50,
) -> (String, String)? {
  if seq.length() >= min_length {
    Some((seq, qual))
  } else {
    None
  }
}

///|
/// Filter reads by GC content range.
///
/// Returns `Some((sequence, quality))` if the GC content of the read falls
/// within the specified range, or `None` if it should be discarded.
pub fn sqt_trim_by_gc(
  seq : String,
  qual : String,
  gc_min? : Double = 0.3,
  gc_max? : Double = 0.7,
) -> (String, String)? {
  let gc = sqt_gc_content(seq)
  if gc >= gc_min && gc <= gc_max {
    Some((seq, qual))
  } else {
    None
  }
}

///|
/// Perform batch trimming on an array of FASTQ reads.
///
/// Applies quality trimming, adapter removal (if adapter_config is provided),
/// and poly-A trimming to each read. The resulting reads are wrapped in
/// TrimResult structures describing what trimming was applied.
pub fn sqt_trim_reads(
  reads : Array[FastqRead],
  qual_config : QualityConfig,
  adapter_config : AdapterConfig?,
) -> Array[TrimResult] {
  let results : Array[TrimResult] = Array::new()

  for read in reads {
    let original_seq = read.sequence
    let original_qual = read.quality
    let mut trimmed_seq = original_seq
    let mut trimmed_qual = original_qual
    let trim_types : Array[String] = Array::new()

    let (q_seq, q_qual) = sqt_quality_trim(
      trimmed_seq,
      trimmed_qual,
      min_qual=qual_config.min_quality,
    )
    if q_seq.length() < trimmed_seq.length() {
      trim_types.push("quality")
    }
    trimmed_seq = q_seq
    trimmed_qual = q_qual

    if adapter_config is Some(_) {
      let ac = adapter_config.unwrap()
      for adapter in ac.adapters {
        let (a_seq, a_qual) = sqt_trim_adapter(
          trimmed_seq,
          trimmed_qual,
          adapter,
          allowed_mismatches=ac.allowed_mismatches,
        )
        if a_seq.length() < trimmed_seq.length() {
          trim_types.push("adapter:" + adapter)
        }
        trimmed_seq = a_seq
        trimmed_qual = a_qual
      }
    }

    let (p_seq, p_qual) = sqt_trim_poly_a(trimmed_seq, trimmed_qual)
    if p_seq.length() < trimmed_seq.length() {
      trim_types.push("poly_A")
    }
    trimmed_seq = p_seq
    trimmed_qual = p_qual

    let trim_type = if trim_types.length() > 0 {
      let mut result = ""
      let mut i = 0
      while i < trim_types.length() {
        if i > 0 {
          result = result + ","
        }
        result = result + trim_types[i]
        i = i + 1
      }
      result
    } else {
      "none"
    }

    results.push(TrimResult::{
      read_id: read.id,
      original_seq,
      trimmed_seq,
      original_qual,
      trimmed_qual,
      trim_type,
    })
  }

  results
}

///|
/// Compute trimming statistics from an array of trim results.
///
/// Compares the original reads with their trimmed counterparts to produce
/// aggregate statistics including read counts, base pair counts, and GC
/// content before and after trimming.
pub fn sqt_compute_stats(
  results : Array[TrimResult],
  original_reads : Array[FastqRead],
) -> TrimStatistics {
  let total_reads = original_reads.length()
  let mut kept_reads = 0
  let mut trimmed_reads = 0
  let mut discarded_reads = 0
  let mut total_bp_before = 0
  let mut total_bp_after = 0
  let mut gc_total_before = 0.0
  let mut gc_total_after = 0.0

  let mut i = 0
  while i < original_reads.length() {
    let read = original_reads[i]
    total_bp_before = total_bp_before + read.sequence.length()
    gc_total_before = gc_total_before + sqt_gc_content(read.sequence)
    i = i + 1
  }

  for result in results {
    total_bp_after = total_bp_after + result.trimmed_seq.length()
    gc_total_after = gc_total_after + sqt_gc_content(result.trimmed_seq)

    if result.trimmed_seq.length() == 0 {
      discarded_reads = discarded_reads + 1
    } else {
      kept_reads = kept_reads + 1
      if result.trim_type != "none" {
        trimmed_reads = trimmed_reads + 1
      }
    }
  }

  let gc_content_before = if total_reads > 0 {
    gc_total_before / total_reads.to_double()
  } else {
    0.0
  }

  let kept_count = kept_reads + trimmed_reads
  let gc_content_after = if kept_count > 0 {
    gc_total_after / kept_count.to_double()
  } else {
    0.0
  }

  TrimStatistics::{
    total_reads,
    kept_reads,
    trimmed_reads,
    discarded_reads,
    total_bp_before,
    total_bp_after,
    gc_content_before,
    gc_content_after,
  }
}

///|
/// Parse a FASTQ formatted string into an array of FastqRead objects.
///
/// Handles Sanger/Phred+33 encoded FASTQ format. Each record consists of
/// four lines: @header, sequence, +separator, and quality string.
pub fn sqt_fastq_parse(content : String) -> Array[FastqRead] {
  let reads : Array[FastqRead] = Array::new()
  let lines = content.split("\n").to_array()
  let n = lines.length()
  let mut i = 0

  while i + 3 < n {
    let line0 = lines[i]
    let line1 = lines[i + 1]
    let line2 = lines[i + 2]
    let line3 = lines[i + 3]

    let len0 = line0.length()
    let mut trimmed0 = len0
    if len0 > 0 && line0.unsafe_get(len0 - 1) == '\r'.to_int().to_uint16() {
      trimmed0 = len0 - 1
    }

    let len2 = line2.length()
    let mut trimmed2 = len2
    if len2 > 0 && line2.unsafe_get(len2 - 1) == '\r'.to_int().to_uint16() {
      trimmed2 = len2 - 1
    }

    if trimmed0 == 0 || trimmed2 == 0 {
      i = i + 4
      continue
    }

    if line0.unsafe_get(0) != '@'.to_int().to_uint16() {
      i = i + 4
      continue
    }

    if line2.unsafe_get(0) != '+'.to_int().to_uint16() {
      i = i + 4
      continue
    }

    let len1 = line1.length()
    let mut trimmed1 = len1
    if len1 > 0 && line1.unsafe_get(len1 - 1) == '\r'.to_int().to_uint16() {
      trimmed1 = len1 - 1
    }

    let len3 = line3.length()
    let mut trimmed3 = len3
    if len3 > 0 && line3.unsafe_get(len3 - 1) == '\r'.to_int().to_uint16() {
      trimmed3 = len3 - 1
    }

    if trimmed1 != trimmed3 {
      i = i + 4
      continue
    }

    let title = if trimmed0 > 1 { line0[1:trimmed0].to_owned() } else { "" }
    let (id, _) = match title.find(" ") {
      Some(pos) => (title[0:pos].to_owned(), title)
      None => (title, title)
    }

    let seq_str = if trimmed1 > 0 { line1[0:trimmed1].to_owned() } else { "" }
    let qual_str = if trimmed3 > 0 { line3[0:trimmed3].to_owned() } else { "" }

    reads.push(FastqRead::{ id, sequence: seq_str, quality: qual_str })
    i = i + 4
  }

  reads
}

///|
/// Serialize an array of FastqRead objects to a FASTQ formatted string.
///
/// Outputs Sanger/Phred+33 encoded FASTQ with one record per read.
pub fn sqt_fastq_serialize(reads : Array[FastqRead]) -> String {
  let mut result = ""
  for read in reads {
    result = result + "@" + read.id + "\n"
    result = result + read.sequence + "\n"
    result = result + "+\n"
    result = result + read.quality + "\n"
  }
  result
}

// ============================================================================
// Inline tests
// ============================================================================

///|
test "sqt_quality_trim_basic" {
  let seq = "ACGTACGT"
  let qual = "IIIIII##"
  let (t_seq, t_qual) = sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq.length(), 6)
  assert_eq(t_seq, "ACGTAC")
  assert_eq(t_qual.length(), 6)
}

///|
test "sqt_quality_trim_all_high_quality" {
  let seq = "ACGTACGT"
  let qual = "IIIIIIII"
  let (t_seq, _) = sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq, "ACGTACGT")
}

///|
test "sqt_quality_trim_all_low_quality" {
  let seq = "ACGTACGT"
  let qual = "########"
  let (t_seq, _) = sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq.length(), 0)
}

///|
test "sqt_quality_trim_5_prime_low" {
  let seq = "ACGTACGT"
  let qual = "##IIIIII"
  let (t_seq, _) = sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq.length(), 6)
  assert_eq(t_seq, "GTACGT")
}

///|
test "sqt_quality_trim_both_ends_low" {
  let seq = "ACGTACGTACGT"
  let qual = "##IIIIIIII##"
  let (t_seq, _) = sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq.length(), 8)
  assert_eq(t_seq, "GTACGTAC")
}

///|
test "sqt_quality_trim_empty" {
  let seq = ""
  let qual = ""
  let (t_seq, _) = sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq.length(), 0)
}

///|
test "sqt_quality_trim_high_threshold" {
  let seq = "ACGTACGT"
  let qual = "IIIIIIII"
  let (t_seq, _) = sqt_quality_trim(seq, qual, min_qual=41)
  assert_eq(t_seq.length(), 0)
}

///|
test "sqt_quality_trim_custom_threshold" {
  let seq = "ACGTACGT"
  let qual = "IIIIIIII"
  let (t_seq, _) = sqt_quality_trim(seq, qual, min_qual=10)
  assert_eq(t_seq, "ACGTACGT")
}

///|
test "sqt_trim_adapter_perfect_match" {
  let seq = "ACGTACGTAGATCGGAA"
  let qual = "IIIIIIIIIIIIIIII"
  let adapter = "AGATCGGAA"
  let (t_seq, _) = sqt_trim_adapter(
    seq, qual, adapter, allowed_mismatches=0,
  )
  assert_eq(t_seq, "ACGTACGT")
}

///|
test "sqt_trim_adapter_no_match" {
  let seq = "ACGTACGTACGTACGT"
  let qual = "IIIIIIIIIIIIIIII"
  let adapter = "GGGGGGGG"
  let (t_seq, _) = sqt_trim_adapter(
    seq, qual, adapter, allowed_mismatches=0,
  )
  assert_eq(t_seq, "ACGTACGTACGTACGT")
}

///|
test "sqt_trim_adapter_with_mismatches" {
  let seq = "ACGTACGTAGATCGGAA"
  let qual = "IIIIIIIIIIIIIIII"
  let adapter = "AGATCGGAA"
  let (t_seq, _) = sqt_trim_adapter(
    seq, qual, adapter, allowed_mismatches=2,
  )
  assert_eq(t_seq, "ACGTACGT")
}

///|
test "sqt_trim_adapter_empty_adapter" {
  let seq = "ACGTACGT"
  let qual = "IIIIIIII"
  let (t_seq, _) = sqt_trim_adapter(seq, qual, "", allowed_mismatches=0)
  assert_eq(t_seq, "ACGTACGT")
}

///|
test "sqt_trim_adapter_short_seq" {
  let seq = "AC"
  let qual = "II"
  let adapter = "AGATCGGAA"
  let (t_seq, _) = sqt_trim_adapter(
    seq, qual, adapter, allowed_mismatches=0,
  )
  assert_eq(t_seq, "AC")
}

///|
test "sqt_trim_adapter_at_beginning" {
  let seq = "AGATCGGAAACGTACGT"
  let qual = "IIIIIIIIIIIIIIII"
  let adapter = "AGATCGGAA"
  let (t_seq, _) = sqt_trim_adapter(
    seq, qual, adapter, allowed_mismatches=0,
  )
  assert_eq(t_seq.length(), 0)
}

///|
test "sqt_trim_poly_a_3_prime" {
  let seq = "ACGTACGGGAAAA"
  let qual = "IIIIIIIIIIIIII"
  let (t_seq, _) = sqt_trim_poly_a(seq, qual, min_tail=3)
  assert_eq(t_seq, "ACGTACGGG")
}

///|
test "sqt_trim_poly_t_3_prime" {
  let seq = "ACGTACGGGTTTT"
  let qual = "IIIIIIIIIIIIII"
  let (t_seq, _) = sqt_trim_poly_a(seq, qual, min_tail=3)
  assert_eq(t_seq, "ACGTACGGG")
}

///|
test "sqt_trim_poly_a_5_prime" {
  let seq = "AAAACGTACGT"
  let qual = "IIIIIIIIIII"
  let (t_seq, _) = sqt_trim_poly_a(seq, qual, min_tail=3)
  assert_eq(t_seq, "CGTACGT")
}

///|
test "sqt_trim_poly_a_no_tail" {
  let seq = "ACGTACGT"
  let qual = "IIIIIIII"
  let (t_seq, _) = sqt_trim_poly_a(seq, qual, min_tail=3)
  assert_eq(t_seq, "ACGTACGT")
}

///|
test "sqt_trim_poly_a_all_poly" {
  let seq = "AAAAAAAA"
  let qual = "IIIIIIII"
  let (t_seq, _) = sqt_trim_poly_a(seq, qual, min_tail=3)
  assert_eq(t_seq.length(), 0)
}

///|
test "sqt_trim_by_length_keep" {
  let seq = "ACGTACGTACGTACGT"
  let qual = "IIIIIIIIIIIIIIII"
  let result = sqt_trim_by_length(seq, qual, min_length=5)
  assert_true(result is Some)
}

///|
test "sqt_trim_by_length_discard" {
  let seq = "ACGT"
  let qual = "IIII"
  let result = sqt_trim_by_length(seq, qual, min_length=10)
  assert_true(result is None)
}

///|
test "sqt_trim_by_length_exact" {
  let seq = "ACGTACGT"
  let qual = "IIIIIIII"
  let result = sqt_trim_by_length(seq, qual, min_length=8)
  assert_true(result is Some)
}

///|
test "sqt_trim_by_gc_within_range" {
  let seq = "ACGTACGT"
  let qual = "IIIIIIII"
  let result = sqt_trim_by_gc(seq, qual, gc_min=0.3, gc_max=0.7)
  assert_true(result is Some)
}

///|
test "sqt_trim_by_gc_out_of_range" {
  let seq = "GGGGGGGG"
  let qual = "IIIIIIII"
  let result = sqt_trim_by_gc(seq, qual, gc_min=0.3, gc_max=0.7)
  assert_true(result is None)
}

///|
test "sqt_trim_by_gc_extreme_at" {
  let seq = "AAAAAAAA"
  let qual = "IIIIIIII"
  let result = sqt_trim_by_gc(seq, qual, gc_min=0.3, gc_max=0.7)
  assert_true(result is None)
}

///|
test "sqt_trim_by_gc_mid_range" {
  let seq = "ACGTACGTACGT"
  let qual = "IIIIIIIIIIII"
  let result = sqt_trim_by_gc(seq, qual, gc_min=0.4, gc_max=0.6)
  assert_true(result is Some)
}

///|
test "sqt_trim_reads_batch" {
  let reads = [
    FastqRead::{
      id: "read1",
      sequence: "ACGTACGT",
      quality: "IIIIIIII",
    },
    FastqRead::{
      id: "read2",
      sequence: "ACGTACGTAGATCGGAA",
      quality: "IIIIIIIIIIIIIIIII",
    },
  ]
  let qual_config = QualityConfig::{
    min_quality: 20,
    min_length: 10,
    window_size: 4,
    avg_quality: 20,
  }
  let adapter_config = AdapterConfig::{
    adapters: ["AGATCGGAA"],
    min_overlap: 4,
    allowed_mismatches: 0,
  }
  let results = sqt_trim_reads(reads, qual_config, Some(adapter_config))
  assert_eq(results.length(), 2)
  assert_eq(results[0].trimmed_seq, "ACGTACGT")
  assert_eq(results[1].trimmed_seq, "ACGTACGT")
}

///|
test "sqt_trim_reads_no_adapter" {
  let reads = [
    FastqRead::{
      id: "read1",
      sequence: "ACGTACGGGAAAA",
      quality: "IIIIIIIIIIIIII",
    },
  ]
  let qual_config = QualityConfig::{
    min_quality: 20,
    min_length: 4,
    window_size: 4,
    avg_quality: 20,
  }
  let results = sqt_trim_reads(reads, qual_config, None)
  assert_eq(results.length(), 1)
  assert_eq(results[0].trimmed_seq, "ACGTACGGG")
  assert_true(
    results[0].trim_type == "poly_A" ||
      results[0].trim_type == "poly_A,quality",
  )
}

///|
test "sqt_compute_stats_basic" {
  let reads = [
    FastqRead::{
      id: "r1",
      sequence: "ACGTACGT",
      quality: "IIIIIIII",
    },
    FastqRead::{
      id: "r2",
      sequence: "TTTTTTTT",
      quality: "IIIIIIII",
    },
  ]
  let results = [
    TrimResult::{
      read_id: "r1",
      original_seq: "ACGTACGT",
      trimmed_seq: "ACGTACGT",
      original_qual: "IIIIIIII",
      trimmed_qual: "IIIIIIII",
      trim_type: "none",
    },
    TrimResult::{
      read_id: "r2",
      original_seq: "TTTTTTTT",
      trimmed_seq: "",
      original_qual: "IIIIIIII",
      trimmed_qual: "",
      trim_type: "poly_A",
    },
  ]
  let stats = sqt_compute_stats(results, reads)
  assert_eq(stats.total_reads, 2)
  assert_eq(stats.kept_reads, 1)
  assert_eq(stats.discarded_reads, 1)
  assert_eq(stats.total_bp_before, 16)
  assert_eq(stats.total_bp_after, 8)
}

///|
test "sqt_compute_stats_empty" {
  let reads : Array[FastqRead] = []
  let results : Array[TrimResult] = []
  let stats = sqt_compute_stats(results, reads)
  assert_eq(stats.total_reads, 0)
  assert_eq(stats.kept_reads, 0)
  assert_eq(stats.total_bp_before, 0)
  assert_eq(stats.gc_content_before, 0.0)
}

///|
test "sqt_fastq_parse_basic" {
  let content = "@read1 desc\nACGTACGT\n+\nIIIIIIII\n"
  let reads = sqt_fastq_parse(content)
  assert_eq(reads.length(), 1)
  assert_eq(reads[0].id, "read1")
  assert_eq(reads[0].sequence, "ACGTACGT")
  assert_eq(reads[0].quality, "IIIIIIII")
}

///|
test "sqt_fastq_parse_multiple" {
  let content = "@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nIIII\n"
  let reads = sqt_fastq_parse(content)
  assert_eq(reads.length(), 2)
  assert_eq(reads[0].id, "r1")
  assert_eq(reads[0].sequence, "ACGT")
  assert_eq(reads[1].id, "r2")
  assert_eq(reads[1].sequence, "TGCA")
}

///|
test "sqt_fastq_parse_empty" {
  let content = ""
  let reads = sqt_fastq_parse(content)
  assert_eq(reads.length(), 0)
}

///|
test "sqt_fastq_parse_malformed" {
  let content = "not a fastq file\n"
  let reads = sqt_fastq_parse(content)
  assert_eq(reads.length(), 0)
}

///|
test "sqt_fastq_serialize_basic" {
  let reads = [
    FastqRead::{
      id: "test1",
      sequence: "ACGT",
      quality: "IIII",
    },
  ]
  let output = sqt_fastq_serialize(reads)
  assert_eq(output, "@test1\nACGT\n+\nIIII\n")
}

///|
test "sqt_fastq_serialize_multiple" {
  let reads = [
    FastqRead::{
      id: "r1",
      sequence: "AAAA",
      quality: "IIII",
    },
    FastqRead::{
      id: "r2",
      sequence: "CCCC",
      quality: "IIII",
    },
  ]
  let output = sqt_fastq_serialize(reads)
  assert_eq(
    output,
    "@r1\nAAAA\n+\nIIII\n@r2\nCCCC\n+\nIIII\n",
  )
}

///|
test "sqt_fastq_roundtrip" {
  let original = "@sample1 desc\nACGTACGT\n+\nIIIIIIII\n@sample2\nTGCA\n+\nIIII\n"
  let reads = sqt_fastq_parse(original)
  assert_eq(reads.length(), 2)
  let serialized = sqt_fastq_serialize(reads)
  let reparsed = sqt_fastq_parse(serialized)
  assert_eq(reparsed.length(), 2)
  assert_eq(reparsed[0].id, "sample1")
  assert_eq(reparsed[0].sequence, "ACGTACGT")
  assert_eq(reparsed[1].id, "sample2")
  assert_eq(reparsed[1].sequence, "TGCA")
}

///|
test "sqt_fastq_parse_with_description" {
  let content = "@read1 some description here\nACGT\n+\nIIII\n"
  let reads = sqt_fastq_parse(content)
  assert_eq(reads.length(), 1)
  assert_eq(reads[0].id, "read1")
}

///|
test "sqt_quality_trim_single_base" {
  let seq = "A"
  let qual = "I"
  let (t_seq, _) = sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq, "A")
}

///|
test "sqt_compute_stats_gc_content" {
  let reads = [
    FastqRead::{
      id: "gc_high",
      sequence: "GGGGCCCC",
      quality: "IIIIIIII",
    },
    FastqRead::{
      id: "at_high",
      sequence: "AAAATTTT",
      quality: "IIIIIIII",
    },
  ]
  let results = [
    TrimResult::{
      read_id: "gc_high",
      original_seq: "GGGGCCCC",
      trimmed_seq: "GGGGCCCC",
      original_qual: "IIIIIIII",
      trimmed_qual: "IIIIIIII",
      trim_type: "none",
    },
    TrimResult::{
      read_id: "at_high",
      original_seq: "AAAATTTT",
      trimmed_seq: "AAAATTTT",
      original_qual: "IIIIIIII",
      trimmed_qual: "IIIIIIII",
      trim_type: "none",
    },
  ]
  let stats = sqt_compute_stats(results, reads)
  assert_true(stats.gc_content_before > 0.4)
  assert_true(stats.gc_content_before < 0.6)
}

///|
test "sqt_trim_reads_quality_trim_only" {
  let reads = [
    FastqRead::{
      id: "r1",
      sequence: "ACGTACGT",
      quality: "IIII####",
    },
  ]
  let qual_config = QualityConfig::{
    min_quality: 20,
    min_length: 4,
    window_size: 4,
    avg_quality: 20,
  }
  let results = sqt_trim_reads(reads, qual_config, None)
  assert_eq(results.length(), 1)
  assert_eq(results[0].trimmed_seq.length(), 4)
}

///|
test "sqt_trim_adapter_partial_overlap" {
  let seq = "ACGTACGTAGATCGGAA"
  let qual = "IIIIIIIIIIIIIIII"
  let adapter = "AGATCGGAA"
  let (t_seq, _) = sqt_trim_adapter(
    seq, qual, adapter, allowed_mismatches=1,
  )
  assert_eq(t_seq, "ACGTACGT")
}

///|
test "sqt_quality_trim_quality_string_preserved" {
  let seq = "ACGTACGTACGT"
  let qual = "IIIIIIIIII##"
  let (t_seq, t_qual) = sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq.length(), t_qual.length())
  assert_eq(t_qual, "IIIIIIIIII")
}

///|
test "sqt_trim_reads_adapter_at_end" {
  let reads = [
    FastqRead::{
      id: "r1",
      sequence: "ACGTACGTACGTAGATCGGAA",
      quality: "IIIIIIIIIIIIIIIIIIIII",
    },
  ]
  let qual_config = QualityConfig::{
    min_quality: 20,
    min_length: 4,
    window_size: 4,
    avg_quality: 20,
  }
  let adapter_config = AdapterConfig::{
    adapters: ["AGATCGGAA"],
    min_overlap: 4,
    allowed_mismatches: 0,
  }
  let results = sqt_trim_reads(reads, qual_config, Some(adapter_config))
  assert_eq(results.length(), 1)
  assert_eq(results[0].trimmed_seq, "ACGTACGTACGT")
}

///|
test "sqt_fastq_parse_crlf" {
  let content = "@r1 desc\r\nACGT\r\n+\r\nIIII\r\n"
  let reads = sqt_fastq_parse(content)
  assert_eq(reads.length(), 1)
  assert_eq(reads[0].id, "r1")
  assert_eq(reads[0].sequence, "ACGT")
}

///|
test "sqt_compute_stats_all_discarded" {
  let reads = [
    FastqRead::{
      id: "r1",
      sequence: "ACGT",
      quality: "IIII",
    },
  ]
  let results = [
    TrimResult::{
      read_id: "r1",
      original_seq: "ACGT",
      trimmed_seq: "",
      original_qual: "IIII",
      trimmed_qual: "",
      trim_type: "quality",
    },
  ]
  let stats = sqt_compute_stats(results, reads)
  assert_eq(stats.total_reads, 1)
  assert_eq(stats.kept_reads, 0)
  assert_eq(stats.discarded_reads, 1)
  assert_eq(stats.total_bp_after, 0)
  assert_eq(stats.gc_content_after, 0.0)
}