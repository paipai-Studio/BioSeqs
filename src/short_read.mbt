///|
pub struct FastqQualityStats {
  n_reads : Int
  read_length : Int
  mut gc_content : Double
  mut mean_quality : Double
  quality_by_cycle : Array[Double]
  n_count_by_cycle : Array[Int]
  quality_distribution : Array[Int]
  mut dup_rate : Double
}

///|
pub fn FastqQualityStats::new(
  n_reads : Int,
  read_length : Int,
) -> FastqQualityStats {
  let quality_by_cycle : Array[Double] = Array::make(read_length, 0.0)
  let n_count_by_cycle : Array[Int] = Array::make(read_length, 0)
  let quality_distribution : Array[Int] = Array::make(42, 0)
  FastqQualityStats::{
    n_reads,
    read_length,
    gc_content: 0.0,
    mean_quality: 0.0,
    quality_by_cycle,
    n_count_by_cycle,
    quality_distribution,
    dup_rate: 0.0,
  }
}

///|
pub fn FastqQualityStats::to_string(self : FastqQualityStats) -> String {
  let mut result = "FastqQualityStats\n"
  result = result + "===================\n"
  result = result + "Total reads: " + self.n_reads.to_string() + "\n"
  result = result + "Read length: " + self.read_length.to_string() + "\n"
  result = result +
    "GC content: " +
    (self.gc_content * 100.0).to_string() +
    "%\n"
  result = result + "Mean quality: " + self.mean_quality.to_string() + "\n"
  result = result +
    "Duplication rate: " +
    (self.dup_rate * 100.0).to_string() +
    "%\n"
  result
}

///|
fn qual_char_to_phred(c : UInt16) -> Int {
  c.to_int() - 33
}

///|
fn is_gc(base : UInt16) -> Bool {
  base == 71 || base == 67 || base == 103 || base == 99
}

///|
fn is_n(base : UInt16) -> Bool {
  base == 78 || base == 110
}

///|
pub fn qa_summary(
  sequences : Array[String],
  qualities : Array[String],
) -> FastqQualityStats {
  let n = sequences.length()
  if n == 0 {
    return FastqQualityStats::new(0, 0)
  }

  let read_length = sequences[0].length()
  let stats = FastqQualityStats::new(n, read_length)

  let cycle_qual_sum : Array[Double] = Array::make(read_length, 0.0)
  let cycle_n_count : Array[Int] = Array::make(read_length, 0)
  let qual_dist : Array[Int] = Array::make(42, 0)

  let mut total_gc = 0.0
  let mut total_bases = 0
  let mut total_qual = 0.0

  let seen : Map[String, Int] = Map([], capacity=n / 2)
  let mut dup_count = 0

  let mut i = 0
  while i < n {
    let seq = sequences[i]
    let qual = qualities[i]

    let count = seen.get(seq)
    if count is Some(_) {
      dup_count = dup_count + 1
      seen[seq] = count.unwrap() + 1
    } else {
      seen[seq] = 1
    }

    let mut j = 0
    while j < read_length && j < seq.length() && j < qual.length() {
      let base = seq.unsafe_get(j)
      let q = qual_char_to_phred(qual.unsafe_get(j))

      cycle_qual_sum[j] = cycle_qual_sum[j] + q.to_double()

      if is_n(base) {
        cycle_n_count[j] = cycle_n_count[j] + 1
      }

      if is_gc(base) {
        total_gc = total_gc + 1.0
      }

      if q >= 0 && q < 42 {
        qual_dist[q] = qual_dist[q] + 1
      }

      total_qual = total_qual + q.to_double()
      total_bases = total_bases + 1

      j = j + 1
    }

    i = i + 1
  }

  let mut k = 0
  while k < read_length {
    stats.quality_by_cycle[k] = cycle_qual_sum[k] / n.to_double()
    stats.n_count_by_cycle[k] = cycle_n_count[k]
    k = k + 1
  }

  let mut m = 0
  while m < 42 {
    stats.quality_distribution[m] = qual_dist[m]
    m = m + 1
  }

  stats.gc_content = if total_bases > 0 {
    total_gc / total_bases.to_double()
  } else {
    0.0
  }
  stats.mean_quality = if total_bases > 0 {
    total_qual / total_bases.to_double()
  } else {
    0.0
  }
  stats.dup_rate = if n > 0 {
    dup_count.to_double() / n.to_double()
  } else {
    0.0
  }

  stats
}

///|
pub struct TrimParams {
  min_quality : Int
  min_length : Int
  trim_left : Int
  trim_right : Int
}

///|
pub fn TrimParams::new() -> TrimParams {
  TrimParams::{ min_quality: 20, min_length: 36, trim_left: 0, trim_right: 0 }
}

///|
pub fn trim_tail(
  sequence : String,
  quality : String,
  min_quality : Int,
) -> (String, String) {
  let len = sequence.length()
  if len == 0 {
    return ("", "")
  }

  let mut trim_pos = 0
  let mut i = len - 1
  while i >= 0 {
    let q = qual_char_to_phred(quality.unsafe_get(i))
    if q >= min_quality {
      trim_pos = i + 1
      break
    }
    i = i - 1
  }

  (sequence[0:trim_pos].to_string(), quality[0:trim_pos].to_string())
}

///|
pub fn trim_leading(
  sequence : String,
  quality : String,
  min_quality : Int,
) -> (String, String) {
  let len = sequence.length()
  if len == 0 {
    return ("", "")
  }

  let mut start = 0
  while start < len {
    let q = qual_char_to_phred(quality.unsafe_get(start))
    if q >= min_quality {
      break
    }
    start = start + 1
  }

  if start >= len {
    return ("", "")
  }

  (sequence[start:len].to_string(), quality[start:len].to_string())
}

///|
pub fn trim_adapter(
  sequence : String,
  quality : String,
  adapter : String,
  max_mismatch : Int,
) -> (String, String) {
  let seq_len = sequence.length()
  let adapter_len = adapter.length()

  if adapter_len == 0 || seq_len < adapter_len {
    return (sequence, quality)
  }

  let mut best_pos = -1
  let mut best_score = -1

  let mut i = 0
  while i <= seq_len - adapter_len / 2 {
    let mut mismatches = 0
    let mut matched = 0
    let mut j = 0
    while j < adapter_len && i + j < seq_len {
      if sequence.unsafe_get(i + j) != adapter.unsafe_get(j) {
        mismatches = mismatches + 1
        if mismatches > max_mismatch {
          break
        }
      } else {
        matched = matched + 1
      }
      j = j + 1
    }

    if matched > best_score && mismatches <= max_mismatch {
      best_score = matched
      best_pos = i
    }

    i = i + 1
  }

  if best_pos >= 0 && best_score >= adapter_len / 2 {
    (sequence[0:best_pos].to_string(), quality[0:best_pos].to_string())
  } else {
    (sequence, quality)
  }
}

///|
pub fn filter_reads(
  sequences : Array[String],
  qualities : Array[String],
  min_quality : Double,
  min_length : Int,
  max_n : Int,
) -> (Array[String], Array[String]) {
  let filtered_seqs : Array[String] = Array::new()
  let filtered_quals : Array[String] = Array::new()

  let mut i = 0
  while i < sequences.length() {
    let seq = sequences[i]
    let qual = qualities[i]

    if seq.length() < min_length {
      i = i + 1
      continue
    }

    let mut total_qual = 0.0
    let mut n_count = 0
    let mut j = 0
    while j < seq.length() && j < qual.length() {
      total_qual = total_qual +
        qual_char_to_phred(qual.unsafe_get(j)).to_double()
      if is_n(seq.unsafe_get(j)) {
        n_count = n_count + 1
      }
      j = j + 1
    }

    let mean_q = if seq.length() > 0 {
      total_qual / seq.length().to_double()
    } else {
      0.0
    }

    if mean_q >= min_quality && n_count <= max_n {
      filtered_seqs.push(seq)
      filtered_quals.push(qual)
    }

    i = i + 1
  }

  (filtered_seqs, filtered_quals)
}

///|
pub struct FastQCReport {
  basic_stats : FastqQualityStats
  per_base_quality : Array[Double]
  per_sequence_quality : Array[Int]
  per_base_n_content : Array[Double]
  gc_distribution : Array[Int]
  status : String
}

///|
pub fn FastQCReport::new(stats : FastqQualityStats) -> FastQCReport {
  let n = stats.n_reads
  let read_len = stats.read_length

  let per_base_n : Array[Double] = Array::make(read_len, 0.0)
  let mut i = 0
  while i < read_len {
    per_base_n[i] = if n > 0 {
      stats.n_count_by_cycle[i].to_double() / n.to_double() * 100.0
    } else {
      0.0
    }
    i = i + 1
  }

  let gc_dist : Array[Int] = Array::make(101, 0)

  let mut status = "PASS"
  if stats.mean_quality < 28.0 {
    status = "WARN"
  }
  if stats.mean_quality < 20.0 {
    status = "FAIL"
  }

  FastQCReport::{
    basic_stats: stats,
    per_base_quality: stats.quality_by_cycle,
    per_sequence_quality: stats.quality_distribution,
    per_base_n_content: per_base_n,
    gc_distribution: gc_dist,
    status,
  }
}

///|
pub fn FastQCReport::to_string(self : FastQCReport) -> String {
  let mut result = "FastQC Report\n"
  result = result + "=============\n"
  result = result + "Status: " + self.status + "\n\n"
  result = result + self.basic_stats.to_string() + "\n"

  result = result + "Per-base quality (first 5 cycles):\n"
  let mut i = 0
  while i < 5 && i < self.per_base_quality.length() {
    result = result +
      "  Cycle " +
      (i + 1).to_string() +
      ": " +
      self.per_base_quality[i].to_string() +
      "\n"
    i = i + 1
  }

  result
}

///|
pub fn generate_fastqc_report(
  sequences : Array[String],
  qualities : Array[String],
) -> FastQCReport {
  let stats = qa_summary(sequences, qualities)
  FastQCReport::new(stats)
}

///|
pub fn reads_to_sequences(
  records : Array[SeqRecord],
) -> (Array[String], Array[String]) {
  let seqs : Array[String] = Array::new()
  let quals : Array[String] = Array::new()

  let mut i = 0
  while i < records.length() {
    let record = records[i]
    seqs.push(record.seq.to_string())

    let qual = record.letter_annotations.get("phred_quality")
    if qual is Some(_) {
      let qual_arr = qual.unwrap()
      let mut qs = ""
      let mut j = 0
      while j < qual_arr.length() {
        qs = qs + (qual_arr[j] + 33).unsafe_to_char().to_string()
        j = j + 1
      }
      quals.push(qs)
    } else {
      let mut qs = ""
      let mut j = 0
      while j < record.seq.length() {
        qs = qs + 'I'.to_string()
        j = j + 1
      }
      quals.push(qs)
    }

    i = i + 1
  }

  (seqs, quals)
}
