///|
/// SeqKit `stats` port — the `seqkit stats -a` sequence statistics. Computes
/// length stats (num/sum/min/avg/max/quartiles), N50/L50, GC%, gap count, N
/// count, and (for FASTQ) Q20%/Q30%/average quality, from already-parsed
/// `SeqRecord`s. Semantics cross-checked against seqkit `cmd/stat.go` + the
/// `bio/util.LengthStats` helpers.

///|
pub struct SeqkitStats {
  num_seqs : Int
  sum_len : Int
  min_len : Int
  avg_len : Double
  max_len : Int
  q1 : Int
  q2 : Int
  q3 : Int
  sum_gap : Int
  n50 : Int
  l50 : Int
  q20_pct : Double
  q30_pct : Double
  avg_qual : Double
  gc_pct : Double
  sum_n : Int
}

///|
fn stat_sort_asc(arr : Array[Int]) -> Array[Int] {
  let out = arr.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && out[j] > key {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
/// Round to `n` decimal places (half-up, values here are non-negative).
fn round_n(x : Double, n : Int) -> Double {
  let f = @math.pow(10.0, n.to_double())
  @math.floor(x * f + 0.5) / f
}

///|
fn log10(x : Double) -> Double {
  @math.ln(x) / @math.ln(10.0)
}

///|
/// Median of `sorted[lo:hi)` (sorted ascending); 0 for an empty slice.
fn stat_median_range(sorted : Array[Int], lo : Int, hi : Int) -> Int {
  let l = hi - lo
  if l == 0 {
    return 0
  }
  if l % 2 == 0 {
    (sorted[lo + l / 2 - 1] + sorted[lo + l / 2]) / 2
  } else {
    sorted[lo + l / 2]
  }
}

///|
/// SeqKit `stats -a` for one set of records. `is_nucleotide` selects N/n
/// (DNA/RNA) vs X/x (protein) for the `sum_n` column; `is_fastq` enables the
/// Q20/Q30/AvgQual quality columns (read from the `phred_quality` letter
/// annotation).
pub fn seqkit_stats(
  records : Array[SeqRecord],
  is_nucleotide : Bool,
  is_fastq : Bool,
) -> SeqkitStats {
  let lengths : Array[Int] = Array(capacity=0)
  let mut sum_len = 0
  let mut gc = 0
  let mut gap = 0
  let mut n = 0
  let mut q20 = 0
  let mut q30 = 0
  let mut err_sum = 0.0
  for rec in records {
    let s = rec.seq.to_string()
    let len = s.length()
    lengths.push(len)
    sum_len = sum_len + len
    let mut i = 0
    while i < len {
      let c = char_at(s, i)
      if c == "G" || c == "C" || c == "g" || c == "c" {
        gc = gc + 1
      }
      if c == "-" || c == " " || c == "." {
        gap = gap + 1
      }
      if is_nucleotide {
        if c == "N" || c == "n" {
          n = n + 1
        }
      } else {
        if c == "X" || c == "x" {
          n = n + 1
        }
      }
      i = i + 1
    }
    if is_fastq {
      match rec.letter_annotations.get("phred_quality") {
        None => ()
        Some(quals) => {
          for q in quals {
            if q >= 20 {
              q20 = q20 + 1
              if q >= 30 {
                q30 = q30 + 1
              }
            }
            err_sum = err_sum + @math.pow(10.0, (0.0 - q.to_double()) / 10.0)
          }
        }
      }
    }
  }

  let num = lengths.length()
  if num == 0 {
    return SeqkitStats::{
      num_seqs: 0,
      sum_len: 0,
      min_len: 0,
      avg_len: 0.0,
      max_len: 0,
      q1: 0,
      q2: 0,
      q3: 0,
      sum_gap: 0,
      n50: 0,
      l50: 0,
      q20_pct: 0.0,
      q30_pct: 0.0,
      avg_qual: 0.0,
      gc_pct: 0.0,
      sum_n: 0,
    }
  }

  let sorted = stat_sort_asc(lengths)
  let min_len = sorted[0]
  let max_len = sorted[num - 1]
  let avg_len = round_n(sum_len.to_double() / num.to_double(), 1)

  // quartiles (integer, truncated), matching seqkit's `quartile`
  let c1 = if num % 2 == 0 { num / 2 } else { (num - 1) / 2 }
  let c2 = if num % 2 == 0 { num / 2 } else { c1 + 1 }
  let q1 = stat_median_range(sorted, 0, c1)
  let q2 = stat_median_range(sorted, 0, num)
  let q3 = stat_median_range(sorted, c2, num)

  // N50 / L50 over lengths descending
  let target = sum_len / 2 + sum_len % 2
  let mut cum = 0
  let mut n50 = 0
  let mut l50 = 0
  let mut i = num - 1
  while i >= 0 {
    cum = cum + sorted[i]
    l50 = l50 + 1
    if cum >= target {
      n50 = sorted[i]
      break
    }
    i = i - 1
  }

  let gc_pct = round_n(gc.to_double() / sum_len.to_double() * 100.0, 2)
  let q20_pct = round_n(q20.to_double() / sum_len.to_double() * 100.0, 2)
  let q30_pct = round_n(q30.to_double() / sum_len.to_double() * 100.0, 2)
  let avg_qual = if err_sum == 0.0 {
    0.0
  } else {
    round_n(-10.0 * log10(err_sum / sum_len.to_double()), 2)
  }

  SeqkitStats::{
    num_seqs: num,
    sum_len,
    min_len,
    avg_len,
    max_len,
    q1,
    q2,
    q3,
    sum_gap: gap,
    n50,
    l50,
    q20_pct,
    q30_pct,
    avg_qual,
    gc_pct,
    sum_n: n,
  }
}