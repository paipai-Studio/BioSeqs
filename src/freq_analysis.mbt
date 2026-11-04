///|
/// Bio.FreqAnalysis module - Frequency analysis for sequences.
/// Inspired by Biopython's Bio.FreqAnalysis module.
/// Provides frequency counting, pattern matching, and statistical analysis.

///|
/// Represents a frequency analysis result.
pub struct FreqAnalysisResult {
  frequencies : Map[String, Double]
  counts : Map[String, Int32]
  total_count : Int32
}

///|
/// Create a new frequency analysis result.
pub fn FreqAnalysisResult::new(
  frequencies : Map[String, Double],
  counts : Map[String, Int32],
  total_count : Int32,
) -> FreqAnalysisResult {
  FreqAnalysisResult::{ frequencies, counts, total_count, }
}

///|
/// Get frequency for a given pattern.
pub fn FreqAnalysisResult::fa_get_frequency(
  self : FreqAnalysisResult,
  pattern : String,
) -> Double {
  self.frequencies.get_or_default(pattern, 0.0)
}

///|
/// Get count for a given pattern.
pub fn FreqAnalysisResult::fa_get_count(
  self : FreqAnalysisResult,
  pattern : String,
) -> Int32 {
  self.counts.get_or_default(pattern, 0)
}

///|
/// Get total count.
pub fn FreqAnalysisResult::fa_get_total_count(
  self : FreqAnalysisResult,
) -> Int32 {
  self.total_count
}

///|
/// Get all patterns.
pub fn FreqAnalysisResult::fa_get_patterns(
  self : FreqAnalysisResult,
) -> Array[String] {
  let patterns : Array[String] = Array(capacity=0)
  for key in self.frequencies.keys() {
    patterns.push(key)
  }
  patterns
}

///|
/// Count occurrences of a pattern in a sequence string.
pub fn fa_count_pattern(sequence : String, pattern : String) -> Int32 {
  let pattern_len = pattern.length()
  if pattern_len == 0 || sequence.length() < pattern_len {
    return 0
  }

  let mut count = 0
  let mut pos = 0
  while pos <= sequence.length() - pattern_len {
    let substr = substring(sequence, pos, pattern_len)
    if substr == pattern {
      count = count + 1
    }
    pos = pos + 1
  }
  count
}

///|
/// Count k-mer frequencies in a sequence string.
pub fn fa_count_kmers(sequence : String, k : Int32) -> FreqAnalysisResult {
  let counts : Map[String, Int32] = Map([], capacity=100)
  let seq_len = sequence.length()

  if seq_len < k {
    return FreqAnalysisResult::new(Map([], capacity=0), Map([], capacity=0), 0)
  }

  let mut total = 0
  let mut pos = 0
  while pos <= seq_len - k {
    let kmer = substring(sequence, pos, k)
    let current = counts.get_or_default(kmer, 0)
    counts.set(kmer, current + 1)
    total = total + 1
    pos = pos + 1
  }

  let frequencies : Map[String, Double] = Map([], capacity=counts.length())
  for key in counts.keys() {
    let count = counts.get_or_default(key, 0)
    frequencies.set(key, count.to_double() / total.to_double())
  }

  FreqAnalysisResult::new(frequencies, counts, total)
}

///|
/// Count nucleotide frequencies in a DNA sequence.
pub fn fa_nucleotide_frequency(sequence : String) -> FreqAnalysisResult {
  fa_count_kmers(sequence, 1)
}

///|
/// Count dinucleotide frequencies.
pub fn fa_dinucleotide_frequency(sequence : String) -> FreqAnalysisResult {
  fa_count_kmers(sequence, 2)
}

///|
/// Calculate expected frequency from mononucleotide frequencies.
/// Expected dinucleotide frequency = freq(A) * freq(B)
pub fn fa_expected_dinucleotide_freq(
  mono_freqs : FreqAnalysisResult,
  dinuc : String,
) -> Double {
  if dinuc.length() != 2 {
    return 0.0
  }
  let base1 = substring(dinuc, 0, 1)
  let base2 = substring(dinuc, 1, 1)
  let f1 = mono_freqs.fa_get_frequency(base1)
  let f2 = mono_freqs.fa_get_frequency(base2)
  f1 * f2
}

///|
/// Count codon usage frequencies in a coding sequence.
pub fn fa_codon_usage(sequence : String) -> FreqAnalysisResult {
  let counts : Map[String, Int32] = Map([], capacity=64)
  let seq_len = sequence.length()

  if seq_len < 3 {
    return FreqAnalysisResult::new(Map([], capacity=0), Map([], capacity=0), 0)
  }

  let mut total = 0
  let mut pos = 0
  while pos <= seq_len - 3 {
    let codon = substring(sequence, pos, 3)
    let current = counts.get_or_default(codon, 0)
    counts.set(codon, current + 1)
    total = total + 1
    pos = pos + 3
  }

  let frequencies : Map[String, Double] = Map([], capacity=counts.length())
  for key in counts.keys() {
    let count = counts.get_or_default(key, 0)
    frequencies.set(key, count.to_double() / total.to_double())
  }

  FreqAnalysisResult::new(frequencies, counts, total)
}

///|
/// Find all occurrences of a motif in a sequence.
pub fn fa_find_motif(sequence : String, motif : String) -> Array[Int32] {
  let positions : Array[Int32] = Array(capacity=0)
  let motif_len = motif.length()
  let seq_len = sequence.length()

  if motif_len == 0 || seq_len < motif_len {
    return positions
  }

  let mut pos = 0
  while pos <= seq_len - motif_len {
    let substr = substring(sequence, pos, motif_len)
    if substr == motif {
      positions.push(pos)
    }
    pos = pos + 1
  }

  positions
}

///|
/// Calculate GC content (guanine and cytosine frequency).
pub fn fa_gc_content(sequence : String) -> Double {
  let counts = fa_nucleotide_frequency(sequence)
  let freq_g = counts.fa_get_frequency("G")
  let freq_c = counts.fa_get_frequency("C")
  freq_g + freq_c
}

///|
/// Calculate AT content (adenine and thymine frequency).
pub fn fa_at_content(sequence : String) -> Double {
  let counts = fa_nucleotide_frequency(sequence)
  let freq_a = counts.fa_get_frequency("A")
  let freq_t = counts.fa_get_frequency("T")
  freq_a + freq_t
}

///|
/// Calculate sequence complexity (Shannon entropy).
pub fn fa_sequence_complexity(sequence : String) -> Double {
  let counts = fa_count_kmers(sequence, 1)
  let mut entropy = 0.0
  let patterns = counts.fa_get_patterns()

  let mut i = 0
  while i < patterns.length() {
    let freq = counts.fa_get_frequency(patterns[i])
    if freq > 0.0 {
      entropy = entropy - freq * @math.ln(freq)
    }
    i = i + 1
  }
  entropy
}
