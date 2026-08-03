///|
/// FIMO-like (Find Individual Motif Occurrences) motif scanning module.
///
/// Inspired by the MEME suite's FIMO tool. Scans biological sequences for
/// occurrences of known motifs represented as Position Weight Matrices (PWMs).
/// Computes log-likelihood scores and p-values for each match position using a
/// dynamic-programming null distribution.
///
/// Base indexing convention: A=0, C=1, G=2, T=3.
/// Matrix orientation: 4 rows (A, C, G, T) x L columns (positions).
///
/// Note: This module defines `MotifScanPSSM` (rather than `PSSM`) and
/// `motif_scan_pwm_from_counts` (rather than `pwm_from_counts`) to avoid name
/// collisions with the pre-existing `PSSM` type alias and `pwm_from_counts`
/// function in `motifs.mbt`.

///|
/// A Position Weight Matrix for a DNA motif.
pub struct PositionWeightMatrix {
  motif_id : String
  motif_name : String
  matrix : Array[Array[Double]] // 4 rows (A, C, G, T), L columns
  background : Array[Double] // background frequencies [A, C, G, T]
  strand : String // "+" or "+-" (both strands)
}

///|
/// A Position-Specific Scoring Matrix (log-likelihood scores).
pub struct MotifScanPSSM {
  motif_id : String
  motif_name : String
  scores : Array[Array[Double]] // 4 rows (A, C, G, T), L columns
  min_score : Double
  max_score : Double
  background : Array[Double] // background used for p-value null distribution
  strand : String // strand mode ("+" or "+-")
}

///|
/// A single motif occurrence reported by the scanner.
pub struct MotifMatch {
  motif_id : String
  motif_name : String
  sequence_id : String
  start : Int // 0-based start position
  end : Int // end position (exclusive)
  strand : String // "+" or "-"
  score : Double // log-likelihood score
  p_value : Double
  matched_sequence : String
}

///|
/// Construct a PositionWeightMatrix with default background and strand.
pub fn PositionWeightMatrix::new(
  motif_id~ : String,
  motif_name~ : String,
  matrix~ : Array[Array[Double]],
  background? : Array[Double] = [0.25, 0.25, 0.25, 0.25],
  strand? : String = "+-",
) -> PositionWeightMatrix {
  PositionWeightMatrix::{ motif_id, motif_name, matrix, background, strand }
}

///|
/// Construct a MotifMatch.
pub fn MotifMatch::new(
  motif_id~ : String,
  motif_name~ : String,
  sequence_id~ : String,
  start~ : Int,
  end~ : Int,
  strand~ : String,
  score~ : Double,
  p_value~ : Double,
  matched_sequence~ : String,
) -> MotifMatch {
  MotifMatch::{
    motif_id,
    motif_name,
    sequence_id,
    start,
    end,
    strand,
    score,
    p_value,
    matched_sequence,
  }
}

///|
/// Build a PositionWeightMatrix from a count matrix.
/// `counts` is 4 rows (A, C, G, T) x L columns of integer frequencies.
/// A pseudocount is added to avoid zero probabilities.
pub fn motif_scan_pwm_from_counts(
  motif_id : String,
  motif_name : String,
  counts : Array[Array[Int]],
  pseudocount : Double,
) -> PositionWeightMatrix {
  let n_bases = counts.length()
  if n_bases == 0 {
    return PositionWeightMatrix::new(
      motif_id=motif_id,
      motif_name=motif_name,
      matrix=[],
      background=[0.25, 0.25, 0.25, 0.25],
      strand="+-",
    )
  }
  let motif_len = counts[0].length()
  let matrix = ms_make_2d(n_bases, motif_len, 0.0)
  for i in 0..<motif_len {
    let mut total = 0.0
    for b in 0..<n_bases {
      total = total + counts[b][i].to_double()
    }
    let norm = total + pseudocount * n_bases.to_double()
    for b in 0..<n_bases {
      matrix[b][i] = (counts[b][i].to_double() + pseudocount) / norm
    }
  }
  PositionWeightMatrix::new(
    motif_id=motif_id,
    motif_name=motif_name,
    matrix=matrix,
    background=[0.25, 0.25, 0.25, 0.25],
    strand="+-",
  )
}

///|
/// Convert a PWM to a PSSM of log-likelihood scores.
/// For each base b and position i: PSSM[b][i] = log2(PWM[b][i] / background[b]).
/// Zero PWM entries use a small pseudocount (0.001).
/// min_score = sum of per-position minima; max_score = sum of per-position maxima.
pub fn build_pssm(pwm : PositionWeightMatrix) -> MotifScanPSSM {
  let n_bases = pwm.matrix.length()
  if n_bases == 0 {
    return MotifScanPSSM::{
      motif_id: pwm.motif_id,
      motif_name: pwm.motif_name,
      scores: [],
      min_score: 0.0,
      max_score: 0.0,
      background: pwm.background,
      strand: pwm.strand,
    }
  }
  let motif_len = pwm.matrix[0].length()
  let scores = ms_make_2d(n_bases, motif_len, 0.0)
  let pseudocount = 0.001
  let mut min_score = 0.0
  let mut max_score = 0.0
  for i in 0..<motif_len {
    let mut col_min = 0.0
    let mut col_max = 0.0
    for b in 0..<n_bases {
      let p = pwm.matrix[b][i]
      let p_eff = if p > 0.0 { p } else { pseudocount }
      let bg = pwm.background[b]
      let bg_eff = if bg > 0.0 { bg } else { 0.25 }
      let ll = ms_log2(p_eff / bg_eff)
      scores[b][i] = ll
      if b == 0 {
        col_min = ll
        col_max = ll
      } else {
        if ll < col_min {
          col_min = ll
        }
        if ll > col_max {
          col_max = ll
        }
      }
    }
    min_score = min_score + col_min
    max_score = max_score + col_max
  }
  MotifScanPSSM::{
    motif_id: pwm.motif_id,
    motif_name: pwm.motif_name,
    scores,
    min_score,
    max_score,
    background: pwm.background,
    strand: pwm.strand,
  }
}

///|
/// Score a sequence against a PSSM by summing per-position log-likelihoods.
/// Scores min(sequence.length, L) positions. Non-ACGT characters are skipped.
pub fn compute_score(pssm : MotifScanPSSM, sequence : String) -> Double {
  let n_bases = pssm.scores.length()
  if n_bases == 0 {
    return 0.0
  }
  let motif_len = pssm.scores[0].length()
  if motif_len == 0 {
    return 0.0
  }
  let n = if sequence.length() < motif_len { sequence.length() } else { motif_len }
  let mut score = 0.0
  for i in 0..<n {
    let c = sequence[i:i + 1].to_owned()
    let b = ms_base_index(c)
    if b >= 0 && b < n_bases {
      score = score + pssm.scores[b][i]
    }
  }
  score
}

///|
/// Compute the p-value P(score >= observed) under the null model (background
/// frequencies) using a dynamic-programming convolution of per-position score
/// distributions. Scores are discretized into 1000 bins over [min_score, max_score].
pub fn compute_pvalue(pssm : MotifScanPSSM, score : Double) -> Double {
  let n_bases = pssm.scores.length()
  if n_bases == 0 {
    return 1.0
  }
  let motif_len = pssm.scores[0].length()
  if motif_len == 0 {
    return 1.0
  }
  let min_s = pssm.min_score
  let max_s = pssm.max_score
  if score > max_s {
    return 0.0
  }
  if score <= min_s {
    return 1.0
  }
  let range = max_s - min_s
  if range <= 0.0 {
    return if score <= min_s { 1.0 } else { 0.0 }
  }
  let nbins = 1000
  let bin_width = range / nbins.to_double()
  let mut dist = Array::make(nbins, 0.0)
  // Initialize: before any position the cumulative score is 0.0 with probability 1.
  let init_bin = ms_clamp_bin(((0.0 - min_s) / bin_width).floor().to_int(), nbins)
  dist[init_bin] = 1.0
  // Convolve the per-position score distributions.
  for i in 0..<motif_len {
    let new_dist = Array::make(nbins, 0.0)
    for b_idx in 0..<nbins {
      let p = dist[b_idx]
      if p > 0.0 {
        let bin_score = min_s + (b_idx.to_double() + 0.5) * bin_width
        for base in 0..<n_bases {
          let bg = pssm.background[base]
          if bg > 0.0 {
            let new_score = bin_score + pssm.scores[base][i]
            let new_bin = ms_clamp_bin(
              ((new_score - min_s) / bin_width).floor().to_int(),
              nbins,
            )
            new_dist[new_bin] = new_dist[new_bin] + p * bg
          }
        }
      }
    }
    dist = new_dist
  }
  // Sum probabilities for scores >= observed (from observed's bin upward).
  let obs_bin = ms_clamp_bin(((score - min_s) / bin_width).floor().to_int(), nbins)
  let mut pval = 0.0
  for b_idx in obs_bin..<nbins {
    pval = pval + dist[b_idx]
  }
  if pval > 1.0 {
    pval = 1.0
  }
  if pval < 0.0 {
    pval = 0.0
  }
  pval
}

///|
/// Scan a single sequence (forward and, if strand is "+-", reverse strand) for
/// motif occurrences with p-value at or below the threshold.
pub fn scan_sequence(
  pssm : MotifScanPSSM,
  sequence : String,
  seq_id : String,
  pvalue_threshold : Double,
) -> Array[MotifMatch] {
  let matches : Array[MotifMatch] = Array::new()
  let n_bases = pssm.scores.length()
  if n_bases == 0 {
    return matches
  }
  let motif_len = pssm.scores[0].length()
  if motif_len == 0 || sequence.length() < motif_len {
    return matches
  }
  let n_pos = sequence.length() - motif_len + 1
  // Forward strand
  for start in 0..<n_pos {
    let subseq = sequence[start:start + motif_len].to_owned()
    let score_val = compute_score(pssm, subseq)
    let pval = compute_pvalue(pssm, score_val)
    if pval <= pvalue_threshold {
      matches.push(MotifMatch::new(
        motif_id=pssm.motif_id,
        motif_name=pssm.motif_name,
        sequence_id=seq_id,
        start=start,
        end=start + motif_len,
        strand="+",
        score=score_val,
        p_value=pval,
        matched_sequence=subseq,
      ))
    }
  }
  // Reverse strand
  if pssm.strand == "+-" {
    let rev_seq = reverse_complement_sequence(sequence)
    for start in 0..<n_pos {
      let subseq = rev_seq[start:start + motif_len].to_owned()
      let score_val = compute_score(pssm, subseq)
      let pval = compute_pvalue(pssm, score_val)
      if pval <= pvalue_threshold {
        let orig_start = sequence.length() - start - motif_len
        matches.push(MotifMatch::new(
          motif_id=pssm.motif_id,
          motif_name=pssm.motif_name,
          sequence_id=seq_id,
          start=orig_start,
          end=orig_start + motif_len,
          strand="-",
          score=score_val,
          p_value=pval,
          matched_sequence=subseq,
        ))
      }
    }
  }
  matches
}

///|
/// Scan multiple sequences. `sequences` is an array of (sequence_id, sequence) tuples.
pub fn scan_sequences(
  pssm : MotifScanPSSM,
  sequences : Array[(String, String)],
  pvalue_threshold : Double,
) -> Array[MotifMatch] {
  let result : Array[MotifMatch] = Array::new()
  for i in 0..<sequences.length() {
    let pair = sequences[i]
    let seq_id = pair.0
    let seq = pair.1
    let matches = scan_sequence(pssm, seq, seq_id, pvalue_threshold)
    for j in 0..<matches.length() {
      result.push(matches[j])
    }
  }
  result
}

///|
/// Compute the reverse complement of a DNA sequence.
pub fn reverse_complement_sequence(seq : String) -> String {
  let mut result = ""
  for i = seq.length() - 1; i >= 0; i = i - 1 {
    let c = seq[i:i + 1].to_owned().to_upper()
    let comp = match c {
      "A" => "T"
      "T" => "A"
      "C" => "G"
      "G" => "C"
      _ => c
    }
    result = result + comp
  }
  result
}

///|
/// Filter matches by p-value threshold (keeps matches with p_value <= threshold).
pub fn filter_matches(
  matches : Array[MotifMatch],
  pvalue_threshold : Double,
) -> Array[MotifMatch] {
  let result : Array[MotifMatch] = Array::new()
  for i in 0..<matches.length() {
    if matches[i].p_value <= pvalue_threshold {
      result.push(matches[i])
    }
  }
  result
}

///|
/// Produce a human-readable summary of motif scan results.
pub fn motif_scan_summary(matches : Array[MotifMatch]) -> String {
  let mut s = "=== FIMO Motif Scan Summary ===\n"
  s = s + "Total matches: " + matches.length().to_string() + "\n"
  if matches.length() > 0 {
    s = s + "\nMatches:\n"
    s = s + "motif_id\tsequence_id\tstart\tend\tstrand\tscore\tp-value\tmatched_sequence\n"
    for i in 0..<matches.length() {
      let m = matches[i]
      s = s +
        m.motif_id +
        "\t" +
        m.sequence_id +
        "\t" +
        m.start.to_string() +
        "\t" +
        m.end.to_string() +
        "\t" +
        m.strand +
        "\t" +
        ms_fmt_double(m.score) +
        "\t" +
        ms_fmt_double(m.p_value) +
        "\t" +
        m.matched_sequence +
        "\n"
    }
  }
  s
}

///|
/// A sample PWM (6-position motif, consensus ATGCAA) for testing and demos.
pub fn motif_scan_sample_pwm() -> PositionWeightMatrix {
  // 4 rows (A, C, G, T) x 6 columns. Consensus: A T G C A A
  let matrix = [
    [0.9, 0.05, 0.05, 0.05, 0.9, 0.9],
    [0.05, 0.05, 0.05, 0.9, 0.05, 0.05],
    [0.05, 0.05, 0.9, 0.05, 0.05, 0.05],
    [0.05, 0.9, 0.05, 0.05, 0.05, 0.05],
  ]
  PositionWeightMatrix::new(
    motif_id="M0001",
    motif_name="SampleTF",
    matrix=matrix,
    background=[0.25, 0.25, 0.25, 0.25],
    strand="+-",
  )
}

///|
/// Sample sequences for testing and demos.
pub fn motif_scan_sample_sequences() -> Array[(String, String)] {
  [
    ("seq1", "GGGGATGCAAGGGG"), // forward ATGCAA at position 4
    ("seq2", "GGGGTTGCATGGGG"), // reverse-strand match (TTGCAT = revcomp of ATGCAA)
    ("seq3", "ATGCAAATGCAA"), // two forward matches at 0 and 6
  ]
}

// --- Internal helpers ---

///|
/// Compute log2(x) using natural log.
fn ms_log2(x : Double) -> Double {
  @math.ln(x) / @math.ln(2.0)
}

///|
/// Create a 2D array with distinct inner arrays (Array::make shares references
/// for mutable element types, so each row must be allocated separately).
fn ms_make_2d(rows : Int, cols : Int, init : Double) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  for _ in 0..<rows {
    result.push(Array::make(cols, init))
  }
  result
}

///|
/// Map a DNA base character to its index (A=0, C=1, G=2, T=3). Returns -1 for non-ACGT.
fn ms_base_index(c : String) -> Int {
  let upper = c.to_upper()
  match upper {
    "A" => 0
    "C" => 1
    "G" => 2
    "T" => 3
    _ => -1
  }
}

///|
/// Clamp a bin index to [0, nbins - 1].
fn ms_clamp_bin(idx : Int, nbins : Int) -> Int {
  if idx < 0 {
    0
  } else if idx >= nbins {
    nbins - 1
  } else {
    idx
  }
}

///|
/// Format a double to a fixed number of decimal places (6) for display.
fn ms_fmt_double(x : Double) -> String {
  let neg = x < 0.0
  let abs_x = if neg { -x } else { x }
  let scaled = (abs_x * 1000000.0 + 0.5).floor()
  let int_part = (scaled / 1000000.0).floor().to_int()
  let frac_int = (scaled - int_part.to_double() * 1000000.0).to_int()
  let frac_str = if frac_int < 10 {
    "00000" + frac_int.to_string()
  } else if frac_int < 100 {
    "0000" + frac_int.to_string()
  } else if frac_int < 1000 {
    "000" + frac_int.to_string()
  } else if frac_int < 10000 {
    "00" + frac_int.to_string()
  } else if frac_int < 100000 {
    "0" + frac_int.to_string()
  } else {
    frac_int.to_string()
  }
  let sign = if neg { "-" } else { "" }
  sign + int_part.to_string() + "." + frac_str
}
