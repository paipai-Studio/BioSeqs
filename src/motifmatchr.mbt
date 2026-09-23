///|
/// PWM motif matching against a set of DNA sequences, ported from
/// Bioconductor's `motifmatchr` package.
///
/// The reference package exposes `matchMotifs(motifs, subjects, ...)`
/// returning a `MotifMatches` object with `motifMatches()`,
/// `motifCounts()`, `motifScores()` and `maxMotifScores()` accessors. This
/// port keeps the same vocabulary: PFM → PWM → PSSM conversion, forward
/// and reverse-complement scanning, and per-position hit/score matrices.

///|
/// A named motif stored as a 4 x L position weight matrix (rows ordered
/// A, C, G, T). When the input is a position frequency matrix (PFM),
/// `motif_pwm_from_counts` converts it to a PWM using pseudocounts and a
/// background nucleotide distribution.
pub struct MatchrMotif {
  name : String
  pwm : Array[Array[Double]]
} derive(Eq)
pub extend MatchrMotif with Eq::{equal, not_equal}

///|
pub fn MatchrMotif::new(
  name~ : String,
  pwm~ : Array[Array[Double]],
) -> MatchrMotif {
  MatchrMotif::{ name, pwm, }
}

///|
/// Length of the motif (number of columns).
pub fn MatchrMotif::len(self : MatchrMotif) -> Int {
  if self.pwm.is_empty() {
    0
  } else {
    self.pwm[0].length()
  }
}

///|
/// Result of matching one motif against a set of sequences.
///
/// `hits[seq][pos]` is true when the PSSM score at `pos` reaches the
/// cutoff; `scores[seq][pos]` is that score (forward vs reverse-complement
/// maximum); `counts[seq]` is the number of hits in the sequence;
/// `max_scores[seq]` is the best score over all positions.
pub struct MotifMatches {
  motif : MatchrMotif
  seq_names : Array[String]
  hits : Array[Array[Bool]]
  scores : Array[Array[Double]]
  counts : Array[Int]
  max_scores : Array[Double]
} derive(Eq)
pub extend MotifMatches with Eq::{equal, not_equal}

///|
/// Default uniform background nucleotide frequencies.
pub fn motifmatchr_default_bg() -> Array[Double] {
  [0.25, 0.25, 0.25, 0.25]
}

///|
/// Index of a nucleotide character in the A, C, G, T ordering, or `-1`
/// for non-ACGT characters (including N).
fn motifmatchr_base_index(c : Int) -> Int {
  match c {
    65 => 0
    67 => 1
    71 => 2
    84 => 3
    _ => -1
  }
}

///|
/// Convert a position frequency matrix (counts) into a PWM.
///
/// `pseudocount` is distributed across bases according to `bg`. The
/// resulting columns each sum to 1.
pub fn motif_pwm_from_counts(
  counts : Array[Array[Int]],
  bg~ : Array[Double] = motifmatchr_default_bg(),
  pseudocount~ : Double = 0.8,
) -> Array[Array[Double]] {
  let n_bases = counts.length()
  let width = if n_bases == 0 { 0 } else { counts[0].length() }
  let pwm : Array[Array[Double]] = Array(capacity=0)
  for _i in 0..<n_bases {
    let row : Array[Double] = Array(capacity=0)
    for _j in 0..<width {
      row.push(0.0)
    }
    pwm.push(row)
  }
  for j in 0..<width {
    let mut col_total = 0.0
    for i in 0..<n_bases {
      col_total += counts[i][j].to_double()
    }
    let denom = col_total + pseudocount
    for i in 0..<n_bases {
      pwm[i][j] = (counts[i][j].to_double() + pseudocount * bg[i]) / denom
    }
  }
  pwm
}

///|
/// Build a PSSM (log-odds) matrix from a PWM and background.
///
/// `pssm[i][j] = log2(pwm[i][j] / bg[i])`. A tiny floor is applied so a
/// zero PWM entry yields a large negative score rather than -infinity.
pub fn motif_pssm_from_pwm(
  pwm : Array[Array[Double]],
  bg~ : Array[Double] = motifmatchr_default_bg(),
) -> Array[Array[Double]] {
  let n_bases = pwm.length()
  let width = if n_bases == 0 { 0 } else { pwm[0].length() }
  let pssm : Array[Array[Double]] = Array(capacity=0)
  for i in 0..<n_bases {
    let row : Array[Double] = Array(capacity=0)
    for j in 0..<width {
      let p = pwm[i][j]
      let ratio = p / bg[i]
      row.push(motifmatchr_log2(ratio))
    }
    pssm.push(row)
  }
  pssm
}

///|
/// log2(x) computed by bisection on 2^k, accurate enough for PSSM scores.
fn motifmatchr_log2(x : Double) -> Double {
  if x <= 0.0 {
    return -20.0
  }
  let mut lo = -20.0
  let mut hi = 20.0
  for _i in 0..<60 {
    let mid = (lo + hi) / 2.0
    let pow = motifmatchr_pow2(mid)
    if pow < x {
      lo = mid
    } else {
      hi = mid
    }
  }
  (lo + hi) / 2.0
}

///|
fn motifmatchr_pow2(x : Double) -> Double {
  // 2^x = exp(x * ln2); use exp via series after range reduction.
  let ln2 = 0.6931471805599453
  motifmatchr_exp(x * ln2)
}

///|
/// exp(x) via Taylor series with argument reduction for stability.
fn motifmatchr_exp(x : Double) -> Double {
  if x > 40.0 {
    return motifmatchr_exp(40.0) * motifmatchr_exp(x - 40.0)
  }
  if x < -40.0 {
    return 0.0
  }
  let mut sum = 1.0
  let mut term = 1.0
  for i in 1..<30 {
    term = term * x / i.to_double()
    sum += term
    let abs_term = if term < 0.0 { 0.0 - term } else { term }
    if abs_term < 1.0e-12 {
      break
    }
  }
  sum
}

///|
/// Reverse complement of a DNA string (UPPERCASE ACGT only; other chars
/// pass through unchanged).
pub fn motifmatchr_revcomp(seq : String) -> String {
  let mut result = ""
  for i = seq.length() - 1; i >= 0; i = i - 1 {
    let c = seq[i:i + 1].to_owned()
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
/// Score a window of `width` characters starting at `start` against the
/// PSSM. Returns a very negative number when the window contains a
/// non-ACGT character.
fn motifmatchr_score_window(
  pssm : Array[Array[Double]],
  seq : String,
  start : Int,
  width : Int,
) -> Double {
  let mut total = 0.0
  for j in 0..<width {
    let idx = motifmatchr_base_index(seq[start + j].to_int())
    if idx < 0 {
      return -1.0e9
    }
    total += pssm[idx][j]
  }
  total
}

///|
/// Match one motif against one sequence, populating per-position hit and
/// score arrays (forward and reverse-complement maximum).
fn motifmatchr_scan_one(
  pssm : Array[Array[Double]],
  seq : String,
  width : Int,
  cutoff : Double,
) -> (Array[Bool], Array[Double], Int, Double) {
  let n = seq.length()
  let hits : Array[Bool] = Array(capacity=0)
  let scores : Array[Double] = Array(capacity=0)
  let mut count = 0
  let mut max_score = -1.0e9
  if width == 0 || n < width {
    return (hits, scores, 0, -1.0e9)
  }
  let rc = motifmatchr_revcomp(seq)
  for pos in 0..<(n - width + 1) {
    let fwd = motifmatchr_score_window(pssm, seq, pos, width)
    let rev = motifmatchr_score_window(pssm, rc, n - pos - width, width)
    let best = if fwd > rev { fwd } else { rev }
    scores.push(best)
    let hit = best >= cutoff
    hits.push(hit)
    if hit {
      count += 1
    }
    if best > max_score {
      max_score = best
    }
  }
  (hits, scores, count, max_score)
}

///|
/// Match a collection of motifs against a collection of sequences.
///
/// Returns one `MotifMatches` per motif. `cutoff` is the minimum PSSM
/// score for a position to count as a hit; the reference package derives
/// it from `p.cutoff` but we expose it directly for transparency.
pub fn match_motifs(
  motifs : Array[MatchrMotif],
  sequences : Array[String],
  seq_names~ : Array[String] = [],
  bg~ : Array[Double] = motifmatchr_default_bg(),
  cutoff~ : Double = 0.0,
) -> Array[MotifMatches] {
  let names : Array[String] = Array(capacity=0)
  for i in 0..<sequences.length() {
    if i < seq_names.length() {
      names.push(seq_names[i])
    } else {
      names.push("seq" + i.to_string())
    }
  }
  let results : Array[MotifMatches] = Array(capacity=0)
  for motif in motifs {
    let pssm = motif_pssm_from_pwm(motif.pwm, bg~)
    let width = motif.len()
    let hits_all : Array[Array[Bool]] = Array(capacity=0)
    let scores_all : Array[Array[Double]] = Array(capacity=0)
    let counts : Array[Int] = Array(capacity=0)
    let max_scores : Array[Double] = Array(capacity=0)
    for seq in sequences {
      let (hits, scores, count, best) = motifmatchr_scan_one(
        pssm, seq, width, cutoff,
      )
      hits_all.push(hits)
      scores_all.push(scores)
      counts.push(count)
      max_scores.push(best)
    }
    results.push(
      MotifMatches::{
        motif,
        seq_names: names,
        hits: hits_all,
        scores: scores_all,
        counts,
        max_scores,
      },
    )
  }
  results
}

///|
/// Total hit count for a motif across all sequences.
pub fn MotifMatches::total_hits(self : MotifMatches) -> Int {
  let mut total = 0
  for c in self.counts {
    total += c
  }
  total
}

///|
/// Best score for a motif across all sequences.
pub fn MotifMatches::best_score(self : MotifMatches) -> Double {
  let mut best = -1.0e9
  for s in self.max_scores {
    if s > best {
      best = s
    }
  }
  best
}

///|
/// Hit matrix as a flat list of `(seq_index, position)` pairs.
pub fn MotifMatches::hit_positions(self : MotifMatches) -> Array[(Int, Int)] {
  let positions : Array[(Int, Int)] = Array(capacity=0)
  for seq_idx in 0..<self.hits.length() {
    for pos in 0..<self.hits[seq_idx].length() {
      if self.hits[seq_idx][pos] {
        positions.push((seq_idx, pos))
      }
    }
  }
  positions
}

///|
/// A small example motif (a 4-mer "ACGT" palindrome PWM) plus sequences,
/// used by the demo and tests.
pub fn motifmatchr_sample_motif() -> MatchrMotif {
  // PFM for ACGT: strong diagonal counts.
  let pfm = [
    [10, 0, 0, 0],
    [0, 10, 0, 0],
    [0, 0, 10, 0],
    [0, 0, 0, 10],
  ]
  let pwm = motif_pwm_from_counts(pfm)
  MatchrMotif::new(name="ACGT", pwm~)
}

///|
/// Sample sequences containing (and not containing) the ACGT motif.
pub fn motifmatchr_sample_sequences() -> Array[String] {
  ["GGGACGTGGG", "TTTTACGTTTT", "ATATATAT", "ACGTACGT"]
}