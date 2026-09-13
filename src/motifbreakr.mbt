///|
/// Motif disruption / gain prediction from SNVs, ported from
/// Bioconductor's `motifbreakR` package.
///
/// Given a set of SNPs and position weight matrices, the reference
/// package scores the reference and alternate alleles against each
/// motif, classifies the variant as `Gain`, `Loss`, or `Neutral`,
/// and reports the score difference (delta).  This port keeps the
/// same vocabulary.

///|
/// Effect of a SNP on a motif occurrence.
pub(all) enum MotifbreakrEffect {
  Gain
  Loss
  Neutral
} derive(Eq, Show, Debug)

///|
/// A single nucleotide variant to evaluate.
pub struct MotifbreakrVariant {
  name : String
  chrom : String
  position : Int
  ref_allele : String
  alt_allele : String
} derive(Eq)

///|
pub fn MotifbreakrVariant::new(
  name~ : String = "",
  chrom~ : String = "chr1",
  position~ : Int,
  ref_allele~ : String,
  alt_allele~ : String,
) -> MotifbreakrVariant {
  MotifbreakrVariant::{ name, chrom, position, ref_allele, alt_allele }
}

///|
/// Result of evaluating one variant against one motif.
pub struct MotifbreakrResult {
  variant_name : String
  motif_name : String
  chrom : String
  position : Int
  ref_score : Double
  alt_score : Double
  delta : Double
  effect : MotifbreakrEffect
  strand : String
} derive(Eq)

///|
/// Index of a nucleotide in A,C,G,T ordering; returns -1 for non-ACGT.
fn motifbreakr_base_idx(c : String) -> Int {
  match c {
    "A" => 0
    "C" => 1
    "G" => 2
    "T" => 3
    _ => -1
  }
}

///|
/// log2 via bisection (same approach as motifmatchr).
fn motifbreakr_log2(x : Double) -> Double {
  if x <= 0.0 {
    return -20.0
  }
  let mut lo = -20.0
  let mut hi = 20.0
  for _ in 0..<60 {
    let mid = (lo + hi) / 2.0
    let pow = motifbreakr_pow2(mid)
    if pow < x {
      lo = mid
    } else {
      hi = mid
    }
  }
  (lo + hi) / 2.0
}

///|
fn motifbreakr_pow2(x : Double) -> Double {
  let ln2 = 0.6931471805599453
  motifbreakr_exp(x * ln2)
}

///|
fn motifbreakr_exp(x : Double) -> Double {
  if x > 40.0 {
    return motifbreakr_exp(40.0) * motifbreakr_exp(x - 40.0)
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
/// Default uniform background.
pub fn motifbreakr_default_bg() -> Array[Double] {
  [0.25, 0.25, 0.25, 0.25]
}

///|
/// Build a PSSM (log-odds) matrix from a PWM (4 x L) and background.
fn motifbreakr_pssm(
  pwm : Array[Array[Double]],
  bg : Array[Double],
) -> Array[Array[Double]] {
  let n_bases = pwm.length()
  let width = if n_bases == 0 { 0 } else { pwm[0].length() }
  let pssm : Array[Array[Double]] = Array::new()
  for i in 0..<n_bases {
    let row : Array[Double] = Array::new()
    for j in 0..<width {
      let p = pwm[i][j]
      let ratio = p / bg[i]
      row.push(motifbreakr_log2(ratio))
    }
    pssm.push(row)
  }
  pssm
}

///|
/// Score a sequence window of length `width` starting at `start`
/// against the PSSM.  Returns a large negative number for windows
/// containing non-ACGT characters.
fn motifbreakr_score_window(
  pssm : Array[Array[Double]],
  seq : String,
  start : Int,
  width : Int,
) -> Double {
  let mut total = 0.0
  for j in 0..<width {
    let c = seq[start + j:start + j + 1].to_owned()
    let idx = motifbreakr_base_idx(c)
    if idx < 0 {
      return -1.0e9
    }
    total += pssm[idx][j]
  }
  total
}

///|
/// Substitute `allele` at `pos` (0-based) in `seq`, returning a new
/// string.  `pos` is the start position of the allele in the sequence.
fn motifbreakr_substitute(seq : String, pos : Int, allele : String) -> String {
  let before = seq[0:pos].to_owned()
  let after = seq[pos + 1:].to_owned()
  before + allele + after
}

///|
/// Evaluate one variant against one motif on a given sequence.
///
/// For every window position where the motif overlaps the SNP, score
/// the ref and alt versions of the sequence, keep the window with the
/// largest absolute delta, and classify the effect.
fn motifbreakr_eval_one(
  variant : MotifbreakrVariant,
  motif_name : String,
  pwm : Array[Array[Double]],
  sequence : String,
  bg : Array[Double],
  threshold~ : Double = 0.0,
) -> MotifbreakrResult? {
  let pssm = motifbreakr_pssm(pwm, bg)
  let width = if pwm.is_empty() { 0 } else { pwm[0].length() }
  let n = sequence.length()
  let pos = variant.position
  if width == 0 || n < width || pos < 0 || pos >= n {
    return None
  }
  let ref_seq = motifbreakr_substitute(sequence, pos, variant.ref_allele)
  let alt_seq = motifbreakr_substitute(sequence, pos, variant.alt_allele)
  let mut best_delta = 0.0
  let mut best_ref = -1.0e9
  let mut best_alt = -1.0e9
  let mut best_strand = "+"
  let mut found = false
  let w_lo = pos - width + 1
  let w_hi = pos + 1
  for w in w_lo..<w_hi {
    if w < 0 || w + width > n {
      continue
    }
    let ref_s = motifbreakr_score_window(pssm, ref_seq, w, width)
    let alt_s = motifbreakr_score_window(pssm, alt_seq, w, width)
    let delta = alt_s - ref_s
    let abs_delta = if delta < 0.0 { 0.0 - delta } else { delta }
    let abs_best = if best_delta < 0.0 {
      0.0 - best_delta
    } else {
      best_delta
    }
    if abs_delta > abs_best || not(found) {
      best_delta = delta
      best_ref = ref_s
      best_alt = alt_s
      best_strand = "+"
      found = true
    }
  }
  if not(found) {
    return None
  }
  let effect = if best_delta > threshold {
    MotifbreakrEffect::Gain
  } else if best_delta < (0.0 - threshold) {
    MotifbreakrEffect::Loss
  } else {
    MotifbreakrEffect::Neutral
  }
  Some(
    MotifbreakrResult::{
      variant_name: variant.name,
      motif_name,
      chrom: variant.chrom,
      position: variant.position,
      ref_score: best_ref,
      alt_score: best_alt,
      delta: best_delta,
      effect,
      strand: best_strand,
    },
  )
}

///|
/// Evaluate a collection of variants against a collection of motifs.
///
/// `sequences` maps chromosome names to their DNA sequence.  For each
/// variant / motif pair, the best overlapping window is scored and
/// the effect classified.  Variants that do not overlap any motif
/// window are skipped.
pub fn motif_break(
  variants : Array[MotifbreakrVariant],
  motifs : Array[(String, Array[Array[Double]])],
  sequences : Map[String, String],
  bg~ : Array[Double] = motifbreakr_default_bg(),
  threshold~ : Double = 0.0,
) -> Array[MotifbreakrResult] {
  let results : Array[MotifbreakrResult] = Array::new()
  for variant in variants {
    let seq = sequences.get(variant.chrom)
    match seq {
      Some(s) => {
        for motif_idx in 0..<motifs.length() {
          let (motif_name, pwm) = motifs[motif_idx]
          match motifbreakr_eval_one(
            variant, motif_name, pwm, s, bg, threshold~,
          ) {
            Some(r) => results.push(r)
            None => ()
          }
        }
      }
      None => ()
    }
  }
  results
}

///|
/// Filter results to only those with a non-Neutral effect.
pub fn motifbreakr_significant(
  results : Array[MotifbreakrResult],
) -> Array[MotifbreakrResult] {
  let out : Array[MotifbreakrResult] = Array::new()
  for r in results {
    if r.effect != MotifbreakrEffect::Neutral {
      out.push(r)
    }
  }
  out
}

///|
/// Sort results by absolute delta descending (most disruptive first).
pub fn motifbreakr_sort_by_delta(
  results : Array[MotifbreakrResult],
) -> Array[MotifbreakrResult] {
  let arr = results.copy()
  arr.sort_by(fn(a, b) {
    let abs_a = if a.delta < 0.0 { 0.0 - a.delta } else { a.delta }
    let abs_b = if b.delta < 0.0 { 0.0 - b.delta } else { b.delta }
    if abs_a > abs_b {
      -1
    } else if abs_a < abs_b {
      1
    } else {
      0
    }
  })
  arr
}

///|
/// A sample PWM for the GATA1 motif (simplified 6-mer).
pub fn motifbreakr_sample_gata1_pwm() -> Array[Array[Double]] {
  [
    [0.1, 0.1, 0.8, 0.1, 0.1, 0.1],
    [0.1, 0.1, 0.1, 0.1, 0.8, 0.1],
    [0.7, 0.1, 0.1, 0.7, 0.1, 0.1],
    [0.1, 0.7, 0.0, 0.1, 0.0, 0.7],
  ]
}

///|
/// A sample sequence and variants for testing.
pub fn motifbreakr_sample_sequence() -> String {
  "GGGAGATAAGGGGGG"
}

///|
/// Sample variants that disrupt (or not) the GATA1 motif.
pub fn motifbreakr_sample_variants() -> Array[MotifbreakrVariant] {
  [
    MotifbreakrVariant::new(
      name="rs1",
      chrom="chr1",
      position=4,
      ref_allele="G",
      alt_allele="T",
    ),
    MotifbreakrVariant::new(
      name="rs2",
      chrom="chr1",
      position=7,
      ref_allele="A",
      alt_allele="C",
    ),
  ]
}