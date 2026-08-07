///|
/// Codon alignment advanced selection-pressure statistics, inspired by
/// Biopython's `Bio.codonalign.test_selection` / `test_neutrality` and the
/// Nei–Gojobori (1986) large-sample variance framework.
///
/// Builds on the base `codon_align.mbt` NG86 dN/dS engine and adds:
///   - Z-test for selection (positive / purifying / neutrality) using the
///     NG86 approximate variance V(dN-dS) = V(dN)+V(dS).
///   - Fisher's exact test for neutrality on the 2×2 contingency table
///     {synonymous, nonsynonymous} × {differences, identities}, evaluated via
///     the hypergeometric distribution with log-factorial numerics.
///   - Codon alignment builder from a protein alignment + coding sequences.
///   - Sliding-window dN/dS scan for selection hotspots.
///   - Pairwise Ka/Ks table with BH-FDR correction.
///   - Benjamini–Hochberg FDR helper.

// ===========================================================================
// Error type
// ===========================================================================

///|
pub suberror CodonAlignAdvancedError {
  CodonAlignAdvancedError(String)
}

///|
fn caa_fail(message : String) -> Unit raise CodonAlignAdvancedError {
  raise CodonAlignAdvancedError::CodonAlignAdvancedError(message)
}

// ===========================================================================
// NG86 intermediate components (private)
// ===========================================================================

///|
priv struct CaaNg86Components {
  dn : Double // Jukes-Cantor corrected dN
  ds : Double // Jukes-Cantor corrected dS
  p_n : Double // pN = Nd / Nn (raw nonsyn difference proportion)
  p_s : Double // pS = Sd / Sn (raw syn difference proportion)
  n_diff : Int // Nd — nonsynonymous differences
  s_diff : Int // Sd — synonymous differences
  n_sites : Double // Nn — nonsynonymous sites (averaged over both seqs)
  s_sites : Double // Sn — synonymous sites (averaged over both seqs)
  n_codons : Int // number of valid codon pairs analysed
}

///|
/// Compute Nei–Gojobori (1986) intermediate values for a codon-aligned pair.
/// Reuses `codon_to_aa`, `is_stop_codon`, `count_codon_sites`,
/// `classify_codon_substitution`, and `jukes_cantor_correction` from the base
/// module so that site/ substitution definitions stay identical.
fn caa_ng86_components(
  seq1 : String,
  seq2 : String,
) -> CaaNg86Components raise CodonAlignAdvancedError {
  let len1 = seq1.length()
  let len2 = seq2.length()
  if len1 < 3 || len2 < 3 {
    caa_fail("CodonAlign advanced requires sequences of at least 3 nt")
  }
  if len1 != len2 {
    caa_fail("CodonAlign advanced requires equal-length sequences")
  }
  if len1 % 3 != 0 {
    caa_fail("CodonAlign advanced requires codon-aligned (multiple-of-3) sequences")
  }
  let n_codons_total = len1 / 3
  let mut s_sites = 0.0
  let mut ns_sites = 0.0
  let mut s_diff = 0
  let mut ns_diff = 0
  let mut valid_codons = 0
  let mut i = 0
  while i < n_codons_total {
    let start = i * 3
    let codon1 = seq1[start:start + 3].to_owned()
    let codon2 = seq2[start:start + 3].to_owned()
    // Skip codon pairs containing gaps or stop codons.
    if !is_stop_codon(codon1) && !is_stop_codon(codon2) && codon1.length() == 3 && codon2.length() == 3 {
      let (s1, ns1) = count_codon_sites(codon1)
      let (s2, ns2) = count_codon_sites(codon2)
      s_sites = s_sites + (s1 + s2) / 2.0
      ns_sites = ns_sites + (ns1 + ns2) / 2.0
      match classify_codon_substitution(codon1, codon2) {
        Synonymous => if codon1 != codon2 {
          s_diff = s_diff + 1
        }
        NonSynonymous => ns_diff = ns_diff + 1
        Unknown => ()
      }
      valid_codons = valid_codons + 1
    }
    i = i + 1
  }
  if valid_codons == 0 {
    caa_fail("CodonAlign advanced found no analysable codon pairs")
  }
  let p_s = if s_sites > 0.0 { s_diff.to_double() / s_sites } else { 0.0 }
  let p_n = if ns_sites > 0.0 {
    ns_diff.to_double() / ns_sites
  } else {
    0.0
  }
  let d_s = jukes_cantor_correction(p_s)
  let d_n = jukes_cantor_correction(p_n)
  CaaNg86Components::{
    dn: d_n,
    ds: d_s,
    p_n,
    p_s,
    n_diff: ns_diff,
    s_diff,
    n_sites: ns_sites,
    s_sites,
    n_codons: valid_codons,
  }
}

// ===========================================================================
// Normal CDF (Abramowitz–Stegun 7.1.26)
// ===========================================================================

///|
fn caa_normal_cdf(x : Double) -> Double {
  let abs_x = if x < 0.0 { -x } else { x }
  let t = 1.0 / (1.0 + 0.2316419 * abs_x)
  let d = 0.3989422804014327 * @math.exp(-x * x / 2.0)
  let p = d *
    t *
    (0.31938153 +
    t *
      (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429))))
  if x > 0.0 {
    1.0 - p
  } else {
    p
  }
}

// ===========================================================================
// Z-test for selection (Nei–Gojobori approximate variance)
// ===========================================================================

///|
/// Result of the Z-test for selection.
pub struct SelectionTestResult {
  dn : Double // corrected nonsynonymous rate
  ds : Double // corrected synonymous rate
  dnds : Double // dn/ds
  variance : Double // V(dN - dS)
  z_score : Double
  p_value : Double
  test_type : String // "positive" | "purifying" | "neutrality"
  conclusion : String // human-readable conclusion
} derive(Debug)

///|
pub fn SelectionTestResult::dn(self : SelectionTestResult) -> Double {
  self.dn
}

///|
pub fn SelectionTestResult::ds(self : SelectionTestResult) -> Double {
  self.ds
}

///|
pub fn SelectionTestResult::dnds(self : SelectionTestResult) -> Double {
  self.dnds
}

///|
pub fn SelectionTestResult::variance(self : SelectionTestResult) -> Double {
  self.variance
}

///|
pub fn SelectionTestResult::z_score(self : SelectionTestResult) -> Double {
  self.z_score
}

///|
pub fn SelectionTestResult::p_value(self : SelectionTestResult) -> Double {
  self.p_value
}

///|
pub fn SelectionTestResult::test_type(self : SelectionTestResult) -> String {
  self.test_type
}

///|
pub fn SelectionTestResult::conclusion(self : SelectionTestResult) -> String {
  self.conclusion
}

///|
/// Z-test for selection using the Nei–Gojobori (1986) large-sample variance.
///
/// `test_type` selects the alternative hypothesis:
///   - "positive"   — H₁: dN > dS  (one-tailed, p = 1 − Φ(Z))
///   - "purifying"  — H₁: dN < dS  (one-tailed, p = Φ(Z))
///   - "neutrality" — H₁: dN ≠ dS  (two-tailed, p = 2(1 − Φ(|Z|)))
///
/// Variance (approximate, covariance ignored):
///   V(dN) = pN(1−pN) / [Nn · (1 − 4pN/3)²]
///   V(dS) = pS(1−pS) / [Sn · (1 − 4pS/3)²]
///   V(dN−dS) = V(dN) + V(dS)
///   Z = (dN − dS) / √V(dN−dS)
pub fn codon_test_selection(
  seq1 : String,
  seq2 : String,
  test_type? : String = "positive",
) -> SelectionTestResult raise CodonAlignAdvancedError {
  let comps = caa_ng86_components(seq1, seq2)
  let p_n = comps.p_n
  let p_s = comps.p_s
  let nn = comps.n_sites
  let sn = comps.s_sites
  // Variance under Jukes-Cantor correction (Nei & Gojobori 1986, eq. 6).
  let v_dn = if nn > 0.0 && (1.0 - 4.0 * p_n / 3.0) > 0.0 {
    p_n * (1.0 - p_n) / (nn * (1.0 - 4.0 * p_n / 3.0) * (1.0 - 4.0 * p_n / 3.0))
  } else {
    0.0
  }
  let v_ds = if sn > 0.0 && (1.0 - 4.0 * p_s / 3.0) > 0.0 {
    p_s * (1.0 - p_s) / (sn * (1.0 - 4.0 * p_s / 3.0) * (1.0 - 4.0 * p_s / 3.0))
  } else {
    0.0
  }
  let variance = v_dn + v_ds
  let dn = comps.dn
  let ds = comps.ds
  let dnds = if ds == 0.0 {
    if dn > 0.0 {
      999.0
    } else {
      0.0
    }
  } else {
    dn / ds
  }
  // Z-score.
  let (z_score, p_value) = if variance <= 0.0 {
    (0.0, 1.0)
  } else {
    let z = (dn - ds) / variance.sqrt()
    let p = match test_type {
      "positive" => 1.0 - caa_normal_cdf(z) // H1: dN > dS
      "purifying" => caa_normal_cdf(z) // H1: dN < dS
      _ => {
        let abs_z = if z < 0.0 { -z } else { z }
        2.0 * (1.0 - caa_normal_cdf(abs_z))
      }
    }
    (z, p)
  }
  // Conclusion at α = 0.05.
  let conclusion = if variance <= 0.0 {
    "inconclusive (zero variance)"
  } else if p_value < 0.05 {
    match test_type {
      "positive" => "significant positive selection (dN > dS)"
      "purifying" => "significant purifying selection (dN < dS)"
      _ => if z_score > 0.0 {
        "significant positive selection (dN > dS)"
      } else {
        "significant purifying selection (dN < dS)"
      }
    }
  } else {
    "not significant (consistent with neutrality)"
  }
  SelectionTestResult::{
    dn,
    ds,
    dnds,
    variance,
    z_score,
    p_value,
    test_type,
    conclusion,
  }
}

// ===========================================================================
// Fisher's exact test for neutrality
// ===========================================================================

///|
/// Log-factorial computed by summation of natural logs (no overflow).
fn caa_log_factorial(n : Int) -> Double {
  if n <= 1 {
    return 0.0
  }
  let mut result = 0.0
  let mut i = 2
  while i <= n {
    result = result + @math.ln(i.to_double())
    i = i + 1
  }
  result
}

///|
/// Log binomial coefficient ln C(n, k) = ln(n!) − ln(k!) − ln((n−k)!).
fn caa_log_binom(n : Int, k : Int) -> Double {
  if k < 0 || k > n {
    return -1.0e300
  }
  caa_log_factorial(n) - caa_log_factorial(k) - caa_log_factorial(n - k)
}

///|
/// Fisher's exact test on a 2×2 contingency table.
///
/// Table layout:
/// ```
///   | a | b |
///   | c | d |
/// ```
/// Returns `(p_value, odds_ratio)`.
///
/// `test_type`:
///   - "two-sided"  — sum of probabilities ≤ observed
///   - "greater"    — P(X ≥ a)
///   - "less"       — P(X ≤ a)
fn caa_fisher_exact(
  a : Int,
  b : Int,
  c : Int,
  d : Int,
  test_type : String,
) -> (Double, Double) {
  let n = a + b + c + d
  let row1 = a + b // Nn (nonsynonymous sites)
  let row2 = c + d // Sn (synonymous sites)
  let col1 = a + c // D  (total differences)
  // Haldane–Anscombe correction for odds ratio when a cell is zero.
  let odds_ratio = if b == 0 || c == 0 || d == 0 || a == 0 {
    let a2 = (a.to_double() + 0.5)
    let b2 = (b.to_double() + 0.5)
    let c2 = (c.to_double() + 0.5)
    let d2 = (d.to_double() + 0.5)
    a2 * d2 / (b2 * c2)
  } else {
    a.to_double() * d.to_double() / (b.to_double() * c.to_double())
  }
  // Range of possible values for the (1,1) cell.
  let k_min = if col1 - row2 > 0 { col1 - row2 } else { 0 }
  let k_max = if col1 < row1 { col1 } else { row1 }
  // Log probability of observed table.
  let log_denom = caa_log_binom(n, col1)
  let log_p_observed = caa_log_binom(row1, a) +
    caa_log_binom(row2, c) -
    log_denom
  let p_observed = @math.exp(log_p_observed)
  let mut p_value = 0.0
  let mut k = k_min
  while k <= k_max {
    let log_p_k = caa_log_binom(row1, k) +
      caa_log_binom(row2, col1 - k) -
      log_denom
    let p_k = @math.exp(log_p_k)
    match test_type {
      "greater" => if k >= a { p_value = p_value + p_k }
      "less" => if k <= a { p_value = p_value + p_k }
      _ => if p_k <= p_observed + 1.0e-15 { p_value = p_value + p_k }
    }
    k = k + 1
  }
  // Clamp to [0, 1].
  if p_value > 1.0 {
    p_value = 1.0
  }
  if p_value < 0.0 {
    p_value = 0.0
  }
  (p_value, odds_ratio)
}

///|
/// Result of Fisher's exact test for neutrality.
pub struct CodonNeutralityTestResult {
  p_value : Double
  odds_ratio : Double
  n_diff : Int // nonsynonymous differences (Nd)
  s_diff : Int // synonymous differences (Sd)
  n_sites : Double // nonsynonymous sites (Nn)
  s_sites : Double // synonymous sites (Sn)
  test_type : String
  conclusion : String
} derive(Debug)

///|
pub fn CodonNeutralityTestResult::p_value(self : CodonNeutralityTestResult) -> Double {
  self.p_value
}

///|
pub fn CodonNeutralityTestResult::odds_ratio(self : CodonNeutralityTestResult) -> Double {
  self.odds_ratio
}

///|
pub fn CodonNeutralityTestResult::n_diff(self : CodonNeutralityTestResult) -> Int {
  self.n_diff
}

///|
pub fn CodonNeutralityTestResult::s_diff(self : CodonNeutralityTestResult) -> Int {
  self.s_diff
}

///|
pub fn CodonNeutralityTestResult::n_sites(self : CodonNeutralityTestResult) -> Double {
  self.n_sites
}

///|
pub fn CodonNeutralityTestResult::s_sites(self : CodonNeutralityTestResult) -> Double {
  self.s_sites
}

///|
pub fn CodonNeutralityTestResult::test_type(self : CodonNeutralityTestResult) -> String {
  self.test_type
}

///|
pub fn CodonNeutralityTestResult::conclusion(
  self : CodonNeutralityTestResult,
) -> String {
  self.conclusion
}

///|
/// Fisher's exact test for neutrality on the 2×2 table:
/// ```
///                | differences | identities |
/// nonsynonymous  | Nd          | Nn − Nd    |
/// synonymous     | Sd          | Sn − Sd    |
/// ```
/// `test_type`: "two-sided" (default), "greater" (positive selection),
/// "less" (purifying selection).
pub fn codon_test_neutrality(
  seq1 : String,
  seq2 : String,
  test_type? : String = "two-sided",
) -> CodonNeutralityTestResult raise CodonAlignAdvancedError {
  let comps = caa_ng86_components(seq1, seq2)
  let nd = comps.n_diff
  let sd = comps.s_diff
  let nn_d = comps.n_sites
  let sn_d = comps.s_sites
  let nn_int = nn_d.round().to_int()
  let sn_int = sn_d.round().to_int()
  // 2×2 table:
  //   a = Nd           b = Nn − Nd
  //   c = Sd           d = Sn − Sd
  let a = nd
  let b = if nn_int - nd > 0 { nn_int - nd } else { 0 }
  let c = sd
  let d = if sn_int - sd > 0 { sn_int - sd } else { 0 }
  let (p_value, odds_ratio) = caa_fisher_exact(a, b, c, d, test_type)
  let conclusion = if p_value < 0.05 {
    match test_type {
      "greater" => "significant positive selection (reject neutrality)"
      "less" => "significant purifying selection (reject neutrality)"
      _ => if odds_ratio > 1.0 {
        "significant positive selection (reject neutrality)"
      } else {
        "significant purifying selection (reject neutrality)"
      }
    }
  } else {
    "not significant (consistent with neutrality)"
  }
  CodonNeutralityTestResult::{
    p_value,
    odds_ratio,
    n_diff: nd,
    s_diff: sd,
    n_sites: nn_d,
    s_sites: sn_d,
    test_type,
    conclusion,
  }
}

// ===========================================================================
// Codon alignment builder
// ===========================================================================

///|
/// A codon alignment: aligned nucleotide sequences built from a protein
/// alignment and un-aligned coding sequences.
pub struct CodonAlignment {
  sequences : Array[String]
  names : Array[String]
  n_codons : Int // number of codon columns (including gap columns)
} derive(Debug)

///|
pub fn CodonAlignment::sequences(self : CodonAlignment) -> Array[String] {
  self.sequences
}

///|
pub fn CodonAlignment::names(self : CodonAlignment) -> Array[String] {
  self.names
}

///|
pub fn CodonAlignment::n_codons(self : CodonAlignment) -> Int {
  self.n_codons
}

///|
/// Build a codon alignment from a protein alignment and coding sequences.
///
/// For each amino-acid position in the protein alignment:
///   - If the residue is `-` (gap), insert `---` in the codon alignment.
///   - Otherwise consume the next codon from the coding sequence and verify
///     that it translates to the expected amino acid.
///
/// `protein_aln` and `nucleotide_seqs` must have the same length (number of
/// sequences).  Each protein string must have identical length across all
/// entries.
pub fn build_codon_alignment(
  protein_aln : Array[String],
  nucleotide_seqs : Array[String],
  names? : Array[String] = [],
) -> CodonAlignment raise CodonAlignAdvancedError {
  let n_seqs = protein_aln.length()
  if n_seqs == 0 {
    caa_fail("CodonAlign advanced builder requires at least one sequence")
  }
  if nucleotide_seqs.length() != n_seqs {
    caa_fail("CodonAlign advanced builder: protein / nucleotide count mismatch")
  }
  let aln_len = protein_aln[0].length()
  for i in 0..<n_seqs {
    if protein_aln[i].length() != aln_len {
      caa_fail("CodonAlign advanced builder: protein alignment rows have different lengths")
    }
  }
  let seq_names : Array[String] = []
  if names.length() == n_seqs {
    for n in names {
      seq_names.push(n)
    }
  } else {
    for i in 0..<n_seqs {
      seq_names.push("seq" + i.to_string())
    }
  }
  let aligned : Array[String] = []
  for i in 0..<n_seqs {
    let protein = protein_aln[i]
    let coding = nucleotide_seqs[i]
    let mut codon_seq = ""
    let mut nt_idx = 0
    let coding_len = coding.length()
    let mut pos = 0
    while pos < aln_len {
      let aa = protein[pos:pos + 1].to_owned()
      if aa == "-" || aa == "." {
        codon_seq = codon_seq + "---"
      } else {
        if nt_idx + 3 > coding_len {
          caa_fail(
            "CodonAlign advanced builder: coding sequence " +
            i.to_string() +
            " exhausted before protein alignment end",
          )
        }
        let codon = coding[nt_idx:nt_idx + 3].to_owned().to_upper()
        let translated = codon_to_aa(codon)
        let aa_upper = aa.to_upper()
        if translated != aa_upper && translated != "?" {
          caa_fail(
            "CodonAlign advanced builder: codon " +
            codon +
            " translates to " +
            translated +
            " but protein alignment expects " +
            aa_upper +
            " at position " +
            pos.to_string() +
            " of sequence " +
            i.to_string(),
          )
        }
        codon_seq = codon_seq + codon
        nt_idx = nt_idx + 3
      }
      pos = pos + 1
    }
    if nt_idx != coding_len {
      caa_fail(
        "CodonAlign advanced builder: coding sequence " +
        i.to_string() +
        " has " +
        (coding_len - nt_idx).to_string() +
        " leftover nucleotide(s) after protein alignment",
      )
    }
    aligned.push(codon_seq)
  }
  CodonAlignment::{ sequences: aligned, names: seq_names, n_codons: aln_len }
}

// ===========================================================================
// Sliding-window dN/dS
// ===========================================================================

///|
/// One window of a sliding-window dN/dS scan.
pub struct WindowDnDs {
  start_codon : Int // 0-based start codon index
  end_codon : Int // 0-based end codon index (exclusive)
  dn : Double
  ds : Double
  dnds : Double
  n_codons : Int // number of valid codon pairs in the window
} derive(Debug)

///|
pub fn WindowDnDs::start_codon(self : WindowDnDs) -> Int {
  self.start_codon
}

///|
pub fn WindowDnDs::end_codon(self : WindowDnDs) -> Int {
  self.end_codon
}

///|
pub fn WindowDnDs::dn(self : WindowDnDs) -> Double {
  self.dn
}

///|
pub fn WindowDnDs::ds(self : WindowDnDs) -> Double {
  self.ds
}

///|
pub fn WindowDnDs::dnds(self : WindowDnDs) -> Double {
  self.dnds
}

///|
pub fn WindowDnDs::n_codons(self : WindowDnDs) -> Int {
  self.n_codons
}

///|
/// Result of a sliding-window dN/dS scan.
pub struct SlidingWindowResult {
  windows : Array[WindowDnDs]
} derive(Debug)

///|
pub fn SlidingWindowResult::windows(self : SlidingWindowResult) -> Array[WindowDnDs] {
  self.windows
}

///|
/// Sliding-window dN/dS scan for detecting selection hotspots.
///
/// Walks along the codon alignment in steps of `step_size` codons, computing
/// NG86 dN/dS for each window of `window_size` codons.  Windows with fewer
/// than 2 valid codon pairs are skipped.
pub fn sliding_window_dnds(
  seq1 : String,
  seq2 : String,
  window_size? : Int = 10,
  step_size? : Int = 5,
) -> SlidingWindowResult raise CodonAlignAdvancedError {
  if window_size < 3 {
    caa_fail("CodonAlign advanced sliding window requires window_size >= 3")
  }
  if step_size < 1 {
    caa_fail("CodonAlign advanced sliding window requires step_size >= 1")
  }
  let len1 = seq1.length()
  let len2 = seq2.length()
  if len1 != len2 || len1 % 3 != 0 {
    caa_fail("CodonAlign advanced sliding window requires equal codon-aligned sequences")
  }
  let n_codons_total = len1 / 3
  let windows : Array[WindowDnDs] = []
  let mut start = 0
  while start < n_codons_total {
    let end = if start + window_size < n_codons_total {
      start + window_size
    } else {
      n_codons_total
    }
    if end - start >= 3 {
      let sub1 = seq1[start * 3:end * 3].to_owned()
      let sub2 = seq2[start * 3:end * 3].to_owned()
      // Use NG86 components directly (does not raise for valid sub-sequences).
      let comps = caa_ng86_components(sub1, sub2) catch {
        CodonAlignAdvancedError(_) => {
          start = start + step_size
          continue
        }
      }
      let dnds = if comps.ds == 0.0 {
        if comps.dn > 0.0 {
          999.0
        } else {
          0.0
        }
      } else {
        comps.dn / comps.ds
      }
      windows.push(
        WindowDnDs::{
          start_codon: start,
          end_codon: end,
          dn: comps.dn,
          ds: comps.ds,
          dnds,
          n_codons: comps.n_codons,
        },
      )
    }
    start = start + step_size
  }
  SlidingWindowResult::{ windows }
}

// ===========================================================================
// Benjamini–Hochberg FDR
// ===========================================================================

///|
/// Benjamini–Hochberg FDR correction.
pub fn codon_align_advanced_bh_fdr(p_values : Array[Double]) -> Array[Double] {
  let count = p_values.length()
  let output : Array[Double] = Array::make(count, 1.0)
  if count == 0 {
    return output
  }
  let indexed : Array[(Double, Int)] = []
  for index in 0..<count {
    let pv = p_values[index]
    let clamped = if pv < 0.0 {
      0.0
    } else if pv > 1.0 {
      1.0
    } else {
      pv
    }
    indexed.push((clamped, index))
  }
  indexed.sort_by(fn(left : (Double, Int), right : (Double, Int)) -> Int {
    if left.0 < right.0 {
      -1
    } else if left.0 > right.0 {
      1
    } else {
      left.1 - right.1
    }
  })
  let mut running = 1.0
  let mut position = count - 1
  while position >= 0 {
    let candidate = (indexed[position].0 *
    count.to_double() /
    (position + 1).to_double()).min(1.0)
    running = running.min(candidate)
    output[indexed[position].1] = running
    position = position - 1
  }
  output
}

// ===========================================================================
// Pairwise Ka/Ks table
// ===========================================================================

///|
/// One row of the pairwise Ka/Ks table with Z-test and FDR.
pub struct KaksTableRow {
  seq1_name : String
  seq2_name : String
  dn : Double
  ds : Double
  dnds : Double
  z_score : Double
  p_value : Double
  fdr : Double
  conclusion : String
} derive(Debug)

///|
pub fn KaksTableRow::seq1_name(self : KaksTableRow) -> String {
  self.seq1_name
}

///|
pub fn KaksTableRow::seq2_name(self : KaksTableRow) -> String {
  self.seq2_name
}

///|
pub fn KaksTableRow::dn(self : KaksTableRow) -> Double {
  self.dn
}

///|
pub fn KaksTableRow::ds(self : KaksTableRow) -> Double {
  self.ds
}

///|
pub fn KaksTableRow::dnds(self : KaksTableRow) -> Double {
  self.dnds
}

///|
pub fn KaksTableRow::z_score(self : KaksTableRow) -> Double {
  self.z_score
}

///|
pub fn KaksTableRow::p_value(self : KaksTableRow) -> Double {
  self.p_value
}

///|
pub fn KaksTableRow::fdr(self : KaksTableRow) -> Double {
  self.fdr
}

///|
pub fn KaksTableRow::conclusion(self : KaksTableRow) -> String {
  self.conclusion
}

///|
/// Compute a pairwise Ka/Ks table for multiple sequences with Z-test for
/// positive selection and BH-FDR correction.
///
/// Each pair of sequences is tested with the one-tailed Z-test
/// (H₁: dN > dS).  The resulting p-values are BH-adjusted across all pairs.
/// `fdr_threshold` (default 0.05) determines the "significant" conclusion.
pub fn pairwise_kaks_table(
  sequences : Array[String],
  names? : Array[String] = [],
  fdr_threshold? : Double = 0.05,
) -> Array[KaksTableRow] raise CodonAlignAdvancedError {
  let n = sequences.length()
  if n < 2 {
    caa_fail("CodonAlign advanced pairwise table requires at least 2 sequences")
  }
  let seq_names : Array[String] = []
  if names.length() == n {
    for nm in names {
      seq_names.push(nm)
    }
  } else {
    for i in 0..<n {
      seq_names.push("seq" + i.to_string())
    }
  }
  // Collect raw results and p-values.
  let raw : Array[(Int, Int, SelectionTestResult)] = []
  let p_values : Array[Double] = []
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let result = codon_test_selection(
        sequences[i],
        sequences[j],
        test_type="positive",
      )
      raw.push((i, j, result))
      p_values.push(result.p_value)
      j = j + 1
    }
    i = i + 1
  }
  let fdr = codon_align_advanced_bh_fdr(p_values)
  let rows : Array[KaksTableRow] = []
  for idx in 0..<raw.length() {
    let (si, sj, result) = raw[idx]
    let f = fdr[idx]
    let conclusion = if f < fdr_threshold {
      "significant positive selection (FDR < " + fdr_threshold.to_string() + ")"
    } else {
      "not significant"
    }
    rows.push(
      KaksTableRow::{
        seq1_name: seq_names[si],
        seq2_name: seq_names[sj],
        dn: result.dn,
        ds: result.ds,
        dnds: result.dnds,
        z_score: result.z_score,
        p_value: result.p_value,
        fdr: f,
        conclusion,
      },
    )
  }
  rows
}

// ===========================================================================
// Demo data
// ===========================================================================

///|
/// Demo: a codon alignment under purifying selection (mostly synonymous changes).
pub fn create_purifying_selection_alignment() -> Array[String] {
  // 9-codon (27 nt) sequences: M-A-L-K-W-Q-Q-P-W
  // Only synonymous codon switches across the variants.
  [
    "ATGGCCCTGAAATGGCAGCAGCCATGG",
    "ATGGCGCTCAAATGGCAGCAGCCATGG",
    "ATGGCCCTGAAATGGCAGCAGCGATGG",
    "ATGGCCCTCAAATGGCAGCAGCCATGG",
  ]
}

///|
/// Demo: a codon alignment under positive selection (many nonsynonymous changes).
pub fn create_positive_selection_alignment() -> Array[String] {
  // 9-codon (27 nt) sequences.
  // Seq1: M A L K W Q Q P W
  // Seq2: M G L E R H P S M  (7 nonsyn, 2 syn)
  [
    "ATGGCCCTGAAATGGCAGCAGCCATGG",
    "ATGGGTCTGGAGAGACACCCCAGCATG",
    "ATGGCCCTGAAATGGCAGCAGCCATGG",
    "ATGGGTCTGGAGAGACACCCCAGCATG",
    "ATGGCCCTGAAATGGCAGCAGCCATGG",
    "ATGGGTCTGGAGAGACACCCCAGCATG",
    "ATGGCCCTGAAATGGCAGCAGCCATGG",
    "ATGGGTCTGGAGAGACACCCCAGCATG",
    "ATGGCCCTGAAATGGCAGCAGCCATGG",
    "ATGGGTCTGGAGAGACACCCCAGCATG",
  ]
}

///|
/// Demo: a small protein alignment and coding sequences for the builder.
pub fn create_demo_protein_alignment() -> (Array[String], Array[String]) {
  // 5-column protein alignment with gaps.
  let protein_aln = ["M-A-K", "MTASK", "M-ASK"]
  let coding_seqs = [
    "ATGGCCAAG", // M-A-K → 3 codons (gaps at positions 2, 4)
    "ATGACCGCCAGCAAG", // MTASK → 5 codons
    "ATGGCCAGCAAG", // M-ASK → 4 codons (gap at position 2)
  ]
  (protein_aln, coding_seqs)
}
