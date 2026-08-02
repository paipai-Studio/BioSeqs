///|
/// MACS2-like ChIP-seq peak calling module.
///
/// Inspired by the MACS2 (Model-based Analysis of ChIP-Seq) algorithm
/// (Zhang et al. 2008, Genome Biol.), this module implements the core
/// peak calling pipeline: sliding window scanning, local background
/// (lambda) estimation from a control sample, Poisson significance
/// testing, candidate peak merging, Benjamini-Hochberg FDR correction,
/// and fold-enrichment based filtering.
///
/// All public identifiers are prefixed with `peak_calling_` or live
/// under the `ChipSeqRead`, `CandidatePeak`, and `PeakCallingParams`
/// types. Internal helpers use the `pc_` prefix to avoid collisions.

// ============================================================================
// Data structures
// ============================================================================

///|
/// A single ChIP-seq read (5' end position + strand).
pub struct ChipSeqRead {
  chr : String
  position : Int
  strand : String
}

///|
/// A candidate peak called by the algorithm. The `fdr` and
/// `significant` fields are mutable so that FDR correction and
/// filtering can update them in place.
pub struct CandidatePeak {
  chr : String
  start : Int
  end : Int
  summit : Int
  read_count : Int
  control_count : Int
  fold_enrichment : Double
  p_value : Double
  mut fdr : Double
  mut significant : Bool
}

///|
/// Parameters controlling peak calling.
pub struct PeakCallingParams {
  window_size : Int
  step_size : Int
  local_lambda_size : Int
  pvalue_threshold : Double
  fdr_threshold : Double
  min_fold_enrichment : Double
}

// ============================================================================
// ChipSeqRead constructors and accessors
// ============================================================================

///|
/// Construct a `ChipSeqRead` from a chromosome name, genomic
/// position, and strand ("+" or "-").
pub fn ChipSeqRead::new(chr~ : String, position~ : Int, strand~ : String) -> ChipSeqRead {
  ChipSeqRead::{ chr, position, strand }
}

///|
/// Chromosome accessor.
pub fn ChipSeqRead::chr(self : ChipSeqRead) -> String {
  self.chr
}

///|
/// Position accessor.
pub fn ChipSeqRead::position(self : ChipSeqRead) -> Int {
  self.position
}

///|
/// Strand accessor.
pub fn ChipSeqRead::strand(self : ChipSeqRead) -> String {
  self.strand
}

// ============================================================================
// CandidatePeak constructors and accessors
// ============================================================================

///|
/// Construct a `CandidatePeak`. `fdr` defaults to 1.0 and `significant`
/// defaults to `false` so the caller can populate them later via FDR
/// correction and filtering.
pub fn CandidatePeak::new(
  chr~ : String,
  start~ : Int,
  end~ : Int,
  summit~ : Int,
  read_count~ : Int,
  control_count~ : Int,
  fold_enrichment~ : Double,
  p_value~ : Double,
) -> CandidatePeak {
  CandidatePeak::{
    chr,
    start,
    end,
    summit,
    read_count,
    control_count,
    fold_enrichment,
    p_value,
    fdr: 1.0,
    significant: false,
  }
}

///|
/// Chromosome accessor.
pub fn CandidatePeak::chr(self : CandidatePeak) -> String {
  self.chr
}

///|
/// Start position accessor.
pub fn CandidatePeak::start(self : CandidatePeak) -> Int {
  self.start
}

///|
/// End position accessor.
pub fn CandidatePeak::end(self : CandidatePeak) -> Int {
  self.end
}

///|
/// Summit position accessor.
pub fn CandidatePeak::summit(self : CandidatePeak) -> Int {
  self.summit
}

///|
/// Treatment read count accessor.
pub fn CandidatePeak::read_count(self : CandidatePeak) -> Int {
  self.read_count
}

///|
/// Control read count accessor.
pub fn CandidatePeak::control_count(self : CandidatePeak) -> Int {
  self.control_count
}

///|
/// Fold enrichment accessor.
pub fn CandidatePeak::fold_enrichment(self : CandidatePeak) -> Double {
  self.fold_enrichment
}

///|
/// P-value accessor.
pub fn CandidatePeak::p_value(self : CandidatePeak) -> Double {
  self.p_value
}

///|
/// Peak width (end - start).
pub fn CandidatePeak::width(self : CandidatePeak) -> Int {
  self.end - self.start
}

// ============================================================================
// PeakCallingParams constructors
// ============================================================================

///|
/// Default parameters (window 500, step 100, local lambda 10kb,
/// p-value 1e-5, FDR 0.05, min fold 2.0).
pub fn PeakCallingParams::default() -> PeakCallingParams {
  PeakCallingParams::{
    window_size: 500,
    step_size: 100,
    local_lambda_size: 10000,
    pvalue_threshold: 1.0e-5,
    fdr_threshold: 0.05,
    min_fold_enrichment: 2.0,
  }
}

///|
/// Construct `PeakCallingParams` with optional overrides for any
/// field. Unspecified fields take the default values.
pub fn PeakCallingParams::new(
  window_size? : Int = 500,
  step_size? : Int = 100,
  local_lambda_size? : Int = 10000,
  pvalue_threshold? : Double = 1.0e-5,
  fdr_threshold? : Double = 0.05,
  min_fold_enrichment? : Double = 2.0,
) -> PeakCallingParams {
  PeakCallingParams::{
    window_size,
    step_size,
    local_lambda_size,
    pvalue_threshold,
    fdr_threshold,
    min_fold_enrichment,
  }
}

// ============================================================================
// Math helpers: incomplete gamma functions
// ============================================================================
//
// Poisson survival P(X >= k; lambda) equals the regularized lower
// incomplete gamma function P(k, lambda) = gamma(k, lambda) / Gamma(k).
// We compute it via the standard Numerical Recipes series / continued
// fraction pair, with a Lanczos approximation for ln(Gamma).

///|
/// Lanczos approximation for ln(Gamma(x)).
fn pc_gammln(xx : Double) -> Double {
  let cof = [
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.1208650973866179e-2,
    -0.5395239384953e-5,
  ]
  let x = xx
  let mut y = xx
  let mut tmp = x + 5.5
  tmp = tmp - (x + 0.5) * @math.ln(tmp)
  let mut ser = 1.000000000190015
  let mut j = 0
  while j < 6 {
    y = y + 1.0
    ser = ser + cof[j] / y
    j = j + 1
  }
  -tmp + @math.ln(2.5066282746310005 * ser / x)
}

///|
/// Regularized lower incomplete gamma function P(a, x).
/// Uses series representation when x < a + 1, otherwise uses
/// the continued fraction for Q(a, x) and returns 1 - Q.
fn pc_gammp(a : Double, x : Double) -> Double {
  if x < 0.0 || a <= 0.0 {
    return 0.0
  }
  if x == 0.0 {
    return 0.0
  }
  if x < a + 1.0 {
    pc_gser(a, x)
  } else {
    1.0 - pc_gcf(a, x)
  }
}

///|
/// Series representation of the lower regularized incomplete gamma.
fn pc_gser(a : Double, x : Double) -> Double {
  let gln = pc_gammln(a)
  let mut ap = a
  let mut sum = 1.0 / a
  let mut del = sum
  let mut n = 0
  while n < 300 {
    ap = ap + 1.0
    del = del * x / ap
    sum = sum + del
    if del.abs() < sum.abs() * 1.0e-14 {
      break
    }
    n = n + 1
  }
  sum * @math.exp(-x + a * @math.ln(x) - gln)
}

///|
/// Continued fraction representation of the upper regularized
// incomplete gamma Q(a, x).
fn pc_gcf(a : Double, x : Double) -> Double {
  let gln = pc_gammln(a)
  let mut b = x + 1.0 - a
  let mut c = 1.0e30
  let mut d = 1.0 / b
  let mut h = d
  let mut i = 1
  while i <= 300 {
    let ii = i.to_double()
    let an = -ii * (ii - a)
    b = b + 2.0
    d = an * d + b
    if d.abs() < 1.0e-30 {
      d = 1.0e-30
    }
    c = b + an / c
    if c.abs() < 1.0e-30 {
      c = 1.0e-30
    }
    d = 1.0 / d
    let del = d * c
    h = h * del
    if (del - 1.0).abs() < 1.0e-14 {
      break
    }
    i = i + 1
  }
  @math.exp(-x + a * @math.ln(x) - gln) * h
}

///|
/// Compute the Poisson survival function P(X >= observed) for
/// X ~ Poisson(expected). For `observed <= 0` the survival is 1.0
/// and for `expected <= 0` with `observed > 0` it is 0.0 (theoretical
/// extreme), but a tiny floor is returned to keep downstream math sane.
pub fn poisson_pvalue(observed : Int, expected : Double) -> Double {
  if observed <= 0 {
    return 1.0
  }
  if expected <= 0.0 {
    return 1.0e-300
  }
  let p = pc_gammp(observed.to_double(), expected)
  if p < 0.0 {
    0.0
  } else if p > 1.0 {
    1.0
  } else {
    p
  }
}

///|
/// Compute fold enrichment = observed / expected. Returns 0.0 when
/// expected is non-positive, and a large value (1.0e6) when observed
/// is positive but expected is effectively zero.
pub fn compute_fold_enrichment(observed : Int, expected : Double) -> Double {
  if expected <= 0.0 {
    return if observed > 0 { 1.0e6 } else { 0.0 }
  }
  observed.to_double() / expected
}

// ============================================================================
// Read counting and local lambda estimation
// ============================================================================

///|
/// Count reads whose 5' position falls in `[start, end)` on `chr`.
pub fn count_reads_in_window(
  reads : Array[ChipSeqRead],
  chr : String,
  start : Int,
  end : Int,
) -> Int {
  let mut count = 0
  for r in reads {
    if r.chr == chr && r.position >= start && r.position < end {
      count = count + 1
    }
  }
  count
}

///|
/// Estimate local background rate (reads per bp) from the control
/// sample in a window of size `window_size` centered on `center`.
pub fn estimate_local_lambda(
  control : Array[ChipSeqRead],
  chr : String,
  center : Int,
  window_size : Int,
) -> Double {
  if window_size <= 0 {
    return 0.0
  }
  let half = window_size / 2
  let start = center - half
  let end = center + half
  let count = count_reads_in_window(control, chr, start, end)
  count.to_double() / window_size.to_double()
}

// ============================================================================
// Internal helpers
// ============================================================================

///|
/// Collect the unique set of chromosomes that appear in either the
/// treatment or control read set, preserving first-seen order.
fn pc_get_chromosomes(
  treatment : Array[ChipSeqRead],
  control : Array[ChipSeqRead],
) -> Array[String] {
  let chrs : Array[String] = Array::new()
  let push_if_new = fn(c : String) -> Unit {
    let mut found = false
    for existing in chrs {
      if existing == c {
        found = true
        break
      }
    }
    if !found {
      chrs.push(c)
    }
  }
  for r in treatment {
    push_if_new(r.chr)
  }
  for r in control {
    push_if_new(r.chr)
  }
  chrs
}

///|
/// Maximum read position on a given chromosome across both libraries.
/// Returns -1 when no reads are present on the chromosome.
fn pc_max_position(
  treatment : Array[ChipSeqRead],
  control : Array[ChipSeqRead],
  chr : String,
) -> Int {
  let mut max_pos = -1
  for r in treatment {
    if r.chr == chr && r.position > max_pos {
      max_pos = r.position
    }
  }
  for r in control {
    if r.chr == chr && r.position > max_pos {
      max_pos = r.position
    }
  }
  max_pos
}

// ============================================================================
// Main peak calling pipeline
// ============================================================================

///|
/// Main MACS2-like peak calling entry point.
///
/// Steps:
/// 1. Build a chromosome list from both libraries.
/// 2. Compute a global scaling factor as treatment_total/control_total.
/// 3. Slide a fixed-size window across each chromosome.
/// 4. For each window with treatment reads, estimate the local
///    background from the control sample and compute the expected
///    count as `lambda_rate * window_size * scaling`.
/// 5. Compute a Poisson survival p-value and fold enrichment.
/// 6. Keep windows passing the p-value and fold enrichment thresholds
///    as candidate peaks.
/// 7. Merge nearby candidate peaks.
/// 8. Apply BH-FDR correction and filter.
pub fn call_peaks(
  treatment : Array[ChipSeqRead],
  control : Array[ChipSeqRead],
  params : PeakCallingParams,
) -> Array[CandidatePeak] {
  let chrs = pc_get_chromosomes(treatment, control)
  let treatment_total = treatment.length().to_double()
  let control_total = control.length().to_double()
  let scaling = if control_total > 0.0 {
    treatment_total / control_total
  } else {
    1.0
  }
  let win_size_d = params.window_size.to_double()

  let candidates : Array[CandidatePeak] = Array::new()

  for chr in chrs {
    let max_pos = pc_max_position(treatment, control, chr)
    if max_pos < 0 {
      continue
    }
    let mut pos = 0
    while pos <= max_pos {
      let win_start = pos
      let win_end = pos + params.window_size
      let center = pos + params.window_size / 2

      let treat_count = count_reads_in_window(
        treatment, chr, win_start, win_end,
      )
      if treat_count == 0 {
        pos = pos + params.step_size
        continue
      }

      // Local background from control, scaled to treatment depth.
      let lambda_rate = estimate_local_lambda(
        control, chr, center, params.local_lambda_size,
      )
      let expected_raw = lambda_rate * win_size_d * scaling
      // Floor expected at a small positive value to avoid div-by-zero.
      let expected = if expected_raw > 0.0 { expected_raw } else { 0.5 }

      let pval = poisson_pvalue(treat_count, expected)
      if pval >= params.pvalue_threshold {
        pos = pos + params.step_size
        continue
      }

      let fold = compute_fold_enrichment(treat_count, expected)
      if fold < params.min_fold_enrichment {
        pos = pos + params.step_size
        continue
      }

      let ctrl_count = count_reads_in_window(
        control, chr, win_start, win_end,
      )

      candidates.push(
        CandidatePeak::new(
          chr~,
          start=win_start,
          end=win_end,
          summit=center,
          read_count=treat_count,
          control_count=ctrl_count,
          fold_enrichment=fold,
          p_value=pval,
        ),
      )

      pos = pos + params.step_size
    }
  }

  // Merge overlapping / nearby peaks. Use a gap equal to twice the
  // step size as the merge tolerance.
  let merge_gap = params.step_size * 2
  let merged = merge_peaks(candidates, merge_gap)

  // Apply BH-FDR and mark significant peaks.
  let with_fdr = apply_fdr(merged, params.fdr_threshold)

  // Final filtering by p-value, fold, and FDR.
  filter_peaks(with_fdr, params)
}

// ============================================================================
// Peak merging
// ============================================================================

///|
/// Merge candidate peaks on the same chromosome whose start/end
/// coordinates overlap or are within `max_gap` bp of each other.
///
/// The merged peak takes the union span, the maximum read count
/// position (summit), summed treatment/control counts, the minimum
/// p-value, and the maximum fold enrichment across the merged set.
pub fn merge_peaks(
  peaks : Array[CandidatePeak],
  max_gap : Int,
) -> Array[CandidatePeak] {
  let n = peaks.length()
  if n == 0 {
    return Array::new()
  }
  // Sort by (chr, start) using a simple stable comparator.
  let sorted = peaks.copy()
  sorted.sort_by(fn(a : CandidatePeak, b : CandidatePeak) -> Int {
    if a.chr < b.chr {
      -1
    } else if a.chr > b.chr {
      1
    } else if a.start < b.start {
      -1
    } else if a.start > b.start {
      1
    } else {
      0
    }
  })

  let merged : Array[CandidatePeak] = Array::new()
  let mut cur = sorted[0]
  let mut i = 1
  while i < n {
    let nxt = sorted[i]
    let same_chr = nxt.chr == cur.chr
    let adjacent = same_chr && nxt.start <= cur.end + max_gap
    if same_chr && adjacent {
      // Merge nxt into cur.
      let new_start = if cur.start < nxt.start { cur.start } else { nxt.start }
      let new_end = if cur.end > nxt.end { cur.end } else { nxt.end }
      // Summit: position of the peak with the higher read count.
      let new_summit = if cur.read_count >= nxt.read_count {
        cur.summit
      } else {
        nxt.summit
      }
      let new_read_count = cur.read_count + nxt.read_count
      let new_control_count = cur.control_count + nxt.control_count
      let new_pval = if cur.p_value < nxt.p_value {
        cur.p_value
      } else {
        nxt.p_value
      }
      let new_fold = if cur.fold_enrichment > nxt.fold_enrichment {
        cur.fold_enrichment
      } else {
        nxt.fold_enrichment
      }
      let cur_chr = cur.chr
      cur = CandidatePeak::new(
        chr=cur_chr,
        start=new_start,
        end=new_end,
        summit=new_summit,
        read_count=new_read_count,
        control_count=new_control_count,
        fold_enrichment=new_fold,
        p_value=new_pval,
      )
    } else {
      merged.push(cur)
      cur = nxt
    }
    i = i + 1
  }
  merged.push(cur)
  merged
}

// ============================================================================
// FDR correction (Benjamini-Hochberg)
// ============================================================================

///|
/// Apply Benjamini-Hochberg FDR correction to a list of peaks.
/// Each peak's `fdr` field is set to its adjusted p-value, and
/// `significant` is set to `true` when `fdr <= fdr_threshold`.
pub fn apply_fdr(
  peaks : Array[CandidatePeak],
  fdr_threshold : Double,
) -> Array[CandidatePeak] {
  let n = peaks.length()
  // Initialize FDR to 1.0 (no significance).
  let mut i = 0
  while i < n {
    peaks[i].fdr = 1.0
    peaks[i].significant = false
    i = i + 1
  }
  if n == 0 {
    return peaks
  }

  // Build (p_value, original_index) pairs and sort by p_value ascending.
  let indexed : Array[(Double, Int)] = Array::new()
  let mut j = 0
  while j < n {
    indexed.push((peaks[j].p_value, j))
    j = j + 1
  }
  indexed.sort_by(fn(a : (Double, Int), b : (Double, Int)) -> Int {
    if a.0 < b.0 {
      -1
    } else if a.0 > b.0 {
      1
    } else {
      0
    }
  })

  // Compute BH-adjusted p-values: p_i * n / rank.
  let adj : Array[Double] = Array::make(n, 1.0)
  let mut k = 0
  while k < n {
    let rank = k + 1
    let val = indexed[k].0 * n.to_double() / rank.to_double()
    adj[k] = if val < 1.0 { val } else { 1.0 }
    k = k + 1
  }
  // Enforce monotonicity from the largest rank down.
  let mut m = n - 2
  while m >= 0 {
    if adj[m] > adj[m + 1] {
      adj[m] = adj[m + 1]
    }
    m = m - 1
  }

  // Map back to original order and update peaks in place.
  let mut p = 0
  while p < n {
    let orig_idx = indexed[p].1
    peaks[orig_idx].fdr = adj[p]
    peaks[orig_idx].significant = adj[p] <= fdr_threshold
    p = p + 1
  }
  peaks
}

// ============================================================================
// Filtering
// ============================================================================

///|
/// Filter peaks by the thresholds stored in `params`. A peak is
/// retained when it satisfies all of:
/// - p-value <= `pvalue_threshold`
/// - fold_enrichment >= `min_fold_enrichment`
/// - `fdr <= fdr_threshold` (significant flag set by `apply_fdr`)
pub fn filter_peaks(
  peaks : Array[CandidatePeak],
  params : PeakCallingParams,
) -> Array[CandidatePeak] {
  let out : Array[CandidatePeak] = Array::new()
  for p in peaks {
    if p.p_value <= params.pvalue_threshold &&
       p.fold_enrichment >= params.min_fold_enrichment &&
       p.fdr <= params.fdr_threshold {
      out.push(p)
    }
  }
  out
}

// ============================================================================
// Sample data and summary
// ============================================================================

///|
/// Generate a small synthetic ChIP-seq dataset with three enriched
/// regions on chr1 and a matched control sample with background
/// reads. Returns `(treatment, control)`.
pub fn peak_calling_sample_data() -> (Array[ChipSeqRead], Array[ChipSeqRead]) {
  let treatment : Array[ChipSeqRead] = Array::new()
  let control : Array[ChipSeqRead] = Array::new()

  // Three enriched regions on chr1.
  // Region A: ~1000-1500 (strong, 30 reads)
  // Region B: ~5000-5500 (medium, 20 reads)
  // Region C: ~9000-9500 (weak, 12 reads)
  let region_starts = [1000, 5000, 9000]
  let region_counts = [30, 20, 12]
  let mut r = 0
  while r < region_starts.length() {
    let base = region_starts[r]
    let n = region_counts[r]
    let mut k = 0
    while k < n {
      // Spread reads across a ~500bp window.
      treatment.push(
        ChipSeqRead::new(chr="chr1", position=base + (k * 17) % 500, strand="+"),
      )
      k = k + 1
    }
    r = r + 1
  }

  // Background treatment reads scattered across the chromosome.
  let mut b = 0
  while b < 40 {
    treatment.push(
      ChipSeqRead::new(chr="chr1", position=b * 250, strand="+"),
    )
    b = b + 1
  }

  // Control: background reads across the chromosome with no strong peaks.
  let mut c = 0
  while c < 80 {
    control.push(
      ChipSeqRead::new(chr="chr1", position=c * 130, strand="+"),
    )
    c = c + 1
  }

  (treatment, control)
}

///|
/// Build a human-readable summary string for a list of peaks:
/// total peaks, number of significant peaks, and the top 5 by p-value.
pub fn peak_calling_summary(peaks : Array[CandidatePeak]) -> String {
  let n = peaks.length()
  if n == 0 {
    return "Peak calling summary: 0 peaks called."
  }
  let mut n_sig = 0
  for p in peaks {
    if p.significant {
      n_sig = n_sig + 1
    }
  }
  let mut buf = "Peak calling summary:\n"
  buf = buf + "  Total peaks: " + n.to_string() + "\n"
  buf = buf + "  Significant peaks: " + n_sig.to_string() + "\n"
  // Sort a copy by p-value ascending for "top 5".
  let copy = peaks.copy()
  copy.sort_by(fn(a : CandidatePeak, b : CandidatePeak) -> Int {
    if a.p_value < b.p_value {
      -1
    } else if a.p_value > b.p_value {
      1
    } else {
      0
    }
  })
  let limit = if n < 5 { n } else { 5 }
  buf = buf + "  Top " + limit.to_string() + " peaks (by p-value):\n"
  let mut i = 0
  while i < limit {
    let p = copy[i]
    buf = buf +
      "    " + p.chr + ":" + p.start.to_string() + "-" + p.end.to_string() +
      " summit=" + p.summit.to_string() +
      " reads=" + p.read_count.to_string() +
      " ctrl=" + p.control_count.to_string() +
      " fold=" + p.fold_enrichment.to_string() +
      " p=" + p.p_value.to_string() +
      " fdr=" + p.fdr.to_string() + "\n"
    i = i + 1
  }
  buf
}
