///|
/// Bioconductor stageR - two-stage hypothesis testing with OFDR control.
///
/// Mirrors the Bioconductor `stageR` package (Van den Berge et al., 2017,
/// Genome Biology 18:151). The stage-wise testing paradigm controls the
/// Overall False Discovery Rate (OFDR) at the gene level while allowing
/// post-hoc assessment of individual confirmation hypotheses (e.g.
/// transcript-level or contrast-level tests).
///
/// Two-stage procedure:
///   1. **Screening stage**: gene-level p-values (aggregated from
///      transcript-level p-values via Simes' method) are BH-adjusted
///      and tested at level alpha.
///   2. **Confirmation stage**: for genes passing screening, confirmation
///      p-values are adjusted using a Holm-like step-down procedure with
///      a modified multiplier vector (accounting for the fact that the
///      screening rejection implies at least one null is false).
///
/// The confirmation-stage p-values are rescaled by G/R (total genes /
/// significant genes) so they can be directly compared to `alpha`.

// ============================================================================
// Configuration and result types
// ============================================================================

///|
/// FWER correction method for the confirmation stage.
pub enum StageRMethod {
  None // no correction within gene
  Holm // standard Holm with screening adjustment
  Dte // differential transcript expression (same as Holm)
  Dtu // differential transcript usage (screening implies >=2 false)
}

///|
pub fn stage_r_method_none() -> StageRMethod {
  StageRMethod::None
}

///|
pub fn stage_r_method_holm() -> StageRMethod {
  StageRMethod::Holm
}

///|
pub fn stage_r_method_dte() -> StageRMethod {
  StageRMethod::Dte
}

///|
pub fn stage_r_method_dtu() -> StageRMethod {
  StageRMethod::Dtu
}

///|
/// Configuration for stage-wise testing.
pub struct StageRConfig {
  /// Target OFDR level (e.g. 0.05).
  alpha : Double
  /// Whether screening p-values are already BH-adjusted.
  pScreenAdjusted : Bool
  /// FWER correction method for confirmation stage.
  stageMethod : StageRMethod
}

///|
pub fn StageRConfig::new(
  alpha? : Double = 0.05,
  pScreenAdjusted? : Bool = false,
  stageMethod? : StageRMethod = Holm,
) -> StageRConfig {
  StageRConfig::{ alpha, pScreenAdjusted, stageMethod }
}

///|
/// Result of stage-wise adjustment.
pub struct StageRResult {
  /// Adjusted screening p-values (length G).
  pAdjScreen : Array[Double]
  /// Adjusted confirmation p-values (G × H matrix).
  /// Rows for non-significant genes are filled with -1.0 (NA sentinel).
  pAdjConfirmation : Array[Array[Double]]
  /// FWER correction method used.
  stageMethod : StageRMethod
  /// Target OFDR level.
  alpha : Double
  /// Adjusted alpha level for confirmation stage = (R/G) * alpha.
  alphaAdjusted : Double
  /// Number of significant genes (passing screening).
  nSignificantGenes : Int
  /// Total number of genes.
  nGenes : Int
  /// Number of confirmation hypotheses per gene.
  nHypotheses : Int
}

// ============================================================================
// BH (Benjamini-Hochberg) FDR adjustment
// ============================================================================

///|
/// Benjamini-Hochberg FDR adjustment.
/// Returns adjusted p-values in the original order.
pub fn stage_r_bh_adjust(pvalues : Array[Double]) -> Array[Double] {
  let n = pvalues.length()
  if n == 0 {
    return []
  }
  // Create index-value pairs and sort by p-value ascending.
  let pairs : Array[(Int, Double)] = []
  let mut i = 0
  while i < n {
    pairs.push((i, pvalues[i]))
    i = i + 1
  }
  // Sort by p-value ascending (bubble sort).
  let mut a = 0
  while a < n {
    let mut b = 0
    while b < n - 1 - a {
      if pairs[b].1 > pairs[b + 1].1 {
        let tmp = pairs[b]
        pairs[b] = pairs[b + 1]
        pairs[b + 1] = tmp
      }
      b = b + 1
    }
    a = a + 1
  }
  // Compute adjusted p-values from largest to smallest.
  let adj_sorted : Array[Double] = Array::make(n, 0.0)
  adj_sorted[n - 1] = pairs[n - 1].1 * n.to_double()
  if adj_sorted[n - 1] > 1.0 {
    adj_sorted[n - 1] = 1.0
  }
  let mut i = n - 2
  while i >= 0 {
    let rank = (i + 1).to_double()
    let raw_adj = pairs[i].1 * n.to_double() / rank
    adj_sorted[i] = if raw_adj < adj_sorted[i + 1] { raw_adj } else { adj_sorted[i + 1] }
    if adj_sorted[i] > 1.0 {
      adj_sorted[i] = 1.0
    }
    i = i - 1
  }
  // Reorder to original indices.
  let result : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    result[pairs[i].0] = adj_sorted[i]
    i = i + 1
  }
  result
}

// ============================================================================
// Simes' method for p-value aggregation
// ============================================================================

///|
/// Simes' method: aggregate H p-values into a single gene-level p-value.
/// p_gene = min_i ( H * p_(i) / i ) where p_(i) is the i-th sorted p-value.
pub fn simes_aggregate(pvalues : Array[Double]) -> Double {
  let n = pvalues.length()
  if n == 0 {
    return 1.0
  }
  let sorted = sort_ascending(pvalues)
  let mut min_val = 1.0
  let mut i = 0
  while i < n {
    let val = n.to_double() * sorted[i] / (i + 1).to_double()
    if val < min_val {
      min_val = val
    }
    i = i + 1
  }
  if min_val > 1.0 {
    1.0
  } else {
    min_val
  }
}

///|
/// Sort an array of doubles in ascending order (returns a new array).
fn sort_ascending(arr : Array[Double]) -> Array[Double] {
  let result : Array[Double] = []
  let mut i = 0
  while i < arr.length() {
    result.push(arr[i])
    i = i + 1
  }
  let n = result.length()
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n - 1 - i {
      if result[j] > result[j + 1] {
        let tmp = result[j]
        result[j] = result[j + 1]
        result[j + 1] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

// ============================================================================
// Holm step-down adjustment (with screening adjustment)
// ============================================================================

///|
/// Holm-like step-down adjustment for confirmation p-values.
/// `adjustment` is the multiplier vector (length n).
/// Returns adjusted p-values in the original order.
fn holm_stepdown(
  pvalues : Array[Double],
  adjustment : Array[Double],
) -> Array[Double] {
  let n = pvalues.length()
  if n == 0 {
    return []
  }
  // Create index-value pairs and sort by p-value ascending.
  let pairs : Array[(Int, Double)] = []
  let mut i = 0
  while i < n {
    pairs.push((i, pvalues[i]))
    i = i + 1
  }
  let mut a = 0
  while a < n {
    let mut b = 0
    while b < n - 1 - a {
      if pairs[b].1 > pairs[b + 1].1 {
        let tmp = pairs[b]
        pairs[b] = pairs[b + 1]
        pairs[b + 1] = tmp
      }
      b = b + 1
    }
    a = a + 1
  }
  // Multiply by adjustment vector and take cumulative max.
  let adj_sorted : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let val = pairs[i].1 * adjustment[i]
    adj_sorted[i] = if val > 1.0 { 1.0 } else { val }
    if i > 0 && adj_sorted[i] < adj_sorted[i - 1] {
      adj_sorted[i] = adj_sorted[i - 1]
    }
    i = i + 1
  }
  // Reorder to original indices.
  let result : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    result[pairs[i].0] = adj_sorted[i]
    i = i + 1
  }
  result
}

///|
/// Compute the adjustment multiplier vector for a given method.
fn get_adjustment_vector(method : StageRMethod, n : Int) -> Array[Double] {
  if method is None {
    let adj : Array[Double] = []
    let mut i = 0
    while i < n {
      adj.push(1.0)
      i = i + 1
    }
    return adj
  }
  if method is Holm || method is Dte {
    // adjustment = [n-1, n-1, n-2, ..., 1] (but at least 1 for each)
    // Actually: c(n-1, (n-1):1) = [n-1, n-1, n-2, ..., 2, 1]
    let adj : Array[Double] = []
    if n == 1 {
      adj.push(0.0) // screening implies the single hypothesis is false
      return adj
    }
    adj.push((n - 1).to_double())
    let mut i = n - 1
    while i >= 1 {
      adj.push(i.to_double())
      i = i - 1
    }
    return adj
  }
  if method is Dtu {
    // DTU: screening implies >=2 hypotheses are false.
    // adjustment = c(n-2, n-2, (n-2):1) for n >= 3
    let adj : Array[Double] = []
    if n <= 2 {
      // For n=2: both hypotheses false, adjustment = [0, 0]
      adj.push(0.0)
      adj.push(0.0)
      return adj
    }
    adj.push((n - 2).to_double())
    adj.push((n - 2).to_double())
    let mut i = n - 2
    while i >= 1 {
      adj.push(i.to_double())
      i = i - 1
    }
    return adj
  }
  // Default: no adjustment
  let adj : Array[Double] = []
  let mut i = 0
  while i < n {
    adj.push(1.0)
    i = i + 1
  }
  adj
}

// ============================================================================
// Main stage-wise adjustment
// ============================================================================

///|
/// Perform the two-stage adjustment.
///
/// Parameters:
///   - `pScreen`: screening p-values (length G, one per gene).
///   - `pConfirmation`: confirmation p-value matrix (G × H).
///   - `config`: configuration (alpha, method, pScreenAdjusted).
pub fn stage_wise_adjustment(
  pScreen : Array[Double],
  pConfirmation : Array[Array[Double]],
  config : StageRConfig,
) -> StageRResult {
  let g = pScreen.length()
  let h = if g > 0 { pConfirmation[0].length() } else { 0 }

  // === SCREENING STAGE ===
  let pAdjScreen = if config.pScreenAdjusted {
    let copy : Array[Double] = []
    let mut i = 0
    while i < g {
      copy.push(pScreen[i])
      i = i + 1
    }
    copy
  } else {
    stage_r_bh_adjust(pScreen)
  }

  // Determine significant genes.
  let mut nSig = 0
  let sig_mask : Array[Bool] = Array::make(g, false)
  let mut i = 0
  while i < g {
    if pAdjScreen[i] <= config.alpha {
      sig_mask[i] = true
      nSig = nSig + 1
    }
    i = i + 1
  }

  // Adjusted alpha level.
  let alphaAdjusted = if g > 0 {
    nSig.to_double() / g.to_double() * config.alpha
  } else {
    0.0
  }

  // === CONFIRMATION STAGE ===
  // Adjustment multiplier vector (same for all genes with same H).
  let adj_vec = get_adjustment_vector(config.stageMethod, h)
  let rescale = if nSig > 0 { g.to_double() / nSig.to_double() } else { 1.0 }

  let pAdjConfirmation : Array[Array[Double]] = []
  let mut gi = 0
  while gi < g {
    let row : Array[Double] = []
    if sig_mask[gi] {
      // Adjust confirmation p-values for this gene.
      let adjusted = holm_stepdown(pConfirmation[gi], adj_vec)
      // Rescale so they can be compared to alpha directly.
      let mut hi = 0
      while hi < h {
        let val = adjusted[hi] * rescale
        row.push(if val > 1.0 { 1.0 } else { val })
        hi = hi + 1
      }
    } else {
      // Non-significant gene: fill with -1.0 (NA sentinel).
      let mut hi = 0
      while hi < h {
        row.push(-1.0)
        hi = hi + 1
      }
    }
    pAdjConfirmation.push(row)
    gi = gi + 1
  }

  StageRResult::{
    pAdjScreen,
    pAdjConfirmation,
    stageMethod: config.stageMethod,
    alpha: config.alpha,
    alphaAdjusted,
    nSignificantGenes: nSig,
    nGenes: g,
    nHypotheses: h,
  }
}

// ============================================================================
// Result accessors
// ============================================================================

///|
/// Get the binary significance matrix (G × (H+1)).
/// Column 0 = gene-level significance, columns 1..H = confirmation significance.
pub fn get_results(result : StageRResult) -> Array[Array[Bool]] {
  let g = result.nGenes
  let h = result.nHypotheses
  let mat : Array[Array[Bool]] = []
  let mut i = 0
  while i < g {
    let row : Array[Bool] = []
    row.push(result.pAdjScreen[i] <= result.alpha)
    let mut j = 0
    while j < h {
      row.push(result.pAdjConfirmation[i][j] >= 0.0 && result.pAdjConfirmation[i][j] <= result.alpha)
      j = j + 1
    }
    mat.push(row)
    i = i + 1
  }
  mat
}

///|
/// Get indices of significant genes (passing screening).
pub fn get_significant_genes(result : StageRResult) -> Array[Int] {
  let genes : Array[Int] = []
  let mut i = 0
  while i < result.nGenes {
    if result.pAdjScreen[i] <= result.alpha {
      genes.push(i)
    }
    i = i + 1
  }
  genes
}

///|
/// Get indices of significant confirmation hypotheses across all genes.
/// Returns (gene_index, hypothesis_index) pairs.
pub fn get_significant_hypotheses(result : StageRResult) -> Array[(Int, Int)] {
  let sigs : Array[(Int, Int)] = []
  let mut i = 0
  while i < result.nGenes {
    if result.pAdjScreen[i] <= result.alpha {
      let mut j = 0
      while j < result.nHypotheses {
        if result.pAdjConfirmation[i][j] >= 0.0 && result.pAdjConfirmation[i][j] <= result.alpha {
          sigs.push((i, j))
        }
        j = j + 1
      }
    }
    i = i + 1
  }
  sigs
}

///|
/// Generate a summary string of the stageR result.
pub fn stage_r_summary(result : StageRResult) -> String {
  let buf = StringBuilder::new()
  buf.write_string("StageR Two-Stage Testing Summary\n")
  buf.write_string("  Method: " + stage_r_method_name(result.stageMethod) + "\n")
  buf.write_string("  Target OFDR (alpha): " + result.alpha.to_string() + "\n")
  buf.write_string("  Adjusted alpha: " + result.alphaAdjusted.to_string() + "\n")
  buf.write_string("  Total genes: " + result.nGenes.to_string() + "\n")
  buf.write_string("  Hypotheses per gene: " + result.nHypotheses.to_string() + "\n")
  buf.write_string("  Significant genes (screening): " + result.nSignificantGenes.to_string() + "\n")
  let sig_hyps = get_significant_hypotheses(result)
  buf.write_string("  Significant hypotheses (confirmation): " + sig_hyps.length().to_string() + "\n")
  buf.to_string()
}

///|
fn stage_r_method_name(method : StageRMethod) -> String {
  match method {
    None => "none"
    Holm => "holm"
    Dte => "dte"
    Dtu => "dtu"
  }
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Generate sample p-values for stageR demo.
/// Returns (pScreen, pConfirmation) where pScreen has 10 genes
/// and pConfirmation has 10×3 (3 hypotheses per gene).
pub fn stage_r_sample() -> (Array[Double], Array[Array[Double]]) {
  let pScreen : Array[Double] = [
    0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.5, 0.8,
  ]
  let pConfirmation : Array[Array[Double]] = []
  // Gene 0: very significant screening, all confirmation significant.
  pConfirmation.push([0.001, 0.002, 0.003])
  // Gene 1: significant screening, 2/3 confirmation significant.
  pConfirmation.push([0.005, 0.01, 0.06])
  // Gene 2: marginally significant, 1/3 confirmation.
  pConfirmation.push([0.01, 0.04, 0.08])
  // Gene 3: significant screening, none confirmation.
  pConfirmation.push([0.03, 0.04, 0.05])
  // Gene 4-9: non-significant or borderline.
  pConfirmation.push([0.1, 0.2, 0.3])
  pConfirmation.push([0.15, 0.25, 0.35])
  pConfirmation.push([0.3, 0.4, 0.5])
  pConfirmation.push([0.4, 0.5, 0.6])
  pConfirmation.push([0.6, 0.7, 0.8])
  pConfirmation.push([0.8, 0.9, 1.0])
  (pScreen, pConfirmation)
}
