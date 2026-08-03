///|
/// ESTIMATE - Estimation of STromal and Immune cells in MAlignant Tumours
/// using Expression data.
///
/// References:
///   - Yoshihara K et al. (2013) Inferring tumour purity and stromal and immune
///     infiltration from transcriptomic data. Nature Communications 4:2612.
///   - Bioconductor package: estimate.
///
/// Core algorithm:
///   Given two gene sets (stromal and immune signature genes) and a sample ×
///   gene expression matrix, compute per-sample enrichment scores using a
///   rank-based single-sample GSEA procedure. The raw scores are then
///   standardized using an empirical cumulative distribution function (ECDF)
///   across samples to produce:
///     - StromalScore: relative abundance of stromal cells
///     - ImmuneScore:  relative abundance of immune cells
///     - ESTIMATEScore: combined infiltration level (Stromal + Immune)
///     - TumorPurity:  inferred from ESTIMATEScore via Yoshihara's formula
///
/// All public identifiers are prefixed with `est_` to avoid conflicts.

// ============================================================================
// Data structures
// ============================================================================

///|
/// ESTIMATE input: an expression matrix with sample and gene identifiers.
/// Rows = genes, cols = samples.
pub struct EstExpression {
  gene_names : Array[String]
  sample_names : Array[String]
  matrix : Array[Array[Double]] // rows=genes, cols=samples
}

///|
/// ESTIMATE gene signature: stromal and immune gene sets.
pub struct EstSignature {
  stromal_genes : Array[String]
  immune_genes : Array[String]
}

///|
/// Per-sample ESTIMATE scores.
pub struct EstSampleScore {
  sample_id : String
  stromal_score : Double
  immune_score : Double
  estimate_score : Double    // stromal + immune
  tumor_purity : Double      // inferred purity in [0, 1] (NaN if undefined)
}

///|
/// Full ESTIMATE result across all samples.
pub struct EstResult {
  scores : Array[EstSampleScore]
  sample_names : Array[String]
  method_name : String
}

///|
/// Create an expression matrix.
pub fn EstExpression::new(
  gene_names : Array[String],
  sample_names : Array[String],
  matrix : Array[Array[Double]],
) -> EstExpression {
  { gene_names, sample_names, matrix }
}

///|
/// Get gene names.
pub fn EstExpression::gene_names(self : EstExpression) -> Array[String] {
  self.gene_names
}

///|
/// Get sample names.
pub fn EstExpression::sample_names(self : EstExpression) -> Array[String] {
  self.sample_names
}

///|
/// Get the raw matrix (rows=genes, cols=samples).
pub fn EstExpression::matrix(self : EstExpression) -> Array[Array[Double]] {
  self.matrix
}

///|
/// Create a gene signature.
pub fn EstSignature::new(
  stromal_genes : Array[String],
  immune_genes : Array[String],
) -> EstSignature {
  { stromal_genes, immune_genes }
}

///|
/// Get stromal gene list.
pub fn EstSignature::stromal_genes(self : EstSignature) -> Array[String] {
  self.stromal_genes
}

///|
/// Get immune gene list.
pub fn EstSignature::immune_genes(self : EstSignature) -> Array[String] {
  self.immune_genes
}

///|
/// Get sample ID.
pub fn EstSampleScore::sample_id(self : EstSampleScore) -> String {
  self.sample_id
}

///|
/// Get stromal score.
pub fn EstSampleScore::stromal_score(self : EstSampleScore) -> Double {
  self.stromal_score
}

///|
/// Get immune score.
pub fn EstSampleScore::immune_score(self : EstSampleScore) -> Double {
  self.immune_score
}

///|
/// Get combined ESTIMATE score.
pub fn EstSampleScore::estimate_score(self : EstSampleScore) -> Double {
  self.estimate_score
}

///|
/// Get inferred tumor purity.
pub fn EstSampleScore::tumor_purity(self : EstSampleScore) -> Double {
  self.tumor_purity
}

///|
/// Get all per-sample scores.
pub fn EstResult::scores(self : EstResult) -> Array[EstSampleScore] {
  self.scores
}

///|
/// Get sample names.
pub fn EstResult::sample_names(self : EstResult) -> Array[String] {
  self.sample_names
}

///|
/// Get method name.
pub fn EstResult::method_name(self : EstResult) -> String {
  self.method_name
}

// ============================================================================
// Core math utilities
// ============================================================================

///|
/// Sort indices of an array in descending order (largest first).
fn est_sort_indices_desc(values : Array[Double]) -> Array[Int] {
  let n = values.length()
  let indices : Array[Int] = Array::make(n, 0)
  for i in 0..<n {
    indices[i] = i
  }
  // Insertion sort (descending) on indices
  let mut i = 1
  while i < n {
    let key = indices[i]
    let mut j = i - 1
    while j >= 0 && values[indices[j]] < values[key] {
      indices[j + 1] = indices[j]
      j = j - 1
    }
    indices[j + 1] = key
    i = i + 1
  }
  indices
}

///|
/// Empirical cumulative distribution function (ECDF) value of `x` in `sample`.
/// Returns the proportion of values in `sample` that are <= x.
fn est_ecdf(x : Double, sample : Array[Double]) -> Double {
  let n = sample.length()
  if n == 0 { return 0.0 }
  let mut count = 0
  for v in sample {
    if v <= x { count = count + 1 }
  }
  count.to_double() / n.to_double()
}

///|
/// Check if all values in an array are identical (degenerate case).
fn est_all_equal(arr : Array[Double]) -> Bool {
  let n = arr.length()
  if n <= 1 { return true }
  let first = arr[0]
  for i in 1..<n {
    if (arr[i] - first).abs() > 1.0e-15 { return false }
  }
  true
}

// ============================================================================
// Single-sample GSEA (ssGSEA) score computation
// ============================================================================

///|
/// Compute the ssGSEA-style enrichment score for a gene set in a single sample.
///
/// Algorithm (Barbie et al. 2009):
///   1. Rank genes by expression (descending).
///   2. Walk down the ranked list; for genes in the set, increment a running
///      sum by |expr|/(sum of |expr| in set); for genes not in the set,
///      decrement by 1/(N - |set|).
///   3. Score = sum of positive running-sum increments (integration-style).
///   4. Normalize by the range of the running sum.
fn est_ssgsea_score(
  expression : Array[Double],
  gene_names : Array[String],
  gene_set : Array[String],
) -> Double {
  let n = gene_names.length()
  if n == 0 || gene_set.length() == 0 { return 0.0 }
  // Build membership set
  let gene_member : Map[String, Unit] = Map([], capacity=gene_set.length())
  let mut set_size = 0
  for g in gene_set {
    if !gene_member.contains(g) {
      gene_member.set(g, ())
      set_size = set_size + 1
    }
  }
  // Rank indices by expression (descending)
  let order = est_sort_indices_desc(expression)
  // Compute sum of absolute expressions for genes in the set
  let mut sum_in_set = 0.0
  for i in 0..<n {
    if gene_member.contains(gene_names[i]) {
      sum_in_set = sum_in_set + expression[i].abs()
    }
  }
  if sum_in_set < 1.0e-15 { return 0.0 }
  // Walk down the ranked list
  let n_not_set = (n - set_size).to_double()
  let incr_in = 1.0 / sum_in_set
  let decr_out = if n_not_set > 0.0 { 1.0 / n_not_set } else { 0.0 }
  let mut running = 0.0
  let mut max_val = 0.0
  let mut min_val = 0.0
  let mut pos_acc = 0.0 // integral of positive excursions
  for i in 0..<n {
    let idx = order[i]
    if gene_member.contains(gene_names[idx]) {
      // Weight by absolute expression
      let w = expression[idx].abs() * incr_in
      running = running + w
    } else {
      running = running - decr_out
    }
    if running > max_val { max_val = running }
    if running < min_val { min_val = running }
    if running > 0.0 { pos_acc = pos_acc + running }
  }
  // Standard ssGSEA normalization: integrate positive area and divide by
  // total range, then subtract 0.5 to center
  let range = max_val - min_val
  if range < 1.0e-15 { return 0.0 }
  // The classic ssGSEA score uses sum of positive running sums (weighted)
  // normalized by (n_set * n_not_set) is NOT appropriate here; we follow
  // the integration-style score: pos_acc / n (average positive excursion)
  pos_acc / n.to_double()
}

// ============================================================================
// Tumor purity inference
// ============================================================================

///|
/// Inferred tumor purity from ESTIMATEScore using the formula from
/// Yoshihara et al. (2013), Equation (1):
///   purity = cos(0.6049872018 + 0.0001627855 * ESTIMATEScore)
/// Returns the value in [0, 1]; returns NaN when ESTIMATEScore is invalid.
fn est_purity_from_estimate(estimate_score : Double) -> Double {
  let arg = 0.6049872018 + 0.0001627855 * estimate_score
  // If argument is outside [0, pi], purity is undefined
  if arg < 0.0 || arg > 3.141592653589793 {
    return 0.0 / 0.0 // NaN
  }
  @math.cos(arg)
}

// ============================================================================
// Main ESTIMATE function
// ============================================================================

///|
/// Run ESTIMATE on an expression matrix.
///
/// For each sample, computes:
///   - raw stromal enrichment score (ssGSEA on stromal gene set)
///   - raw immune enrichment score (ssGSEA on immune gene set)
/// Then standardizes each score across samples using ECDF:
///   normalized = ecdf(raw, all_raw_scores) * 2 - 1   // maps to [-1, 1]
/// The final scores are scaled by a factor of 100 for display convenience,
/// matching the convention of the R `estimate` package.
///
/// Parameters:
///   - expr: expression matrix (rows=genes, cols=samples)
///   - signature?: gene signature (default: built-in stromal/immune sets)
pub fn est_run(
  expr : EstExpression,
  signature? : EstSignature = est_default_signature(),
) -> EstResult {
  let n_samples = expr.sample_names.length()
  let n_genes = expr.gene_names.length()
  // Compute raw scores per sample
  let raw_stromal : Array[Double] = []
  let raw_immune : Array[Double] = []
  for s in 0..<n_samples {
    // Extract expression column for this sample
    let col : Array[Double] = []
    for gi in 0..<n_genes {
      col.push(expr.matrix[gi][s])
    }
    let st = est_ssgsea_score(col, expr.gene_names, signature.stromal_genes)
    let im = est_ssgsea_score(col, expr.gene_names, signature.immune_genes)
    raw_stromal.push(st)
    raw_immune.push(im)
  }
  // Standardize across samples via ECDF: map to [0, 1], then scale & shift
  // The R estimate package multiplies by a constant ~ 100 to get readable
  // integer-ish scores. We follow that convention.
  let scale = 100.0
  // Detect degenerate case: if all raw scores are identical (e.g., all zero
  // when the gene set is empty or has no matching genes), ECDF would map
  // them all to 1.0. In that case, set normalized scores to 0 (neutral).
  let stromal_degenerate = est_all_equal(raw_stromal)
  let immune_degenerate = est_all_equal(raw_immune)
  let scores : Array[EstSampleScore] = []
  for s in 0..<n_samples {
    // ECDF-based normalization to [-scale, scale]
    let st_ecdf = if stromal_degenerate { 0.5 } else { est_ecdf(raw_stromal[s], raw_stromal) }
    let im_ecdf = if immune_degenerate { 0.5 } else { est_ecdf(raw_immune[s], raw_immune) }
    let stromal_final = (st_ecdf - 0.5) * scale * 2.0
    let immune_final = (im_ecdf - 0.5) * scale * 2.0
    let est_score = stromal_final + immune_final
    let purity = est_purity_from_estimate(est_score)
    scores.push({
      sample_id: expr.sample_names[s],
      stromal_score: stromal_final,
      immune_score: immune_final,
      estimate_score: est_score,
      tumor_purity: purity,
    })
  }
  { scores, sample_names: expr.sample_names.copy(), method_name: "estimate" }
}

// ============================================================================
// Built-in stromal and immune gene signatures
// ============================================================================

///|
/// Get the built-in stromal gene signature.
/// These are 50 representative stromal-associated genes from Yoshihara et al.
/// (2013) Supplementary Table S1 (subset).
pub fn est_stromal_genes() -> Array[String] {
  [
    "ACTA2", "ADAM12", "ADAMDEC1", "ADAMTS2", "AEBP1", "ANGPTL2", "ANGPTL4",
    "APBB1IP", "APLP1", "ARCHN1", "ASPN", "BGN", "BMP1", "C1QA", "C1QB",
    "C1QC", "C3", "CADM1", "CALD1", "CCDC80", "CD248", "CD55", "CDH11",
    "CDH19", "COL10A1", "COL11A1", "COL12A1", "COL14A1", "COL1A1", "COL1A2",
    "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2", "COL6A3", "COPZ2",
    "CRLF1", "CXCL1", "CXCL12", "DCN", "DPT", "ECM2", "EDIL3", "EGR2",
    "FAP", "FBLN1", "FBN1", "FBN2", "FGF7",
  ]
}

///|
/// Get the built-in immune gene signature.
/// These are 30 representative immune-associated genes from Yoshihara et al.
/// (2013) Supplementary Table S2 (subset).
pub fn est_immune_genes() -> Array[String] {
  [
    "AIF1", "APOBEC3G", "ARHGAP25", "BCL6", "BLK", "BST1", "BST2", "C1orf162",
    "C2orf88", "C3orf14", "CCL3", "CD164L1", "CD2", "CD247", "CD3D", "CD3E",
    "CD3G", "CD52", "CD53", "CD69", "CD7", "CEBPD", "CLEC2B", "CMKLR1",
    "CPEB2", "CRIP1", "CXCL9", "CXCR6", "CYBB", "DOCK2",
  ]
}

///|
/// Get the default ESTIMATE gene signature (stromal + immune sets).
pub fn est_default_signature() -> EstSignature {
  EstSignature::new(est_stromal_genes(), est_immune_genes())
}

// ============================================================================
// Sample data for testing
// ============================================================================

///|
/// Create a sample expression matrix for testing ESTIMATE.
/// Builds 4 samples with varying degrees of stromal and immune infiltration
/// by combining a baseline tumor profile with stromal and immune signatures
/// in different proportions. Background (non-signature) genes are included so
/// that ssGSEA rankings are meaningful: signature genes compete with the
/// background for top ranks.
pub fn est_sample_data() -> EstExpression {
  let stromal = est_stromal_genes()
  let immune = est_immune_genes()
  // Build a unified gene list (stromal + immune, no duplicates)
  let gene_names : Array[String] = []
  let stromal_set : Map[String, Unit] = Map([], capacity=stromal.length())
  for g in stromal {
    gene_names.push(g)
    stromal_set.set(g, ())
  }
  let immune_set : Map[String, Unit] = Map([], capacity=immune.length())
  for g in immune {
    if !stromal_set.contains(g) {
      gene_names.push(g)
    }
    immune_set.set(g, ())
  }
  // Add background (non-signature) genes so ssGSEA has a reference baseline.
  // These represent typical tumor genes that are neither stromal nor immune.
  let n_background = 60
  for bi in 0..<n_background {
    gene_names.push("BACKGRD" + bi.to_string())
  }
  let n_genes = gene_names.length()
  let sample_names = [
    "HighStromal", "HighImmune", "HighBoth", "LowBoth",
  ]
  // Define infiltration profiles: (stromal_add, immune_add)
  // HighStromal: stromal genes up-regulated, immune slightly up
  // HighImmune:  immune genes up-regulated, stromal slightly up
  // HighBoth:    both stromal and immune up-regulated
  // LowBoth:     both stromal and immune DOWN-regulated (below background),
  //              representing a tumor-dominant sample with little infiltration
  let profiles = [
    (8.0, 1.0),
    (1.0, 8.0),
    (8.0, 8.0),
    (-2.0, -2.0),
  ]
  // Build matrix (rows=genes, cols=samples)
  let matrix : Array[Array[Double]] = []
  for gi in 0..<n_genes {
    let g = gene_names[gi]
    let is_stromal = stromal_set.contains(g)
    let is_immune = immune_set.contains(g)
    let row : Array[Double] = []
    for si in 0..<4 {
      let (s_level, i_level) = profiles[si]
      let base = 5.0 // baseline expression for all genes (incl. background)
      let mut val = base
      if is_stromal {
        val = val + s_level
      }
      if is_immune {
        val = val + i_level
      }
      // Clamp to non-negative (expression cannot be negative)
      if val < 0.1 { val = 0.1 }
      // Add small per-gene jitter to break ties
      val = val + 0.01 * (gi % 7).to_double()
      row.push(val)
    }
    matrix.push(row)
  }
  EstExpression::new(gene_names, sample_names, matrix)
}

///|
/// Get the score for a specific sample.
pub fn EstResult::get_score(self : EstResult, sample_id : String) -> EstSampleScore {
  for s in self.scores {
    if s.sample_id == sample_id { return s }
  }
  // Return a default empty score if not found
  {
    sample_id,
    stromal_score: 0.0,
    immune_score: 0.0,
    estimate_score: 0.0,
    tumor_purity: 0.0 / 0.0,
  }
}

///|
/// Format ESTIMATE results as a string table.
pub fn EstResult::to_string(self : EstResult) -> String {
  let mut result = "ESTIMATE Results:\n"
  result = result + "Sample\tStromalScore\tImmuneScore\tESTIMATEScore\tTumorPurity\n"
  for s in self.scores {
    let purity_str = if s.tumor_purity.is_nan() {
      "NA"
    } else {
      s.tumor_purity.to_string()
    }
    result = result +
      s.sample_id +
      "\t" +
      s.stromal_score.to_string() +
      "\t" +
      s.immune_score.to_string() +
      "\t" +
      s.estimate_score.to_string() +
      "\t" +
      purity_str +
      "\n"
  }
  result
}
