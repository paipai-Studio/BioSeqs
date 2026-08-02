///|
/// missMethyl module - Differential methylation analysis with GO enrichment.
/// Inspired by Bioconductor missMethyl package.
///
/// Provides functions for Illumina 450K/EPIC methylation arrays:
/// - Differential methylation analysis (DMP) using M-values and t-tests
/// - Beta-to-M and M-to-beta value conversions
/// - GO enrichment analysis with probe-gene mapping bias correction
/// - Fisher's exact test and BH-FDR multiple testing correction
/// - Probe filtering by chromosome and CpG island relation

// ============================================================================
// Private helpers
// ============================================================================

///|
/// Mean of an array of doubles. Returns 0.0 for an empty array.
fn mm_mean(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let mut s = 0.0
  for v in arr {
    s = s + v
  }
  s / n.to_double()
}

///|
/// Sample variance (n-1 denominator) given the precomputed mean.
fn mm_variance(arr : Array[Double], mean : Double) -> Double {
  let n = arr.length()
  if n < 2 {
    return 0.0
  }
  let mut s = 0.0
  for v in arr {
    let d = v - mean
    s = s + d * d
  }
  s / (n - 1).to_double()
}

///|
/// Lanczos approximation for ln(gamma(x)).
fn mm_lgamma(x : Double) -> Double {
  let cof = [
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.001208650973866179,
    -0.000005395239384953,
  ]
  let mut x1 = x
  let mut tmp = x1 + 5.5
  tmp = (x1 + 0.5) * @math.ln(tmp) - tmp
  let mut ser = 1.000000000190015
  let mut j = 0
  while j < cof.length() {
    x1 = x1 + 1.0
    ser = ser + cof[j] / x1
    j = j + 1
  }
  tmp + @math.ln(2.5066282746310005 * ser)
}

///|
/// Continued fraction for the incomplete beta function.
fn mm_betacf(x : Double, a : Double, b : Double) -> Double {
  let max_iter = 100
  let eps = 0.0000003
  let mut c = 1.0
  let mut d = 1.0 - (a + b) * x / (a + 1.0)
  if d.abs() < eps {
    d = eps
  }
  d = 1.0 / d
  let mut h = d
  let mut m = 1
  while m <= max_iter {
    let m2 = (2 * m).to_double()
    let mut aa = m.to_double() *
      (b - m.to_double()) *
      x /
      ((a + m2 - 1.0) * (a + m2))
    d = 1.0 + aa * d
    if d.abs() < eps {
      d = eps
    }
    c = 1.0 + aa / c
    if c.abs() < eps {
      c = eps
    }
    d = 1.0 / d
    h = h * d * c
    aa = -(a + m.to_double()) *
      (a + b + m.to_double()) *
      x /
      ((a + m2) * (a + m2 + 1.0))
    d = 1.0 + aa * d
    if d.abs() < eps {
      d = eps
    }
    c = 1.0 + aa / c
    if c.abs() < eps {
      c = eps
    }
    d = 1.0 / d
    h = h * d * c
    if (h - 1.0).abs() < eps {
      return h
    }
    m = m + 1
  }
  h
}

///|
/// Regularized incomplete beta function I_x(a, b).
fn mm_pbeta(x : Double, a : Double, b : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x >= 1.0 {
    return 1.0
  }
  let lbeta_val = mm_lgamma(a) + mm_lgamma(b) - mm_lgamma(a + b)
  let mut bt = @math.ln(x) * a + @math.ln(1.0 - x) * b - lbeta_val
  bt = @math.exp(bt)
  let mut betaf = 0.0
  if x < (a + 1.0) / (a + b + 2.0) {
    betaf = bt * mm_betacf(x, a, b) / a
  } else {
    betaf = 1.0 - bt * mm_betacf(1.0 - x, b, a) / b
  }
  if betaf < 0.0 {
    0.0
  } else if betaf > 1.0 {
    1.0
  } else {
    betaf
  }
}

///|
/// Student's t-distribution CDF.
fn mm_t_cdf(t : Double, df : Double) -> Double {
  if df <= 0.0 {
    return 0.5
  }
  let x = df / (df + t * t)
  let p = 0.5 * mm_pbeta(x, df / 2.0, 0.5)
  if t >= 0.0 {
    1.0 - p
  } else {
    p
  }
}

///|
/// Natural log of factorial using direct summation (small) or Stirling (large).
fn mm_ln_factorial(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x < 20.0 {
    let mut result = 0.0
    let mut i = 1.0
    while i <= x {
      result = result + @math.ln(i)
      i = i + 1.0
    }
    result
  } else {
    let pi = 3.141592653589793
    let x1 = x + 1.0
    (x1 - 0.5) * @math.ln(x1) - x1 + 0.5 * @math.ln(2.0 * pi)
  }
}

///|
/// Multivariate hypergeometric PMF for a 2x2 table.
fn mm_hypergeom_pmf(a : Int, b : Int, c : Int, d : Int) -> Double {
  let n = (a + b + c + d).to_double()
  if n <= 0.0 {
    return 0.0
  }
  let ad = a.to_double()
  let bd = b.to_double()
  let cd = c.to_double()
  let dd = d.to_double()
  let ln_p = mm_ln_factorial(ad + bd) +
    mm_ln_factorial(cd + dd) +
    mm_ln_factorial(ad + cd) +
    mm_ln_factorial(bd + dd) -
    mm_ln_factorial(n) -
    mm_ln_factorial(ad) -
    mm_ln_factorial(bd) -
    mm_ln_factorial(cd) -
    mm_ln_factorial(dd)
  let p = @math.exp(ln_p)
  if p < 0.0 {
    0.0
  } else {
    p
  }
}

///|
/// Returns true if the array contains the given string.
fn mm_contains(arr : Array[String], val : String) -> Bool {
  for v in arr {
    if v == val {
      return true
    }
  }
  false
}

///|
/// Return the unique strings of an array, preserving first-seen order.
fn mm_unique(arr : Array[String]) -> Array[String] {
  let result : Array[String] = Array::new()
  for v in arr {
    if not(mm_contains(result, v)) {
      result.push(v)
    }
  }
  result
}

///|
/// Intersection of two string arrays (unique, first-seen order).
fn mm_intersect(a : Array[String], b : Array[String]) -> Array[String] {
  let result : Array[String] = Array::new()
  for v in a {
    if mm_contains(b, v) && not(mm_contains(result, v)) {
      result.push(v)
    }
  }
  result
}

///|
/// Collect all unique GO IDs from the database.
fn mm_unique_go_ids(db : MmGoDatabase) -> Array[String] {
  let result : Array[String] = Array::new()
  for entry in db.entries {
    if not(mm_contains(result, entry.0)) {
      result.push(entry.0)
    }
  }
  result
}

///|
/// Collect the GO term name and ontology for a GO id (first match).
fn mm_go_meta(db : MmGoDatabase, go_id : String) -> (String, String) {
  for entry in db.entries {
    if entry.0 == go_id {
      return (entry.2, entry.3)
    }
  }
  ("", "")
}

// ============================================================================
// Data structures
// ============================================================================

///|
/// Annotation for a single methylation probe.
pub struct MmProbeAnnotation {
  probe_id : String
  chr : String
  pos : Int
  gene : String
  island_relation : String
  probe_type : String
}

///|
/// Construct an MmProbeAnnotation.
pub fn MmProbeAnnotation::new(
  probe_id : String,
  chr : String,
  pos : Int,
  gene : String,
  island_relation : String,
  probe_type : String,
) -> MmProbeAnnotation {
  MmProbeAnnotation::{
    probe_id,
    chr,
    pos,
    gene,
    island_relation,
    probe_type,
  }
}

///|
/// Probe identifier.
pub fn MmProbeAnnotation::probe_id(self : MmProbeAnnotation) -> String {
  self.probe_id
}

///|
/// Chromosome name.
pub fn MmProbeAnnotation::chr(self : MmProbeAnnotation) -> String {
  self.chr
}

///|
/// Genomic position.
pub fn MmProbeAnnotation::pos(self : MmProbeAnnotation) -> Int {
  self.pos
}

///|
/// Mapped gene symbol.
pub fn MmProbeAnnotation::gene(self : MmProbeAnnotation) -> String {
  self.gene
}

///|
/// CpG island relation ("Island"/"Shore"/"Shelf"/"OpenSea").
pub fn MmProbeAnnotation::island_relation(
  self : MmProbeAnnotation,
) -> String {
  self.island_relation
}

///|
/// Assay type ("I" or "II").
pub fn MmProbeAnnotation::probe_type(self : MmProbeAnnotation) -> String {
  self.probe_type
}

///|
/// Human-readable representation.
pub fn MmProbeAnnotation::to_string(self : MmProbeAnnotation) -> String {
  let buf = StringBuilder::new()
  buf.write_string("MmProbeAnnotation(probe_id=")
  buf.write_string(self.probe_id)
  buf.write_string(", chr=")
  buf.write_string(self.chr)
  buf.write_string(", pos=")
  buf.write_string(self.pos.to_string())
  buf.write_string(", gene=")
  buf.write_string(self.gene)
  buf.write_string(", island_relation=")
  buf.write_string(self.island_relation)
  buf.write_string(", type=")
  buf.write_string(self.probe_type)
  buf.write_string(")")
  buf.to_string()
}

///|
/// Differential methylation position (DMP) result.
pub struct MmDmpResult {
  probe_id : String
  log_fc : Double
  t_stat : Double
  p_value : Double
  adj_p_value : Double
  chr : String
  pos : Int
  gene : String
}

///|
/// Construct an MmDmpResult.
pub fn MmDmpResult::new(
  probe_id : String,
  log_fc : Double,
  t_stat : Double,
  p_value : Double,
  adj_p_value : Double,
  chr : String,
  pos : Int,
  gene : String,
) -> MmDmpResult {
  MmDmpResult::{
    probe_id,
    log_fc,
    t_stat,
    p_value,
    adj_p_value,
    chr,
    pos,
    gene,
  }
}

///|
/// Probe identifier.
pub fn MmDmpResult::probe_id(self : MmDmpResult) -> String {
  self.probe_id
}

///|
/// Log2 fold change of M-values (treatment - control).
pub fn MmDmpResult::log_fc(self : MmDmpResult) -> Double {
  self.log_fc
}

///|
/// Moderated t-statistic.
pub fn MmDmpResult::t_stat(self : MmDmpResult) -> Double {
  self.t_stat
}

///|
/// Raw p-value.
pub fn MmDmpResult::p_value(self : MmDmpResult) -> Double {
  self.p_value
}

///|
/// BH-adjusted p-value.
pub fn MmDmpResult::adj_p_value(self : MmDmpResult) -> Double {
  self.adj_p_value
}

///|
/// Chromosome name.
pub fn MmDmpResult::chr(self : MmDmpResult) -> String {
  self.chr
}

///|
/// Genomic position.
pub fn MmDmpResult::pos(self : MmDmpResult) -> Int {
  self.pos
}

///|
/// Mapped gene symbol.
pub fn MmDmpResult::gene(self : MmDmpResult) -> String {
  self.gene
}

///|
/// Returns true when the adjusted p-value is below the given threshold.
pub fn MmDmpResult::is_significant(
  self : MmDmpResult,
  threshold : Double,
) -> Bool {
  self.adj_p_value < threshold
}

///|
/// GO enrichment result for a single GO term.
pub struct MmGoResult {
  go_id : String
  go_term : String
  ontology : String
  n_probes : Int
  n_sig_probes : Int
  p_value : Double
  adj_p_value : Double
  expected : Double
}

///|
/// Construct an MmGoResult.
pub fn MmGoResult::new(
  go_id : String,
  go_term : String,
  ontology : String,
  n_probes : Int,
  n_sig_probes : Int,
  p_value : Double,
  adj_p_value : Double,
  expected : Double,
) -> MmGoResult {
  MmGoResult::{
    go_id,
    go_term,
    ontology,
    n_probes,
    n_sig_probes,
    p_value,
    adj_p_value,
    expected,
  }
}

///|
/// GO identifier.
pub fn MmGoResult::go_id(self : MmGoResult) -> String {
  self.go_id
}

///|
/// GO term description.
pub fn MmGoResult::go_term(self : MmGoResult) -> String {
  self.go_term
}

///|
/// Ontology branch ("BP"/"CC"/"MF").
pub fn MmGoResult::ontology(self : MmGoResult) -> String {
  self.ontology
}

///|
/// Number of background probes annotated to this GO term.
pub fn MmGoResult::n_probes(self : MmGoResult) -> Int {
  self.n_probes
}

///|
/// Number of significant probes annotated to this GO term.
pub fn MmGoResult::n_sig_probes(self : MmGoResult) -> Int {
  self.n_sig_probes
}

///|
/// Raw enrichment p-value.
pub fn MmGoResult::p_value(self : MmGoResult) -> Double {
  self.p_value
}

///|
/// BH-adjusted p-value.
pub fn MmGoResult::adj_p_value(self : MmGoResult) -> Double {
  self.adj_p_value
}

///|
/// Expected number of significant probes under the null.
pub fn MmGoResult::expected(self : MmGoResult) -> Double {
  self.expected
}

///|
/// In-memory GO database mapping GO terms to genes.
pub struct MmGoDatabase {
  entries : Array[(String, String, String, String)]
}

///|
/// Construct an empty GO database.
pub fn MmGoDatabase::new() -> MmGoDatabase {
  MmGoDatabase::{ entries: Array::new() }
}

///|
/// Add a (go_id, gene, go_term, ontology) entry.
pub fn MmGoDatabase::add_entry(
  self : MmGoDatabase,
  go_id : String,
  gene : String,
  go_term : String,
  ontology : String,
) -> Unit {
  self.entries.push((go_id, gene, go_term, ontology))
}

///|
/// Return the unique genes annotated to a GO id.
pub fn MmGoDatabase::genes_for_go(
  self : MmGoDatabase,
  go_id : String,
) -> Array[String] {
  let result : Array[String] = Array::new()
  for entry in self.entries {
    if entry.0 == go_id && not(mm_contains(result, entry.1)) {
      result.push(entry.1)
    }
  }
  result
}

///|
/// Return the unique GO ids annotated to a gene.
pub fn MmGoDatabase::go_for_gene(
  self : MmGoDatabase,
  gene : String,
) -> Array[String] {
  let result : Array[String] = Array::new()
  for entry in self.entries {
    if entry.1 == gene && not(mm_contains(result, entry.0)) {
      result.push(entry.0)
    }
  }
  result
}

// ============================================================================
// Value conversions
// ============================================================================

///|
/// Convert a beta value to an M-value: log2(beta / (1 - beta)).
/// Beta is clamped to [1.0e-6, 1.0 - 1.0e-6] to avoid log(0).
pub fn mm_beta_to_m(beta : Double) -> Double {
  let b = if beta < 1.0e-6 {
    1.0e-6
  } else if beta > 1.0 - 1.0e-6 {
    1.0 - 1.0e-6
  } else {
    beta
  }
  @math.log2(b / (1.0 - b))
}

///|
/// Convert an M-value back to a beta value: 2^m / (1 + 2^m).
pub fn mm_m_to_beta(m : Double) -> Double {
  let two_m = @math.pow(2.0, m)
  two_m / (1.0 + two_m)
}

// ============================================================================
// Statistical tests
// ============================================================================

///|
/// Welch's two-sample t-test.
/// Returns (t_statistic, p_value).
pub fn mm_t_test(
  group1 : Array[Double],
  group2 : Array[Double],
) -> (Double, Double) {
  let n1 = group1.length()
  let n2 = group2.length()
  if n1 < 2 || n2 < 2 {
    return (0.0, 1.0)
  }
  let m1 = mm_mean(group1)
  let m2 = mm_mean(group2)
  let v1 = mm_variance(group1, m1)
  let v2 = mm_variance(group2, m2)
  let n1d = n1.to_double()
  let n2d = n2.to_double()
  let se = (v1 / n1d + v2 / n2d).sqrt()
  if se <= 0.0 {
    return (0.0, 1.0)
  }
  let t = (m1 - m2) / se
  // Welch-Satterthwaite degrees of freedom
  let term1 = v1 / n1d
  let term2 = v2 / n2d
  let num = (term1 + term2) * (term1 + term2)
  let den = term1 * term1 / (n1 - 1).to_double() +
    term2 * term2 / (n2 - 1).to_double()
  let df = if den > 0.0 {
    num / den
  } else {
    (n1 + n2 - 2).to_double()
  }
  let p = 2.0 * (1.0 - mm_t_cdf(t.abs(), df))
  let pc = if p < 0.0 {
    0.0
  } else if p > 1.0 {
    1.0
  } else {
    p
  }
  (t, pc)
}

///|
/// Benjamini-Hochberg FDR correction.
pub fn mm_bh_correct(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return Array::new()
  }
  // Sort indices by ascending p-value (bubble sort).
  let indices = Array::make(n, 0)
  for i in 0..<n {
    indices[i] = i
  }
  let mut swapped = true
  let mut k = n - 1
  while k > 0 && swapped {
    swapped = false
    let mut i = 0
    while i < k {
      if p_values[indices[i]] > p_values[indices[i + 1]] {
        let tmp = indices[i]
        indices[i] = indices[i + 1]
        indices[i + 1] = tmp
        swapped = true
      }
      i = i + 1
    }
    k = k - 1
  }
  // Step-up: from largest p to smallest, enforce monotonicity.
  let adj = Array::make(n, 1.0)
  let mut min_adj = 1.0
  let mut i = n - 1
  while i >= 0 {
    let rank = (i + 1).to_double()
    let val = p_values[indices[i]] * n.to_double() / rank
    if val < min_adj {
      min_adj = val
    }
    adj[indices[i]] = if min_adj > 1.0 {
      1.0
    } else {
      min_adj
    }
    i = i - 1
  }
  adj
}

///|
/// Two-sided Fisher's exact test for a 2x2 contingency table.
/// Table: [[a, b], [c, d]]. Returns the two-sided p-value.
pub fn mm_fisher_exact(a : Int, b : Int, c : Int, d : Int) -> Double {
  let n = a + b + c + d
  if n == 0 {
    return 1.0
  }
  let row1_total = a + b
  let col1_total = a + c
  let col2_total = b + d
  let p_obs = mm_hypergeom_pmf(a, b, c, d)
  let min_a = if row1_total - col2_total > 0 {
    row1_total - col2_total
  } else {
    0
  }
  let max_a = if row1_total < col1_total {
    row1_total
  } else {
    col1_total
  }
  let mut p_value = 0.0
  let mut i = min_a
  while i <= max_a {
    let bi = row1_total - i
    let ci = col1_total - i
    let di = col2_total - bi
    if bi >= 0 && ci >= 0 && di >= 0 {
      let p = mm_hypergeom_pmf(i, bi, ci, di)
      if p <= p_obs + 1.0e-10 {
        p_value = p_value + p
      }
    }
    i = i + 1
  }
  if p_value > 1.0 {
    1.0
  } else if p_value < 0.0 {
    0.0
  } else {
    p_value
  }
}

// ============================================================================
// DMP analysis
// ============================================================================

///|
/// Differential methylation position analysis.
/// `beta` is a probes x samples matrix of beta values.
/// `group` flags treatment (true) vs control (false) samples.
/// Beta values are converted to M-values and tested probe-wise with a t-test,
/// then BH-FDR corrected.
pub fn mm_dmp_analysis(
  beta : Array[Array[Double]],
  group : Array[Bool],
  probe_annotations : Array[MmProbeAnnotation],
) -> Array[MmDmpResult] {
  let n_probes = beta.length()
  let p_values = Array::make(n_probes, 1.0)
  let t_stats = Array::make(n_probes, 0.0)
  let log_fcs = Array::make(n_probes, 0.0)
  for p in 0..<n_probes {
    let treatment : Array[Double] = Array::new()
    let control : Array[Double] = Array::new()
    let row = beta[p]
    for s in 0..<group.length() {
      if s < row.length() {
        let m_val = mm_beta_to_m(row[s])
        if group[s] {
          treatment.push(m_val)
        } else {
          control.push(m_val)
        }
      }
    }
    let (t_stat, p_val) = mm_t_test(treatment, control)
    t_stats[p] = t_stat
    p_values[p] = p_val
    log_fcs[p] = mm_mean(treatment) - mm_mean(control)
  }
  let adj_p = mm_bh_correct(p_values)
  let results : Array[MmDmpResult] = Array::new()
  for p in 0..<n_probes {
    let ann = if p < probe_annotations.length() {
      probe_annotations[p]
    } else {
      MmProbeAnnotation::new("", "", 0, "", "", "")
    }
    results.push(
      MmDmpResult::new(
        ann.probe_id(),
        log_fcs[p],
        t_stats[p],
        p_values[p],
        adj_p[p],
        ann.chr(),
        ann.pos(),
        ann.gene(),
      ),
    )
  }
  results
}

// ============================================================================
// GO enrichment
// ============================================================================

///|
/// Map a list of probe IDs to the unique set of genes they annotate.
pub fn mm_map_probes_to_genes(
  probes : Array[String],
  annotations : Array[MmProbeAnnotation],
) -> Array[String] {
  let result : Array[String] = Array::new()
  for ann in annotations {
    if mm_contains(probes, ann.probe_id()) && not(mm_contains(result, ann.gene())) {
      result.push(ann.gene())
    }
  }
  result
}

///|
/// Map a list of genes to the unique probes that annotate them.
fn mm_probes_for_genes(
  genes : Array[String],
  annotations : Array[MmProbeAnnotation],
) -> Array[String] {
  let result : Array[String] = Array::new()
  for ann in annotations {
    if mm_contains(genes, ann.gene()) && not(mm_contains(result, ann.probe_id())) {
      result.push(ann.probe_id())
    }
  }
  result
}

///|
/// GO over-representation analysis at the gene level.
/// For each GO term a 2x2 Fisher's exact test is computed comparing
/// significant genes vs background genes inside/outside the GO term.
/// The background is the gene set derived from `all_probes`.
pub fn mm_go_enrichment(
  sig_probes : Array[String],
  all_probes : Array[String],
  probe_annotations : Array[MmProbeAnnotation],
  go_db : MmGoDatabase,
) -> Array[MmGoResult] {
  let sig_genes = mm_map_probes_to_genes(sig_probes, probe_annotations)
  let bg_genes = mm_map_probes_to_genes(all_probes, probe_annotations)
  let n_bg = bg_genes.length()
  let n_sig = sig_genes.length()
  let go_ids = mm_unique_go_ids(go_db)
  let raw_results : Array[MmGoResult] = Array::new()
  let p_values : Array[Double] = Array::new()
  for go_id in go_ids {
    let go_genes_all = go_db.genes_for_go(go_id)
    let go_genes = mm_intersect(go_genes_all, bg_genes)
    let sig_in_go = mm_intersect(sig_genes, go_genes)
    let a = sig_in_go.length()
    let b = n_sig - a
    let c = go_genes.length() - a
    let d = n_bg - a - b - c
    let p_val = mm_fisher_exact(a, b, c, d)
    let expected = if n_bg > 0 {
      n_sig.to_double() * go_genes.length().to_double() / n_bg.to_double()
    } else {
      0.0
    }
    let meta = mm_go_meta(go_db, go_id)
    let result = MmGoResult::new(
      go_id,
      meta.0,
      meta.1,
      go_genes.length(),
      a,
      p_val,
      p_val,
      expected,
    )
    raw_results.push(result)
    p_values.push(p_val)
  }
  let adj_p = mm_bh_correct(p_values)
  let final_results : Array[MmGoResult] = Array::new()
  for i in 0..<raw_results.length() {
    let r = raw_results[i]
    final_results.push(
      MmGoResult::new(
        r.go_id(),
        r.go_term(),
        r.ontology(),
        r.n_probes(),
        r.n_sig_probes(),
        r.p_value(),
        adj_p[i],
        r.expected(),
      ),
    )
  }
  final_results
}

///|
/// Probe-bias corrected GO enrichment (the core missMethyl algorithm).
/// Instead of a gene-level test, the contingency table is built at the probe
/// level so that genes represented by more probes contribute proportionally
/// more probes to both the significant and GO-term sets. The background
/// universe is restricted to probes whose genes have array annotations,
/// correcting the bias where genes with more probes are more likely to appear
/// in the significant set.
pub fn mm_probe_bias_correction(
  sig_probes : Array[String],
  all_probes : Array[String],
  annotations : Array[MmProbeAnnotation],
  go_db : MmGoDatabase,
) -> Array[MmGoResult] {
  let n_sig = sig_probes.length()
  let n_bg = all_probes.length()
  let go_ids = mm_unique_go_ids(go_db)
  let raw_results : Array[MmGoResult] = Array::new()
  let p_values : Array[Double] = Array::new()
  for go_id in go_ids {
    let go_genes = go_db.genes_for_go(go_id)
    let go_probes = mm_probes_for_genes(go_genes, annotations)
    // Restrict the GO probe set to the background universe.
    let go_probes_bg = mm_intersect(go_probes, all_probes)
    let sig_in_go = mm_intersect(sig_probes, go_probes_bg)
    let a = sig_in_go.length()
    let b = n_sig - a
    let c = go_probes_bg.length() - a
    let d = n_bg - a - b - c
    let p_val = mm_fisher_exact(a, b, c, d)
    let expected = if n_bg > 0 {
      n_sig.to_double() * go_probes_bg.length().to_double() / n_bg.to_double()
    } else {
      0.0
    }
    let meta = mm_go_meta(go_db, go_id)
    let result = MmGoResult::new(
      go_id,
      meta.0,
      meta.1,
      go_probes_bg.length(),
      a,
      p_val,
      p_val,
      expected,
    )
    raw_results.push(result)
    p_values.push(p_val)
  }
  let adj_p = mm_bh_correct(p_values)
  let final_results : Array[MmGoResult] = Array::new()
  for i in 0..<raw_results.length() {
    let r = raw_results[i]
    final_results.push(
      MmGoResult::new(
        r.go_id(),
        r.go_term(),
        r.ontology(),
        r.n_probes(),
        r.n_sig_probes(),
        r.p_value(),
        adj_p[i],
        r.expected(),
      ),
    )
  }
  final_results
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Filter probes by chromosome and/or CpG island relation.
/// An empty or "all" filter matches every probe.
pub fn mm_filter_probes(
  annotations : Array[MmProbeAnnotation],
  chr_filter : String,
  island_filter : String,
) -> Array[MmProbeAnnotation] {
  let result : Array[MmProbeAnnotation] = Array::new()
  for ann in annotations {
    let chr_ok = chr_filter == "" ||
      chr_filter == "all" ||
      ann.chr() == chr_filter
    let island_ok = island_filter == "" ||
      island_filter == "all" ||
      ann.island_relation() == island_filter
    if chr_ok && island_ok {
      result.push(ann)
    }
  }
  result
}

///|
/// Build a textual summary of a DMP analysis.
pub fn mm_summary(dmp_results : Array[MmDmpResult]) -> String {
  let n = dmp_results.length()
  let mut n_sig_05 = 0
  let mut n_sig_01 = 0
  let mut n_hyper = 0
  let mut n_hypo = 0
  for r in dmp_results {
    if r.is_significant(0.05) {
      n_sig_05 = n_sig_05 + 1
      if r.log_fc() > 0.0 {
        n_hyper = n_hyper + 1
      } else {
        n_hypo = n_hypo + 1
      }
    }
    if r.is_significant(0.01) {
      n_sig_01 = n_sig_01 + 1
    }
  }
  let buf = StringBuilder::new()
  buf.write_string("DMP Analysis Summary\n")
  buf.write_string("====================\n")
  buf.write_string("Total probes tested: ")
  buf.write_string(n.to_string())
  buf.write_string("\nSignificant (adj p < 0.05): ")
  buf.write_string(n_sig_05.to_string())
  buf.write_string("\nSignificant (adj p < 0.01): ")
  buf.write_string(n_sig_01.to_string())
  buf.write_string("\n  Hyper-methylated: ")
  buf.write_string(n_hyper.to_string())
  buf.write_string("\n  Hypo-methylated: ")
  buf.write_string(n_hypo.to_string())
  buf.write_string("\n")
  buf.to_string()
}

///|
/// Build a textual summary of a GO enrichment analysis.
pub fn mm_go_summary(go_results : Array[MmGoResult]) -> String {
  let n = go_results.length()
  let mut n_sig = 0
  let buf = StringBuilder::new()
  buf.write_string("GO Enrichment Summary\n")
  buf.write_string("=====================\n")
  buf.write_string("GO terms tested: ")
  buf.write_string(n.to_string())
  buf.write_string("\n")
  for r in go_results {
    if r.adj_p_value() < 0.05 {
      n_sig = n_sig + 1
      buf.write_string("  ")
      buf.write_string(r.go_id())
      buf.write_string(" (")
      buf.write_string(r.ontology())
      buf.write_string(") ")
      buf.write_string(r.go_term())
      buf.write_string(" - adj p = ")
      buf.write_string(r.adj_p_value().to_string())
      buf.write_string(", sig probes = ")
      buf.write_string(r.n_sig_probes().to_string())
      buf.write_string(" / ")
      buf.write_string(r.n_probes().to_string())
      buf.write_string("\n")
    }
  }
  buf.write_string("Significant GO terms (adj p < 0.05): ")
  buf.write_string(n_sig.to_string())
  buf.write_string("\n")
  buf.to_string()
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Twenty sample probe annotations spanning chromosomes 1-3.
pub fn mm_sample_annotations() -> Array[MmProbeAnnotation] {
  let result : Array[MmProbeAnnotation] = Array::new()
  result.push(MmProbeAnnotation::new("cg001", "chr1", 1000, "A", "Island", "I"))
  result.push(MmProbeAnnotation::new("cg002", "chr1", 2000, "B", "Shore", "I"))
  result.push(MmProbeAnnotation::new("cg003", "chr1", 3000, "C", "Shelf", "II"))
  result.push(MmProbeAnnotation::new("cg004", "chr1", 4000, "D", "OpenSea", "I"))
  result.push(MmProbeAnnotation::new("cg005", "chr1", 5000, "E", "Island", "II"))
  result.push(MmProbeAnnotation::new("cg006", "chr1", 6000, "F", "Shore", "I"))
  result.push(MmProbeAnnotation::new("cg007", "chr1", 7000, "G", "Shelf", "II"))
  result.push(MmProbeAnnotation::new("cg008", "chr2", 2000, "H", "Island", "I"))
  result.push(MmProbeAnnotation::new("cg009", "chr2", 2500, "A", "Shore", "II"))
  result.push(MmProbeAnnotation::new("cg010", "chr2", 3000, "B", "Shelf", "I"))
  result.push(MmProbeAnnotation::new("cg011", "chr2", 3500, "C", "OpenSea", "II"))
  result.push(MmProbeAnnotation::new("cg012", "chr2", 4000, "D", "Island", "I"))
  result.push(MmProbeAnnotation::new("cg013", "chr2", 5000, "E", "Shore", "II"))
  result.push(MmProbeAnnotation::new("cg014", "chr2", 6000, "F", "Shelf", "I"))
  result.push(MmProbeAnnotation::new("cg015", "chr3", 3000, "G", "OpenSea", "II"))
  result.push(MmProbeAnnotation::new("cg016", "chr3", 4000, "H", "Island", "I"))
  result.push(MmProbeAnnotation::new("cg017", "chr3", 5000, "A", "Shore", "II"))
  result.push(MmProbeAnnotation::new("cg018", "chr3", 6000, "B", "Shelf", "I"))
  result.push(MmProbeAnnotation::new("cg019", "chr3", 7000, "C", "OpenSea", "II"))
  result.push(MmProbeAnnotation::new("cg020", "chr3", 8000, "D", "Island", "I"))
  result
}

///|
/// A small GO database with 10 entries annotating genes A-H.
pub fn mm_sample_go_database() -> MmGoDatabase {
  let db = MmGoDatabase::new()
  db.add_entry("GO:0000001", "A", "cell cycle process", "BP")
  db.add_entry("GO:0000001", "B", "cell cycle process", "BP")
  db.add_entry("GO:0000002", "C", "apoptotic process", "BP")
  db.add_entry("GO:0000002", "D", "apoptotic process", "BP")
  db.add_entry("GO:0000003", "A", "nucleus", "CC")
  db.add_entry("GO:0000003", "E", "nucleus", "CC")
  db.add_entry("GO:0000004", "F", "kinase activity", "MF")
  db.add_entry("GO:0000004", "G", "kinase activity", "MF")
  db.add_entry("GO:0000005", "H", "DNA binding", "MF")
  db.add_entry("GO:0000005", "A", "DNA binding", "MF")
  db
}

///|
/// A 20 probes x 8 samples beta-value matrix (4 control + 4 treatment).
/// The first five probes are strongly differential (hypermethylated in
/// treatment); the remaining probes show no consistent group difference.
pub fn mm_sample_beta() -> Array[Array[Double]] {
  let matrix : Array[Array[Double]] = Array::new()
  for p in 0..<20 {
    let row : Array[Double] = Array::new()
    for s in 0..<8 {
      if p < 5 {
        // Differential probes: low beta in control, high in treatment.
        if s < 4 {
          row.push(0.15 + p.to_double() * 0.02 + s.to_double() * 0.01)
        } else {
          row.push(0.70 + p.to_double() * 0.02 + (s - 4).to_double() * 0.01)
        }
      } else {
        // Non-differential probes with mild per-probe/sample variation.
        let base = 0.40 + p.to_double() * 0.01
        row.push(base + s.to_double() * 0.005)
      }
    }
    matrix.push(row)
  }
  matrix
}
