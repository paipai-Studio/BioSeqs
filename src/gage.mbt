///|
/// GAGE module - Generally Applicable Gene-set Enrichment analysis.
/// Inspired by Bioconductor GAGE package.
///
/// References:
///   - https://bioconductor.org/packages/release/bioc/html/gage.html
///   - Luo W et al. (2009) GAGE: generally applicable gene set enrichment
///     for pathway analysis. BMC Bioinformatics 10:161.

///|
/// Gene set test result.
pub struct GageResult {
  gene_set_id : String
  gene_set_name : String
  /// Number of genes in the set
  n_genes : Int
  /// Mean fold change of the gene set
  mean_fc : Double
  /// Test statistic (t-value or z-score)
  stat : Double
  /// Raw p-value
  p_val : Double
  /// Adjusted p-value (BH)
  mut p_adj : Double
  /// Is the set up-regulated?
  is_up : Bool
  /// Is the set down-regulated?
  is_down : Bool
  /// Is the set significant?
  mut significant : Bool
}

///|
/// Container for GAGE results.
pub struct GageResults {
  results : Array[GageResult]
  n_gene_sets : Int
  n_significant : Int
  fdr_threshold : Double
  /// Whether paired test was used
  paired : Bool
}

///|
/// A gene set definition for GAGE analysis.
pub struct GageGeneSet {
  id : String
  name : String
  gene_ids : Array[String]
}

///|
/// Create a GageGeneSet.
pub fn GageGeneSet::new(id : String, name : String, gene_ids : Array[String]) -> GageGeneSet {
  GageGeneSet::{ id, name, gene_ids }
}

///|
/// Run GAGE gene set enrichment analysis.
/// Uses fold-change-based t-test for each gene set.
///
/// Parameters:
///   - gene_fcs: Map from gene_id to fold change (log2)
///   - gene_sets: List of gene sets to test
///   - fdr_threshold: FDR threshold for significance
///   - paired: Whether to use paired test (changes variance calculation)
pub fn run_gage(
  gene_fcs : Map[String, Double],
  gene_sets : Array[GageGeneSet],
  fdr_threshold? : Double = 0.05,
  paired? : Bool = false,
) -> GageResults {
  let all_results : Array[GageResult] = Array::new()
  let mut gsi = 0
  while gsi < gene_sets.length() {
    let gs = gene_sets[gsi]
    // Collect fold changes for genes in this set
    let fcs : Array[Double] = Array::new()
    let mut gi = 0
    while gi < gs.gene_ids.length() {
      let gid = gs.gene_ids[gi]
      match gene_fcs.get(gid) {
        Some(fc) => fcs.push(fc)
        None => ()
      }
      gi = gi + 1
    }
    if fcs.length() > 0 {
      let mean_fc = gage_mean(fcs)
      // Compute test statistic
      let se = gage_se(fcs, paired)
      let stat = if se > 1.0e-10 { mean_fc / se } else { 0.0 }
      // Compute p-value using normal approximation
      let p_val = gage_pvalue_from_z(stat)
      let is_up = mean_fc > 0.0
      let is_down = mean_fc < 0.0
      all_results.push(GageResult::{
        gene_set_id: gs.id,
        gene_set_name: gs.name,
        n_genes: fcs.length(),
        mean_fc,
        stat,
        p_val,
        p_adj: p_val,
        is_up,
        is_down,
        significant: false,
      })
    }
    gsi = gsi + 1
  }
  // Apply BH correction
  let p_values : Array[Double] = Array::new()
  let mut ri = 0
  while ri < all_results.length() {
    p_values.push(all_results[ri].p_val)
    ri = ri + 1
  }
  let adjusted = gage_bh_adjust(p_values)
  let mut n_sig = 0
  let mut si = 0
  while si < all_results.length() {
    all_results[si].p_adj = adjusted[si]
    all_results[si].significant = adjusted[si] < fdr_threshold
    if all_results[si].significant {
      n_sig = n_sig + 1
    }
    si = si + 1
  }
  GageResults::{
    results: all_results,
    n_gene_sets: all_results.length(),
    n_significant: n_sig,
    fdr_threshold,
    paired,
  }
}

///|
/// Compute mean of an array.
fn gage_mean(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < values.length() {
    sum = sum + values[i]
    i = i + 1
  }
  sum / values.length().to_double()
}

///|
/// Compute standard error for GAGE test.
fn gage_se(values : Array[Double], paired : Bool) -> Double {
  let n = values.length()
  if n < 2 {
    return 1.0
  }
  let mean = gage_mean(values)
  let mut sum_sq = 0.0
  let mut i = 0
  while i < n {
    let diff = values[i] - mean
    sum_sq = sum_sq + diff * diff
    i = i + 1
  }
  let variance = sum_sq / (n - 1).to_double()
  // For paired test, reduce variance by factor of 2 (correlated samples)
  let adj_variance = if paired { variance / 2.0 } else { variance }
  (adj_variance / n.to_double()).sqrt()
}

///|
/// Compute two-sided p-value from z-score using normal approximation.
fn gage_pvalue_from_z(z : Double) -> Double {
  let abs_z = if z < 0.0 { -z } else { z }
  // Use complementary error function approximation
  // p = 2 * (1 - Phi(|z|)) = erfc(|z| / sqrt(2))
  let x = abs_z / 1.4142135623730951
  let p = erfc_approx(x)
  if p > 1.0 { 1.0 } else if p < 0.0 { 0.0 } else { p }
}

///|
/// Approximate complementary error function erfc(x).
/// Uses Abramowitz and Stegun approximation 7.1.26.
fn erfc_approx(x : Double) -> Double {
  if x < 0.0 {
    return 2.0 - erfc_approx(-x)
  }
  let t = 1.0 / (1.0 + 0.3275911 * x)
  let t2 = t * t
  let t3 = t2 * t
  let t4 = t3 * t
  let t5 = t4 * t
  let poly = 1.0 + 0.254829592 * t + -0.284496736 * t2 + 1.421413741 * t3 + -1.453152027 * t4 + 1.061405429 * t5
  @math.exp(-x * x) * poly
}

///|
/// Benjamini-Hochberg FDR adjustment.
fn gage_bh_adjust(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return Array::new()
  }
  // Create index-value pairs and sort by p-value
  let indexed : Array[(Int, Double)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((i, p_values[i]))
    i = i + 1
  }
  // Sort by p-value (bubble sort)
  let mut j = 0
  while j < n - 1 {
    let mut k = 0
    while k < n - 1 - j {
      if indexed[k].1 > indexed[k + 1].1 {
        let tmp = indexed[k]
        indexed.set(k, indexed[k + 1])
        indexed.set(k + 1, tmp)
      }
      k = k + 1
    }
    j = j + 1
  }
  // Compute adjusted p-values
  let adjusted : Array[Double] = Array::new()
  let mut ai = 0
  while ai < n {
    adjusted.push(0.0)
    ai = ai + 1
  }
  let mut prev_adj = 1.0
  let mut ri = n - 1
  while ri >= 0 {
    let rank = ri + 1
    let raw_p = indexed[ri].1
    let bh_val = raw_p * n.to_double() / rank.to_double()
    if bh_val < prev_adj {
      prev_adj = bh_val
    }
    if prev_adj > 1.0 {
      adjusted.set(indexed[ri].0, 1.0)
    } else {
      adjusted.set(indexed[ri].0, prev_adj)
    }
    ri = ri - 1
  }
  adjusted
}

///|
/// Get significant results only.
pub fn GageResults::get_significant(self : GageResults) -> Array[GageResult] {
  let sig : Array[GageResult] = Array::new()
  let mut i = 0
  while i < self.results.length() {
    if self.results[i].significant {
      sig.push(self.results[i])
    }
    i = i + 1
  }
  sig
}

///|
/// Get up-regulated gene sets.
pub fn GageResults::get_upregulated(self : GageResults) -> Array[GageResult] {
  let up : Array[GageResult] = Array::new()
  let mut i = 0
  while i < self.results.length() {
    if self.results[i].is_up {
      up.push(self.results[i])
    }
    i = i + 1
  }
  up
}

///|
/// Get down-regulated gene sets.
pub fn GageResults::get_downregulated(self : GageResults) -> Array[GageResult] {
  let down : Array[GageResult] = Array::new()
  let mut i = 0
  while i < self.results.length() {
    if self.results[i].is_down {
      down.push(self.results[i])
    }
    i = i + 1
  }
  down
}

///|
/// Get top N gene sets by adjusted p-value.
pub fn GageResults::get_top_sets(self : GageResults, n : Int) -> Array[GageResult] {
  let sorted = self.results.copy()
  let len = sorted.length()
  let mut i = 0
  while i < len - 1 {
    let mut j = 0
    while j < len - 1 - i {
      if sorted[j].p_adj > sorted[j + 1].p_adj {
        let tmp = sorted[j]
        sorted.set(j, sorted[j + 1])
        sorted.set(j + 1, tmp)
      }
      j = j + 1
    }
    i = i + 1
  }
  let result : Array[GageResult] = Array::new()
  let take = if n < len { n } else { len }
  let mut k = 0
  while k < take {
    result.push(sorted[k])
    k = k + 1
  }
  result
}

///|
/// Generate a summary string.
pub fn GageResults::summary(self : GageResults) -> String {
  let mut s = "=== GAGE Analysis Summary ===\n".to_string()
  s = s + "Gene sets tested: " + self.n_gene_sets.to_string() + "\n"
  s = s + "Significant (FDR < " + self.fdr_threshold.to_string() + "): " + self.n_significant.to_string() + "\n"
  s = s + "Paired test: " + self.paired.to_string() + "\n"
  s
}

///|
/// Create sample gene sets for testing.
pub fn gage_sample_gene_sets() -> Array[GageGeneSet] {
  let sets : Array[GageGeneSet] = Array::new()
  sets.push(GageGeneSet::new("GS1", "Pathway A", ["Gene1", "Gene2", "Gene3"]))
  sets.push(GageGeneSet::new("GS2", "Pathway B", ["Gene4", "Gene5", "Gene6"]))
  sets.push(GageGeneSet::new("GS3", "Pathway C", ["Gene1", "Gene4", "Gene7"]))
  sets
}

///|
/// Create sample fold changes for testing.
pub fn gage_sample_fold_changes() -> Map[String, Double] {
  let fcs : Map[String, Double] = Map::new()
  fcs.set("Gene1", 2.0)
  fcs.set("Gene2", 1.5)
  fcs.set("Gene3", 1.8)
  fcs.set("Gene4", -1.2)
  fcs.set("Gene5", -0.8)
  fcs.set("Gene6", -1.5)
  fcs.set("Gene7", 0.5)
  fcs
}
