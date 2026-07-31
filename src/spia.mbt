///|
/// SPIA module - Signaling Pathway Impact Analysis.
/// Inspired by Bioconductor SPIA package.
///
/// References:
///   - https://bioconductor.org/packages/release/bioc/html/SPIA.html
///   - Tarca AL et al. (2009) A novel signaling pathway impact analysis (SPIA)
///     for microarray data. Bioinformatics 25:75-82.

///|
/// Pathway graph node (gene).
pub struct PathwayNode {
  gene_id : String
  /// Gene name/symbol
  gene_name : String
  /// List of downstream gene IDs (targets)
  downstream : Array[String]
  /// List of upstream gene IDs (regulators)
  upstream : Array[String]
  /// Activation (1) or inhibition (-1) for each downstream interaction
  activation : Array[Int]
}

///|
/// A signaling pathway.
pub struct SignalingPathway {
  id : String
  name : String
  nodes : Array[PathwayNode]
}

///|
/// SPIA result for a single pathway.
pub struct SpiaResult {
  pathway_id : String
  pathway_name : String
  /// Number of genes in the pathway
  n_genes : Int
  /// Number of DE genes in the pathway
  n_de : Int
  /// Over-representation p-value (from hypergeometric test)
  p_or : Double
  /// Perturbation accumulation factor
  pert_factor : Double
  /// Perturbation p-value
  p_pert : Double
  /// Combined p-value (using Fisher's method)
  p_combined : Double
  /// Adjusted combined p-value (BH)
  mut p_adj : Double
  /// Is the pathway significant?
  mut significant : Bool
  /// Activation status: 1 = activated, -1 = inhibited, 0 = neutral
  activation_status : Int
}

///|
/// Container for SPIA results.
pub struct SpiaResults {
  results : Array[SpiaResult]
  n_pathways : Int
  n_significant : Int
  fdr_threshold : Double
}

///|
/// Create a PathwayNode.
pub fn PathwayNode::new(gene_id : String, gene_name : String) -> PathwayNode {
  PathwayNode::{
    gene_id,
    gene_name,
    downstream: Array::new(),
    upstream: Array::new(),
    activation: Array::new(),
  }
}

///|
/// Add a downstream interaction.
pub fn PathwayNode::add_downstream(
  self : PathwayNode,
  target_id : String,
  is_activation : Bool,
) -> Unit {
  self.downstream.push(target_id)
  if is_activation {
    self.activation.push(1)
  } else {
    self.activation.push(-1)
  }
}

///|
/// Create a SignalingPathway.
pub fn SignalingPathway::new(id : String, name : String) -> SignalingPathway {
  SignalingPathway::{ id, name, nodes: Array::new() }
}

///|
/// Add a node to the pathway.
pub fn SignalingPathway::add_node(self : SignalingPathway, node : PathwayNode) -> Unit {
  self.nodes.push(node)
}

///|
/// Get a node by gene ID.
pub fn SignalingPathway::get_node(self : SignalingPathway, gene_id : String) -> PathwayNode? {
  let mut i = 0
  while i < self.nodes.length() {
    if self.nodes[i].gene_id == gene_id {
      return Some(self.nodes[i])
    }
    i = i + 1
  }
  None
}

///|
/// Run SPIA on a set of pathways.
///
/// Parameters:
///   - pathways: Array of signaling pathways
///   - de_genes: Map from gene_id to log2 fold change (only DE genes)
///   - total_genes: Total number of genes measured
///   - fdr_threshold: FDR threshold for significance
pub fn run_spia(
  pathways : Array[SignalingPathway],
  de_genes : Map[String, Double],
  total_genes : Int,
  fdr_threshold? : Double = 0.05,
) -> SpiaResults {
  let n_total = if total_genes > 0 { total_genes } else { 1 }
  let n_de = de_genes.size()
  let all_results : Array[SpiaResult] = Array::new()
  let mut pi = 0
  while pi < pathways.length() {
    let pathway = pathways[pi]
    // Count DE genes in pathway
    let n_pathway_genes = pathway.nodes.length()
    let mut n_de_in_pathway = 0
    let mut ni = 0
    while ni < n_pathway_genes {
      if de_genes.contains(pathway.nodes[ni].gene_id) {
        n_de_in_pathway = n_de_in_pathway + 1
      }
      ni = ni + 1
    }
    // Over-representation p-value (hypergeometric test)
    let p_or = hypergeometric_pvalue(
      n_de_in_pathway,
      n_pathway_genes,
      n_de,
      n_total,
    )
    // Perturbation accumulation
    let pert_factor = compute_perturbation(pathway, de_genes)
    // Perturbation p-value (normal approximation)
    let p_pert = perturbation_pvalue(pert_factor, n_de_in_pathway)
    // Combine p-values using Fisher's method
    let p_combined = fisher_combine(p_or, p_pert)
    // Activation status
    let activation_status = if pert_factor > 0.0 { 1 } else if pert_factor < 0.0 { -1 } else { 0 }
    all_results.push(SpiaResult::{
      pathway_id: pathway.id,
      pathway_name: pathway.name,
      n_genes: n_pathway_genes,
      n_de: n_de_in_pathway,
      p_or,
      pert_factor,
      p_pert,
      p_combined,
      p_adj: p_combined,
      significant: false,
      activation_status,
    })
    pi = pi + 1
  }
  // Apply BH correction
  let p_values : Array[Double] = Array::new()
  let mut ri = 0
  while ri < all_results.length() {
    p_values.push(all_results[ri].p_combined)
    ri = ri + 1
  }
  let adjusted = spia_bh_adjust(p_values)
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
  SpiaResults::{
    results: all_results,
    n_pathways: all_results.length(),
    n_significant: n_sig,
    fdr_threshold,
  }
}

///|
/// Compute perturbation accumulation factor for a pathway.
/// This measures the net flow of perturbation through the pathway.
fn compute_perturbation(
  pathway : SignalingPathway,
  de_genes : Map[String, Double],
) -> Double {
  let n_nodes = pathway.nodes.length()
  if n_nodes == 0 {
    return 0.0
  }
  // Compute perturbation for each gene in the pathway
  // PF(g) = sum over upstream DE genes of (logFC * activation_sign) + logFC(g) if g is DE
  let perturbations : Map[String, Double] = Map([], capacity=64)
  let mut i = 0
  while i < n_nodes {
    let node = pathway.nodes[i]
    let mut pert = 0.0
    // Start with own logFC if DE
    match de_genes.get(node.gene_id) {
      Some(lfc) => pert = pert + lfc
      None => ()
    }
    // Add perturbations from upstream DE genes
    // We need to find all nodes that have this node as downstream
    let mut j = 0
    while j < n_nodes {
      let upstream_node = pathway.nodes[j]
      // Check if upstream_node has node as downstream
      let mut dk = 0
      while dk < upstream_node.downstream.length() {
        if upstream_node.downstream[dk] == node.gene_id {
          // Found an upstream interaction
          match de_genes.get(upstream_node.gene_id) {
            Some(lfc) => {
              let sign = upstream_node.activation[dk].to_double()
              pert = pert + lfc * sign
            }
            None => ()
          }
        }
        dk = dk + 1
      }
      j = j + 1
    }
    perturbations.set(node.gene_id, pert)
    i = i + 1
  }
  // Sum all perturbations and normalize
  let mut total_pert = 0.0
  let keys = perturbations.keys().collect()
  let mut ki = 0
  while ki < keys.length() {
    total_pert = total_pert + perturbations[keys[ki]]
    ki = ki + 1
  }
  // Normalize by number of DE genes (avoid division by zero)
  let n_de_in_pathway = {
    let mut cnt = 0
    let mut ci = 0
    while ci < n_nodes {
      if de_genes.contains(pathway.nodes[ci].gene_id) {
        cnt = cnt + 1
      }
      ci = ci + 1
    }
    cnt
  }
  if n_de_in_pathway > 0 {
    total_pert / n_de_in_pathway.to_double()
  } else {
    0.0
  }
}

///|
/// Compute hypergeometric p-value (one-sided, over-representation).
fn hypergeometric_pvalue(
  k : Int, /// DE genes in pathway
  m : Int, /// genes in pathway
  n : Int, /// total DE genes
  n_total : Int, /// total genes
) -> Double {
  // Using a simplified chi-squared approximation
  // Actually, let's use a simple Fisher's exact test approach
  if k == 0 || m == 0 || n == 0 {
    return 1.0
  }
  // Expected value
  let expected = n.to_double() * m.to_double() / n_total.to_double()
  if expected < 1.0e-10 {
    return 1.0
  }
  // Chi-squared statistic
  let observed = k.to_double()
  let chi_sq = (observed - expected) * (observed - expected) / expected
  // One-sided p-value from chi-squared (df=1)
  // p = exp(-chi_sq / 2) / sqrt(2 * pi * chi_sq) approx
  if chi_sq < 1.0e-10 {
    return 1.0
  }
  let p = @math.exp(-chi_sq / 2.0) / (2.0 * 3.14159265358979 * chi_sq).sqrt()
  if p > 1.0 { 1.0 } else if p < 0.0 { 0.0 } else { p }
}

///|
/// Compute perturbation p-value using normal approximation.
fn perturbation_pvalue(pert_factor : Double, n_de : Int) -> Double {
  if n_de == 0 {
    return 1.0
  }
  // The perturbation factor under null follows approximately N(0, sigma^2)
  // sigma ≈ sqrt(n_de) * mean_sigma
  let sigma = n_de.to_double().sqrt()
  if sigma < 1.0e-10 {
    return 1.0
  }
  let z = pert_factor / sigma
  let abs_z = if z < 0.0 { -z } else { z }
  // Two-sided p-value using normal CDF approximation
  let p = spia_erfc_approx(abs_z / 1.4142135623730951)
  if p > 1.0 { 1.0 } else if p < 0.0 { 0.0 } else { p }
}

///|
/// Combine two p-values using Fisher's method.
fn fisher_combine(p1 : Double, p2 : Double) -> Double {
  if p1 <= 0.0 || p2 <= 0.0 {
    return 0.0
  }
  if p1 >= 1.0 && p2 >= 1.0 {
    return 1.0
  }
  let stat = -2.0 * (@math.ln(p1) + @math.ln(p2))
  if stat <= 0.0 {
    return 1.0
  }
  // p-value from chi-squared with df=4 (since 2 p-values)
  // Use the formula: p = exp(-stat/2) * (1 + stat/2)
  let p = @math.exp(-stat / 2.0) * (1.0 + stat / 2.0)
  if p > 1.0 { 1.0 } else if p < 0.0 { 0.0 } else { p }
}

///|
/// Approximate complementary error function erfc(x) for SPIA module.
fn spia_erfc_approx(x : Double) -> Double {
  if x < 0.0 {
    return 2.0 - spia_erfc_approx(-x)
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
fn spia_bh_adjust(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return Array::new()
  }
  let indexed : Array[(Int, Double)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((i, p_values[i]))
    i = i + 1
  }
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
/// Get significant pathways.
pub fn SpiaResults::get_significant(self : SpiaResults) -> Array[SpiaResult] {
  let sig : Array[SpiaResult] = Array::new()
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
/// Get activated pathways.
pub fn SpiaResults::get_activated(self : SpiaResults) -> Array[SpiaResult] {
  let act : Array[SpiaResult] = Array::new()
  let mut i = 0
  while i < self.results.length() {
    if self.results[i].activation_status > 0 {
      act.push(self.results[i])
    }
    i = i + 1
  }
  act
}

///|
/// Get inhibited pathways.
pub fn SpiaResults::get_inhibited(self : SpiaResults) -> Array[SpiaResult] {
  let inh : Array[SpiaResult] = Array::new()
  let mut i = 0
  while i < self.results.length() {
    if self.results[i].activation_status < 0 {
      inh.push(self.results[i])
    }
    i = i + 1
  }
  inh
}

///|
/// Get top N pathways by combined p-value.
pub fn SpiaResults::get_top_pathways(self : SpiaResults, n : Int) -> Array[SpiaResult] {
  let sorted = self.results.copy()
  let len = sorted.length()
  let mut i = 0
  while i < len - 1 {
    let mut j = 0
    while j < len - 1 - i {
      if sorted[j].p_combined > sorted[j + 1].p_combined {
        let tmp = sorted[j]
        sorted.set(j, sorted[j + 1])
        sorted.set(j + 1, tmp)
      }
      j = j + 1
    }
    i = i + 1
  }
  let result : Array[SpiaResult] = Array::new()
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
pub fn SpiaResults::summary(self : SpiaResults) -> String {
  let mut s = "=== SPIA Analysis Summary ===\n".to_string()
  s = s + "Pathways tested: " + self.n_pathways.to_string() + "\n"
  s = s + "Significant (FDR < " + self.fdr_threshold.to_string() + "): " + self.n_significant.to_string() + "\n"
  let mut i = 0
  while i < self.results.length() {
    let r = self.results[i]
    if r.significant {
      let status = if r.activation_status > 0 { "Activated" } else if r.activation_status < 0 { "Inhibited" } else { "Neutral" }
      s = s + "  " + r.pathway_name + " (" + status + ", p_adj=" + r.p_adj.to_string() + ")\n"
    }
    i = i + 1
  }
  s
}

///|
/// Create sample pathways for testing.
pub fn spia_sample_pathways() -> Array[SignalingPathway] {
  let pathways : Array[SignalingPathway] = Array::new()
  // Pathway 1: Simple activation cascade
  let p1 = SignalingPathway::new("PW1", "MAPK Cascade")
  let n1 = PathwayNode::new("Gene1", "EGFR")
  n1.add_downstream("Gene2", true) // activates
  p1.add_node(n1)
  let n2 = PathwayNode::new("Gene2", "RAS")
  n2.add_downstream("Gene3", true)
  p1.add_node(n2)
  let n3 = PathwayNode::new("Gene3", "RAF")
  n3.add_downstream("Gene4", true)
  p1.add_node(n3)
  let n4 = PathwayNode::new("Gene4", "MEK")
  n4.add_downstream("Gene5", true)
  p1.add_node(n4)
  let n5 = PathwayNode::new("Gene5", "ERK")
  p1.add_node(n5)
  pathways.push(p1)
  // Pathway 2: Apoptosis pathway with inhibition
  let p2 = SignalingPathway::new("PW2", "Apoptosis")
  let a1 = PathwayNode::new("Gene6", "BAX")
  a1.add_downstream("Gene8", true)
  p2.add_node(a1)
  let a2 = PathwayNode::new("Gene7", "BCL2")
  a2.add_downstream("Gene8", false) // inhibits
  p2.add_node(a2)
  let a3 = PathwayNode::new("Gene8", "CASP3")
  p2.add_node(a3)
  pathways.push(p2)
  pathways
}

///|
/// Create sample DE genes for testing.
pub fn spia_sample_de_genes() -> Map[String, Double] {
  let de : Map[String, Double] = Map::new()
  de.set("Gene1", 2.5) // upregulated
  de.set("Gene2", 1.8)
  de.set("Gene5", -1.5) // downregulated
  de.set("Gene6", 2.0)
  de
}
