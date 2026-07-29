///|
/// Singscore - Rank-based single-sample gene set scoring.
/// Inspired by the singscore R package (Foroutan et al., 2020).
/// Computes up/down regulated gene set scores using median-centered ranks.

///|
/// A single sample expression profile.
pub struct SampleExpression {
  pub sample_id : String
  pub gene_names : Array[String]
  pub expression : Array[Double]
}

///|
/// Get sample ID.
pub fn SampleExpression::get_sample_id(self : SampleExpression) -> String {
  self.sample_id
}

///|
/// Get gene names.
pub fn SampleExpression::get_gene_names(self : SampleExpression) -> Array[String] {
  self.gene_names
}

///|
/// Get expression values.
pub fn SampleExpression::get_expression(self : SampleExpression) -> Array[Double] {
  self.expression
}

///|
/// Create a sample expression profile.
pub fn SampleExpression::new(
  sample_id : String,
  gene_names : Array[String],
  expression : Array[Double],
) -> SampleExpression {
  SampleExpression::{ sample_id, gene_names, expression }
}

///|
/// A gene set specification with up-regulated and down-regulated genes.
pub struct GeneSetSpec {
  pub name : String
  pub up_genes : Array[String]
  pub down_genes : Array[String]
}

///|
/// Get gene set name.
pub fn GeneSetSpec::get_name(self : GeneSetSpec) -> String {
  self.name
}

///|
/// Get up-regulated genes.
pub fn GeneSetSpec::get_up_genes(self : GeneSetSpec) -> Array[String] {
  self.up_genes
}

///|
/// Get down-regulated genes.
pub fn GeneSetSpec::get_down_genes(self : GeneSetSpec) -> Array[String] {
  self.down_genes
}

///|
/// Create a gene set specification.
pub fn GeneSetSpec::new(
  name : String,
  up_genes : Array[String],
  down_genes? : Array[String] = [],
) -> GeneSetSpec {
  GeneSetSpec::{ name, up_genes, down_genes }
}

///|
/// Result of singscore computation.
pub struct SingscoreResult {
  pub sample_id : String
  pub score : Double
  pub up_score : Double
  pub down_score : Double
  pub up_dispersion : Double
  pub down_dispersion : Double
  pub n_up_found : Int
  pub n_down_found : Int
}

///|
/// Get sample ID.
pub fn SingscoreResult::get_sample_id(self : SingscoreResult) -> String {
  self.sample_id
}

///|
/// Get overall score.
pub fn SingscoreResult::get_score(self : SingscoreResult) -> Double {
  self.score
}

///|
/// Get up-regulated score.
pub fn SingscoreResult::get_up_score(self : SingscoreResult) -> Double {
  self.up_score
}

///|
/// Get down-regulated score.
pub fn SingscoreResult::get_down_score(self : SingscoreResult) -> Double {
  self.down_score
}

///|
/// Get up-regulated dispersion.
pub fn SingscoreResult::get_up_dispersion(self : SingscoreResult) -> Double {
  self.up_dispersion
}

///|
/// Get down-regulated dispersion.
pub fn SingscoreResult::get_down_dispersion(self : SingscoreResult) -> Double {
  self.down_dispersion
}

///|
/// Get number of up-regulated genes found.
pub fn SingscoreResult::get_n_up_found(self : SingscoreResult) -> Int {
  self.n_up_found
}

///|
/// Get number of down-regulated genes found.
pub fn SingscoreResult::get_n_down_found(self : SingscoreResult) -> Int {
  self.n_down_found
}

///|
/// Score a single sample against a gene set specification.
/// Score is bounded in [-0.5, 0.5] range.
pub fn singscore_score(
  sample : SampleExpression,
  spec : GeneSetSpec,
  sample_id : String,
) -> SingscoreResult {
  let (up_score, up_disp, n_up) = singscore_directional(sample, spec.up_genes, true)
  let (down_score, down_disp, n_down) = singscore_directional(sample, spec.down_genes, false)
  
  let n_dir = (if n_up > 0 { 1 } else { 0 }) + (if n_down > 0 { 1 } else { 0 })
  let total = if n_dir > 0 { (up_score + down_score) / n_dir.to_double() } else { 0.0 }
  
  SingscoreResult::{
    sample_id,
    score: total,
    up_score,
    down_score,
    up_dispersion: up_disp,
    down_dispersion: down_disp,
    n_up_found: n_up,
    n_down_found: n_down,
  }
}

///|
/// Score a sample against multiple gene sets.
pub fn singscore_score_multiple(
  sample : SampleExpression,
  specs : Array[GeneSetSpec],
) -> Array[SingscoreResult] {
  let results : Array[SingscoreResult] = Array::new()
  for spec in specs {
    results.push(singscore_score(sample, spec, sample.sample_id))
  }
  results
}

///|
/// Compute directional score (up or down regulated).
fn singscore_directional(
  sample : SampleExpression,
  genes : Array[String],
  is_up : Bool,
) -> (Double, Double, Int) {
  if genes.length() == 0 {
    return (0.0, 0.0, 0)
  }
  
  let n_genes = sample.gene_names.length()
  if n_genes == 0 {
    return (0.0, 0.0, 0)
  }
  
  // Build rank array (centered around 0)
  let ranks : Array[Double] = Array::make(n_genes, 0.0)
  // Create index array and sort by expression
  let indices : Array[Int] = Array::make(n_genes, 0)
  {
    let mut i = 0
    while i < n_genes {
      indices[i] = i
      i = i + 1
    }
  }
  
  // Simple insertion sort for ranking (handles ties by averaging)
  let sorted_indices = singscore_sort_indices(sample.expression)
  
  // Assign ranks (1-based, then center)
  {
    let mut i = 0
    while i < n_genes {
      ranks[sorted_indices[i]] = (i + 1).to_double()
      i = i + 1
    }
  }
  
  // Center ranks around 0
  let n = n_genes.to_double()
  {
    let mut i = 0
    while i < n_genes {
      ranks[i] = ranks[i] - (n + 1.0) / 2.0
      i = i + 1
    }
  }
  
  // Build gene name -> index map
  let gene_idx : Map[String, Int] = Map([], capacity=n_genes)
  {
    let mut i = 0
    while i < n_genes {
      gene_idx.set(sample.gene_names[i], i)
      i = i + 1
    }
  }
  
  // Sum ranks for specified genes
  let mut sum = 0.0
  let mut n_found = 0
  let mut i = 0
  while i < genes.length() {
    let gene = genes[i]
    match gene_idx.get(gene) {
      Some(idx) => {
        let rank_val = ranks[idx]
        sum = sum + (if is_up { rank_val } else { -rank_val })
        n_found = n_found + 1
      }
      None => ()
    }
    i = i + 1
  }
  
  if n_found == 0 {
    return (0.0, 0.0, 0)
  }
  
  // Normalize to [-0.5, 0.5] using (n-1) bound
  let nf = n_found.to_double()
  let nf1 = (n_found - 1).to_double()
  let score = if nf1 > 0.0 {
    sum / nf1
  } else {
    sum / nf
  }
  
  // Clip to [-0.5, 0.5]
  let clamped = if score > 0.5 {
    0.5
  } else if score < -0.5 {
    -0.5
  } else {
    score
  }
  
  // Compute dispersion (variance of ranks for found genes)
  let dispersion = singscore_compute_dispersion(sample, genes, is_up, gene_idx)
  
  (clamped, dispersion, n_found)
}

///|
/// Compute dispersion of ranks for found genes.
fn singscore_compute_dispersion(
  sample : SampleExpression,
  genes : Array[String],
  is_up : Bool,
  gene_idx : Map[String, Int],
) -> Double {
  let ranks_found : Array[Double] = Array::new()
  let mut i = 0
  while i < genes.length() {
    let gene = genes[i]
    match gene_idx.get(gene) {
      Some(idx) => {
        let rank = (idx + 1).to_double() - (sample.gene_names.length().to_double() + 1.0) / 2.0
        ranks_found.push(if is_up { rank } else { -rank })
      }
      None => ()
    }
    i = i + 1
  }
  
  let n = ranks_found.length()
  if n <= 1 {
    return 0.0
  }
  
  let (mean, var) = stat_mean_var(ranks_found)
  let _ = mean
  var
}

///|
/// Sort indices by expression value (ascending), return sorted indices.
fn singscore_sort_indices(values : Array[Double]) -> Array[Int] {
  let n = values.length()
  let indices : Array[Int] = Array::make(n, 0)
  {
    let mut i = 0
    while i < n {
      indices[i] = i
      i = i + 1
    }
  }
  
  // Quicksort
  singscore_quicksort(values, indices, 0, n - 1)
  indices
}

fn singscore_quicksort(
  values : Array[Double],
  indices : Array[Int],
  lo : Int,
  hi : Int,
) -> Unit {
  if lo < hi {
    let pivot = values[indices[hi]]
    let mut i = lo - 1
    let mut j = lo
    while j < hi {
      if values[indices[j]] <= pivot {
        i = i + 1
        let tmp = indices[i]
        indices[i] = indices[j]
        indices[j] = tmp
      }
      j = j + 1
    }
    let tmp = indices[i + 1]
    indices[i + 1] = indices[hi]
    indices[hi] = tmp
    
    let p = i + 1
    singscore_quicksort(values, indices, lo, p - 1)
    singscore_quicksort(values, indices, p + 1, hi)
  }
}

///|
/// Create example data for singscore testing.
pub fn singscore_create_example() -> (SampleExpression, GeneSetSpec) {
  let genes = [
    "TP53", "BRCA1", "BRCA2", "EGFR", "MYC", "KRAS",
    "PTEN", "AKT1", "MTOR", "JAK2", "STAT3", "NFKB1",
    "CDK4", "CDK6", "RB1", "E2F1", "APC", "CTNNB1",
    "SMAD4", "TGFBR1",
  ]
  let expression = [
    10.0, 8.5, 6.2, 15.0, 20.0, 5.0,
    3.2, 12.0, 18.0, 7.0, 9.5, 11.0,
    14.0, 13.0, 4.0, 16.0, 2.0, 1.5,
    2.5, 3.0,
  ]
  
  let sample = SampleExpression::new("sample_1", genes, expression)
  
  let up_genes = ["TP53", "BRCA1", "BRCA2", "EGFR", "MYC", "KRAS", "AKT1", "MTOR"]
  let down_genes = ["PTEN", "RB1", "APC", "CTNNB1"]
  
  let spec = GeneSetSpec::new("cell_cycle", up_genes, down_genes=down_genes)
  (sample, spec)
}
