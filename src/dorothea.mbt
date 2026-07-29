///|
/// dorothea - TF activity prediction using VIPER.
/// Inspired by Bioconductor dorothea package.
///
/// dorothea is a method for predicting transcription factor (TF) activity
/// from gene expression data using regulon enrichment analysis.
///
/// Key concepts:
/// - Regulon: set of target genes for a TF with direction (activation/repression)
/// - TF activity scoring: based on collective expression of target genes
/// - VIPER: statistical framework for TF activity assessment
/// - Permutation testing: for significance estimation

///|
/// A target gene in a regulon.
pub struct TargetGene {
  gene : String
  direction : String
}

///|
/// A regulon (TF and its target genes).
pub struct Regulon {
  tf : String
  targets : Array[TargetGene]
}

///|
/// TF activity result for a single TF.
pub struct TFActivityResult {
  tf : String
  activity : Double
  normalized_activity : Double
  p_value : Double
  z_score : Double
  n_targets : Int
}

///|
/// Parameters for dorothea TF activity scoring.
pub struct DorotheaParams {
  n_permutations : Int
  min_targets : Int
  significance : Double
}

///|
/// Create default dorothea parameters.
pub fn DorotheaParams::new() -> DorotheaParams {
  DorotheaParams::{
    n_permutations: 1000,
    min_targets: 5,
    significance: 0.05,
  }
}

///|
/// Create custom dorothea parameters.
pub fn DorotheaParams::with_params(n_perms : Int, min_targets : Int, alpha : Double) -> DorotheaParams {
  DorotheaParams::{
    n_permutations: n_perms,
    min_targets,
    significance: alpha,
  }
}

///|
/// Get a curated set of regulons for common human TFs.
/// Based on the dorothea package's regulon database.
pub fn dorothea_get_regulons() -> Array[Regulon] {
  [
    Regulon::{
      tf: "TP53",
      targets: [
        TargetGene::{ gene: "CDKN1A", direction: "activation" },
        TargetGene::{ gene: "BAX", direction: "activation" },
        TargetGene::{ gene: "BBC3", direction: "activation" },
        TargetGene::{ gene: "PMAIP1", direction: "activation" },
        TargetGene::{ gene: "SFN", direction: "activation" },
        TargetGene::{ gene: "CCNG1", direction: "activation" },
        TargetGene::{ gene: "MDM2", direction: "activation" },
        TargetGene::{ gene: "RRM2B", direction: "activation" },
        TargetGene::{ gene: "SESN1", direction: "activation" },
        TargetGene::{ gene: "SESN2", direction: "activation" },
        TargetGene::{ gene: "CDK4", direction: "repression" },
        TargetGene::{ gene: "CDK6", direction: "repression" },
      ],
    },
    Regulon::{
      tf: "MYC",
      targets: [
        TargetGene::{ gene: "CCND1", direction: "activation" },
        TargetGene::{ gene: "CCNE1", direction: "activation" },
        TargetGene::{ gene: "E2F1", direction: "activation" },
        TargetGene::{ gene: "RRM1", direction: "activation" },
        TargetGene::{ gene: "RRM2", direction: "activation" },
        TargetGene::{ gene: "TYMS", direction: "activation" },
        TargetGene::{ gene: "LDHA", direction: "activation" },
        TargetGene::{ gene: "HK2", direction: "activation" },
        TargetGene::{ gene: "PFKL", direction: "activation" },
        TargetGene::{ gene: "CDKN1A", direction: "repression" },
        TargetGene::{ gene: "CDKN2A", direction: "repression" },
      ],
    },
    Regulon::{
      tf: "JUN",
      targets: [
        TargetGene::{ gene: "CCL2", direction: "activation" },
        TargetGene::{ gene: "CXCL8", direction: "activation" },
        TargetGene::{ gene: "MMP1", direction: "activation" },
        TargetGene::{ gene: "MMP9", direction: "activation" },
        TargetGene::{ gene: "VEGFA", direction: "activation" },
        TargetGene::{ gene: "EGFR", direction: "activation" },
        TargetGene::{ gene: "FOS", direction: "activation" },
        TargetGene::{ gene: "NFKB1", direction: "activation" },
      ],
    },
    Regulon::{
      tf: "FOS",
      targets: [
        TargetGene::{ gene: "AP-1", direction: "activation" },
        TargetGene::{ gene: "MMP1", direction: "activation" },
        TargetGene::{ gene: "MMP3", direction: "activation" },
        TargetGene::{ gene: "MMP9", direction: "activation" },
        TargetGene::{ gene: "VEGFA", direction: "activation" },
        TargetGene::{ gene: "IL6", direction: "activation" },
        TargetGene::{ gene: "CXCL8", direction: "activation" },
      ],
    },
    Regulon::{
      tf: "NFKB1",
      targets: [
        TargetGene::{ gene: "IL6", direction: "activation" },
        TargetGene::{ gene: "IL1B", direction: "activation" },
        TargetGene::{ gene: "TNF", direction: "activation" },
        TargetGene::{ gene: "CXCL8", direction: "activation" },
        TargetGene::{ gene: "CCL2", direction: "activation" },
        TargetGene::{ gene: "VEGFA", direction: "activation" },
        TargetGene::{ gene: "BCL2", direction: "activation" },
        TargetGene::{ gene: "XIAP", direction: "activation" },
      ],
    },
    Regulon::{
      tf: "STAT3",
      targets: [
        TargetGene::{ gene: "SOCS3", direction: "activation" },
        TargetGene::{ gene: "BCL2", direction: "activation" },
        TargetGene::{ gene: "BCL2L1", direction: "activation" },
        TargetGene::{ gene: "MCL1", direction: "activation" },
        TargetGene::{ gene: "CCND1", direction: "activation" },
        TargetGene::{ gene: "MYC", direction: "activation" },
        TargetGene::{ gene: "IL17", direction: "activation" },
      ],
    },
    Regulon::{
      tf: "HIF1A",
      targets: [
        TargetGene::{ gene: "VEGFA", direction: "activation" },
        TargetGene::{ gene: "EPO", direction: "activation" },
        TargetGene::{ gene: "LDHA", direction: "activation" },
        TargetGene::{ gene: "HK2", direction: "activation" },
        TargetGene::{ gene: "PFKL", direction: "activation" },
        TargetGene::{ gene: "PGK1", direction: "activation" },
        TargetGene::{ gene: "CA9", direction: "activation" },
      ],
    },
    Regulon::{
      tf: "E2F1",
      targets: [
        TargetGene::{ gene: "CCNE1", direction: "activation" },
        TargetGene::{ gene: "CCND1", direction: "activation" },
        TargetGene::{ gene: "RRM1", direction: "activation" },
        TargetGene::{ gene: "RRM2", direction: "activation" },
        TargetGene::{ gene: "TYMS", direction: "activation" },
        TargetGene::{ gene: "POLA1", direction: "activation" },
        TargetGene::{ gene: "POLE", direction: "activation" },
        TargetGene::{ gene: "PCNA", direction: "activation" },
        TargetGene::{ gene: "CDKN1A", direction: "repression" },
      ],
    },
    Regulon::{
      tf: "SMAD4",
      targets: [
        TargetGene::{ gene: "CDKN1A", direction: "activation" },
        TargetGene::{ gene: "CDKN2B", direction: "activation" },
        TargetGene::{ gene: "FOS", direction: "activation" },
        TargetGene::{ gene: "JUN", direction: "activation" },
        TargetGene::{ gene: "MYC", direction: "repression" },
        TargetGene::{ gene: "CCND1", direction: "repression" },
      ],
    },
    Regulon::{
      tf: "NFE2L2",
      targets: [
        TargetGene::{ gene: "HMOX1", direction: "activation" },
        TargetGene::{ gene: "NQO1", direction: "activation" },
        TargetGene::{ gene: "GCLC", direction: "activation" },
        TargetGene::{ gene: "GCLM", direction: "activation" },
        TargetGene::{ gene: "TXNRD1", direction: "activation" },
        TargetGene::{ gene: "GSR", direction: "activation" },
        TargetGene::{ gene: "KEAP1", direction: "repression" },
      ],
    },
    Regulon::{
      tf: "SOX2",
      targets: [
        TargetGene::{ gene: "NANOG", direction: "activation" },
        TargetGene::{ gene: "OCT4", direction: "activation" },
        TargetGene::{ gene: "KLF4", direction: "activation" },
        TargetGene::{ gene: "MYC", direction: "activation" },
        TargetGene::{ gene: "LIN28A", direction: "activation" },
      ],
    },
    Regulon::{
      tf: "OCT4",
      targets: [
        TargetGene::{ gene: "SOX2", direction: "activation" },
        TargetGene::{ gene: "NANOG", direction: "activation" },
        TargetGene::{ gene: "KLF4", direction: "activation" },
        TargetGene::{ gene: "MYC", direction: "activation" },
        TargetGene::{ gene: "CD45", direction: "activation" },
      ],
    },
    Regulon::{
      tf: "GATA3",
      targets: [
        TargetGene::{ gene: "IL4", direction: "activation" },
        TargetGene::{ gene: "IL5", direction: "activation" },
        TargetGene::{ gene: "IL13", direction: "activation" },
        TargetGene::{ gene: "CCR8", direction: "activation" },
        TargetGene::{ gene: "CXCR4", direction: "activation" },
      ],
    },
  ]
}

///|
/// Find the index of a gene in the gene names array.
fn dorothea_find_gene_index(gene_names : Array[String], gene : String) -> Int {
  let n = gene_names.length()
  let mut i = 0
  while i < n {
    if gene_names[i] == gene {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Get the target gene expression values for a regulon.
pub fn dorothea_get_target_expression(
  cell_expression : Array[Double],
  gene_names : Array[String],
  regulon : Regulon,
) -> Array[(String, Double)] {
  let result : Array[(String, Double)] = Array::new()
  for target in regulon.targets {
    let idx = dorothea_find_gene_index(gene_names, target.gene)
    if idx >= 0 {
      result.push((target.direction, cell_expression[idx]))
    }
  }
  result
}

///|
/// Compute simple TF activity score based on target gene expression.
/// The score is the weighted mean of target gene expression,
/// where activation targets get positive weight and repression targets get negative weight.
pub fn dorothea_compute_activity_simple(
  cell_expression : Array[Double],
  gene_names : Array[String],
  regulon : Regulon,
) -> Double {
  let mut sum = 0.0
  let mut n = 0
  for target in regulon.targets {
    let idx = dorothea_find_gene_index(gene_names, target.gene)
    if idx >= 0 {
      let expr = cell_expression[idx]
      if target.direction == "activation" {
        sum = sum + expr
      } else {
        sum = sum - expr
      }
      n = n + 1
    }
  }
  if n == 0 {
    0.0
  } else {
    sum / n.to_double()
  }
}

///|
/// Compute VIPER-like TF activity score.
/// Uses a more sophisticated weighting based on target gene importance.
pub fn dorothea_compute_viper_activity(
  cell_expression : Array[Double],
  gene_names : Array[String],
  regulon : Regulon,
) -> Double {
  let mut activation_sum = 0.0
  let mut repression_sum = 0.0
  let mut n_activation = 0
  let mut n_repression = 0

  for target in regulon.targets {
    let idx = dorothea_find_gene_index(gene_names, target.gene)
    if idx >= 0 {
      let expr = cell_expression[idx]
      if target.direction == "activation" {
        activation_sum = activation_sum + expr
        n_activation = n_activation + 1
      } else {
        repression_sum = repression_sum + expr
        n_repression = n_repression + 1
      }
    }
  }

  let activation_mean = if n_activation > 0 { activation_sum / n_activation.to_double() } else { 0.0 }
  let repression_mean = if n_repression > 0 { repression_sum / n_repression.to_double() } else { 0.0 }

  // VIPER score: difference between activation and repression means
  activation_mean - repression_mean
}

///|
/// Perform permutation test for TF activity significance.
pub fn dorothea_permutation_test(
  cell_expression : Array[Double],
  gene_names : Array[String],
  regulon : Regulon,
  n_permutations : Int,
) -> (Double, Double) {
  let observed = dorothea_compute_viper_activity(cell_expression, gene_names, regulon)

  // Generate null distribution by shuffling gene labels
  let null_scores : Array[Double] = Array::new()
  let n_genes = cell_expression.length()

  let mut perm = 0
  while perm < n_permutations {
    // Create shuffled expression by randomly reordering
    let shuffled : Array[Double] = Array::make(n_genes, 0.0)
    // Simple shuffle: swap random pairs
    let mut i = 0
    while i < n_genes {
      shuffled[i] = cell_expression[i]
      i = i + 1
    }
    // Fisher-Yates shuffle
    let mut j = n_genes - 1
    while j > 0 {
      let k = (j.to_double() * dorothea_pseudo_random(perm * n_genes + j)).to_int()
      let temp = shuffled[j]
      shuffled[j] = shuffled[k]
      shuffled[k] = temp
      j = j - 1
    }

    let null_score = dorothea_compute_viper_activity(shuffled, gene_names, regulon)
    null_scores.push(null_score)
    perm = perm + 1
  }

  // Compute p-value as proportion of null scores >= observed
  let mut count_extreme = 0
  for null_score in null_scores {
    if null_score >= observed {
      count_extreme = count_extreme + 1
    }
  }

  let p_value = count_extreme.to_double() / n_permutations.to_double()

  // Compute z-score
  let null_mean = dorothea_array_mean(null_scores)
  let null_std = dorothea_array_std(null_scores)
  let z_score = if null_std > 0.0 {
    (observed - null_mean) / null_std
  } else {
    0.0
  }

  (p_value, z_score)
}

///|
/// Simple pseudo-random number generator (deterministic for reproducibility).
fn dorothea_pseudo_random(seed : Int) -> Double {
  let x = (seed * 1103515245 + 12345) & 0x7FFFFFFF
  (x.to_double() / 2147483647.0).abs()
}

///|
/// Compute mean of an array.
pub fn dorothea_array_mean(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for val in arr {
    sum = sum + val
  }
  sum / n.to_double()
}

///|
/// Compute standard deviation of an array.
pub fn dorothea_array_std(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n < 2 {
    return 0.0
  }
  let mean = dorothea_array_mean(arr)
  let mut sum_sq = 0.0
  for val in arr {
    let diff = val - mean
    sum_sq = sum_sq + diff * diff
  }
  (sum_sq / (n - 1).to_double()).sqrt()
}

///|
/// Run dorothea TF activity analysis on a single cell.
pub fn dorothea_analyze_cell(
  cell_expression : Array[Double],
  gene_names : Array[String],
  regulons : Array[Regulon],
  params : DorotheaParams,
) -> Array[TFActivityResult] {
  let results : Array[TFActivityResult] = Array::new()

  for regulon in regulons {
    if regulon.targets.length() < params.min_targets {
      continue
    }

    let activity = dorothea_compute_viper_activity(cell_expression, gene_names, regulon)
    let (p_value, z_score) = dorothea_permutation_test(
      cell_expression,
      gene_names,
      regulon,
      params.n_permutations,
    )

    results.push(TFActivityResult::{
      tf: regulon.tf,
      activity,
      normalized_activity: z_score,
      p_value,
      z_score,
      n_targets: regulon.targets.length(),
    })
  }

  results
}

///|
/// Run dorothea analysis on multiple cells and return averaged results.
pub fn dorothea_analyze_cells(
  counts : Array[Array[Double]],
  cell_ids : Array[String],
  gene_names : Array[String],
  regulons : Array[Regulon],
  params : DorotheaParams,
) -> Array[TFActivityResult] {
  let n_cells = cell_ids.length()
  let n_regulons = regulons.length()

  // Accumulate activities across cells
  let tf_activities : Map[String, Array[Double]] = Map([], capacity=n_regulons)

  for regulon in regulons {
    if regulon.targets.length() >= params.min_targets {
      tf_activities.set(regulon.tf, [])
    }
  }

  // Process each cell
  let mut c = 0
  while c < n_cells {
    let cell_expr : Array[Double] = Array::new()
    let n_genes = counts.length()
    let mut g = 0
    while g < n_genes {
      cell_expr.push(counts[g][c])
      g = g + 1
    }

    for regulon in regulons {
      if regulon.targets.length() >= params.min_targets {
        let activity = dorothea_compute_viper_activity(cell_expr, gene_names, regulon)
        let current = tf_activities.get(regulon.tf).unwrap_or([])
        let updated = current.copy()
        updated.push(activity)
        tf_activities.set(regulon.tf, updated)
      }
    }
    c = c + 1
  }

  // Compute summary statistics
  let results : Array[TFActivityResult] = Array::new()
  let keys = tf_activities.keys()
  for tf in keys {
    let activities = tf_activities.get(tf).unwrap_or([])
    let mean_activity = dorothea_array_mean(activities)
    let std_activity = dorothea_array_std(activities)

    // Find the corresponding regulon for n_targets
    let mut n_targets = 0
    for regulon in regulons {
      if regulon.tf == tf {
        n_targets = regulon.targets.length()
        break
      }
    }

    // Simple significance: z-score of mean activity
    let z_score = if std_activity > 0.0 {
      mean_activity / std_activity
    } else {
      0.0
    }

    results.push(TFActivityResult::{
      tf,
      activity: mean_activity,
      normalized_activity: z_score,
      p_value: 0.0,
      z_score,
      n_targets,
    })
  }

  results
}

///|
/// Filter TF activity results by significance.
pub fn dorothea_filter_results(
  results : Array[TFActivityResult],
  alpha : Double,
) -> Array[TFActivityResult] {
  let filtered : Array[TFActivityResult] = Array::new()
  for result in results {
    if result.p_value <= alpha {
      filtered.push(result)
    }
  }
  filtered
}

///|
/// Sort TF activity results by absolute z-score (descending).
pub fn dorothea_sort_by_activity(results : Array[TFActivityResult]) -> Array[TFActivityResult] {
  // Simple bubble sort by abs(z_score) descending
  let n = results.length()
  let sorted = results.copy()
  let mut i = 0
  while i < n - 1 {
    let mut j = 0
    while j < n - i - 1 {
      if (sorted[j].z_score).abs() < (sorted[j + 1].z_score).abs() {
        let temp = sorted[j]
        sorted[j] = sorted[j + 1]
        sorted[j + 1] = temp
      }
      j = j + 1
    }
    i = i + 1
  }
  sorted
}

///|
/// Create synthetic test data for dorothea testing.
pub fn dorothea_create_test_data() -> (Array[Array[Double]], Array[String], Array[String]) {
  let regulons = dorothea_get_regulons()

  // Collect all unique genes from regulons
  let gene_names : Array[String] = Array::new()
  for regulon in regulons {
    for target in regulon.targets {
      if !gene_names.contains(target.gene) {
        gene_names.push(target.gene)
      }
    }
  }

  let n_genes = gene_names.length()
  let n_cells = 10

  // Create expression data
  let gene_expr : Array[Array[Double]] = Array::new()
  let cell_ids : Array[String] = Array::new()

  let mut g = 0
  while g < n_genes {
    let col : Array[Double] = Array::new()
    let mut c = 0
    while c < n_cells {
      // Vary expression across cells with some structure
      let base = 5.0
      let variation = (c.to_double() + 1.0) * 0.5
      let expr = base + variation + (g.to_double() * 0.01)
      col.push(expr)
      c = c + 1
    }
    gene_expr.push(col)
    g = g + 1
  }

  // Generate cell IDs
  let mut idx = 0
  while idx < n_cells {
    cell_ids.push("cell_\{(idx + 1).to_string()}")
    idx = idx + 1
  }

  (gene_expr, cell_ids, gene_names)
}

///|
/// Get top TF results by absolute activity.
pub fn dorothea_get_top_tfs(results : Array[TFActivityResult], n : Int) -> Array[TFActivityResult] {
  let sorted = dorothea_sort_by_activity(results)
  let top : Array[TFActivityResult] = Array::new()
  let limit = if n < sorted.length() { n } else { sorted.length() }
  let mut i = 0
  while i < limit {
    top.push(sorted[i])
    i = i + 1
  }
  top
}