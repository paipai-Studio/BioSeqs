///|
/// Cyclone - Cell cycle scoring using gene pairs.
/// Inspired by Bioconductor cyclone package.
///
/// Cyclone is a method for assigning cell cycle phase scores to
/// single-cell RNA-seq data using known gene pair markers.
///
/// Key concepts:
/// - Gene pairs: sets of gene pairs where one gene is up-regulated
///   and another is down-regulated in a specific phase
/// - Scoring: proportion of gene pairs matching expected expression
/// - Phase assignment: the phase with the highest score for each cell

///|
/// A gene pair for cell cycle scoring.
pub struct GenePair {
  phase : String
  up_gene : String
  down_gene : String
}

///|
/// Result of cell cycle scoring for a single cell.
pub struct CycloneResult {
  cell_id : String
  scores : Map[String, Double]
  phases : Array[String]
  assigned_phase : String
  g1_score : Double
  s_score : Double
  g2_score : Double
  m_score : Double
}

///|
/// Parameters for Cyclone scoring.
pub struct CycloneParams {
  min_score : Double
  n_top_pairs : Int
}

///|
/// Create default Cyclone parameters.
pub fn CycloneParams::new() -> CycloneParams {
  CycloneParams::{
    min_score: 0.0,
    n_top_pairs: 1,
  }
}

///|
/// Create custom Cyclone parameters.
pub fn CycloneParams::with_params(min_score : Double, n_top_pairs : Int) -> CycloneParams {
  CycloneParams::{ min_score, n_top_pairs }
}

///|
/// Get the standard cell cycle gene pairs for human cells.
/// Based on the cyclone package's marker gene pairs.
pub fn cyclone_get_gene_pairs() -> Array[GenePair] {
  [
    // G1 phase pairs
    GenePair::{ phase: "G1", up_gene: "CCND1", down_gene: "CCNE1" },
    GenePair::{ phase: "G1", up_gene: "CCND2", down_gene: "CCNE2" },
    GenePair::{ phase: "G1", up_gene: "CCND3", down_gene: "CCNE1" },
    GenePair::{ phase: "G1", up_gene: "RB1", down_gene: "E2F1" },
    GenePair::{ phase: "G1", up_gene: "CDKN1A", down_gene: "CDK2" },
    GenePair::{ phase: "G1", up_gene: "CDKN2A", down_gene: "CDK4" },
    GenePair::{ phase: "G1", up_gene: "ARF", down_gene: "MDM2" },
    GenePair::{ phase: "G1", up_gene: "CD44", down_gene: "CD24" },
    GenePair::{ phase: "G1", up_gene: "KLF4", down_gene: "LIN28A" },
    GenePair::{ phase: "G1", up_gene: "MYC", down_gene: "KRAS" },
    // S phase pairs
    GenePair::{ phase: "S", up_gene: "PCNA", down_gene: "CDK2" },
    GenePair::{ phase: "S", up_gene: "MCM2", down_gene: "CDK1" },
    GenePair::{ phase: "S", up_gene: "MCM3", down_gene: "CDK2" },
    GenePair::{ phase: "S", up_gene: "MCM4", down_gene: "CDK1" },
    GenePair::{ phase: "S", up_gene: "MCM5", down_gene: "CDK2" },
    GenePair::{ phase: "S", up_gene: "MCM6", down_gene: "CDK1" },
    GenePair::{ phase: "S", up_gene: "MCM7", down_gene: "CDK2" },
    GenePair::{ phase: "S", up_gene: "RRM1", down_gene: "TYMS" },
    GenePair::{ phase: "S", up_gene: "RRM2", down_gene: "TYMK" },
    GenePair::{ phase: "S", up_gene: "POLA1", down_gene: "POLA2" },
    GenePair::{ phase: "S", up_gene: "POLE", down_gene: "POLH" },
    GenePair::{ phase: "S", up_gene: "FEN1", down_gene: "RAD51" },
    // G2 phase pairs
    GenePair::{ phase: "G2", up_gene: "CCN1", down_gene: "CDK2" },
    GenePair::{ phase: "G2", up_gene: "CCN2", down_gene: "CDK1" },
    GenePair::{ phase: "G2", up_gene: "PLK1", down_gene: "AURKA" },
    GenePair::{ phase: "G2", up_gene: "WEE1", down_gene: "CDC25A" },
    GenePair::{ phase: "G2", up_gene: "CHEK1", down_gene: "CHEK2" },
    GenePair::{ phase: "G2", up_gene: "BRCA1", down_gene: "BRCA2" },
    GenePair::{ phase: "G2", up_gene: "RAD50", down_gene: "RAD51" },
    GenePair::{ phase: "G2", up_gene: "MRE11", down_gene: "NBN" },
    GenePair::{ phase: "G2", up_gene: "ATM", down_gene: "ATR" },
    // M phase pairs
    GenePair::{ phase: "M", up_gene: "CCN3", down_gene: "CDK1" },
    GenePair::{ phase: "M", up_gene: "CCN4", down_gene: "CDK2" },
    GenePair::{ phase: "M", up_gene: "BUB1", down_gene: "BUB1B" },
    GenePair::{ phase: "M", up_gene: "MAD2L1", down_gene: "MAD2L2" },
    GenePair::{ phase: "M", up_gene: "CDC20", down_gene: "FZR1" },
    GenePair::{ phase: "M", up_gene: "AURKB", down_gene: "AURKA" },
    GenePair::{ phase: "M", up_gene: "TTK", down_gene: "NEK2" },
    GenePair::{ phase: "M", up_gene: "KIF11", down_gene: "KIF20A" },
    GenePair::{ phase: "M", up_gene: "KIF2C", down_gene: "KIF2A" },
    GenePair::{ phase: "M", up_gene: "DLGAP5", down_gene: "ASPM" },
    GenePair::{ phase: "M", up_gene: "TPX2", down_gene: "XKL2" },
  ]
}

///|
/// Get the unique phases from gene pairs.
pub fn cyclone_get_phases(gene_pairs : Array[GenePair]) -> Array[String] {
  let phases : Array[String] = Array::new()
  for pair in gene_pairs {
    if !phases.contains(pair.phase) {
      phases.push(pair.phase)
    }
  }
  phases
}

///|
/// Count gene pairs for each phase.
pub fn cyclone_count_pairs_by_phase(gene_pairs : Array[GenePair]) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([])
  for pair in gene_pairs {
    let current = counts.get(pair.phase).unwrap_or(0)
    counts.set(pair.phase, current + 1)
  }
  counts
}

///|
/// Score a single cell for a specific phase using gene pairs.
/// The score is the proportion of gene pairs where the up_gene
/// has higher expression than the down_gene.
pub fn cyclone_score_phase(
  cell_expression : Array[Double],
  gene_names : Array[String],
  gene_pairs : Array[GenePair],
  phase : String,
) -> Double {
  let mut score = 0.0
  let mut n_pairs = 0

  for pair in gene_pairs {
    if pair.phase == phase {
      let up_idx = cyclone_find_gene_index(gene_names, pair.up_gene)
      let down_idx = cyclone_find_gene_index(gene_names, pair.down_gene)

      if up_idx >= 0 && down_idx >= 0 {
        let up_expr = cell_expression[up_idx]
        let down_expr = cell_expression[down_idx]
        if up_expr > down_expr {
          score = score + 1.0
        }
        n_pairs = n_pairs + 1
      }
    }
  }

  if n_pairs == 0 {
    0.0
  } else {
    score / n_pairs.to_double()
  }
}

///|
/// Score a single cell across all phases.
pub fn cyclone_score_cell(
  cell_expression : Array[Double],
  gene_names : Array[String],
  gene_pairs : Array[GenePair],
) -> Map[String, Double] {
  let phases = cyclone_get_phases(gene_pairs)
  let scores : Map[String, Double] = Map([], capacity=phases.length())

  for phase in phases {
    let score = cyclone_score_phase(cell_expression, gene_names, gene_pairs, phase)
    scores.set(phase, score)
  }

  scores
}

///|
/// Find the index of a gene in the gene names array.
fn cyclone_find_gene_index(gene_names : Array[String], gene : String) -> Int {
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
/// Assign cell cycle phase based on scores.
pub fn cyclone_assign_phase(scores : Map[String, Double]) -> String {
  let keys = scores.keys()
  let mut best_phase = ""
  let mut best_score = -1.0

  for key in keys {
    let score = scores.get(key).unwrap_or(0.0)
    if score > best_score {
      best_score = score
      best_phase = key
    }
  }

  best_phase
}

///|
/// Run Cyclone scoring on a single cell.
pub fn cyclone_score_single_cell(
  cell_expression : Array[Double],
  cell_id : String,
  gene_names : Array[String],
  gene_pairs : Array[GenePair],
) -> CycloneResult {
  let scores = cyclone_score_cell(cell_expression, gene_names, gene_pairs)
  let phases = cyclone_get_phases(gene_pairs)
  let assigned = cyclone_assign_phase(scores)

  CycloneResult::{
    cell_id,
    scores,
    phases,
    assigned_phase: assigned,
    g1_score: scores.get("G1").unwrap_or(0.0),
    s_score: scores.get("S").unwrap_or(0.0),
    g2_score: scores.get("G2").unwrap_or(0.0),
    m_score: scores.get("M").unwrap_or(0.0),
  }
}

///|
/// Run Cyclone scoring on multiple cells.
pub fn cyclone_score_cells(
  counts : Array[Array[Double]],
  cell_ids : Array[String],
  gene_names : Array[String],
) -> Array[CycloneResult] {
  let gene_pairs = cyclone_get_gene_pairs()
  let results : Array[CycloneResult] = Array::new()
  let n_cells = cell_ids.length()

  let mut i = 0
  while i < n_cells {
    let cell_expr : Array[Double] = Array::new()
    let n_genes = counts.length()
    let mut j = 0
    while j < n_genes {
      cell_expr.push(counts[j][i])
      j = j + 1
    }
    let result = cyclone_score_single_cell(cell_expr, cell_ids[i], gene_names, gene_pairs)
    results.push(result)
    i = i + 1
  }

  results
}

///|
/// Get a summary of cell cycle phase distribution.
pub fn cyclone_phase_distribution(results : Array[CycloneResult]) -> Map[String, Double] {
  let n = results.length()
  if n == 0 {
    return Map([])
  }

  let phase_counts : Map[String, Int] = Map([])
  for result in results {
    let current = phase_counts.get(result.assigned_phase).unwrap_or(0)
    phase_counts.set(result.assigned_phase, current + 1)
  }

  let distribution : Map[String, Double] = Map([])
  let keys = phase_counts.keys()
  for key in keys {
    let count = phase_counts.get(key).unwrap_or(0)
    distribution.set(key, count.to_double() / n.to_double())
  }

  distribution
}

///|
/// Create synthetic test data for Cyclone testing.
pub fn cyclone_create_test_data() -> (Array[Array[Double]], Array[String], Array[String]) {
  let gene_pairs = cyclone_get_gene_pairs()

  // Collect all unique gene names from pairs
  let gene_names : Array[String] = Array::new()
  for pair in gene_pairs {
    if !gene_names.contains(pair.up_gene) {
      gene_names.push(pair.up_gene)
    }
    if !gene_names.contains(pair.down_gene) {
      gene_names.push(pair.down_gene)
    }
  }

  let n_genes = gene_names.length()
  let n_cells = 20

  // Create expression data with phase-specific patterns
  let counts : Array[Array[Double]] = Array::new()
  let cell_ids : Array[String] = Array::new()

  // Create cells for each phase: 5 G1, 5 S, 5 G2, 5 M
  let gene_expr : Array[Array[Double]] = Array::new()
  let mut g = 0
  while g < n_genes {
    let col : Array[Double] = Array::new()
    let mut c = 0
    while c < n_cells {
      // Determine which phase this cell belongs to
      let cell_phase_idx = c / 5
      let cell_phase = if cell_phase_idx == 0 {
        "G1"
      } else if cell_phase_idx == 1 {
        "S"
      } else if cell_phase_idx == 2 {
        "G2"
      } else {
        "M"
      }

      // Check if this gene is an up_gene or down_gene for this phase
      let gene_name = gene_names[g]
      let mut expr = 0.5

      for pair in gene_pairs {
        if pair.phase == cell_phase {
          if pair.up_gene == gene_name {
            expr = 10.0 // High expression for up genes in active phase
          }
          if pair.down_gene == gene_name {
            expr = 0.1 // Low expression for down genes in active phase
          }
        }
      }

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
/// Get the average scores for each phase across all cells.
pub fn cyclone_average_scores(results : Array[CycloneResult]) -> Map[String, Double] {
  if results.length() == 0 {
    return Map([])
  }

  let phases = results[0].phases
  let avg_scores : Map[String, Double] = Map([])

  for phase in phases {
    let mut sum = 0.0
    for result in results {
      sum = sum + result.scores.get(phase).unwrap_or(0.0)
    }
    avg_scores.set(phase, sum / results.length().to_double())
  }

  avg_scores
}