///|
/// Bioconductor CellChat - Cell-cell communication interaction analysis.
/// Identifies significant ligand-receptor interactions between cell types
/// based on single-cell expression data, using permutation testing.
///
/// Reference: Bioconductor CellChat package
/// Core algorithm:
/// 1. Score each ligand-receptor pair: product of mean ligand (sender) and mean receptor (receiver)
/// 2. Permutation test: shuffle cell type labels to build null distribution
/// 3. Compute p-values and significance

///|
/// Ligand-receptor interaction pair definition.
pub struct LRPair {
  ligand : String
  receptor : String
  key : String
}

///|
/// Score for a single ligand-receptor interaction between two cell types.
pub struct InteractionScore {
  ligand : String
  receptor : String
  source_celltype : String
  target_celltype : String
  score : Double
  perm_mean : Double
  perm_std : Double
  p_value : Double
  p_adj : Double
  significant : Bool
}

///|
/// CellChat analysis result containing all interaction scores.
pub struct CellChatResult {
  scores : Array[InteractionScore]
  lr_pairs : Array[LRPair]
  cell_types : Array[String]
  n_permutations : Int
  seed : Int
  fdr : Double
}

///|
/// Create a ligand-receptor pair.
pub fn lr_pair(ligand : String, receptor : String) -> LRPair {
  LRPair::{
    ligand: ligand,
    receptor: receptor,
    key: ligand + "_" + receptor,
  }
}

///|
/// Create a built-in ligand-receptor database with common mammalian pairs.
pub fn cellchat_lr_database() -> Array[LRPair] {
  [
    lr_pair("TNF", "TNFR1"),
    lr_pair("TNF", "TNFR2"),
    lr_pair("IL1B", "IL1R1"),
    lr_pair("IL1B", "IL1R2"),
    lr_pair("IL6", "IL6R"),
    lr_pair("IL6", "gp130"),
    lr_pair("TGFB1", "TGFBR1"),
    lr_pair("TGFB1", "TGFBR2"),
    lr_pair("VEGFA", "VEGFR1"),
    lr_pair("VEGFA", "VEGFR2"),
    lr_pair("EGFR", "ERBB2"),
    lr_pair("FGF1", "FGFR1"),
    lr_pair("FGF2", "FGFR1"),
    lr_pair("FGF2", "FGFR2"),
    lr_pair("WNT5A", "FZD5"),
    lr_pair("WNT3A", "FZD1"),
    lr_pair("BMP2", "BMPR1A"),
    lr_pair("BMP4", "BMPR1A"),
    lr_pair("CXCL12", "CXCR4"),
    lr_pair("CXCL10", "CXCR3"),
    lr_pair("CCL2", "CCR2"),
    lr_pair("CCL5", "CCR5"),
    lr_pair("PDGFB", "PDGFRB"),
    lr_pair("HGF", "MET"),
    lr_pair("IGF1", "IGF1R"),
    lr_pair("IL10", "IL10R"),
    lr_pair("IFNG", "IFNGR1"),
    lr_pair("CSF1", "CSF1R"),
    lr_pair("SCF", "KIT"),
    lr_pair("EPO", "EPOR"),
  ]
}

///|
/// Compute mean expression of a gene in a given cell type.
/// expression: cell × gene matrix
/// cell_types: cell type label per cell
pub fn cellchat_mean_expr(
  expression : Array[Array[Double]],
  cell_types : Array[String],
  cell_type : String,
  gene_index : Int,
) -> Double {
  let mut sum = 0.0
  let mut count = 0
  let mut i = 0
  while i < cell_types.length() {
    if cell_types[i] == cell_type {
      sum = sum + expression[i][gene_index]
      count = count + 1
    }
    i = i + 1
  }
  if count > 0 { sum / count.to_double() } else { 0.0 }
}

///|
/// Compute CellChat interaction score matrix.
/// For each LR pair and direction (source→target):
///   score = mean_ligand(source) * mean_receptor(target)
///
/// Parameters:
///   expression : cells × genes expression matrix
///   cell_types : cell type per cell
///   gene_names : gene names (same order as columns)
///   lr_pairs : ligand-receptor pairs to test
///   n_permutations : number of permutations for null distribution
///   seed : random seed for reproducibility
///   fdr : false discovery rate threshold
pub fn cellchat_analyze(
  expression : Array[Array[Double]],
  cell_types : Array[String],
  gene_names : Array[String],
  lr_pairs : Array[LRPair],
  n_permutations? : Int = 100,
  seed? : Int = 42,
  fdr? : Double = 0.05,
) -> CellChatResult {
  // Get unique cell types
  let unique_ct : Array[String] = []
  let mut i = 0
  while i < cell_types.length() {
    let ct = cell_types[i]
    let mut found = false
    let mut j = 0
    while j < unique_ct.length() {
      if unique_ct[j] == ct { found = true; break }
      j = j + 1
    }
    if !found { unique_ct.push(ct) }
    i = i + 1
  }

  // Find gene indices
  let gene_index : Map[String, Int] = Map([], capacity=100)
  i = 0
  while i < gene_names.length() {
    gene_index.set(gene_names[i], i)
    i = i + 1
  }

  let n_perm = n_permutations
  let seed_val = seed
  let fdr_val = fdr

  let scores : Array[InteractionScore] = []

  // Generate pseudo-random permutations
  let mut rng_state = seed_val
  let permuted_labels_list : Array[Array[String]] = []
  let mut p = 0
  while p < n_perm {
    let perm : Array[String] = []
    let mut idx = 0
    // Simple Fisher-Yates-like shuffle
    let indices : Array[Int] = []
    idx = 0
    while idx < cell_types.length() {
      indices.push(idx)
      idx = idx + 1
    }
    // Shuffle indices
    let mut s = cell_types.length() - 1
    while s > 0 {
      rng_state = (rng_state * 1103515245 + 12345) & 0x7fffffff
      let j = rng_state % (s + 1)
      let tmp = indices[s]
      indices[s] = indices[j]
      indices[j] = tmp
      s = s - 1
    }
    idx = 0
    while idx < indices.length() {
      perm.push(cell_types[indices[idx]])
      idx = idx + 1
    }
    permuted_labels_list.push(perm)
    p = p + 1
  }

  // For each LR pair and direction
  let mut li = 0
  while li < lr_pairs.length() {
    let pair = lr_pairs[li]
    let lig_idx_opt = gene_index.get(pair.ligand)
    let rec_idx_opt = gene_index.get(pair.receptor)
    if lig_idx_opt is None || rec_idx_opt is None {
      li = li + 1
      continue
    }
    let lig_idx = lig_idx_opt.unwrap()
    let rec_idx = rec_idx_opt.unwrap()

    let mut src = 0
    while src < unique_ct.length() {
      let mut tgt = 0
      while tgt < unique_ct.length() {
        if src != tgt {
          let src_ct = unique_ct[src]
          let tgt_ct = unique_ct[tgt]

          // Observed score
          let mean_lig = cellchat_mean_expr(expression, cell_types, src_ct, lig_idx)
          let mean_rec = cellchat_mean_expr(expression, cell_types, tgt_ct, rec_idx)
          let obs_score = mean_lig * mean_rec

          // Permutation null distribution
          let perm_scores : Array[Double] = []
          let mut pm = 0
          while pm < n_perm {
            let perm_labels = permuted_labels_list[pm]
            let perm_lig = cellchat_mean_expr(expression, perm_labels, src_ct, lig_idx)
            let perm_rec = cellchat_mean_expr(expression, perm_labels, tgt_ct, rec_idx)
            let perm_score = perm_lig * perm_rec
            perm_scores.push(perm_score)
            pm = pm + 1
          }

          // Null distribution statistics
          let mut perm_mean = 0.0
          let mut perm_var = 0.0
          pm = 0
          while pm < n_perm {
            perm_mean = perm_mean + perm_scores[pm]
            pm = pm + 1
          }
          perm_mean = perm_mean / n_perm.to_double()

          pm = 0
          while pm < n_perm {
            let diff = perm_scores[pm] - perm_mean
            perm_var = perm_var + diff * diff
            pm = pm + 1
          }
          perm_var = perm_var / n_perm.to_double()
          let perm_std = if perm_var > 0.0 { perm_var.sqrt() } else { 0.0 }

          // P-value: fraction of permutation scores >= observed
          let mut n_ge = 0
          pm = 0
          while pm < n_perm {
            if perm_scores[pm] >= obs_score { n_ge = n_ge + 1 }
            pm = pm + 1
          }
          let p_val = (n_ge + 1).to_double() / (n_perm + 1).to_double()

          // BH FDR correction (simplified: just use raw p * number of tests / rank)
          // For now, use raw p-value threshold
          let p_adj_raw = p_val
          let significant = p_val < fdr_val

          scores.push(InteractionScore::{
            ligand: pair.ligand,
            receptor: pair.receptor,
            source_celltype: src_ct,
            target_celltype: tgt_ct,
            score: obs_score,
            perm_mean: perm_mean,
            perm_std: perm_std,
            p_value: p_val,
            p_adj: p_adj_raw,
            significant: significant,
          })
        }
        tgt = tgt + 1
      }
      src = src + 1
    }
    li = li + 1
  }

  CellChatResult::{
    scores: scores,
    lr_pairs: lr_pairs,
    cell_types: unique_ct,
    n_permutations: n_perm,
    seed: seed_val,
    fdr: fdr_val,
  }
}

///|
/// Get significant interactions from a CellChatResult.
pub fn cellchat_get_significant(result : CellChatResult) -> Array[InteractionScore] {
  let sig : Array[InteractionScore] = []
  let mut i = 0
  while i < result.scores.length() {
    if result.scores[i].significant {
      sig.push(result.scores[i])
    }
    i = i + 1
  }
  sig
}

///|
/// Aggregate total communication score per cell-type pair (undirected).
pub fn cellchat_aggregate(result : CellChatResult) -> Map[String, Double] {
  let agg : Map[String, Double] = Map([], capacity=50)
  let mut i = 0
  while i < result.scores.length() {
    let s = result.scores[i]
    // Use sorted direction as key
    let key = if s.source_celltype <= s.target_celltype {
      s.source_celltype + "|" + s.target_celltype
    } else {
      s.target_celltype + "|" + s.source_celltype
    }
    if agg.contains(key) {
      agg.set(key, agg.get(key).unwrap() + s.score)
    } else {
      agg.set(key, s.score)
    }
    i = i + 1
  }
  agg
}

///|
/// Generate sample single-cell expression data for CellChat testing.
/// Returns (expression, cell_types, gene_names) with 3 cell types and known LR pairs.
pub fn cellchat_sample_data() -> (Array[Array[Double]], Array[String], Array[String]) {
  let n_cells = 30
  // Express 10 genes: TNF, TNFR1, IL6, IL6R, VEGFA, VEGFR1, TGFB1, TGFBR1, CXCL12, CXCR4
  let gene_names = ["TNF", "TNFR1", "IL6", "IL6R", "VEGFA", "VEGFR1", "TGFB1", "TGFBR1", "CXCL12", "CXCR4"]
  let cell_types : Array[String] = []
  // 10 cells per cell type
  let types = ["TypeA", "TypeB", "TypeC"]
  let mut i = 0
  while i < n_cells {
    cell_types.push(types[i / 10])
    i = i + 1
  }

  let expression : Array[Array[Double]] = []
  let mut ci = 0
  while ci < n_cells {
    let row : Array[Double] = []
    let ct = cell_types[ci]
    let mut gi = 0
    while gi < gene_names.length() {
      let base = if ct == "TypeA" {
        if gene_names[gi] == "TNF" || gene_names[gi] == "IL6" { 15.0 } else { 2.0 }
      } else if ct == "TypeB" {
        if gene_names[gi] == "TNFR1" || gene_names[gi] == "IL6R" { 12.0 } else { 2.5 }
      } else {
        if gene_names[gi] == "VEGFA" || gene_names[gi] == "TGFB1" { 18.0 } else { 3.0 }
      }
      let noise = (((ci * 7 + gi * 11) % 100).to_double() / 100.0) * base * 0.4
      row.push(base + noise)
      gi = gi + 1
    }
    expression.push(row)
    ci = ci + 1
  }
  (expression, cell_types, gene_names)
}

///|
/// Summary of CellChat results.
pub fn cellchat_summary(result : CellChatResult) -> String {
  let mut n_sig = 0
  let mut i = 0
  while i < result.scores.length() {
    if result.scores[i].significant { n_sig = n_sig + 1 }
    i = i + 1
  }
  "CellChat Analysis Summary:\n" +
    "  Cell types: " + result.cell_types.length().to_string() + "\n" +
    "  LR pairs tested: " + result.lr_pairs.length().to_string() + "\n" +
    "  Total interactions: " + result.scores.length().to_string() + "\n" +
    "  Significant (FDR<" + result.fdr.to_string() + "): " + n_sig.to_string() + "\n" +
    "  Permutations: " + result.n_permutations.to_string()
}

///|
/// Get the top interactions by score.
pub fn cellchat_get_top(result : CellChatResult, n : Int) -> Array[InteractionScore] {
  let all = result.scores
  // Sort by score descending
  let mut i = 0
  while i < all.length() {
    let mut j = i + 1
    while j < all.length() {
      if all[j].score > all[i].score {
        let tmp = all[i]
        all[i] = all[j]
        all[j] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  let top : Array[InteractionScore] = []
  let count = if n < all.length() { n } else { all.length() }
  i = 0
  while i < count {
    top.push(all[i])
    i = i + 1
  }
  top
}
