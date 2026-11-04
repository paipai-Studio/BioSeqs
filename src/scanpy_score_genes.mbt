/// Scanpy-inspired: tl.score_genes, tl.score_genes_cell_cycle — gene set scoring

///|
/// Score cells for a set of genes.
/// The score is the mean expression of genes in the set minus the mean
/// expression of a reference (control) gene set.
///
/// Parameters:
/// - `x`: cells × genes matrix (log-normalized)
/// - `var_names`: gene names
/// - `gene_list`: genes to score
/// - `n_bins`: number of bins for control gene selection (default 25)
/// - `ctrl_size`: control set size as fraction of gene_list (default 1.0)
pub fn scanpy_score_genes(
  x : Array[Array[Double]],
  var_names : Array[String],
  gene_list : Array[String],
  n_bins~ : Int32 = 25,
  ctrl_size~ : Double = 1.0,
) -> Array[Double] {
  let n_cells = x.length()
  if n_cells == 0 {
    return []
  }
  let n_genes = x[0].length()
  // Map gene names to indices
  let name_to_idx : Map[String, Int32] = Map([], capacity=n_genes)
  for j in 0..<n_genes {
    name_to_idx[var_names[j]] = j
  }
  // Find target genes
  let target_indices = Array(capacity=0)
  for g in gene_list {
    match name_to_idx.get(g) {
      Some(idx) => target_indices.push(idx)
      None => ()
    }
  }
  if target_indices.length() == 0 {
    return {
      let a : Array[Double] = []
      for _ in 0..<n_cells {
        a.push(0.0)
      }
      a
    }
  }
  // Compute mean expression per gene for binning
  let gene_means = Array(capacity=n_genes)
  for j in 0..<n_genes {
    let mut s = 0.0
    for i in 0..<n_cells {
      s += x[i][j]
    }
    gene_means.push(s / n_cells.to_double())
  }
  // Bin genes by mean expression
  let mut min_mean = 1.0e100
  let mut max_mean = -1.0e100
  for v in gene_means {
    if v < min_mean {
      min_mean = v
    }
    if v > max_mean {
      max_mean = v
    }
  }
  let bin_width = if max_mean > min_mean {
    (max_mean - min_mean) / n_bins.to_double()
  } else {
    1.0
  }
  // For each target gene, select control genes from the same bin
  let ctrl_indices = Array(capacity=0)
  let _ctrl_count = (target_indices.length().to_double() * ctrl_size)
    .floor()
    .to_int()
  let mut rng : Int32 = 42
  for t in target_indices {
    let bin_idx = ((gene_means[t] - min_mean) / bin_width).floor().to_int()
    let lo = min_mean + bin_idx.to_double() * bin_width
    let hi = lo + bin_width
    // Collect genes in same bin (excluding target genes)
    let bin_genes = Array(capacity=0)
    for j in 0..<n_genes {
      if gene_means[j] >= lo && gene_means[j] < hi {
        let mut is_target = false
        for ti in target_indices {
          if j == ti {
            is_target = true
          }
        }
        if !is_target {
          bin_genes.push(j)
        }
      }
    }
    if bin_genes.length() > 0 {
      rng = (rng * 1103515245 + 12345) & 0x7FFFFFFF
      let pick = rng.reinterpret_as_uint() %
        bin_genes.length().reinterpret_as_uint()
      ctrl_indices.push(bin_genes[pick.reinterpret_as_int()])
    }
  }
  // Compute scores: mean(target) - mean(control) per cell
  let scores = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let mut target_sum = 0.0
    for t in target_indices {
      target_sum += x[i][t]
    }
    let target_mean = target_sum / target_indices.length().to_double()
    let mut ctrl_mean = 0.0
    if ctrl_indices.length() > 0 {
      let mut ctrl_sum = 0.0
      for c in ctrl_indices {
        ctrl_sum += x[i][c]
      }
      ctrl_mean = ctrl_sum / ctrl_indices.length().to_double()
    }
    scores.push(target_mean - ctrl_mean)
  }
  scores
}

///|
/// Score cell cycle phases (S and G2M) based on gene lists.
///
/// Parameters:
/// - `x`: cells × genes matrix (log-normalized)
/// - `var_names`: gene names
/// - `s_genes`: S-phase marker genes
/// - `g2m_genes`: G2/M-phase marker genes
///
/// Returns: (s_scores, g2m_scores, phase) where phase is "S", "G2M", or "G1"
pub fn scanpy_score_genes_cell_cycle(
  x : Array[Array[Double]],
  var_names : Array[String],
  s_genes : Array[String],
  g2m_genes : Array[String],
) -> (Array[Double], Array[Double], Array[String]) {
  let s_scores = scanpy_score_genes(x, var_names, s_genes)
  let g2m_scores = scanpy_score_genes(x, var_names, g2m_genes)
  // Assign phase: S if S > G2M, G2M if G2M > S, else G1
  let phases = Array(capacity=s_scores.length())
  for i in 0..<s_scores.length() {
    if s_scores[i] > g2m_scores[i] && s_scores[i] > 0.0 {
      phases.push("S")
    } else if g2m_scores[i] > s_scores[i] && g2m_scores[i] > 0.0 {
      phases.push("G2M")
    } else {
      phases.push("G1")
    }
  }
  (s_scores, g2m_scores, phases)
}
