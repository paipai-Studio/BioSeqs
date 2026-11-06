/// Scanpy second-round: Gene set scoring extensions
/// Inspired by scanpy.tl.score_genes (random_state, n_bins, ctrl_size),
/// score_genes_cell_cycle (advanced), gene_set_enrichment

///|
/// Score genes with control gene set.
/// Selects control genes matched by expression bins, computes score
/// as mean(gene_list) - mean(control_set).
pub fn scanpy_score_genes_with_control(
  x : Array[Array[Double]],
  gene_names : Array[String],
  gene_list : Array[String],
  n_bins~ : Int32 = 25,
  ctrl_size~ : Int32 = 50,
  random_state~ : Int32 = 0,
) -> Array[Double] {
  let n_cells = x.length()
  let n_genes = gene_names.length()
  // Find indices of gene_list
  let gene_set = Set([])
  for g in gene_list {
    gene_set.add(g)
  }
  let target_indices = Array(capacity=0)
  for g in 0..<n_genes {
    if gene_set.contains(gene_names[g]) {
      target_indices.push(g)
    }
  }
  if target_indices.length() == 0 {
    return Array(capacity=n_cells)
  }
  // Compute mean expression per gene for binning
  let gene_means = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let mut sum = 0.0
    for i in 0..<n_cells {
      sum = sum + x[i][g]
    }
    gene_means.push(if n_cells > 0 { sum / n_cells.to_double() } else { 0.0 })
  }
  // Bin genes by mean expression
  let sorted_means = Array(capacity=n_genes)
  for g in 0..<n_genes {
    sorted_means.push((g, gene_means[g]))
  }
  sorted_means.sort_by(fn(a, b) {
    if a.1 < b.1 {
      -1
    } else if a.1 > b.1 {
      1
    } else {
      0
    }
  })
  let bin_size = if n_bins > 0 { n_genes / n_bins } else { n_genes }
  let bins = Array(capacity=n_bins)
  for b in 0..<n_bins {
    let bin_genes = Array(capacity=bin_size)
    let start = b * bin_size
    let end = if b == n_bins - 1 { n_genes } else { (b + 1) * bin_size }
    for k in start..<end {
      if k < sorted_means.length() {
        bin_genes.push(sorted_means[k].0)
      }
    }
    bins.push(bin_genes)
  }
  // For each target gene, select control genes from same bin
  let ctrl_indices = Array(capacity=0)
  let mut rng_state = random_state
  for t in target_indices {
    // Find which bin this gene belongs to
    let t_mean = gene_means[t]
    let mut bin_idx = 0
    for b in 0..<n_bins {
      if bins[b].length() > 0 {
        let last_gene = bins[b][bins[b].length() - 1]
        if t_mean <= gene_means[last_gene] {
          bin_idx = b
          break
        }
        bin_idx = b
      }
    }
    // Select ctrl_size genes from same bin (deterministic with random_state)
    let bin = bins[bin_idx]
    for _c in 0..<ctrl_size {
      if bin.length() > 0 {
        rng_state = (rng_state * 1103515245 + 12345) % 2147483647
        let idx = rng_state % bin.length()
        if idx >= 0 {
          ctrl_indices.push(bin[idx])
        }
      }
    }
  }
  // Compute scores
  let scores = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let mut target_sum = 0.0
    for t in target_indices {
      target_sum = target_sum + x[i][t]
    }
    let target_mean = target_sum / target_indices.length().to_double()
    let mut ctrl_sum = 0.0
    for c in ctrl_indices {
      ctrl_sum = ctrl_sum + x[i][c]
    }
    let ctrl_mean = if ctrl_indices.length() > 0 {
      ctrl_sum / ctrl_indices.length().to_double()
    } else {
      0.0
    }
    scores.push(target_mean - ctrl_mean)
  }
  scores
}

///|
/// Advanced cell cycle scoring with S-phase and G2/M-phase gene sets.
/// Returns (S_score, G2M_score, phase) per cell.
pub fn scanpy_score_cell_cycle_advanced(
  x : Array[Array[Double]],
  gene_names : Array[String],
  s_genes : Array[String],
  g2m_genes : Array[String],
) -> Array[(Double, Double, String)] {
  let n_cells = x.length()
  let gene_set_s = Set([])
  for g in s_genes {
    gene_set_s.add(g)
  }
  let gene_set_g2m = Set([])
  for g in g2m_genes {
    gene_set_g2m.add(g)
  }
  let result = Array(capacity=n_cells)
  for i in 0..<n_cells {
    // S-phase score
    let mut s_sum = 0.0
    let mut s_count = 0
    for g in 0..<gene_names.length() {
      if gene_set_s.contains(gene_names[g]) {
        s_sum = s_sum + x[i][g]
        s_count = s_count + 1
      }
    }
    let s_score = if s_count > 0 { s_sum / s_count.to_double() } else { 0.0 }
    // G2M-phase score
    let mut g2m_sum = 0.0
    let mut g2m_count = 0
    for g in 0..<gene_names.length() {
      if gene_set_g2m.contains(gene_names[g]) {
        g2m_sum = g2m_sum + x[i][g]
        g2m_count = g2m_count + 1
      }
    }
    let g2m_score = if g2m_count > 0 {
      g2m_sum / g2m_count.to_double()
    } else {
      0.0
    }
    // Assign phase
    let phase = if s_score > g2m_score && s_score > 0.0 {
      "S"
    } else if g2m_score > s_score && g2m_score > 0.0 {
      "G2M"
    } else {
      "G1"
    }
    result.push((s_score, g2m_score, phase))
  }
  result
}

///|
/// Gene set enrichment analysis (GSEA-like).
/// Given a ranked gene list and gene sets, compute enrichment scores.
pub fn scanpy_gene_set_enrichment(
  ranked_genes : Array[String], // sorted by significance
  gene_sets : Map[String, Array[String]], // set_name -> gene list
  n_perm~ : Int32 = 100,
) -> Array[(String, Double, Double)] {
  let n = ranked_genes.length()
  let results = Array(capacity=0)
  // Build rank lookup
  let rank_map = Map([])
  for i in 0..<n {
    rank_map[ranked_genes[i]] = i
  }
  for set_name in gene_sets.keys() {
    let gene_set = gene_sets[set_name]
    let set_size = gene_set.length()
    if set_size == 0 {
      results.push((set_name, 0.0, 1.0))
      continue
    }
    // Compute running enrichment score (simplified KS statistic)
    let mut es = 0.0
    let mut max_es = 0.0
    let mut min_es = 0.0
    let hit_weight = 1.0 / set_size.to_double()
    let miss_weight = 1.0 / (n - set_size).to_double()
    for i in 0..<n {
      let in_set = gene_set.contains(ranked_genes[i])
      if in_set {
        es = es + hit_weight
      } else {
        es = es - miss_weight
      }
      if es > max_es {
        max_es = es
      }
      if es < min_es {
        min_es = es
      }
    }
    let score = if -min_es > max_es { min_es } else { max_es }
    // Simplified p-value: based on permutation
    let mut perm_count = 0
    let mut rng = 42
    for _p in 0..<n_perm {
      // Shuffle gene set membership
      let mut perm_es = 0.0
      let mut perm_max = 0.0
      for _i in 0..<n {
        rng = (rng * 1103515245 + 12345) % 2147483647
        let hit = rng % n < set_size
        if hit {
          perm_es = perm_es + hit_weight
        } else {
          perm_es = perm_es - miss_weight
        }
        if perm_es > perm_max {
          perm_max = perm_es
        }
      }
      if perm_max >= score {
        perm_count = perm_count + 1
      }
    }
    let pval = (perm_count + 1).to_double() / (n_perm + 1).to_double()
    results.push((set_name, score, pval))
  }
  results.sort_by(fn(a, b) {
    if a.2 < b.2 {
      -1
    } else if a.2 > b.2 {
      1
    } else {
      0
    }
  })
  results
}
