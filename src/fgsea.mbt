///|
/// FGSEA - Fast Gene Set Enrichment Analysis.
/// Inspired by Bioconductor fgsea package.
/// Implements GSEA with pre-ranked gene lists and permutation-based p-values.

///|
/// Represents a gene set with name and member genes.
pub struct GeneSet {
  name : String
  genes : Array[String]
}

///|
/// Create a new gene set.
pub fn GeneSet::new(name : String, genes : Array[String]) -> GeneSet {
  GeneSet::{ name, genes }
}

///|
/// Result of GSEA analysis for a single gene set.
pub struct FgseaResult {
  pathway : String
  pval : Double
  padj : Double
  es : Double
  nes : Double
  nMoreExtreme : Int
  size : Int
  leadingEdge : Array[String]
}

///|
/// Create a new FGSEA result.
pub fn FgseaResult::new(
  pathway : String,
  pval : Double,
  padj : Double,
  es : Double,
  nes : Double,
  nMoreExtreme : Int,
  size : Int,
  leadingEdge : Array[String],
) -> FgseaResult {
  FgseaResult::{ pathway, pval, padj, es, nes, nMoreExtreme, size, leadingEdge }
}

///|
/// Convert array of (gene, stat) pairs to a map for quick lookup.
fn gene_stat_map(ranks : Array[(String, Double)]) -> Map[String, Double] {
  Map::from_array(ranks)
}

///|
/// Calculate enrichment score (es) for a gene set using the running sum approach.
fn calc_enrichment_score(
  ranks : Array[(String, Double)],
  gene_set : Array[String],
  gene_map : Map[String, Double],
) -> (Double, Int) {
  let n = ranks.length()
  let gene_in_set : Array[Bool] = []
  let mut set_size = 0
  for i = 0; i < n; i = i + 1 {
    let (gene, _) = ranks[i]
    let mut found = false
    for g in gene_set {
      if g == gene {
        found = true
        set_size = set_size + 1
        break
      }
    }
    gene_in_set.push(found)
  }

  if set_size == 0 || set_size == n {
    return (0.0, 0)
  }

  let mut sum_correl = 0.0
  for i = 0; i < n; i = i + 1 {
    if gene_in_set[i] {
      let (_, stat) = ranks[i]
      sum_correl = sum_correl + stat.abs()
    }
  }

  if sum_correl == 0.0 {
    return (0.0, 0)
  }

  let p_hit = 1.0 / sum_correl
  let p_miss = 1.0 / (n - set_size).to_double()

  let mut max_es = 0.0
  let mut max_pos = 0
  let mut running_sum = 0.0

  for i = 0; i < n; i = i + 1 {
    if gene_in_set[i] {
      let (_, stat) = ranks[i]
      running_sum = running_sum + stat.abs() * p_hit
    } else {
      running_sum = running_sum - p_miss
    }
    if running_sum.abs() > max_es.abs() {
      max_es = running_sum
      max_pos = i
    }
  }

  (max_es, max_pos)
}

///|
/// Extract leading edge genes from ranked list and gene set.
fn extract_leading_edge(
  ranks : Array[(String, Double)],
  gene_set : Array[String],
  peak_pos : Int,
  es_positive : Bool,
) -> Array[String] {
  let result : Array[String] = []
  let start = if es_positive { 0 } else { peak_pos }
  let end = if es_positive { peak_pos + 1 } else { ranks.length() }

  for i = start; i < end; i = i + 1 {
    let (gene, _) = ranks[i]
    for g in gene_set {
      if g == gene {
        result.push(gene)
        break
      }
    }
  }
  result
}

///|
/// Perform a single permutation by shuffling gene labels.
fn permute_es(
  ranks : Array[(String, Double)],
  gene_set : Array[String],
  gene_map : Map[String, Double],
  seed : Int,
) -> Double {
  let n = ranks.length()
  let perm_ranks : Array[(String, Double)] = []
  let indices : Array[Int] = []
  for i = 0; i < n; i = i + 1 {
    indices.push(i)
  }

  let mut s = seed
  for i = n - 1; i > 0; i = i - 1 {
    s = (s * 1103515245 + 12345) % 2147483647
    if s < 0 {
      s = s + 2147483647
    }
    let j = s % (i + 1)
    let tmp = indices[i]
    indices[i] = indices[j]
    indices[j] = tmp
  }

  for i = 0; i < n; i = i + 1 {
    let (gene, _) = ranks[indices[i]]
    let (_, stat) = ranks[i]
    perm_ranks.push((gene, stat))
  }

  let perm_map = gene_stat_map(perm_ranks)
  let (es, _) = calc_enrichment_score(perm_ranks, gene_set, perm_map)
  es
}

///|
/// Run FGSEA analysis on multiple gene sets with permutation-based p-values.
pub fn fgsea(
  pathways : Array[GeneSet],
  stats : Array[(String, Double)],
  nperm? : Int = 100,
  min_size? : Int = 1,
  max_size? : Int = 500,
) -> Array[FgseaResult] {
  let gene_map = gene_stat_map(stats)
  let results : Array[FgseaResult] = []

  for pathway in pathways {
    let gene_set = pathway.genes
    let set_size = gene_set.length()

    if set_size < min_size || set_size > max_size {
      continue
    }

    let (es, peak_pos) = calc_enrichment_score(stats, gene_set, gene_map)
    let es_positive = es > 0.0

    let mut n_more_extreme = 0
    let mut perm_es_sum_pos = 0.0
    let mut perm_es_sum_neg = 0.0
    let mut perm_count_pos = 0
    let mut perm_count_neg = 0
    let perm_es_pos_list : Array[Double] = []
    let perm_es_neg_list : Array[Double] = []

    for p = 0; p < nperm; p = p + 1 {
      let seed = p * 10007 + 12345
      let perm_es = permute_es(stats, gene_set, gene_map, seed)
      if perm_es >= 0.0 {
        perm_es_sum_pos = perm_es_sum_pos + perm_es
        perm_count_pos = perm_count_pos + 1
        perm_es_pos_list.push(perm_es)
      } else {
        perm_es_sum_neg = perm_es_sum_neg + perm_es
        perm_count_neg = perm_count_neg + 1
        perm_es_neg_list.push(perm_es)
      }
      if (es >= 0.0 && perm_es >= es) || (es < 0.0 && perm_es <= es) {
        n_more_extreme = n_more_extreme + 1
      }
    }

    let pval = (n_more_extreme + 1).to_double() / (nperm + 1).to_double()

    let nes = if es >= 0.0 {
      if perm_count_pos > 0 {
        let mean_pos = perm_es_sum_pos / perm_count_pos.to_double()
        if mean_pos > 0.0 {
          es / mean_pos
        } else {
          0.0
        }
      } else {
        0.0
      }
    } else if perm_count_neg > 0 {
      let mean_neg = perm_es_sum_neg / perm_count_neg.to_double()
      if mean_neg < 0.0 {
        es / mean_neg.abs()
      } else {
        0.0
      }
    } else {
      0.0
    }

    let leading_edge = extract_leading_edge(
      stats, gene_set, peak_pos, es_positive,
    )

    let result = FgseaResult::new(
      pathway.name,
      pval,
      pval,
      es,
      nes,
      n_more_extreme,
      set_size,
      leading_edge,
    )
    results.push(result)
  }

  results
}

///|
/// Rank genes by statistic value (descending order).
pub fn rank_genes(
  genes : Array[String],
  stats : Array[Double],
) -> Array[(String, Double)] {
  let n = genes.length()
  let pairs : Array[(String, Double)] = []
  for i = 0; i < n; i = i + 1 {
    pairs.push((genes[i], stats[i]))
  }

  for i = 0; i < n - 1; i = i + 1 {
    for j = 0; j < n - 1 - i; j = j + 1 {
      let (_, s1) = pairs[j]
      let (_, s2) = pairs[j + 1]
      if s1 < s2 {
        let tmp = pairs[j]
        pairs[j] = pairs[j + 1]
        pairs[j + 1] = tmp
      }
    }
  }

  pairs
}

///|
/// Filter FGSEA results by adjusted p-value threshold.
pub fn fgsea_filter(
  results : Array[FgseaResult],
  padj_cutoff? : Double = 0.05,
) -> Array[FgseaResult] {
  let filtered : Array[FgseaResult] = []
  for r in results {
    if r.padj <= padj_cutoff {
      filtered.push(r)
    }
  }
  filtered
}

///|
/// Get top up-regulated gene sets (positive nes).
pub fn fgsea_top_up(
  results : Array[FgseaResult],
  n? : Int = 10,
) -> Array[FgseaResult] {
  let up : Array[FgseaResult] = []
  for r in results {
    if r.nes > 0.0 {
      up.push(r)
    }
  }

  for i = 0; i < up.length() - 1; i = i + 1 {
    for j = 0; j < up.length() - 1 - i; j = j + 1 {
      if up[j].nes < up[j + 1].nes {
        let tmp = up[j]
        up[j] = up[j + 1]
        up[j + 1] = tmp
      }
    }
  }

  let top : Array[FgseaResult] = []
  for i = 0; i < up.length() && i < n; i = i + 1 {
    top.push(up[i])
  }
  top
}

///|
/// Get top down-regulated gene sets (negative nes).
pub fn fgsea_top_down(
  results : Array[FgseaResult],
  n? : Int = 10,
) -> Array[FgseaResult] {
  let down : Array[FgseaResult] = []
  for r in results {
    if r.nes < 0.0 {
      down.push(r)
    }
  }

  for i = 0; i < down.length() - 1; i = i + 1 {
    for j = 0; j < down.length() - 1 - i; j = j + 1 {
      if down[j].nes > down[j + 1].nes {
        let tmp = down[j]
        down[j] = down[j + 1]
        down[j + 1] = tmp
      }
    }
  }

  let top : Array[FgseaResult] = []
  for i = 0; i < down.length() && i < n; i = i + 1 {
    top.push(down[i])
  }
  top
}

///|
/// Collapse (average) redundant gene sets based on leading edge similarity.
pub fn fgsea_collapse_pathways(
  results : Array[FgseaResult],
  pathways : Array[GeneSet],
  threshold? : Double = 0.5,
) -> Array[Array[Int]] {
  let n = results.length()
  let clusters : Array[Array[Int]] = []
  let assigned : Array[Bool] = []
  for i = 0; i < n; i = i + 1 {
    assigned.push(false)
  }

  for i = 0; i < n; i = i + 1 {
    if assigned[i] {
      continue
    }
    assigned[i] = true
    let cluster : Array[Int] = [i]
    let le_i = results[i].leadingEdge

    for j = i + 1; j < n; j = j + 1 {
      if assigned[j] {
        continue
      }
      let le_j = results[j].leadingEdge
      let mut overlap = 0
      for a in le_i {
        for b in le_j {
          if a == b {
            overlap = overlap + 1
            break
          }
        }
      }
      let min_size = if le_i.length() < le_j.length() {
        le_i.length()
      } else {
        le_j.length()
      }
      if min_size > 0 && overlap.to_double() / min_size.to_double() >= threshold {
        assigned[j] = true
        cluster.push(j)
      }
    }
    clusters.push(cluster)
  }

  clusters
}

///|
/// Sample gene sets for testing/demo purposes.
pub fn sample_gene_sets() -> Array[GeneSet] {
  [
    GeneSet::new("CELL_CYCLE", [
      "CDK1", "CCNA2", "CCNB1", "CDC20", "CDKN1A", "TP53", "RB1", "E2F1", "MYC",
      "PCNA", "MCM2", "MCM3", "ORC1", "CDC6", "CDT1",
    ]),
    GeneSet::new("APOPTOSIS", [
      "BCL2", "BAX", "BAK1", "CASP3", "CASP8", "CASP9", "FAS", "TNFRSF1A", "TP53",
      "BAD", "BID", "BIM", "XIAP", "DIABLO", "APAF1",
    ]),
    GeneSet::new("METABOLISM", [
      "GAPDH", "LDHA", "PFKL", "PKM", "HK2", "PGK1", "ENO1", "ALDOA", "GLUT1", "IDH1",
      "SDHA", "MDH2", "CS", "ACLY", "FASN",
    ]),
    GeneSet::new("IMMUNE_ResPONSE", [
      "IL6", "TNF", "IL1B", "CXCL8", "CCL2", "IL10", "IFNG", "IL4", "IL12A", "STAT1",
      "NFKB1", "RELA", "IRF3", "TLR4", "CD40",
    ]),
    GeneSet::new("DNA_REPAIR", [
      "ATM", "ATR", "BRCA1", "BRCA2", "RAD51", "TP53", "CHEK1", "CHEK2", "MSH2",
      "MLH1", "PMS2", "MSH6", "XPA", "ERCC1", "PCNA",
    ]),
  ]
}

///|
/// Sample gene statistics for testing/demo purposes.
pub fn sample_gene_stats() -> Array[(String, Double)] {
  let genes = [
    "CDK1", "CCNA2", "CCNB1", "CDC20", "CDKN1A", "TP53", "RB1", "E2F1", "MYC", "PCNA",
    "MCM2", "MCM3", "ORC1", "CDC6", "CDT1", "BCL2", "BAX", "BAK1", "CASP3", "CASP8",
    "CASP9", "FAS", "TNFRSF1A", "BAD", "BID", "BIM", "XIAP", "DIABLO", "APAF1", "GAPDH",
    "LDHA", "PFKL", "PKM", "HK2", "PGK1", "ENO1", "ALDOA", "GLUT1", "IDH1", "SDHA",
    "MDH2", "CS", "ACLY", "FASN", "IL6", "TNF", "IL1B", "CXCL8", "CCL2", "IL10",
    "IFNG", "IL4", "IL12A", "STAT1", "NFKB1", "RELA", "IRF3", "TLR4", "CD40", "ATM",
    "ATR", "BRCA1", "BRCA2", "RAD51", "CHEK1", "CHEK2", "MSH2", "MLH1", "PMS2", "MSH6",
    "XPA", "ERCC1", "ACTB", "GAPDH", "TUBB", "HPRT1", "B2M", "YWHAZ",
  ]
  let stats = [
    2.5, 2.3, 2.1, 1.9, -1.5, -2.0, -1.8, 1.6, 1.7, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9,
    -1.2, 1.8, 1.5, 2.0, 1.7, 1.6, -0.8, -0.6, 1.3, 1.1, 1.4, -1.0, 0.9, 0.8, 0.5,
    0.3, 0.4, 0.2, 0.6, 0.7, -0.3, -0.2, 0.1, -0.4, -0.5, -0.1, 0.0, -0.6, -0.7,
    2.2, 1.9, 1.8, 1.5, 1.3, -0.9, 1.6, -1.1, 1.4, 1.7, 1.2, 1.0, 0.8, 1.1, 0.9,
    -1.3, -1.4, -1.6, -1.7, -1.5, -0.5, -0.7, -0.8, -0.6, -0.4, -0.3, -0.2, -0.1,
    0.05, 0.03, 0.04, 0.02, 0.06, 0.01,
  ]
  rank_genes(genes, stats)
}
