///|
/// Bioconductor HTSFilter - RNA-seq count filtering based on CPM (Counts Per Million) threshold.
/// Filters low-expressed genes prior to differential expression analysis.
///
/// Reference: Bioconductor HTSFilter package
/// Rows = genes, Columns = samples.
/// Groups are used to define the minimum number of samples that must pass the threshold.

///|
/// HTSFilter configuration and results.
pub struct HTSFilterResult {
  cpm_threshold : Double
  min_samples_per_group : Int
  groups : Array[Int]
  cpm_matrix : Array[Array[Double]]
  keep : Array[Bool]
  n_genes_input : Int
  n_genes_kept : Int
  n_genes_removed : Int
}

///|
/// Compute CPM (Counts Per Million) matrix from raw counts and library sizes.
/// counts: genes × samples
/// Returns: genes × samples CPM values
pub fn hts_filter_cpm(counts : Array[Array[Double]], library_sizes : Array[Double]) -> Array[Array[Double]] {
  let n_genes = counts.length()
  if n_genes == 0 {
    return []
  }
  let n_samples = counts[0].length()
  let result : Array[Array[Double]] = []
  let mut g = 0
  while g < n_genes {
    let row : Array[Double] = []
    let mut s = 0
    while s < n_samples {
      let lib = library_sizes[s]
      if lib > 0.0 {
        let cpm = counts[g][s] / lib * 1_000_000.0
        row.push(cpm)
      } else {
        row.push(0.0)
      }
      s = s + 1
    }
    result.push(row)
    g = g + 1
  }
  result
}

///|
/// Compute library sizes (column sums) from a count matrix.
pub fn hts_filter_library_sizes(counts : Array[Array[Double]]) -> Array[Double] {
  if counts.length() == 0 {
    return []
  }
  let n_samples = counts[0].length()
  let lib : Array[Double] = []
  let mut s = 0
  while s < n_samples {
    let mut total = 0.0
    let mut g = 0
    while g < counts.length() {
      total = total + counts[g][s]
      g = g + 1
    }
    lib.push(total)
    s = s + 1
  }
  lib
}

///|
/// Filter a count matrix by CPM threshold.
/// For each gene, keep it if at least min_samples_per_group samples in any group
/// have CPM >= cpm_threshold.
///
/// Parameters:
///   counts : genes × samples raw count matrix
///   groups : group assignment per sample (e.g., [0,0,0,1,1,1])
///   cpm_threshold : minimum CPM to consider a gene expressed
///   min_samples_per_group : minimum number of samples per group that must pass
///
/// Returns HTSFilterResult with keep mask and statistics.
pub fn hts_filter(
  counts : Array[Array[Double]],
  groups : Array[Int],
  cpm_threshold : Double,
  min_samples_per_group : Int,
) -> HTSFilterResult {
  let n_genes = counts.length()
  if n_genes == 0 {
    return HTSFilterResult::{
      cpm_threshold: cpm_threshold,
      min_samples_per_group: min_samples_per_group,
      groups: groups,
      cpm_matrix: [],
      keep: [],
      n_genes_input: 0,
      n_genes_kept: 0,
      n_genes_removed: 0,
    }
  }
  let n_samples = counts[0].length()
  let lib_sizes = hts_filter_library_sizes(counts)
  let cpm_mat = hts_filter_cpm(counts, lib_sizes)

  // Find unique groups
  let unique_groups : Array[Int] = []
  let mut i = 0
  while i < groups.length() {
    let g = groups[i]
    let mut found = false
    let mut j = 0
    while j < unique_groups.length() {
      if unique_groups[j] == g { found = true; break }
      j = j + 1
    }
    if !found { unique_groups.push(g) }
    i = i + 1
  }

  let keep : Array[Bool] = []
  let mut kept = 0
  let mut g = 0
  while g < n_genes {
    let mut gene_kept = false
    let mut ug = 0
    while ug < unique_groups.length() {
      let group_id = unique_groups[ug]
      let mut count_passing = 0
      let mut s = 0
      while s < n_samples {
        if groups[s] == group_id && cpm_mat[g][s] >= cpm_threshold {
          count_passing = count_passing + 1
        }
        s = s + 1
      }
      if count_passing >= min_samples_per_group {
        gene_kept = true
        break
      }
      ug = ug + 1
    }
    keep.push(gene_kept)
    if gene_kept { kept = kept + 1 }
    g = g + 1
  }

  HTSFilterResult::{
    cpm_threshold: cpm_threshold,
    min_samples_per_group: min_samples_per_group,
    groups: groups,
    cpm_matrix: cpm_mat,
    keep: keep,
    n_genes_input: n_genes,
    n_genes_kept: kept,
    n_genes_removed: n_genes - kept,
  }
}

///|
/// Apply a keep mask to filter a count matrix.
/// Returns only the rows (genes) where keep[i] is true.
pub fn hts_filter_apply(counts : Array[Array[Double]], keep : Array[Bool]) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  let mut i = 0
  while i < counts.length() {
    if keep[i] {
      result.push(counts[i])
    }
    i = i + 1
  }
  result
}

///|
/// Apply keep mask to a gene name list.
pub fn hts_filter_apply_names(gene_names : Array[String], keep : Array[Bool]) -> Array[String] {
  let result : Array[String] = []
  let mut i = 0
  while i < gene_names.length() {
    if keep[i] {
      result.push(gene_names[i])
    }
    i = i + 1
  }
  result
}

///|
/// Summarize the filtering result as a string.
pub fn hts_filter_summary(result : HTSFilterResult) -> String {
  "HTSFilter Summary:\n" +
    "  CPM threshold: " + result.cpm_threshold.to_string() + "\n" +
    "  Min samples per group: " + result.min_samples_per_group.to_string() + "\n" +
    "  Input genes: " + result.n_genes_input.to_string() + "\n" +
    "  Genes kept: " + result.n_genes_kept.to_string() + "\n" +
    "  Genes removed: " + result.n_genes_removed.to_string() + "\n" +
    "  Retention rate: " +
    (if result.n_genes_input > 0 {
      (result.n_genes_kept.to_double() / result.n_genes_input.to_double() * 100.0).to_string()
    } else { "0.0" }) + "%"
}

///|
/// Generate sample RNA-seq count data for testing.
/// Returns tuple of (counts, groups, gene_names).
pub fn hts_filter_sample_data() -> (Array[Array[Double]], Array[Int], Array[String]) {
  // 20 genes x 6 samples (3 control + 3 treatment)
  let n_genes = 20
  let n_samples = 6
  let counts : Array[Array[Double]] = []
  let gene_names : Array[String] = []
  let groups = [0, 0, 0, 1, 1, 1]

  let mut g = 0
  while g < n_genes {
    let row : Array[Double] = []
    let mut s = 0
    while s < n_samples {
      // First 10 genes: high expression (~1000-5000 CPM)
      // Next 5 genes: moderate expression (~50-200 CPM)
      // Last 5 genes: low expression (~1-10 CPM)
      let base = if g < 10 { 2000.0 } else if g < 15 { 100.0 } else { 5.0 }
      let noise = (((g * 7 + s * 13) % 100).to_double() / 100.0) * base * 0.5
      let val = base + noise
      row.push(val)
      s = s + 1
    }
    counts.push(row)
    gene_names.push("Gene" + (g + 1).to_string())
    g = g + 1
  }
  (counts, groups, gene_names)
}

///|
/// Compute CPM for a single count value given library size.
pub fn hts_filter_single_cpm(count : Double, library_size : Double) -> Double {
  if library_size > 0.0 {
    count / library_size * 1_000_000.0
  } else {
    0.0
  }
}

///|
/// Get the keep mask from a filtering result.
pub fn hts_filter_get_keep(result : HTSFilterResult) -> Array[Bool] {
  result.keep
}

///|
/// Calculate the retention rate (fraction of genes kept).
pub fn hts_filter_retention_rate(result : HTSFilterResult) -> Double {
  if result.n_genes_input > 0 {
    result.n_genes_kept.to_double() / result.n_genes_input.to_double()
  } else {
    0.0
  }
}
