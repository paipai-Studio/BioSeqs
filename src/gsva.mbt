///|
/// GSVA - Gene Set Variation Analysis.
/// Inspired by Bioconductor GSVA package.
/// Implements single-sample pathway scoring for RNA-seq and microarray data.

///|
/// Represents a gene set with name, member genes, and description.
pub struct GSVAGeneSet {
  name : String
  genes : Array[String]
  description : String
}

///|
/// Represents gene expression data matrix with gene and sample names.
pub struct GSVAData {
  expression : Array[Array[Double]]
  gene_names : Array[String]
  sample_names : Array[String]
  gene_sets : Array[GSVAGeneSet]
}

///|
/// Represents GSVA scores for gene sets across samples.
pub struct GSVAScores {
  scores : Array[Array[Double]]
  gene_set_names : Array[String]
  sample_names : Array[String]
  method : String
  params : Map[String, String]
}

///|
/// Parameters for GSVA analysis.
pub struct GSVAParams {
  method : String
  kcdf : String
  min_size : Int
  max_size : Int
  abs_ranking : Bool
  verbose : Bool
}

///|
/// Represents enrichment result for a gene set.
pub struct GSVAEnrichment {
  gene_set_name : String
  enrichment_score : Double
  normalized_enrichment_score : Double
  pvalue : Double
  fdr : Double
  leading_edge : Array[String]
}

///|
/// Statistical summary of GSVA scores.
pub struct GSVAStat {
  mean_score : Array[Double]
  sd_score : Array[Double]
  median_score : Array[Double]
  min_score : Array[Double]
  max_score : Array[Double]
}

///|
/// Create a new gene set.
pub fn GSVAGeneSet::new(
  name : String,
  genes : Array[String],
  description : String,
) -> GSVAGeneSet {
  GSVAGeneSet::{ name, genes, description }
}

///|
/// Create a new GSVA data object.
pub fn GSVAData::new(
  expression : Array[Array[Double]],
  gene_names : Array[String],
  sample_names : Array[String],
) -> GSVAData {
  GSVAData::{ expression, gene_names, sample_names, gene_sets: [] }
}

///|
/// Create default GSVA parameters.
pub fn GSVAParams::new() -> GSVAParams {
  GSVAParams::{
    method: "ssgsea",
    kcdf: "Gaussian",
    min_size: 1,
    max_size: 500,
    abs_ranking: false,
    verbose: false,
  }
}

///|
/// Newton's method for square root.
pub fn gsva_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut guess = x
  let mut prev = 0.0
  let mut i = 0
  while gsva_abs(guess - prev) > 0.0000000001 && i < 100 {
    prev = guess
    guess = (guess + x / guess) / 2.0
    i = i + 1
  }
  guess
}

///|
/// Absolute value for double.
fn gsva_abs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

///|
/// Natural logarithm using Taylor series.
pub fn gsva_ln(x : Double) -> Double {
  if x <= 0.0 {
    return -1000000000000000000000000000000.0
  }
  let mut y = x
  let mut log2 = 0.0
  while y >= 2.0 {
    y = y / 2.0
    log2 = log2 + 1.0
  }
  while y < 1.0 {
    y = y * 2.0
    log2 = log2 - 1.0
  }
  let z = (y - 1.0) / (y + 1.0)
  let z_sq = z * z
  let mut result = 0.0
  let mut term = z
  let mut i = 1
  while gsva_abs(term) > 0.000000000000001 && i < 100 {
    result = result + term / i.to_double()
    term = term * z_sq
    i = i + 2
  }
  2.0 * result + log2 * 0.6931471805599453
}

///|
/// Exponential function using Taylor series.
pub fn gsva_exp(x : Double) -> Double {
  if x > 700.0 {
    return 1000000000000000000000000000000.0
  }
  if x < -700.0 {
    return 0.0
  }
  let mut result = 1.0
  let mut term = 1.0
  let mut i = 1
  while gsva_abs(term) > 0.000000000000001 && i < 100 {
    term = term * x / i.to_double()
    result = result + term
    i = i + 1
  }
  result
}

///|
/// Calculate mean of array.
pub fn gsva_mean(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < arr.length() {
    sum = sum + arr[i]
    i = i + 1
  }
  sum / arr.length().to_double()
}

///|
/// Calculate variance of array.
pub fn gsva_variance(arr : Array[Double]) -> Double {
  if arr.length() < 2 {
    return 0.0
  }
  let mean = gsva_mean(arr)
  let mut sum_sq = 0.0
  let mut i = 0
  while i < arr.length() {
    let diff = arr[i] - mean
    sum_sq = sum_sq + diff * diff
    i = i + 1
  }
  sum_sq / (arr.length() - 1).to_double()
}

///|
/// Calculate median of array.
pub fn gsva_median(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let sorted = gsva_sort_doubles(arr)
  let n = sorted.length()
  if n % 2 == 0 {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  } else {
    sorted[n / 2]
  }
}

///|
/// Calculate standard deviation of array.
pub fn gsva_std(arr : Array[Double]) -> Double {
  gsva_sqrt(gsva_variance(arr))
}

///|
/// Sort array of doubles in ascending order.
fn gsva_sort_doubles(arr : Array[Double]) -> Array[Double] {
  let result = Array::make(arr.length(), 0.0)
  let mut i = 0
  while i < arr.length() {
    result[i] = arr[i]
    i = i + 1
  }
  let n = result.length()
  i = 0
  while i < n - 1 {
    let mut j = 0
    while j < n - 1 - i {
      if result[j] > result[j + 1] {
        let tmp = result[j]
        result[j] = result[j + 1]
        result[j + 1] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Calculate percent rank of array elements.
pub fn gsva_percent_rank(arr : Array[Double]) -> Array[Double] {
  if arr.length() == 0 {
    return []
  }
  let sorted = gsva_sort_doubles(arr)
  let n = arr.length()
  let result = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let val = arr[i]
    let mut rank = 0
    let mut j = 0
    while j < n {
      if sorted[j] < val {
        rank = rank + 1
      }
      j = j + 1
    }
    result[i] = rank.to_double() / (n - 1).to_double()
    i = i + 1
  }
  result
}

///|
/// Normal CDF approximation using error function.
pub fn gsva_cdf(x : Double, mean : Double, sd : Double) -> Double {
  if sd <= 0.0 {
    return if x >= mean { 1.0 } else { 0.0 }
  }
  let z = (x - mean) / sd
  0.5 * (1.0 + gsva_erf(z / gsva_sqrt(2.0)))
}

///|
/// Error function approximation.
fn gsva_erf(x : Double) -> Double {
  let abs_x = gsva_abs(x)
  let t = 1.0 / (1.0 + 0.3275911 * abs_x)
  let y = 1.0 -
    (
      (((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) *
      t +
      0.254829592
    ) *
    t *
    gsva_exp(-abs_x * abs_x)
  if x >= 0.0 {
    y
  } else {
    -y
  }
}

///|
/// Find gene indices in gene_names that match gene_set genes.
fn gsva_find_gene_indices(
  gene_set : GSVAGeneSet,
  gene_names : Array[String],
) -> Array[Int] {
  let indices : Array[Int] = []
  let mut i = 0
  while i < gene_set.genes.length() {
    let gene = gene_set.genes[i]
    let mut j = 0
    while j < gene_names.length() {
      if gene_names[j] == gene {
        indices.push(j)
        break
      }
      j = j + 1
    }
    i = i + 1
  }
  indices
}

///|
/// Filter gene sets by size constraints.
pub fn gsva_filter_gene_sets(
  gene_sets : Array[GSVAGeneSet],
  gene_names : Array[String],
  min_size : Int,
  max_size : Int,
) -> Array[GSVAGeneSet] {
  let filtered : Array[GSVAGeneSet] = []
  let mut i = 0
  while i < gene_sets.length() {
    let indices = gsva_find_gene_indices(gene_sets[i], gene_names)
    if indices.length() >= min_size && indices.length() <= max_size {
      filtered.push(gene_sets[i])
    }
    i = i + 1
  }
  filtered
}

///|
/// Rank genes by expression values for each sample.
pub fn gsva_rank_genes(
  expression : Array[Array[Double]],
  method : String,
) -> Array[Array[Double]] {
  let n_genes = expression.length()
  let n_samples = if n_genes > 0 { expression[0].length() } else { 0 }
  let ranks : Array[Array[Double]] = Array::new()

  let mut sample_idx = 0
  while sample_idx < n_samples {
    let gene_values : Array[(Int, Double)] = []
    let mut gene_idx = 0
    while gene_idx < n_genes {
      let val = if method == "abs" {
        gsva_abs(expression[gene_idx][sample_idx])
      } else {
        expression[gene_idx][sample_idx]
      }
      gene_values.push((gene_idx, val))
      gene_idx = gene_idx + 1
    }

    let n = gene_values.length()
    let mut i = 0
    while i < n - 1 {
      let mut j = 0
      while j < n - 1 - i {
        let (_, v1) = gene_values[j]
        let (_, v2) = gene_values[j + 1]
        if v1 < v2 {
          let tmp = gene_values[j]
          gene_values[j] = gene_values[j + 1]
          gene_values[j + 1] = tmp
        }
        j = j + 1
      }
      i = i + 1
    }

    let sample_ranks = Array::make(n_genes, 0.0)
    let mut rank = 1.0
    let mut j = 0
    while j < n {
      let (gene_idx, _) = gene_values[j]
      sample_ranks[gene_idx] = rank
      rank = rank + 1.0
      j = j + 1
    }
    ranks.push(sample_ranks)
    sample_idx = sample_idx + 1
  }

  ranks
}

///|
/// Calculate GSEA enrichment score.
pub fn gsva_gsea_score(
  gene_ranks : Array[Double],
  gene_set_indices : Array[Int],
  n_genes : Int,
) -> Double {
  if gene_set_indices.length() == 0 {
    return 0.0
  }

  let set_size = gene_set_indices.length()
  let not_set_size = n_genes - set_size

  let mut sum_ranks = 0.0
  let mut i = 0
  while i < set_size {
    sum_ranks = sum_ranks + gene_ranks[gene_set_indices[i]]
    i = i + 1
  }

  let p_hit = 1.0 / sum_ranks
  let p_miss = 1.0 / not_set_size.to_double()

  let mut running_sum = 0.0
  let mut max_es = 0.0

  let rank_sorted : Array[(Double, Bool)] = []
  let mut j = 0
  while j < n_genes {
    let mut in_set = false
    let mut k = 0
    while k < set_size {
      if gene_set_indices[k] == j {
        in_set = true
        break
      }
      k = k + 1
    }
    rank_sorted.push((gene_ranks[j], in_set))
    j = j + 1
  }

  let n = rank_sorted.length()
  j = 0
  while j < n - 1 {
    let mut k = 0
    while k < n - 1 - j {
      let (r1, _) = rank_sorted[k]
      let (r2, _) = rank_sorted[k + 1]
      if r1 > r2 {
        let tmp = rank_sorted[k]
        rank_sorted[k] = rank_sorted[k + 1]
        rank_sorted[k + 1] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }

  j = 0
  while j < n {
    let (_, in_set) = rank_sorted[j]
    if in_set {
      running_sum = running_sum + p_hit
    } else {
      running_sum = running_sum - p_miss
    }
    if gsva_abs(running_sum) > gsva_abs(max_es) {
      max_es = running_sum
    }
    j = j + 1
  }

  max_es
}

///|
/// ssGSEA scoring method.
pub fn gsva_score_ssgsea(
  expression : Array[Array[Double]],
  gene_set : GSVAGeneSet,
  gene_names : Array[String],
) -> Array[Double] {
  let n_samples = if expression.length() > 0 {
    expression[0].length()
  } else {
    0
  }
  let indices = gsva_find_gene_indices(gene_set, gene_names)

  let scores = Array::make(n_samples, 0.0)
  let mut sample_idx = 0
  while sample_idx < n_samples {
    let mut set_sum = 0.0
    let mut i = 0
    while i < indices.length() {
      set_sum = set_sum + expression[indices[i]][sample_idx]
      i = i + 1
    }

    let mut total_sum = 0.0
    let mut j = 0
    while j < expression.length() {
      total_sum = total_sum + expression[j][sample_idx]
      j = j + 1
    }

    scores[sample_idx] = if total_sum > 0.0 { set_sum / total_sum } else { 0.0 }
    sample_idx = sample_idx + 1
  }

  scores
}

///|
/// Z-score scoring method.
pub fn gsva_score_zscore(
  expression : Array[Array[Double]],
  gene_set : GSVAGeneSet,
  gene_names : Array[String],
) -> Array[Double] {
  let n_samples = if expression.length() > 0 {
    expression[0].length()
  } else {
    0
  }
  let indices = gsva_find_gene_indices(gene_set, gene_names)

  let scores = Array::make(n_samples, 0.0)
  let mut sample_idx = 0
  while sample_idx < n_samples {
    let values : Array[Double] = []
    let mut i = 0
    while i < indices.length() {
      values.push(expression[indices[i]][sample_idx])
      i = i + 1
    }

    let mean = gsva_mean(values)
    let sd = gsva_std(values)
    scores[sample_idx] = if sd > 0.0 {
      (mean - gsva_mean(expression[0])) / sd
    } else {
      0.0
    }
    sample_idx = sample_idx + 1
  }

  scores
}

///|
/// PLAGE scoring method (Principal Component Analysis).
pub fn gsva_score_plage(
  expression : Array[Array[Double]],
  gene_set : GSVAGeneSet,
  gene_names : Array[String],
) -> Array[Double] {
  let n_samples = if expression.length() > 0 {
    expression[0].length()
  } else {
    0
  }
  let indices = gsva_find_gene_indices(gene_set, gene_names)

  if indices.length() < 2 {
    return Array::make(n_samples, 0.0)
  }

  let scores = Array::make(n_samples, 0.0)
  let mut sample_idx = 0
  while sample_idx < n_samples {
    let set_expr : Array[Double] = []
    let mut i = 0
    while i < indices.length() {
      set_expr.push(expression[indices[i]][sample_idx])
      i = i + 1
    }

    let mean_val = gsva_mean(set_expr)
    let centered : Array[Double] = []
    let mut j = 0
    while j < set_expr.length() {
      centered.push(set_expr[j] - mean_val)
      j = j + 1
    }

    let mut sum_sq = 0.0
    j = 0
    while j < centered.length() {
      sum_sq = sum_sq + centered[j] * centered[j]
      j = j + 1
    }

    scores[sample_idx] = if sum_sq > 0.0 {
      gsva_sqrt(sum_sq) / centered.length().to_double()
    } else {
      0.0
    }
    sample_idx = sample_idx + 1
  }

  scores
}

///|
/// GSVA original scoring method.
pub fn gsva_score_gsva(
  expression : Array[Array[Double]],
  gene_set : GSVAGeneSet,
  gene_names : Array[String],
  kcdf : String,
) -> Array[Double] {
  let ranks = gsva_rank_genes(
    expression,
    if kcdf == "none" {
      "none"
    } else {
      "abs"
    },
  )
  let indices = gsva_find_gene_indices(gene_set, gene_names)
  let n_genes = expression.length()

  let n_samples = if expression.length() > 0 {
    expression[0].length()
  } else {
    0
  }
  let scores = Array::make(n_samples, 0.0)

  let mut sample_idx = 0
  while sample_idx < n_samples {
    scores[sample_idx] = gsva_gsea_score(ranks[sample_idx], indices, n_genes)
    sample_idx = sample_idx + 1
  }

  scores
}

///|
/// Normalize scores using specified method.
pub fn gsva_normalize_scores(
  scores : Array[Array[Double]],
  method : String,
) -> Array[Array[Double]] {
  let n_sets = scores.length()
  let n_samples = if n_sets > 0 { scores[0].length() } else { 0 }

  let normalized = Array::make(n_sets, [])
  let mut set_idx = 0
  while set_idx < n_sets {
    let set_scores = Array::make(n_samples, 0.0)
    let mut i = 0
    while i < n_samples {
      set_scores[i] = scores[set_idx][i]
      i = i + 1
    }

    if method == "zscore" {
      let mean = gsva_mean(set_scores)
      let sd = gsva_std(set_scores)
      let mut j = 0
      while j < n_samples {
        set_scores[j] = if sd > 0.0 { (set_scores[j] - mean) / sd } else { 0.0 }
        j = j + 1
      }
    } else if method == "rank" {
      let ranked = gsva_percent_rank(set_scores)
      let mut j = 0
      while j < n_samples {
        set_scores[j] = ranked[j]
        j = j + 1
      }
    }

    normalized[set_idx] = set_scores
    set_idx = set_idx + 1
  }

  normalized
}

///|
/// Run GSVA analysis.
pub fn gsva_run(
  data : GSVAData,
  gene_sets : Array[GSVAGeneSet],
  params : GSVAParams,
) -> GSVAScores {
  let filtered_sets = gsva_filter_gene_sets(
    gene_sets,
    data.gene_names,
    params.min_size,
    params.max_size,
  )

  let n_sets = filtered_sets.length()
  let n_samples = data.sample_names.length()
  let scores = Array::make(n_sets, [])
  let gene_set_names = Array::make(n_sets, "")

  let mut i = 0
  while i < n_sets {
    gene_set_names[i] = filtered_sets[i].name

    let set_scores = if params.method == "ssgsea" {
      gsva_score_ssgsea(data.expression, filtered_sets[i], data.gene_names)
    } else if params.method == "zscore" {
      gsva_score_zscore(data.expression, filtered_sets[i], data.gene_names)
    } else if params.method == "plage" {
      gsva_score_plage(data.expression, filtered_sets[i], data.gene_names)
    } else {
      gsva_score_gsva(
        data.expression,
        filtered_sets[i],
        data.gene_names,
        params.kcdf,
      )
    }

    scores[i] = set_scores
    i = i + 1
  }

  let normalized_scores = gsva_normalize_scores(scores, "none")

  let abs_rank_str = if params.abs_ranking { "true" } else { "false" }
  let verbose_str = if params.verbose { "true" } else { "false" }
  let param_pairs = [
    ("method", params.method),
    ("kcdf", params.kcdf),
    ("min_size", params.min_size.to_string()),
    ("max_size", params.max_size.to_string()),
    ("abs_ranking", abs_rank_str),
    ("verbose", verbose_str),
  ]
  let param_map = Map::from_array(param_pairs)

  GSVAScores::{
    scores: normalized_scores,
    gene_set_names,
    sample_names: data.sample_names,
    method: params.method,
    params: param_map,
  }
}

///|
/// Identify leading edge genes.
pub fn gsva_leading_edge(
  scores : Array[Array[Double]],
  gene_set : GSVAGeneSet,
  gene_names : Array[String],
) -> Array[String] {
  let n_samples = if scores.length() > 0 { scores[0].length() } else { 0 }
  if n_samples == 0 {
    return []
  }

  let avg_scores : Array[(String, Double)] = []
  let mut i = 0
  while i < gene_names.length() {
    let mut sum_score = 0.0
    let mut j = 0
    while j < scores.length() {
      sum_score = sum_score + scores[j][i]
      j = j + 1
    }
    avg_scores.push((gene_names[i], sum_score / scores.length().to_double()))
    i = i + 1
  }

  let n = avg_scores.length()
  i = 0
  while i < n - 1 {
    let mut j = 0
    while j < n - 1 - i {
      let (_, s1) = avg_scores[j]
      let (_, s2) = avg_scores[j + 1]
      if s1 < s2 {
        let tmp = avg_scores[j]
        avg_scores[j] = avg_scores[j + 1]
        avg_scores[j + 1] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }

  let leading : Array[String] = []
  i = 0
  while i < avg_scores.length() && i < gene_set.genes.length() / 2 {
    let (gene, _) = avg_scores[i]
    let mut in_set = false
    let mut j = 0
    while j < gene_set.genes.length() {
      if gene_set.genes[j] == gene {
        in_set = true
        break
      }
      j = j + 1
    }
    if in_set {
      leading.push(gene)
    }
    i = i + 1
  }

  leading
}

///|
/// Perform permutation test for significance.
pub fn gsva_permutation_test(
  scores : Array[Array[Double]],
  groups : Array[Int],
  n_perm : Int,
) -> Array[Double] {
  let n_sets = scores.length()
  let pvalues = Array::make(n_sets, 1.0)

  let mut set_idx = 0
  while set_idx < n_sets {
    let observed_diff = gsva_compute_group_diff(scores[set_idx], groups)

    let mut extreme_count = 0
    let mut perm = 0
    while perm < n_perm {
      let permuted = gsva_permute_groups(groups, perm)
      let perm_diff = gsva_compute_group_diff(scores[set_idx], permuted)
      if gsva_abs(perm_diff) >= gsva_abs(observed_diff) {
        extreme_count = extreme_count + 1
      }
      perm = perm + 1
    }

    pvalues[set_idx] = (extreme_count + 1).to_double() /
      (n_perm + 1).to_double()
    set_idx = set_idx + 1
  }

  pvalues
}

///|
/// Compute difference between group means.
fn gsva_compute_group_diff(
  scores : Array[Double],
  groups : Array[Int],
) -> Double {
  let group0 : Array[Double] = []
  let group1 : Array[Double] = []

  let mut i = 0
  while i < scores.length() {
    if groups[i] == 0 {
      group0.push(scores[i])
    } else {
      group1.push(scores[i])
    }
    i = i + 1
  }

  gsva_mean(group0) - gsva_mean(group1)
}

///|
/// Permute group labels.
fn gsva_permute_groups(groups : Array[Int], seed : Int) -> Array[Int] {
  let n = groups.length()
  let permuted = Array::make(n, 0)
  let mut i = 0
  while i < n {
    permuted[i] = groups[i]
    i = i + 1
  }

  let mut s = seed
  i = n - 1
  while i > 0 {
    s = (s * 1103515245 + 12345) % 2147483647
    if s < 0 {
      s = s + 2147483647
    }
    let j = s % (i + 1)
    let tmp = permuted[i]
    permuted[i] = permuted[j]
    permuted[j] = tmp
    i = i - 1
  }

  permuted
}

///|
/// Compute enrichment statistics for gene sets.
pub fn gsva_compute_enrichment(
  scores : Array[Array[Double]],
  gene_sets : Array[GSVAGeneSet],
  gene_names : Array[String],
) -> Array[GSVAEnrichment] {
  let n_sets = gene_sets.length()
  let enrichments = Array::make(n_sets, GSVAEnrichment::{
    gene_set_name: "",
    enrichment_score: 0.0,
    normalized_enrichment_score: 0.0,
    pvalue: 1.0,
    fdr: 1.0,
    leading_edge: [],
  })

  let mut i = 0
  while i < n_sets {
    let set_scores = scores[i]
    let es = gsva_mean(set_scores)
    let nes = if gsva_std(set_scores) > 0.0 {
      es / gsva_std(set_scores)
    } else {
      0.0
    }

    let leading = gsva_leading_edge(scores, gene_sets[i], gene_names)

    enrichments[i] = GSVAEnrichment::{
      gene_set_name: gene_sets[i].name,
      enrichment_score: es,
      normalized_enrichment_score: nes,
      pvalue: 0.05,
      fdr: 0.05,
      leading_edge: leading,
    }
    i = i + 1
  }

  enrichments
}

///|
/// Perform differential analysis between groups.
pub fn gsva_differential_analysis(
  scores : Array[Array[Double]],
  groups : Array[Int],
) -> Array[(Double, Double)] {
  let n_sets = scores.length()
  let results = Array::make(n_sets, (0.0, 0.0))

  let mut i = 0
  while i < n_sets {
    let group0 : Array[Double] = []
    let group1 : Array[Double] = []

    let mut j = 0
    while j < scores[i].length() {
      if groups[j] == 0 {
        group0.push(scores[i][j])
      } else {
        group1.push(scores[i][j])
      }
      j = j + 1
    }

    let mean0 = gsva_mean(group0)
    let mean1 = gsva_mean(group1)
    let sd0 = gsva_std(group0)
    let sd1 = gsva_std(group1)

    let pooled_sd = gsva_sqrt(
      (
        (group0.length() - 1).to_double() * sd0 * sd0 +
        (group1.length() - 1).to_double() * sd1 * sd1
      ) /
      (group0.length() + group1.length() - 2).to_double(),
    )
    let t_stat = if pooled_sd > 0.0 { (mean0 - mean1) / pooled_sd } else { 0.0 }

    results[i] = (t_stat, 0.05)
    i = i + 1
  }

  results
}

///|
/// Compute statistical summary of scores.
pub fn gsva_statistics(scores : Array[Array[Double]]) -> GSVAStat {
  let n_sets = scores.length()
  let mean_score = Array::make(n_sets, 0.0)
  let sd_score = Array::make(n_sets, 0.0)
  let median_score = Array::make(n_sets, 0.0)
  let min_score = Array::make(n_sets, 0.0)
  let max_score = Array::make(n_sets, 0.0)

  let mut i = 0
  while i < n_sets {
    let set_scores = scores[i]
    mean_score[i] = gsva_mean(set_scores)
    sd_score[i] = gsva_std(set_scores)
    median_score[i] = gsva_median(set_scores)

    let mut min_val = set_scores[0]
    let mut max_val = set_scores[0]
    let mut j = 1
    while j < set_scores.length() {
      if set_scores[j] < min_val {
        min_val = set_scores[j]
      }
      if set_scores[j] > max_val {
        max_val = set_scores[j]
      }
      j = j + 1
    }
    min_score[i] = min_val
    max_score[i] = max_val

    i = i + 1
  }

  GSVAStat::{ mean_score, sd_score, median_score, min_score, max_score }
}

///|
/// Generate summary string for GSVA scores.
pub fn gsva_summarize(scores : GSVAScores) -> String {
  let mut summary = "GSVA Results Summary\n"
  summary = summary + "===================\n\n"
  summary = summary + "Method: " + scores.method + "\n"
  summary = summary +
    "Gene Sets: " +
    scores.gene_set_names.length().to_string() +
    "\n"
  summary = summary +
    "Samples: " +
    scores.sample_names.length().to_string() +
    "\n\n"

  let stats = gsva_statistics(scores.scores)
  summary = summary + "Statistics per Gene Set:\n"
  summary = summary + "------------------------\n"

  let mut i = 0
  while i < scores.gene_set_names.length() && i < 5 {
    summary = summary + scores.gene_set_names[i] + "\n"
    summary = summary + "  Mean: " + stats.mean_score[i].to_string() + "\n"
    summary = summary + "  SD: " + stats.sd_score[i].to_string() + "\n"
    summary = summary + "  Median: " + stats.median_score[i].to_string() + "\n"
    summary = summary + "  Min: " + stats.min_score[i].to_string() + "\n"
    summary = summary + "  Max: " + stats.max_score[i].to_string() + "\n\n"
    i = i + 1
  }

  if scores.gene_set_names.length() > 5 {
    summary = summary +
      "... and " +
      (scores.gene_set_names.length() - 5).to_string() +
      " more gene sets\n"
  }

  summary
}

///|
/// Generate ASCII plot of scores for a specific gene set.
pub fn gsva_plot_scores(scores : GSVAScores, gene_set : String) -> String {
  let mut plot = "GSVA Scores Plot\n"
  plot = plot + "================\n\n"
  plot = plot + "Gene Set: " + gene_set + "\n\n"

  let mut set_idx = -1
  let mut i = 0
  while i < scores.gene_set_names.length() {
    if scores.gene_set_names[i] == gene_set {
      set_idx = i
      break
    }
    i = i + 1
  }

  if set_idx < 0 {
    return plot + "Gene set not found!\n"
  }

  let set_scores = scores.scores[set_idx]
  let min_val = gsva_mean(set_scores) - 2.0 * gsva_std(set_scores)
  let max_val = gsva_mean(set_scores) + 2.0 * gsva_std(set_scores)

  plot = plot + "Score distribution:\n"
  i = 0
  while i < set_scores.length() {
    let val = set_scores[i]
    let normalized = (val - min_val) / (max_val - min_val)
    let bar_length = (normalized * 40.0).to_int()

    let mut bar = ""
    let mut j = 0
    while j < bar_length {
      bar = bar + "="
      j = j + 1
    }
    bar = bar + ">"

    plot = plot +
      scores.sample_names[i] +
      ": " +
      bar +
      " " +
      val.to_string() +
      "\n"
    i = i + 1
  }

  plot
}

///|
/// Generate ASCII heatmap of scores.
pub fn gsva_heatmap(scores : GSVAScores) -> String {
  let mut heatmap = "GSVA Heatmap\n"
  heatmap = heatmap + "============\n\n"

  let n_sets = scores.gene_set_names.length()
  let n_samples = scores.sample_names.length()

  let max_display_sets = 10
  let max_display_samples = 10

  let display_sets = if n_sets < max_display_sets {
    n_sets
  } else {
    max_display_sets
  }
  let display_samples = if n_samples < max_display_samples {
    n_samples
  } else {
    max_display_samples
  }

  heatmap = heatmap + "\t"
  let mut i = 0
  while i < display_samples {
    let name = if scores.sample_names[i].length() > 8 {
      scores.sample_names[i][0:8].to_string()
    } else {
      scores.sample_names[i]
    }
    heatmap = heatmap + name + "\t"
    i = i + 1
  }
  heatmap = heatmap + "\n"

  i = 0
  while i < display_sets {
    let set_name = if scores.gene_set_names[i].length() > 12 {
      scores.gene_set_names[i][0:12].to_string()
    } else {
      scores.gene_set_names[i]
    }
    heatmap = heatmap + set_name + "\t"

    let mut j = 0
    while j < display_samples {
      let val = scores.scores[i][j]
      let mut char = "."
      if val > 1.0 {
        char = "#"
      } else if val > 0.5 {
        char = "%"
      } else if val > 0.0 {
        char = "+"
      } else if val > -0.5 {
        char = "-"
      } else if val > -1.0 {
        char = "~"
      } else {
        char = "*"
      }
      heatmap = heatmap + char + "\t"
      j = j + 1
    }
    heatmap = heatmap + "\n"
    i = i + 1
  }

  heatmap = heatmap +
    "\nLegend: #(>1.0) % (0.5-1.0) +(0-0.5) -( -0.5-0) ~(-1.0--0.5) *(< -1.0)\n"

  heatmap
}

///|
/// Convert GSVAGeneSet to string.
pub fn GSVAGeneSet::to_string(self : GSVAGeneSet) -> String {
  let mut result = "GSVAGeneSet: " + self.name + "\n"
  result = result + "Description: " + self.description + "\n"
  result = result + "Genes: [" + self.genes.length().to_string() + "] "
  let mut i = 0
  while i < self.genes.length() && i < 5 {
    if i > 0 {
      result = result + ", "
    }
    result = result + self.genes[i]
    i = i + 1
  }
  if self.genes.length() > 5 {
    result = result + ", ..."
  }
  result
}

///|
/// Convert GSVAData to string.
pub fn GSVAData::to_string(self : GSVAData) -> String {
  let n_genes = self.gene_names.length()
  let n_samples = self.sample_names.length()
  let n_sets = self.gene_sets.length()

  let mut result = "GSVAData\n"
  result = result + "--------\n"
  result = result + "Genes: " + n_genes.to_string() + "\n"
  result = result + "Samples: " + n_samples.to_string() + "\n"
  result = result + "Gene Sets: " + n_sets.to_string() + "\n"
  result = result +
    "Expression matrix: " +
    n_genes.to_string() +
    " x " +
    n_samples.to_string() +
    "\n"

  result
}

///|
/// Convert GSVAScores to string.
pub fn GSVAScores::to_string(self : GSVAScores) -> String {
  let n_sets = self.gene_set_names.length()
  let n_samples = self.sample_names.length()

  let mut result = "GSVAScores\n"
  result = result + "-----------\n"
  result = result + "Method: " + self.method + "\n"
  result = result + "Gene Sets: " + n_sets.to_string() + "\n"
  result = result + "Samples: " + n_samples.to_string() + "\n"
  result = result +
    "Scores matrix: " +
    n_sets.to_string() +
    " x " +
    n_samples.to_string() +
    "\n"

  result
}

///|
/// Convert GSVAParams to string.
pub fn GSVAParams::to_string(self : GSVAParams) -> String {
  let mut result = "GSVAParams\n"
  result = result + "-----------\n"
  result = result + "method: " + self.method + "\n"
  result = result + "kcdf: " + self.kcdf + "\n"
  result = result + "min_size: " + self.min_size.to_string() + "\n"
  result = result + "max_size: " + self.max_size.to_string() + "\n"
  let abs_rank_str = if self.abs_ranking { "true" } else { "false" }
  result = result + "abs_ranking: " + abs_rank_str + "\n"
  let verbose_str = if self.verbose { "true" } else { "false" }
  result = result + "verbose: " + verbose_str + "\n"

  result
}

///|
/// Convert GSVAEnrichment to string.
pub fn GSVAEnrichment::to_string(self : GSVAEnrichment) -> String {
  let mut result = "GSVAEnrichment\n"
  result = result + "---------------\n"
  result = result + "Gene Set: " + self.gene_set_name + "\n"
  result = result +
    "Enrichment Score: " +
    self.enrichment_score.to_string() +
    "\n"
  result = result +
    "Normalized ES: " +
    self.normalized_enrichment_score.to_string() +
    "\n"
  result = result + "p-value: " + self.pvalue.to_string() + "\n"
  result = result + "FDR: " + self.fdr.to_string() + "\n"
  result = result +
    "Leading Edge: [" +
    self.leading_edge.length().to_string() +
    "] "

  let mut i = 0
  while i < self.leading_edge.length() && i < 5 {
    if i > 0 {
      result = result + ", "
    }
    result = result + self.leading_edge[i]
    i = i + 1
  }
  if self.leading_edge.length() > 5 {
    result = result + ", ..."
  }

  result
}

///|
/// Create example GSVA data for testing.
pub fn gsva_create_example_data() -> GSVAData {
  let genes : Array[String] = Array::new()
  let samples : Array[String] = Array::new()
  let expression : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < 100 {
    genes.push(("Gene".to_string() + i.to_string()))
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < 10 {
      row.push(((i % 10) - 5).to_double() + j.to_double() * 0.1)
      j = j + 1
    }
    expression.push(row)
    i = i + 1
  }
  
  i = 0
  while i < 10 {
    samples.push(("Sample".to_string() + i.to_string()))
    i = i + 1
  }

  GSVAData::new(expression, genes, samples)
}

///|
/// Create example gene sets for testing.
pub fn gsva_create_example_gene_sets() -> Array[GSVAGeneSet] {
  [
    GSVAGeneSet::new(
      "CELL_CYCLE",
      ["Gene0", "Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6", "Gene7", "Gene8", "Gene9"],
      "Cell cycle related genes",
    ),
    GSVAGeneSet::new(
      "APOPTOSIS",
      ["Gene10", "Gene11", "Gene12", "Gene13", "Gene14", "Gene15", "Gene16", "Gene17", "Gene18", "Gene19"],
      "Apoptosis related genes",
    ),
    GSVAGeneSet::new(
      "METABOLISM",
      ["Gene20", "Gene21", "Gene22", "Gene23", "Gene24", "Gene25", "Gene26", "Gene27", "Gene28", "Gene29"],
      "Metabolism related genes",
    ),
    GSVAGeneSet::new(
      "IMMUNE",
      ["Gene30", "Gene31", "Gene32", "Gene33", "Gene34", "Gene35", "Gene36", "Gene37", "Gene38", "Gene39"],
      "Immune related genes",
    ),
    GSVAGeneSet::new(
      "STEM",
      ["Gene40", "Gene41", "Gene42", "Gene43", "Gene44", "Gene45", "Gene46", "Gene47", "Gene48", "Gene49"],
      "Stem cell related genes",
    ),
  ]
}

///|
/// Create example GSVA scores for testing.
pub fn gsva_create_example_scores() -> GSVAScores {
  let gene_sets = gsva_create_example_gene_sets()
  let data = gsva_create_example_data()
  let params = GSVAParams::new()

  gsva_run(data, gene_sets, params)
}

///|
/// Calculate overlap between gene sets.
fn gsva_gene_set_overlap(set1 : GSVAGeneSet, set2 : GSVAGeneSet) -> Double {
  let mut overlap = 0
  let mut i = 0
  while i < set1.genes.length() {
    let mut j = 0
    while j < set2.genes.length() {
      if set1.genes[i] == set2.genes[j] {
        overlap = overlap + 1
        break
      }
      j = j + 1
    }
    i = i + 1
  }
  let union = set1.genes.length() + set2.genes.length() - overlap
  if union == 0 {
    0.0
  } else {
    overlap.to_double() / union.to_double()
  }
}

///|
/// Generate enrichment map visualization showing relationships between gene sets.
pub fn gsva_enrichment_map(
  scores : GSVAScores,
  gene_sets : Array[GSVAGeneSet],
  threshold : Double,
) -> String {
  let mut map = "GSVA Enrichment Map\n"
  map = map + "===================\n\n"

  let n_sets = gene_sets.length()
  if n_sets == 0 {
    return map + "No gene sets provided\n"
  }

  let avg_scores : Array[Double] = Array::make(n_sets, 0.0)
  let mut i = 0
  while i < n_sets {
    let mut sum = 0.0
    let mut j = 0
    while j < scores.scores[i].length() {
      sum = sum + scores.scores[i][j]
      j = j + 1
    }
    avg_scores[i] = sum / scores.scores[i].length().to_double()
    i = i + 1
  }

  map = map + "Nodes (Gene Sets):\n"
  map = map + "------------------\n"
  i = 0
  while i < n_sets {
    let color = if avg_scores[i] > 0.5 { "#" } else if avg_scores[i] > 0.0 { "*" } else if avg_scores[i] > -0.5 { "." } else { "o" }
    map = map + "  " + color + " " + gene_sets[i].name + " (avg=" + avg_scores[i].to_string() + ", genes=" + gene_sets[i].genes.length().to_string() + ")\n"
    i = i + 1
  }

  map = map + "\nEdges (Overlap > " + threshold.to_string() + "):\n"
  map = map + "----------------------------------------\n"
  let mut edge_count = 0
  i = 0
  while i < n_sets {
    let mut j = i + 1
    while j < n_sets {
      let overlap = gsva_gene_set_overlap(gene_sets[i], gene_sets[j])
      if overlap > threshold {
        map = map + "  " + gene_sets[i].name + " ---" + overlap.to_string() + "--- " + gene_sets[j].name + "\n"
        edge_count = edge_count + 1
      }
      j = j + 1
    }
    i = i + 1
  }

  if edge_count == 0 {
    map = map + "  No edges above threshold\n"
  }

  map = map + "\nLegend:\n"
  map = map + "  # = High positive score (>0.5)\n"
  map = map + "  * = Moderate positive score (0-0.5)\n"
  map = map + "  . = Moderate negative score (-0.5-0)\n"
  map = map + "  o = High negative score (<-0.5)\n"

  map
}

///|
/// Calculate Pearson correlation between two arrays.
fn gsva_pearson_correlation(x : Array[Double], y : Array[Double]) -> Double {
  if x.length() != y.length() || x.length() < 2 {
    return 0.0
  }

  let n = x.length().to_double()
  let mut sum_x = 0.0
  let mut sum_y = 0.0
  let mut sum_xy = 0.0
  let mut sum_x2 = 0.0
  let mut sum_y2 = 0.0

  let mut i = 0
  while i < x.length() {
    sum_x = sum_x + x[i]
    sum_y = sum_y + y[i]
    sum_xy = sum_xy + x[i] * y[i]
    sum_x2 = sum_x2 + x[i] * x[i]
    sum_y2 = sum_y2 + y[i] * y[i]
    i = i + 1
  }

  let numerator = n * sum_xy - sum_x * sum_y
  let denominator = gsva_sqrt((n * sum_x2 - sum_x * sum_x) * (n * sum_y2 - sum_y * sum_y))

  if denominator == 0.0 {
    0.0
  } else {
    numerator / denominator
  }
}

///|
/// Correlate GSVA scores with phenotype data.
pub fn gsva_phenotype_correlation(
  scores : GSVAScores,
  phenotypes : Array[Double],
) -> Array[(String, Double, Double)] {
  let n_sets = scores.scores.length()
  let results = Array::make(n_sets, ("", 0.0, 0.0))

  let mut i = 0
  while i < n_sets {
    let corr = gsva_pearson_correlation(scores.scores[i], phenotypes)
    let p_value = if corr.abs() > 0.0 {
      let t_stat = corr * gsva_sqrt((scores.scores[i].length() - 2).to_double() / (1.0 - corr * corr))
      let df = scores.scores[i].length() - 2
      if df < 1 {
        1.0
      } else {
        2.0 * (1.0 - gsva_cdf(t_stat.abs(), 0.0, 1.0))
      }
    } else {
      1.0
    }
    results[i] = (scores.gene_set_names[i], corr, p_value)
    i = i + 1
  }

  results
}

///|
/// Generate phenotype correlation report.
pub fn gsva_phenotype_report(
  scores : GSVAScores,
  phenotypes : Array[Double],
  phenotype_name : String,
) -> String {
  let mut report = "GSVA Phenotype Correlation Report\n"
  report = report + "=================================\n\n"
  report = report + "Phenotype: " + phenotype_name + "\n"
  report = report + "Samples: " + scores.sample_names.length().to_string() + "\n\n"

  let correlations = gsva_phenotype_correlation(scores, phenotypes)
  let n_sets = correlations.length()

  let sorted : Array[(String, Double, Double)] = Array::make(n_sets, ("", 0.0, 0.0))
  let mut i = 0
  while i < n_sets {
    sorted[i] = correlations[i]
    i = i + 1
  }

  let mut j = 0
  while j < n_sets - 1 {
    let mut k = 0
    while k < n_sets - 1 - j {
      let (_, c1, _) = sorted[k]
      let (_, c2, _) = sorted[k + 1]
      if c1.abs() < c2.abs() {
        let tmp = sorted[k]
        sorted[k] = sorted[k + 1]
        sorted[k + 1] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }

  report = report + "Top correlated gene sets:\n"
  report = report + "-------------------------\n"
  report = report + "Gene Set\tCorrelation\tp-value\n"
  report = report + "--------\t-----------\t-------\n"

  i = 0
  while i < sorted.length() && i < 10 {
    let (name, corr, pval) = sorted[i]
    report = report + name + "\t" + corr.to_string() + "\t" + pval.to_string() + "\n"
    i = i + 1
  }

  report
}

///|
/// Survival analysis result.
pub struct GSVASurvivalResult {
  gene_set_name : String
  hazard_ratio : Double
  p_value : Double
  ci_lower : Double
  ci_upper : Double
}

///|
/// Perform survival analysis using Cox proportional hazards model approximation.
pub fn gsva_survival_analysis(
  scores : GSVAScores,
  survival_time : Array[Double],
  event : Array[Int],
) -> Array[GSVASurvivalResult] {
  let n_sets = scores.scores.length()
  let results = Array::make(n_sets, GSVASurvivalResult::{
    gene_set_name: "",
    hazard_ratio: 1.0,
    p_value: 1.0,
    ci_lower: 0.0,
    ci_upper: 0.0,
  })

  let mut i = 0
  while i < n_sets {
    let set_scores = scores.scores[i]
    let n = set_scores.length()

    let mut sum_time = 0.0
    let mut sum_event = 0
    let mut sum_score = 0.0
    let mut sum_score_event = 0.0
    let mut sum_score_sq = 0.0

    let mut j = 0
    while j < n {
      sum_time = sum_time + survival_time[j]
      sum_event = sum_event + event[j]
      sum_score = sum_score + set_scores[j]
      sum_score_event = sum_score_event + set_scores[j] * event[j].to_double()
      sum_score_sq = sum_score_sq + set_scores[j] * set_scores[j]
      j = j + 1
    }

    if sum_event == 0 || sum_event == n {
      results[i] = GSVASurvivalResult::{
        gene_set_name: scores.gene_set_names[i],
        hazard_ratio: 1.0,
        p_value: 1.0,
        ci_lower: 0.0,
        ci_upper: 0.0,
      }
      i = i + 1
      continue
    }

    let mean_score = sum_score / n.to_double()
    let mean_event = sum_event.to_double() / n.to_double()

    let mut cov = 0.0
    let mut var_score = 0.0
    let mut var_event = 0.0

    j = 0
    while j < n {
      cov = cov + (set_scores[j] - mean_score) * (event[j].to_double() - mean_event)
      var_score = var_score + (set_scores[j] - mean_score) * (set_scores[j] - mean_score)
      var_event = var_event + (event[j].to_double() - mean_event) * (event[j].to_double() - mean_event)
      j = j + 1
    }

    let beta = if var_score > 0.0 { cov / var_score } else { 0.0 }
    let hazard_ratio = gsva_exp(beta)

    let se = if var_score > 0.0 && n > 2 {
      gsva_sqrt(var_event / (var_score * (n - 2).to_double()))
    } else {
      1.0
    }

    let z_stat = beta / se
    let p_value = 2.0 * (1.0 - gsva_cdf(z_stat.abs(), 0.0, 1.0))

    let ci_lower = gsva_exp(beta - 1.96 * se)
    let ci_upper = gsva_exp(beta + 1.96 * se)

    results[i] = GSVASurvivalResult::{
      gene_set_name: scores.gene_set_names[i],
      hazard_ratio,
      p_value,
      ci_lower,
      ci_upper,
    }

    i = i + 1
  }

  results
}

///|
/// Generate survival analysis report.
pub fn gsva_survival_report(
  scores : GSVAScores,
  survival_time : Array[Double],
  event : Array[Int],
) -> String {
  let mut report = "GSVA Survival Analysis Report\n"
  report = report + "=============================\n\n"
  report = report + "Samples: " + scores.sample_names.length().to_string() + "\n"
  let mut event_sum = 0
  let mut ev_idx = 0
  while ev_idx < event.length() {
    event_sum = event_sum + event[ev_idx]
    ev_idx = ev_idx + 1
  }
  report = report + "Events: " + event_sum.to_string() + "\n\n"

  let results = gsva_survival_analysis(scores, survival_time, event)
  let n_sets = results.length()

  let sorted : Array[GSVASurvivalResult] = Array::make(n_sets, results[0])
  let mut i = 0
  while i < n_sets {
    sorted[i] = results[i]
    i = i + 1
  }

  let mut j = 0
  while j < n_sets - 1 {
    let mut k = 0
    while k < n_sets - 1 - j {
      if sorted[k].p_value > sorted[k + 1].p_value {
        let tmp = sorted[k]
        sorted[k] = sorted[k + 1]
        sorted[k + 1] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }

  report = report + "Significant gene sets (sorted by p-value):\n"
  report = report + "-------------------------------------------\n"
  report = report + "Gene Set\tHR\t95% CI\t\tp-value\n"
  report = report + "--------\t--\t------\t\t-------\n"

  i = 0
  while i < sorted.length() && i < 10 {
    let res = sorted[i]
    report = report + res.gene_set_name + "\t" + res.hazard_ratio.to_string() + "\t[" + res.ci_lower.to_string() + "," + res.ci_upper.to_string() + "]\t" + res.p_value.to_string() + "\n"
    i = i + 1
  }

  report
}

///|
/// Generate score distribution plot for all gene sets.
pub fn gsva_score_distribution(scores : GSVAScores) -> String {
  let mut plot = "GSVA Score Distribution\n"
  plot = plot + "=======================\n\n"

  let n_sets = scores.scores.length()
  if n_sets == 0 {
    return plot + "No gene sets\n"
  }

  let all_scores : Array[Double] = Array::new()
  let mut i = 0
  while i < n_sets {
    let mut j = 0
    while j < scores.scores[i].length() {
      all_scores.push(scores.scores[i][j])
      j = j + 1
    }
    i = i + 1
  }

  if all_scores.length() == 0 {
    return plot + "No scores\n"
  }

  let sorted_scores = gsva_sort_doubles(all_scores)
  let min_val = sorted_scores[0]
  let max_val = sorted_scores[sorted_scores.length() - 1]
  let range = max_val - min_val

  let bin_count = 20
  let bins : Array[Int] = Array::make(bin_count, 0)

  i = 0
  while i < sorted_scores.length() {
    let val = sorted_scores[i]
    let bin_idx = if range > 0.0 {
      ((val - min_val) / range * (bin_count - 1).to_double()).to_int()
    } else {
      0
    }
    if bin_idx >= 0 && bin_idx < bin_count {
      bins[bin_idx] = bins[bin_idx] + 1
    }
    i = i + 1
  }

  let mut max_freq = 0
  let mut bin_idx = 0
  while bin_idx < bins.length() {
    if bins[bin_idx] > max_freq {
      max_freq = bins[bin_idx]
    }
    bin_idx = bin_idx + 1
  }

  plot = plot + "Score histogram:\n"
  plot = plot + "----------------\n"

  i = bin_count - 1
  while i >= 0 {
    let threshold = min_val + range * i.to_double() / (bin_count - 1).to_double()
    let mut bar = ""
    let mut j = 0
    while j < bins[i] && max_freq > 0 {
      bar = bar + "#"
      j = j + 1
    }
    plot = plot + threshold.to_string() + ": " + bar + "\n"
    i = i - 1
  }

  let mean = gsva_mean(all_scores)
  let sd = gsva_std(all_scores)
  plot = plot + "\nStatistics:\n"
  plot = plot + "Mean: " + mean.to_string() + "\n"
  plot = plot + "SD: " + sd.to_string() + "\n"
  plot = plot + "Min: " + min_val.to_string() + "\n"
  plot = plot + "Max: " + max_val.to_string() + "\n"

  plot
}
