///|
/// infercnv module - Copy Number Variation inference from single-cell RNA-seq.
/// Inspired by Broad Institute inferCNV / Bioconductor infercnv package and
/// infercnvpy (Python reimplementation).
///
/// References:
///   - https://bioconductor.org/packages/release/bioc/html/infercnv.html
///   - Tirosh I et al. (2016) Dissecting the multicellular ecosystem of metastatic
///     melanoma by single-cell RNA-seq. Science 352(6282):189-196.
///   - infercnvpy: https://icbi-lab.github.io/infercnvpy/
///
/// Core algorithm workflow:
///   1. Order genes by chromosome and genomic position, exclude chrM/chrX/chrY.
///   2. Select reference expression profile:
///        - default: mean across all cells
///        - cell-type: mean of specified "normal" reference categories
///        - custom: user-provided reference array
///   3. Per-cell chunk processing:
///        a. Compute log fold change (single ref: subtraction; multiple refs: bounded)
///        b. Clip values to +/- lfc_cap (default 1.5)
///        c. Smooth across genomic position with pyramidal-weighted running window
///        d. Center each cell by its median
///        e. Apply noise threshold (set |value| < noise to 0)
///   4. Produce CNV matrix (cells x ordered genes) and per-cell CNV scores.

// ============================================================================
// Enumerations and helper constructors
// ============================================================================

///|
/// Method used to build the reference expression profile.
pub enum ReferenceMethod {
  /// Use the mean of all cells as reference (no prior knowledge required).
  GlobalMean
  /// Use the mean of specified reference cell categories (e.g. immune / stromal
  /// cells annotated in the cell metadata).
  ReferenceCategories
  /// User-provided custom reference vector (length = number of ordered genes).
  Custom(Array[Double])
}

///|
/// Global-mean reference method.
pub fn ref_method_global_mean() -> ReferenceMethod {
  ReferenceMethod::GlobalMean
}

///|
/// Reference-categories reference method.
pub fn ref_method_reference_categories() -> ReferenceMethod {
  ReferenceMethod::ReferenceCategories
}

///|
/// Custom reference method wrapping a user-provided per-gene vector.
pub fn ref_method_custom(ref_profile : Array[Double]) -> ReferenceMethod {
  ReferenceMethod::Custom(ref_profile)
}

// ============================================================================
// Core data structures
// ============================================================================

///|
/// Genomic position of a gene.
pub struct GenePosition {
  gene_id : String
  chromosome : String
  start : Int
  end : Int
} derive(Debug, Default)

///|
/// Create a GenePosition.
pub fn gene_position(
  gene_id : String,
  chromosome : String,
  start : Int,
  end : Int,
) -> GenePosition {
  {
    gene_id: gene_id,
    chromosome: chromosome,
    start: start,
    end: end,
  }
}

///|
/// An ordered gene set together with per-gene genomic positions.
pub struct OrderedGenes {
  gene_order : Array[GenePosition]        // sorted by chr, then start position
  n_genes : Int
} derive(Debug, Default)

///|
/// Build an OrderedGenes struct from an unordered list of GenePositions.
/// Chromosomes are sorted naturally by number/letter; within each chromosome
/// genes are sorted by start position.
pub fn ordered_genes(raw : Array[GenePosition]) -> OrderedGenes {
  let arr = raw.copy()
  // bubble sort by natural chromosome key then start position
  let n = arr.length()
  for i in 0..<n {
    for j in 0..<(n - 1 - i) {
      let ca = natural_chrom_key(arr[j].chromosome)
      let cb = natural_chrom_key(arr[j + 1].chromosome)
      let should_swap = if ca != cb { ca > cb } else { arr[j].start > arr[j + 1].start }
      if should_swap {
        let tmp = arr[j]
        arr[j] = arr[j + 1]
        arr[j + 1] = tmp
      }
    }
  }
  {
    gene_order: arr,
    n_genes: arr.length(),
  }
}

fn natural_chrom_key(chr : String) -> String {
  // Trim optional "chr" prefix, then:
  //   numeric chromosomes 1..22 -> "01" .. "22"
  //   X -> "A0", Y -> "A1", M/MT -> "A2" (excluded by default anyway)
  //   other -> "B" + original
  let s = chr.to_lower()
  let base = if s.starts_with("chr") { s.substring(start=3, end=s.length()) } else { s }
  match parse_double(base) {
    Some(v) => {
      let n = v.to_int()
      if n < 10 { "0" + n.to_string() } else { n.to_string() }
    }
    None => {
      if base == "x" { "A0" }
      else if base == "y" { "A1" }
      else if base == "m" || base == "mt" { "A2" }
      else { "B" + base }
    }
  }
}

///|
/// Exclude genes on chromosomes such as chrM, chrX, chrY.
pub fn exclude_chromosomes(
  og : OrderedGenes,
  excluded : Array[String],
) -> OrderedGenes {
  let lower_excluded = excluded.map(fn (s) -> String { s.to_lower() })
  let kept = og.gene_order.filter(fn (gp : GenePosition) -> Bool {
    let c = gp.chromosome.to_lower()
    not(lower_excluded.contains(c))
  })
  ordered_genes(kept)
}

// ============================================================================
// Input data: expression matrix + cell annotations
// ============================================================================

///|
/// Cell metadata record.
pub struct CellAnnotation {
  cell_id : String
  category : String   // e.g. cell type / sample label
} derive(Debug, Default)

///|
/// Create a cell annotation.
pub fn cell_annotation(cell_id : String, category : String) -> CellAnnotation {
  { cell_id: cell_id, category: category }
}

///|
/// The input expression matrix together with cell annotations and ordered genes.
/// Rows of `expression` are cells; columns are genes in the order of
/// `ordered_genes.gene_order`. Values are expected to be log-normalized
/// (typically log2(TP100K+1)). Use `log_normalize_counts` if you have raw counts.
pub struct InferCNVInput {
  n_cells : Int
  n_genes : Int
  expression : Array[Array[Double]]     // cell x gene (ordered)
  annotations : Array[CellAnnotation]   // length n_cells
  ordered_genes : OrderedGenes
} derive(Debug, Default)

// ============================================================================
// Utility: raw count normalization (TP100K -> log)
// ============================================================================

///|
/// Log-normalize a raw counts matrix (cell x gene) to log2(TP100K + 1).
/// Cells with zero total counts are preserved as zeros.
pub fn log_normalize_counts(
  raw_counts : Array[Array[Double]],
  target_sum? : Double = 100000.0,
) -> Array[Array[Double]] {
  let n_cells = raw_counts.length()
  let out : Array[Array[Double]] = []
  for i in 0..<n_cells {
    let row = raw_counts[i]
    let total = infercnv_array_sum(row)
    let n_genes = row.length()
    let norm_row : Array[Double] = []
    if total <= 0.0 {
      for _j in 0..<n_genes { norm_row.push(0.0) }
    } else {
      let scale = target_sum / total
      for j in 0..<n_genes {
        norm_row.push(@math.log2(row[j] * scale + 1.0))
      }
    }
    out.push(norm_row)
  }
  out
}

fn infercnv_array_sum(arr : Array[Double]) -> Double {
  let mut sum = 0.0
  for v in arr { sum += v }
  sum
}

// ============================================================================
// Result structure
// ============================================================================

///|
/// CNV inference result.
pub struct CNVResult {
  n_cells : Int
  n_genes : Int
  cnv_matrix : Array[Array[Double]]     // cell x gene (ordered)
  cell_ids : Array[String]
  gene_ids : Array[String]
  chromosome_boundaries : Array[(String, Int, Int)]  // (chr, start_idx, end_idx) per chromosome
  per_cell_cnv_score : Array[Double]    // mean |CNV| of each cell across all genes
  per_cluster_cnv_score : Array[(String, Double)] // (category, mean score)
}

///|
/// Parameters exposed to the user for controlling the infercnv pipeline.
pub struct CNVParams {
  window_size : Int           // number of adjacent genes in smoothing window
  lfc_cap : Double            // clip fold changes to +/- lfc_cap
  noise_threshold : Double    // values with |x| < noise are set to 0
  reference_method : ReferenceMethod
  reference_categories : Array[String]  // used with ReferenceCategories method
}

///|
/// Default parameters used in infercnvpy.
pub fn default_cnv_params() -> CNVParams {
  {
    window_size: 100,
    lfc_cap: 1.5,
    noise_threshold: 0.2,
    reference_method: ReferenceMethod::GlobalMean,
    reference_categories: [],
  }
}

///|
/// Create custom CNVParams from outside the package.
pub fn make_cnv_params(
  window_size : Int,
  lfc_cap : Double,
  noise_threshold : Double,
  reference_method : ReferenceMethod,
  reference_categories : Array[String],
) -> CNVParams {
  {
    window_size: window_size,
    lfc_cap: lfc_cap,
    noise_threshold: noise_threshold,
    reference_method: reference_method,
    reference_categories: reference_categories,
  }
}

// ============================================================================
// Reference profile computation
// ============================================================================

fn compute_reference_profile(
  input : InferCNVInput,
  params : CNVParams,
) -> Array[Array[Double]] {  // returns array of reference profiles; single ref = 1 profile
  match params.reference_method {
    ReferenceMethod::Custom(ref_profile) => {
      // ensure length matches n_genes
      let n = input.n_genes
      let arr = ref_profile.copy()
      while arr.length() < n { arr.push(0.0) }
      // truncate to n genes
      let truncated : Array[Double] = []
      for i in 0..<n { truncated.push(arr[i]) }
      [truncated]
    }
    ReferenceMethod::GlobalMean => {
      [mean_per_gene(input.expression)]
    }
    ReferenceMethod::ReferenceCategories => {
      let profiles : Array[Array[Double]] = []
      for cat in params.reference_categories {
        let indices = cells_in_category(input.annotations, cat)
        if indices.length() == 0 { continue }
        let subset : Array[Array[Double]] = []
        for idx in indices { subset.push(input.expression[idx]) }
        profiles.push(mean_per_gene(subset))
      }
      if profiles.length() == 0 {
        // fallback to global mean if no reference category matches
        [mean_per_gene(input.expression)]
      } else {
        profiles
      }
    }
  }
}

fn mean_per_gene(matrix : Array[Array[Double]]) -> Array[Double] {
  let n_cells = matrix.length()
  if n_cells == 0 { return [] }
  let n_genes = matrix[0].length()
  let out : Array[Double] = []
  for j in 0..<n_genes {
    let mut sum = 0.0
    for i in 0..<n_cells { sum += matrix[i][j] }
    out.push(sum / (n_cells.to_double()))
  }
  out
}

fn cells_in_category(
  ann : Array[CellAnnotation],
  cat : String,
) -> Array[Int] {
  let out : Array[Int] = []
  for i in 0..<ann.length() {
    if ann[i].category == cat { out.push(i) }
  }
  out
}

// ============================================================================
// Core step: compute log fold change for a single cell
// ============================================================================

fn compute_lfc(
  cell_expr : Array[Double],
  references : Array[Array[Double]],
  lfc_cap : Double,
) -> Array[Double] {
  let n = cell_expr.length()
  let out : Array[Double] = []
  if references.length() == 1 {
    let ref1 = references[0]
    for j in 0..<n {
      let v = cell_expr[j] - ref1[j]
      let clipped = if v > lfc_cap { lfc_cap } else if v < -lfc_cap { -lfc_cap } else { v }
      out.push(clipped)
    }
  } else {
    // Multiple reference categories -> bounded log fold change.
    // Values in [min_ref, max_ref] -> 0
    // Values above max -> cell - max
    // Values below min -> cell - min
    let n_refs = references.length()
    for j in 0..<n {
      let mut min_r = references[0][j]
      let mut max_r = references[0][j]
      for r in 1..<n_refs {
        let v = references[r][j]
        if v < min_r { min_r = v }
        if v > max_r { max_r = v }
      }
      let c = cell_expr[j]
      let lfc = if c > max_r { c - max_r }
                else if c < min_r { c - min_r }
                else { 0.0 }
      let clipped = if lfc > lfc_cap { lfc_cap }
                    else if lfc < -lfc_cap { -lfc_cap }
                    else { lfc }
      out.push(clipped)
    }
  }
  out
}

// ============================================================================
// Genomic smoothing per chromosome (pyramidal-weighted running mean)
// ============================================================================

fn chromosome_boundaries(
  og : OrderedGenes,
) -> Array[(String, Int, Int)] {
  let result : Array[(String, Int, Int)] = []
  if og.n_genes == 0 { return result }
  let mut current_chr = og.gene_order[0].chromosome
  let mut start_idx = 0
  for i in 1..<og.n_genes {
    if og.gene_order[i].chromosome != current_chr {
      result.push((current_chr, start_idx, i - 1))
      current_chr = og.gene_order[i].chromosome
      start_idx = i
    }
  }
  result.push((current_chr, start_idx, og.n_genes - 1))
  result
}

/// Pyramidal weights of length 2*half + 1:
///   w[i] = (half + 1 - |i - half|) / (half + 1)
/// Example half=2: [0.33, 0.66, 1.0, 0.66, 0.33]
fn pyramidal_weights(half : Int) -> Array[Double] {
  let total = 2 * half + 1
  let w : Array[Double] = []
  let denom = (half + 1).to_double()
  for i in 0..<total {
    let d = (i - half).abs()
    w.push(((half + 1 - d).to_double()) / denom)
  }
  w
}

fn array_sum_weighted(vec : Array[Double], w : Array[Double]) -> Double {
  let mut sum = 0.0
  for i in 0..<w.length() { sum += vec[i] * w[i] }
  sum
}

fn smooth_vector(vec : Array[Double], half_window : Int) -> Array[Double] {
  let n = vec.length()
  if half_window <= 0 || n <= 1 { return vec.copy() }
  let w = pyramidal_weights(half_window)
  let w_total = infercnv_array_sum(w)
  let out : Array[Double] = []
  for i in 0..<n {
    let mut sum = 0.0
    let mut sum_w = 0.0
    for k in 0..<w.length() {
      let pos = i - half_window + k
      if pos >= 0 && pos < n {
        sum += vec[pos] * w[k]
        sum_w += w[k]
      }
    }
    if sum_w > 0.0 {
      out.push(sum / sum_w)
    } else {
      out.push(0.0)
    }
    let _ = w_total
  }
  out
}

fn smooth_per_chromosome(
  vec : Array[Double],
  boundaries : Array[(String, Int, Int)],
  window_size : Int,
) -> Array[Double] {
  let half = window_size / 2
  let result = vec.copy()
  for entry in boundaries {
    let (chr, s, e) = entry
    let len = e - s + 1
    let segment : Array[Double] = []
    for i in s..=e { segment.push(result[i]) }
    let smoothed = smooth_vector(segment, half)
    for i in 0..<len { result[s + i] = smoothed[i] }
    let _ = chr  // silence unused warning
  }
  result
}

// ============================================================================
// Median of a vector (for cell-level centering)
// ============================================================================

fn median_of(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 { return 0.0 }
  let sorted = arr.copy()
  // simple insertion sort for Double
  for i in 1..<n {
    let key = sorted[i]
    let mut j = i - 1
    while j >= 0 && sorted[j] > key {
      sorted[j + 1] = sorted[j]
      j = j - 1
    }
    sorted[j + 1] = key
  }
  if n % 2 == 1 { sorted[n / 2] }
  else {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  }
}

fn center_by_median(vec : Array[Double]) -> Array[Double] {
  let m = median_of(vec)
  vec.map(fn (x) -> Double { x - m })
}

fn apply_noise_threshold(vec : Array[Double], noise : Double) -> Array[Double] {
  vec.map(fn (x) -> Double { if x.abs() < noise { 0.0 } else { x } })
}

// ============================================================================
// Mean absolute per-cell CNV score
// ============================================================================

fn per_cell_cnv_scores(matrix : Array[Array[Double]]) -> Array[Double] {
  let n_cells = matrix.length()
  let out : Array[Double] = []
  for i in 0..<n_cells {
    let row = matrix[i]
    if row.length() == 0 {
      out.push(0.0)
    } else {
      let mut sum_abs = 0.0
      for x in row { sum_abs += x.abs() }
      out.push(sum_abs / (row.length().to_double()))
    }
  }
  out
}

fn per_category_scores(
  scores : Array[Double],
  ann : Array[CellAnnotation],
) -> Array[(String, Double)] {
  // group indices by category using parallel arrays
  let cats : Array[String] = []
  let sums : Array[Double] = []
  let counts : Array[Int] = []
  for i in 0..<ann.length() {
    let cat = ann[i].category
    let mut found = false
    for k in 0..<cats.length() {
      if cats[k] == cat {
        sums[k] += scores[i]
        counts[k] += 1
        found = true
        break
      }
    }
    if not(found) {
      cats.push(cat)
      sums.push(scores[i])
      counts.push(1)
    }
  }
  let out : Array[(String, Double)] = []
  for k in 0..<cats.length() {
    if counts[k] == 0 {
      out.push((cats[k], 0.0))
    } else {
      out.push((cats[k], sums[k] / counts[k].to_double()))
    }
  }
  out
}

// ============================================================================
// Top-level entrypoint
// ============================================================================

///|
/// Run the inferCNV pipeline end-to-end.
///
/// Parameters:
///   - input: expression matrix + annotations + ordered genes
///   - params: window_size, lfc_cap, noise_threshold, reference method
///   - excluded_chromosomes: chr names to drop (default [chrM, chrX, chrY])
pub fn run_infercnv(
  input : InferCNVInput,
  params? : CNVParams = default_cnv_params(),
  excluded_chromosomes? : Array[String] = ["chrM", "chrX", "chrY"],
) -> CNVResult {
  // Step 1: filter out excluded chromosomes and reorder
  let filtered = exclude_chromosomes(input.ordered_genes, excluded_chromosomes)

  // Step 1b: keep only expression columns whose genes survived the filter.
  // Reorder input expression to match the filtered gene order by gene_id lookup.
  let keep_gene_ids = filtered.gene_order.map(fn (g) -> String { g.gene_id })
  let original_gene_ids = input.ordered_genes.gene_order.map(fn (g) -> String { g.gene_id })
  // build index map
  let orig_keys : Array[String] = []
  let orig_vals : Array[Int] = []
  for j in 0..<original_gene_ids.length() {
    orig_keys.push(original_gene_ids[j])
    orig_vals.push(j)
  }
  let remapped_expr : Array[Array[Double]] = []
  for i in 0..<input.n_cells {
    let new_row : Array[Double] = []
    for gid in keep_gene_ids {
      let mut found_idx = -1
      for k in 0..<orig_keys.length() {
        if orig_keys[k] == gid { found_idx = orig_vals[k]; break }
      }
      if found_idx >= 0 {
        new_row.push(input.expression[i][found_idx])
      } else {
        new_row.push(0.0)
      }
    }
    remapped_expr.push(new_row)
  }

  // Build remapped input
  let effective_input : InferCNVInput = {
    n_cells: input.n_cells,
    n_genes: filtered.n_genes,
    expression: remapped_expr,
    annotations: input.annotations,
    ordered_genes: filtered,
  }

  // Step 2: reference profile
  let references = compute_reference_profile(effective_input, params)
  let boundaries = chromosome_boundaries(effective_input.ordered_genes)

  // Step 3: per-cell processing
  let cnv_matrix : Array[Array[Double]] = []
  for i in 0..<effective_input.n_cells {
    let cell_expr = effective_input.expression[i]
    let lfc = compute_lfc(cell_expr, references, params.lfc_cap)
    let smoothed = smooth_per_chromosome(lfc, boundaries, params.window_size)
    let centered = center_by_median(smoothed)
    let denoised = apply_noise_threshold(centered, params.noise_threshold)
    cnv_matrix.push(denoised)
  }

  let cell_ids = effective_input.annotations.map(fn (a) -> String { a.cell_id })
  let gene_ids = effective_input.ordered_genes.gene_order.map(fn (g) -> String { g.gene_id })

  let cell_scores = per_cell_cnv_scores(cnv_matrix)
  let cluster_scores = per_category_scores(cell_scores, effective_input.annotations)

  {
    n_cells: effective_input.n_cells,
    n_genes: effective_input.n_genes,
    cnv_matrix: cnv_matrix,
    cell_ids: cell_ids,
    gene_ids: gene_ids,
    chromosome_boundaries: boundaries,
    per_cell_cnv_score: cell_scores,
    per_cluster_cnv_score: cluster_scores,
  }
}

// ============================================================================
// Example / synthetic data generator
// ============================================================================

///|
/// Generate a synthetic infercnv example mimicking a tumour with normal cells.
///
/// Creates `n_tumour + n_normal` cells over `n_genes_per_chr * n_chr` genes
/// spread across `n_chr` autosomes. Tumour cells have a simulated deletion on
/// the first half of chromosome 1 (low expression) and amplification on the
/// second half of chromosome 2 (high expression). Returns a tuple
/// `(input, tumour_categories, normal_categories)`.
pub fn infercnv_sample_data(
  n_tumour? : Int = 40,
  n_normal? : Int = 20,
  n_chr? : Int = 4,
  n_genes_per_chr? : Int = 40,
  seed? : Int = 42,
) -> (InferCNVInput, Array[String], Array[String]) {
  let mut rng_state = seed
  let next_rand = fn () -> Double {
    // simple LCG
    rng_state = (rng_state * 1103515245 + 12345) & 0x7fffffff
    ((rng_state & 0xffffff).to_double()) / 16777215.0
  }
  let normal_mean = 3.0
  let n_genes = n_chr * n_genes_per_chr
  // Build ordered genes
  let positions : Array[GenePosition] = []
  for c in 0..<n_chr {
    let chr = "chr" + (c + 1).to_string()
    for g in 0..<n_genes_per_chr {
      let start = g * 10000
      let gid = "Gene_" + (c + 1).to_string() + "_" + (g + 1).to_string()
      positions.push(gene_position(gid, chr, start, start + 2000))
    }
  }
  let og = ordered_genes(positions)

  // helper: generate one cell with optional CNV effect
  let make_cell = fn (is_tumour : Bool) -> Array[Double] {
    let expr : Array[Double] = []
    for c in 0..<n_chr {
      for g in 0..<n_genes_per_chr {
        // base expression from normal distribution-ish
        let base = normal_mean + (next_rand() - 0.5) * 1.5
        let mut val = base
        if is_tumour {
          if c == 0 && g < n_genes_per_chr / 2 {
            // deletion on chr1 first half
            val -= 1.2
          }
          if c == 1 && g >= n_genes_per_chr / 2 {
            // amplification on chr2 second half
            val += 1.2
          }
        }
        if val < 0.0 { val = 0.0 }
        expr.push(val)
      }
    }
    expr
  }

  // Build expression and annotations
  let expressions : Array[Array[Double]] = []
  let annotations : Array[CellAnnotation] = []
  for i in 0..<n_tumour {
    expressions.push(make_cell(true))
    annotations.push(cell_annotation("Tumour_" + (i + 1).to_string(), "Tumour"))
  }
  for i in 0..<n_normal {
    expressions.push(make_cell(false))
    annotations.push(cell_annotation("Normal_" + (i + 1).to_string(), "Normal"))
  }

  let n_cells = n_tumour + n_normal
  let input : InferCNVInput = {
    n_cells: n_cells,
    n_genes: n_genes,
    expression: expressions,
    annotations: annotations,
    ordered_genes: og,
  }
  (input, ["Tumour"], ["Normal"])
}

// ============================================================================
// Convenience query methods on CNVResult
// ============================================================================

///|
/// Return the mean CNV score of a category (cluster / cell type), or 0.0 if unknown.
pub fn CNVResult::cluster_score(self : CNVResult, category : String) -> Double {
  for entry in self.per_cluster_cnv_score {
    let (cat, score) = entry
    if cat == category { return score }
  }
  0.0
}

///|
/// Return the list of chromosome boundaries as triples (chr, start, end).
pub fn CNVResult::chromosome_list(self : CNVResult) -> Array[(String, Int, Int)] {
  self.chromosome_boundaries.copy()
}

///|
/// Return a short one-line summary string describing dimensions and scores.
pub fn CNVResult::summary(self : CNVResult) -> String {
  let total = self.per_cell_cnv_score.length()
  let mut sum = 0.0
  for v in self.per_cell_cnv_score { sum += v }
  let mean = if total == 0 { 0.0 } else { sum / (total.to_double()) }
  let mut max_cat = ""
  let mut max_score = 0.0
  for entry in self.per_cluster_cnv_score {
    let (cat, score) = entry
    if score > max_score { max_cat = cat; max_score = score }
  }
  "CNVResult(cells=" + total.to_string() +
    ", genes=" + self.n_genes.to_string() +
    ", mean_score=" + mean.to_string() +
    ", top_cluster=" + max_cat + "=" + max_score.to_string() + ")"
}

///|
/// Return CNV matrix value at cell index i, gene index j.
pub fn CNVResult::at(self : CNVResult, i : Int, j : Int) -> Double {
  if i < 0 || i >= self.n_cells || j < 0 || j >= self.n_genes { 0.0 }
  else { self.cnv_matrix[i][j] }
}

///|
/// Classify cells as putatively tumour-like when their per-cell CNV score
/// exceeds `threshold * normal_cluster_score` (1.5x by default).
/// Returns indices of cells flagged as malignant.
pub fn CNVResult::predict_tumour_cells(
  self : CNVResult,
  normal_category : String,
  threshold_factor? : Double = 1.5,
) -> Array[Int] {
  let base = self.cluster_score(normal_category)
  let threshold = base * threshold_factor
  let out : Array[Int] = []
  for i in 0..<self.per_cell_cnv_score.length() {
    if self.per_cell_cnv_score[i] > threshold { out.push(i) }
  }
  out
}
