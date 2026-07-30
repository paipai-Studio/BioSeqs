///|
/// muscat module - Differential state analysis for single-cell RNA-seq.
/// Inspired by Bioconductor muscat package.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/muscat.html
///   - Crowell HL et al. (2020) muscat detects subpopulation-specific state transitions
///     in multi-sample single-cell transcriptomics data. Nature Communications 11:6077.

///|
/// Aggregation method for pseudo-bulk counts.
pub enum AggregationMethod {
  Sum
  Mean
  Median
}

///|
/// Create a Sum aggregation method.
pub fn aggregation_sum() -> AggregationMethod {
  AggregationMethod::Sum
}

///|
/// Create a Mean aggregation method.
pub fn aggregation_mean() -> AggregationMethod {
  AggregationMethod::Mean
}

///|
/// Create a Median aggregation method.
pub fn aggregation_median() -> AggregationMethod {
  AggregationMethod::Median
}

///|
/// DS test method type.
pub enum DSMethod {
  EdgeR
  DESeq2
  Limma
}

///|
/// Create edgeR DS method.
pub fn ds_method_edger() -> DSMethod {
  DSMethod::EdgeR
}

///|
/// Create DESeq2 DS method.
pub fn ds_method_deseq2() -> DSMethod {
  DSMethod::DESeq2
}

///|
/// Create limma DS method.
pub fn ds_method_limma() -> DSMethod {
  DSMethod::Limma
}

///|
/// A single cell in the experiment.
pub struct SingleCell {
  cell_id : String
  sample_id : String
  cluster_id : String
  group_id : String
  gene_counts : Map[String, Double]
}

///|
/// Pseudo-bulk sample aggregated from single cells.
pub struct PseudoBulk {
  sample_id : String
  cluster_id : String
  group_id : String
  n_cells : Int
  gene_counts : Map[String, Double]
}

///|
/// Differential state test result for a single gene.
pub struct DSResult {
  gene_id : String
  cluster_id : String
  log2fc : Double
  p_val : Double
  mut p_adj : Double
  mean_ctrl : Double
  mean_stim : Double
  mut significant : Bool
}

///|
/// Container for DS results across all genes/clusters.
pub struct DSResults {
  results : Array[DSResult]
  ds_method : DSMethod
  n_genes : Int
  n_clusters : Int
  n_sig : Int
  fdr_threshold : Double
}

///|
/// QC summary for a sample-cluster combination.
pub struct SampleQC {
  sample_id : String
  cluster_id : String
  n_cells : Int
  n_genes : Int
  total_counts : Double
  median_genes_per_cell : Double
}

///|
/// AggregationMethod to string.
pub fn AggregationMethod::to_string(self : AggregationMethod) -> String {
  match self {
    AggregationMethod::Sum => "sum"
    AggregationMethod::Mean => "mean"
    AggregationMethod::Median => "median"
  }
}

///|
/// DSMethod to string.
pub fn DSMethod::to_string(self : DSMethod) -> String {
  match self {
    DSMethod::EdgeR => "edgeR"
    DSMethod::DESeq2 => "DESeq2"
    DSMethod::Limma => "limma"
  }
}

///|
/// Create a SingleCell.
pub fn SingleCell::new(
  cell_id : String,
  sample_id : String,
  cluster_id : String,
  group_id : String,
) -> SingleCell {
  SingleCell::{
    cell_id,
    sample_id,
    cluster_id,
    group_id,
    gene_counts: Map::new(),
  }
}

///|
/// Set gene count for a cell.
pub fn SingleCell::set_count(self : SingleCell, gene_id : String, count : Double) -> Unit {
  self.gene_counts.set(gene_id, count)
}

///|
/// Get gene count for a cell.
pub fn SingleCell::get_count(self : SingleCell, gene_id : String) -> Double {
  self.gene_counts.get(gene_id).unwrap_or(0.0)
}

///|
/// Get total counts for a cell.
pub fn SingleCell::total_counts(self : SingleCell) -> Double {
  let mut total = 0.0
  let keys = self.gene_counts.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let v = self.gene_counts.get(keys[i]).unwrap_or(0.0)
    total = total + v
    i = i + 1
  }
  total
}

///|
/// Get number of expressed genes for a cell.
pub fn SingleCell::n_expressed(self : SingleCell) -> Int {
  let mut n = 0
  let keys = self.gene_counts.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let v = self.gene_counts.get(keys[i]).unwrap_or(0.0)
    if v > 0.0 {
      n = n + 1
    }
    i = i + 1
  }
  n
}

///|
/// Aggregate single cells into pseudo-bulk samples.
pub fn aggregate_cells(
  cells : Array[SingleCell],
  method? : AggregationMethod = AggregationMethod::Sum,
) -> Array[PseudoBulk] {
  // Group cells by (sample_id, cluster_id)
  let groups : Map[String, Array[SingleCell]] = Map::new()
  let mut i = 0
  while i < cells.length() {
    let cell = cells[i]
    let key = cell.sample_id + "::" + cell.cluster_id
    let arr = groups.get(key).unwrap_or(Array::new())
    arr.push(cell)
    groups.set(key, arr)
    i = i + 1
  }

  let result : Array[PseudoBulk] = Array::new()
  let group_keys = groups.keys().collect()
  let mut gk = 0
  while gk < group_keys.length() {
    let k = group_keys[gk]
    let v = groups.get(k).unwrap()
    let parts = k.split("::")
    let sample_id = parts.next().unwrap_or("").to_string()
    let cluster_id = parts.next().unwrap_or("").to_string()

    // Determine group_id from first cell
    let group_id = if v.length() > 0 { v[0].group_id } else { "" }

    // Collect all gene ids
    let gene_ids : Map[String, Unit] = Map::new()
    let mut idx = 0
    while idx < v.length() {
      let cell_keys = v[idx].gene_counts.keys().collect()
      let mut cki = 0
      while cki < cell_keys.length() {
        gene_ids.set(cell_keys[cki], ())
        cki = cki + 1
      }
      idx = idx + 1
    }

    let pb_counts : Map[String, Double] = Map::new()
    let gene_id_list = gene_ids.keys().collect()
    let mut ggi = 0
    while ggi < gene_id_list.length() {
      let gid = gene_id_list[ggi]
      let values : Array[Double] = Array::new()
      let mut ci = 0
      while ci < v.length() {
        let c = v[ci].gene_counts.get(gid).unwrap_or(0.0)
        values.push(c)
        ci = ci + 1
      }
      let agg = match method {
        AggregationMethod::Sum => {
          let mut s = 0.0
          let mut vi = 0
          while vi < values.length() {
            s = s + values[vi]
            vi = vi + 1
          }
          s
        }
        AggregationMethod::Mean => {
          let mut s = 0.0
          let mut vi = 0
          while vi < values.length() {
            s = s + values[vi]
            vi = vi + 1
          }
          if values.length() > 0 { s / values.length().to_double() } else { 0.0 }
        }
        AggregationMethod::Median => {
          let sorted = values.copy()
          let slen = sorted.length()
          let mut si = 0
          while si < slen - 1 {
            let mut sj = 0
            while sj < slen - 1 - si {
              if sorted[sj] > sorted[sj + 1] {
                let tmp = sorted[sj]
                sorted.set(sj, sorted[sj + 1])
                sorted.set(sj + 1, tmp)
              }
              sj = sj + 1
            }
            si = si + 1
          }
          let mid = slen / 2
          if slen % 2 == 0 {
            (sorted[mid - 1] + sorted[mid]) / 2.0
          } else {
            sorted[mid]
          }
        }
      }
      pb_counts.set(gid, agg)
      ggi = ggi + 1
    }

    result.push(PseudoBulk::{
      sample_id,
      cluster_id,
      group_id,
      n_cells: v.length(),
      gene_counts: pb_counts,
    })
    gk = gk + 1
  }

  result
}

///|
/// Get gene count from pseudo-bulk.
pub fn PseudoBulk::get_count(self : PseudoBulk, gene_id : String) -> Double {
  self.gene_counts.get(gene_id).unwrap_or(0.0)
}

///|
/// Get total counts from pseudo-bulk.
pub fn PseudoBulk::total_counts(self : PseudoBulk) -> Double {
  let mut total = 0.0
  let keys = self.gene_counts.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let v = self.gene_counts.get(keys[i]).unwrap_or(0.0)
    total = total + v
    i = i + 1
  }
  total
}

///|
/// Get number of expressed genes from pseudo-bulk.
pub fn PseudoBulk::n_expressed(self : PseudoBulk) -> Int {
  let mut n = 0
  let keys = self.gene_counts.keys().collect()
  let mut i = 0
  while i < keys.length() {
    let v = self.gene_counts.get(keys[i]).unwrap_or(0.0)
    if v > 0.0 {
      n = n + 1
    }
    i = i + 1
  }
  n
}

///|
/// Run differential state analysis on pseudo-bulk data.
pub fn run_ds_analysis(
  pb_samples : Array[PseudoBulk],
  method? : DSMethod = DSMethod::EdgeR,
  fdr_threshold? : Double = 0.05,
  log2fc_threshold? : Double = 1.0,
) -> DSResults {
  // Group by cluster_id
  let cluster_groups : Map[String, Array[PseudoBulk]] = Map::new()
  let mut i = 0
  while i < pb_samples.length() {
    let pb = pb_samples[i]
    let arr = cluster_groups.get(pb.cluster_id).unwrap_or(Array::new())
    arr.push(pb)
    cluster_groups.set(pb.cluster_id, arr)
    i = i + 1
  }

  let all_results : Array[DSResult] = Array::new()
  let n_clusters = cluster_groups.size()

  let cluster_keys = cluster_groups.keys().collect()
  let mut cki = 0
  while cki < cluster_keys.length() {
    let cluster_id = cluster_keys[cki]
    let samples = cluster_groups.get(cluster_id).unwrap()
    // Split into control and stimulation groups
    let ctrl : Array[PseudoBulk] = Array::new()
    let stim : Array[PseudoBulk] = Array::new()
    let mut si = 0
    while si < samples.length() {
      if samples[si].group_id == "ctrl" || samples[si].group_id == "control" {
        ctrl.push(samples[si])
      } else {
        stim.push(samples[si])
      }
      si = si + 1
    }

    // Collect all gene ids
    let gene_ids : Map[String, Unit] = Map::new()
    let mut ci = 0
    while ci < samples.length() {
      let sck = samples[ci].gene_counts.keys().collect()
      let mut scki = 0
      while scki < sck.length() {
        gene_ids.set(sck[scki], ())
        scki = scki + 1
      }
      ci = ci + 1
    }

    let gene_id_list = gene_ids.keys().collect()
    let mut ggi = 0
    while ggi < gene_id_list.length() {
      let gid = gene_id_list[ggi]
      // Get counts for control and stim
      let ctrl_counts : Array[Double] = Array::new()
      let mut cci = 0
      while cci < ctrl.length() {
        ctrl_counts.push(ctrl[cci].gene_counts.get(gid).unwrap_or(0.0))
        cci = cci + 1
      }
      let stim_counts : Array[Double] = Array::new()
      let mut sci = 0
      while sci < stim.length() {
        stim_counts.push(stim[sci].gene_counts.get(gid).unwrap_or(0.0))
        sci = sci + 1
      }

      let mean_ctrl = compute_mean(ctrl_counts)
      let mean_stim = compute_mean(stim_counts)

      // Compute log2 fold change
      let log2fc = if mean_ctrl > 0.0 && mean_stim > 0.0 {
        @math.log2(mean_stim) - @math.log2(mean_ctrl)
      } else if mean_ctrl > 0.0 {
        -@math.log2(mean_ctrl + 1.0)
      } else if mean_stim > 0.0 {
        @math.log2(mean_stim + 1.0)
      } else {
        0.0
      }

      // Simplified p-value computation using t-test approximation
      let p_val = compute_ttest_pvalue(ctrl_counts, stim_counts)

      // BH correction will be applied globally later
      all_results.push(DSResult::{
        gene_id: gid,
        cluster_id,
        log2fc,
        p_val,
        p_adj: p_val,
        mean_ctrl,
        mean_stim,
        significant: false,
      })
      ggi = ggi + 1
    }
    cki = cki + 1
  }

  // Apply Benjamini-Hochberg correction
  let p_values : Array[Double] = Array::new()
  let mut ri = 0
  while ri < all_results.length() {
    p_values.push(all_results[ri].p_val)
    ri = ri + 1
  }
  let adjusted = muscat_bh_adjust(p_values)

  let n_sig = 0
  let mut si2 = 0
  while si2 < all_results.length() {
    all_results[si2].p_adj = adjusted[si2]
    let abs_lfc = if all_results[si2].log2fc < 0.0 { -all_results[si2].log2fc } else { all_results[si2].log2fc }
    let sig = adjusted[si2] < fdr_threshold && abs_lfc >= log2fc_threshold
    all_results[si2].significant = sig
    si2 = si2 + 1
  }

  let mut n_sig_count = 0
  let mut si3 = 0
  while si3 < all_results.length() {
    if all_results[si3].significant {
      n_sig_count = n_sig_count + 1
    }
    si3 = si3 + 1
  }

  let n_genes = if all_results.length() > 0 {
    let gene_set : Map[String, Unit] = Map::new()
    let mut gi = 0
    while gi < all_results.length() {
      gene_set.set(all_results[gi].gene_id, ())
      gi = gi + 1
    }
    gene_set.size()
  } else {
    0
  }

  DSResults::{
    results: all_results,
    ds_method: method,
    n_genes,
    n_clusters,
    n_sig: n_sig_count,
    fdr_threshold,
  }
}

///|
/// Compute mean of an array.
fn compute_mean(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < values.length() {
    sum = sum + values[i]
    i = i + 1
  }
  sum / values.length().to_double()
}

///|
/// Compute variance of an array.
fn compute_variance(values : Array[Double]) -> Double {
  if values.length() < 2 {
    return 0.0
  }
  let mean = compute_mean(values)
  let mut sum_sq = 0.0
  let mut i = 0
  while i < values.length() {
    let diff = values[i] - mean
    sum_sq = sum_sq + diff * diff
    i = i + 1
  }
  sum_sq / (values.length() - 1).to_double()
}

///|
/// Compute a simplified t-test p-value (two-sided).
fn compute_ttest_pvalue(ctrl : Array[Double], stim : Array[Double]) -> Double {
  let n1 = ctrl.length()
  let n2 = stim.length()
  if n1 < 2 || n2 < 2 {
    return 1.0
  }

  let m1 = compute_mean(ctrl)
  let m2 = compute_mean(stim)
  let v1 = compute_variance(ctrl)
  let v2 = compute_variance(stim)

  // Pooled standard error
  let pooled_se = (v1 / n1.to_double() + v2 / n2.to_double()).sqrt()
  if pooled_se < 1.0e-10 {
    return if (m1 - m2).abs() < 1.0e-10 { 1.0 } else { 0.001 }
  }

  let t_stat = (m2 - m1) / pooled_se

  // Degrees of freedom (Welch's)
  let df_num = (v1 / n1.to_double() + v2 / n2.to_double()) * (v1 / n1.to_double() + v2 / n2.to_double())
  let df_den = (v1 / n1.to_double()) * (v1 / n1.to_double()) / (n1 - 1).to_double() + (v2 / n2.to_double()) * (v2 / n2.to_double()) / (n2 - 1).to_double()
  let df = if df_den < 1.0e-10 { (n1 + n2 - 2).to_double() } else { df_num / df_den }

  // Approximate two-sided p-value from t-distribution
  // Using a simple approximation: p = 2 * (1 - cdf(|t|, df))
  // Approximation: p ≈ 2 * exp(-df * ln(1 + t^2/df) / 2) / sqrt(df * pi)
  let abs_t = if t_stat < 0.0 { -t_stat } else { t_stat }
  let t_sq = abs_t * abs_t
  let p_approx = 2.0 * @math.exp(-df * @math.ln(1.0 + t_sq / df) / 2.0) / (df * 3.14159265358979).sqrt()

  // Clamp p-value to [0, 1]
  if p_approx > 1.0 { 1.0 } else if p_approx < 0.0 { 0.0 } else { p_approx }
}

///|
/// Benjamini-Hochberg FDR adjustment.
fn muscat_bh_adjust(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return Array::new()
  }

  // Create index-value pairs and sort by p-value
  let indexed : Array[(Int, Double)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((i, p_values[i]))
    i = i + 1
  }

  // Sort by p-value (bubble sort)
  let mut j = 0
  while j < n - 1 {
    let mut k = 0
    while k < n - 1 - j {
      if indexed[k].1 > indexed[k + 1].1 {
        let tmp = indexed[k]
        indexed.set(k, indexed[k + 1])
        indexed.set(k + 1, tmp)
      }
      k = k + 1
    }
    j = j + 1
  }

  // Compute adjusted p-values
  let adjusted : Array[Double] = Array::new()
  let mut ai = 0
  while ai < n {
    adjusted.push(0.0)
    ai = ai + 1
  }

  let mut prev_adj = 1.0
  let mut ri = n - 1
  while ri >= 0 {
    let rank = ri + 1
    let raw_p = indexed[ri].1
    let bh_val = raw_p * n.to_double() / rank.to_double()
    if bh_val < prev_adj {
      prev_adj = bh_val
    }
    if prev_adj > 1.0 {
      adjusted.set(indexed[ri].0, 1.0)
    } else {
      adjusted.set(indexed[ri].0, prev_adj)
    }
    ri = ri - 1
  }

  adjusted
}

///|
/// Get number of results.
pub fn DSResults::get_n_results(self : DSResults) -> Int {
  self.results.length()
}

///|
/// Get significant results only.
pub fn DSResults::get_significant(self : DSResults) -> Array[DSResult] {
  let sig : Array[DSResult] = Array::new()
  let mut i = 0
  while i < self.results.length() {
    if self.results[i].significant {
      sig.push(self.results[i])
    }
    i = i + 1
  }
  sig
}

///|
/// Get results for a specific cluster.
pub fn DSResults::get_cluster_results(self : DSResults, cluster_id : String) -> Array[DSResult] {
  let filtered : Array[DSResult] = Array::new()
  let mut i = 0
  while i < self.results.length() {
    if self.results[i].cluster_id == cluster_id {
      filtered.push(self.results[i])
    }
    i = i + 1
  }
  filtered
}

///|
/// Get top N genes by adjusted p-value.
pub fn DSResults::get_top_genes(self : DSResults, n : Int) -> Array[DSResult] {
  let sorted = self.results.copy()
  let len = sorted.length()
  let mut i = 0
  while i < len - 1 {
    let mut j = 0
    while j < len - 1 - i {
      if sorted[j].p_adj > sorted[j + 1].p_adj {
        let tmp = sorted[j]
        sorted.set(j, sorted[j + 1])
        sorted.set(j + 1, tmp)
      }
      j = j + 1
    }
    i = i + 1
  }

  let result : Array[DSResult] = Array::new()
  let take = if n < len { n } else { len }
  let mut k = 0
  while k < take {
    result.push(sorted[k])
    k = k + 1
  }
  result
}

///|
/// Generate a summary string.
pub fn DSResults::summary(self : DSResults) -> String {
  let mut s = "=== muscat DS Analysis Summary ===\n".to_string()
  s = s + "Method: " + self.ds_method.to_string() + "\n"
  s = s + "Total genes tested: " + self.n_genes.to_string() + "\n"
  s = s + "Clusters: " + self.n_clusters.to_string() + "\n"
  s = s + "Significant (FDR < " + self.fdr_threshold.to_string() + "): " + self.n_sig.to_string() + "\n"
  s = s + "Total results: " + self.results.length().to_string() + "\n"
  s
}

///|
/// Compute QC summary for pseudo-bulk samples.
pub fn compute_qc(pb_samples : Array[PseudoBulk]) -> Array[SampleQC] {
  let qc_list : Array[SampleQC] = Array::new()
  let mut i = 0
  while i < pb_samples.length() {
    let pb = pb_samples[i]
    let n_genes = pb.n_expressed()
    let total = pb.total_counts()
    let median_genes = if pb.n_cells > 0 {
      n_genes.to_double()
    } else {
      0.0
    }
    qc_list.push(SampleQC::{
      sample_id: pb.sample_id,
      cluster_id: pb.cluster_id,
      n_cells: pb.n_cells,
      n_genes,
      total_counts: total,
      median_genes_per_cell: median_genes,
    })
    i = i + 1
  }
  qc_list
}

///|
/// Create sample single-cell data for testing.
pub fn muscat_sample_data() -> Array[SingleCell] {
  let cells : Array[SingleCell] = Array::new()
  let genes = ["Gene1", "Gene2", "Gene3", "Gene4", "Gene5"]
  let samples = ["S1", "S2", "S3", "S4"]
  let clusters = ["C1", "C2"]
  let groups = ["ctrl", "ctrl", "stim", "stim"]

  let mut si = 0
  while si < samples.length() {
    let mut ci = 0
    while ci < clusters.length() {
      let mut cell_idx = 0
      while cell_idx < 5 {
        let cell = SingleCell::new(
          samples[si] + "_" + clusters[ci] + "_" + cell_idx.to_string(),
          samples[si],
          clusters[ci],
          groups[si],
        )
        let mut gi = 0
        while gi < genes.length() {
          let count = if groups[si] == "stim" {
            (cell_idx + 1).to_double() * (gi + 1).to_double() * 2.0
          } else {
            (cell_idx + 1).to_double() * (gi + 1).to_double()
          }
          cell.set_count(genes[gi], count)
          gi = gi + 1
        }
        cells.push(cell)
        cell_idx = cell_idx + 1
      }
      ci = ci + 1
    }
    si = si + 1
  }
  cells
}
