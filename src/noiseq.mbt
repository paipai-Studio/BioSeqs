///|
/// NOISeq module - Noise-robust differential expression for RNA-seq.
/// Inspired by Bioconductor NOISeq package.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/NOISeq.html
///   - Tarazona S et al. (2015) Data quality aware analysis of differential expression
///     in RNA-seq with NOISeq Bioconductor package. Nucleic Acids Research 43(21):e140.

///|
/// Normalization method for NOISeq.
pub enum NOISeqNorm {
  RPKM
  TMM
  UpperQuartile
  None
}

///|
/// Create RPKM normalization.
pub fn noiseq_norm_rpkm() -> NOISeqNorm {
  NOISeqNorm::RPKM
}

///|
/// Create TMM normalization.
pub fn noiseq_norm_tmm() -> NOISeqNorm {
  NOISeqNorm::TMM
}

///|
/// Create Upper Quartile normalization.
pub fn noiseq_norm_upper() -> NOISeqNorm {
  NOISeqNorm::UpperQuartile
}

///|
/// Create no normalization.
pub fn noiseq_norm_none() -> NOISeqNorm {
  NOISeqNorm::None
}

///|
/// NOISeq method type.
pub enum NOISeqMethod {
  NOISeqBio
  NOISeqSim
}

///|
/// Create NOISeqBio method.
pub fn noiseq_method_bio() -> NOISeqMethod {
  NOISeqMethod::NOISeqBio
}

///|
/// Create NOISeqSim method.
pub fn noiseq_method_sim() -> NOISeqMethod {
  NOISeqMethod::NOISeqSim
}

///|
/// A gene expression sample.
pub struct NOISeqSample {
  sample_id : String
  condition : String
  counts : Map[String, Double]
}

///|
/// Differential expression result for a gene.
pub struct NOISeqResult {
  gene_id : String
  mean_control : Double
  mean_treatment : Double
  log2fc : Double
  divergence : Double
  prob : Double
  significant : Bool
  mut ranking : Int
}

///|
/// Container for NOISeq results.
pub struct NOISeqResults {
  results : Array[NOISeqResult]
  ds_method : NOISeqMethod
  n_genes : Int
  n_sig : Int
  prob_threshold : Double
  norm_method : NOISeqNorm
}

///|
/// QC diagnostic result.
pub struct NOISeqQC {
  biotype_counts : Map[String, Int]
  gc_content : Double
  length_distribution : Map[String, Int]
  bias_detected : Bool
}

///|
/// NOISeqNorm to string.
pub fn NOISeqNorm::to_string(self : NOISeqNorm) -> String {
  match self {
    NOISeqNorm::RPKM => "RPKM"
    NOISeqNorm::TMM => "TMM"
    NOISeqNorm::UpperQuartile => "UpperQuartile"
    NOISeqNorm::None => "none"
  }
}

///|
/// NOISeqMethod to string.
pub fn NOISeqMethod::to_string(self : NOISeqMethod) -> String {
  match self {
    NOISeqMethod::NOISeqBio => "NOISeqBio"
    NOISeqMethod::NOISeqSim => "NOISeqSim"
  }
}

///|
/// Create a NOISeqSample.
pub fn NOISeqSample::new(sample_id : String, condition : String) -> NOISeqSample {
  NOISeqSample::{ sample_id, condition, counts: Map::new() }
}

///|
/// Set gene count for a sample.
pub fn NOISeqSample::set_count(self : NOISeqSample, gene_id : String, count : Double) -> Unit {
  self.counts.set(gene_id, count)
}

///|
/// Get gene count from a sample.
pub fn NOISeqSample::get_count(self : NOISeqSample, gene_id : String) -> Double {
  self.counts.get(gene_id).unwrap_or(0.0)
}

///|
/// Get total counts (library size) for a sample.
pub fn NOISeqSample::library_size(self : NOISeqSample) -> Double {
  let mut total = 0.0
  let key_list = self.counts.keys().collect()
  let mut i = 0
  while i < key_list.length() {
    let v = self.counts.get(key_list[i]).unwrap_or(0.0)
    total = total + v
    i = i + 1
  }
  total
}

///|
/// Get number of expressed genes.
pub fn NOISeqSample::n_expressed(self : NOISeqSample) -> Int {
  let mut n = 0
  let key_list = self.counts.keys().collect()
  let mut i = 0
  while i < key_list.length() {
    let v = self.counts.get(key_list[i]).unwrap_or(0.0)
    if v > 0.0 {
      n = n + 1
    }
    i = i + 1
  }
  n
}

///|
/// Normalize counts using specified method.
pub fn noiseq_normalize(
  samples : Array[NOISeqSample],
  norm_method? : NOISeqNorm = NOISeqNorm::TMM,
) -> Array[NOISeqSample] {
  match norm_method {
    NOISeqNorm::None => samples
    NOISeqNorm::RPKM => {
      // RPKM: Reads Per Kilobase Million (simplified: just scale to library size)
      let normalized : Array[NOISeqSample] = Array::new()
      let mut i = 0
      while i < samples.length() {
        let sample = samples[i]
        let norm_sample = NOISeqSample::new(sample.sample_id, sample.condition)
        let lib_size = sample.library_size()
        let scale = if lib_size > 0.0 { 1.0e6 / lib_size } else { 1.0 }
        let gid_list = sample.counts.keys().collect()
        let mut gi = 0
        while gi < gid_list.length() {
          let gid = gid_list[gi]
          let count = sample.counts.get(gid).unwrap_or(0.0)
          // Simplified RPKM: CPM (counts per million)
          norm_sample.set_count(gid, count * scale)
          gi = gi + 1
        }
        normalized.push(norm_sample)
        i = i + 1
      }
      normalized
    }
    NOISeqNorm::TMM => {
      // TMM: Trimmed Mean of M-values (simplified version)
      // Compute normalization factors relative to first sample
      let normalized : Array[NOISeqSample] = Array::new()
      if samples.length() == 0 {
        return normalized
      }

      let ref_sample = samples[0]
      let ref_lib = ref_sample.library_size()

      let mut i = 0
      while i < samples.length() {
        let sample = samples[i]
        let norm_sample = NOISeqSample::new(sample.sample_id, sample.condition)
        let lib_size = sample.library_size()

        // Compute simple scaling factor
        let scale = if lib_size > 0.0 && ref_lib > 0.0 { ref_lib / lib_size } else { 1.0 }

        let gid_list = sample.counts.keys().collect()
        let mut gi = 0
        while gi < gid_list.length() {
          let gid = gid_list[gi]
          let count = sample.counts.get(gid).unwrap_or(0.0)
          norm_sample.set_count(gid, count * scale)
          gi = gi + 1
        }
        normalized.push(norm_sample)
        i = i + 1
      }
      normalized
    }
    NOISeqNorm::UpperQuartile => {
      // Upper Quartile normalization
      let normalized : Array[NOISeqSample] = Array::new()
      let mut i = 0
      while i < samples.length() {
        let sample = samples[i]
        let norm_sample = NOISeqSample::new(sample.sample_id, sample.condition)

        // Collect non-zero counts
        let counts : Array[Double] = Array::new()
        let key_list = sample.counts.keys().collect()
        let mut ki = 0
        while ki < key_list.length() {
          let v = sample.counts.get(key_list[ki]).unwrap_or(0.0)
          if v > 0.0 {
            counts.push(v)
          }
          ki = ki + 1
        }

        // Sort and find upper quartile
        let slen = counts.length()
        let mut si = 0
        while si < slen - 1 {
          let mut sj = 0
          while sj < slen - 1 - si {
            if counts[sj] > counts[sj + 1] {
              let tmp = counts[sj]
              counts.set(sj, counts[sj + 1])
              counts.set(sj + 1, tmp)
            }
            sj = sj + 1
          }
          si = si + 1
        }

        let uq_pos = (slen * 3) / 4
        let upper_quartile = if slen > 0 && uq_pos < slen {
          counts[uq_pos]
        } else {
          1.0
        }

        let scale = if upper_quartile > 0.0 { 1.0e6 / upper_quartile } else { 1.0 }

        let gid_list = sample.counts.keys().collect()
        let mut gi = 0
        while gi < gid_list.length() {
          let gid = gid_list[gi]
          let count = sample.counts.get(gid).unwrap_or(0.0)
          norm_sample.set_count(gid, count * scale)
          gi = gi + 1
        }
        normalized.push(norm_sample)
        i = i + 1
      }
      normalized
    }
  }
}

///|
/// Run NOISeq differential expression analysis.
pub fn noiseq_run(
  control_samples : Array[NOISeqSample],
  treatment_samples : Array[NOISeqSample],
  ds_method? : NOISeqMethod = NOISeqMethod::NOISeqBio,
  norm? : NOISeqNorm = NOISeqNorm::TMM,
  prob_threshold? : Double = 0.8,
  log2fc_threshold? : Double = 0.0,
) -> NOISeqResults {
  // Normalize
  let norm_ctrl = noiseq_normalize(control_samples, norm_method=norm)
  let norm_trt = noiseq_normalize(treatment_samples, norm_method=norm)

  // Collect all gene ids
  let gene_ids : Map[String, Unit] = Map::new()
  let mut ci = 0
  while ci < norm_ctrl.length() {
    let ctrl_gid_list = norm_ctrl[ci].counts.keys().collect()
    let mut cgi = 0
    while cgi < ctrl_gid_list.length() {
      gene_ids.set(ctrl_gid_list[cgi], ())
      cgi = cgi + 1
    }
    ci = ci + 1
  }
  let mut ti = 0
  while ti < norm_trt.length() {
    let trt_gid_list = norm_trt[ti].counts.keys().collect()
    let mut tgi = 0
    while tgi < trt_gid_list.length() {
      gene_ids.set(trt_gid_list[tgi], ())
      tgi = tgi + 1
    }
    ti = ti + 1
  }

  let results : Array[NOISeqResult] = Array::new()

  let gid_list = gene_ids.keys().collect()
  let mut gi = 0
  while gi < gid_list.length() {
    let gid = gid_list[gi]
    // Compute mean expression per condition
    let ctrl_vals : Array[Double] = Array::new()
    let mut ci2 = 0
    while ci2 < norm_ctrl.length() {
      ctrl_vals.push(norm_ctrl[ci2].get_count(gid))
      ci2 = ci2 + 1
    }
    let trt_vals : Array[Double] = Array::new()
    let mut ti2 = 0
    while ti2 < norm_trt.length() {
      trt_vals.push(norm_trt[ti2].get_count(gid))
      ti2 = ti2 + 1
    }

    let mean_ctrl = noiseq_mean(ctrl_vals)
    let mean_trt = noiseq_mean(trt_vals)

    // Compute M (log2 fold change) and D (divergence)
    let m_val = if mean_ctrl > 0.0 && mean_trt > 0.0 {
      @math.log2(mean_trt) - @math.log2(mean_ctrl)
    } else if mean_ctrl > 0.0 {
      -@math.log2(mean_ctrl + 1.0)
    } else if mean_trt > 0.0 {
      @math.log2(mean_trt + 1.0)
    } else {
      0.0
    }

    // Compute D: average of absolute differences within conditions
    let d_ctrl = noiseq_mean_abs_diff(ctrl_vals)
    let d_trt = noiseq_mean_abs_diff(trt_vals)
    let d_val = (d_ctrl + d_trt) / 2.0

    // Compute probability using noise distribution
    // In NOISeq, the probability is based on the proportion of noise
    // genes with (M, D) values less extreme than the observed value
    let prob = noiseq_compute_prob(m_val, d_val, ds_method)

    let abs_m = if m_val < 0.0 { -m_val } else { m_val }
    let significant = prob >= prob_threshold && abs_m >= log2fc_threshold

    results.push(NOISeqResult::{
      gene_id: gid,
      mean_control: mean_ctrl,
      mean_treatment: mean_trt,
      log2fc: m_val,
      divergence: d_val,
      prob,
      significant,
      ranking: 0,
    })
    gi = gi + 1
  }

  // Sort results by probability (descending) and assign rankings
  let sorted = results.copy()
  let slen = sorted.length()
  let mut si = 0
  while si < slen - 1 {
    let mut sj = 0
    while sj < slen - 1 - si {
      if sorted[sj].prob < sorted[sj + 1].prob {
        let tmp = sorted[sj]
        sorted.set(sj, sorted[sj + 1])
        sorted.set(sj + 1, tmp)
      }
      sj = sj + 1
    }
    si = si + 1
  }

  // Assign rankings
  let mut ri = 0
  while ri < slen {
    sorted[ri].ranking = ri + 1
    ri = ri + 1
  }

  // Count significant
  let mut n_sig = 0
  let mut si2 = 0
  while si2 < slen {
    if sorted[si2].significant {
      n_sig = n_sig + 1
    }
    si2 = si2 + 1
  }

  NOISeqResults::{
    results: sorted,
    ds_method,
    n_genes: gene_ids.length(),
    n_sig,
    prob_threshold,
    norm_method: norm,
  }
}

///|
/// Compute mean.
fn noiseq_mean(values : Array[Double]) -> Double {
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
/// Compute mean absolute difference (D statistic).
fn noiseq_mean_abs_diff(values : Array[Double]) -> Double {
  if values.length() < 2 {
    return 0.0
  }
  let mean = noiseq_mean(values)
  let mut sum_abs = 0.0
  let mut i = 0
  while i < values.length() {
    let diff = values[i] - mean
    sum_abs = sum_abs + (if diff < 0.0 { -diff } else { diff })
    i = i + 1
  }
  sum_abs / values.length().to_double()
}

///|
/// Compute probability from noise distribution.
/// In NOISeq, this is the probability that a random gene from the noise
/// distribution has a (M, D) value more extreme than the observed.
fn noiseq_compute_prob(m : Double, d : Double, ds_method : NOISeqMethod) -> Double {
  // The magnitude of the (M, D) vector indicates how far the gene is from the noise
  let magnitude = (m * m + d * d).sqrt()

  // Higher magnitude -> higher probability of being truly differentially expressed
  // This is a simplified version of the actual NOISeq algorithm
  // The actual algorithm uses bootstrapping/resampling to build the noise distribution
  match ds_method {
    NOISeqMethod::NOISeqBio => {
      // NOISeqBio uses a more conservative estimate
      // probability ≈ 1 - exp(-magnitude / sigma)
      let sigma = 1.0
      let prob = 1.0 - @math.exp(-magnitude / sigma)
      if prob > 1.0 { 1.0 } else { prob }
    }
    NOISeqMethod::NOISeqSim => {
      // NOISeqSim is less conservative
      let sigma = 0.8
      let prob = 1.0 - @math.exp(-magnitude / sigma)
      if prob > 1.0 { 1.0 } else { prob }
    }
  }
}

///|
/// Get number of results.
pub fn NOISeqResults::get_n_results(self : NOISeqResults) -> Int {
  self.results.length()
}

///|
/// Get significant results.
pub fn NOISeqResults::get_significant(self : NOISeqResults) -> Array[NOISeqResult] {
  let sig : Array[NOISeqResult] = Array::new()
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
/// Get top N genes by probability.
pub fn NOISeqResults::get_top_genes(self : NOISeqResults, n : Int) -> Array[NOISeqResult] {
  let result : Array[NOISeqResult] = Array::new()
  let take = if n < self.results.length() { n } else { self.results.length() }
  let mut i = 0
  while i < take {
    result.push(self.results[i])
    i = i + 1
  }
  result
}

///|
/// Get up-regulated genes.
pub fn NOISeqResults::get_up_regulated(self : NOISeqResults) -> Array[NOISeqResult] {
  let result : Array[NOISeqResult] = Array::new()
  let mut i = 0
  while i < self.results.length() {
    if self.results[i].significant && self.results[i].log2fc > 0.0 {
      result.push(self.results[i])
    }
    i = i + 1
  }
  result
}

///|
/// Get down-regulated genes.
pub fn NOISeqResults::get_down_regulated(self : NOISeqResults) -> Array[NOISeqResult] {
  let result : Array[NOISeqResult] = Array::new()
  let mut i = 0
  while i < self.results.length() {
    if self.results[i].significant && self.results[i].log2fc < 0.0 {
      result.push(self.results[i])
    }
    i = i + 1
  }
  result
}

///|
/// Generate summary string.
pub fn NOISeqResults::summary(self : NOISeqResults) -> String {
  let mut s = "=== NOISeq Results Summary ===\n".to_string()
  s = s + "Method: " + self.ds_method.to_string() + "\n"
  s = s + "Normalization: " + self.norm_method.to_string() + "\n"
  s = s + "Genes tested: " + self.n_genes.to_string() + "\n"
  s = s + "Significant (prob >= " + self.prob_threshold.to_string() + "): " + self.n_sig.to_string() + "\n"
  let n_up = self.get_up_regulated().length()
  let n_down = self.get_down_regulated().length()
  s = s + "  Up-regulated: " + n_up.to_string() + "\n"
  s = s + "  Down-regulated: " + n_down.to_string() + "\n"
  s
}

///|
/// Compute QC diagnostics for the data.
pub fn noiseq_qc(samples : Array[NOISeqSample]) -> NOISeqQC {
  let biotype_counts : Map[String, Int] = Map::new()
  let gc_content = 50.0 // Simplified: assume 50% GC on average
  let length_dist : Map[String, Int] = Map::new()

  // Count genes per biotype (simplified: all genes are "protein_coding")
  let mut total_genes = 0
  let mut i = 0
  while i < samples.length() {
    total_genes = total_genes + samples[i].n_expressed()
    i = i + 1
  }

  biotype_counts.set("protein_coding", total_genes)

  // Length distribution (simplified)
  length_dist.set("0-500", 30)
  length_dist.set("500-1000", 25)
  length_dist.set("1000-2000", 30)
  length_dist.set("2000-5000", 10)
  length_dist.set("5000+", 5)

  // Bias detection (simplified: check if library sizes are very different)
  let lib_sizes : Array[Double] = Array::new()
  let mut li = 0
  while li < samples.length() {
    lib_sizes.push(samples[li].library_size())
    li = li + 1
  }

  let mean_lib = noiseq_mean(lib_sizes)
  let max_diff = {
    let mut max_d = 0.0
    let mut mi = 0
    while mi < lib_sizes.length() {
      let diff = (lib_sizes[mi] - mean_lib).abs() / (mean_lib + 1.0e-10)
      if diff > max_d { max_d = diff }
      mi = mi + 1
    }
    max_d
  }

  let bias_detected = max_diff > 0.5 // More than 50% deviation

  NOISeqQC::{
    biotype_counts,
    gc_content,
    length_distribution: length_dist,
    bias_detected,
  }
}

///|
/// Create sample data for testing.
pub fn noiseq_sample_data() -> (Array[NOISeqSample], Array[NOISeqSample]) {
  let genes = ["Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6", "Gene7", "Gene8", "Gene9", "Gene10"]

  // Control samples
  let ctrl_samples : Array[NOISeqSample] = Array::new()
  let mut si = 0
  while si < 3 {
    let sample = NOISeqSample::new("ctrl_" + si.to_string(), "control")
    let mut gi = 0
    while gi < genes.length() {
      let count = (gi + 1).to_double() * 100.0 * (1.0 + si.to_double() * 0.1)
      sample.set_count(genes[gi], count)
      gi = gi + 1
    }
    ctrl_samples.push(sample)
    si = si + 1
  }

  // Treatment samples (with differential expression)
  let trt_samples : Array[NOISeqSample] = Array::new()
  let mut si2 = 0
  while si2 < 3 {
    let sample = NOISeqSample::new("trt_" + si2.to_string(), "treatment")
    let mut gi = 0
    while gi < genes.length() {
      // First 3 genes are upregulated, next 3 down, rest unchanged
      let count = if gi < 3 {
        (gi + 1).to_double() * 100.0 * 4.0 * (1.0 + si2.to_double() * 0.1)
      } else if gi < 6 {
        (gi + 1).to_double() * 100.0 * 0.2 * (1.0 + si2.to_double() * 0.1)
      } else {
        (gi + 1).to_double() * 100.0 * (1.0 + si2.to_double() * 0.1)
      }
      sample.set_count(genes[gi], count)
      gi = gi + 1
    }
    trt_samples.push(sample)
    si2 = si2 + 1
  }

  (ctrl_samples, trt_samples)
}
