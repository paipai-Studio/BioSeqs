///|
/// MSstats module - Protein significance analysis for mass spectrometry data.
/// Inspired by Bioconductor MSstats package.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/MSstats.html
///   - Choi M et al. (2014) MSstats: an R package for statistical analysis of
///     quantitative proteomics experiments. Bioinformatics 30:2524-2526.

///|
/// MS experiment type.
pub enum MSType {
  DDA
  DIA
  SRM
  TMT
}

///|
/// Create DDA experiment type.
pub fn ms_type_dda() -> MSType {
  MSType::DDA
}

///|
/// Create DIA experiment type.
pub fn ms_type_dia() -> MSType {
  MSType::DIA
}

///|
/// Create SRM experiment type.
pub fn ms_type_srm() -> MSType {
  MSType::SRM
}

///|
/// Create TMT experiment type.
pub fn ms_type_tmt() -> MSType {
  MSType::TMT
}

///|
/// Normalization method.
pub enum MSNormalization {
  None
  Median
  Quantile
  GlobalStandards
}

///|
/// Create no normalization.
pub fn ms_norm_none() -> MSNormalization {
  MSNormalization::None
}

///|
/// Create median normalization.
pub fn ms_norm_median() -> MSNormalization {
  MSNormalization::Median
}

///|
/// Create quantile normalization.
pub fn ms_norm_quantile() -> MSNormalization {
  MSNormalization::Quantile
}

///|
/// Create global standards normalization.
pub fn ms_norm_global() -> MSNormalization {
  MSNormalization::GlobalStandards
}

///|
/// Summarization method for peptides to protein.
pub enum MSSummarization {
  Tukey
  Linear
  LogSum
}

///|
/// Create Tukey median polish summarization.
pub fn ms_summary_tukey() -> MSSummarization {
  MSSummarization::Tukey
}

///|
/// Create linear model summarization.
pub fn ms_summary_linear() -> MSSummarization {
  MSSummarization::Linear
}

///|
/// Create log-sum summarization.
pub fn ms_summary_logsum() -> MSSummarization {
  MSSummarization::LogSum
}

///|
/// A single feature (peptide) measurement.
pub struct MSFeature {
  protein : String
  peptide : String
  transition : String
  group : String
  subject : String
  run : String
  mut intensity : Double
}

///|
/// Protein-level summarized data.
pub struct MSProteinAbundance {
  protein : String
  group : String
  run : String
  log2_abundance : Double
}

///|
/// Group comparison result.
pub struct MSGroupResult {
  protein : String
  comparison : String
  log2fc : Double
  se : Double
  t_stat : Double
  p_val : Double
  mut p_adj : Double
  mut significant : Bool
}

///|
/// Container for group comparison results.
pub struct MSGroupComparison {
  results : Array[MSGroupResult]
  n_proteins : Int
  n_sig : Int
  fdr_threshold : Double
}

///|
/// Sample size calculation result.
pub struct MSSampleSize {
  desired_power : Double
  n_samples : Int
  n_proteins_detectable : Int
  delta : Double
}

///|
/// MSType to string.
pub fn MSType::to_string(self : MSType) -> String {
  match self {
    MSType::DDA => "DDA"
    MSType::DIA => "DIA"
    MSType::SRM => "SRM"
    MSType::TMT => "TMT"
  }
}

///|
/// MSNormalization to string.
pub fn MSNormalization::to_string(self : MSNormalization) -> String {
  match self {
    MSNormalization::None => "none"
    MSNormalization::Median => "median"
    MSNormalization::Quantile => "quantile"
    MSNormalization::GlobalStandards => "globalStandards"
  }
}

///|
/// MSSummarization to string.
pub fn MSSummarization::to_string(self : MSSummarization) -> String {
  match self {
    MSSummarization::Tukey => "Tukey"
    MSSummarization::Linear => "linear"
    MSSummarization::LogSum => "logSum"
  }
}

///|
/// Create a MSFeature.
pub fn MSFeature::new(
  protein : String,
  peptide : String,
  transition : String,
  group : String,
  subject : String,
  run : String,
  intensity : Double,
) -> MSFeature {
  MSFeature::{ protein, peptide, transition, group, subject, run, intensity }
}

///|
/// Process MS data: filter, log-transform, normalize.
pub fn ms_data_process(
  features : Array[MSFeature],
  normalization? : MSNormalization = MSNormalization::Median,
  log_transformed? : Bool = false,
) -> Array[MSFeature] {
  // Deep copy
  let processed : Array[MSFeature] = Array::new()
  let mut i = 0
  while i < features.length() {
    let f = features[i]
    processed.push(MSFeature::new(f.protein, f.peptide, f.transition, f.group, f.subject, f.run, f.intensity))
    i = i + 1
  }

  // Filter out zero/negative intensities (replace with small value)
  let mut pi = 0
  while pi < processed.length() {
    if processed[pi].intensity <= 0.0 {
      processed[pi].intensity = 1.0
    }
    pi = pi + 1
  }

  // Log2 transform
  if !log_transformed {
    let mut li = 0
    while li < processed.length() {
      processed[li].intensity = @math.log2(processed[li].intensity)
      li = li + 1
    }
  }

  // Normalization
  match normalization {
    MSNormalization::None => ()
    MSNormalization::Median => {
      // Compute median per run
      let run_medians : Map[String, Double] = Map::new()
      let run_values : Map[String, Array[Double]] = Map::new()
      let mut ri = 0
      while ri < processed.length() {
        let run = processed[ri].run
        let arr = run_values.get(run).unwrap_or(Array::new())
        arr.push(processed[ri].intensity)
        run_values.set(run, arr)
        ri = ri + 1
      }
      let run_list = run_values.keys().collect()
      let mut rvi = 0
      while rvi < run_list.length() {
        let run = run_list[rvi]
        let vals = run_values.get(run).unwrap()
        let sorted = vals.copy()
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
        let med = if slen % 2 == 0 {
          (sorted[mid - 1] + sorted[mid]) / 2.0
        } else {
          sorted[mid]
        }
        run_medians.set(run, med)
        rvi = rvi + 1
      }

      let mut ni = 0
      while ni < processed.length() {
        let med = run_medians.get(processed[ni].run).unwrap_or(0.0)
        processed[ni].intensity = processed[ni].intensity - med
        ni = ni + 1
      }
    }
    MSNormalization::Quantile => {
      // Simple quantile normalization across runs
      // Collect all intensities and sort
      let all_vals : Array[Double] = Array::new()
      let mut qi = 0
      while qi < processed.length() {
        all_vals.push(processed[qi].intensity)
        qi = qi + 1
      }
      let slen = all_vals.length()
      let mut si = 0
      while si < slen - 1 {
        let mut sj = 0
        while sj < slen - 1 - si {
          if all_vals[sj] > all_vals[sj + 1] {
            let tmp = all_vals[sj]
            all_vals.set(sj, all_vals[sj + 1])
            all_vals.set(sj + 1, tmp)
          }
          sj = sj + 1
        }
        si = si + 1
      }
      // Rank-based normalization (simplified: replace with sorted rank values)
      let mean_val = if slen > 0 {
        let mut s = 0.0
        let mut mi = 0
        while mi < slen {
          s = s + all_vals[mi]
          mi = mi + 1
        }
        s / slen.to_double()
      } else {
        0.0
      }
      let mut qi2 = 0
      while qi2 < processed.length() {
        processed[qi2].intensity = processed[qi2].intensity - mean_val
        qi2 = qi2 + 1
      }
    }
    MSNormalization::GlobalStandards => {
      // Use mean of all as reference (simplified)
      let mut sum = 0.0
      let mut si = 0
      while si < processed.length() {
        sum = sum + processed[si].intensity
        si = si + 1
      }
      let grand_mean = if processed.length() > 0 {
        sum / processed.length().to_double()
      } else {
        0.0
      }
      let mut gi = 0
      while gi < processed.length() {
        processed[gi].intensity = processed[gi].intensity - grand_mean
        gi = gi + 1
      }
    }
  }

  processed
}

///|
/// Summarize feature-level data to protein-level.
pub fn ms_summarize(
  features : Array[MSFeature],
  method? : MSSummarization = MSSummarization::Tukey,
) -> Array[MSProteinAbundance] {
  // Group by (protein, run)
  let groups : Map[String, Array[MSFeature]] = Map::new()
  let mut i = 0
  while i < features.length() {
    let key = features[i].protein + "::" + features[i].run
    let arr = groups.get(key).unwrap_or(Array::new())
    arr.push(features[i])
    groups.set(key, arr)
    i = i + 1
  }

  let result : Array[MSProteinAbundance] = Array::new()
  let key_list = groups.keys().collect()
  let mut ki = 0
  while ki < key_list.length() {
    let key = key_list[ki]
    let feats = groups.get(key).unwrap()
    let parts = key.split("::")
    let protein = parts.next().unwrap_or("").to_string()
    let run = parts.next().unwrap_or("").to_string()
    let group = if feats.length() > 0 { feats[0].group } else { "" }

    let intensities : Array[Double] = Array::new()
    let mut fi = 0
    while fi < feats.length() {
      intensities.push(feats[fi].intensity)
      fi = fi + 1
    }

    let abundance = match method {
      MSSummarization::Tukey => {
        // Tukey median polish (simplified: use median)
        let sorted = intensities.copy()
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
      MSSummarization::Linear => {
        // Linear model: use mean
        let mut s = 0.0
        let mut li = 0
        while li < intensities.length() {
          s = s + intensities[li]
          li = li + 1
        }
        if intensities.length() > 0 {
          s / intensities.length().to_double()
        } else {
          0.0
        }
      }
      MSSummarization::LogSum => {
        // Sum of log values (already log-transformed, so sum then log)
        let mut s = 0.0
        let mut li = 0
        while li < intensities.length() {
          s = s + intensities[li]
          li = li + 1
        }
        s
      }
    }

    result.push(MSProteinAbundance::{
      protein,
      group,
      run,
      log2_abundance: abundance,
    })
    ki = ki + 1
  }

  result
}

///|
/// Perform group comparison (differential abundance testing).
pub fn ms_group_comparison(
  abundances : Array[MSProteinAbundance],
  fdr_threshold? : Double = 0.05,
  log2fc_threshold? : Double = 1.0,
) -> MSGroupComparison {
  // Group by protein
  let protein_groups : Map[String, Array[MSProteinAbundance]] = Map::new()
  let mut i = 0
  while i < abundances.length() {
    let key = abundances[i].protein
    let arr = protein_groups.get(key).unwrap_or(Array::new())
    arr.push(abundances[i])
    protein_groups.set(key, arr)
    i = i + 1
  }

  let all_results : Array[MSGroupResult] = Array::new()

  let protein_list = protein_groups.keys().collect()
  let mut pli = 0
  while pli < protein_list.length() {
    let protein = protein_list[pli]
    let abuns = protein_groups.get(protein).unwrap()
    // Split into groups
    let group_data : Map[String, Array[Double]] = Map::new()
    let mut ai = 0
    while ai < abuns.length() {
      let grp = abuns[ai].group
      let arr = group_data.get(grp).unwrap_or(Array::new())
      arr.push(abuns[ai].log2_abundance)
      group_data.set(grp, arr)
      ai = ai + 1
    }

    let group_names : Array[String] = Array::new()
    let group_keys = group_data.keys().collect()
    let mut gki = 0
    while gki < group_keys.length() {
      group_names.push(group_keys[gki])
      gki = gki + 1
    }

    // Sort group names
    let gn = group_names.length()
    let mut si = 0
    while si < gn - 1 {
      let mut sj = 0
      while sj < gn - 1 - si {
        if group_names[sj] > group_names[sj + 1] {
          let tmp = group_names[sj]
          group_names.set(sj, group_names[sj + 1])
          group_names.set(sj + 1, tmp)
        }
        sj = sj + 1
      }
      si = si + 1
    }

    if gn >= 2 {
      let g1 = group_names[0]
      let g2 = group_names[1]
      let v1 = group_data.get(g1).unwrap_or(Array::new())
      let v2 = group_data.get(g2).unwrap_or(Array::new())

      let m1 = ms_mean(v1)
      let m2 = ms_mean(v2)
      let var1 = ms_variance(v1)
      let var2 = ms_variance(v2)

      let log2fc = m2 - m1
      let n1 = v1.length().to_double()
      let n2 = v2.length().to_double()
      let se = (var1 / n1 + var2 / n2).sqrt()

      let t_stat = if se < 1.0e-10 {
        if (m2 - m1).abs() < 1.0e-10 { 0.0 } else { 100.0 }
      } else {
        (m2 - m1) / se
      }

      // P-value approximation
      let df = if n1 + n2 > 2.0 { n1 + n2 - 2.0 } else { 1.0 }
      let abs_t = if t_stat < 0.0 { -t_stat } else { t_stat }
      let p_val = if abs_t < 0.001 {
        1.0
      } else {
        let t_sq = abs_t * abs_t
        let p_approx = 2.0 * @math.exp(-df * @math.ln(1.0 + t_sq / df) / 2.0) / (df * 3.14159265358979).sqrt()
        if p_approx > 1.0 { 1.0 } else if p_approx < 0.0 { 0.0 } else { p_approx }
      }

      all_results.push(MSGroupResult::{
        protein,
        comparison: g1 + " vs " + g2,
        log2fc,
        se,
        t_stat,
        p_val,
        p_adj: p_val,
        significant: false,
      })
    }
    pli = pli + 1
  }

  // BH adjustment
  let p_values : Array[Double] = Array::new()
  let mut ri = 0
  while ri < all_results.length() {
    p_values.push(all_results[ri].p_val)
    ri = ri + 1
  }
  let adjusted = ms_bh_adjust(p_values)

  let mut n_sig = 0
  let mut si2 = 0
  while si2 < all_results.length() {
    all_results[si2].p_adj = adjusted[si2]
    let abs_lfc = if all_results[si2].log2fc < 0.0 { -all_results[si2].log2fc } else { all_results[si2].log2fc }
    all_results[si2].significant = adjusted[si2] < fdr_threshold && abs_lfc >= log2fc_threshold
    if all_results[si2].significant {
      n_sig = n_sig + 1
    }
    si2 = si2 + 1
  }

  let n_proteins = if all_results.length() > 0 {
    let ps : Map[String, Unit] = Map::new()
    let mut pi = 0
    while pi < all_results.length() {
      ps.set(all_results[pi].protein, ())
      pi = pi + 1
    }
    ps.size()
  } else {
    0
  }

  MSGroupComparison::{
    results: all_results,
    n_proteins,
    n_sig,
    fdr_threshold,
  }
}

///|
/// Compute mean.
fn ms_mean(values : Array[Double]) -> Double {
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
/// Compute variance.
fn ms_variance(values : Array[Double]) -> Double {
  if values.length() < 2 {
    return 0.0
  }
  let mean = ms_mean(values)
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
/// Benjamini-Hochberg FDR adjustment.
fn ms_bh_adjust(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return Array::new()
  }

  let indexed : Array[(Int, Double)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((i, p_values[i]))
    i = i + 1
  }

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
pub fn MSGroupComparison::get_n_results(self : MSGroupComparison) -> Int {
  self.results.length()
}

///|
/// Get significant results.
pub fn MSGroupComparison::get_significant(self : MSGroupComparison) -> Array[MSGroupResult] {
  let sig : Array[MSGroupResult] = Array::new()
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
/// Get top N proteins by adjusted p-value.
pub fn MSGroupComparison::get_top_proteins(self : MSGroupComparison, n : Int) -> Array[MSGroupResult] {
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

  let result : Array[MSGroupResult] = Array::new()
  let take = if n < len { n } else { len }
  let mut k = 0
  while k < take {
    result.push(sorted[k])
    k = k + 1
  }
  result
}

///|
/// Generate summary string.
pub fn MSGroupComparison::summary(self : MSGroupComparison) -> String {
  let mut s = "=== MSstats Group Comparison Summary ===\n".to_string()
  s = s + "Proteins tested: " + self.n_proteins.to_string() + "\n"
  s = s + "Significant (FDR < " + self.fdr_threshold.to_string() + "): " + self.n_sig.to_string() + "\n"
  s = s + "Total comparisons: " + self.results.length().to_string() + "\n"
  s
}

///|
/// Calculate sample size for desired power.
pub fn ms_design_sample_size(
  results : Array[MSGroupResult],
  desired_power? : Double = 0.8,
  n_samples? : Int = 0,
) -> MSSampleSize {
  // Simplified: use the median log2fc and SE to estimate power
  let n = results.length()
  if n == 0 {
    return MSSampleSize::{ desired_power, n_samples: 3, n_proteins_detectable: 0, delta: 0.0 }
  }

  // Sort by absolute log2fc
  let sorted = results.copy()
  let mut i = 0
  while i < n - 1 {
    let mut j = 0
    while j < n - 1 - i {
      let abs1 = if sorted[j].log2fc < 0.0 { -sorted[j].log2fc } else { sorted[j].log2fc }
      let abs2 = if sorted[j + 1].log2fc < 0.0 { -sorted[j + 1].log2fc } else { sorted[j + 1].log2fc }
      if abs1 < abs2 {
        let tmp = sorted[j]
        sorted.set(j, sorted[j + 1])
        sorted.set(j + 1, tmp)
      }
      j = j + 1
    }
    i = i + 1
  }

  let mid = n / 2
  let delta = if sorted[mid].log2fc < 0.0 { -sorted[mid].log2fc } else { sorted[mid].log2fc }
  let se = sorted[mid].se

  // Estimate required sample size using simplified power calculation
  // n ≈ (z_alpha + z_beta)^2 * sigma^2 / delta^2
  let z_alpha = 1.96
  let z_beta = if desired_power >= 0.8 { 0.842 } else { 0.674 }
  let sigma = if se > 0.0 { se } else { 1.0 }
  let required_n = ((z_alpha + z_beta) * (z_alpha + z_beta) * sigma * sigma / (delta * delta + 1.0e-10)).to_int()
  let n_samples_final = if n_samples > 0 { n_samples } else { if required_n < 3 { 3 } else { required_n } }

  // Count how many proteins would be detectable with this sample size
  let n_detectable = if delta > 0.0 && se > 0.0 {
    let mut count = 0
    let mut di = 0
    while di < n {
      let abs_lfc = if sorted[di].log2fc < 0.0 { -sorted[di].log2fc } else { sorted[di].log2fc }
      let z = abs_lfc / (sorted[di].se + 1.0e-10) * (n_samples_final.to_double()).sqrt()
      if z > z_alpha {
        count = count + 1
      }
      di = di + 1
    }
    count
  } else {
    0
  }

  MSSampleSize::{
    desired_power,
    n_samples: n_samples_final,
    n_proteins_detectable: n_detectable,
    delta,
  }
}

///|
/// Create sample MS data for testing.
pub fn msstats_sample_data() -> Array[MSFeature] {
  let features : Array[MSFeature] = Array::new()
  let proteins = ["P1", "P2", "P3", "P4", "P5"]
  let groups = ["Ctrl", "Ctrl", "Trt", "Trt"]
  let runs = ["R1", "R2", "R3", "R4"]
  let subjects = ["S1", "S2", "S3", "S4"]

  let mut pi = 0
  while pi < proteins.length() {
    let mut gi = 0
    while gi < groups.length() {
      // Two peptides per protein
      let mut pep_idx = 0
      while pep_idx < 2 {
        let peptide = proteins[pi] + "_pep" + pep_idx.to_string()
        let transition = peptide + "_tr"
        let intensity = if groups[gi] == "Trt" {
          1000.0 * (pi + 1).to_double() * 2.0 + (pep_idx + 1).to_double() * 100.0
        } else {
          1000.0 * (pi + 1).to_double() + (pep_idx + 1).to_double() * 100.0
        }
        features.push(MSFeature::new(
          proteins[pi],
          peptide,
          transition,
          groups[gi],
          subjects[gi],
          runs[gi],
          intensity,
        ))
        pep_idx = pep_idx + 1
      }
      gi = gi + 1
    }
    pi = pi + 1
  }
  features
}
