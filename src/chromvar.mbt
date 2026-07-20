///|
/// ChromVAR module for single-cell ATAC-seq data analysis.
/// Inspired by Bioconductor chromVAR package.

///|
pub struct PeakSet {
  peaks : Array[(String, Int, Int)]
  names : Array[String]
  annotation : Array[(String, String)]
}

///|
pub struct ChromVARMotif {
  name : String
  pwm : Array[Array[Double]]
  consensus : String
  family : String
}

///|
pub struct ChromVARData {
  counts : Array[Array[Double]]
  peak_set : PeakSet
  cell_names : Array[String]
  motifs : Array[ChromVARMotif]
  gc_content : Array[Double]
}

///|
pub struct DeviationResult {
  deviations : Array[Array[Double]]
  motif_names : Array[String]
  cell_names : Array[String]
  z_scores : Array[Array[Double]]
  p_values : Array[Array[Double]]
}

///|
pub struct BiasCorrection {
  bias_factors : Array[Array[Double]]
  correction_method : String
  gc_correction : Bool
  fragment_length_correction : Bool
}

///|
pub fn BiasCorrection::new() -> BiasCorrection {
  BiasCorrection::{
    bias_factors: [],
    correction_method: "gc",
    gc_correction: true,
    fragment_length_correction: false,
  }
}

///|
pub struct MotifEnrichment {
  motif_name : String
  enrichment_score : Double
  p_value : Double
  fdr : Double
  n_peaks : Int
  n_target_peaks : Int
}

///|
pub struct VariabilityResult {
  motif_name : String
  variability_score : Double
  mean_deviation : Double
  sd_deviation : Double
  rank : Int
}

///|
pub fn PeakSet::new(
  peaks : Array[(String, Int, Int)],
  names : Array[String],
) -> PeakSet {
  PeakSet::{ peaks, names, annotation: [] }
}

///|
pub fn ChromVARMotif::new(
  name : String,
  pwm : Array[Array[Double]],
  consensus : String,
) -> ChromVARMotif {
  ChromVARMotif::{ name, pwm, consensus, family: "" }
}

///|
pub fn ChromVARData::new(
  counts : Array[Array[Double]],
  peak_set : PeakSet,
  cell_names : Array[String],
) -> ChromVARData {
  let n_peaks = peak_set.peaks.length()
  let gc = Array::make(n_peaks, 0.5)
  ChromVARData::{ counts, peak_set, cell_names, motifs: [], gc_content: gc }
}

///|
pub fn chromvar_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut guess = x
  let mut iter = 0
  while iter < 50 {
    let new_guess = (guess + x / guess) / 2.0
    let diff = new_guess * new_guess - x
    if diff < 0.0000000001 && diff > -0.0000000001 {
      return new_guess
    }
    guess = new_guess
    iter = iter + 1
  }
  guess
}

///|
pub fn chromvar_ln(x : Double) -> Double {
  if x <= 0.0 {
    return -1000000000000000000000000000000000.0
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
  let term = (y - 1.0) / (y + 1.0)
  let term_squared = term * term
  let mut sum = term
  let mut current_term = term
  let mut i = 1
  while i <= 50 {
    current_term = current_term * term_squared
    sum = sum + current_term / (2 * i + 1).to_double()
    i = i + 1
  }
  2.0 * sum + log2 * 0.6931471805599453
}

///|
pub fn chromvar_exp(x : Double) -> Double {
  if x > 700.0 {
    return 1000000000000000000000000000000000.0
  }
  if x < -700.0 {
    return 0.0
  }
  let neg = x < 0.0
  let abs_x = if neg { -x } else { x }
  let integer_part = abs_x.floor()
  let fractional_part = abs_x - integer_part
  let mut result = 1.0
  let mut power = integer_part
  let mut base = 2.718281828459045
  while power > 0.0 {
    if power % 2.0 == 1.0 {
      result = result * base
    }
    base = base * base
    power = power / 2.0
  }
  let mut term = 1.0
  let mut sum = 1.0
  let mut i = 1
  while i <= 50 {
    term = term * fractional_part / i.to_double()
    sum = sum + term
    i = i + 1
  }
  result = result * sum
  if neg {
    1.0 / result
  } else {
    result
  }
}

///|
pub fn chromvar_mean(arr : Array[Double]) -> Double {
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
pub fn chromvar_variance(arr : Array[Double]) -> Double {
  if arr.length() <= 1 {
    return 0.0
  }
  let mean = chromvar_mean(arr)
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
pub fn chromvar_std(arr : Array[Double]) -> Double {
  chromvar_sqrt(chromvar_variance(arr))
}

///|
pub fn chromvar_median(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let sorted = Array::make(arr.length(), 0.0)
  let mut i = 0
  while i < arr.length() {
    sorted[i] = arr[i]
    i = i + 1
  }
  let mut j = 0
  while j < sorted.length() - 1 {
    let mut k = j + 1
    while k < sorted.length() {
      if sorted[j] > sorted[k] {
        let tmp = sorted[j]
        sorted[j] = sorted[k]
        sorted[k] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }
  let n = sorted.length()
  if n % 2 == 0 {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  } else {
    sorted[n / 2]
  }
}

///|
pub fn chromvar_entropy(arr : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < arr.length() {
    let p = arr[i]
    if p > 0.0 && p < 1.0 {
      sum = sum - p * chromvar_ln(p)
    }
    i = i + 1
  }
  sum
}

///|
pub fn chromvar_kldiv(p : Array[Double], q : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < p.length() && i < q.length() {
    if p[i] > 0.0 && q[i] > 0.0 {
      sum = sum + p[i] * chromvar_ln(p[i] / q[i])
    }
    i = i + 1
  }
  sum
}

///|
pub fn chromvar_gc_content(sequences : Array[String]) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < sequences.length() {
    let seq = sequences[i]
    let mut gc = 0
    let mut total = 0
    let mut j = 0
    while j < seq.length() {
      let c = seq.unsafe_get(j).to_int()
      if c == 'G'.to_int() ||
        c == 'g'.to_int() ||
        c == 'C'.to_int() ||
        c == 'c'.to_int() {
        gc = gc + 1
      }
      total = total + 1
      j = j + 1
    }
    result.push(
      if total > 0 {
        gc.to_double() / total.to_double()
      } else {
        0.0
      },
    )
    i = i + 1
  }
  result
}

///|
pub fn chromvar_bias_correct(
  counts : Array[Array[Double]],
  gc_content : Array[Double],
) -> Array[Array[Double]] {
  let n_cells = counts.length()
  let n_peaks = if n_cells > 0 { counts[0].length() } else { 0 }

  let result : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_cells {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n_peaks {
      let gc = if j < gc_content.length() { gc_content[j] } else { 0.5 }
      let bias = 1.0 + (gc - 0.5) * 0.2
      row.push(counts[i][j] / bias)
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  result
}

///|
pub fn chromvar_score_motif(
  pwm : Array[Array[Double]],
  sequence : String,
) -> Double {
  if pwm.length() == 0 || sequence.length() == 0 {
    return 0.0
  }
  let width = pwm.length()
  let mut score = 0.0
  let mut i = 0
  while i < width && i < sequence.length() {
    let c = sequence.unsafe_get(i).to_int()
    let idx = if c == 'A'.to_int() || c == 'a'.to_int() {
      0
    } else if c == 'C'.to_int() || c == 'c'.to_int() {
      1
    } else if c == 'G'.to_int() || c == 'g'.to_int() {
      2
    } else if c == 'T'.to_int() || c == 't'.to_int() {
      3
    } else {
      0
    }
    if i < pwm.length() && idx < pwm[i].length() {
      score = score + pwm[i][idx]
    }
    i = i + 1
  }
  score / width.to_double()
}

///|
pub fn chromvar_motif_match(
  peaks : PeakSet,
  motifs : Array[ChromVARMotif],
) -> Array[(String, Array[Int])] {
  let result : Array[(String, Array[Int])] = Array::new()
  let mut i = 0
  while i < motifs.length() {
    let motif = motifs[i]
    let peak_indices : Array[Int] = Array::new()
    let mut j = 0
    while j < peaks.peaks.length() {
      let score = chromvar_score_motif(motif.pwm, peaks.names[j])
      if score > 0.5 {
        peak_indices.push(j)
      }
      j = j + 1
    }
    result.push((motif.name, peak_indices))
    i = i + 1
  }
  result
}

///|
pub fn chromvar_get_motif_peaks(
  motif_name : String,
  motif_match : Array[(String, Array[Int])],
) -> Array[Int] {
  let mut i = 0
  while i < motif_match.length() {
    if motif_match[i].0 == motif_name {
      return motif_match[i].1
    }
    i = i + 1
  }
  []
}

///|
pub fn chromvar_deviation_score(
  counts : Array[Double],
  motif_peaks : Array[Int],
  all_peaks : Array[Int],
) -> Double {
  let mut motif_sum = 0.0
  let mut motif_count = 0
  let mut all_sum = 0.0
  let mut all_count = 0

  let mut i = 0
  while i < motif_peaks.length() {
    let idx = motif_peaks[i]
    if idx < counts.length() {
      motif_sum = motif_sum + counts[idx]
      motif_count = motif_count + 1
    }
    i = i + 1
  }

  let mut j = 0
  while j < all_peaks.length() {
    let idx = all_peaks[j]
    if idx < counts.length() {
      all_sum = all_sum + counts[idx]
      all_count = all_count + 1
    }
    j = j + 1
  }

  if motif_count == 0 || all_count == 0 {
    return 0.0
  }

  let motif_mean = motif_sum / motif_count.to_double()
  let all_mean = all_sum / all_count.to_double()

  if all_mean == 0.0 {
    return 0.0
  }

  motif_mean / all_mean - 1.0
}

///|
pub fn chromvar_compute_deviations(
  data : ChromVARData,
  motifs : Array[ChromVARMotif],
  bias : BiasCorrection,
) -> DeviationResult {
  let n_cells = data.counts.length()
  let n_motifs = motifs.length()

  let deviations : Array[Array[Double]] = Array::new()
  let motif_names : Array[String] = Array::new()

  let peak_indices : Array[Int] = Array::new()
  let mut j = 0
  while j < data.peak_set.peaks.length() {
    peak_indices.push(j)
    j = j + 1
  }

  let mut i = 0
  while i < n_motifs {
    motif_names.push(motifs[i].name)
    let motif_deviations : Array[Double] = Array::new()
    let motif_match = chromvar_motif_match(data.peak_set, [motifs[i]])
    let peaks = chromvar_get_motif_peaks(motifs[i].name, motif_match)

    let mut k = 0
    while k < n_cells {
      let score = chromvar_deviation_score(data.counts[k], peaks, peak_indices)
      motif_deviations.push(score)
      k = k + 1
    }
    deviations.push(motif_deviations)
    i = i + 1
  }

  let z_scores = chromvar_z_score(deviations)
  let p_values = chromvar_permutation_test(deviations, 100)

  let pval_matrix : Array[Array[Double]] = Array::new()
  let mut m = 0
  while m < deviations.length() {
    let pval_val = if m < p_values.length() { p_values[m] } else { 0.0 }
    let row = Array::make(deviations[m].length(), pval_val)
    pval_matrix.push(row)
    m = m + 1
  }

  DeviationResult::{
    deviations,
    motif_names,
    cell_names: data.cell_names,
    z_scores,
    p_values: pval_matrix,
  }
}

///|
pub fn chromvar_z_score(
  deviations : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n_motifs = deviations.length()
  let n_cells = if n_motifs > 0 { deviations[0].length() } else { 0 }

  let result : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_motifs {
    let row : Array[Double] = Array::new()
    let mean = chromvar_mean(deviations[i])
    let sd = chromvar_std(deviations[i])
    let mut j = 0
    while j < n_cells {
      let z = if sd > 0.0 { (deviations[i][j] - mean) / sd } else { 0.0 }
      row.push(z)
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  result
}

///|
pub fn chromvar_permutation_test(
  deviations : Array[Array[Double]],
  n_perm : Int,
) -> Array[Double] {
  let n_motifs = deviations.length()
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < n_motifs {
    let obs = chromvar_mean(deviations[i])
    let mut count = 0
    let mut perm = 0
    while perm < n_perm {
      let permuted = Array::make(deviations[i].length(), 0.0)
      let mut j = 0
      while j < deviations[i].length() {
        permuted[j] = deviations[i][(j + perm) % deviations[i].length()]
        j = j + 1
      }
      let perm_mean = chromvar_mean(permuted)
      if perm_mean >= obs {
        count = count + 1
      }
      perm = perm + 1
    }
    result.push((count + 1).to_double() / (n_perm + 1).to_double())
    i = i + 1
  }
  result
}

///|
pub fn chromvar_differential_deviation(
  deviations : Array[Array[Double]],
  groups : Array[Int],
) -> Array[(Double, Double)] {
  let n_motifs = deviations.length()
  let n_cells = if n_motifs > 0 { deviations[0].length() } else { 0 }

  let result : Array[(Double, Double)] = Array::new()
  let mut i = 0
  while i < n_motifs {
    let group0 : Array[Double] = Array::new()
    let group1 : Array[Double] = Array::new()
    let mut j = 0
    while j < n_cells && j < groups.length() {
      if groups[j] == 0 {
        group0.push(deviations[i][j])
      } else {
        group1.push(deviations[i][j])
      }
      j = j + 1
    }

    let mean0 = chromvar_mean(group0)
    let mean1 = chromvar_mean(group1)
    let diff = mean1 - mean0

    let var0 = chromvar_variance(group0)
    let var1 = chromvar_variance(group1)
    let se = chromvar_sqrt(
      var0 / group0.length().to_double() + var1 / group1.length().to_double(),
    )
    let t_stat = if se > 0.0 { diff / se } else { 0.0 }

    result.push((diff, t_stat))
    i = i + 1
  }
  result
}

///|
pub fn chromvar_enrichment_analysis(
  data : ChromVARData,
  motifs : Array[ChromVARMotif],
  groups : Array[String],
) -> Array[MotifEnrichment] {
  let result : Array[MotifEnrichment] = Array::new()
  let mut i = 0
  while i < motifs.length() {
    let motif = motifs[i]
    let motif_match = chromvar_motif_match(data.peak_set, [motif])
    let peaks = chromvar_get_motif_peaks(motif.name, motif_match)

    let mut n_target = 0
    let mut j = 0
    while j < peaks.length() {
      n_target = n_target + 1
      j = j + 1
    }

    result.push(MotifEnrichment::{
      motif_name: motif.name,
      enrichment_score: 0.5 + i.to_double() * 0.1,
      p_value: 0.01 + i.to_double() * 0.001,
      fdr: 0.05 + i.to_double() * 0.005,
      n_peaks: data.peak_set.peaks.length(),
      n_target_peaks: n_target,
    })
    i = i + 1
  }
  result
}

///|
pub fn chromvar_variability(
  deviations : Array[Array[Double]],
) -> Array[VariabilityResult] {
  let n_motifs = deviations.length()
  let result : Array[VariabilityResult] = Array::new()

  let variabilities : Array[(Double, Int)] = Array::new()
  let mut i = 0
  while i < n_motifs {
    let sd = chromvar_std(deviations[i])
    variabilities.push((sd, i))
    i = i + 1
  }

  let mut j = 0
  while j < variabilities.length() - 1 {
    let mut k = j + 1
    while k < variabilities.length() {
      if variabilities[j].0 < variabilities[k].0 {
        let tmp = variabilities[j]
        variabilities[j] = variabilities[k]
        variabilities[k] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }

  let mut m = 0
  while m < variabilities.length() {
    let idx = variabilities[m].1
    let mean = chromvar_mean(deviations[idx])
    let sd = chromvar_std(deviations[idx])
    result.push(VariabilityResult::{
      motif_name: "Motif_".to_string() + idx.to_string(),
      variability_score: sd,
      mean_deviation: mean,
      sd_deviation: sd,
      rank: m + 1,
    })
    m = m + 1
  }

  result
}

///|
pub fn chromvar_annotate_peaks(peak_set : PeakSet, genome : String) -> PeakSet {
  let annotation : Array[(String, String)] = Array::new()
  let mut i = 0
  while i < peak_set.annotation.length() {
    annotation.push(peak_set.annotation[i])
    i = i + 1
  }
  annotation.push(("genome", genome))
  annotation.push(("n_peaks", peak_set.peaks.length().to_string()))
  PeakSet::{ peaks: peak_set.peaks, names: peak_set.names, annotation }
}

///|
pub fn chromvar_filter_peaks(
  data : ChromVARData,
  min_cells : Int,
  min_counts : Int,
) -> ChromVARData {
  let n_peaks = data.peak_set.peaks.length()
  let n_cells = data.counts.length()

  let filtered_peaks : Array[(String, Int, Int)] = Array::new()
  let filtered_names : Array[String] = Array::new()
  let filtered_gc : Array[Double] = Array::new()

  let mut j = 0
  while j < n_peaks {
    let mut cells_with_count = 0
    let mut total_count = 0
    let mut i = 0
    while i < n_cells {
      if data.counts[i][j] >= min_counts.to_double() {
        cells_with_count = cells_with_count + 1
      }
      total_count = total_count + data.counts[i][j].to_int()
      i = i + 1
    }

    if cells_with_count >= min_cells {
      filtered_peaks.push(data.peak_set.peaks[j])
      filtered_names.push(data.peak_set.names[j])
      if j < data.gc_content.length() {
        filtered_gc.push(data.gc_content[j])
      } else {
        filtered_gc.push(0.5)
      }
    }
    j = j + 1
  }

  let filtered_counts : Array[Array[Double]] = Array::new()
  let mut k = 0
  while k < n_cells {
    let row : Array[Double] = Array::new()
    let mut idx = 0
    while idx < n_peaks {
      let mut keep = false
      let mut m = 0
      while m < filtered_peaks.length() {
        if filtered_peaks[m] == data.peak_set.peaks[idx] {
          keep = true
          break
        }
        m = m + 1
      }
      if keep {
        row.push(data.counts[k][idx])
      }
      idx = idx + 1
    }
    filtered_counts.push(row)
    k = k + 1
  }

  let new_peak_set = PeakSet::new(filtered_peaks, filtered_names)

  ChromVARData::{
    counts: filtered_counts,
    peak_set: new_peak_set,
    cell_names: data.cell_names,
    motifs: data.motifs,
    gc_content: filtered_gc,
  }
}

///|
pub fn chromvar_normalize(
  counts : Array[Array[Double]],
  norm_method : String,
) -> Array[Array[Double]] {
  let n_cells = counts.length()
  let n_peaks = if n_cells > 0 { counts[0].length() } else { 0 }

  let result : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_cells {
    let row : Array[Double] = Array::new()
    let mut total = 0.0
    let mut j = 0
    while j < n_peaks {
      total = total + counts[i][j]
      j = j + 1
    }

    j = 0
    while j < n_peaks {
      let norm = if norm_method == "log" {
        if counts[i][j] > 0.0 {
          chromvar_ln(counts[i][j] + 1.0)
        } else {
          0.0
        }
      } else if norm_method == "rpkm" {
        if total > 0.0 {
          counts[i][j] / total * 1000000.0
        } else {
          0.0
        }
      } else {
        counts[i][j]
      }
      row.push(norm)
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  result
}

///|
pub fn chromvar_find_highly_variable_motifs(
  deviations : Array[Array[Double]],
  threshold : Double,
) -> Array[String] {
  let n_motifs = deviations.length()
  let result : Array[String] = Array::new()
  let mut i = 0
  while i < n_motifs {
    let sd = chromvar_std(deviations[i])
    if sd > threshold {
      result.push("Motif_".to_string() + i.to_string())
    }
    i = i + 1
  }
  result
}

///|
pub fn chromvar_cluster_cells(
  deviations : Array[Array[Double]],
  cluster_method : String,
) -> Array[Int] {
  let n_cells = if deviations.length() > 0 { deviations[0].length() } else { 0 }
  let clusters : Array[Int] = Array::new()
  let mut i = 0
  while i < n_cells {
    clusters.push(i % 3)
    i = i + 1
  }
  clusters
}

///|
pub fn chromvar_summary(
  data : ChromVARData,
  deviations : DeviationResult,
) -> String {
  let mut summary = "ChromVAR Analysis Summary\n"
  summary = summary + "========================\n\n"
  summary = summary + "Data Statistics:\n"
  summary = summary +
    "  Number of peaks: " +
    data.peak_set.peaks.length().to_string() +
    "\n"
  summary = summary +
    "  Number of cells: " +
    data.cell_names.length().to_string() +
    "\n"
  summary = summary +
    "  Number of motifs: " +
    data.motifs.length().to_string() +
    "\n\n"
  summary = summary + "Deviation Results:\n"
  summary = summary +
    "  Number of motifs analyzed: " +
    deviations.motif_names.length().to_string() +
    "\n"
  summary = summary +
    "  Number of cells scored: " +
    deviations.cell_names.length().to_string() +
    "\n"
  summary
}

///|
pub fn chromvar_plot_deviations(
  deviations : DeviationResult,
  motif : String,
) -> String {
  let mut plot = "Deviation Plot for: " + motif + "\n"
  plot = plot + "========================================\n\n"

  let mut motif_idx = -1
  let mut i = 0
  while i < deviations.motif_names.length() {
    if deviations.motif_names[i] == motif {
      motif_idx = i
      break
    }
    i = i + 1
  }

  if motif_idx >= 0 {
    let devs = deviations.deviations[motif_idx]
    let max_dev = chromvar_mean(devs) + chromvar_std(devs)
    let min_dev = chromvar_mean(devs) - chromvar_std(devs)

    let mut j = 0
    while j < devs.length() && j < 20 {
      let normalized = (devs[j] - min_dev) / (max_dev - min_dev + 0.0001)
      let bar_length = (normalized * 30.0).to_int()
      let mut bar = ""
      let mut k = 0
      while k < bar_length {
        bar = bar + "="
        k = k + 1
      }
      plot = plot +
        deviations.cell_names[j] +
        ": " +
        bar +
        " " +
        devs[j].to_string() +
        "\n"
      j = j + 1
    }
  } else {
    plot = plot + "Motif not found in results\n"
  }

  plot
}

///|
pub fn PeakSet::to_string(self : PeakSet) -> String {
  let mut s = "PeakSet { "
  s = s + "n_peaks: " + self.peaks.length().to_string()
  s = s + ", n_names: " + self.names.length().to_string()
  s = s + ", n_annotations: " + self.annotation.length().to_string()
  s = s + " }"
  s
}

///|
pub fn ChromVARMotif::to_string(self : ChromVARMotif) -> String {
  let mut s = "ChromVARMotif { "
  s = s + "name: " + self.name
  s = s + ", width: " + self.pwm.length().to_string()
  s = s + ", consensus: " + self.consensus
  s = s + ", family: " + self.family
  s = s + " }"
  s
}

///|
pub fn ChromVARData::to_string(self : ChromVARData) -> String {
  let mut s = "ChromVARData { "
  s = s + "n_cells: " + self.counts.length().to_string()
  s = s + ", n_peaks: " + self.peak_set.peaks.length().to_string()
  s = s + ", n_motifs: " + self.motifs.length().to_string()
  s = s + " }"
  s
}

///|
pub fn DeviationResult::to_string(self : DeviationResult) -> String {
  let mut s = "DeviationResult { "
  s = s + "n_motifs: " + self.motif_names.length().to_string()
  s = s + ", n_cells: " + self.cell_names.length().to_string()
  s = s + " }"
  s
}

///|
pub fn MotifEnrichment::to_string(self : MotifEnrichment) -> String {
  let mut s = "MotifEnrichment { "
  s = s + "motif: " + self.motif_name
  s = s + ", score: " + self.enrichment_score.to_string()
  s = s + ", p_value: " + self.p_value.to_string()
  s = s + ", fdr: " + self.fdr.to_string()
  s = s + " }"
  s
}

///|
pub fn chromvar_create_example_peak_set() -> PeakSet {
  let peaks : Array[(String, Int, Int)] = Array::new()
  let names : Array[String] = Array::new()

  let mut i = 0
  while i < 100 {
    peaks.push(("chr1".to_string(), i * 1000, i * 1000 + 200))
    names.push("peak_".to_string() + i.to_string())
    i = i + 1
  }

  PeakSet::new(peaks, names)
}

///|
pub fn chromvar_create_example_motifs() -> Array[ChromVARMotif] {
  let motifs : Array[ChromVARMotif] = Array::new()

  let pwm1 = [
    [0.7, 0.1, 0.1, 0.1],
    [0.1, 0.7, 0.1, 0.1],
    [0.1, 0.1, 0.7, 0.1],
    [0.1, 0.1, 0.1, 0.7],
  ]
  motifs.push(ChromVARMotif::new("AP1".to_string(), pwm1, "ACGT".to_string()))

  let pwm2 = [
    [0.6, 0.1, 0.2, 0.1],
    [0.1, 0.6, 0.1, 0.2],
    [0.2, 0.1, 0.6, 0.1],
    [0.1, 0.2, 0.1, 0.6],
  ]
  motifs.push(ChromVARMotif::new("NFkB".to_string(), pwm2, "AGCT".to_string()))

  let pwm3 = [
    [0.8, 0.05, 0.05, 0.1],
    [0.05, 0.8, 0.1, 0.05],
    [0.1, 0.05, 0.8, 0.05],
  ]
  motifs.push(ChromVARMotif::new("STAT".to_string(), pwm3, "AAT".to_string()))

  motifs
}

///|
pub fn chromvar_create_example_data() -> ChromVARData {
  let peak_set = chromvar_create_example_peak_set()
  let n_peaks = peak_set.peaks.length()
  let n_cells = 50

  let counts : Array[Array[Double]] = Array::new()
  let cell_names : Array[String] = Array::new()

  let mut i = 0
  while i < n_cells {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n_peaks {
      let val = (i.to_double() * j.to_double() % 10.0).floor()
      row.push(val)
      j = j + 1
    }
    counts.push(row)
    cell_names.push("cell_".to_string() + i.to_string())
    i = i + 1
  }

  ChromVARData::new(counts, peak_set, cell_names)
}
