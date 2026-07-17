///|
/// Affymetrix microarray data analysis module, inspired by Biopython Bio.Affy.

///|
pub struct ProbeSet {
  probe_set_id : String
  pm_intensities : Array[Array[Double]]
  mm_intensities : Array[Array[Double]]
}

///|
pub fn ProbeSet::new(
  probe_set_id : String,
  pm_intensities : Array[Array[Double]],
  mm_intensities : Array[Array[Double]],
) -> ProbeSet {
  ProbeSet::{ probe_set_id, pm_intensities, mm_intensities }
}

///|
pub struct AffyBatch {
  probe_set_ids : Array[String]
  probe_sets : Array[ProbeSet]
  sample_names : Array[String]
  n_samples : Int
  n_probe_sets : Int
}

///|
pub fn AffyBatch::new(
  probe_set_ids : Array[String],
  probe_sets : Array[ProbeSet],
  sample_names : Array[String],
) -> AffyBatch {
  let n_samples = if probe_sets.length() > 0 &&
    probe_sets[0].pm_intensities.length() > 0 {
    probe_sets[0].pm_intensities[0].length()
  } else {
    0
  }
  AffyBatch::{
    probe_set_ids,
    probe_sets,
    sample_names,
    n_samples,
    n_probe_sets: probe_set_ids.length(),
  }
}

///|
pub fn AffyBatch::n_probes_per_set(self : AffyBatch) -> Int {
  if self.probe_sets.length() > 0 {
    self.probe_sets[0].pm_intensities.length()
  } else {
    0
  }
}

///|
fn log2(x : Double) -> Double {
  @math.log2(x)
}

///|
fn sort_doubles(arr : Array[Double]) -> Array[Double] {
  let sorted = arr.copy()
  let n = sorted.length()
  let mut i = 0
  while i < n - 1 {
    let mut j = i + 1
    while j < n {
      if sorted[i] > sorted[j] {
        let tmp = sorted[i]
        sorted[i] = sorted[j]
        sorted[j] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  sorted
}

///|
fn median(arr : Array[Double]) -> Double {
  let sorted = sort_doubles(arr)
  let n = sorted.length()
  if n % 2 == 0 {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  } else {
    sorted[n / 2]
  }
}

///|
fn quantile_normalize(data : Array[Array[Double]]) -> Array[Array[Double]] {
  let n_rows = data.length()
  if n_rows == 0 {
    return data
  }
  let n_cols = data[0].length()

  let sorted_data : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_rows {
    sorted_data.push(sort_doubles(data[i]))
    i = i + 1
  }

  let quantiles : Array[Double] = Array::new()
  let mut j = 0
  while j < n_cols {
    let col_vals : Array[Double] = Array::new()
    let mut k = 0
    while k < n_rows {
      col_vals.push(sorted_data[k][j])
      k = k + 1
    }
    quantiles.push(median(col_vals))
    j = j + 1
  }

  let result : Array[Array[Double]] = Array::new()
  i = 0
  while i < n_rows {
    let row = data[i]
    let sorted_row = sort_doubles(row)
    let ranks : Array[Int] = Array::new()
    let mut m = 0
    while m < row.length() {
      let val = row[m]
      let mut rank = 0
      let mut n = 0
      while n < sorted_row.length() {
        if sorted_row[n] == val {
          rank = n
          break
        }
        n = n + 1
      }
      ranks.push(rank)
      m = m + 1
    }

    let normalized_row : Array[Double] = Array::new()
    m = 0
    while m < ranks.length() {
      normalized_row.push(quantiles[ranks[m]])
      m = m + 1
    }
    result.push(normalized_row)
    i = i + 1
  }

  result
}

///|
pub fn background_correct_rma(
  pm_intensities : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n_probes = pm_intensities.length()
  if n_probes == 0 {
    return pm_intensities
  }
  let n_samples = pm_intensities[0].length()

  let corrected : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_probes {
    let probe_vals : Array[Double] = Array::new()
    let mut j = 0
    while j < n_samples {
      let val = pm_intensities[i][j]
      let bg = 50.0
      let corrected_val = if val > bg { val - bg } else { 1.0 }
      probe_vals.push(corrected_val)
      j = j + 1
    }
    corrected.push(probe_vals)
    i = i + 1
  }

  corrected
}

///|
pub fn rma_normalize(affy_batch : AffyBatch) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()

  let mut i = 0
  while i < affy_batch.probe_sets.length() {
    let probe_set = affy_batch.probe_sets[i]
    let pm = probe_set.pm_intensities

    let bg_corrected = background_correct_rma(pm)

    let log_transformed : Array[Array[Double]] = Array::new()
    let mut j = 0
    while j < bg_corrected.length() {
      let log_row : Array[Double] = Array::new()
      let mut k = 0
      while k < bg_corrected[j].length() {
        log_row.push(log2(bg_corrected[j][k]))
        k = k + 1
      }
      log_transformed.push(log_row)
      j = j + 1
    }

    let normalized = quantile_normalize(log_transformed)

    let summarized : Array[Double] = Array::new()
    j = 0
    while j < affy_batch.n_samples {
      let col_vals : Array[Double] = Array::new()
      let mut k = 0
      while k < normalized.length() {
        col_vals.push(normalized[k][j])
        k = k + 1
      }
      summarized.push(median(col_vals))
      j = j + 1
    }

    result.push(summarized)
    i = i + 1
  }

  result
}

///|
pub fn affy_batch_from_matrix(
  probe_set_ids : Array[String],
  pm_matrix : Array[Array[Array[Double]]],
  sample_names : Array[String],
) -> AffyBatch {
  let probe_sets : Array[ProbeSet] = Array::new()
  let mut i = 0
  while i < probe_set_ids.length() {
    let mm : Array[Array[Double]] = Array::new()
    probe_sets.push(ProbeSet::new(probe_set_ids[i], pm_matrix[i], mm))
    i = i + 1
  }
  AffyBatch::new(probe_set_ids, probe_sets, sample_names)
}

///|
pub fn compute_pm_mm_difference(probe_set : ProbeSet) -> Array[Array[Double]] {
  let diff : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < probe_set.pm_intensities.length() {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < probe_set.pm_intensities[i].length() {
      let mm_val = if i < probe_set.mm_intensities.length() &&
        j < probe_set.mm_intensities[i].length() {
        probe_set.mm_intensities[i][j]
      } else {
        0.0
      }
      row.push(probe_set.pm_intensities[i][j] - mm_val)
      j = j + 1
    }
    diff.push(row)
    i = i + 1
  }
  diff
}

///|
pub fn summarize_probeset(
  probe_set : ProbeSet,
  method : String,
) -> Array[Double] {
  let pm = probe_set.pm_intensities
  let n_samples = if pm.length() > 0 { pm[0].length() } else { 0 }

  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < n_samples {
    let vals : Array[Double] = Array::new()
    let mut j = 0
    while j < pm.length() {
      vals.push(pm[j][i])
      j = j + 1
    }

    let val = if method == "median" {
      median(vals)
    } else if method == "mean" {
      let mut sum = 0.0
      let mut k = 0
      while k < vals.length() {
        sum = sum + vals[k]
        k = k + 1
      }
      sum / vals.length().to_double()
    } else {
      median(vals)
    }

    result.push(val)
    i = i + 1
  }

  result
}

///|
pub fn affy_batch_summarize(
  affy_batch : AffyBatch,
  method : String,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < affy_batch.probe_sets.length() {
    result.push(summarize_probeset(affy_batch.probe_sets[i], method))
    i = i + 1
  }
  result
}
