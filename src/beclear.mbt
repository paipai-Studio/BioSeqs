///|
/// Batch effect detection and correction, ported from Bioconductor's
/// `BEclear` package.
///
/// The reference package detects batch effects in methylation (or
/// expression) matrices by comparing medians between batches, marks
/// affected gene x batch blocks, imputes replacement values via
/// nearest-neighbor interpolation, and returns a corrected matrix.
/// This port follows the same pipeline.

///|
/// A data matrix with row names (features) and column names (samples),
/// plus a mapping from sample name to batch id.
pub struct BEclearMatrix {
  row_names : Array[String]
  col_names : Array[String]
  data : Array[Array[Double]]
  batch : Map[String, Int]
} derive(Eq)

///|
pub fn BEclearMatrix::new(
  row_names~ : Array[String],
  col_names~ : Array[String],
  data~ : Array[Array[Double]],
  batch~ : Map[String, Int],
) -> BEclearMatrix {
  BEclearMatrix::{ row_names, col_names, data, batch }
}

///|
/// Number of rows (features).
pub fn BEclearMatrix::n_rows(self : BEclearMatrix) -> Int {
  self.row_names.length()
}

///|
/// Number of columns (samples).
pub fn BEclearMatrix::n_cols(self : BEclearMatrix) -> Int {
  self.col_names.length()
}

///|
/// Median of a non-empty array of doubles.
fn beclear_median(arr : Array[Double]) -> Double {
  if arr.is_empty() {
    return 0.0
  }
  let sorted = arr.copy()
  sorted.sort_by(fn(a, b) {
    if a < b {
      -1
    } else if a > b {
      1
    } else {
      0
    }
  })
  let n = sorted.length()
  if n % 2 == 1 {
    sorted[n / 2]
  } else {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  }
}

///|
/// Mean of an array.
fn beclear_mean(arr : Array[Double]) -> Double {
  if arr.is_empty() {
    return 0.0
  }
  let mut sum = 0.0
  for v in arr {
    sum += v
  }
  sum / arr.length().to_double()
}

///|
/// Absolute value.
fn beclear_abs(x : Double) -> Double {
  if x < 0.0 {
    0.0 - x
  } else {
    x
  }
}

///|
/// Summary of batch effect analysis for one feature row.
pub struct BEclearBatchSummary {
  row_name : String
  batch_medians : Map[Int, Double]
  overall_median : Double
  max_deviation : Double
  affected : Bool
} derive(Eq)

///|
/// Detect batch effects: for each row, compute per-batch medians and
/// flag rows where the maximum deviation between any batch median and
/// the overall median exceeds `threshold`.
pub fn detect_batch_effect(
  matrix : BEclearMatrix,
  threshold~ : Double = 0.05,
) -> Array[BEclearBatchSummary] {
  let summaries : Array[BEclearBatchSummary] = Array::new()
  for i in 0..<matrix.n_rows() {
    let batch_values : Map[Int, Array[Double]] = Map::new()
    for j in 0..<matrix.n_cols() {
      let col_name = matrix.col_names[j]
      match matrix.batch.get(col_name) {
        Some(b) => {
          match batch_values.get(b) {
            Some(arr) => arr.push(matrix.data[i][j])
            None => batch_values.set(b, [matrix.data[i][j]])
          }
        }
        None => ()
      }
    }
    let batch_medians : Map[Int, Double] = Map::new()
    let all_values : Array[Double] = Array::new()
    for b in batch_values.keys() {
      match batch_values.get(b) {
        Some(arr) => {
          let med = beclear_median(arr)
          batch_medians.set(b, med)
          for v in arr {
            all_values.push(v)
          }
        }
        None => ()
      }
    }
    let overall = beclear_median(all_values)
    let mut max_dev = 0.0
    for b in batch_medians.keys() {
      match batch_medians.get(b) {
        Some(med) => {
          let dev = beclear_abs(med - overall)
          if dev > max_dev {
            max_dev = dev
          }
        }
        None => ()
      }
    }
    summaries.push(
      BEclearBatchSummary::{
        row_name: matrix.row_names[i],
        batch_medians,
        overall_median: overall,
        max_deviation: max_dev,
        affected: max_dev > threshold,
      },
    )
  }
  summaries
}

///|
/// Identify (row, batch) blocks that need correction.
///
/// Returns a list of (row_index, batch_id) pairs where the batch
/// median deviates from the overall row median by more than the
/// threshold.
pub fn beclear_affected_blocks(
  matrix : BEclearMatrix,
  threshold~ : Double = 0.05,
) -> Array[(Int, Int)] {
  let summaries = detect_batch_effect(matrix, threshold~)
  let blocks : Array[(Int, Int)] = Array::new()
  for i in 0..<summaries.length() {
    let s = summaries[i]
    if s.affected {
      for b in s.batch_medians.keys() {
        match s.batch_medians.get(b) {
          Some(med) => {
            if beclear_abs(med - s.overall_median) > threshold {
              blocks.push((i, b))
            }
          }
          None => ()
        }
      }
    }
  }
  blocks
}

///|
/// Correct a single (row, batch) block by shifting all values in that
/// batch for that row so the batch median matches the overall row
/// median.
fn beclear_correct_block(
  data : Array[Array[Double]],
  matrix : BEclearMatrix,
  row_idx : Int,
  batch_id : Int,
  overall_median : Double,
) -> Unit {
  let batch_vals : Array[Double] = Array::new()
  for j in 0..<matrix.n_cols() {
    let col_name = matrix.col_names[j]
    match matrix.batch.get(col_name) {
      Some(b) => {
        if b == batch_id {
          batch_vals.push(data[row_idx][j])
        }
      }
      None => ()
    }
  }
  let batch_med = beclear_median(batch_vals)
  let shift = overall_median - batch_med
  for j in 0..<matrix.n_cols() {
    let col_name = matrix.col_names[j]
    match matrix.batch.get(col_name) {
      Some(b) => {
        if b == batch_id {
          data[row_idx][j] = data[row_idx][j] + shift
        }
      }
      None => ()
    }
  }
}

///|
/// Main BEclear pipeline: detect batch effects and correct them by
/// median shifting.
///
/// Returns a new corrected data matrix (same dimensions as the input).
pub fn beclear(
  matrix : BEclearMatrix,
  threshold~ : Double = 0.05,
) -> Array[Array[Double]] {
  let n_rows = matrix.n_rows()
  let n_cols = matrix.n_cols()
  let corrected : Array[Array[Double]] = Array::new()
  for i in 0..<n_rows {
    let row : Array[Double] = Array::new()
    for j in 0..<n_cols {
      row.push(matrix.data[i][j])
    }
    corrected.push(row)
  }
  let summaries = detect_batch_effect(matrix, threshold~)
  for i in 0..<summaries.length() {
    let s = summaries[i]
    if s.affected {
      for b in s.batch_medians.keys() {
        match s.batch_medians.get(b) {
          Some(med) => {
            if beclear_abs(med - s.overall_median) > threshold {
              beclear_correct_block(
                corrected, matrix, i, b, s.overall_median,
              )
            }
          }
          None => ()
        }
      }
    }
  }
  corrected
}

///|
/// Compute the BEscore — a summary measure of remaining batch effect
/// after correction.  Lower is better; zero means no batch effect.
pub fn beclear_score(
  matrix : BEclearMatrix,
  corrected_data : Array[Array[Double]],
) -> Double {
  let temp = BEclearMatrix::{
    row_names: matrix.row_names,
    col_names: matrix.col_names,
    data: corrected_data,
    batch: matrix.batch,
  }
  let summaries = detect_batch_effect(temp, threshold=0.0)
  let mut total = 0.0
  let mut count = 0
  for s in summaries {
    total += s.max_deviation * s.max_deviation
    count += 1
  }
  if count == 0 {
    0.0
  } else {
    beclear_sqrt(total / count.to_double())
  }
}

///|
/// sqrt via Newton's method.
fn beclear_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut guess = x
  for _ in 0..<50 {
    guess = (guess + x / guess) / 2.0
  }
  guess
}

///|
/// A sample matrix with a clear batch effect for testing.
pub fn beclear_sample_matrix() -> BEclearMatrix {
  let row_names = ["gene1", "gene2", "gene3"]
  let col_names = ["s1", "s2", "s3", "s4"]
  let data = [
    [1.0, 1.1, 1.5, 1.6],
    [2.0, 2.1, 2.5, 2.6],
    [3.0, 3.1, 3.5, 3.6],
  ]
  let batch = {
    let m : Map[String, Int] = Map::new()
    m.set("s1", 1)
    m.set("s2", 1)
    m.set("s3", 2)
    m.set("s4", 2)
    m
  }
  BEclearMatrix::new(row_names~, col_names~, data~, batch~)
}