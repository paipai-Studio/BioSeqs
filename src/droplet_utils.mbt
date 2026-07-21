///|
/// DropletUtils module, inspired by Bioconductor DropletUtils.
/// Provides functions for single-cell RNA-seq data analysis, including
/// empty droplet detection and barcode processing.

///|
pub struct DropletStats {
  total_counts : Array[Double]
  n_genes : Array[Int]
  total_features_by_barcode : Array[(String, Double, Int)]
}

///|
pub struct EmptyDropsResult {
  barcode : Array[String]
  total : Array[Double]
  log_prob_empty : Array[Double]
  p_value : Array[Double]
  fdr : Array[Double]
  is_cell : Array[Bool]
}

///|
pub fn calculate_droplet_stats(
  counts : Array[Array[Int]],
  barcodes : Array[String],
) -> DropletStats {
  let n_barcodes = counts.length()
  let n_genes = if n_barcodes > 0 { counts[0].length() } else { 0 }

  let total_counts : Array[Double] = Array::new()
  let gene_counts : Array[Int] = Array::new()
  let features : Array[(String, Double, Int)] = Array::new()

  let mut i = 0
  while i < n_barcodes {
    let mut total = 0.0
    let mut detected = 0

    let mut j = 0
    while j < n_genes {
      let count = counts[i][j]
      total = total + count.to_double()
      if count > 0 {
        detected = detected + 1
      }
      j = j + 1
    }

    total_counts.push(total)
    gene_counts.push(detected)

    if i < barcodes.length() {
      features.push((barcodes[i], total, detected))
    }

    i = i + 1
  }

  DropletStats::{
    total_counts,
    n_genes: gene_counts,
    total_features_by_barcode: features,
  }
}

///|
pub fn barcode_ranking(total_counts : Array[Double]) -> Array[(Double, Int)] {
  let ranked : Array[(Double, Int)] = Array::new()

  let mut i = 0
  while i < total_counts.length() {
    ranked.push((total_counts[i], i))
    i = i + 1
  }

  let mut j = 0
  while j < ranked.length() - 1 {
    let mut k = j + 1
    while k < ranked.length() {
      if ranked[j].0 < ranked[k].0 {
        let tmp = ranked[j]
        ranked[j] = ranked[k]
        ranked[k] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }

  ranked
}

///|
pub fn find_knee_point(ranked_counts : Array[(Double, Int)]) -> Int {
  if ranked_counts.length() < 3 {
    return 0
  }

  let max_count = ranked_counts[0].0
  let min_count = ranked_counts[ranked_counts.length() - 1].0

  let mut knee_index = 0
  let mut max_curvature = 0.0

  let mut i = 1
  while i < ranked_counts.length() - 1 {
    let x1 = (i - 1).to_double() / (ranked_counts.length() - 1).to_double()
    let y1 = (ranked_counts[i - 1].0 - min_count) / (max_count - min_count)

    let x2 = i.to_double() / (ranked_counts.length() - 1).to_double()
    let y2 = (ranked_counts[i].0 - min_count) / (max_count - min_count)

    let x3 = (i + 1).to_double() / (ranked_counts.length() - 1).to_double()
    let y3 = (ranked_counts[i + 1].0 - min_count) / (max_count - min_count)

    let a = ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)).sqrt()
    let b = ((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2)).sqrt()
    let c = ((x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1)).sqrt()

    let curvature = if a > 0.0 && b > 0.0 {
      (a * a + b * b - c * c) / (2.0 * a * b)
    } else {
      0.0
    }

    if curvature > max_curvature {
      max_curvature = curvature
      knee_index = i
    }

    i = i + 1
  }

  knee_index
}

///|
pub fn empty_drops(
  counts : Array[Array[Int]],
  barcodes : Array[String],
  lower : Int,
  niters : Int,
) -> EmptyDropsResult {
  let n_barcodes = counts.length()
  let n_genes = if n_barcodes > 0 { counts[0].length() } else { 0 }

  let stats = calculate_droplet_stats(counts, barcodes)

  let ambient_barcodes : Array[Array[Int]] = Array::new()
  let ambient_indices : Array[Int] = Array::new()

  let mut i = 0
  while i < n_barcodes {
    if stats.total_counts[i] <= lower.to_double() {
      ambient_barcodes.push(counts[i])
      ambient_indices.push(i)
    }
    i = i + 1
  }

  let ambient_profile : Array[Double] = Array::make(n_genes, 0.0)
  let mut total_ambient = 0.0

  i = 0
  while i < ambient_barcodes.length() {
    let mut j = 0
    while j < n_genes {
      ambient_profile[j] = ambient_profile[j] +
        ambient_barcodes[i][j].to_double()
      total_ambient = total_ambient + ambient_barcodes[i][j].to_double()
      j = j + 1
    }
    i = i + 1
  }

  i = 0
  while i < n_genes {
    if total_ambient > 0.0 {
      ambient_profile[i] = ambient_profile[i] / total_ambient
    }
    i = i + 1
  }

  let barcodes_result : Array[String] = Array::new()
  let total_result : Array[Double] = Array::new()
  let log_prob_empty_result : Array[Double] = Array::new()
  let p_value_result : Array[Double] = Array::new()

  i = 0
  while i < n_barcodes {
    if stats.total_counts[i] > lower.to_double() {
      let barcode = if i < barcodes.length() {
        barcodes[i]
      } else {
        i.to_string()
      }
      barcodes_result.push(barcode)
      total_result.push(stats.total_counts[i])

      let total_umi = stats.total_counts[i].to_int()

      let mut log_prob = 0.0
      let mut j = 0
      while j < n_genes {
        let count = counts[i][j]
        if count > 0 && ambient_profile[j] > 0.0 {
          log_prob = log_prob + count.to_double() * @math.ln(ambient_profile[j])
        }
        j = j + 1
      }

      let mut fact = 1.0
      let mut f = 1
      while f <= total_umi && f <= 20 {
        fact = fact * f.to_double()
        f = f + 1
      }
      log_prob = log_prob - total_umi.to_double() + @math.ln(fact)

      let mut n_extreme = 0
      let mut iter = 0
      while iter < niters {
        let mut simulated_log_prob = 0.0
        let mut remaining = total_umi

        let mut j = 0
        while j < n_genes && remaining > 0 {
          let prob = ambient_profile[j]
          if prob > 0.0 && remaining > 0 {
            let expected = prob * remaining.to_double()
            let sampled = if expected < 100.0 {
              let mut sum = 0
              let mut k = 0
              while k < remaining {
                let u = (k * 7 + i * 13).to_double() % 10.0 / 10.0
                if u < prob {
                  sum = sum + 1
                }
                k = k + 1
              }
              sum
            } else {
              expected.to_int()
            }

            if sampled > 0 {
              simulated_log_prob = simulated_log_prob +
                sampled.to_double() * @math.ln(prob)
            }
            remaining = remaining - sampled
          }
          j = j + 1
        }

        if simulated_log_prob <= log_prob {
          n_extreme = n_extreme + 1
        }

        iter = iter + 1
      }

      log_prob_empty_result.push(log_prob)

      let p_val = (n_extreme + 1).to_double() / (niters + 1).to_double()
      p_value_result.push(p_val)
    }
    i = i + 1
  }

  let fdr_result = benjamini_hochberg(p_value_result)

  let is_cell_result : Array[Bool] = Array::new()
  i = 0
  while i < fdr_result.length() {
    is_cell_result.push(fdr_result[i] <= 0.01)
    i = i + 1
  }

  EmptyDropsResult::{
    barcode: barcodes_result,
    total: total_result,
    log_prob_empty: log_prob_empty_result,
    p_value: p_value_result,
    fdr: fdr_result,
    is_cell: is_cell_result,
  }
}

///|
fn benjamini_hochberg(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()

  let indexed : Array[(Double, Int)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((p_values[i], i))
    i = i + 1
  }

  let mut j = 0
  while j < indexed.length() - 1 {
    let mut k = j + 1
    while k < indexed.length() {
      if indexed[j].0 > indexed[k].0 {
        let tmp = indexed[j]
        indexed[j] = indexed[k]
        indexed[k] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }

  let fdr : Array[Double] = Array::make(n, 1.0)
  let mut min_fdr = 1.0

  i = n - 1
  while i >= 0 {
    let p = indexed[i].0
    let rank = (i + 1).to_double()
    let current_fdr = p * n.to_double() / rank

    if current_fdr < min_fdr {
      min_fdr = current_fdr
    } else {
      min_fdr = min_fdr
    }

    fdr[indexed[i].1] = min_fdr

    i = i - 1
  }

  fdr
}

///|
pub fn filter_cells_by_empty_drops(
  counts : Array[Array[Int]],
  barcodes : Array[String],
  result : EmptyDropsResult,
) -> (Array[Array[Int]], Array[String]) {
  let filtered_counts : Array[Array[Int]] = Array::new()
  let filtered_barcodes : Array[String] = Array::new()

  let barcode_map : Map[String, Bool] = Map(
    [],
    capacity=result.barcode.length(),
  )
  let mut i = 0
  while i < result.barcode.length() {
    if result.is_cell[i] {
      barcode_map[result.barcode[i]] = true
    }
    i = i + 1
  }

  i = 0
  while i < counts.length() {
    let barcode = if i < barcodes.length() {
      barcodes[i]
    } else {
      i.to_string()
    }
    let is_cell = barcode_map.get(barcode)
    if is_cell.is_some() && is_cell.unwrap() {
      filtered_counts.push(counts[i])
      filtered_barcodes.push(barcode)
    }
    i = i + 1
  }

  (filtered_counts, filtered_barcodes)
}

///|
pub fn estimate_cell_number(total_counts : Array[Double]) -> Int {
  let ranked = barcode_ranking(total_counts)
  let knee = find_knee_point(ranked)
  knee + 1
}

///|
pub fn create_example_droplet_data() -> (Array[Array[Int]], Array[String]) {
  let n_barcodes = 100
  let n_genes = 500

  let counts : Array[Array[Int]] = Array::new()
  let barcodes : Array[String] = Array::new()

  let mut i = 0
  while i < n_barcodes {
    let row : Array[Int] = Array::new()
    let s = i.to_string()
    let pad = "00000000"
    let padded = if s.length() < 8 {
      pad[0:8 - s.length()].to_string() + s
    } else {
      s
    }
    let barcode = "BC" + padded
    barcodes.push(barcode)

    let mut j = 0
    while j < n_genes {
      let mut count = 0
      if i < 30 {
        if (j * 7 + i * 13).to_double() % 10.0 / 10.0 < 0.3 {
          count = ((j * 7 + i * 13).to_double() % 10.0 / 10.0 * 20.0).to_int() + 1
        }
      } else if (j * 7 + i * 13).to_double() % 10.0 / 10.0 < 0.05 {
        count = ((j * 7 + i * 13).to_double() % 10.0 / 10.0 * 5.0).to_int()
      }
      row.push(count)
      j = j + 1
    }

    counts.push(row)
    i = i + 1
  }

  (counts, barcodes)
}
