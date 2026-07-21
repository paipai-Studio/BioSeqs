///|
pub struct SingleCellData {
  counts : Array[Array[Double]]
  cell_names : Array[String]
  gene_names : Array[String]
}

///|
pub fn SingleCellData::new(
  counts : Array[Array[Double]],
  cell_names : Array[String],
  gene_names : Array[String],
) -> SingleCellData {
  SingleCellData::{ counts, cell_names, gene_names }
}

///|
pub struct DoubletScore {
  cell_name : String
  score : Double
  is_doublet : Bool
}

///|
pub fn DoubletScore::new(
  cell_name : String,
  score : Double,
  is_doublet : Bool,
) -> DoubletScore {
  DoubletScore::{ cell_name, score, is_doublet }
}

///|
pub fn scdf_compute_doublet_score(
  data : SingleCellData,
  k : Int,
) -> Array[DoubletScore] {
  let n_cells = data.cell_names.length()
  let n_genes = data.gene_names.length()

  let scores : Array[DoubletScore] = Array::new()

  let mut i = 0
  while i < n_cells {
    let mut score = 0.0
    let mut j = 0
    while j < n_cells && j != i {
      let distance = scdf_compute_distance(data.counts[i], data.counts[j])
      if distance < 1.0 {
        score = score + (1.0 - distance)
      }
      j = j + 1
    }

    score = score / (n_cells - 1).to_double()
    let is_doublet = score > 0.3

    scores.push(DoubletScore::new(data.cell_names[i], score, is_doublet))
    i = i + 1
  }

  scores
}

///|
pub fn scdf_compute_distance(
  vec1 : Array[Double],
  vec2 : Array[Double],
) -> Double {
  let mut sum_sq = 0.0
  let mut i = 0
  while i < vec1.length() && i < vec2.length() {
    let diff = vec1[i] - vec2[i]
    sum_sq = sum_sq + diff * diff
    i = i + 1
  }
  sum_sq.sqrt()
}

///|
pub fn scdf_find_nearest_neighbors(
  data : SingleCellData,
  cell_idx : Int,
  k : Int,
) -> Array[(String, Double)] {
  let n_cells = data.cell_names.length()
  let distances : Array[(String, Double)] = Array::new()

  let mut i = 0
  while i < n_cells {
    if i != cell_idx {
      let dist = scdf_compute_distance(data.counts[cell_idx], data.counts[i])
      distances.push((data.cell_names[i], dist))
    }
    i = i + 1
  }

  let sorted = scdf_sort_by_distance(distances)
  let result : Array[(String, Double)] = Array::new()
  let mut j = 0
  while j < k && j < sorted.length() {
    result.push(sorted[j])
    j = j + 1
  }
  result
}

///|
pub fn scdf_sort_by_distance(
  distances : Array[(String, Double)],
) -> Array[(String, Double)] {
  let sorted = distances.copy()
  let n = sorted.length()
  if n <= 1 {
    return sorted
  }
  
  let mut i = 0
  while i < n - 1 {
    let min_idx = scdf_find_min(sorted, i, n)
    if min_idx != i {
      let temp = sorted[i]
      sorted[i] = sorted[min_idx]
      sorted[min_idx] = temp
    }
    i = i + 1
  }
  sorted
}

///|
fn scdf_find_min(
  arr : Array[(String, Double)],
  start : Int,
  end : Int,
) -> Int {
  let mut min_idx = start
  let mut i = start + 1
  while i < end {
    if arr[i].1 < arr[min_idx].1 {
      min_idx = i
    }
    i = i + 1
  }
  min_idx
}

///|
pub fn scdf_detect_doublets(
  scores : Array[DoubletScore],
  threshold : Double,
) -> Array[String] {
  let doublets : Array[String] = Array::new()

  let mut i = 0
  while i < scores.length() {
    if scores[i].score > threshold {
      doublets.push(scores[i].cell_name)
    }
    i = i + 1
  }

  doublets
}

///|
pub fn scdf_doublet_summary(scores : Array[DoubletScore]) -> String {
  let mut result = "Doublet Detection Summary\n"
  result = result + "========================\n"
  result = result + "Total cells: " + scores.length().to_string() + "\n"

  let mut n_doublets = 0
  let mut sum_score = 0.0
  let mut min_score = 1.0
  let mut max_score = 0.0

  let mut i = 0
  while i < scores.length() {
    if scores[i].is_doublet {
      n_doublets = n_doublets + 1
    }
    sum_score = sum_score + scores[i].score
    if scores[i].score < min_score {
      min_score = scores[i].score
    }
    if scores[i].score > max_score {
      max_score = scores[i].score
    }
    i = i + 1
  }

  result = result + "Detected doublets: " + n_doublets.to_string() + "\n"
  result = result +
    "Doublet rate: " +
    (n_doublets.to_double() / scores.length().to_double() * 100.0).to_string() +
    "%\n"
  result = result +
    "Mean score: " +
    (sum_score / scores.length().to_double()).to_string() +
    "\n"
  result = result + "Min score: " + min_score.to_string() + "\n"
  result = result + "Max score: " + max_score.to_string() + "\n"

  result
}

///|
pub fn scdf_filter_doublets(
  data : SingleCellData,
  scores : Array[DoubletScore],
) -> SingleCellData {
  let new_counts : Array[Array[Double]] = Array::new()
  let new_cell_names : Array[String] = Array::new()

  let mut i = 0
  while i < scores.length() {
    if !scores[i].is_doublet {
      new_counts.push(data.counts[i])
      new_cell_names.push(data.cell_names[i])
    }
    i = i + 1
  }

  SingleCellData::new(new_counts, new_cell_names, data.gene_names)
}

///|
pub fn scdf_compute_pca(
  data : SingleCellData,
  n_components : Int,
) -> Array[Array[Double]] {
  let n_cells = data.cell_names.length()
  let n_genes = data.gene_names.length()

  let col_means : Array[Double] = Array::new()
  let mut j = 0
  while j < n_genes {
    let mut sum = 0.0
    let mut i = 0
    while i < n_cells {
      sum = sum + data.counts[i][j]
      i = i + 1
    }
    col_means.push(sum / n_cells.to_double())
    j = j + 1
  }

  let centered : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_cells {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n_genes {
      row.push(data.counts[i][j] - col_means[j])
      j = j + 1
    }
    centered.push(row)
    i = i + 1
  }

  let cov_matrix : Array[Array[Double]] = Array::new()
  let mut k = 0
  while k < n_genes {
    let row : Array[Double] = Array::new()
    let mut l = 0
    while l < n_genes {
      let mut sum = 0.0
      let mut m = 0
      while m < n_cells {
        sum = sum + centered[m][k] * centered[m][l]
        m = m + 1
      }
      row.push(sum / (n_cells - 1).to_double())
      l = l + 1
    }
    cov_matrix.push(row)
    k = k + 1
  }

  let eigenvalues = scdf_compute_eigenvalues(cov_matrix)
  let pca_result : Array[Array[Double]] = Array::new()

  let mut i = 0
  while i < n_cells {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n_components && j < eigenvalues.length() {
      let mut sum = 0.0
      let mut k = 0
      while k < n_genes {
        sum = sum + centered[i][k] * cov_matrix[k][j]
        k = k + 1
      }
      row.push(sum)
      j = j + 1
    }
    pca_result.push(row)
    i = i + 1
  }

  pca_result
}

///|
pub fn scdf_compute_eigenvalues(matrix : Array[Array[Double]]) -> Array[Double] {
  let n = matrix.length()
  let eigenvalues : Array[Double] = Array::new()

  let mut i = 0
  while i < n && i < 10 {
    let v : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      v.push(1.0)
      j = j + 1
    }

    let mut iteration = 0
    while iteration < 50 {
      let new_v : Array[Double] = Array::new()
      let mut j = 0
      while j < n {
        let mut sum = 0.0
        let mut k = 0
        while k < n {
          sum = sum + matrix[j][k] * v[k]
          k = k + 1
        }
        new_v.push(sum)
        j = j + 1
      }

      let mut norm = 0.0
      let mut j = 0
      while j < n {
        norm = norm + new_v[j] * new_v[j]
        j = j + 1
      }
      norm = norm.sqrt()

      let mut j = 0
      while j < n {
        v[j] = new_v[j] / norm
        j = j + 1
      }

      iteration = iteration + 1
    }

    let mut eigenvalue = 0.0
    let mut j = 0
    while j < n {
      let mut sum = 0.0
      let mut k = 0
      while k < n {
        sum = sum + matrix[j][k] * v[k]
        k = k + 1
      }
      eigenvalue = eigenvalue + v[j] * sum
      j = j + 1
    }
    eigenvalues.push(eigenvalue)

    let mut j = 0
    while j < n {
      let mut k = 0
      while k < n {
        matrix[j][k] = matrix[j][k] - eigenvalue * v[j] * v[k]
        k = k + 1
      }
      j = j + 1
    }

    i = i + 1
  }

  eigenvalues
}

///|
pub fn scdf_create_example_data() -> SingleCellData {
  let cell_names : Array[String] = Array::new()
  let mut i = 0
  while i < 50 {
    cell_names.push("Cell_" + i.to_string())
    i = i + 1
  }

  let gene_names = [
    "TP53", "BRCA1", "BRCA2", "ATM", "CDK4", "MYC", "EGFR", "ACTB", "GAPDH", "RB1",
  ]

  let counts : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < 50 {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < 10 {
      if i < 10 {
        row.push((j * 10 + i * 2).to_double())
      } else if i < 20 {
        row.push((j * 15 + i * 3).to_double())
      } else if i < 30 {
        row.push((j * 8 + i * 2).to_double())
      } else if i < 40 {
        row.push((j * 12 + i * 4).to_double())
      } else {
        row.push((j * 5 + i * 1).to_double())
      }
      j = j + 1
    }
    counts.push(row)
    i = i + 1
  }

  SingleCellData::new(counts, cell_names, gene_names)
}
