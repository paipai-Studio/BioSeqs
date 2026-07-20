///|
/// ComplexHeatmap - MoonBit implementation of Bioconductor ComplexHeatmap package

///|
pub enum AnnotationType {
  RowAnnotation
  ColumnAnnotation
  CellAnnotation
}

///|
pub struct HeatmapData {
  matrix : Array[Array[Double]]
  row_names : Array[String]
  col_names : Array[String]
  row_dend : Array[Array[Double]]
  col_dend : Array[Array[Double]]
}

///|
pub struct HeatmapAnnotation {
  annot_type : AnnotationType
  labels : Array[String]
  colors : Array[String]
  values : Array[Double]
}

///|
pub struct ColorMap {
  colors : Array[String]
  breaks : Array[Double]
  na_color : String
}

///|
pub struct HeatmapResult {
  text : String
  row_order : Array[Int]
  col_order : Array[Int]
  row_clusters : Array[Int]
  col_clusters : Array[Int]
}

///|
pub struct HeatmapAnnotationResult {
  row_annot : String
  col_annot : String
  legend : String
}

///|
pub struct HeatmapOptions {
  row_title : String
  col_title : String
  show_row_names : Bool
  show_col_names : Bool
  show_legend : Bool
  row_cluster : Bool
  col_cluster : Bool
  cluster_method : String
  scale : String
  center : Bool
  na_color : String
  border_color : String
  row_gap : Int
  col_gap : Int
}

///|
pub fn HeatmapData::new(
  matrix : Array[Array[Double]],
  row_names : Array[String],
  col_names : Array[String],
) -> HeatmapData {
  HeatmapData::{ matrix, row_names, col_names, row_dend: [], col_dend: [] }
}

///|
pub fn HeatmapAnnotation::new(
  annot_type : AnnotationType,
  labels : Array[String],
  colors : Array[String],
) -> HeatmapAnnotation {
  HeatmapAnnotation::{ annot_type, labels, colors, values: [] }
}

///|
pub fn ColorMap::new(
  colors : Array[String],
  breaks : Array[Double],
) -> ColorMap {
  ColorMap::{ colors, breaks, na_color: "#FFFFFF" }
}

///|
pub fn HeatmapOptions::new() -> HeatmapOptions {
  HeatmapOptions::{
    row_title: "",
    col_title: "",
    show_row_names: true,
    show_col_names: true,
    show_legend: true,
    row_cluster: true,
    col_cluster: true,
    cluster_method: "complete",
    scale: "none",
    center: false,
    na_color: "#FFFFFF",
    border_color: "#000000",
    row_gap: 0,
    col_gap: 0,
  }
}

///|
pub fn ch_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else {
    let mut guess = x / 2.0
    let mut prev = 0.0
    while (guess - prev).abs() > 0.0000001 {
      prev = guess
      guess = (guess + x / guess) / 2.0
    }
    guess
  }
}

///|
pub fn ch_ln(x : Double) -> Double {
  if x <= 0.0 {
    -1000000.0
  } else {
    let y = (x - 1.0) / (x + 1.0)
    let mut sum = y
    let mut y_power = y
    let mut n = 1
    while n < 100 {
      y_power = y_power * y * y
      sum = sum + y_power / (2.0 * n.to_double() + 1.0)
      n = n + 1
    }
    2.0 * sum
  }
}

///|
pub fn ch_exp(x : Double) -> Double {
  if x > 700.0 {
    1000000000000000000000000000000000.0
  } else if x < -700.0 {
    0.0
  } else {
    let mut sum = 1.0
    let mut term = 1.0
    let mut n = 1
    while n < 100 {
      term = term * x / n.to_double()
      sum = sum + term
      n = n + 1
    }
    sum
  }
}

///|
pub fn ch_mean(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    0.0
  } else {
    let mut sum = 0.0
    let mut i = 0
    while i < arr.length() {
      sum = sum + arr[i]
      i = i + 1
    }
    sum / arr.length().to_double()
  }
}

///|
pub fn ch_variance(arr : Array[Double]) -> Double {
  if arr.length() < 2 {
    0.0
  } else {
    let mean_val = ch_mean(arr)
    let mut sum_sq = 0.0
    let mut i = 0
    while i < arr.length() {
      let diff = arr[i] - mean_val
      sum_sq = sum_sq + diff * diff
      i = i + 1
    }
    sum_sq / (arr.length() - 1).to_double()
  }
}

///|
pub fn ch_std(arr : Array[Double]) -> Double {
  ch_sqrt(ch_variance(arr))
}

///|
pub fn ch_min(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    0.0
  } else {
    let mut min_val = arr[0]
    let mut i = 1
    while i < arr.length() {
      if arr[i] < min_val {
        min_val = arr[i]
      }
      i = i + 1
    }
    min_val
  }
}

///|
pub fn ch_max(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    0.0
  } else {
    let mut max_val = arr[0]
    let mut i = 1
    while i < arr.length() {
      if arr[i] > max_val {
        max_val = arr[i]
      }
      i = i + 1
    }
    max_val
  }
}

///|
pub fn AnnotationType::to_string(self : AnnotationType) -> String {
  match self {
    RowAnnotation => "RowAnnotation"
    ColumnAnnotation => "ColumnAnnotation"
    CellAnnotation => "CellAnnotation"
  }
}

///|
pub fn HeatmapData::to_string(self : HeatmapData) -> String {
  let mut s = "HeatmapData {\n"
  s = s + "  Rows: " + self.row_names.length().to_string() + "\n"
  s = s + "  Cols: " + self.col_names.length().to_string() + "\n"
  s = s + "}"
  s
}

///|
pub fn HeatmapResult::to_string(self : HeatmapResult) -> String {
  let mut s = "HeatmapResult {\n"
  s = s + "  Row order: " + self.row_order.length().to_string() + " elements\n"
  s = s + "  Col order: " + self.col_order.length().to_string() + " elements\n"
  s = s + "}"
  s
}

///|
pub fn HeatmapAnnotation::to_string(self : HeatmapAnnotation) -> String {
  let mut s = "HeatmapAnnotation {\n"
  s = s + "  Type: " + self.annot_type.to_string() + "\n"
  s = s + "  Labels: " + self.labels.length().to_string() + "\n"
  s = s + "}"
  s
}

///|
pub fn ColorMap::to_string(self : ColorMap) -> String {
  let mut s = "ColorMap {\n"
  s = s + "  Colors: " + self.colors.length().to_string() + "\n"
  s = s + "  Breaks: " + self.breaks.length().to_string() + "\n"
  s = s + "}"
  s
}

///|
pub fn ch_euclidean_distance(v1 : Array[Double], v2 : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < v1.length() && i < v2.length() {
    let diff = v1[i] - v2[i]
    sum = sum + diff * diff
    i = i + 1
  }
  ch_sqrt(sum)
}

///|
pub fn ch_pearson_distance(v1 : Array[Double], v2 : Array[Double]) -> Double {
  let n = if v1.length() < v2.length() { v1.length() } else { v2.length() }
  if n < 2 {
    1.0
  } else {
    let mean1 = ch_mean(v1)
    let mean2 = ch_mean(v2)
    let mut cov = 0.0
    let mut var1 = 0.0
    let mut var2 = 0.0
    let mut i = 0
    while i < n {
      let diff1 = v1[i] - mean1
      let diff2 = v2[i] - mean2
      cov = cov + diff1 * diff2
      var1 = var1 + diff1 * diff1
      var2 = var2 + diff2 * diff2
      i = i + 1
    }
    let denom = ch_sqrt(var1) * ch_sqrt(var2)
    if denom == 0.0 {
      1.0
    } else {
      1.0 - cov / denom
    }
  }
}

///|
pub fn ch_cosine_distance(v1 : Array[Double], v2 : Array[Double]) -> Double {
  let mut dot = 0.0
  let mut norm1 = 0.0
  let mut norm2 = 0.0
  let mut i = 0
  while i < v1.length() && i < v2.length() {
    dot = dot + v1[i] * v2[i]
    norm1 = norm1 + v1[i] * v1[i]
    norm2 = norm2 + v2[i] * v2[i]
    i = i + 1
  }
  let denom = ch_sqrt(norm1) * ch_sqrt(norm2)
  if denom == 0.0 {
    1.0
  } else {
    1.0 - dot / denom
  }
}

///|
pub fn ch_distance_matrix(
  matrix : Array[Array[Double]],
  dist_method : String,
) -> Array[Array[Double]] {
  let n = matrix.length()
  let dist_matrix : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let dist = if dist_method == "pearson" {
        ch_pearson_distance(matrix[i], matrix[j])
      } else if dist_method == "cosine" {
        ch_cosine_distance(matrix[i], matrix[j])
      } else {
        ch_euclidean_distance(matrix[i], matrix[j])
      }
      dist_matrix[i][j] = dist
      dist_matrix[j][i] = dist
      j = j + 1
    }
    i = i + 1
  }
  dist_matrix
}

///|
pub fn ch_hclust(
  dist_matrix : Array[Array[Double]],
  clust_method : String,
) -> Array[Int] {
  let n = dist_matrix.length()
  let mut clusters : Array[Array[Int]] = Array::make(n, Array::new())
  let mut i = 0
  while i < n {
    clusters[i] = [i]
    i = i + 1
  }

  while clusters.length() > 1 {
    let mut min_dist = 1000000000000000000000000000000000.0
    let mut min_i = 0
    let mut min_j = 1

    let m = clusters.length()
    let mut a = 0
    while a < m {
      let mut b = a + 1
      while b < m {
        let mut dist = 0.0
        if clust_method == "single" {
          dist = 1000000000000000000000000000000000.0
          let mut x = 0
          while x < clusters[a].length() {
            let mut y = 0
            while y < clusters[b].length() {
              let d = dist_matrix[clusters[a][x]][clusters[b][y]]
              if d < dist {
                dist = d
              }
              y = y + 1
            }
            x = x + 1
          }
        } else if clust_method == "average" {
          let mut sum = 0.0
          let mut count = 0
          let mut x = 0
          while x < clusters[a].length() {
            let mut y = 0
            while y < clusters[b].length() {
              sum = sum + dist_matrix[clusters[a][x]][clusters[b][y]]
              count = count + 1
              y = y + 1
            }
            x = x + 1
          }
          dist = sum / count.to_double()
        } else {
          dist = 0.0
          let mut x = 0
          while x < clusters[a].length() {
            let mut y = 0
            while y < clusters[b].length() {
              let d = dist_matrix[clusters[a][x]][clusters[b][y]]
              if d > dist {
                dist = d
              }
              y = y + 1
            }
            x = x + 1
          }
        }

        if dist < min_dist {
          min_dist = dist
          min_i = a
          min_j = b
        }
        b = b + 1
      }
      a = a + 1
    }

    let merged = clusters[min_i] + clusters[min_j]
    let new_clusters : Array[Array[Int]] = Array::make(
      clusters.length() - 1,
      Array::new(),
    )
    let mut idx = 0
    let mut k = 0
    while k < clusters.length() {
      if k != min_i && k != min_j {
        new_clusters[idx] = clusters[k]
        idx = idx + 1
      }
      k = k + 1
    }
    new_clusters[idx] = merged
    clusters = new_clusters
  }

  if clusters.length() > 0 {
    clusters[0]
  } else {
    []
  }
}

///|
pub fn ch_cluster_rows(
  matrix : Array[Array[Double]],
  clust_method : String,
) -> Array[Int] {
  let dist_matrix = ch_distance_matrix(matrix, "euclidean")
  ch_hclust(dist_matrix, clust_method)
}

///|
pub fn ch_cluster_cols(
  matrix : Array[Array[Double]],
  clust_method : String,
) -> Array[Int] {
  let n_rows = matrix.length()
  let n_cols = if n_rows > 0 { matrix[0].length() } else { 0 }
  let t_matrix : Array[Array[Double]] = Array::make(
    n_cols,
    Array::make(n_rows, 0.0),
  )
  let mut i = 0
  while i < n_rows {
    let mut j = 0
    while j < n_cols {
      t_matrix[j][i] = matrix[i][j]
      j = j + 1
    }
    i = i + 1
  }
  ch_cluster_rows(t_matrix, clust_method)
}

///|
pub fn ch_row_scale(
  matrix : Array[Array[Double]],
  center : Bool,
) -> Array[Array[Double]] {
  let n_rows = matrix.length()
  let scaled : Array[Array[Double]] = Array::make(n_rows, Array::new())
  let mut i = 0
  while i < n_rows {
    let row = matrix[i]
    let mean_val = if center { ch_mean(row) } else { 0.0 }
    let std_val = ch_std(row)
    let scaled_row : Array[Double] = Array::make(row.length(), 0.0)
    let mut j = 0
    while j < row.length() {
      scaled_row[j] = if std_val > 0.0 {
        (row[j] - mean_val) / std_val
      } else {
        0.0
      }
      j = j + 1
    }
    scaled[i] = scaled_row
    i = i + 1
  }
  scaled
}

///|
pub fn ch_col_scale(
  matrix : Array[Array[Double]],
  center : Bool,
) -> Array[Array[Double]] {
  let n_rows = matrix.length()
  let n_cols = if n_rows > 0 { matrix[0].length() } else { 0 }
  let scaled : Array[Array[Double]] = Array::make(
    n_rows,
    Array::make(n_cols, 0.0),
  )
  let mut j = 0
  while j < n_cols {
    let col : Array[Double] = Array::make(n_rows, 0.0)
    let mut i = 0
    while i < n_rows {
      col[i] = matrix[i][j]
      i = i + 1
    }
    let mean_val = if center { ch_mean(col) } else { 0.0 }
    let std_val = ch_std(col)
    let mut i = 0
    while i < n_rows {
      scaled[i][j] = if std_val > 0.0 {
        (matrix[i][j] - mean_val) / std_val
      } else {
        0.0
      }
      i = i + 1
    }
    j = j + 1
  }
  scaled
}

///|
pub fn ch_scale_matrix(
  matrix : Array[Array[Double]],
  scale : String,
  center : Bool,
) -> Array[Array[Double]] {
  if scale == "row" {
    ch_row_scale(matrix, center)
  } else if scale == "column" {
    ch_col_scale(matrix, center)
  } else {
    matrix
  }
}

///|
pub fn ch_colorscale(
  value : Double,
  min_val : Double,
  max_val : Double,
  cmap : ColorMap,
) -> String {
  if cmap.colors.length() == 0 {
    return "#808080"
  }

  let range = max_val - min_val
  if range <= 0.0 {
    return cmap.colors[0]
  }

  let normalized = (value - min_val) / range
  let idx = (normalized * (cmap.colors.length() - 1).to_double()).to_int()

  if idx < 0 {
    cmap.colors[0]
  } else if idx >= cmap.colors.length() {
    cmap.colors[cmap.colors.length() - 1]
  } else {
    cmap.colors[idx]
  }
}

///|
pub fn ch_pad_right(s : String, len : Int) -> String {
  let mut result = s
  while result.length() < len {
    result = result + " "
  }
  result
}

///|
pub fn ch_draw_heatmap(
  data : HeatmapData,
  options : HeatmapOptions,
) -> HeatmapResult {
  let mut matrix = data.matrix

  if options.scale != "none" {
    matrix = ch_scale_matrix(matrix, options.scale, options.center)
  }

  let mut row_order : Array[Int] = Array::make(matrix.length(), 0)
  let mut i = 0
  while i < matrix.length() {
    row_order[i] = i
    i = i + 1
  }

  let mut col_order : Array[Int] = Array::make(
    if matrix.length() > 0 {
      matrix[0].length()
    } else {
      0
    },
    0,
  )
  let mut j = 0
  while j < col_order.length() {
    col_order[j] = j
    j = j + 1
  }

  if options.row_cluster {
    row_order = ch_cluster_rows(matrix, options.cluster_method)
  }

  if options.col_cluster {
    col_order = ch_cluster_cols(matrix, options.cluster_method)
  }

  let all_values : Array[Double] = Array::new()
  let mut x = 0
  while x < matrix.length() {
    let mut y = 0
    while y < matrix[x].length() {
      all_values.push(matrix[x][y])
      y = y + 1
    }
    x = x + 1
  }

  let min_val = ch_min(all_values)
  let max_val = ch_max(all_values)

  let cmap = ColorMap::new(
    [
      "#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090",
      "#FDAE61", "#F46D43", "#D73027", "#A50026",
    ],
    [],
  )

  let mut text = ""

  if options.show_col_names {
    text = text + "     "
    let mut k = 0
    while k < col_order.length() {
      let col_name = data.col_names[col_order[k]]
      text = text + ch_pad_right(col_name, 6)
      k = k + 1
    }
    text = text + "\n"
  }

  let row_clusters : Array[Int] = Array::make(row_order.length(), 0)
  let col_clusters : Array[Int] = Array::make(col_order.length(), 0)

  let mut m = 0
  while m < row_order.length() {
    let row_idx = row_order[m]
    if options.show_row_names {
      let row_name = data.row_names[row_idx]
      text = text + ch_pad_right(row_name, 5)
    }

    let mut n = 0
    while n < col_order.length() {
      let col_idx = col_order[n]
      let value = matrix[row_idx][col_idx]

      let cell_char = if value < min_val + (max_val - min_val) * 0.2 {
        "█"
      } else if value < min_val + (max_val - min_val) * 0.4 {
        "▓"
      } else if value < min_val + (max_val - min_val) * 0.6 {
        "▒"
      } else if value < min_val + (max_val - min_val) * 0.8 {
        "░"
      } else {
        " "
      }

      text = text + cell_char + cell_char
      n = n + 1
    }
    text = text + "\n"
    m = m + 1
  }

  HeatmapResult::{ text, row_order, col_order, row_clusters, col_clusters }
}

///|
pub fn ch_create_annotation(
  annot : HeatmapAnnotation,
) -> HeatmapAnnotationResult {
  let mut row_annot = ""
  let mut col_annot = ""
  let mut legend = ""

  match annot.annot_type {
    RowAnnotation => {
      let mut i = 0
      while i < annot.labels.length() {
        row_annot = row_annot + annot.labels[i] + ": " + annot.colors[i] + "\n"
        i = i + 1
      }
    }
    ColumnAnnotation => {
      let mut i = 0
      while i < annot.labels.length() {
        col_annot = col_annot + annot.labels[i] + ": " + annot.colors[i] + "\n"
        i = i + 1
      }
    }
    CellAnnotation =>
      legend = "Cell annotations: " + annot.labels.length().to_string()
  }

  legend = "Legend:\n"
  let mut j = 0
  while j < annot.labels.length() && j < annot.colors.length() {
    legend = legend + annot.colors[j] + " " + annot.labels[j] + "\n"
    j = j + 1
  }

  HeatmapAnnotationResult::{ row_annot, col_annot, legend }
}

///|
pub fn ch_legend(cmap : ColorMap, title : String) -> String {
  let mut legend = title + "\n"
  legend = legend + "Min"
  let mut i = 0
  while i < cmap.colors.length() {
    legend = legend + " ██"
    i = i + 1
  }
  legend = legend + " Max\n"

  if cmap.breaks.length() > 0 {
    let mut j = 0
    while j < cmap.breaks.length() {
      legend = legend + cmap.breaks[j].to_string() + " "
      j = j + 1
    }
    legend = legend + "\n"
  }

  legend
}

///|
pub fn ch_merge_heatmaps(heatmaps : Array[HeatmapResult]) -> String {
  if heatmaps.length() == 0 {
    return ""
  }

  let mut merged = heatmaps[0].text
  let mut i = 1
  while i < heatmaps.length() {
    merged = merged + "\n---\n" + heatmaps[i].text
    i = i + 1
  }

  merged
}

///|
pub fn ch_row_split(
  data : HeatmapData,
  groups : Array[String],
) -> HeatmapResult {
  let group_map : Map[String, Array[Int]] = Map([], capacity=groups.length())
  let mut i = 0
  while i < groups.length() {
    let group = groups[i]
    match group_map.get(group) {
      Some(indices) => {
        let new_indices = indices + [i]
        group_map[group] = new_indices
      }
      None => group_map[group] = [i]
    }
    i = i + 1
  }

  let row_order : Array[Int] = Array::new()
  for group in group_map.keys().collect() {
    let indices = group_map[group]
    let mut j = 0
    while j < indices.length() {
      row_order.push(indices[j])
      j = j + 1
    }
  }

  let col_order : Array[Int] = Array::make(data.col_names.length(), 0)
  let mut k = 0
  while k < col_order.length() {
    col_order[k] = k
    k = k + 1
  }

  HeatmapResult::{
    text: "",
    row_order,
    col_order,
    row_clusters: [],
    col_clusters: [],
  }
}

///|
pub fn ch_col_split(
  data : HeatmapData,
  groups : Array[String],
) -> HeatmapResult {
  let group_map : Map[String, Array[Int]] = Map([], capacity=groups.length())
  let mut i = 0
  while i < groups.length() {
    let group = groups[i]
    match group_map.get(group) {
      Some(indices) => {
        let new_indices = indices + [i]
        group_map[group] = new_indices
      }
      None => group_map[group] = [i]
    }
    i = i + 1
  }

  let col_order : Array[Int] = Array::new()
  for group in group_map.keys().collect() {
    let indices = group_map[group]
    let mut j = 0
    while j < indices.length() {
      col_order.push(indices[j])
      j = j + 1
    }
  }

  let row_order : Array[Int] = Array::make(data.row_names.length(), 0)
  let mut k = 0
  while k < row_order.length() {
    row_order[k] = k
    k = k + 1
  }

  HeatmapResult::{
    text: "",
    row_order,
    col_order,
    row_clusters: [],
    col_clusters: [],
  }
}

///|
pub fn ch_add_annotation(
  heatmap : HeatmapResult,
  annot : HeatmapAnnotation,
) -> HeatmapResult {
  let annot_result = ch_create_annotation(annot)
  let new_text = heatmap.text + "\nAnnotations:\n" + annot_result.legend

  HeatmapResult::{
    text: new_text,
    row_order: heatmap.row_order,
    col_order: heatmap.col_order,
    row_clusters: heatmap.row_clusters,
    col_clusters: heatmap.col_clusters,
  }
}

///|
pub fn ch_annotate_rows(
  data : HeatmapData,
  annotations : Array[HeatmapAnnotation],
) -> String {
  let mut result = "Row Annotations:\n"
  let mut i = 0
  while i < annotations.length() {
    let annot = annotations[i]
    let annot_result = ch_create_annotation(annot)
    result = result + annot_result.row_annot
    i = i + 1
  }
  result
}

///|
pub fn ch_annotate_cols(
  data : HeatmapData,
  annotations : Array[HeatmapAnnotation],
) -> String {
  let mut result = "Column Annotations:\n"
  let mut i = 0
  while i < annotations.length() {
    let annot = annotations[i]
    let annot_result = ch_create_annotation(annot)
    result = result + annot_result.col_annot
    i = i + 1
  }
  result
}

///|
pub fn ch_summary(heatmap : HeatmapResult) -> String {
  let mut summary = "Heatmap Summary:\n"
  summary = summary + "Rows: " + heatmap.row_order.length().to_string() + "\n"
  summary = summary +
    "Columns: " +
    heatmap.col_order.length().to_string() +
    "\n"
  summary = summary +
    "Row clusters: " +
    heatmap.row_clusters.length().to_string() +
    "\n"
  summary = summary +
    "Column clusters: " +
    heatmap.col_clusters.length().to_string() +
    "\n"
  summary = summary +
    "Heatmap text length: " +
    heatmap.text.length().to_string() +
    " characters\n"
  summary
}

///|
pub fn ch_create_example_data() -> HeatmapData {
  let matrix = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [2.0, 3.0, 4.0, 5.0, 6.0],
    [3.0, 4.0, 5.0, 6.0, 7.0],
    [4.0, 5.0, 6.0, 7.0, 8.0],
    [5.0, 6.0, 7.0, 8.0, 9.0],
  ]
  let row_names = ["Gene1", "Gene2", "Gene3", "Gene4", "Gene5"]
  let col_names = ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"]

  HeatmapData::new(matrix, row_names, col_names)
}

///|
pub fn ch_create_example_annotation() -> HeatmapAnnotation {
  HeatmapAnnotation::new(
    RowAnnotation,
    ["GroupA", "GroupB", "GroupA", "GroupB", "GroupA"],
    ["#FF0000", "#00FF00", "#FF0000", "#00FF00", "#FF0000"],
  )
}

///|
pub fn ch_create_example_colormap() -> ColorMap {
  ColorMap::new(
    [
      "#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090",
      "#FDAE61", "#F46D43", "#D73027", "#A50026",
    ],
    [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0],
  )
}

///|
pub fn ch_create_example_options() -> HeatmapOptions {
  let options = HeatmapOptions::new()
  options
}
