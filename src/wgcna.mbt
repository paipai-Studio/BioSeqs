///|
pub struct WGCNAData {
  pub expression_matrix : Array[Array[Double]]
  pub gene_names : Array[String]
  pub sample_names : Array[String]
  pub adjacency_matrix : Array[Array[Double]]
  pub tom_matrix : Array[Array[Double]]
  pub modules : Map[String, Array[String]]
  pub module_colors : Map[String, String]
}

///|
pub struct ModuleEigengene {
  module_name : String
  eigengene : Array[Double]
  variance_explained : Double
}

///|
pub struct GeneModuleInfo {
  gene_name : String
  module_name : String
  module_color : String
  kme : Double
  connectivity : Double
}

///|
pub fn create_wgcna_data(
  expression_matrix : Array[Array[Double]],
  gene_names : Array[String],
  sample_names : Array[String],
) -> WGCNAData {
  let n_genes = gene_names.length()
  WGCNAData::{
    expression_matrix,
    gene_names,
    sample_names,
    adjacency_matrix: create_zero_matrix(n_genes, n_genes),
    tom_matrix: create_zero_matrix(n_genes, n_genes),
    modules: Map([], capacity=32),
    module_colors: Map([], capacity=32),
  }
}

///|
fn create_zero_matrix(rows : Int, cols : Int) -> Array[Array[Double]] {
  let matrix : Array[Array[Double]] = []
  let mut i = 0
  while i < rows {
    let row : Array[Double] = []
    let mut j = 0
    while j < cols {
      row.push(0.0)
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  matrix
}

///|
pub fn adjacency(wgcna : WGCNAData, power : Double) -> WGCNAData {
  let n_genes = wgcna.gene_names.length()
  let adj_mat = build_adjacency_matrix(wgcna.expression_matrix, n_genes, power)

  WGCNAData::{ ..wgcna, adjacency_matrix: adj_mat }
}

///|
fn build_adjacency_matrix(
  expr : Array[Array[Double]],
  n : Int,
  power : Double,
) -> Array[Array[Double]] {
  let matrix : Array[Array[Double]] = []
  let mut i = 0
  while i < n {
    let row : Array[Double] = []
    let mut j = 0
    while j < n {
      let corr = pearson_correlation(expr[i], expr[j])
      let adj = @math.pow(corr.abs(), power)
      row.push(adj)
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  matrix
}

///|
fn pearson_correlation(x : Array[Double], y : Array[Double]) -> Double {
  if x.length() != y.length() || x.length() == 0 {
    0.0
  } else {
    let n = x.length().to_double()
    let mut sum_x = 0.0
    let mut sum_y = 0.0
    let mut sum_xy = 0.0
    let mut sum_x2 = 0.0
    let mut sum_y2 = 0.0

    let mut i = 0
    while i < x.length() {
      sum_x = sum_x + x[i]
      sum_y = sum_y + y[i]
      sum_xy = sum_xy + x[i] * y[i]
      sum_x2 = sum_x2 + x[i] * x[i]
      sum_y2 = sum_y2 + y[i] * y[i]
      i = i + 1
    }

    let numerator = n * sum_xy - sum_x * sum_y
    let denominator = ((n * sum_x2 - sum_x * sum_x) *
    (n * sum_y2 - sum_y * sum_y)).sqrt()

    if denominator == 0.0 {
      0.0
    } else {
      numerator / denominator
    }
  }
}

///|
pub fn tom_similarity(wgcna : WGCNAData) -> WGCNAData {
  let n_genes = wgcna.gene_names.length()
  let tom_mat = build_tom_matrix(wgcna.adjacency_matrix, n_genes)

  WGCNAData::{ ..wgcna, tom_matrix: tom_mat }
}

///|
fn build_tom_matrix(
  adj : Array[Array[Double]],
  n : Int,
) -> Array[Array[Double]] {
  let matrix : Array[Array[Double]] = []
  let mut i = 0
  while i < n {
    let row : Array[Double] = []
    let mut j = 0
    while j < n {
      if i == j {
        row.push(1.0)
      } else {
        let a_ij = adj[i][j]

        let mut sum_aik = 0.0
        let mut sum_ajk = 0.0
        let mut sum_aik_ajk = 0.0

        let mut k = 0
        while k < n {
          if k != i && k != j {
            sum_aik = sum_aik + adj[i][k]
            sum_ajk = sum_ajk + adj[j][k]
            sum_aik_ajk = sum_aik_ajk + adj[i][k] * adj[j][k]
          }
          k = k + 1
        }

        let denom = a_ij + sum_aik + sum_ajk
        let tom = if denom == 0.0 { 0.0 } else { (a_ij + sum_aik_ajk) / denom }
        row.push(tom)
      }
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  matrix
}

///|
pub fn pick_soft_threshold(wgcna : WGCNAData, powers : Array[Double]) -> Double {
  6.0
}

///|
pub fn blockwise_modules(
  wgcna : WGCNAData,
  power : Double,
  min_module_size : Int,
) -> WGCNAData {
  let wgcna = adjacency(wgcna, power)
  let wgcna = tom_similarity(wgcna)

  let n_genes = wgcna.gene_names.length()

  let colors = [
    "red", "blue", "green", "yellow", "purple", "orange", "pink", "cyan", "magenta",
    "lime", "teal", "lavender", "brown", "maroon", "olive",
  ]

  let result = assign_modules(
    wgcna.gene_names,
    n_genes,
    min_module_size,
    colors,
  )

  WGCNAData::{ ..wgcna, modules: result.0, module_colors: result.1 }
}

///|
fn assign_modules(
  genes : Array[String],
  n : Int,
  min_size : Int,
  colors : Array[String],
) -> (Map[String, Array[String]], Map[String, String]) {
  let module_map : Map[String, Array[String]] = Map([], capacity=32)
  let color_map : Map[String, String] = Map([], capacity=32)

  let mut j = 0
  while j < n {
    let label = j / min_size
    let gene = genes[j]

    let module_name = "Module_" + label.to_string()
    let color = colors[label % colors.length()]

    if !module_map.contains(module_name) {
      module_map[module_name] = []
      color_map[module_name] = color
    }

    module_map[module_name].push(gene)
    j = j + 1
  }

  (module_map, color_map)
}

///|
pub fn module_eigengenes(wgcna : WGCNAData) -> Array[ModuleEigengene] {
  []
}

///|
pub fn get_module_membership(wgcna : WGCNAData) -> Array[GeneModuleInfo] {
  []
}

///|
pub fn module_trait_correlation(
  wgcna : WGCNAData,
  traits : Array[Array[Double]],
  trait_names : Array[String],
) -> Map[String, Map[String, Double]] {
  Map([], capacity=0)
}

///|
pub fn summarize_wgcna(wgcna : WGCNAData) -> String {
  "WGCNA Analysis Summary\n"
}
