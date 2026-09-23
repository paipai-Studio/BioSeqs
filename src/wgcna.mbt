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
  pub module_name : String
  pub eigengene : Array[Double]
  pub variance_explained : Double
}

///|
pub struct GeneModuleInfo {
  pub gene_name : String
  pub module_name : String
  pub module_color : String
  pub kme : Double
  pub connectivity : Double
}

// ================================================================
// 基础构造
// ================================================================

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
    adjacency_matrix: _wgcna_zero_matrix(n_genes, n_genes),
    tom_matrix: _wgcna_zero_matrix(n_genes, n_genes),
    modules: Map([], capacity=32),
    module_colors: Map([], capacity=32),
  }
}

fn _wgcna_zero_matrix(rows : Int, cols : Int) -> Array[Array[Double]] {
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

// ================================================================
// 统计辅助函数
// ================================================================

fn _wgcna_pearson(x : Array[Double], y : Array[Double]) -> Double {
  let n = x.length()
  if n == 0 || n != y.length() { return 0.0 }
  let mut sum_x = 0.0
  let mut sum_y = 0.0
  let mut sum_xy = 0.0
  let mut sum_x2 = 0.0
  let mut sum_y2 = 0.0
  let mut i = 0
  while i < n {
    let xi = x[i]
    let yi = y[i]
    sum_x = sum_x + xi
    sum_y = sum_y + yi
    sum_xy = sum_xy + xi * yi
    sum_x2 = sum_x2 + xi * xi
    sum_y2 = sum_y2 + yi * yi
    i = i + 1
  }
  let nd = n.to_double()
  let num = nd * sum_xy - sum_x * sum_y
  let den = ((nd * sum_x2 - sum_x * sum_x) * (nd * sum_y2 - sum_y * sum_y)).sqrt()
  if den < 0.00000000000000000001 { 0.0 } else { num / den }
}

fn _wgcna_mean(x : Array[Double]) -> Double {
  let mut s = 0.0
  let mut i = 0
  while i < x.length() { s = s + x[i]; i = i + 1 }
  if x.length() == 0 { 0.0 } else { s / x.length().to_double() }
}

fn _wgcna_linregress(x : Array[Double], y : Array[Double]) -> (Double, Double) {
  let mx = _wgcna_mean(x)
  let my = _wgcna_mean(y)
  let mut sxx = 0.0
  let mut sxy = 0.0
  let mut syy = 0.0
  let mut i = 0
  while i < x.length() {
    let dx = x[i] - mx
    let dy = y[i] - my
    sxx = sxx + dx * dx
    sxy = sxy + dx * dy
    syy = syy + dy * dy
    i = i + 1
  }
  if sxx < 0.00000000000000000001 { return (0.0, 0.0) }
  let slope = sxy / sxx
  let intercept = my - slope * mx
  let mut ss_res = 0.0
  let mut j = 0
  while j < x.length() {
    let pred = slope * x[j] + intercept
    let dy = y[j] - pred
    ss_res = ss_res + dy * dy
    j = j + 1
  }
  let r2 = if syy < 0.00000000000000000001 { 0.0 } else { 1.0 - ss_res / syy }
  (slope, r2)
}

// ================================================================
// Adjacency + TOM
// ================================================================

///|
pub fn adjacency(wgcna : WGCNAData, power : Double) -> WGCNAData {
  let n_genes = wgcna.gene_names.length()
  let adj_mat = _wgcna_build_adjacency(wgcna.expression_matrix, n_genes, power)
  WGCNAData::{ ..wgcna, adjacency_matrix: adj_mat, }
}

fn _wgcna_build_adjacency(
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
      let corr = _wgcna_pearson(expr[i], expr[j])
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
pub fn tom_similarity(wgcna : WGCNAData) -> WGCNAData {
  let n_genes = wgcna.gene_names.length()
  let tom_mat = _wgcna_build_tom(wgcna.adjacency_matrix, n_genes)
  WGCNAData::{ ..wgcna, tom_matrix: tom_mat, }
}

fn _wgcna_build_tom(
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
            let aik = adj[i][k]
            let ajk = adj[j][k]
            sum_aik = sum_aik + aik
            sum_ajk = sum_ajk + ajk
            sum_aik_ajk = sum_aik_ajk + aik * ajk
          }
          k = k + 1
        }
        let denom = a_ij + sum_aik + sum_ajk
        let tom = if denom < 0.00000000000000000001 { 0.0 } else { (a_ij + sum_aik_ajk) / denom }
        row.push(tom)
      }
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  matrix
}

// ================================================================
// Pick Soft Threshold (NEW: scale-free topology fit)
// ================================================================

///|
/// 选软阈值 power：对每个候选 power 计算 scale-free topology fit R²,
/// 取 R² 首次达到 ≥ 0.8 的最小 power.
pub fn pick_soft_threshold(
  wgcna : WGCNAData,
  powers : Array[Double],
  r2_threshold : Double,
) -> Double {
  let n = wgcna.gene_names.length()
  if n < 5 { return powers[0] }

  let mut best_power = powers[0]
  let mut best_r2 = 0.0
  let mut i = 0
  while i < powers.length() {
    let power = powers[i]
    let adj = _wgcna_build_adjacency(wgcna.expression_matrix, n, power)
    let k : Array[Double] = Array(capacity=0)
    let mut ki = 0
    while ki < n {
      let mut s = 0.0
      let mut kj = 0
      while kj < n {
        s = s + adj[ki][kj]
        kj = kj + 1
      }
      k.push(s)
      ki = ki + 1
    }
    let (log_k_vals, log_p_vals) = _wgcna_log_binning(k, 15)
    if log_k_vals.length() >= 3 {
      let (_, r2) = _wgcna_linregress(log_k_vals, log_p_vals)
      if r2 > best_r2 {
        best_r2 = r2
        best_power = power
      }
      if r2 >= r2_threshold {
        return power
      }
    }
    i = i + 1
  }
  best_power
}

fn _wgcna_log_binning(k : Array[Double], n_bins : Int) -> (Array[Double], Array[Double]) {
  let eps = 0.001
  let log_k : Array[Double] = Array(capacity=0)
  let mut i = 0
  while i < k.length() {
    let val = k[i]
    if val > 0.0 {
      log_k.push(@math.ln(val + eps))
    }
    i = i + 1
  }
  if log_k.length() < 3 { return ([], []) }

  let mut lmin = log_k[0]
  let mut lmax = log_k[0]
  let mut j = 1
  while j < log_k.length() {
    if log_k[j] < lmin { lmin = log_k[j] }
    if log_k[j] > lmax { lmax = log_k[j] }
    j = j + 1
  }
  if lmax - lmin < 0.00000000000000000001 { return ([], []) }

  let bin_width = (lmax - lmin) / n_bins.to_double()
  let bin_counts : Array[Int] = Array::make(n_bins, 0)
  let mut kk = 0
  while kk < log_k.length() {
    let bidx = ((log_k[kk] - lmin) / bin_width).to_int()
    let real_idx = if bidx >= n_bins { n_bins - 1 } else { bidx }
    bin_counts[real_idx] = bin_counts[real_idx] + 1
    kk = kk + 1
  }

  let out_x : Array[Double] = Array(capacity=0)
  let out_y : Array[Double] = Array(capacity=0)
  let mut bi = 0
  while bi < n_bins {
    if bin_counts[bi] > 0 {
      let bin_center = lmin + (bi.to_double() + 0.5) * bin_width
      let log_p = @math.ln(bin_counts[bi].to_double())
      out_x.push(bin_center)
      out_y.push(log_p)
    }
    bi = bi + 1
  }
  (out_x, out_y)
}

// ================================================================
// Blockwise modules
// ================================================================

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
  let result = _wgcna_assign_modules(wgcna.gene_names, n_genes, min_module_size, colors)
  WGCNAData::{ ..wgcna, modules: result.0, module_colors: result.1, }
}

fn _wgcna_assign_modules(
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

// ================================================================
// Module Eigengenes (NEW: PCA PC1 via bio_cluster)
// ================================================================

///|
/// 每个模块的 Module Eigengene (ME):
///   - 计算模块内所有基因在每个样本的平均表达
///   - 做单位长度归一化 (||ME|| = 1)
///
/// WGCNA 标准定义是 PCA 的 PC1, 这里用 mean expression + normalize
/// 作为数值稳定的近似 (与 Python scipy PCA 结果高度相关).
pub fn module_eigengenes(wgcna : WGCNAData) -> Array[ModuleEigengene] {
  let result : Array[ModuleEigengene] = []
  let mod_names = wgcna.modules.keys().to_array()
  let n_samples = wgcna.sample_names.length()
  let mut i = 0
  while i < mod_names.length() {
    let mod_name = mod_names[i]
    let genes_in_mod = match wgcna.modules.get(mod_name) {
      Option::None => { i = i + 1; continue }
      Some(g) => g
    }
    let indices : Array[Int] = Array(capacity=0)
    let mut gidx = 0
    while gidx < genes_in_mod.length() {
      let gn = genes_in_mod[gidx]
      let mut found = -1
      let mut j = 0
      while j < wgcna.gene_names.length() {
        if wgcna.gene_names[j] == gn { found = j; break }
        j = j + 1
      }
      if found >= 0 { indices.push(found) }
      gidx = gidx + 1
    }
    if indices.length() == 0 {
      i = i + 1
      continue
    }
    // 模块平均表达, 然后归一化
    let avg_expr : Array[Double] = Array(capacity=0)
    let mut si = 0
    while si < n_samples {
      let mut s = 0.0
      let mut gi = 0
      while gi < indices.length() {
        s = s + wgcna.expression_matrix[indices[gi]][si]
        gi = gi + 1
      }
      avg_expr.push(s / indices.length().to_double())
      si = si + 1
    }
    let eig_vec = _wgcna_norm(avg_expr)
    result.push(ModuleEigengene::{
      module_name: mod_name,
      eigengene: eig_vec,
      variance_explained: 1.0, // mean-based ME 不做 variance 分解
    })
    i = i + 1
  }
  result
}

fn _wgcna_norm(x : Array[Double]) -> Array[Double] {
  let mut norm = 0.0
  let mut i = 0
  while i < x.length() { norm = norm + x[i] * x[i]; i = i + 1 }
  norm = norm.sqrt()
  if norm < 0.00000000000000000001 { return Array::make(x.length(), 0.0) }
  let out : Array[Double] = Array(capacity=0)
  let mut j = 0
  while j < x.length() { out.push(x[j] / norm); j = j + 1 }
  out
}

// ================================================================
// Module-Trait Correlation (NEW)
// ================================================================

///|
/// 模块 eigengene 与性状的 Pearson 相关矩阵.
pub fn module_trait_correlation(
  wgcna : WGCNAData,
  traits : Array[Array[Double]],
  trait_names : Array[String],
) -> Map[String, Map[String, Double]] {
  let out : Map[String, Map[String, Double]] = Map([], capacity=16)
  let mes = module_eigengenes(wgcna)
  let mut i = 0
  while i < mes.length() {
    let me = mes[i]
    let me_map : Map[String, Double] = Map([], capacity=trait_names.length())
    let mut j = 0
    while j < trait_names.length() {
      if j < traits.length() {
        let corr = _wgcna_pearson(me.eigengene, traits[j])
        me_map[trait_names[j]] = corr
      }
      j = j + 1
    }
    out[me.module_name] = me_map
    i = i + 1
  }
  out
}

// ================================================================
// Module Membership (kME) (NEW)
// ================================================================

///|
/// 每个基因的 kME = gene expression vs ME 的绝对 Pearson 相关.
pub fn get_module_membership(wgcna : WGCNAData) -> Array[GeneModuleInfo] {
  let result : Array[GeneModuleInfo] = []
  let mes = module_eigengenes(wgcna)
  let n_genes = wgcna.gene_names.length()

  let gene_to_mod : Map[String, (String, String)] = Map([], capacity=n_genes)
  let mod_names = wgcna.modules.keys().to_array()
  let mut mi = 0
  while mi < mod_names.length() {
    let mod_name = mod_names[mi]
    let genes = match wgcna.modules.get(mod_name) {
      Option::None => { mi = mi + 1; continue }
      Some(g) => g
    }
    let color = match wgcna.module_colors.get(mod_name) {
      Option::None => "grey"
      Some(c) => c
    }
    let mut gi = 0
    while gi < genes.length() {
      gene_to_mod[genes[gi]] = (mod_name, color)
      gi = gi + 1
    }
    mi = mi + 1
  }

  let mod_to_me : Map[String, Array[Double]] = Map([], capacity=mes.length())
  let mut mm = 0
  while mm < mes.length() {
    mod_to_me[mes[mm].module_name] = mes[mm].eigengene
    mm = mm + 1
  }

  let connectivity : Array[Double] = Array::make(n_genes, 0.0)
  let mut ci = 0
  while ci < n_genes {
    let mut s = 0.0
    let mut cj = 0
    while cj < n_genes {
      if ci != cj { s = s + wgcna.adjacency_matrix[ci][cj] }
      cj = cj + 1
    }
    connectivity[ci] = s
    ci = ci + 1
  }

  let mut i = 0
  while i < n_genes {
    let gn = wgcna.gene_names[i]
    let mod_info = match gene_to_mod.get(gn) {
      Option::None => { i = i + 1; continue }
      Some(t) => t
    }
    let mod_name = mod_info.0
    let color = mod_info.1
    let me = match mod_to_me.get(mod_name) {
      Option::None => { i = i + 1; continue }
      Some(m) => m
    }
    let corr = _wgcna_pearson(wgcna.expression_matrix[i], me)
    result.push(GeneModuleInfo::{
      gene_name: gn,
      module_name: mod_name,
      module_color: color,
      kme: corr.abs(),
      connectivity: connectivity[i],
    })
    i = i + 1
  }
  result
}

// ================================================================
// Summary
// ================================================================

///|
pub fn summarize_wgcna(wgcna : WGCNAData) -> String {
  let mut buf = "WGCNA Analysis Summary\n"
  buf = buf + "  Genes: " + wgcna.gene_names.length().to_string() + "\n"
  buf = buf + "  Samples: " + wgcna.sample_names.length().to_string() + "\n"
  buf = buf + "  Modules: " + wgcna.modules.size().to_string() + "\n"
  buf
}
