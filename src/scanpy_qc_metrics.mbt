/// Scanpy second-round: QC quality control metrics
/// Inspired by scanpy.pp.calculate_qc_metrics

///|
/// QC metrics result for a single cell or gene.
pub struct QcMetrics {
  total_counts : Array[Double] // total counts per cell
  n_genes_by_counts : Array[Int32] // number of genes with counts > 0 per cell
  total_counts_mt : Array[Double] // total mitochondrial counts per cell
  pct_counts_mt : Array[Double] // percent mitochondrial counts per cell
  total_counts_ribo : Array[Double] // total ribosomal counts per cell
  pct_counts_ribo : Array[Double] // percent ribosomal counts per cell
  total_counts_hemo : Array[Double] // total hemoglobin counts per cell
  pct_counts_hemo : Array[Double] // percent hemoglobin counts per cell
}

///|
/// Calculate standard QC metrics for each cell.
/// Computes total counts, number of detected genes, and
/// percentages of mitochondrial / ribosomal / hemoglobin counts.
pub fn scanpy_calculate_qc_metrics(
  x : Array[Array[Double]],
  var_names : Array[String],
  mt_prefix~ : String = "MT-",
  ribo_prefix~ : String = "RPS|RPL",
  hemo_prefix~ : String = "HBA|HBB|HBD",
) -> QcMetrics {
  let n_cells = x.length()
  let n_genes = if n_cells > 0 { x[0].length() } else { 0 }
  let total_counts = Array(capacity=n_cells)
  let n_genes_by_counts = Array(capacity=n_cells)
  let total_counts_mt = Array(capacity=n_cells)
  let pct_counts_mt = Array(capacity=n_cells)
  let total_counts_ribo = Array(capacity=n_cells)
  let pct_counts_ribo = Array(capacity=n_cells)
  let total_counts_hemo = Array(capacity=n_cells)
  let pct_counts_hemo = Array(capacity=n_cells)
  // Precompute masks
  let mt_mask = Array(capacity=n_genes)
  let ribo_mask = Array(capacity=n_genes)
  let hemo_mask = Array(capacity=n_genes)
  for j in 0..<n_genes {
    let name = var_names[j].to_upper()
    mt_mask.push(name.has_prefix(mt_prefix.to_upper()))
    ribo_mask.push(name.has_prefix("RPS") || name.has_prefix("RPL"))
    hemo_mask.push(
      name.has_prefix("HBA") || name.has_prefix("HBB") || name.has_prefix("HBD"),
    )
    ignore(ribo_prefix)
    ignore(hemo_prefix)
  }
  for i in 0..<n_cells {
    let mut tc = 0.0
    let mut ng = 0
    let mut tc_mt = 0.0
    let mut tc_ribo = 0.0
    let mut tc_hemo = 0.0
    for j in 0..<n_genes {
      let v = x[i][j]
      tc = tc + v
      if v > 0.0 {
        ng = ng + 1
      }
      if mt_mask[j] {
        tc_mt = tc_mt + v
      }
      if ribo_mask[j] {
        tc_ribo = tc_ribo + v
      }
      if hemo_mask[j] {
        tc_hemo = tc_hemo + v
      }
    }
    total_counts.push(tc)
    n_genes_by_counts.push(ng)
    total_counts_mt.push(tc_mt)
    total_counts_ribo.push(tc_ribo)
    total_counts_hemo.push(tc_hemo)
    if tc > 0.0 {
      pct_counts_mt.push(tc_mt / tc * 100.0)
      pct_counts_ribo.push(tc_ribo / tc * 100.0)
      pct_counts_hemo.push(tc_hemo / tc * 100.0)
    } else {
      pct_counts_mt.push(0.0)
      pct_counts_ribo.push(0.0)
      pct_counts_hemo.push(0.0)
    }
  }
  QcMetrics::{
    total_counts,
    n_genes_by_counts,
    total_counts_mt,
    pct_counts_mt,
    total_counts_ribo,
    pct_counts_ribo,
    total_counts_hemo,
    pct_counts_hemo,
  }
}

///|
/// Compute per-gene QC metrics: total counts, number of cells expressing, mean.
pub fn scanpy_calculate_qc_metrics_genes(
  x : Array[Array[Double]],
) -> Array[(Double, Int32, Double)] {
  let n_cells = x.length()
  let n_genes = if n_cells > 0 { x[0].length() } else { 0 }
  let result = Array(capacity=n_genes)
  for j in 0..<n_genes {
    let mut total = 0.0
    let mut n_cells_expr = 0
    for i in 0..<n_cells {
      let v = x[i][j]
      total = total + v
      if v > 0.0 {
        n_cells_expr = n_cells_expr + 1
      }
    }
    let mean = if n_cells > 0 { total / n_cells.to_double() } else { 0.0 }
    result.push((total, n_cells_expr, mean))
  }
  result
}

///|
/// Recommend filtering thresholds based on QC metric percentiles.
/// Returns (min_total_counts, max_total_counts, min_n_genes, max_pct_mt).
pub fn scanpy_qc_filter_threshold(
  qc : QcMetrics,
  low_quantile~ : Double = 0.01,
  high_quantile~ : Double = 0.99,
) -> (Double, Double, Int32, Double) {
  let n = qc.total_counts.length()
  let sorted_tc = Array(capacity=n)
  let sorted_ng = Array(capacity=n)
  let sorted_mt = Array(capacity=n)
  for i in 0..<n {
    sorted_tc.push(qc.total_counts[i])
    sorted_ng.push(qc.n_genes_by_counts[i].to_double())
    sorted_mt.push(qc.pct_counts_mt[i])
  }
  sorted_tc.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  sorted_ng.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  sorted_mt.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  let lo_idx = (low_quantile * n.to_double()).floor().to_int()
  let hi_idx = (high_quantile * n.to_double()).floor().to_int()
  let lo = if lo_idx >= n { n - 1 } else { lo_idx }
  let hi = if hi_idx >= n { n - 1 } else { hi_idx }
  let lo2 = if lo < 0 { 0 } else { lo }
  let hi2 = if hi < 0 { 0 } else { hi }
  (
    sorted_tc[lo2],
    sorted_tc[hi2],
    sorted_ng[lo2].floor().to_int(),
    sorted_mt[hi2],
  )
}
