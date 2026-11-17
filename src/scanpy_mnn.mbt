///|
/// Batch correction via Mutual Nearest Neighbors (MNN, Haghverdi et al. 2018 —
/// the `scanpy.external.pp.mnn_correct` counterpart). Batch 0 is the reference;
/// for every other batch we find mutually nearest cell pairs against the
/// reference, derive Gaussian-weighted mean difference vectors and add them to
/// every cell of the batch.

///|
/// MNN correction result: `corrected[b]` is the adjusted expression matrix of
/// batch `b` (reference batch unchanged), `n_mnn[b]` the number of mutual
/// nearest-neighbor pairs found for batch `b`.
pub struct ScanpyMnnResult {
  corrected : Array[Array[Array[Double]]]
  n_mnn : Array[Int32]
}

///|
/// Squared Euclidean distance between two vectors (truncated to the shorter one).
fn scanpy_mnn_d2(a : Array[Double], b : Array[Double]) -> Double {
  let m = if a.length() < b.length() { a.length() } else { b.length() }
  let mut s = 0.0
  for g = 0; g < m; g = g + 1 {
    let diff = a[g] - b[g]
    s = s + diff * diff
  }
  s
}

///|
/// Indices of the `kk` smallest entries of `d2` (the k nearest neighbors).
fn scanpy_mnn_nearest(d2 : Array[Double], kk : Int32) -> Set[Int32] {
  let order : Array[Int32] = Array(capacity=0)
  for i = 0; i < d2.length(); i = i + 1 {
    order.push(i)
  }
  order.sort_by(fn(a, b) {
    if d2[a] < d2[b] {
      -1
    } else if d2[a] > d2[b] {
      1
    } else {
      0
    }
  })
  let s : Set[Int32] = Set([], capacity=0)
  let t_max = if kk < order.length() { kk } else { order.length() }
  for t = 0; t < t_max; t = t + 1 {
    s.add(order[t])
  }
  s
}

///|
/// MNN batch correction over plain cell × gene matrices, one per batch.
/// `k` neighbors per direction define mutuality; `sigma` is the Gaussian
/// kernel width weighting difference vectors by batch-internal distance.
/// The reference batch (index 0) is returned unchanged.
pub fn scanpy_mnn_correct(
  mats : Array[Array[Array[Double]]],
  k~ : Int32 = 20,
  sigma~ : Double = 1.0,
) -> ScanpyMnnResult {
  let nb = mats.length()
  let corrected : Array[Array[Array[Double]]] = Array(capacity=nb)
  let n_mnn : Array[Int32] = Array(capacity=nb)
  if nb == 0 {
    return ScanpyMnnResult::{ corrected, n_mnn, }
  }
  corrected.push(mats[0].map(fn(row) { row.copy() }))
  n_mnn.push(0)
  let ref_mat = mats[0]
  let n_ref = ref_mat.length()
  let clamp_k = fn(n : Int32) -> Int32 {
    let v = if k < 1 { 1 } else { k }
    if v > n - 1 {
      n - 1
    } else {
      v
    }
  }
  let k_ref = clamp_k(n_ref)
  for b = 1; b < nb; b = b + 1 {
    let bat_mat = mats[b]
    let n_bat = bat_mat.length()
    let k_bat = clamp_k(n_bat)
    // cross distances ref[i] <-> batch[j]
    let d2 : Array[Array[Double]] = Array(capacity=n_ref)
    for i = 0; i < n_ref; i = i + 1 {
      let row : Array[Double] = Array(capacity=n_bat)
      for j = 0; j < n_bat; j = j + 1 {
        row.push(scanpy_mnn_d2(ref_mat[i], bat_mat[j]))
      }
      d2.push(row)
    }
    // MNN pairs: (i, j) mutually among the k nearest of each other
    let is_mnn : Array[Array[Bool]] = Array(capacity=n_ref)
    for _i in 0..<n_ref {
      is_mnn.push(Array::make(n_bat, false))
    }
    // precompute nearest-neighbor candidate sets per side
    let ref_nearest : Array[Set[Int32]] = Array(capacity=n_bat)
    for j = 0; j < n_bat; j = j + 1 {
      let col : Array[Double] = Array(capacity=n_ref)
      for i = 0; i < n_ref; i = i + 1 {
        col.push(d2[i][j])
      }
      ref_nearest.push(scanpy_mnn_nearest(col, k_ref))
    }
    let bat_nearest : Array[Set[Int32]] = Array(capacity=n_ref)
    for i = 0; i < n_ref; i = i + 1 {
      bat_nearest.push(scanpy_mnn_nearest(d2[i], k_bat))
    }
    let mut pairs = 0
    let mnn_i : Array[Int32] = Array(capacity=0)
    let mnn_j : Array[Int32] = Array(capacity=0)
    for i = 0; i < n_ref; i = i + 1 {
      for j in bat_nearest[i] {
        if ref_nearest[j].contains(i) {
          is_mnn[i][j] = true
          mnn_i.push(i)
          mnn_j.push(j)
          pairs = pairs + 1
        }
      }
    }
    n_mnn.push(pairs)
    // correction: every batch cell gets the Gaussian-weighted mean difference
    // over all MNN pairs, weighted by its batch-internal distance to the
    // paired batch cell
    let corrected_b : Array[Array[Double]] = Array(capacity=n_bat)
    for j = 0; j < n_bat; j = j + 1 {
      let m = bat_mat[j].length()
      let out = bat_mat[j].copy()
      let mut w_sum = 0.0
      let corr = Array::make(m, 0.0)
      for p = 0; p < mnn_i.length(); p = p + 1 {
        let i = mnn_i[p]
        let jj = mnn_j[p]
        let w = @math.pow(
          2.718281828459045,
          -scanpy_mnn_d2(bat_mat[j], bat_mat[jj]) / (2.0 * sigma * sigma),
        )
        w_sum = w_sum + w
        for g = 0
            g < m && g < ref_mat[i].length() && g < bat_mat[jj].length()
            g = g + 1 {
          corr[g] = corr[g] + w * (ref_mat[i][g] - bat_mat[jj][g])
        }
      }
      if w_sum > 0.0 {
        for g = 0; g < m; g = g + 1 {
          out[g] = out[g] + corr[g] / w_sum
        }
      }
      corrected_b.push(out)
    }
    corrected.push(corrected_b)
  }
  ScanpyMnnResult::{ corrected, n_mnn, }
}

///|
/// AnnData convenience: correct every batch in `adatas` against batch 0 in
/// place (replaces each `adata.x`). All batches must share the same variables.
pub fn scanpy_mnn_correct_adatas(
  adatas : Array[AnnData],
  k~ : Int32 = 20,
  sigma~ : Double = 1.0,
) -> ScanpyMnnResult {
  let mats : Array[Array[Array[Double]]] = Array(capacity=adatas.length())
  for ad in adatas {
    mats.push(ad.x_matrix())
  }
  let res = scanpy_mnn_correct(mats, k~, sigma~)
  for b = 0; b < adatas.length() && b < res.corrected.length(); b = b + 1 {
    adatas[b].x = res.corrected[b]
  }
  res
}
