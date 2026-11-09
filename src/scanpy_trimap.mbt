///|
/// scanpy `tl.trimap` — TriMAP: large-scale dimensionality reduction using
/// triplets (Amid & Warmuth 2019).  Each point is anchored by sampled
/// triplets (point, inlier, outlier): the embedding must place `point` closer
/// to its inlier (a kNN neighbour, weight `1/n_inliers`) than to randomly
/// drawn outliers (weight `1/n_outliers`).  The triplet loss
/// `log(1 + (d_ij^2 + gamma) / (d_ik^2 + gamma))` is minimized by normalized
/// SGD with a linearly decaying learning rate.

///|
/// TriMAP result: coordinates (`n x n_components`) and the triplet sampling
/// configuration actually used.
pub struct ScanpyTriMap {
  embedding : Array[Array[Double]]
  n_inliers : Int32
  n_outliers : Int32
  n_random : Int32
  n_iter : Int32
}

///|
/// TriMAP embedding of an `n x m` matrix.  `n_random` is accepted for API
/// compatibility with the reference implementation (random outliers already
/// realize the global-structure sampling) and reported back unchanged.
pub fn scanpy_trimap(
  x : Array[Array[Double]],
  n_components~ : Int32 = 2,
  n_inliers~ : Int32 = 10,
  n_outliers~ : Int32 = 5,
  n_random~ : Int32 = 5,
  n_iter~ : Int32 = 400,
  weight_adj~ : Double = 500.0,
  lr~ : Double = 0.2,
  seed~ : Int32 = 0,
) -> ScanpyTriMap {
  let n = x.length()
  if n == 0 {
    return ScanpyTriMap::{
      embedding: [],
      n_inliers: 0,
      n_outliers: 0,
      n_random: 0,
      n_iter: 0,
    }
  }
  let mut rng = if seed == 0 { 3 } else { seed }
  // kNN (inlier pool)
  let inl = if n_inliers < 1 { 1 } else { n_inliers }
  let k = if inl + 1 < n { inl + 1 } else { n }
  let (knn_idx, _knn_dists) = scanpy_knn_distances(x, k)
  // gamma: small loss-softening constant relative to the embedding scale
  // (TriMAP default 1e-2).  Scaling it with the data-space max distance would
  // swamp the triplet gradient, since embedding distances start near 1.
  let gamma = if weight_adj > 0.0 {
    1.0e-2 * weight_adj / 500.0
  } else {
    1.0e-2
  }
  // random init in [-1, 1]^n_components
  let pos : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for _c = 0; _c < n_components; _c = _c + 1 {
      rng = scanpy_trimap_next(rng)
      row.push(rng.to_double() / 1073741823.5 - 1.0)
    }
    pos.push(row)
  }
  let n_out = if n_outliers < 1 { 1 } else { n_outliers }
  // degenerate inputs: not enough points to form triplets
  if n < 3 || k < 2 {
    return ScanpyTriMap::{
      embedding: pos,
      n_inliers: inl,
      n_outliers: n_out,
      n_random,
      n_iter,
    }
  }
  for iter = 0; iter < n_iter; iter = iter + 1 {
    let cur_lr = lr * (1.0 - iter.to_double() / n_iter.to_double()) + 0.01 * lr
    let inlier_slot = iter % (k - 1) + 1 // rotate through the kNN list
    for i = 0; i < n; i = i + 1 {
      let j = knn_idx[i][inlier_slot]
      if j == i {
        continue
      }
      for _o = 0; _o < n_out; _o = _o + 1 {
        rng = scanpy_trimap_next(rng)
        let kk = (rng % n).abs()
        if kk == i || kk == j {
          continue
        }
        // embedding-space distances
        let mut dij2 = 0.0
        let mut dik2 = 0.0
        for c = 0; c < n_components; c = c + 1 {
          let diff_ij = pos[i][c] - pos[j][c]
          let diff_ik = pos[i][c] - pos[kk][c]
          dij2 = dij2 + diff_ij * diff_ij
          dik2 = dik2 + diff_ik * diff_ik
        }
        let a = dij2 + gamma
        let b = dik2 + gamma
        if a / b < 0.01 {
          continue // triplet already satisfied
        }
        // gradients of log(1 + a/b); step normalized per triplet
        let denom = b * (a + b)
        let step = cur_lr
        for c = 0; c < n_components; c = c + 1 {
          let vij = pos[i][c] - pos[j][c]
          let vik = pos[i][c] - pos[kk][c]
          let grad_i = 2.0 * (b * vij - a * vik) / denom
          let grad_j = -2.0 * vij / (a + b)
          let grad_k = 2.0 * a * vik / denom
          pos[i][c] = pos[i][c] - step * grad_i / (grad_i.abs() + 1.0)
          pos[j][c] = pos[j][c] - step * grad_j / (grad_j.abs() + 1.0)
          pos[kk][c] = pos[kk][c] - step * grad_k / (grad_k.abs() + 1.0)
        }
      }
    }
  }
  ScanpyTriMap::{
    embedding: pos,
    n_inliers: inl,
    n_outliers: n_out,
    n_random,
    n_iter,
  }
}

///|
/// LCG step for the TriMAP sampler.
fn scanpy_trimap_next(state : Int32) -> Int32 {
  (state * 1103515245 + 12345) & 0x7FFFFFFF
}
