///|
/// scanpy `tl.umap` — Uniform Manifold Approximation and Projection, 2D
/// embedding of a cell-connectivity graph. Port of umap-learn's
/// `simplicial_set_embedding` + `optimize_layout_euclidean` (spectral init +
/// stochastic gradient descent with negative sampling).

///|
fn scanpy_umap_lcg(seed : Int) -> Int {
  (seed * 1103515245 + 12345) & 0x7FFFFFFF
}

///|
fn scanpy_umap_clip(x : Double) -> Double {
  if x > 4.0 {
    4.0
  } else if x < -4.0 {
    -4.0
  } else {
    x
  }
}

///|
/// Uniform Manifold Approximation and Projection (2D). `connectivities` is the
/// symmetric fuzzy graph from `scanpy_neighbors`; returns an `n × 2` embedding.
pub fn scanpy_umap(
  connectivities : Array[Array[Double]],
  n_epochs : Int,
) -> Array[Array[Double]] {
  let n = connectivities.length()
  if n < 3 {
    // too few points to embed meaningfully: return a trivial layout
    let out : Array[Array[Double]] = Array(capacity=0)
    for i = 0; i < n; i = i + 1 {
      out.push([i.to_double(), 0.0])
    }
    return out
  }

  // collect positive-weight edges (upper triangle of the symmetric graph)
  let head : Array[Int] = Array(capacity=0)
  let tail : Array[Int] = Array(capacity=0)
  let weights : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let w = connectivities[i][j]
      if w > 0.0 {
        head.push(i)
        tail.push(j)
        weights.push(w)
      }
    }
  }
  let n_edges = head.length()

  // spectral init: eigenvectors 1 and 2 of the diffmap (skip the constant 0)
  let dm = scanpy_diffmap(connectivities, 3)
  let embedding : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    embedding.push([
      dm.eigen_basis[i][1],
      if dm.eigen_values.length() > 2 { dm.eigen_basis[i][2] } else { 0.0 },
    ])
  }
  if n_edges == 0 {
    return embedding
  }

  // UMAP curve params for min_dist = 0.1 (spread = 1)
  let a = 1.5769434603113077
  let b = 0.8950608779109733
  let gamma = 1.0
  let neg_rate = 5.0

  // epochs per sample = max_weight / weight
  let mut max_w = 0.0
  for w in weights {
    if w > max_w {
      max_w = w
    }
  }
  let eps = Array::make(n_edges, 0.0)
  for i = 0; i < n_edges; i = i + 1 {
    eps[i] = max_w / weights[i]
  }
  let next_sample = eps.copy()
  let eps_neg = Array::make(n_edges, 0.0)
  for i = 0; i < n_edges; i = i + 1 {
    eps_neg[i] = eps[i] / neg_rate
  }
  let next_neg = eps_neg.copy()

  let mut seed = 42
  let mut alpha = 1.0
  for e = 0; e < n_epochs; e = e + 1 {
    for i = 0; i < n_edges; i = i + 1 {
      if next_sample[i] <= e.to_double() {
        let j = head[i]
        let k = tail[i]
        let mut dist2 = 0.0
        for d = 0; d < 2; d = d + 1 {
          let diff = embedding[j][d] - embedding[k][d]
          dist2 = dist2 + diff * diff
        }
        let grad_coeff = if dist2 > 0.0 {
          -2.0 * a * b *
            @math.pow(dist2, b - 1.0) /
            (a * @math.pow(dist2, b) + 1.0)
        } else {
          0.0
        }
        for d = 0; d < 2; d = d + 1 {
          let gd = scanpy_umap_clip(grad_coeff * (embedding[j][d] - embedding[k][d]))
          embedding[j][d] = embedding[j][d] + gd * alpha
          embedding[k][d] = embedding[k][d] - gd * alpha
        }
        next_sample[i] = next_sample[i] + eps[i]

        let n_neg_f = (e.to_double() - next_neg[i]) / eps_neg[i]
        if n_neg_f > 0.0 {
          let n_neg = n_neg_f.to_int()
          for _p = 0; _p < n_neg; _p = _p + 1 {
            seed = scanpy_umap_lcg(seed)
            let kk = seed % n
            if kk == j {
              continue
            }
            let mut dn2 = 0.0
            for d = 0; d < 2; d = d + 1 {
              let diff = embedding[j][d] - embedding[kk][d]
              dn2 = dn2 + diff * diff
            }
            let gcoeff = if dn2 > 0.0 {
              2.0 * gamma * b /
                ((0.001 + dn2) * (a * @math.pow(dn2, b) + 1.0))
            } else {
              0.0
            }
            if gcoeff > 0.0 {
              for d = 0; d < 2; d = d + 1 {
                let gd = scanpy_umap_clip(gcoeff * (embedding[j][d] - embedding[kk][d]))
                embedding[j][d] = embedding[j][d] + gd * alpha
              }
            }
          }
        }
        next_neg[i] = next_neg[i] + n_neg_f.to_int().to_double() * eps_neg[i]
      }
    }
    alpha = 1.0 * (1.0 - e.to_double() / n_epochs.to_double())
  }

  embedding
}