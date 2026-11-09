///|
/// scanpy `tl.phate` — PHATE (Potential of Heat-diffusion for Affinity-based
/// Trajectory Embedding, Moon et al. 2019).  Pipeline: kNN graph with a
/// local-scale Gaussian kernel, row-normalized diffusion operator, diffusion
/// time chosen by matching the von Neumann entropy decay (80% of VNE(t=1)),
/// diffusion distances computed from the symmetrized eigenbasis
/// `coord_k(i) = lambda_k^t * phi_k(i) / sqrt(deg(i))`, and a SMACOF
/// metric multidimensional scaling onto `n_components`.

///|
/// PHATE embedding result: coordinates (`n x n_components`), the diffusion
/// time used (auto-selected when `t=0`) and the von Neumann entropy of the
/// diffusion operator at that time.
pub struct ScanpyPhate {
  embedding : Array[Array[Double]]
  t : Int32
  vne : Double
}

///|
/// Small deterministic LCG (same parameters as `scanpy_draw_graph`).
fn scanpy_phate_rand(state : Int32) -> Int32 {
  (state * 1103515245 + 12345) & 0x7FFFFFFF
}

///|
/// Row-normalized diffusion operator of the local-scale Gaussian kNN kernel
/// as a dense matrix (diagonal self-similarity 1, edge weights
/// `exp(-d^2 / (sigma_i * sigma_j))` averaged over the two directions).
fn scanpy_phate_diff_op(
  knn_idx : Array[Array[Int32]],
  knn_dists : Array[Array[Double]],
  k : Int32,
  n : Int32,
) -> Array[Array[Double]] {
  let sigma = Array::make(n, 1.0)
  for i = 0; i < n; i = i + 1 {
    let sk = knn_dists[i][k - 1]
    sigma[i] = if sk > 1.0e-12 { sk } else { 1.0 }
  }
  // accumulate directed kernel values, then average both directions
  let sums : Map[Int32, Double] = Map([], capacity=n * k)
  let counts : Map[Int32, Int32] = Map([], capacity=n * k)
  for i = 0; i < n; i = i + 1 {
    for j = 1; j < k; j = j + 1 {
      let nb = knn_idx[i][j]
      if nb == i {
        continue
      }
      let d = knn_dists[i][j]
      let denom = sigma[i] * sigma[nb]
      let w = if denom > 1.0e-300 { @math.exp(-(d * d / denom)) } else { 0.0 }
      let a = if i < nb { i } else { nb }
      let b = if i < nb { nb } else { i }
      let key = a * n + b
      match sums.get(key) {
        None => sums.set(key, w)
        Some(x) => sums.set(key, x + w)
      }
      match counts.get(key) {
        None => counts.set(key, 1)
        Some(c) => counts.set(key, c + 1)
      }
    }
  }
  let a_mat : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    row[i] = 1.0
    a_mat.push(row)
  }
  for kv in sums {
    let key = kv.0
    let w = kv.1 / counts.get(key).unwrap().to_double()
    let a = key / n
    let b = key % n
    a_mat[a][b] = w
    a_mat[b][a] = w
  }
  // row-normalize
  for i = 0; i < n; i = i + 1 {
    let mut s = 0.0
    for j = 0; j < n; j = j + 1 {
      s = s + a_mat[i][j]
    }
    if s > 0.0 {
      for j = 0; j < n; j = j + 1 {
        a_mat[i][j] = a_mat[i][j] / s
      }
    }
  }
  a_mat
}

///|
/// Von Neumann entropy of the t-step diffusion operator restricted to the
/// top-k eigenvalues: `-sum(lambda^t * log(lambda^t))`.
fn scanpy_phate_vne(evals : Array[Double], t : Int32) -> Double {
  let mut s = 0.0
  let tf = t.to_double()
  for l in evals {
    let lam = if l > 0.0 { l } else { 0.0 }
    let lt = @math.pow(lam, tf)
    if lt > 1.0e-300 {
      s = s - lt * @math.ln(lt)
    }
  }
  s
}

///|
/// SMACOF metric multidimensional scaling: stress-minimizing Guttman
/// transforms on the (precomputed) target distance matrix.  Returns the
/// final configuration.
fn scanpy_phate_smacof(
  target : Array[Array[Double]],
  n_components : Int32,
  n_iter : Int32,
  seed : Int32,
) -> Array[Array[Double]] {
  let n = target.length()
  if n == 0 {
    return []
  }
  let mut rng = if seed == 0 { 7 } else { seed }
  let pos : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for _c = 0; _c < n_components; _c = _c + 1 {
      rng = scanpy_phate_rand(rng)
      row.push(rng.to_double() / 1073741823.5 - 1.0)
    }
    pos.push(row)
  }
  for _iter = 0; _iter < n_iter; _iter = _iter + 1 {
    let new_pos : Array[Array[Double]] = Array(capacity=n)
    for i = 0; i < n; i = i + 1 {
      let row : Array[Double] = Array(capacity=n_components)
      for c = 0; c < n_components; c = c + 1 {
        row.push(0.0)
      }
      new_pos.push(row)
    }
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        if i == j {
          continue
        }
        // current embedding distance (clipped)
        let mut d2 = 0.0
        for c = 0; c < n_components; c = c + 1 {
          let diff = pos[i][c] - pos[j][c]
          d2 = d2 + diff * diff
        }
        let dz = if d2 > 1.0e-24 { d2.sqrt() } else { 1.0e-12 }
        let ratio = target[i][j] / dz
        for c = 0; c < n_components; c = c + 1 {
          new_pos[i][c] = new_pos[i][c] + ratio * (pos[i][c] - pos[j][c])
        }
      }
    }
    for i = 0; i < n; i = i + 1 {
      for c = 0; c < n_components; c = c + 1 {
        pos[i][c] = new_pos[i][c] / n.to_double()
      }
    }
  }
  pos
}

///|
/// PHATE embedding of an `n x m` matrix.  `t = 0` selects the diffusion time
/// automatically as the smallest `t` whose von Neumann entropy falls to
/// `0.2 * VNE(1)` (80% entropy decay).  The embedding is the SMACOF MDS of
/// the diffusion distances.
pub fn scanpy_phate(
  x : Array[Array[Double]],
  n_components~ : Int32 = 2,
  n_neighbors~ : Int32 = 5,
  t~ : Int32 = 0,
  mds_iter~ : Int32 = 200,
  seed~ : Int32 = 0,
) -> ScanpyPhate {
  let n = x.length()
  if n == 0 {
    return ScanpyPhate::{ embedding: [], t: 0, vne: 0.0, }
  }
  let k = if n_neighbors < n { n_neighbors } else { n }
  let (knn_idx, knn_dists) = scanpy_knn_distances(x, k)
  let diff_op = scanpy_phate_diff_op(knn_idx, knn_dists, k, n)
  // symmetric similar operator D^-1/2 A D^-1/2 for the eigenvalues
  let deg = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let mut s = 0.0
    for j = 0; j < n; j = j + 1 {
      s = s + diff_op[i][j]
    }
    deg[i] = if s > 1.0e-12 { s } else { 1.0 }
  }
  let mat_vec = fn(v : Array[Double]) -> Array[Double] {
    let out = Array::make(n, 0.0)
    for i = 0; i < n; i = i + 1 {
      let mut acc = 0.0
      for j = 0; j < n; j = j + 1 {
        acc = acc + diff_op[i][j] / (deg[i] * deg[j]).sqrt() * v[j]
      }
      out[i] = acc
    }
    out
  }
  let nc = if n < 32 { n } else { 32 }
  let (evals, evecs) = linalg_lanczos_topk(mat_vec, n, nc)
  // auto diffusion time by von Neumann entropy decay
  let vne1 = scanpy_phate_vne(evals, 1)
  let t_sel = if t > 0 {
    t
  } else {
    let target = 0.2 * vne1
    let mut tt = 1
    while tt < 100 {
      if scanpy_phate_vne(evals, tt) <= target {
        break
      }
      tt = tt + 1
    }
    tt
  }
  // diffusion coordinates: coord_k(i) = lambda^t * phi_k(i) / sqrt(deg_i)
  let coords : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=nc)
    for c = 0; c < nc; c = c + 1 {
      let lam = if evals[c] > 0.0 { evals[c] } else { 0.0 }
      let lt = @math.pow(lam, t_sel.to_double())
      row.push(lt * evecs[i][c] / deg[i].sqrt())
    }
    coords.push(row)
  }
  // diffusion distance matrix
  let dist : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    dist.push(row)
  }
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let mut d2 = 0.0
      for c = 0; c < nc; c = c + 1 {
        let diff = coords[i][c] - coords[j][c]
        d2 = d2 + diff * diff
      }
      let d = d2.sqrt()
      dist[i][j] = d
      dist[j][i] = d
    }
  }
  let embedding = scanpy_phate_smacof(dist, n_components, mds_iter, seed)
  ScanpyPhate::{ embedding, t: t_sel, vne: scanpy_phate_vne(evals, t_sel), }
}
