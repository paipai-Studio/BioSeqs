/// AnnData obsp/varp pairwise distance storage and operations.
///
/// `adata.obsp["distances"]` and `adata.obsp["connectivities"]` equivalents.
///
/// Stores pairwise matrices for obs (n_obs × n_obs) and var (n_vars × n_vars).

///|
/// Pairwise matrix type tag.
pub enum AnnDataPairwiseType {
  PairwiseDistances // distance matrix
  PairwiseConnectivities // connectivity/adjacency matrix
  PairwiseAffinity // affinity/similarity matrix
  PairwiseCustom(String) // custom named type
}

///|
/// A named pairwise matrix.
pub struct AnnDataPairwise {
  name : String
  matrix : Array[Array[Double]] // n × n
  kind : AnnDataPairwiseType
}

///|
pub fn AnnDataPairwise::name(self : AnnDataPairwise) -> String {
  self.name
}

///|
pub fn AnnDataPairwise::matrix(self : AnnDataPairwise) -> Array[Array[Double]] {
  self.matrix
}

///|
/// Get element (i, j) from pairwise matrix.
pub fn AnnDataPairwise::get(
  self : AnnDataPairwise,
  i : Int32,
  j : Int32,
) -> Double {
  if i < self.matrix.length() && j < self.matrix[i].length() {
    self.matrix[i][j]
  } else {
    0.0
  }
}

///|
/// Set element (i, j) in pairwise matrix.
pub fn AnnDataPairwise::set(
  self : AnnDataPairwise,
  i : Int32,
  j : Int32,
  v : Double,
) -> Unit {
  if i < self.matrix.length() && j < self.matrix[i].length() {
    self.matrix[i][j] = v
  }
}

///|
/// Size (n × n).
pub fn AnnDataPairwise::size(self : AnnDataPairwise) -> Int32 {
  self.matrix.length()
}

///|
/// Check if matrix is symmetric (M[i][j] == M[j][i]).
pub fn AnnDataPairwise::is_symmetric(self : AnnDataPairwise) -> Bool {
  let n = self.matrix.length()
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if i < self.matrix.length() && j < self.matrix[i].length() {
        if j < self.matrix.length() && i < self.matrix[j].length() {
          if (self.matrix[i][j] - self.matrix[j][i]).abs() > 0.000000000001 {
            return false
          }
        }
      }
    }
  }
  true
}

///|
/// Check if diagonal is all zero.
pub fn AnnDataPairwise::is_zero_diagonal(self : AnnDataPairwise) -> Bool {
  let n = self.matrix.length()
  for i = 0; i < n; i = i + 1 {
    if i < self.matrix.length() && i < self.matrix[i].length() {
      if self.matrix[i][i].abs() > 0.000000000001 {
        return false
      }
    }
  }
  true
}

///|
/// Create a distance matrix from obs coordinates (obsm).
/// Uses Euclidean distance.
pub fn anndata_pairwise_distances_euclidean(
  coords : Array[Array[Double]],
) -> AnnDataPairwise {
  let n = coords.length()
  let mat : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let mut sum_sq = 0.0
      let ni = coords[i].length()
      let nj = coords[j].length()
      let d = if ni < nj { ni } else { nj }
      for k = 0; k < d; k = k + 1 {
        let diff = coords[i][k] - coords[j][k]
        sum_sq = sum_sq + diff * diff
      }
      let dist = sum_sq.sqrt()
      mat[i][j] = dist
      mat[j][i] = dist
    }
  }
  { name: "distances", matrix: mat, kind: PairwiseDistances, }
}

///|
/// Create k-NN connectivity matrix from a distance matrix.
/// Each row has 1.0 for its k nearest neighbors.
pub fn anndata_pairwise_connectivities_knn(
  distances : AnnDataPairwise,
  k : Int32,
) -> AnnDataPairwise {
  let n = distances.size()
  let mat : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  for i = 0; i < n; i = i + 1 {
    // Collect (distance, index) pairs
    let pairs : Array[(Double, Int32)] = Array(capacity=0)
    for j = 0; j < n; j = j + 1 {
      if j != i {
        pairs.push((distances.get(i, j), j))
      }
    }
    // Sort by distance
    pairs.sort_by(fn(a, b) {
      if a.0 < b.0 {
        -1
      } else if a.0 > b.0 {
        1
      } else {
        0
      }
    })
    // Take top k
    let kk = if k < pairs.length() { k } else { pairs.length() }
    for idx = 0; idx < kk; idx = idx + 1 {
      let j = pairs[idx].1
      mat[i][j] = 1.0
      mat[j][i] = 1.0 // symmetric
    }
  }
  { name: "connectivities", matrix: mat, kind: PairwiseConnectivities, }
}

///|
/// Convert distance matrix to affinity (Gaussian kernel).
/// affinity = exp(-d^2 / (2 * sigma^2))
pub fn anndata_pairwise_to_affinity(
  distances : AnnDataPairwise,
  sigma : Double,
) -> AnnDataPairwise {
  let n = distances.size()
  let mat : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  let denom = 2.0 * sigma * sigma
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      let d = distances.get(i, j)
      mat[i][j] = @math.exp(-(d * d) / denom)
    }
  }
  { name: "affinity", matrix: mat, kind: PairwiseAffinity, }
}

///|
/// Extract upper triangle as flat array (for condensed distance).
pub fn anndata_pairwise_condensed(distances : AnnDataPairwise) -> Array[Double] {
  let n = distances.size()
  let out : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      out.push(distances.get(i, j))
    }
  }
  out
}

///|
/// Subset pairwise matrix by indices.
pub fn AnnDataPairwise::subset(
  self : AnnDataPairwise,
  idx : Array[Int32],
) -> AnnDataPairwise {
  let n = idx.length()
  let mat : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      mat[i][j] = self.get(idx[i], idx[j])
    }
  }
  { name: self.name, matrix: mat, kind: self.kind, }
}
