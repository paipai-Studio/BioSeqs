/// Scanpy-inspired: tl.embedding_density — compute density in embedding space

///|
/// Compute cell density in an embedding (e.g., UMAP, t-SNE).
/// Uses Gaussian kernel density estimation.
///
/// Parameters:
/// - `embedding`: n_cells × n_dims (e.g., UMAP coordinates)
/// - `groupby`: optional group labels per cell (if empty, compute global density)
/// - `bandwidth`: KDE bandwidth (default 0.5; if 0, auto-estimate)
pub fn scanpy_embedding_density(
  embedding : Array[Array[Double]],
  _groupby~ : Array[String] = [],
  bandwidth~ : Double = 0.5,
) -> Array[Double] {
  let n = embedding.length()
  if n == 0 {
    return []
  }
  let dim = embedding[0].length()
  // Auto-estimate bandwidth if 0: use mean nearest-neighbor distance
  let bw = if bandwidth > 0.0 {
    bandwidth
  } else {
    let mut total_nn_dist = 0.0
    for i in 0..<n {
      let mut min_dist = 1.0e100
      for j in 0..<n {
        if i != j {
          let mut d = 0.0
          for k in 0..<dim {
            let diff = embedding[i][k] - embedding[j][k]
            d += diff * diff
          }
          if d < min_dist {
            min_dist = d
          }
        }
      }
      total_nn_dist += min_dist.sqrt()
    }
    total_nn_dist / n.to_double()
  }
  if bw <= 0.0 {
    return {
      let a = Array(capacity=0)
      for _ in 0..<n {
        a.push(0.0)
      }
      a
    }
  }
  // Compute Gaussian KDE: density_i = sum_j exp(-||x_i - x_j||^2 / (2*bw^2))
  let two_bw_sq = 2.0 * bw * bw
  let densities = Array(capacity=n)
  for i in 0..<n {
    let mut density = 0.0
    for j in 0..<n {
      let mut d_sq = 0.0
      for k in 0..<dim {
        let diff = embedding[i][k] - embedding[j][k]
        d_sq += diff * diff
      }
      density += @math.exp(-d_sq / two_bw_sq)
    }
    // Normalize by n and bandwidth
    densities.push(
      density /
      (
        n.to_double() *
        @math.pow(2.0 * 3.141592653589793, dim.to_double() / 2.0) *
        @math.pow(bw, dim.to_double())
      ),
    )
  }
  densities
}

///|
/// Compute embedding density per group.
/// Returns a map from group label to density array.
pub fn scanpy_embedding_density_by_group(
  embedding : Array[Array[Double]],
  groupby : Array[String],
  bandwidth~ : Double = 0.5,
) -> Map[String, Array[Double]] {
  let n = embedding.length()
  let result : Map[String, Array[Double]] = Map([], capacity=0)
  if n == 0 {
    return result
  }
  // Get unique groups
  let group_set : Set[String] = Set::from_iter(groupby.iter())
  for group in group_set {
    // Collect embedding for this group
    let group_emb = Array(capacity=0)
    for i in 0..<n {
      if groupby[i] == group {
        group_emb.push(embedding[i])
      }
    }
    let d = scanpy_embedding_density(group_emb, bandwidth~)
    result[group] = d
  }
  result
}
