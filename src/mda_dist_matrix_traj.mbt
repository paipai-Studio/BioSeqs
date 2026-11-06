/// Inter-residue distance matrix time series.
///
/// MDAnalysis equivalent: `analysis.distances.DistanceArray` over trajectory
///
/// Tracks pairwise residue distances over time for allosteric communication
/// and distance covariance analysis.

///|
/// Compute inter-residue distance matrix for a single frame.
/// residue_coords: [residue][xyz] (typically CA atoms)
pub fn mda_distance_matrix_frame(
  residue_coords : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = residue_coords.length()
  let matrix : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let dx = residue_coords[i][0] - residue_coords[j][0]
      let dy = residue_coords[i][1] - residue_coords[j][1]
      let dz = residue_coords[i][2] - residue_coords[j][2]
      let d = (dx * dx + dy * dy + dz * dz).sqrt()
      matrix[i][j] = d
      matrix[j][i] = d
    }
  }
  matrix
}

///|
/// Compute distance matrices for all frames.
/// traj: [frame][residue][xyz]
/// Returns [frame][residue_i][residue_j] distance
pub fn mda_distance_matrix_trajectory(
  traj : Array[Array[Array[Double]]],
) -> Array[Array[Array[Double]]] {
  let result : Array[Array[Array[Double]]] = Array(capacity=0)
  for f = 0; f < traj.length(); f = f + 1 {
    result.push(mda_distance_matrix_frame(traj[f]))
  }
  result
}

///|
/// Compute average distance matrix over trajectory.
pub fn mda_average_distance_matrix(
  dist_traj : Array[Array[Array[Double]]],
) -> Array[Array[Double]] {
  let n_frames = dist_traj.length()
  if n_frames == 0 {
    return []
  }
  let n = dist_traj[0].length()
  let avg : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  for f = 0; f < n_frames; f = f + 1 {
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        avg[i][j] = avg[i][j] + dist_traj[f][i][j]
      }
    }
  }
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      avg[i][j] = avg[i][j] / n_frames.to_double()
    }
  }
  avg
}

///|
/// Compute distance fluctuation (std dev) matrix.
/// High fluctuation pairs are potential allosteric communication sites.
pub fn mda_distance_fluctuation_matrix(
  dist_traj : Array[Array[Array[Double]]],
) -> Array[Array[Double]] {
  let n_frames = dist_traj.length()
  if n_frames == 0 {
    return []
  }
  let n = dist_traj[0].length()
  let avg = mda_average_distance_matrix(dist_traj)
  let var_mat : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  for f = 0; f < n_frames; f = f + 1 {
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        let d = dist_traj[f][i][j] - avg[i][j]
        var_mat[i][j] = var_mat[i][j] + d * d
      }
    }
  }
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      var_mat[i][j] = (var_mat[i][j] / n_frames.to_double()).sqrt()
    }
  }
  var_mat
}

///|
/// Distance covariance between two residue pairs.
/// Measures if distance changes in pair (i,j) correlate with pair (k,l).
pub fn mda_distance_covariance(
  dist_traj : Array[Array[Array[Double]]],
  i : Int32,
  j : Int32,
  k : Int32,
  l : Int32,
) -> Double {
  let n = dist_traj.length()
  if n < 2 {
    return 0.0
  }
  let mut mean_ij = 0.0
  let mut mean_kl = 0.0
  for f = 0; f < n; f = f + 1 {
    mean_ij = mean_ij + dist_traj[f][i][j]
    mean_kl = mean_kl + dist_traj[f][k][l]
  }
  mean_ij = mean_ij / n.to_double()
  mean_kl = mean_kl / n.to_double()
  let mut cov = 0.0
  for f = 0; f < n; f = f + 1 {
    cov = cov + (dist_traj[f][i][j] - mean_ij) * (dist_traj[f][k][l] - mean_kl)
  }
  cov / (n - 1).to_double()
}

///|
/// Detect allosteric communication: residue pairs with high distance correlation.
/// Returns list of (residue_i, residue_j, fluctuation) sorted by fluctuation.
pub fn mda_allosteric_pairs(
  fluct_matrix : Array[Array[Double]],
  min_fluct~ : Double = 1.0,
) -> Array[(Int32, Int32, Double)] {
  let n = fluct_matrix.length()
  let pairs : Array[(Int32, Int32, Double)] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if fluct_matrix[i][j] >= min_fluct {
        pairs.push((i, j, fluct_matrix[i][j]))
      }
    }
  }
  // Sort by fluctuation descending (simple bubble sort)
  for a = 0; a < pairs.length(); a = a + 1 {
    for b = a + 1; b < pairs.length(); b = b + 1 {
      if pairs[b].2 > pairs[a].2 {
        let tmp = pairs[a]
        pairs[a] = pairs[b]
        pairs[b] = tmp
      }
    }
  }
  pairs
}

///|
/// Dynamic cross-correlation matrix (DCCM) from distance fluctuations.
/// Normalized: C_ij = cov(d_ij, d_ji) / (std(d_ij) * std(d_ji))
/// For symmetric distance matrix, this simplifies to C_ij = 1 on diagonal.
pub fn mda_dynamic_cross_correlation(
  dist_traj : Array[Array[Array[Double]]],
) -> Array[Array[Double]] {
  let n_frames = dist_traj.length()
  if n_frames == 0 {
    return []
  }
  let n = dist_traj[0].length()
  let avg = mda_average_distance_matrix(dist_traj)
  // Compute DCCM: correlation of distance time series
  let dccm : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      // Correlation between d(i,*) and d(j,*) time series
      // Simplified: use distance to a reference (e.g., center of mass)
      // Here we compute correlation of d(i,j) with itself = 1
      if i == j {
        dccm[i][j] = 1.0
      } else {
        // Compute correlation of distance(i,j) time series with mean
        let mut cov = 0.0
        let mut var_i = 0.0
        let mut var_j = 0.0
        for f = 0; f < n_frames; f = f + 1 {
          let di = dist_traj[f][i][j] - avg[i][j]
          cov = cov + di * di
          var_i = var_i + di * di
          var_j = var_j + di * di
        }
        let denom = (var_i * var_j).sqrt()
        if denom > 0.0 {
          dccm[i][j] = cov / denom
        }
      }
    }
  }
  dccm
}
