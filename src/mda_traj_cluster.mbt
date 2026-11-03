/// Trajectory clustering by conformational similarity.
///
/// MDAnalysis equivalent: `encore.clustering.Cluster`
///
/// Groups trajectory frames into clusters by RMSD similarity,
/// identifying distinct conformational states.

///|
/// Cluster result for trajectory frames
pub struct MdaTrajCluster {
  assignments : Array[Int32] // cluster ID per frame
  centroids : Array[Int32] // frame index of each centroid
  sizes : Array[Int32] // number of frames per cluster
  n_clusters : Int32
}

///|
pub fn MdaTrajCluster::assignments(self : MdaTrajCluster) -> Array[Int32] {
  self.assignments
}

///|
pub fn MdaTrajCluster::centroids(self : MdaTrajCluster) -> Array[Int32] {
  self.centroids
}

///|
pub fn MdaTrajCluster::sizes(self : MdaTrajCluster) -> Array[Int32] {
  self.sizes
}

///|
pub fn MdaTrajCluster::n_clusters(self : MdaTrajCluster) -> Int32 {
  self.n_clusters
}

///|
/// Compute RMSD between two frames (each frame is atom×xyz).
fn mda_frame_rmsd(
  frame1 : Array[Array[Double]],
  frame2 : Array[Array[Double]],
) -> Double {
  let n = if frame1.length() < frame2.length() {
    frame1.length()
  } else {
    frame2.length()
  }
  if n == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for a = 0; a < n; a = a + 1 {
    let dx = frame1[a][0] - frame2[a][0]
    let dy = frame1[a][1] - frame2[a][1]
    let dz = frame1[a][2] - frame2[a][2]
    sum = sum + dx * dx + dy * dy + dz * dz
  }
  (sum / n.to_double()).sqrt()
}

///|
/// K-means clustering of trajectory frames by RMSD.
///
/// traj: frames[atom][xyz]
/// k: number of clusters
/// max_iter: maximum iterations
pub fn mda_trajectory_kmeans(
  traj : Array[Array[Array[Double]]],
  k : Int32,
  max_iter~ : Int32 = 100,
) -> MdaTrajCluster {
  let n = traj.length()
  if n == 0 || k <= 0 {
    return { assignments: [], centroids: [], sizes: [], n_clusters: 0, }
  }
  // Initialize centroids: evenly spaced frames
  let centroids : Array[Int32] = Array(capacity=0)
  for i = 0; i < k && i < n; i = i + 1 {
    centroids.push(i * (n - 1) / (k - 1))
  }
  let actual_k = centroids.length()
  let assignments : Array[Int32] = Array::make(n, 0)
  for iter = 0; iter < max_iter; iter = iter + 1 {
    let mut changed = false
    // Assign each frame to nearest centroid
    for f = 0; f < n; f = f + 1 {
      let mut best_c = 0
      let mut best_d = 10000000000.0
      for c = 0; c < actual_k; c = c + 1 {
        let d = mda_frame_rmsd(traj[f], traj[centroids[c]])
        if d < best_d {
          best_d = d
          best_c = c
        }
      }
      if assignments[f] != best_c {
        assignments[f] = best_c
        changed = true
      }
    }
    if !changed {
      break
    }
    // Update centroids: find frame closest to cluster mean
    for c = 0; c < actual_k; c = c + 1 {
      let mut best_f = centroids[c]
      let mut best_avg = 10000000000.0
      for f = 0; f < n; f = f + 1 {
        if assignments[f] == c {
          let mut sum_d = 0.0
          let mut count = 0.0
          for f2 = 0; f2 < n; f2 = f2 + 1 {
            if assignments[f2] == c {
              sum_d = sum_d + mda_frame_rmsd(traj[f], traj[f2])
              count = count + 1.0
            }
          }
          if count > 0.0 && sum_d / count < best_avg {
            best_avg = sum_d / count
            best_f = f
          }
        }
      }
      centroids[c] = best_f
    }
  }
  // Compute cluster sizes
  let sizes : Array[Int32] = Array::make(actual_k, 0)
  for f = 0; f < n; f = f + 1 {
    sizes[assignments[f]] = sizes[assignments[f]] + 1
  }
  { assignments, centroids, sizes, n_clusters: actual_k, }
}

///|
/// Hierarchical clustering (single linkage) by RMSD cutoff.
/// Frames within cutoff RMSD of any cluster member join that cluster.
pub fn mda_trajectory_cluster_by_cutoff(
  traj : Array[Array[Array[Double]]],
  cutoff : Double,
) -> MdaTrajCluster {
  let n = traj.length()
  if n == 0 {
    return { assignments: [], centroids: [], sizes: [], n_clusters: 0, }
  }
  let assignments : Array[Int32] = Array::make(n, -1)
  let centroids : Array[Int32] = Array(capacity=0)
  let mut cluster_id = 0
  for f = 0; f < n; f = f + 1 {
    if assignments[f] >= 0 {
      continue
    }
    // Start new cluster
    assignments[f] = cluster_id
    centroids.push(f)
    // Grow cluster: add frames within cutoff
    let mut added = true
    while added {
      added = false
      for f2 = 0; f2 < n; f2 = f2 + 1 {
        if assignments[f2] >= 0 {
          continue
        }
        // Check if close to any member of this cluster
        for f3 = 0; f3 < n; f3 = f3 + 1 {
          if assignments[f3] == cluster_id {
            let d = mda_frame_rmsd(traj[f2], traj[f3])
            if d <= cutoff {
              assignments[f2] = cluster_id
              added = true
              break
            }
          }
        }
      }
    }
    cluster_id = cluster_id + 1
  }
  let sizes : Array[Int32] = Array::make(cluster_id, 0)
  for f = 0; f < n; f = f + 1 {
    sizes[assignments[f]] = sizes[assignments[f]] + 1
  }
  { assignments, centroids, sizes, n_clusters: cluster_id, }
}

///|
/// Compute cluster populations (fractions).
pub fn mda_cluster_populations(clusters : MdaTrajCluster) -> Array[Double] {
  let total = clusters.assignments.length()
  if total == 0 {
    return []
  }
  let result : Array[Double] = Array::make(clusters.n_clusters, 0.0)
  for c = 0; c < clusters.n_clusters; c = c + 1 {
    result[c] = clusters.sizes[c].to_double() / total.to_double()
  }
  result
}
