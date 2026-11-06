/// mda_clustering.mbt — Trajectory clustering, average structure, density peaks, cluster centroids
/// Item 22-25 (clustering, averaging, density peaks, centroids)

///|
/// Cluster result
pub struct MdaClusterResult {
  n_clusters : Int32
  assignments : Array[Int32] // cluster index for each frame
  centroids : Array[Int32] // frame index of each cluster centroid
  cluster_sizes : Array[Int32]
}

///|
/// K-means clustering on trajectory frames (using RMSD as distance metric).
pub fn mda_cluster_kmeans(
  frames : Array[Array[(Double, Double, Double)]],
  n_clusters~ : Int32 = 3,
  max_iter~ : Int32 = 50,
) -> MdaClusterResult {
  let n_frames = frames.length()
  if n_frames == 0 {
    return MdaClusterResult::{
      n_clusters: 0,
      assignments: Array(capacity=0),
      centroids: Array(capacity=0),
      cluster_sizes: Array(capacity=0),
    }
  }
  let nc = if n_clusters < n_frames { n_clusters } else { n_frames }
  // Initialize centroids: evenly spaced
  let centroids : Array[Int32] = Array(capacity=nc)
  for i = 0; i < nc; i = i + 1 {
    centroids.push(i * (n_frames - 1) / (nc - 1).max(1))
  }
  let assignments : Array[Int32] = Array(capacity=n_frames)
  for i = 0; i < n_frames; i = i + 1 {
    assignments.push(0)
  }
  // Iterate
  for iter = 0; iter < max_iter; iter = iter + 1 {
    let mut changed = false
    // Assign each frame to nearest centroid
    for f = 0; f < n_frames; f = f + 1 {
      let mut best_cluster = 0
      let mut best_dist = 999999.0
      for c = 0; c < nc; c = c + 1 {
        let d = mda_cl_rmsd(frames[f], frames[centroids[c]])
        if d < best_dist {
          best_dist = d
          best_cluster = c
        }
      }
      if assignments[f] != best_cluster {
        assignments[f] = best_cluster
        changed = true
      }
    }
    if !changed {
      break
    }
    // Update centroids: find frame closest to cluster mean
    for c = 0; c < nc; c = c + 1 {
      let mut best_frame = centroids[c]
      let mut best_avg_dist = 999999.0
      for f = 0; f < n_frames; f = f + 1 {
        if assignments[f] != c {
          continue
        }
        let mut avg_dist = 0.0
        let mut count = 0
        for f2 = 0; f2 < n_frames; f2 = f2 + 1 {
          if assignments[f2] != c {
            continue
          }
          avg_dist = avg_dist + mda_cl_rmsd(frames[f], frames[f2])
          count = count + 1
        }
        if count > 0 {
          avg_dist = avg_dist / count.to_double()
        }
        if avg_dist < best_avg_dist {
          best_avg_dist = avg_dist
          best_frame = f
        }
      }
      centroids[c] = best_frame
    }
  }
  // Cluster sizes
  let cluster_sizes : Array[Int32] = Array(capacity=nc)
  for c = 0; c < nc; c = c + 1 {
    cluster_sizes.push(0)
  }
  for f = 0; f < n_frames; f = f + 1 {
    cluster_sizes[assignments[f]] = cluster_sizes[assignments[f]] + 1
  }
  MdaClusterResult::{ n_clusters: nc, assignments, centroids, cluster_sizes, }
}

///|
/// Compute average structure across selected frames.
pub fn mda_average_structure(
  frames : Array[Array[(Double, Double, Double)]],
  frame_indices~ : Array[Int32] = Array(capacity=0),
) -> Array[(Double, Double, Double)] {
  let indices = if frame_indices.length() > 0 {
    frame_indices
  } else {
    let all : Array[Int32] = Array(capacity=frames.length())
    for i = 0; i < frames.length(); i = i + 1 {
      all.push(i)
    }
    all
  }
  if indices.length() == 0 || frames[0].length() == 0 {
    return Array(capacity=0)
  }
  let n_atoms = frames[0].length()
  let result : Array[(Double, Double, Double)] = Array(capacity=n_atoms)
  for a = 0; a < n_atoms; a = a + 1 {
    let mut sx = 0.0
    let mut sy = 0.0
    let mut sz = 0.0
    for idx in indices {
      if idx < frames.length() && a < frames[idx].length() {
        sx = sx + frames[idx][a].0
        sy = sy + frames[idx][a].1
        sz = sz + frames[idx][a].2
      }
    }
    let count = indices.length().to_double()
    result.push((sx / count, sy / count, sz / count))
  }
  result
}

///|
/// Density peak clustering (simplified Rodriguez-Laio method).
pub fn mda_density_peaks(
  frames : Array[Array[(Double, Double, Double)]],
  dc~ : Double = 2.0,
) -> MdaClusterResult {
  let n = frames.length()
  if n == 0 {
    return MdaClusterResult::{
      n_clusters: 0,
      assignments: Array(capacity=0),
      centroids: Array(capacity=0),
      cluster_sizes: Array(capacity=0),
    }
  }
  // Compute local density rho_i = count of frames within dc
  let rho : Array[Int32] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let mut count = 0
    for j = 0; j < n; j = j + 1 {
      if i == j {
        continue
      }
      if mda_cl_rmsd(frames[i], frames[j]) < dc {
        count = count + 1
      }
    }
    rho.push(count)
  }
  // Sort by density descending
  let order : Array[Int32] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    order.push(i)
  }
  order.sort_by(fn(a, b) {
    if rho[b] > rho[a] {
      1
    } else if rho[b] < rho[a] {
      -1
    } else {
      0
    }
  })
  // Pick cluster centers: top frames with high density and far from existing centers
  let centroids : Array[Int32] = Array(capacity=0)
  let assignments : Array[Int32] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    assignments.push(0)
  }
  // Simple: pick top 3 density frames as centers
  let n_centers = if n >= 3 { 3 } else { n }
  for c = 0; c < n_centers; c = c + 1 {
    centroids.push(order[c])
  }
  // Assign each frame to nearest center
  for f = 0; f < n; f = f + 1 {
    let mut best = 0
    let mut best_dist = 999999.0
    for c = 0; c < n_centers; c = c + 1 {
      let d = mda_cl_rmsd(frames[f], frames[centroids[c]])
      if d < best_dist {
        best_dist = d
        best = c
      }
    }
    assignments[f] = best
  }
  let cluster_sizes : Array[Int32] = Array(capacity=n_centers)
  for c = 0; c < n_centers; c = c + 1 {
    cluster_sizes.push(0)
  }
  for f = 0; f < n; f = f + 1 {
    cluster_sizes[assignments[f]] = cluster_sizes[assignments[f]] + 1
  }
  MdaClusterResult::{
    n_clusters: n_centers,
    assignments,
    centroids,
    cluster_sizes,
  }
}

///|
/// Compute cluster centroid coordinates (average structure per cluster).
pub fn mda_cluster_centroids(
  frames : Array[Array[(Double, Double, Double)]],
  cluster_result : MdaClusterResult,
) -> Array[Array[(Double, Double, Double)]] {
  let result : Array[Array[(Double, Double, Double)]] = Array(
    capacity=cluster_result.n_clusters,
  )
  for c = 0; c < cluster_result.n_clusters; c = c + 1 {
    let indices : Array[Int32] = Array(capacity=0)
    for f = 0; f < frames.length(); f = f + 1 {
      if cluster_result.assignments[f] == c {
        indices.push(f)
      }
    }
    result.push(mda_average_structure(frames, frame_indices=indices))
  }
  result
}

///|
/// Hierarchical clustering (simplified: single-linkage agglomerative).
pub fn mda_cluster_hierarchical(
  frames : Array[Array[(Double, Double, Double)]],
  n_clusters~ : Int32 = 3,
) -> MdaClusterResult {
  let n = frames.length()
  if n == 0 {
    return MdaClusterResult::{
      n_clusters: 0,
      assignments: Array(capacity=0),
      centroids: Array(capacity=0),
      cluster_sizes: Array(capacity=0),
    }
  }
  // Start with each frame in its own cluster
  let cluster_id : Array[Int32] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    cluster_id.push(i)
  }
  // Merge closest clusters until we have n_clusters
  let target = if n_clusters < n { n_clusters } else { n }
  let mut n_current = n
  while n_current > target {
    // Find two closest clusters (single linkage: min distance between members)
    let mut best_i = 0
    let mut best_j = 1
    let mut best_dist = 999999.0
    for i = 0; i < n; i = i + 1 {
      for j = i + 1; j < n; j = j + 1 {
        if cluster_id[i] == cluster_id[j] {
          continue
        }
        let d = mda_cl_rmsd(frames[i], frames[j])
        if d < best_dist {
          best_dist = d
          best_i = i
          best_j = j
        }
      }
    }
    // Merge cluster of best_j into cluster of best_i
    let old_id = cluster_id[best_j]
    let new_id = cluster_id[best_i]
    for i = 0; i < n; i = i + 1 {
      if cluster_id[i] == old_id {
        cluster_id[i] = new_id
      }
    }
    n_current = n_current - 1
  }
  // Remap cluster IDs to 0..n_clusters-1
  let id_map : Map[Int32, Int32] = Map([], capacity=32)
  let mut next_id = 0
  let assignments : Array[Int32] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let cid = cluster_id[i]
    match id_map.get(cid) {
      Some(mapped) => assignments.push(mapped)
      None => {
        id_map[cid] = next_id
        assignments.push(next_id)
        next_id = next_id + 1
      }
    }
  }
  // Centroids: frame closest to cluster average
  let centroids : Array[Int32] = Array(capacity=next_id)
  let cluster_sizes : Array[Int32] = Array(capacity=next_id)
  for c = 0; c < next_id; c = c + 1 {
    centroids.push(0)
    cluster_sizes.push(0)
  }
  for f = 0; f < n; f = f + 1 {
    cluster_sizes[assignments[f]] = cluster_sizes[assignments[f]] + 1
  }
  // Find centroid for each cluster
  for c = 0; c < next_id; c = c + 1 {
    let mut best_frame = 0
    let mut best_avg_dist = 999999.0
    for f = 0; f < n; f = f + 1 {
      if assignments[f] != c {
        continue
      }
      let mut avg_dist = 0.0
      let mut count = 0
      for f2 = 0; f2 < n; f2 = f2 + 1 {
        if assignments[f2] != c {
          continue
        }
        avg_dist = avg_dist + mda_cl_rmsd(frames[f], frames[f2])
        count = count + 1
      }
      if count > 0 {
        avg_dist = avg_dist / count.to_double()
      }
      if avg_dist < best_avg_dist {
        best_avg_dist = avg_dist
        best_frame = f
      }
    }
    centroids[c] = best_frame
  }
  MdaClusterResult::{
    n_clusters: next_id,
    assignments,
    centroids,
    cluster_sizes,
  }
}

///|
fn mda_cl_rmsd(
  a : Array[(Double, Double, Double)],
  b : Array[(Double, Double, Double)],
) -> Double {
  let n = a.length()
  if n == 0 || n != b.length() {
    return 999999.0
  }
  let mut sum_sq = 0.0
  for i = 0; i < n; i = i + 1 {
    let dx = a[i].0 - b[i].0
    let dy = a[i].1 - b[i].1
    let dz = a[i].2 - b[i].2
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
  }
  @math.pow(sum_sq / n.to_double(), 0.5)
}

///|
pub fn mda_clustering_dummy() -> Unit {
  ignore(MdaClusterResult::{
    n_clusters: 0,
    assignments: Array(capacity=0),
    centroids: Array(capacity=0),
    cluster_sizes: Array(capacity=0),
  })
}
