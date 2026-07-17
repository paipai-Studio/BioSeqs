///|
/// Cluster - Sequence clustering analysis inspired by Biopython Bio.Cluster.

///|
pub struct ClusterResult {
  clusters : Array[Array[Int]]
  distances : Array[Array[Double]]
  cluster_labels : Array[Int]
}

///|
pub fn cluster_result(
  clusters : Array[Array[Int]],
  distances : Array[Array[Double]],
  cluster_labels : Array[Int],
) -> ClusterResult {
  ClusterResult::{ clusters, distances, cluster_labels }
}

///|
pub fn seq_distance(seq1 : String, seq2 : String) -> Double {
  let len1 = seq1.length()
  let len2 = seq2.length()
  let min_len = if len1 < len2 { len1 } else { len2 }
  let max_len = if len1 > len2 { len1 } else { len2 }

  let mut mismatches = 0.0
  for i = 0; i < min_len; i = i + 1 {
    if seq1.unsafe_get(i) != seq2.unsafe_get(i) {
      mismatches = mismatches + 1.0
    }
  }

  (mismatches + (max_len - min_len).to_double()) / max_len.to_double()
}

///|
pub fn distance_matrix(sequences : Array[String]) -> Array[Array[Double]] {
  let n = sequences.length()
  let matrix : Array[Array[Double]] = Array::new()

  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::new()
    for j = 0; j < n; j = j + 1 {
      if i == j {
        row.push(0.0)
      } else {
        row.push(seq_distance(sequences[i], sequences[j]))
      }
    }
    matrix.push(row)
  }

  matrix
}

///|
pub fn single_linkage_clustering(
  distances : Array[Array[Double]],
  threshold : Double,
) -> ClusterResult {
  let n = distances.length()
  let cluster_labels : Array[Int] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    cluster_labels[i] = i
  }

  let mut n_clusters = n

  let mut changed = true
  while changed {
    changed = false
    let mut min_dist = 2.0
    let mut min_i = 0
    let mut min_j = 0

    for i = 0; i < n; i = i + 1 {
      for j = i + 1; j < n; j = j + 1 {
        if cluster_labels[i] != cluster_labels[j] && distances[i][j] < min_dist {
          min_dist = distances[i][j]
          min_i = i
          min_j = j
        }
      }
    }

    if min_dist <= threshold {
      let old_label = cluster_labels[min_j]
      let new_label = cluster_labels[min_i]

      for i = 0; i < n; i = i + 1 {
        if cluster_labels[i] == old_label {
          cluster_labels[i] = new_label
        }
      }

      n_clusters = n_clusters - 1
      changed = true
    }
  }

  let clusters : Array[Array[Int]] = Array::make(n_clusters, Array::new())
  let label_map : Map[Int, Int] = Map([], capacity=n_clusters)
  let mut current_idx = 0

  for i = 0; i < n; i = i + 1 {
    match label_map.get(cluster_labels[i]) {
      Some(idx) => clusters[idx].push(i)
      None => {
        label_map[cluster_labels[i]] = current_idx
        clusters[current_idx] = [i]
        current_idx = current_idx + 1
      }
    }
  }

  ClusterResult::{ clusters, distances, cluster_labels }
}

///|
pub fn complete_linkage_clustering(
  distances : Array[Array[Double]],
  threshold : Double,
) -> ClusterResult {
  let n = distances.length()
  let cluster_labels : Array[Int] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    cluster_labels[i] = i
  }

  let mut changed = true
  while changed {
    changed = false
    let mut min_max_dist = 2.0
    let mut min_i = 0
    let mut min_j = 0

    for i = 0; i < n; i = i + 1 {
      for j = i + 1; j < n; j = j + 1 {
        if cluster_labels[i] != cluster_labels[j] {
          let mut max_dist = 0.0

          for a = 0; a < n; a = a + 1 {
            for b = 0; b < n; b = b + 1 {
              if cluster_labels[a] == cluster_labels[i] &&
                cluster_labels[b] == cluster_labels[j] {
                if distances[a][b] > max_dist {
                  max_dist = distances[a][b]
                }
              }
            }
          }

          if max_dist < min_max_dist {
            min_max_dist = max_dist
            min_i = i
            min_j = j
          }
        }
      }
    }

    if min_max_dist <= threshold {
      let old_label = cluster_labels[min_j]
      let new_label = cluster_labels[min_i]

      for i = 0; i < n; i = i + 1 {
        if cluster_labels[i] == old_label {
          cluster_labels[i] = new_label
        }
      }

      changed = true
    }
  }

  let labels_seen : Map[Int, Bool] = Map([], capacity=n)
  for i = 0; i < n; i = i + 1 {
    labels_seen[cluster_labels[i]] = true
  }

  let n_clusters = labels_seen.keys().collect().length()
  let clusters : Array[Array[Int]] = Array::make(n_clusters, Array::new())
  let label_map : Map[Int, Int] = Map([], capacity=n_clusters)
  let mut current_idx = 0

  for i = 0; i < n; i = i + 1 {
    match label_map.get(cluster_labels[i]) {
      Some(idx) => clusters[idx].push(i)
      None => {
        label_map[cluster_labels[i]] = current_idx
        clusters[current_idx] = [i]
        current_idx = current_idx + 1
      }
    }
  }

  ClusterResult::{ clusters, distances, cluster_labels }
}

///|
pub fn average_linkage_clustering(
  distances : Array[Array[Double]],
  threshold : Double,
) -> ClusterResult {
  let n = distances.length()
  let cluster_labels : Array[Int] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    cluster_labels[i] = i
  }

  let mut changed = true
  while changed {
    changed = false
    let mut min_avg_dist = 2.0
    let mut min_i = 0
    let mut min_j = 0

    for i = 0; i < n; i = i + 1 {
      for j = i + 1; j < n; j = j + 1 {
        if cluster_labels[i] != cluster_labels[j] {
          let mut sum_dist = 0.0
          let mut count = 0

          for a = 0; a < n; a = a + 1 {
            for b = 0; b < n; b = b + 1 {
              if cluster_labels[a] == cluster_labels[i] &&
                cluster_labels[b] == cluster_labels[j] {
                sum_dist = sum_dist + distances[a][b]
                count = count + 1
              }
            }
          }

          let avg_dist = sum_dist / count.to_double()

          if avg_dist < min_avg_dist {
            min_avg_dist = avg_dist
            min_i = i
            min_j = j
          }
        }
      }
    }

    if min_avg_dist <= threshold {
      let old_label = cluster_labels[min_j]
      let new_label = cluster_labels[min_i]

      for i = 0; i < n; i = i + 1 {
        if cluster_labels[i] == old_label {
          cluster_labels[i] = new_label
        }
      }

      changed = true
    }
  }

  let clusters : Array[Array[Int]] = Array::make(n, Array::new())
  let label_map : Map[Int, Int] = Map([], capacity=n)
  let mut current_idx = 0

  for i = 0; i < n; i = i + 1 {
    match label_map.get(cluster_labels[i]) {
      Some(idx) => clusters[idx].push(i)
      None => {
        label_map[cluster_labels[i]] = current_idx
        clusters[current_idx] = [i]
        current_idx = current_idx + 1
      }
    }
  }

  let final_clusters : Array[Array[Int]] = Array::make(
    current_idx,
    Array::new(),
  )
  for i = 0; i < current_idx; i = i + 1 {
    final_clusters[i] = clusters[i]
  }

  ClusterResult::{ clusters: final_clusters, distances, cluster_labels }
}

///|
pub fn hclust_to_newick(
  result : ClusterResult,
  names : Array[String],
) -> String {
  let clusters = result.clusters

  if clusters.length() == 0 {
    return ""
  }

  let mut tree_str = "("

  for i = 0; i < clusters.length(); i = i + 1 {
    if clusters[i].length() == 1 {
      tree_str = tree_str + names[clusters[i][0]]
    } else {
      tree_str = tree_str + "("
      for j = 0; j < clusters[i].length(); j = j + 1 {
        tree_str = tree_str + names[clusters[i][j]]
        if j < clusters[i].length() - 1 {
          tree_str = tree_str + ","
        }
      }
      tree_str = tree_str + ")"
    }

    if i < clusters.length() - 1 {
      tree_str = tree_str + ","
    }
  }

  tree_str = tree_str + ");"
  tree_str
}

///|
pub fn seq_silhouette_score(
  sequences : Array[String],
  labels : Array[Int],
) -> Double {
  let n = sequences.length()
  if n < 2 {
    return 0.0
  }

  let distances = distance_matrix(sequences)
  let mut sum_silhouette = 0.0

  for i = 0; i < n; i = i + 1 {
    let mut a = 0.0
    let mut a_count = 0
    let mut b = 2.0

    for j = 0; j < n; j = j + 1 {
      if i != j {
        if labels[j] == labels[i] {
          a = a + distances[i][j]
          a_count = a_count + 1
        } else if distances[i][j] < b {
          b = distances[i][j]
        }
      }
    }

    if a_count > 0 {
      a = a / a_count.to_double()
    }

    let s = (b - a) / (if a > b { a } else { b })
    sum_silhouette = sum_silhouette + s
  }

  sum_silhouette / n.to_double()
}

///|
pub fn cluster_size(cluster : Array[Int]) -> Int {
  cluster.length()
}

///|
pub fn total_clusters(result : ClusterResult) -> Int {
  result.clusters.length()
}
