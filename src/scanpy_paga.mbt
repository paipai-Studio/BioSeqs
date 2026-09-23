///|
/// scanpy `tl.paga` — partition-based graph abstraction: connectivity between
/// clusters, normalized by the random-graph null model.

///|
pub fn scanpy_paga(
  connectivities : Array[Array[Double]],
  labels : Array[Int],
) -> Array[Array[Double]] {
  let n = connectivities.length()
  if n == 0 {
    return []
  }
  let mut max_label = 0
  for l in labels {
    if l > max_label {
      max_label = l
    }
  }
  let n_clusters = max_label + 1

  let ns = Array::make(n_clusters, 0)
  for l in labels {
    ns[l] = ns[l] + 1
  }

  // degree + incident weight per cluster
  let deg = Array::make(n, 0.0)
  for u = 0; u < n; u = u + 1 {
    for v = 0; v < n; v = v + 1 {
      deg[u] = deg[u] + connectivities[u][v]
    }
  }
  let es = Array::make(n_clusters, 0.0)
  for u = 0; u < n; u = u + 1 {
    es[labels[u]] = es[labels[u]] + deg[u]
  }

  // inter-cluster weight
  let inter : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n_clusters; i = i + 1 {
    inter.push(Array::make(n_clusters, 0.0))
  }
  for u = 0; u < n; u = u + 1 {
    for v = 0; v < n; v = v + 1 {
      inter[labels[u]][labels[v]] = inter[labels[u]][labels[v]] + connectivities[u][v]
    }
  }

  // connectivity = inter / null-model expectation, clamped to [0,1]
  let conn : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n_clusters; i = i + 1 {
    conn.push(Array::make(n_clusters, 0.0))
  }
  for i = 0; i < n_clusters; i = i + 1 {
    for j = 0; j < n_clusters; j = j + 1 {
      if i == j {
        continue
      }
      let v = inter[i][j] + inter[j][i]
      let expected = (es[i] * ns[j].to_double() + es[j] * ns[i].to_double()) /
        (n - 1).to_double()
      let scaled = if expected != 0.0 { v / expected } else { 1.0 }
      conn[i][j] = if scaled > 1.0 { 1.0 } else if scaled < 0.0 { 0.0 } else { scaled }
    }
  }
  conn
}