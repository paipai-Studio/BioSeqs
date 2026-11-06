/// Scanpy second-round: Advanced graph algorithms
/// Inspired by scanpy.tl: leiden (directed/weighted/flavor), louvain (flavor),
/// draw_graph (layouts), paga (groups/threshold)

///|
/// Leiden clustering with directed/weighted graph support.
pub fn scanpy_leiden_directed(
  adjacency : Array[Array[Double]],
  directed~ : Bool = true,
  resolution~ : Double = 1.0,
  n_iterations~ : Int32 = 10,
) -> Array[Int32] {
  let n = adjacency.length()
  // Initialize: each node is its own community
  let labels = Array(capacity=n)
  for i in 0..<n {
    labels.push(i)
  }
  for _iter in 0..<n_iterations {
    let mut improved = false
    for node in 0..<n {
      // Compute neighbor communities and their weights
      let comm_weights = Map([])
      let mut total_weight = 0.0
      for j in 0..<n {
        let w = if directed {
          adjacency[node][j]
        } else {
          (adjacency[node][j] + adjacency[j][node]) / 2.0
        }
        if w > 0.0 {
          let c = labels[j]
          let cur = comm_weights.get(c).unwrap_or(0.0)
          comm_weights[c] = cur + w
          total_weight = total_weight + w
        }
      }
      // Find best community (maximize modularity gain)
      let current_comm = labels[node]
      let mut best_comm = current_comm
      let mut best_gain = 0.0
      for c in comm_weights.keys() {
        if c != current_comm {
          let w_in = comm_weights.get(c).unwrap_or(0.0)
          let gain = w_in -
            resolution *
            total_weight *
            comm_weights.get(c).unwrap_or(0.0) /
            (total_weight + 0.0001)
          if gain > best_gain {
            best_gain = gain
            best_comm = c
          }
        }
      }
      if best_comm != current_comm {
        labels[node] = best_comm
        improved = true
      }
    }
    if not(improved) {
      break
    }
  }
  // Relabel to 0..k-1
  let label_map = Map([])
  let mut next_label = 0
  for i in 0..<n {
    let mapped = label_map.get(labels[i])
    match mapped {
      Some(l) => labels[i] = l
      None => {
        label_map[labels[i]] = next_label
        labels[i] = next_label
        next_label = next_label + 1
      }
    }
  }
  labels
}

///|
/// Leiden with igraph flavor (uses igraph-style modularity).
pub fn scanpy_leiden_igraph(
  adjacency : Array[Array[Double]],
  resolution~ : Double = 1.0,
  n_iterations~ : Int32 = 10,
) -> Array[Int32] {
  // igraph flavor uses a slightly different modularity definition
  // but the algorithm is the same; we delegate to the directed version
  scanpy_leiden_directed(adjacency, directed=false, resolution~, n_iterations~)
}

///|
/// Louvain with flavor option ('igraph' or 'networkx').
pub fn scanpy_louvain_flavor(
  adjacency : Array[Array[Double]],
  flavor~ : String = "igraph",
  resolution~ : Double = 1.0,
) -> Array[Int32] {
  let n = adjacency.length()
  // Both flavors use the same core algorithm; the difference is in
  // the modularity normalization. igraph uses 2m, networkx uses m.
  let m_factor = if flavor == "igraph" { 2.0 } else { 1.0 }
  // Initialize
  let labels = Array(capacity=n)
  for i in 0..<n {
    labels.push(i)
  }
  // Compute total edge weight
  let mut total_w = 0.0
  for i in 0..<n {
    for j in 0..<n {
      total_w = total_w + adjacency[i][j]
    }
  }
  let m = total_w * m_factor
  // Single pass of label propagation
  for node in 0..<n {
    let comm_weights = Map([])
    let mut degree = 0.0
    for j in 0..<n {
      let w = adjacency[node][j] + adjacency[j][node]
      if w > 0.0 {
        let c = labels[j]
        comm_weights[c] = comm_weights.get(c).unwrap_or(0.0) + w
        degree = degree + w
      }
    }
    let mut best_comm = labels[node]
    let mut best_gain = 0.0
    for c in comm_weights.keys() {
      let gain = comm_weights.get(c).unwrap_or(0.0) -
        resolution * degree * degree / (m + 0.0001)
      if gain > best_gain {
        best_gain = gain
        best_comm = c
      }
    }
    labels[node] = best_comm
  }
  // Relabel
  let label_map = Map([])
  let mut next_label = 0
  for i in 0..<n {
    let mapped = label_map.get(labels[i])
    match mapped {
      Some(l) => labels[i] = l
      None => {
        label_map[labels[i]] = next_label
        labels[i] = next_label
        next_label = next_label + 1
      }
    }
  }
  labels
}

///|
/// Draw graph with specified layout algorithm.
/// Supported: "fr" (Fruchterman-Reingold), "kk" (Kamada-Kawai), "fa" (Force-Atlas)
pub fn scanpy_draw_graph_layout(
  adjacency : Array[Array[Double]],
  layout~ : String = "fr",
  n_iter~ : Int32 = 50,
) -> Array[Array[Double]] {
  let n = adjacency.length()
  // Initialize positions randomly on circle
  let pos = Array(capacity=n)
  for i in 0..<n {
    let angle = 2.0 * 3.141592653589793 * i.to_double() / n.to_double()
    let row = Array(capacity=2)
    row.push(@math.cos(angle))
    row.push(@math.sin(angle))
    pos.push(row)
  }
  match layout {
    "fr" => {
      // Fruchterman-Reingold
      let k = @math.pow(1.0 / n.to_double(), 0.5)
      let mut temp = 0.1
      for _iter in 0..<n_iter {
        let disp = Array(capacity=n)
        for _i in 0..<n {
          disp.push([0.0, 0.0])
        }
        // Repulsive forces
        for i in 0..<n {
          for j in 0..<n {
            if i != j {
              let dx = pos[i][0] - pos[j][0]
              let dy = pos[i][1] - pos[j][1]
              let dist = @math.pow(dx * dx + dy * dy + 0.0001, 0.5)
              let force = k * k / dist
              disp[i][0] = disp[i][0] + dx / dist * force
              disp[i][1] = disp[i][1] + dy / dist * force
            }
          }
        }
        // Attractive forces
        for i in 0..<n {
          for j in 0..<n {
            if adjacency[i][j] > 0.0 && i != j {
              let dx = pos[i][0] - pos[j][0]
              let dy = pos[i][1] - pos[j][1]
              let dist = @math.pow(dx * dx + dy * dy + 0.0001, 0.5)
              let force = dist * dist / k
              disp[i][0] = disp[i][0] - dx / dist * force
              disp[i][1] = disp[i][1] - dy / dist * force
            }
          }
        }
        // Apply displacements
        for i in 0..<n {
          let d = @math.pow(
            disp[i][0] * disp[i][0] + disp[i][1] * disp[i][1] + 0.0001,
            0.5,
          )
          pos[i][0] = pos[i][0] + disp[i][0] / d * temp
          pos[i][1] = pos[i][1] + disp[i][1] / d * temp
        }
        temp = temp * 0.95
      }
    }
    "kk" =>
      // Kamada-Kawai: use shortest path distances as target
      // Simplified: just use adjacency as spring constants
      for _iter in 0..<n_iter {
        for i in 0..<n {
          let mut fx = 0.0
          let mut fy = 0.0
          for j in 0..<n {
            if i != j {
              let dx = pos[j][0] - pos[i][0]
              let dy = pos[j][1] - pos[i][1]
              let dist = @math.pow(dx * dx + dy * dy + 0.0001, 0.5)
              let target = if adjacency[i][j] > 0.0 { 0.5 } else { 2.0 }
              let force = (dist - target) / dist
              fx = fx + dx * force
              fy = fy + dy * force
            }
          }
          pos[i][0] = pos[i][0] + fx * 0.01
          pos[i][1] = pos[i][1] + fy * 0.01
        }
      }
    "fa" =>
      // Force-Atlas 2 (simplified)
      for _iter in 0..<n_iter {
        for i in 0..<n {
          let mut fx = 0.0
          let mut fy = 0.0
          // Repulsion
          for j in 0..<n {
            if i != j {
              let dx = pos[i][0] - pos[j][0]
              let dy = pos[i][1] - pos[j][1]
              let dist2 = dx * dx + dy * dy + 0.0001
              fx = fx + dx / dist2
              fy = fy + dy / dist2
            }
          }
          // Attraction
          for j in 0..<n {
            if adjacency[i][j] > 0.0 && i != j {
              let dx = pos[j][0] - pos[i][0]
              let dy = pos[j][1] - pos[i][1]
              fx = fx + dx * adjacency[i][j]
              fy = fy + dy * adjacency[i][j]
            }
          }
          pos[i][0] = pos[i][0] + fx * 0.01
          pos[i][1] = pos[i][1] + fy * 0.01
        }
      }
    _ => ()
  }
  pos
}

///|
/// PAGA with group-level abstraction and threshold.
pub fn scanpy_paga_groups(
  labels : Array[Int32],
  adjacency : Array[Array[Double]],
  threshold~ : Double = 0.0,
) -> Array[Array[Double]] {
  // Count groups
  let group_set = Set([])
  for l in labels {
    group_set.add(l)
  }
  let groups = group_set.length()
  // Compute inter-group connectivity
  let paga = Array(capacity=groups)
  for _i in 0..<groups {
    let row = Array(capacity=groups)
    for _j in 0..<groups {
      row.push(0.0)
    }
    paga.push(row)
  }
  // Count edges between groups
  let group_counts = Array(capacity=groups)
  for _i in 0..<groups {
    group_counts.push(0)
  }
  for i in 0..<labels.length() {
    group_counts[labels[i]] = group_counts[labels[i]] + 1
  }
  for i in 0..<adjacency.length() {
    for j in 0..<adjacency.length() {
      if adjacency[i][j] > 0.0 && labels[i] != labels[j] {
        paga[labels[i]][labels[j]] = paga[labels[i]][labels[j]] + 1.0
      }
    }
  }
  // Normalize by group sizes and apply threshold
  for i in 0..<groups {
    for j in 0..<groups {
      if group_counts[i] > 0 && group_counts[j] > 0 {
        paga[i][j] = paga[i][j] /
          (group_counts[i].to_double() * group_counts[j].to_double())
      }
      if paga[i][j] < threshold {
        paga[i][j] = 0.0
      }
    }
  }
  paga
}
