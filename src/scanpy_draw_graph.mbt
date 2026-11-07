/// Scanpy-inspired: tl.draw_graph — force-directed graph layout
/// Uses a simplified Fruchterman-Reingold force-directed algorithm.

///|
/// Draw graph using force-directed layout.
///
/// Parameters:
/// - `connectivities`: adjacency/connectivity matrix (n × n)
/// - `n_iter`: number of iterations (default 500)
/// - `init_pos`: initial positions (optional, random if empty)
/// - `seed`: random seed
pub fn scanpy_draw_graph(
  connectivities : Array[Array[Double]],
  n_iter~ : Int32 = 500,
  init_pos~ : Array[Array[Double]] = [],
  seed~ : Int32 = 0,
) -> Array[Array[Double]] {
  let n = connectivities.length()
  if n == 0 {
    return []
  }
  let dim = 2 // 2D layout
  // Initialize positions
  let pos = Array(capacity=n)
  if init_pos.length() == n {
    for i in 0..<n {
      pos.push(Array::from_iter(init_pos[i].iter()))
    }
  } else {
    // Random initialization in unit square
    let mut rng = seed
    for _ in 0..<n {
      let row = Array(capacity=0)
      for _ in 0..<dim {
        rng = (rng * 1103515245 + 12345) & 0x7FFFFFFF
        row.push(rng.to_double() / 2147483647.0)
      }
      pos.push(row)
    }
  }
  // Optimal distance
  let area = 1.0
  let k = (area / n.to_double()).sqrt()
  let k_sq = k * k
  // Iterate
  for iter in 0..<n_iter {
    let temperature = 1.0 * (1.0 - iter.to_double() / n_iter.to_double())
    // Compute repulsive forces (all pairs)
    let disp = Array(capacity=n)
    for _ in 0..<n {
      disp.push(
        {
          let a = Array(capacity=0)
          for _ in 0..<dim {
            a.push(0.0)
          }
          a
        },
      )
    }
    for i in 0..<n {
      for j in 0..<n {
        if i != j {
          let delta : Array[Double] = []
          let mut dist_sq = 0.0
          for d in 0..<dim {
            let diff = pos[i][d] - pos[j][d]
            delta.push(diff)
            dist_sq += diff * diff
          }
          if dist_sq < 0.01 {
            dist_sq = 0.01
          }
          let dist = dist_sq.sqrt()
          // Repulsive force: k^2 / dist
          let force = k_sq / dist
          for d in 0..<dim {
            disp[i][d] += delta[d] / dist * force
          }
        }
      }
    }
    // Compute attractive forces (connected pairs)
    for i in 0..<n {
      for j in 0..<n {
        if i != j && connectivities[i][j] > 0.0 {
          let delta : Array[Double] = []
          let mut dist_sq = 0.0
          for d in 0..<dim {
            let diff = pos[j][d] - pos[i][d]
            delta.push(diff)
            dist_sq += diff * diff
          }
          if dist_sq < 0.01 {
            dist_sq = 0.01
          }
          let dist = dist_sq.sqrt()
          // Attractive force: dist^2 / k * weight
          let force = dist / k * connectivities[i][j]
          for d in 0..<dim {
            disp[i][d] += delta[d] / dist * force
          }
        }
      }
    }
    // Apply forces with temperature cooling
    for i in 0..<n {
      let mut disp_mag = 0.0
      for d in 0..<dim {
        disp_mag += disp[i][d] * disp[i][d]
      }
      disp_mag = disp_mag.sqrt()
      if disp_mag > 0.0 {
        let limited = if disp_mag > temperature {
          temperature / disp_mag
        } else {
          1.0
        }
        for d in 0..<dim {
          pos[i][d] += disp[i][d] / disp_mag * limited
        }
      }
    }
  }
  pos
}
