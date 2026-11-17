///|
/// squidpy-inspired spatial statistics: neighborhood enrichment
/// (`sq.gr.nhood_enrichment`), co-occurrence (`sq.gr.co_occurrence`) and
/// Ripley's K/L/G (`sq.gr.ripley`). Operate on plain label / distance / coordinate
/// arrays; a symmetric adjacency or distance matrix can come from
/// `scanpy_spatial_neighbors`.

///|
/// Neighborhood enrichment result: raw neighbor-pair counts and their
/// z-scores, both indexed by the cluster labels in `clusters` order.
pub struct ScanpyNhoodEnrichment {
  clusters : Array[String]
  counts : Array[Array[Double]]
  zscores : Array[Array[Double]]
}

///|
/// Cluster label index: first-seen order over `labels`.
fn scanpy_spatial_cluster_levels(
  labels : Array[String],
) -> (Array[String], Map[String, Int32]) {
  let levels : Array[String] = Array(capacity=0)
  let idx : Map[String, Int32] = Map([], capacity=4)
  for l in labels {
    if !idx.contains(l) {
      idx[l] = levels.length()
      levels.push(l)
    }
  }
  (levels, idx)
}

///|
/// `sq.gr.nhood_enrichment` — z-scored enrichment of cluster pairs among
/// graph neighbors. `weights` is the symmetric connectivity matrix (e.g.
/// `SpatialNeighborsResult::connectivities`, diagonal ignored). For every
/// cluster pair the observed number of neighbor edges is compared against the
/// expectation under random label permutation (binomial approximation over the
/// total edge count), giving z = (observed - expected) / sd.
pub fn scanpy_nhood_enrichment(
  clusters : Array[String],
  weights : Array[Array[Double]],
) -> ScanpyNhoodEnrichment {
  let (levels, idx) = scanpy_spatial_cluster_levels(clusters)
  let k = levels.length()
  let n = clusters.length()
  let counts : Array[Array[Double]] = Array(capacity=k)
  for _a in 0..<k {
    counts.push(Array::make(k, 0.0))
  }
  // observed ordered neighbor pairs (i != j)
  let mut total_edges = 0.0
  for i = 0; i < n && i < weights.length(); i = i + 1 {
    let row = weights[i]
    let ci = match idx.get(clusters[i]) {
      Some(v) => v
      None => 0
    }
    for j = 0; j < n && j < row.length(); j = j + 1 {
      if i != j && row[j] > 0.0 {
        let cj = match idx.get(clusters[j]) {
          Some(v) => v
          None => 0
        }
        counts[ci][cj] = counts[ci][cj] + 1.0
        total_edges = total_edges + 1.0
      }
    }
  }
  // cluster sizes
  let sizes : Array[Double] = Array(capacity=k)
  for _a in 0..<k {
    sizes.push(0.0)
  }
  for c in clusters {
    let a = match idx.get(c) {
      Some(v) => v
      None => 0
    }
    sizes[a] = sizes[a] + 1.0
  }
  let nf = n.to_double()
  let zscores : Array[Array[Double]] = Array(capacity=k)
  for a = 0; a < k; a = a + 1 {
    let row : Array[Double] = Array(capacity=k)
    for b = 0; b < k; b = b + 1 {
      let na = sizes[a]
      let nb = sizes[b]
      let p = if a == b {
        na * (na - 1.0) / (nf * (nf - 1.0))
      } else {
        2.0 * na * nb / (nf * (nf - 1.0))
      }
      let expected = total_edges * p
      let sd = (total_edges * p * (1.0 - p)).sqrt()
      row.push(if sd > 0.0 { (counts[a][b] - expected) / sd } else { 0.0 })
    }
    zscores.push(row)
  }
  ScanpyNhoodEnrichment::{ clusters: levels, counts, zscores, }
}

///|
/// Co-occurrence result: `ratios[band][a][b]` is the ratio of the probability
/// of finding a cluster-`b` cell within the `band`-th distance interval of a
/// cluster-`a` cell, over the global fraction of cluster `b`. `intervals` holds
/// the band upper bounds (uniform split of the max nonzero distance).
pub struct ScanpyCoOccurrence {
  clusters : Array[String]
  intervals : Array[Double]
  ratios : Array[Array[Array[Double]]]
}

///|
/// `sq.gr.co_occurrence` — cluster co-occurrence probability ratios over
/// `n_bands` uniform distance intervals of the symmetric `distances` matrix.
/// ratio > 1 means cluster `b` appears near cluster `a` more often than expected.
pub fn scanpy_co_occurrence(
  clusters : Array[String],
  distances : Array[Array[Double]],
  n_bands~ : Int32 = 3,
) -> ScanpyCoOccurrence {
  let (levels, idx) = scanpy_spatial_cluster_levels(clusters)
  let k = levels.length()
  let n = clusters.length()
  // max nonzero distance
  let mut dmax = 0.0
  for i = 0; i < n && i < distances.length(); i = i + 1 {
    let row = distances[i]
    for j = 0; j < n && j < row.length(); j = j + 1 {
      if i != j && row[j] > dmax {
        dmax = row[j]
      }
    }
  }
  let bands : Array[Double] = Array(capacity=n_bands)
  let los : Array[Double] = Array(capacity=n_bands)
  let step = dmax / n_bands.to_double()
  for b = 0; b < n_bands; b = b + 1 {
    los.push(step * b.to_double())
    bands.push(step * (b + 1).to_double())
  }
  // cluster sizes
  let sizes : Array[Double] = Array(capacity=k)
  for _a in 0..<k {
    sizes.push(0.0)
  }
  for c in clusters {
    let a = match idx.get(c) {
      Some(v) => v
      None => 0
    }
    sizes[a] = sizes[a] + 1.0
  }
  let baseline : Array[Double] = Array(capacity=k)
  for b = 0; b < k; b = b + 1 {
    baseline.push(if n > 1 { sizes[b] / (n - 1).to_double() } else { 0.0 })
  }
  let ratios : Array[Array[Array[Double]]] = Array(capacity=n_bands)
  for band = 0; band < n_bands; band = band + 1 {
    let lo = los[band]
    let hi = bands[band]
    let mat : Array[Array[Double]] = Array(capacity=k)
    for a = 0; a < k; a = a + 1 {
      let row : Array[Double] = Array(capacity=k)
      for b = 0; b < k; b = b + 1 {
        row.push(0.0)
      }
      mat.push(row)
    }
    // accumulate conditional fractions per anchor cell, then average
    let sums : Array[Array[Double]] = Array(capacity=k)
    let n_a : Array[Int32] = Array(capacity=k)
    for _a in 0..<k {
      sums.push(Array::make(k, 0.0))
      n_a.push(0)
    }
    for i = 0; i < n && i < distances.length(); i = i + 1 {
      let ci = match idx.get(clusters[i]) {
        Some(v) => v
        None => 0
      }
      n_a[ci] = n_a[ci] + 1
      let row = distances[i]
      let in_band = Array::make(k, 0)
      let mut denom = 0
      for j = 0; j < n && j < row.length(); j = j + 1 {
        if i != j && row[j] > lo && row[j] <= hi {
          let cj = match idx.get(clusters[j]) {
            Some(v) => v
            None => 0
          }
          in_band[cj] = in_band[cj] + 1
          denom = denom + 1
        }
      }
      if denom > 0 && baseline[ci] > 0.0 {
        for b = 0; b < k; b = b + 1 {
          // per-anchor fraction of cluster b among in-band cells
          let frac = in_band[b].to_double() / denom.to_double()
          // normalized by baseline expectation for cluster b
          sums[ci][b] = sums[ci][b] + frac / baseline[b]
        }
      }
    }
    for a = 0; a < k; a = a + 1 {
      for b = 0; b < k; b = b + 1 {
        mat[a][b] = if n_a[a] > 0 {
          sums[a][b] / n_a[a].to_double()
        } else {
          0.0
        }
      }
    }
    ratios.push(mat)
  }
  ScanpyCoOccurrence::{ clusters: levels, intervals: bands, ratios, }
}

///|
/// Ripley's K / L result over a sequence of radii.
/// K(t) = area / (n(n-1)) * #{ordered pairs with distance <= t},
/// L(t) = sqrt(K(t)/pi) - t (variance-stabilizing transform, ~0 under CSR).
pub struct ScanpyRipleyK {
  radii : Array[Double]
  k : Array[Double]
  l : Array[Double]
  g : Array[Double] // nearest-neighbor distribution G(t)
}

///|
/// `sq.gr.ripley` — Ripley's K function with the variance-stabilizing L
/// transform, plus the nearest-neighbor distribution G. Coordinates use the
/// first two dimensions; the study area is the bounding box. `radius_max`
/// defaults to half the bounding-box diagonal; radii are `n_steps` uniform steps.
pub fn scanpy_ripley_k(
  coords : Array[Array[Double]],
  radius_max~ : Double = 0.0,
  n_steps~ : Int32 = 10,
) -> ScanpyRipleyK {
  let n = coords.length()
  let radii : Array[Double] = Array(capacity=n_steps)
  let k_vals : Array[Double] = Array(capacity=n_steps)
  let l_vals : Array[Double] = Array(capacity=n_steps)
  let g_vals : Array[Double] = Array(capacity=n_steps)
  if n < 2 || n_steps <= 0 {
    for _s in 0..<n_steps {
      radii.push(0.0)
      k_vals.push(0.0)
      l_vals.push(0.0)
      g_vals.push(0.0)
    }
    return ScanpyRipleyK::{ radii, k: k_vals, l: l_vals, g: g_vals, }
  }
  // bounding box over dims 0/1
  let mut min_x = coords[0][0]
  let mut max_x = coords[0][0]
  let mut min_y = coords[0][1]
  let mut max_y = coords[0][1]
  for i = 0; i < n; i = i + 1 {
    if coords[i][0] < min_x {
      min_x = coords[i][0]
    }
    if coords[i][0] > max_x {
      max_x = coords[i][0]
    }
    if coords[i][1] < min_y {
      min_y = coords[i][1]
    }
    if coords[i][1] > max_y {
      max_y = coords[i][1]
    }
  }
  let area = (max_x - min_x) * (max_y - min_y)
  let rmax = if radius_max > 0.0 {
    radius_max
  } else {
    let dx = max_x - min_x
    let dy = max_y - min_y
    (dx * dx + dy * dy).sqrt() / 2.0
  }
  // nearest-neighbor distances (for G)
  let nn = Array::make(n, @double.infinity)
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      if i != j {
        let dx = coords[i][0] - coords[j][0]
        let dy = coords[i][1] - coords[j][1]
        let d = (dx * dx + dy * dy).sqrt()
        if d < nn[i] {
          nn[i] = d
        }
      }
    }
  }
  // ordered pair counts per radius (single pass over radii)
  for s = 0; s < n_steps; s = s + 1 {
    let t = if n_steps == 1 {
      rmax
    } else {
      rmax * s.to_double() / (n_steps - 1).to_double()
    }
    let mut pairs = 0.0
    let mut g_count = 0
    for i = 0; i < n; i = i + 1 {
      if nn[i] <= t {
        g_count = g_count + 1
      }
      for j = 0; j < n; j = j + 1 {
        if i != j {
          let dx = coords[i][0] - coords[j][0]
          let dy = coords[i][1] - coords[j][1]
          if (dx * dx + dy * dy).sqrt() <= t {
            pairs = pairs + 1.0
          }
        }
      }
    }
    let k_val = if area > 0.0 && n > 1 {
      area / (n.to_double() * (n - 1).to_double()) * pairs
    } else {
      0.0
    }
    radii.push(t)
    k_vals.push(k_val)
    l_vals.push((k_val / @math.PI).sqrt() - t)
    g_vals.push(g_count.to_double() / n.to_double())
  }
  ScanpyRipleyK::{ radii, k: k_vals, l: l_vals, g: g_vals, }
}
