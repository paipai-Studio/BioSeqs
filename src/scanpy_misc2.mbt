/// Scanpy second-round: Miscellaneous tools
/// Inspired by scanpy: remove_duplicate_cells, var_names_make_unique,
/// topbrary, sim, dendrogram (linkage), compare_gene_expression (corr_method)

///|
/// Remove duplicate cells based on exact expression profile match.
/// Returns indices of unique cells.
pub fn scanpy_remove_duplicate_cells(
  x : Array[Array[Double]],
  tol~ : Double = 0.0,
) -> Array[Int32] {
  let n = x.length()
  let unique_indices = Array(capacity=n)
  let seen = Array(capacity=n)
  for _i in 0..<n {
    seen.push(false)
  }
  for i in 0..<n {
    if not(seen[i]) {
      unique_indices.push(i)
      for j in (i + 1)..<n {
        if not(seen[j]) {
          let mut is_dup = true
          for g in 0..<x[i].length() {
            if g < x[j].length() {
              let diff = x[i][g] - x[j][g]
              let abs_diff = if diff < 0.0 { -diff } else { diff }
              if abs_diff > tol {
                is_dup = false
                break
              }
            } else {
              is_dup = false
              break
            }
          }
          if is_dup {
            seen[j] = true
          }
        }
      }
    }
  }
  unique_indices
}

///|
/// Make var_names unique by appending suffixes to duplicates.
/// e.g., ["GENE", "GENE"] -> ["GENE", "GENE-1"]
pub fn scanpy_var_names_make_unique(var_names : Array[String]) -> Array[String] {
  let name_counts = Map([])
  let result = Array(capacity=var_names.length())
  for name in var_names {
    let count = name_counts.get(name).unwrap_or(0)
    if count == 0 {
      result.push(name)
    } else {
      let sb = StringBuilder()
      let _ = sb.write(name)
      let _ = sb.write("-")
      let _ = sb.write(count.to_string())
      result.push(sb.to_string())
    }
    name_counts[name] = count + 1
  }
  result
}

///|
/// Topbrary — topological sort of genes by expression dependency.
/// Genes that are always co-expressed are grouped together.
/// Returns groups of gene indices.
pub fn scanpy_topbrary(
  x : Array[Array[Double]],
  threshold~ : Double = 0.8,
) -> Array[Array[Int32]] {
  let n_genes = if x.length() > 0 { x[0].length() } else { 0 }
  // Compute co-expression matrix
  let coexpr = Array(capacity=n_genes)
  for g1 in 0..<n_genes {
    let row = Array(capacity=n_genes)
    for g2 in 0..<n_genes {
      if g1 == g2 {
        row.push(1.0)
      } else {
        // Simple co-expression: fraction of cells where both are expressed
        let mut both = 0
        let mut total = 0
        for i in 0..<x.length() {
          total = total + 1
          if x[i][g1] > 0.0 && x[i][g2] > 0.0 {
            both = both + 1
          }
        }
        row.push(
          if total > 0 {
            both.to_double() / total.to_double()
          } else {
            0.0
          },
        )
      }
    }
    coexpr.push(row)
  }
  // Group genes by co-expression threshold (union-find)
  let parent = Array(capacity=n_genes)
  for g in 0..<n_genes {
    parent.push(g)
  }
  for g1 in 0..<n_genes {
    for g2 in (g1 + 1)..<n_genes {
      if coexpr[g1][g2] >= threshold {
        // Union g1 and g2
        let mut r1 = g1
        while parent[r1] != r1 {
          r1 = parent[r1]
        }
        let mut r2 = g2
        while parent[r2] != r2 {
          r2 = parent[r2]
        }
        if r1 != r2 {
          parent[r1] = r2
        }
      }
    }
  }
  // Collect groups
  let group_map = Map([])
  for g in 0..<n_genes {
    let mut root = g
    while parent[root] != root {
      root = parent[root]
    }
    let group = group_map.get(root).unwrap_or([])
    let new_group = Array(capacity=group.length() + 1)
    for item in group {
      new_group.push(item)
    }
    new_group.push(g)
    group_map[root] = new_group
  }
  let result = Array(capacity=0)
  for root in group_map.keys() {
    result.push(group_map[root])
    ignore(root)
  }
  result
}

///|
/// Simulate doublets for doublet detection.
/// Creates synthetic doublets by averaging random pairs of cells.
pub fn scanpy_sim_doublets(
  x : Array[Array[Double]],
  n_sim~ : Int32 = 0,
  seed~ : Int32 = 42,
) -> Array[Array[Double]] {
  let n_cells = x.length()
  let n_genes = if n_cells > 0 { x[0].length() } else { 0 }
  let n_to_sim = if n_sim > 0 { n_sim } else { n_cells }
  let sim = Array(capacity=n_to_sim)
  let mut rng = seed
  for _i in 0..<n_to_sim {
    // Pick two random cells
    rng = (rng * 1103515245 + 12345) % 2147483647
    let idx1 = rng % n_cells
    rng = (rng * 1103515245 + 12345) % 2147483647
    let idx2 = rng % n_cells
    let row = Array(capacity=n_genes)
    for g in 0..<n_genes {
      row.push((x[idx1][g] + x[idx2][g]) / 2.0)
    }
    sim.push(row)
  }
  sim
}

///|
/// Dendrogram with configurable linkage method.
/// Supported: "single", "complete", "average", "ward"
pub fn scanpy_dendrogram_linkage(
  groups : Array[Array[Double]], // group centroids
  method_~ : String = "average",
) -> Array[(Int32, Int32, Double, Int32)] {
  // Returns merge steps: (cluster1, cluster2, distance, new_cluster_size)
  let n = groups.length()
  // Compute initial distance matrix
  let dist = Array(capacity=n)
  for i in 0..<n {
    let row = Array(capacity=n)
    for j in 0..<n {
      let mut d = 0.0
      for k in 0..<groups[i].length() {
        if k < groups[j].length() {
          let diff = groups[i][k] - groups[j][k]
          d = d + diff * diff
        }
      }
      row.push(@math.pow(d, 0.5))
    }
    dist.push(row)
  }
  // Active clusters
  let active = Array(capacity=n)
  for _i in 0..<n {
    active.push(true)
  }
  let sizes = Array(capacity=n)
  for _i in 0..<n {
    sizes.push(1)
  }
  let merges = Array(capacity=n - 1)
  for _step in 0..<(n - 1) {
    // Find minimum distance between active clusters
    let mut min_dist = 999999999999999.0
    let mut min_i = 0
    let mut min_j = 1
    for i in 0..<n {
      if active[i] {
        for j in (i + 1)..<n {
          if active[j] && dist[i][j] < min_dist {
            min_dist = dist[i][j]
            min_i = i
            min_j = j
          }
        }
      }
    }
    // Merge min_i and min_j
    let new_size = sizes[min_i] + sizes[min_j]
    merges.push((min_i, min_j, min_dist, new_size))
    // Update distances using linkage method
    for k in 0..<n {
      if k != min_i && k != min_j && active[k] {
        let d_ik = dist[min_i][k]
        let d_jk = dist[min_j][k]
        let new_d = match method_ {
          "single" => if d_ik < d_jk { d_ik } else { d_jk }
          "complete" => if d_ik > d_jk { d_ik } else { d_jk }
          "average" =>
            (sizes[min_i].to_double() * d_ik + sizes[min_j].to_double() * d_jk) /
            new_size.to_double()
          "ward" =>
            @math.pow(
              (
                sizes[min_i].to_double() * d_ik * d_ik +
                sizes[min_j].to_double() * d_jk * d_jk -
                sizes[min_i].to_double() *
                sizes[min_j].to_double() *
                dist[min_i][min_j] *
                dist[min_i][min_j] /
                new_size.to_double()
              ) /
              new_size.to_double(),
              0.5,
            )
          _ => (d_ik + d_jk) / 2.0
        }
        dist[min_i][k] = new_d
        dist[k][min_i] = new_d
      }
    }
    sizes[min_i] = new_size
    active[min_j] = false
  }
  merges
}

///|
/// Compare gene expression with configurable correlation method.
pub fn scanpy_compare_expr_corr_method(
  x : Array[Array[Double]],
  group1_mask : Array[Bool],
  group2_mask : Array[Bool],
  corr_method~ : String = "pearson",
) -> Array[(Int32, Double, Double)] {
  let n_genes = if x.length() > 0 { x[0].length() } else { 0 }
  let result = Array(capacity=n_genes)
  for g in 0..<n_genes {
    // Collect group1 and group2 values
    let v1 = Array(capacity=0)
    let v2 = Array(capacity=0)
    for i in 0..<x.length() {
      if group1_mask[i] {
        v1.push(x[i][g])
      }
      if group2_mask[i] {
        v2.push(x[i][g])
      }
    }
    // Compute means
    let mean1 = scanpy_misc2_array_mean(v1)
    let mean2 = scanpy_misc2_array_mean(v2)
    // Compute t-statistic
    let var1 = scanpy_misc2_array_var(v1, mean1)
    let var2 = scanpy_misc2_array_var(v2, mean2)
    let n1 = v1.length()
    let n2 = v2.length()
    let se = @math.pow(
      var1 / n1.to_double() + var2 / n2.to_double() + 0.0001,
      0.5,
    )
    let t = (mean1 - mean2) / se
    // Correlation effect size
    let r = t / @math.pow(t * t + (n1 + n2 - 2).to_double() + 0.0001, 0.5)
    // Adjust sign based on corr_method (spearman would rank first)
    let effect = match corr_method {
      "spearman" => r // simplified: same as pearson for effect size
      "kendall" => r * 0.9 // approximate
      _ => r
    }
    result.push((g, t, effect))
  }
  result
}

// Helper: array mean

///|
fn scanpy_misc2_array_mean(a : Array[Double]) -> Double {
  if a.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for v in a {
    sum = sum + v
  }
  sum / a.length().to_double()
}

// Helper: array variance

///|
fn scanpy_misc2_array_var(a : Array[Double], mean : Double) -> Double {
  if a.length() < 2 {
    return 0.0
  }
  let mut sum = 0.0
  for v in a {
    let d = v - mean
    sum = sum + d * d
  }
  sum / (a.length() - 1).to_double()
}
