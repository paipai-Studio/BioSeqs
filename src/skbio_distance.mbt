///|
/// scikit-bio distance-based statistical tests.
///
/// Port of `skbio.stats.distance.anosim` (Analysis of Similarities) and
/// `skbio.stats.distance.permanova` (PERMANOVA). Both compute a deterministic
/// test statistic from a distance matrix and a group assignment, plus a
/// permutation-based p-value. The permutation uses a deterministic LCG +
/// Fisher-Yates shuffle, so the statistic matches scikit-bio exactly while the
/// p-value is a valid permutation p-value (but not bit-identical to numpy's
/// RNG).

///|
pub struct AnosimResult {
  statistic : Double
  p_value : Double
  n : Int
  num_groups : Int
}

///|
pub struct PermanovaResult {
  statistic : Double
  p_value : Double
  n : Int
  num_groups : Int
}

///|
fn skbd_nan() -> Double {
  0.0 / 0.0
}

///|
/// Condensed (upper-triangle off-diagonal) form of a square distance matrix.
fn skbd_condensed(matrix : Array[Array[Double]]) -> Array[Double] {
  let n = matrix.length()
  let out : Array[Double] = Array::new()
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      out.push(matrix[i][j])
    }
  }
  out
}

///|
/// Number of distinct group labels.
fn skbd_num_groups(grouping : Array[Int]) -> Int {
  let seen : Array[Int] = Array::new()
  for g in grouping {
    let mut found = false
    for s in seen {
      if s == g {
        found = true
        break
      }
    }
    if !found {
      seen.push(g)
    }
  }
  seen.length()
}

///|
/// Average-rank of a vector (scipy.stats.rankdata with method="average").
fn skbd_rankdata_average(values : Array[Double]) -> Array[Double] {
  let n = values.length()
  let ranks = Array::make(n, 0.0)
  let idx = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    idx[i] = i
  }
  // sort indices by value (insertion sort)
  for i = 1; i < n; i = i + 1 {
    let key = values[idx[i]]
    let ki = idx[i]
    let mut j = i - 1
    while j >= 0 && values[idx[j]] > key {
      idx[j + 1] = idx[j]
      j = j - 1
    }
    idx[j + 1] = ki
  }
  let mut i = 0
  while i < n {
    let mut j = i
    while j + 1 < n && values[idx[j + 1]] == values[idx[i]] {
      j = j + 1
    }
    let avg = (i + j).to_double() / 2.0 + 1.0
    let mut k = i
    while k <= j {
      ranks[idx[k]] = avg
      k = k + 1
    }
    i = j + 1
  }
  ranks
}

///|
fn skbd_lcg(seed : Int) -> Int {
  (seed * 1103515245 + 12345) & 0x7FFFFFFF
}

///|
/// ANOSIM R statistic from pre-ranked condensed distances + grouping.
fn skbd_anosim_r(ranks : Array[Double], grouping : Array[Int]) -> Double {
  let n = grouping.length()
  let mut r_w_sum = 0.0
  let mut r_w_n = 0
  let mut r_b_sum = 0.0
  let mut r_b_n = 0
  let mut p = 0
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if grouping[i] == grouping[j] {
        r_w_sum = r_w_sum + ranks[p]
        r_w_n = r_w_n + 1
      } else {
        r_b_sum = r_b_sum + ranks[p]
        r_b_n = r_b_n + 1
      }
      p = p + 1
    }
  }
  let r_w = r_w_sum / r_w_n.to_double()
  let r_b = r_b_sum / r_b_n.to_double()
  (r_b - r_w) / (n.to_double() * (n - 1).to_double() / 4.0)
}

///|
/// Analysis of Similarities (ANOSIM).
pub fn skbio_anosim(
  dm : Array[Array[Double]],
  grouping : Array[Int],
  permutations : Int,
  seed : Int,
) -> AnosimResult {
  let n = dm.length()
  let ranks = skbd_rankdata_average(skbd_condensed(dm))
  let obs = skbd_anosim_r(ranks, grouping)
  let mut count_better = 0
  let mut rng = seed
  for _p = 0; _p < permutations; _p = _p + 1 {
    let g = grouping.copy()
    let mut i = n - 1
    while i > 0 {
      rng = skbd_lcg(rng)
      let j = rng % (i + 1)
      let t = g[i]
      g[i] = g[j]
      g[j] = t
      i = i - 1
    }
    if skbd_anosim_r(ranks, g) >= obs {
      count_better = count_better + 1
    }
  }
  let p_value = if permutations == 0 {
    skbd_nan()
  } else {
    (count_better + 1).to_double() / (permutations + 1).to_double()
  }
  AnosimResult::{ statistic: obs, p_value, n, num_groups: skbd_num_groups(grouping) }
}

///|
/// PERMANOVA within-group sum of squares (upper-triangle convention). A fixed
/// factor of 2 cancels in the final F ratio, matching scikit-bio.
fn skbd_permanova_s_w(dm : Array[Array[Double]], grouping : Array[Int]) -> Double {
  let n = dm.length()
  let mut s_w = 0.0
  let seen : Array[Int] = Array::new()
  for g in grouping {
    let mut dup = false
    for s in seen {
      if s == g {
        dup = true
        break
      }
    }
    if dup {
      continue
    }
    seen.push(g)
    // count size and accumulate within-group sum of squared distances
    let mut size = 0
    for h in grouping {
      if h == g {
        size = size + 1
      }
    }
    let mut sum = 0.0
    for i = 0; i < n; i = i + 1 {
      for j = i + 1; j < n; j = j + 1 {
        if grouping[i] == g && grouping[j] == g {
          sum = sum + dm[i][j] * dm[i][j]
        }
      }
    }
    if size > 0 {
      s_w = s_w + sum / size.to_double()
    }
  }
  s_w
}

///|
/// PERMANOVA pseudo-F statistic from a distance matrix + grouping.
fn skbd_permanova_f(dm : Array[Array[Double]], grouping : Array[Int]) -> Double {
  let n = dm.length()
  let g = skbd_num_groups(grouping)
  let mut s_t = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      s_t = s_t + dm[i][j] * dm[i][j]
    }
  }
  s_t = s_t / n.to_double()
  let s_w = skbd_permanova_s_w(dm, grouping)
  let s_a = s_t - s_w
  let f = (s_a / (g - 1).to_double()) / (s_w / (n - g).to_double())
  f
}

///|
/// PERMANOVA.
pub fn skbio_permanova(
  dm : Array[Array[Double]],
  grouping : Array[Int],
  permutations : Int,
  seed : Int,
) -> PermanovaResult {
  let n = dm.length()
  let obs = skbd_permanova_f(dm, grouping)
  let g = skbd_num_groups(grouping)
  let mut count_better = 0
  let mut rng = seed
  for _p = 0; _p < permutations; _p = _p + 1 {
    let grp = grouping.copy()
    let mut i = n - 1
    while i > 0 {
      rng = skbd_lcg(rng)
      let j = rng % (i + 1)
      let t = grp[i]
      grp[i] = grp[j]
      grp[j] = t
      i = i - 1
    }
    let f_perm = skbd_permanova_f(dm, grp)
    if f_perm >= obs {
      count_better = count_better + 1
    }
  }
  let p_value = if permutations == 0 {
    skbd_nan()
  } else {
    (count_better + 1).to_double() / (permutations + 1).to_double()
  }
  PermanovaResult::{ statistic: obs, p_value, n, num_groups: g }
}
///|
pub struct PermdispResult {
  statistic : Double
  p_value : Double
  n : Int
  num_groups : Int
}

///|
/// Grouped distance-to-centroid arrays used by the one-way ANOVA.
fn skbd_permdisp_distances(
  pc : Array[Array[Double]],
  grouping : Array[Int],
) -> Array[Array[Double]] {
  let n = pc.length()
  let d = if n > 0 { pc[0].length() } else { 0 }
  let seen : Array[Int] = Array::new()
  let result : Array[Array[Double]] = Array::new()
  for g in grouping {
    let mut dup = false
    for s in seen {
      if s == g {
        dup = true
        break
      }
    }
    if dup {
      continue
    }
    seen.push(g)
    let members : Array[Int] = Array::new()
    for i = 0; i < n; i = i + 1 {
      if grouping[i] == g {
        members.push(i)
      }
    }
    let ng = members.length()
    let centroid = Array::make(d, 0.0)
    for j = 0; j < d; j = j + 1 {
      let mut s = 0.0
      for m in members {
        s = s + pc[m][j]
      }
      centroid[j] = s / ng.to_double()
    }
    let dists : Array[Double] = Array::new()
    for m in members {
      let mut sq = 0.0
      for j = 0; j < d; j = j + 1 {
        let diff = pc[m][j] - centroid[j]
        sq = sq + diff * diff
      }
      dists.push(sq.sqrt())
    }
    result.push(dists)
  }
  result
}

///|
/// One-way ANOVA F-statistic over the grouped distances.
fn skbd_one_way_anova_f(dist_grouped : Array[Array[Double]]) -> Double {
  let a = dist_grouped.length()
  let mut grand_sum = 0.0
  let mut grand_n = 0
  for g in dist_grouped {
    for v in g {
      grand_sum = grand_sum + v
      grand_n = grand_n + 1
    }
  }
  let grand_mean = grand_sum / grand_n.to_double()
  let mut ss_between = 0.0
  let mut ss_within = 0.0
  for g in dist_grouped {
    let ng = g.length()
    let mut s = 0.0
    for v in g {
      s = s + v
    }
    let mean_g = s / ng.to_double()
    ss_between = ss_between + ng.to_double() * (mean_g - grand_mean) * (mean_g - grand_mean)
    for v in g {
      ss_within = ss_within + (v - mean_g) * (v - mean_g)
    }
  }
  let df_between = a - 1
  let df_within = grand_n - a
  (ss_between / df_between.to_double()) / (ss_within / df_within.to_double())
}

///|
fn skbd_permdisp_f(pc : Array[Array[Double]], grouping : Array[Int]) -> Double {
  skbd_one_way_anova_f(skbd_permdisp_distances(pc, grouping))
}

///|
/// PERMDISP (permutational multivariate dispersion, centroid test).
pub fn skbio_permdisp(
  dm : Array[Array[Double]],
  grouping : Array[Int],
  permutations : Int,
  seed : Int,
) -> PermdispResult {
  let n = dm.length()
  // principal coordinates via PCoA (Gower centering + eigen)
  let b = skbio_center_distance_matrix(dm)
  let (raw_evals, raw_evecs) = eigh_symmetric(b, n)
  let (evals, evecs) = linalg_sort_eigen(raw_evals, raw_evecs, n)
  let mut d = 0
  while d < n && evals[d] > 1.0e-12 {
    d = d + 1
  }
  let pc : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(d, 0.0)
    for j = 0; j < d; j = j + 1 {
      row[j] = evecs[i][j] * evals[j].sqrt()
    }
    pc.push(row)
  }
  let obs = skbd_permdisp_f(pc, grouping)
  let mut count_better = 0
  let mut rng = seed
  for _p = 0; _p < permutations; _p = _p + 1 {
    let g = grouping.copy()
    let mut i = n - 1
    while i > 0 {
      rng = skbd_lcg(rng)
      let j = rng % (i + 1)
      let t = g[i]
      g[i] = g[j]
      g[j] = t
      i = i - 1
    }
    if skbd_permdisp_f(pc, g) >= obs {
      count_better = count_better + 1
    }
  }
  let p_value = if permutations == 0 {
    skbd_nan()
  } else {
    (count_better + 1).to_double() / (permutations + 1).to_double()
  }
  PermdispResult::{ statistic: obs, p_value, n, num_groups: skbd_num_groups(grouping) }
}
