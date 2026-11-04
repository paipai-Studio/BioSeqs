/// Pairwise PERMANOVA with multiple testing correction.
///
/// scikit-bio equivalent: `stats.permanova_pairwise()`
///
/// Performs PERMANOVA for all pairs of groups, then applies
/// multiple testing correction (Bonferroni or Benjamini-Hochberg).

///|
/// Result of a single pairwise PERMANOVA test
pub struct SkbioPermanovaPairResult {
  group1 : String
  group2 : String
  statistic : Double
  p_value : Double
  corrected_p : Double
  n_permutations : Int32
}

///|
pub fn SkbioPermanovaPairResult::group1(
  self : SkbioPermanovaPairResult,
) -> String {
  self.group1
}

///|
pub fn SkbioPermanovaPairResult::group2(
  self : SkbioPermanovaPairResult,
) -> String {
  self.group2
}

///|
pub fn SkbioPermanovaPairResult::statistic(
  self : SkbioPermanovaPairResult,
) -> Double {
  self.statistic
}

///|
pub fn SkbioPermanovaPairResult::p_value(
  self : SkbioPermanovaPairResult,
) -> Double {
  self.p_value
}

///|
pub fn SkbioPermanovaPairResult::corrected_p(
  self : SkbioPermanovaPairResult,
) -> Double {
  self.corrected_p
}

///|
/// Compute pseudo-F statistic for two groups.
fn permanova_f_statistic(
  dm : Array[Array[Double]],
  group_labels : Array[Int32],
  n_group1 : Int32,
  n_group2 : Int32,
) -> Double {
  let n = dm.length()
  let total = n_group1 + n_group2
  // Total sum of squares
  let mut sst = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      sst = sst + dm[i][j] * dm[i][j]
    }
  }
  sst = sst / total.to_double()
  // Within-group sum of squares
  let mut ssw = 0.0
  // Group 1
  let mut ss1 = 0.0
  for i = 0; i < n; i = i + 1 {
    if group_labels[i] == 1 {
      for j = 0; j < n; j = j + 1 {
        if group_labels[j] == 1 {
          ss1 = ss1 + dm[i][j] * dm[i][j]
        }
      }
    }
  }
  if n_group1 > 0 {
    ss1 = ss1 / n_group1.to_double()
  }
  // Group 2
  let mut ss2 = 0.0
  for i = 0; i < n; i = i + 1 {
    if group_labels[i] == 2 {
      for j = 0; j < n; j = j + 1 {
        if group_labels[j] == 2 {
          ss2 = ss2 + dm[i][j] * dm[i][j]
        }
      }
    }
  }
  if n_group2 > 0 {
    ss2 = ss2 / n_group2.to_double()
  }
  ssw = ss1 + ss2
  let ssb = sst - ssw
  let df_between = 1.0 // k - 1 = 2 - 1
  let df_within = (total - 2).to_double()
  if df_within > 0.0 && ssw > 0.0 {
    ssb / df_between / (ssw / df_within)
  } else {
    0.0
  }
}

///|
/// Simple permutation test p-value for two groups
fn permanova_pvalue(
  dm : Array[Array[Double]],
  group_labels : Array[Int32],
  n_perm : Int32,
) -> (Double, Double) {
  let n = dm.length()
  let mut n1 : Int32 = 0
  for l in group_labels {
    if l == 1 {
      n1 = n1 + 1
    }
  }
  let mut n2 : Int32 = 0
  for l in group_labels {
    if l == 2 {
      n2 = n2 + 1
    }
  }
  let observed_f = permanova_f_statistic(dm, group_labels, n1, n2)
  let mut count = 0
  let perm = group_labels.copy()
  let mut seed = 12345
  for p = 0; p < n_perm; p = p + 1 {
    // Fisher-Yates shuffle
    for i = n - 1; i > 0; i = i - 1 {
      seed = (seed * 1103515245 + 12345) % 2147483647
      let j = seed % (i + 1)
      let tmp = perm[i]
      perm[i] = perm[j]
      perm[j] = tmp
    }
    let f = permanova_f_statistic(dm, perm, n1, n2)
    if f >= observed_f {
      count = count + 1
    }
  }
  let p = (count.to_double() + 1.0) / (n_perm.to_double() + 1.0)
  (observed_f, p)
}

///|
/// Benjamini-Hochberg correction
fn bh_correct(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return []
  }
  // Sort indices by p-value
  let indices : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    indices[i] = i
  }
  indices.sort_by(fn(a, b) {
    if p_values[a] < p_values[b] {
      -1
    } else if p_values[a] > p_values[b] {
      1
    } else {
      0
    }
  })
  let corrected : Array[Double] = Array::make(n, 0.0)
  let mut prev = 1.0
  for k = n - 1; k >= 0; k = k - 1 {
    let rank = k + 1
    let val = p_values[indices[k]] * n.to_double() / rank.to_double()
    let adj = if val < prev { val } else { prev }
    corrected[indices[k]] = if adj > 1.0 { 1.0 } else { adj }
    prev = adj
  }
  corrected
}

///|
/// Perform pairwise PERMANOVA across all groups.
///
/// dm: distance matrix
/// grouping: group label for each sample
/// n_permutations: number of permutations per test (default 999)
/// correction: "bh" for Benjamini-Hochberg, "bonferroni" for Bonferroni
pub fn skbio_permanova_pairwise(
  dm : Array[Array[Double]],
  grouping : Array[String],
  n_permutations~ : Int32 = 999,
  correction~ : String = "bh",
) -> Array[SkbioPermanovaPairResult] {
  // Get unique groups
  let group_set : Set[String] = Set([])
  for g in grouping {
    group_set.add(g)
  }
  let groups = group_set.to_array()
  groups.sort()
  // Generate all pairs
  let raw_results : Array[(String, String, Double, Double)] = Array(capacity=0)
  for i = 0; i < groups.length(); i = i + 1 {
    for j = i + 1; j < groups.length(); j = j + 1 {
      let g1 = groups[i]
      let g2 = groups[j]
      // Extract sub-matrix and labels
      let sub_indices : Array[Int32] = Array(capacity=0)
      let labels : Array[Int32] = Array(capacity=0)
      for k = 0; k < grouping.length(); k = k + 1 {
        if grouping[k] == g1 {
          sub_indices.push(k)
          labels.push(1)
        } else if grouping[k] == g2 {
          sub_indices.push(k)
          labels.push(2)
        }
      }
      let m = sub_indices.length()
      let sub_dm : Array[Array[Double]] = Array::make(m, [])
      for a = 0; a < m; a = a + 1 {
        sub_dm[a] = Array::make(m, 0.0)
        for b = 0; b < m; b = b + 1 {
          sub_dm[a][b] = dm[sub_indices[a]][sub_indices[b]]
        }
      }
      let (f, p) = permanova_pvalue(sub_dm, labels, n_permutations)
      raw_results.push((g1, g2, f, p))
    }
  }
  // Extract p-values for correction
  let p_vals : Array[Double] = Array(capacity=0)
  for r in raw_results {
    p_vals.push(r.3)
  }
  let corrected = match correction {
    "bonferroni" => {
      let n_tests = p_vals.length().to_double()
      p_vals.map(fn(p) {
        let v = p * n_tests
        if v > 1.0 {
          1.0
        } else {
          v
        }
      })
    }
    _ => bh_correct(p_vals)
  }
  // Build results
  let results : Array[SkbioPermanovaPairResult] = Array(capacity=0)
  for i = 0; i < raw_results.length(); i = i + 1 {
    let r = raw_results[i]
    results.push({
      group1: r.0,
      group2: r.1,
      statistic: r.2,
      p_value: r.3,
      corrected_p: corrected[i],
      n_permutations,
    })
  }
  results
}
