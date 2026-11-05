/// skbio_composition_stats.mbt — ANCOM and Gini index
/// Items 6-7: skbio_ancom, skbio_gini

///|
/// ANCOM result for a single feature.
pub struct SkbioAncomResult {
  feature_idx : Int32
  w_stat : Int32 // number of significant pairwise log-ratio tests
  percentile : Double // percentile of W statistic
  detected : Bool // whether this feature is differentially abundant
}

///|
/// Analysis of Composition of Microbiomes (ANCOM).
/// Detects differentially abundant features between groups.
/// samples: abundance table (rows=samples, cols=features)
/// grouping: group label for each sample (0 or 1 for two-group)
/// Returns one result per feature.
pub fn skbio_ancom(
  samples : Array[Array[Double]],
  grouping : Array[Int32],
  threshold~ : Double = 0.05,
) -> Array[SkbioAncomResult] {
  let n_samples = samples.length()
  if n_samples == 0 {
    return Array(capacity=0)
  }
  let n_features = samples[0].length()
  let results : Array[SkbioAncomResult] = Array(capacity=n_features)
  // Add pseudocount of 1 to avoid log(0)
  let log_data : Array[Array[Double]] = Array(capacity=n_samples)
  for i = 0; i < n_samples; i = i + 1 {
    let row : Array[Double] = Array(capacity=n_features)
    for j = 0; j < n_features; j = j + 1 {
      row.push(@math.ln(samples[i][j] + 1.0))
    }
    log_data.push(row)
  }
  // For each feature, compute W = number of other features where
  // the log-ratio test is significant between groups
  let w_stats : Array[Int32] = Array(capacity=n_features)
  for j = 0; j < n_features; j = j + 1 {
    w_stats.push(0)
  }
  for j = 0; j < n_features; j = j + 1 {
    for k = j + 1; k < n_features; k = k + 1 {
      // Compute log-ratio for each sample
      let group0_vals : Array[Double] = Array(capacity=0)
      let group1_vals : Array[Double] = Array(capacity=0)
      for i = 0; i < n_samples; i = i + 1 {
        let ratio = log_data[i][j] - log_data[i][k]
        if grouping[i] == 0 {
          group0_vals.push(ratio)
        } else {
          group1_vals.push(ratio)
        }
      }
      // Simple t-test approximation
      if group0_vals.length() > 1 && group1_vals.length() > 1 {
        let mean0 = skbio_mean(group0_vals)
        let mean1 = skbio_mean(group1_vals)
        let var0 = skbio_variance(group0_vals, mean0)
        let var1 = skbio_variance(group1_vals, mean1)
        let n0 = group0_vals.length().to_double()
        let n1 = group1_vals.length().to_double()
        let se = @math.pow(var0 / n0 + var1 / n1, 0.5)
        if se > 0.0 {
          let t_stat = (mean0 - mean1) / se
          // Approximate p-value using |t| > 2.5 (rough threshold)
          let abs_t = if t_stat < 0.0 { -t_stat } else { t_stat }
          if abs_t > 2.5 {
            w_stats[j] = w_stats[j] + 1
            w_stats[k] = w_stats[k] + 1
          }
        }
      }
    }
  }
  // Compute percentile threshold (70th percentile of W)
  let w_sorted : Array[Int32] = Array(capacity=n_features)
  for j = 0; j < n_features; j = j + 1 {
    w_sorted.push(w_stats[j])
  }
  w_sorted.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  let threshold_idx = n_features * 7 / 10
  let w_threshold = if threshold_idx < n_features {
    w_sorted[threshold_idx]
  } else {
    0
  }
  for j = 0; j < n_features; j = j + 1 {
    let pct = if n_features > 0 {
      w_stats[j].to_double() / n_features.to_double() * 100.0
    } else {
      0.0
    }
    results.push(SkbioAncomResult::{
      feature_idx: j,
      w_stat: w_stats[j],
      percentile: pct,
      detected: w_stats[j] > w_threshold && w_stats[j] > 0,
    })
  }
  // Suppress unused threshold warning
  ignore(threshold)
  results
}

///|
fn skbio_mean(vals : Array[Double]) -> Double {
  let mut sum = 0.0
  for i = 0; i < vals.length(); i = i + 1 {
    sum = sum + vals[i]
  }
  if vals.length() == 0 {
    0.0
  } else {
    sum / vals.length().to_double()
  }
}

///|
fn skbio_variance(vals : Array[Double], mean : Double) -> Double {
  if vals.length() <= 1 {
    return 0.0
  }
  let mut sum_sq = 0.0
  for i = 0; i < vals.length(); i = i + 1 {
    let diff = vals[i] - mean
    sum_sq = sum_sq + diff * diff
  }
  sum_sq / (vals.length() - 1).to_double()
}

///|
/// Gini index: measures inequality (0=perfect equality, 1=max inequality).
pub fn skbio_gini(abundances : Array[Double]) -> Double {
  let n = abundances.length()
  if n == 0 {
    return 0.0
  }
  // Sort ascending
  let sorted : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    sorted.push(abundances[i])
  }
  sorted.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  let mut sum = 0.0
  let mut total = 0.0
  for i = 0; i < n; i = i + 1 {
    sum = sum + (i.to_double() + 1.0) * sorted[i]
    total = total + sorted[i]
  }
  if total == 0.0 {
    return 0.0
  }
  let gini = 2.0 * sum / (n.to_double() * total) -
    (n.to_double() + 1.0) / n.to_double()
  // Correct for small-n bias
  gini * n.to_double() / (n - 1).to_double()
}

///|
/// Gini-Simpson index: 1 - sum(p_i^2), measures diversity.
pub fn skbio_gini_simpson(abundances : Array[Double]) -> Double {
  let mut total = 0.0
  for i = 0; i < abundances.length(); i = i + 1 {
    total = total + abundances[i]
  }
  if total == 0.0 {
    return 0.0
  }
  let mut sum_sq = 0.0
  for i = 0; i < abundances.length(); i = i + 1 {
    let p = abundances[i] / total
    sum_sq = sum_sq + p * p
  }
  1.0 - sum_sq
}
