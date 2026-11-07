// Scikit-bio third-round: rarefaction and diversity curves (items 24-26)
// rarefaction_single, rarefaction_double, pd_rarefaction

///|
pub fn skbio3_rarefaction_single(
  counts : Array[Int32],
  max_depth : Int32,
  step : Int32,
) -> Array[(Int32, Double)] {
  // Single rarefaction: expected number of OTUs at each sampling depth
  // E[S] = sum_i (1 - C(N, n) * C(N - n_i, n) / C(N, n))
  // where N = total, n = sample depth, n_i = count of OTU i
  let n_otus = counts.length()
  let mut total = 0
  for i in 0..<n_otus {
    total = total + counts[i]
  }
  let result : Array[(Int32, Double)] = []
  let mut depth = step
  while depth <= max_depth && depth <= total {
    if depth > total {
      break
    }
    // Expected OTU count at this depth
    let mut expected_s = 0.0
    for i in 0..<n_otus {
      let ni = counts[i]
      if ni == 0 {
        continue
      }
      // P(OTU i is present) = 1 - C(N-ni, depth) / C(N, depth)
      // Use log to avoid overflow: log(C(N-ni, d) / C(N, d))
      let log_ratio = skbio3_log_choose(total - ni, depth) -
        skbio3_log_choose(total, depth)
      let p_absent = @math.exp(log_ratio)
      expected_s = expected_s + (1.0 - p_absent)
    }
    result.push((depth, expected_s))
    depth = depth + step
  }
  result
}

///|
fn skbio3_log_choose(n : Int32, k : Int32) -> Double {
  // log(C(n, k)) = log(n!) - log(k!) - log((n-k)!)
  if k < 0 || k > n {
    return 0.0 // C(n,k) = 0, log(0) = -inf, but we return 0 for practical purposes
  }
  if k == 0 || k == n {
    return 0.0 // C(n,0) = C(n,n) = 1, log(1) = 0
  }
  // Use smaller k for efficiency: C(n,k) = C(n, n-k)
  let kk = if k > n - k { n - k } else { k }
  let mut result = 0.0
  for i in 0..<kk {
    result = result +
      @math.ln((n - i).to_double()) -
      @math.ln((i + 1).to_double())
  }
  result
}

///|
pub fn skbio3_rarefaction_double(
  counts1 : Array[Int32],
  counts2 : Array[Int32],
  max_depth : Int32,
  step : Int32,
) -> Array[(Int32, Double)] {
  // Double rarefaction: expected shared OTUs between two samples at given depth
  let result : Array[(Int32, Double)] = []
  let mut total1 = 0
  for i in 0..<counts1.length() {
    total1 = total1 + counts1[i]
  }
  let mut total2 = 0
  for i in 0..<counts2.length() {
    total2 = total2 + counts2[i]
  }
  // Get union of OTU IDs (use index as ID)
  let n_otus = if counts1.length() > counts2.length() {
    counts1.length()
  } else {
    counts2.length()
  }
  let mut depth = step
  while depth <= max_depth && depth <= total1 && depth <= total2 {
    let mut expected_shared = 0.0
    for i in 0..<n_otus {
      let n1i = if i < counts1.length() { counts1[i] } else { 0 }
      let n2i = if i < counts2.length() { counts2[i] } else { 0 }
      if n1i == 0 || n2i == 0 {
        continue
      }
      // P(OTU i present in both) = P(present in 1) * P(present in 2)
      let p1_absent = @math.exp(
        skbio3_log_choose(total1 - n1i, depth) -
        skbio3_log_choose(total1, depth),
      )
      let p2_absent = @math.exp(
        skbio3_log_choose(total2 - n2i, depth) -
        skbio3_log_choose(total2, depth),
      )
      let p_both_present = (1.0 - p1_absent) * (1.0 - p2_absent)
      expected_shared = expected_shared + p_both_present
    }
    result.push((depth, expected_shared))
    depth = depth + step
  }
  result
}

///|
pub fn skbio3_pd_rarefaction(
  counts : Array[Int32],
  branch_lengths : Array[Double],
  max_depth : Int32,
  step : Int32,
) -> Array[(Int32, Double)] {
  // Phylogenetic rarefaction: expected Faith PD at each sampling depth
  // E[PD] = sum_b (1 - P(branch b not sampled))
  // P(branch b not sampled) = C(N - n_b, depth) / C(N, depth)
  // where n_b = total count of tips descending from branch b
  let n_branches = branch_lengths.length()
  let mut total = 0
  for i in 0..<counts.length() {
    total = total + counts[i]
  }
  let result : Array[(Int32, Double)] = []
  let mut depth = step
  while depth <= max_depth && depth <= total {
    let mut expected_pd = 0.0
    for b in 0..<n_branches {
      let bl = branch_lengths[b]
      // n_b = number of tips descending from this branch
      // Simplified: use the count at index b (assumes 1:1 mapping)
      let n_b = if b < counts.length() { counts[b] } else { 0 }
      if n_b == 0 {
        continue
      }
      // P(branch present) = 1 - C(N - n_b, depth) / C(N, depth)
      let log_ratio = skbio3_log_choose(total - n_b, depth) -
        skbio3_log_choose(total, depth)
      let p_absent = @math.exp(log_ratio)
      expected_pd = expected_pd + bl * (1.0 - p_absent)
    }
    result.push((depth, expected_pd))
    depth = depth + step
  }
  result
}
