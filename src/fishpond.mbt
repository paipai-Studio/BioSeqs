///|
/// fishpond - Swish non-parametric differential expression analysis.
///
/// References:
///   - Bioconductor fishpond package (Swish method)
///   - Zhu A, Srivastava A, Ibrahim JG, Patro R, Love MI (2019) Nonparametric
///     expression analysis using inferential replicate counts. Nucleic Acids
///     Research 47(18):e105.
///
/// Swish (SAMseq With Inferential Samples Handles) performs non-parametric
/// testing for differential transcript/gene expression. It uses the
/// Mann-Whitney-Wilcoxon rank-sum statistic combined with permutation-based
/// inference, which avoids distributional assumptions about count data.
///
/// All public identifiers are prefixed with `fish_` to avoid conflicts.

// ============================================================================
// Data structures
// ============================================================================

///|
/// Expression counts for a set of transcripts across samples.
/// rows = transcripts, cols = samples.
pub struct FishCounts {
  transcript_names : Array[String]
  sample_names : Array[String]
  condition : Array[String] // condition label per sample (e.g. "case"/"control")
  counts : Array[Array[Double]] // rows=transcripts, cols=samples
}

///|
/// Create a counts object.
pub fn FishCounts::new(
  transcript_names : Array[String],
  sample_names : Array[String],
  condition : Array[String],
  counts : Array[Array[Double]],
) -> FishCounts {
  { transcript_names, sample_names, condition, counts }
}

///|
/// Get transcript names.
pub fn FishCounts::transcript_names(self : FishCounts) -> Array[String] {
  self.transcript_names
}

///|
/// Get sample names.
pub fn FishCounts::sample_names(self : FishCounts) -> Array[String] {
  self.sample_names
}

///|
/// Get condition labels.
pub fn FishCounts::condition(self : FishCounts) -> Array[String] {
  self.condition
}

///|
/// Get counts matrix.
pub fn FishCounts::counts(self : FishCounts) -> Array[Array[Double]] {
  self.counts
}

///|
/// Differential expression result for a single transcript.
pub struct FishResult {
  transcript : String
  statistic : Double // Mann-Whitney-Wilcoxon statistic
  log2_fold_change : Double // mean(case) - mean(control) in log2
  p_value : Double // permutation-based p-value
  fdr : Double // BH-adjusted p-value
  direction : String // "up", "down", or "ns"
}

///|
/// Create a result.
pub fn FishResult::new(
  transcript : String,
  statistic : Double,
  log2_fold_change : Double,
  p_value : Double,
  fdr : Double,
  direction : String,
) -> FishResult {
  { transcript, statistic, log2_fold_change, p_value, fdr, direction }
}

///|
/// Get transcript name.
pub fn FishResult::transcript(self : FishResult) -> String {
  self.transcript
}

///|
/// Get test statistic.
pub fn FishResult::statistic(self : FishResult) -> Double {
  self.statistic
}

///|
/// Get log2 fold change.
pub fn FishResult::log2_fold_change(self : FishResult) -> Double {
  self.log2_fold_change
}

///|
/// Get p-value.
pub fn FishResult::p_value(self : FishResult) -> Double {
  self.p_value
}

///|
/// Get FDR.
pub fn FishResult::fdr(self : FishResult) -> Double {
  self.fdr
}

///|
/// Get direction ("up", "down", or "ns").
pub fn FishResult::direction(self : FishResult) -> String {
  self.direction
}

// ============================================================================
// Core algorithm
// ============================================================================

///|
/// Compute the Mann-Whitney-Wilcoxon rank-sum statistic for a single
/// transcript, comparing case vs control samples.
/// Returns the W statistic (rank sum of the first group).
pub fn fish_mann_whitney_test(
  case_vals : Array[Double],
  control_vals : Array[Double],
) -> Double {
  let n1 = case_vals.length()
  let n2 = control_vals.length()
  if n1 == 0 || n2 == 0 {
    return 0.0
  }
  // Combine and rank
  let combined : Array[(Double, Int)] = [] // (value, group: 0=case, 1=control)
  for v in case_vals {
    combined.push((v, 0))
  }
  for v in control_vals {
    combined.push((v, 1))
  }
  // Sort by value (ascending) using bubble sort
  let n = combined.length()
  let mut i = 0
  while i < n - 1 {
    let mut j = 0
    while j < n - i - 1 {
      if combined[j].0 > combined[j + 1].0 {
        let tmp = combined[j]
        combined[j] = combined[j + 1]
        combined[j + 1] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  // Assign ranks (handle ties with average rank)
  let ranks : Array[Double] = Array::make(n, 0.0)
  let mut idx = 0
  while idx < n {
    let mut j = idx + 1
    while j < n && combined[j].0 == combined[idx].0 {
      j = j + 1
    }
    // Tied from idx to j-1
    let mut avg_rank = 0.0
    for k in idx..<j {
      avg_rank = avg_rank + (k + 1).to_double()
    }
    avg_rank = avg_rank / (j - idx).to_double()
    for k in idx..<j {
      ranks[k] = avg_rank
    }
    idx = j
  }
  // Sum of ranks for case group (group=0)
  let mut w_stat = 0.0
  for k in 0..<n {
    if combined[k].1 == 0 {
      w_stat = w_stat + ranks[k]
    }
  }
  w_stat
}

///|
/// Compute log2 fold change between case and control.
pub fn fish_log2fc_test(
  case_vals : Array[Double],
  control_vals : Array[Double],
) -> Double {
  let case_mean = stat_mean(case_vals)
  let control_mean = stat_mean(control_vals)
  // Add pseudocount of 1 to avoid log(0)
  let case_log = if case_mean > 0.0 { @math.log2(case_mean + 1.0) } else { 0.0 }
  let control_log = if control_mean > 0.0 { @math.log2(control_mean + 1.0) } else { 0.0 }
  case_log - control_log
}

///|
/// Split sample values for a transcript into case and control groups.
fn fish_split_groups(
  row : Array[Double],
  condition : Array[String],
  case_label : String,
) -> (Array[Double], Array[Double]) {
  let case_vals : Array[Double] = []
  let control_vals : Array[Double] = []
  for i in 0..<row.length() {
    if condition[i] == case_label {
      case_vals.push(row[i])
    } else {
      control_vals.push(row[i])
    }
  }
  (case_vals, control_vals)
}

///|
/// Permute condition labels and recompute the Mann-Whitney statistic.
/// Returns the permuted statistic for a single transcript.
fn fish_permuted_statistic(
  row : Array[Double],
  condition : Array[String],
  case_label : String,
  seed : Int,
) -> Double {
  let n = row.length()
  // Count case samples
  let mut n_case = 0
  for c in condition {
    if c == case_label {
      n_case = n_case + 1
    }
  }
  // Shuffle the row values (Fisher-Yates) to simulate permutation
  let shuffled : Array[Double] = []
  for v in row {
    shuffled.push(v)
  }
  let mut rng_state = seed
  let n_minus_1 = n - 1
  for i in 0..<n_minus_1 {
    rng_state = (rng_state * 1103515245 + 12345) & 0x7FFFFFFF
    let r = rng_state % (n - i)
    let swap_idx = i + r
    let tmp = shuffled[i]
    shuffled[i] = shuffled[swap_idx]
    shuffled[swap_idx] = tmp
  }
  // Split shuffled values: first n_case = "case", rest = "control"
  let case_vals : Array[Double] = []
  let control_vals : Array[Double] = []
  for i in 0..<n {
    if i < n_case {
      case_vals.push(shuffled[i])
    } else {
      control_vals.push(shuffled[i])
    }
  }
  fish_mann_whitney_test(case_vals, control_vals)
}

///|
/// Run the Swish non-parametric differential expression analysis.
///
/// Parameters:
/// - counts: transcript × sample count matrix with condition labels
/// - case_label: the condition label for the case group (default "case")
/// - n_permutations: number of permutations for p-value estimation
///
/// Returns differential expression results with permutation p-values
/// and BH-FDR correction.
pub fn fish_swish(
  counts : FishCounts,
  case_label : String,
  n_permutations : Int,
) -> Array[FishResult] {
  let n_transcripts = counts.transcript_names.length()
  let n_samples = counts.sample_names.length()
  if n_transcripts == 0 || n_samples < 2 {
    return []
  }
  let results : Array[FishResult] = []
  let p_values : Array[Double] = Array::make(n_transcripts, 1.0)
  let stats : Array[Double] = Array::make(n_transcripts, 0.0)
  let log2fcs : Array[Double] = Array::make(n_transcripts, 0.0)
  for t in 0..<n_transcripts {
    let row = counts.counts[t]
    let (case_vals, control_vals) = fish_split_groups(
      row,
      counts.condition,
      case_label,
    )
    if case_vals.length() == 0 || control_vals.length() == 0 {
      stats[t] = 0.0
      log2fcs[t] = 0.0
      p_values[t] = 1.0
      continue
    }
    // Observed statistic
    let obs_stat = fish_mann_whitney_test(case_vals, control_vals)
    stats[t] = obs_stat
    log2fcs[t] = fish_log2fc_test(case_vals, control_vals)
    // Expected statistic under null (mean rank sum)
    let n1 = case_vals.length()
    let n2 = control_vals.length()
    let expected = n1.to_double() * (n1 + n2 + 1).to_double() / 2.0
    let obs_dev = fish_abs(obs_stat - expected)
    // Permutation testing
    let mut perm_count = 0
    for p in 0..<n_permutations {
      let seed = 42 + t * 7919 + p * 31
      let perm_stat = fish_permuted_statistic(
        row,
        counts.condition,
        case_label,
        seed,
      )
      let perm_dev = fish_abs(perm_stat - expected)
      if perm_dev >= obs_dev {
        perm_count = perm_count + 1
      }
    }
    p_values[t] = (perm_count.to_double() + 1.0) / (n_permutations.to_double() + 1.0)
  }
  // BH-FDR correction
  let fdrs = stat_bh(p_values)
  // Build results
  for t in 0..<n_transcripts {
    let direction = if fdrs[t] < 0.05 {
      if log2fcs[t] > 0.0 { "up" } else { "down" }
    } else {
      "ns"
    }
    results.push(
      FishResult::new(
        counts.transcript_names[t],
        stats[t],
        log2fcs[t],
        p_values[t],
        fdrs[t],
        direction,
      ),
    )
  }
  results
}

///|
/// Absolute value helper.
fn fish_abs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Generate sample count data with 10 transcripts, 16 samples (8 case + 8 control).
/// Transcripts 0-2 are up-regulated in case, 3-5 down-regulated, 6-9 not significant.
pub fn fish_sample_counts() -> FishCounts {
  let transcript_names = [
    "TX1", "TX2", "TX3", "TX4", "TX5", "TX6", "TX7", "TX8", "TX9", "TX10",
  ]
  let sample_names = [
    "case1", "case2", "case3", "case4", "case5", "case6", "case7", "case8",
    "ctrl1", "ctrl2", "ctrl3", "ctrl4", "ctrl5", "ctrl6", "ctrl7", "ctrl8",
  ]
  let condition = [
    "case", "case", "case", "case", "case", "case", "case", "case",
    "control", "control", "control", "control", "control", "control", "control", "control",
  ]
  // transcript × sample counts (10 × 16)
  let counts : Array[Array[Double]] = []
  // TX1-3: up in case (high in case, low in control)
  for t in 0..<3 {
    let row : Array[Double] = []
    for s in 0..<16 {
      let val = if s < 8 { 50.0 + t.to_double() * 10.0 + s.to_double() * 2.0 } else { 10.0 + (s - 8).to_double() * 1.0 }
      row.push(val)
    }
    counts.push(row)
  }
  // TX4-6: down in case (low in case, high in control)
  for t in 0..<3 {
    let row : Array[Double] = []
    for s in 0..<16 {
      let val = if s < 8 { 10.0 + s.to_double() * 1.0 } else { 50.0 + t.to_double() * 10.0 + (s - 8).to_double() * 2.0 }
      row.push(val)
    }
    counts.push(row)
  }
  // TX7-10: not significant (constant values per transcript; all samples identical
  // so case/control distributions are identical and W == expected, p == 1)
  for t in 0..<4 {
    let row : Array[Double] = []
    let const_val = 25.0 + t.to_double() * 5.0
    for s in 0..<16 {
      ignore(s)
      row.push(const_val)
    }
    counts.push(row)
  }
  FishCounts::new(transcript_names, sample_names, condition, counts)
}

// ============================================================================
// Utilities
// ============================================================================

///|
/// Filter results by FDR threshold.
pub fn fish_significant(
  results : Array[FishResult],
  fdr_threshold : Double,
) -> Array[FishResult] {
  let sig : Array[FishResult] = []
  for r in results {
    if r.fdr < fdr_threshold {
      sig.push(r)
    }
  }
  sig
}

///|
/// Format a FishResult as a readable string.
pub fn FishResult::to_string(self : FishResult) -> String {
  self.transcript + " log2FC=" + self.log2_fold_change.to_string() + " p=" + self.p_value.to_string() + " fdr=" + self.fdr.to_string() + " dir=" + self.direction
}

///|
/// Format all results as a table.
pub fn fish_results_to_string(results : Array[FishResult]) -> String {
  let mut s = "Swish Results (" + results.length().to_string() + " transcripts)\n"
  s = s + "transcript\tstatistic\tlog2FC\tp_value\tfdr\tdirection\n"
  for r in results {
    s = s + r.transcript + "\t" + r.statistic.to_string() + "\t" + r.log2_fold_change.to_string() + "\t" + r.p_value.to_string() + "\t" + r.fdr.to_string() + "\t" + r.direction + "\n"
  }
  s
}
