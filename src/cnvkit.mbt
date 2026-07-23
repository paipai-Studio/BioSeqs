///|
/// CNVkit - Copy number variation detection and CBS segmentation.
/// Inspired by Bioconductor CNVkit.
/// Provides copy number data structures, CBS binary segmentation algorithm,
/// log2 ratio smoothing, copy number state calling, breakpoint detection,
/// and CNV summary statistics.

///|
pub struct CNVProbe {
  probe_id : String
  chromosome : String
  start : Int
  end : Int
  log2_ratio : Double
  weight : Double
}

///|
pub struct CNVSegment {
  chromosome : String
  start : Int
  end : Int
  probes : Int
  mean_log2 : Double
  state : String
}

///|
pub struct CNVDataset {
  probes : Array[CNVProbe]
  segments : Array[CNVSegment]
  chromosomes : Array[String]
}

///|
pub struct CBSResult {
  segments : Array[CNVSegment]
  breakpoints : Array[Int]
  n_segments : Int
}

///|
pub struct CopyNumberCall {
  segment : CNVSegment
  copy_number : Double
  state : String
}

///|
pub struct CNVSummary {
  total_probes : Int
  total_segments : Int
  avg_segment_size : Double
  deletion_count : Int
  amplification_count : Int
  neutral_count : Int
}

///|
pub fn CNVProbe::new(
  probe_id : String,
  chromosome : String,
  start : Int,
  end : Int,
  log2_ratio : Double,
  weight : Double
) -> CNVProbe {
  CNVProbe::{
    probe_id,
    chromosome,
    start,
    end,
    log2_ratio,
    weight
  }
}

///|
pub fn CNVSegment::new(
  chromosome : String,
  start : Int,
  end : Int,
  probes : Int,
  mean_log2 : Double,
  state : String
) -> CNVSegment {
  CNVSegment::{
    chromosome,
    start,
    end,
    probes,
    mean_log2,
    state
  }
}

///|
pub fn CNVDataset::new() -> CNVDataset {
  CNVDataset::{
    probes: Array::new(),
    segments: Array::new(),
    chromosomes: Array::new()
  }
}

///|
pub fn CNVDataset::add_probe(self : CNVDataset, probe : CNVProbe) -> CNVDataset {
  let new_probes = self.probes
  new_probes.push(probe)
  
  let new_chroms = self.chromosomes
  if !cnvkit_contains_string(new_chroms, probe.chromosome) {
    new_chroms.push(probe.chromosome)
  }
  
  CNVDataset::{
    probes: new_probes,
    segments: self.segments,
    chromosomes: new_chroms
  }
}

///|
pub fn CNVDataset::count_probes(self : CNVDataset) -> Int {
  self.probes.length()
}

///|
pub fn CNVDataset::count_segments(self : CNVDataset) -> Int {
  self.segments.length()
}

///|
pub fn CNVDataset::filter_chromosome(self : CNVDataset, chrom : String) -> CNVDataset {
  let filtered : Array[CNVProbe] = Array::new()
  let mut i = 0
  while i < self.probes.length() {
    if self.probes[i].chromosome == chrom {
      filtered.push(self.probes[i])
    }
    i = i + 1
  }
  
  let segments : Array[CNVSegment] = Array::new()
  let mut j = 0
  while j < self.segments.length() {
    if self.segments[j].chromosome == chrom {
      segments.push(self.segments[j])
    }
    j = j + 1
  }
  
  CNVDataset::{
    probes: filtered,
    segments,
    chromosomes: [chrom]
  }
}

///|
pub fn cbs_segment(probes : Array[CNVProbe], alpha : Double) -> CBSResult {
  let n = probes.length()
  if n == 0 {
    return CBSResult::{
      segments: Array::new(),
      breakpoints: Array::new(),
      n_segments: 0
    }
  }
  
  // Extract log2 ratios
  let log2ratios : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    log2ratios.push(probes[i].log2_ratio)
    i = i + 1
  }
  
  // Perform CBS segmentation
  let breakpoints = cbs_find_breakpoints(log2ratios, alpha)
  
  // Create segments from breakpoints
  let segments = cbs_create_segments(probes, breakpoints)
  
  CBSResult::{
    segments,
    breakpoints,
    n_segments: segments.length()
  }
}

///|
fn cbs_find_breakpoints(data : Array[Double], alpha : Double) -> Array[Int] {
  let n = data.length()
  let breakpoints : Array[Int] = Array::new()
  
  if n == 0 {
    return breakpoints
  }
  
  if n < 3 {
    breakpoints.push(0)
    breakpoints.push(n)
    return breakpoints
  }
  
  // Find change points using t-statistic
  let mut i = 1
  while i < n - 1 {
    let t_stat = cbs_calculate_t_statistic(data, i)
    let p_value = cbs_approximate_p_value(t_stat, n)
    
    if p_value < alpha {
      breakpoints.push(i)
    }
    i = i + 1
  }
  
  // Ensure breakpoints include start and end
  if breakpoints.length() > 0 {
    if breakpoints[0] != 0 {
      breakpoints.insert(0, 0)
    }
    if breakpoints[breakpoints.length() - 1] != n {
      breakpoints.push(n)
    }
  } else {
    breakpoints.push(0)
    breakpoints.push(n)
  }
  
  breakpoints
}

///|
fn cbs_calculate_t_statistic(data : Array[Double], breakpoint : Int) -> Double {
  let n = data.length()
  let n_left = breakpoint
  let n_right = n - breakpoint
  
  if n_left < 1 || n_right < 1 {
    return 0.0
  }
  
  // Calculate means
  let mut sum_left = 0.0
  let mut i = 0
  while i < n_left {
    sum_left = sum_left + data[i]
    i = i + 1
  }
  let mean_left = sum_left / n_left.to_double()
  
  let mut sum_right = 0.0
  let mut j = n_left
  while j < n {
    sum_right = sum_right + data[j]
    j = j + 1
  }
  let mean_right = sum_right / n_right.to_double()
  
  // Calculate variances
  let mut var_left = 0.0
  i = 0
  while i < n_left {
    let diff = data[i] - mean_left
    var_left = var_left + diff * diff
    i = i + 1
  }
  var_left = var_left / n_left.to_double()
  
  let mut var_right = 0.0
  j = n_left
  while j < n {
    let diff = data[j] - mean_right
    var_right = var_right + diff * diff
    j = j + 1
  }
  var_right = var_right / n_right.to_double()
  
  let pooled_var = (n_left.to_double() * var_left + n_right.to_double() * var_right) / 
    (n_left.to_double() + n_right.to_double())
  let se = (pooled_var * (1.0 / n_left.to_double() + 1.0 / n_right.to_double())).sqrt()
  
  if se == 0.0 {
    0.0
  } else {
    (mean_left - mean_right) / se
  }
}

///|
fn cbs_approximate_p_value(t : Double, n : Int) -> Double {
  // Approximate p-value using normal distribution
  let abs_t = if t >= 0.0 { t } else { -t }
  let df = n.to_double() - 2.0
  
  // Simple approximation
  let x = abs_t / (1.0 + abs_t / (2.0 * df))
  2.0 * (1.0 - cbs_normal_cdf(x))
}

///|
fn cbs_normal_cdf(x : Double) -> Double {
  // Standard normal CDF approximation
  let abs_x = if x >= 0.0 { x } else { -x }
  let t = 1.0 / (1.0 + 0.2316419 * abs_x)
  let d = 0.3989422804014327 * @math.exp(-abs_x * abs_x / 2.0)
  let p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))))
  if x >= 0.0 {
    1.0 - p
  } else {
    p
  }
}

///|
fn cbs_create_segments(probes : Array[CNVProbe], breakpoints : Array[Int]) -> Array[CNVSegment] {
  let segments : Array[CNVSegment] = Array::new()
  let n = breakpoints.length()
  
  if n < 2 {
    return segments
  }
  
  let mut i = 0
  while i < n - 1 {
    let start_idx = breakpoints[i]
    let end_idx = breakpoints[i + 1]
    
    if start_idx < end_idx && start_idx < probes.length() {
      let chrom = probes[start_idx].chromosome
      let start_pos = probes[start_idx].start
      let end_pos = if end_idx > 0 && end_idx <= probes.length() {
        probes[end_idx - 1].end
      } else {
        probes[probes.length() - 1].end
      }
      let count = end_idx - start_idx
      
      // Calculate mean log2 ratio
      let mut sum = 0.0
      let mut j = start_idx
      while j < end_idx && j < probes.length() {
        sum = sum + probes[j].log2_ratio
        j = j + 1
      }
      let mean_log2 = if count > 0 { sum / count.to_double() } else { 0.0 }
      
      // Determine state
      let state = cbs_call_state(mean_log2)
      
      segments.push(CNVSegment::new(
        chrom,
        start_pos,
        end_pos,
        count,
        mean_log2,
        state
      ))
    }
    i = i + 1
  }
  
  segments
}

///|
pub fn cbs_call_state(mean_log2 : Double) -> String {
  if mean_log2 < -0.3 {
    "deletion"
  } else if mean_log2 > 0.3 {
    "amplification"
  } else {
    "neutral"
  }
}

///|
pub fn smooth_log2_ratios(probes : Array[CNVProbe], window_size : Int) -> Array[CNVProbe] {
  let n = probes.length()
  let smoothed : Array[CNVProbe] = Array::new()
  
  let mut i = 0
  while i < n {
    let start = if i - window_size / 2 >= 0 { i - window_size / 2 } else { 0 }
    let end = if i + window_size / 2 < n { i + window_size / 2 } else { n - 1 }
    
    let mut sum = 0.0
    let mut count = 0
    let mut j = start
    while j <= end {
      sum = sum + probes[j].log2_ratio
      count = count + 1
      j = j + 1
    }
    
    let smoothed_ratio = if count > 0 { sum / count.to_double() } else { probes[i].log2_ratio }
    
    smoothed.push(CNVProbe::{
      ..probes[i],
      log2_ratio: smoothed_ratio
    })
    i = i + 1
  }
  
  smoothed
}

///|
pub fn detect_breakpoints(segments : Array[CNVSegment], threshold : Double) -> Array[(Int, String)] {
  let breakpoints : Array[(Int, String)] = Array::new()
  let n = segments.length()
  
  if n < 2 {
    return breakpoints
  }
  
  let mut i = 1
  while i < n {
    let prev = segments[i - 1]
    let curr = segments[i]
    
    // Check for state change
    if prev.state != curr.state {
      breakpoints.push((curr.start, curr.state))
    }
    
    // Check for significant mean change
    let mean_diff = (prev.mean_log2 - curr.mean_log2).abs()
    if mean_diff > threshold {
      breakpoints.push((curr.start, "change"))
    }
    
    i = i + 1
  }
  
  breakpoints
}

///|
pub fn call_copy_numbers(segments : Array[CNVSegment], ploidy : Int) -> Array[CopyNumberCall] {
  let calls : Array[CopyNumberCall] = Array::new()
  
  let mut i = 0
  while i < segments.length() {
    let segment = segments[i]
    let copy_number = cbs_log2_to_copynumber(segment.mean_log2, ploidy)
    let state = cbs_call_state(segment.mean_log2)
    
    calls.push(CopyNumberCall::{
      segment,
      copy_number,
      state
    })
    i = i + 1
  }
  
  calls
}

///|
fn cbs_log2_to_copynumber(log2_ratio : Double, ploidy : Int) -> Double {
  let ploidy_double = ploidy.to_double()
  let ratio = 2.0.pow(log2_ratio)
  (ploidy_double * ratio).round()
}

///|
pub fn summarize_cnv(dataset : CNVDataset) -> CNVSummary {
  let total_probes = dataset.count_probes()
  let total_segments = dataset.count_segments()
  
  let avg_segment_size = if total_segments > 0 {
    total_probes.to_double() / total_segments.to_double()
  } else {
    0.0
  }
  
  let mut deletion_count = 0
  let mut amplification_count = 0
  let mut neutral_count = 0
  
  let segments = dataset.segments
  let mut i = 0
  while i < segments.length() {
    match segments[i].state {
      "deletion" => deletion_count = deletion_count + 1
      "amplification" => amplification_count = amplification_count + 1
      _ => neutral_count = neutral_count + 1
    }
    i = i + 1
  }
  
  CNVSummary::{
    total_probes,
    total_segments,
    avg_segment_size,
    deletion_count,
    amplification_count,
    neutral_count
  }
}

///|
pub fn create_example_cnv_dataset() -> CNVDataset {
  let mut dataset = CNVDataset::new()
  
  // Create a dataset with 100 probes showing a copy number change
  let mut i = 0
  while i < 100 {
    let chrom = if i < 50 { "chr1" } else { "chr2" }
    let start = if chrom == "chr1" { 1000000 + i * 10000 } else { 2000000 + (i - 50) * 10000 }
    
    // Create a deletion in the middle of chr1
    let log2_ratio = if chrom == "chr1" && i >= 20 && i < 35 {
      -0.8
    } else if chrom == "chr1" && i >= 50 && i < 65 {
      0.8
    } else {
      if i % 3 == 0 { 0.05 } else if i % 3 == 1 { -0.03 } else { 0.02 }
    }
    
    let probe = CNVProbe::new(
      "probe_" + i.to_string(),
      chrom,
      start,
      start + 1000,
      log2_ratio,
      1.0
    )
    
    dataset = dataset.add_probe(probe)
    i = i + 1
  }
  
  dataset
}

///|
pub fn cbs_summarize_dataset(dataset : CNVDataset) -> String {
  let summary = summarize_cnv(dataset)
  let mut result = "\nCNV Dataset Summary\n"
  result = result + "==============================\n"
  result = result + "Total probes: " + summary.total_probes.to_string() + "\n"
  result = result + "Total segments: " + summary.total_segments.to_string() + "\n"
  result = result + "Avg segment size: " + summary.avg_segment_size.to_string() + "\n"
  result = result + "Deletions: " + summary.deletion_count.to_string() + "\n"
  result = result + "Amplifications: " + summary.amplification_count.to_string() + "\n"
  result = result + "Neutral: " + summary.neutral_count.to_string() + "\n"
  result
}

// Helper functions

fn cnvkit_contains_string(arr : Array[String], value : String) -> Bool {
  let mut i = 0
  while i < arr.length() {
    if arr[i] == value {
      return true
    }
    i = i + 1
  }
  false
}
