// HMMcopy - Copy Number Variation detection using Hidden Markov Model.
// Inspired by Bioconductor HMMcopy package.
//
// Calls copy number states from read coverage in genomic bins using a
// Hidden Markov Model with Viterbi decoding. Includes GC bias correction
// and segmentation of consecutive bins with the same state.
//
// References:
//   - Lai D et al. (2018) HMMcopy: Copy number prediction with read
//     depth information. Bioconductor package.
//   - Ha G et al. (2012) Integrative analysis of genome-wide loss of
//     heterozygosity and mono-allelic expression at nucleotide
//     resolution. Bioinformatics 28:2796-2802.

// ===========================================================================
// Data structures
// ===========================================================================

///|
/// A genomic bin with read coverage.
pub struct HMMcopyBin {
  pub chr : String
  pub start : Int
  pub end : Int
  pub reads : Int
  pub gc : Double // GC content (0..1)
  pub mappability : Double // mappability score (0..1)
  pub mut corrected_reads : Double // GC-corrected coverage
  pub mut copy : Int // inferred copy number (after Viterbi)
  pub mut state : Int // HMM state index (0=CN0,1=CN1,2=CN2,...)
}

///|
/// Construct a new HMMcopyBin.
pub fn HMMcopyBin::new(
  chr~ : String,
  start~ : Int,
  end~ : Int,
  reads~ : Int,
  gc~ : Double = 0.5,
  mappability~ : Double = 1.0,
) -> HMMcopyBin {
  HMMcopyBin::{
    chr,
    start,
    end,
    reads,
    gc,
    mappability,
    corrected_reads: reads.to_double(),
    copy: 2,
    state: 2,
  }
}

///|
/// Get bin width.
pub fn HMMcopyBin::width(self : HMMcopyBin) -> Int {
  self.end - self.start
}

///|
/// Set corrected reads.
pub fn HMMcopyBin::set_corrected_reads(self : HMMcopyBin, val : Double) -> Unit {
  self.corrected_reads = val
}

///|
/// Set copy number.
pub fn HMMcopyBin::set_copy(self : HMMcopyBin, val : Int) -> Unit {
  self.copy = val
}

///|
/// Set HMM state index.
pub fn HMMcopyBin::set_state(self : HMMcopyBin, val : Int) -> Unit {
  self.state = val
}

///|
/// Get corrected coverage per bp.
pub fn HMMcopyBin::normalized_coverage(self : HMMcopyBin) -> Double {
  let w = self.width()
  if w == 0 {
    self.corrected_reads
  } else {
    self.corrected_reads / w.to_double()
  }
}

///|
/// HMM model parameters.
pub struct HMMcopyModel {
  pub num_states : Int // e.g. 7: CN0..CN6
  pub transition : Array[Array[Double]] // transition[s_from][s_to]
  pub emissions : Array[Double] // mean coverage per state (in CN2 units)
  pub variances : Array[Double] // variance per state
  pub initial : Array[Double] // initial state distribution
  pub ploidy : Int // expected copy number (default 2)
}

///|
/// Construct a new HMMcopyModel with default parameters.
/// States: CN0, CN1, CN2, CN3, CN4, CN5, CN6 (7 states).
pub fn HMMcopyModel::new(
  num_states~ : Int = 7,
  ploidy~ : Int = 2,
) -> HMMcopyModel {
  // Default emission means: 0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 (in CN2 units)
  let emissions = Array::make(num_states, 0.0)
  let variances = Array::make(num_states, 0.1)
  let mut i = 0
  while i < num_states {
    emissions[i] = i.to_double() / ploidy.to_double()
    variances[i] = 0.1 + i.to_double() * 0.02
    i = i + 1
  }
  // Default transition: high self-transition, low cross-transition
  let transition = Array::make(num_states, Array::make(num_states, 0.0))
  i = 0
  while i < num_states {
    let mut j = 0
    while j < num_states {
      if i == j {
        transition[i][j] = 0.999
      } else {
        transition[i][j] = 0.001 / (num_states - 1).to_double()
      }
      j = j + 1
    }
    i = i + 1
  }
  // Default initial: uniform
  let initial = Array::make(num_states, 1.0 / num_states.to_double())
  HMMcopyModel::{
    num_states,
    transition,
    emissions,
    variances,
    initial,
    ploidy,
  }
}

///|
/// Set a transition probability.
pub fn HMMcopyModel::set_transition(
  self : HMMcopyModel,
  from : Int,
  to : Int,
  prob : Double,
) -> Unit {
  self.transition[from][to] = prob
}

///|
/// Set emission mean for a state (in CN2 units).
pub fn HMMcopyModel::set_emission(
  self : HMMcopyModel,
  state : Int,
  mean : Double,
) -> Unit {
  self.emissions[state] = mean
}

///|
/// Set emission variance for a state.
pub fn HMMcopyModel::set_variance(
  self : HMMcopyModel,
  state : Int,
  var : Double,
) -> Unit {
  self.variances[state] = var
}

///|
/// Set initial probability for a state.
pub fn HMMcopyModel::set_initial(
  self : HMMcopyModel,
  state : Int,
  prob : Double,
) -> Unit {
  self.initial[state] = prob
}

///|
/// Get the copy number for a state index.
pub fn HMMcopyModel::state_to_copy(self : HMMcopyModel, state : Int) -> Int {
  state
}

///|
/// Get the state index for a copy number.
pub fn HMMcopyModel::copy_to_state(self : HMMcopyModel, copy : Int) -> Int {
  if copy < 0 {
    0
  } else if copy >= self.num_states {
    self.num_states - 1
  } else {
    copy
  }
}

///|
/// A segment of consecutive bins with the same copy number.
pub struct HMMcopySegment {
  pub chr : String
  pub start : Int
  pub end : Int
  pub copy : Int
  pub state : Int
  pub num_bins : Int
  pub mean_reads : Double
}

///|
/// Construct a HMMcopySegment.
pub fn HMMcopySegment::new(
  chr~ : String,
  start~ : Int,
  end~ : Int,
  copy~ : Int,
  state~ : Int,
  num_bins~ : Int,
  mean_reads~ : Double,
) -> HMMcopySegment {
  HMMcopySegment::{
    chr,
    start,
    end,
    copy,
    state,
    num_bins,
    mean_reads,
  }
}

// ===========================================================================
// GC bias correction
// ===========================================================================

///|
/// Correct GC bias in read counts using LOESS-like binning.
/// Bins are grouped by GC content and scaled to the global mean.
pub fn hmmcopy_correct_gc(
  bins : Array[HMMcopyBin],
  num_gc_bins~ : Int = 20,
) -> Array[HMMcopyBin] {
  let n = bins.length()
  if n == 0 {
    return []
  }
  // Compute global mean of reads
  let mut total_reads = 0
  for b in bins {
    total_reads = total_reads + b.reads
  }
  let global_mean = total_reads.to_double() / n.to_double()
  // Group bins by GC content
  let gc_bins = Array::make(num_gc_bins, Array::new())
  for b in bins {
    let gc_idx = (b.gc * num_gc_bins.to_double()).to_int()
    let clamped = if gc_idx < 0 { 0 } else if gc_idx >= num_gc_bins { num_gc_bins - 1 } else { gc_idx }
    gc_bins[clamped].push(b)
  }
  // Compute mean reads per GC bin
  let gc_means = Array::make(num_gc_bins, 0.0)
  let mut i = 0
  while i < num_gc_bins {
    let group = gc_bins[i]
    if group.length() > 0 {
      let mut sum = 0.0
      for b in group {
        sum = sum + b.reads.to_double()
      }
      gc_means[i] = sum / group.length().to_double()
    } else {
      gc_means[i] = global_mean
    }
    i = i + 1
  }
  // Correct each bin: scale reads so that GC bin mean equals global mean
  let result = Array::new()
  for b in bins {
    let gc_idx = (b.gc * num_gc_bins.to_double()).to_int()
    let clamped = if gc_idx < 0 {
      0
    } else if gc_idx >= num_gc_bins {
      num_gc_bins - 1
    } else {
      gc_idx
    }
    let gc_mean = gc_means[clamped]
    let corrected = if gc_mean > 0.0 {
      b.reads.to_double() * global_mean / gc_mean
    } else {
      b.reads.to_double()
    }
    let new_bin = HMMcopyBin::new(
      chr=b.chr,
      start=b.start,
      end=b.end,
      reads=b.reads,
      gc=b.gc,
      mappability=b.mappability,
    )
    new_bin.corrected_reads = corrected
    result.push(new_bin)
  }
  result
}

// ===========================================================================
// Viterbi algorithm
// ===========================================================================

///|
/// Compute Gaussian log emission probability.
fn hmmcopy_log_gaussian(x : Double, mean : Double, variance : Double) -> Double {
  if variance <= 0.0 {
    return -1.0e30
  }
  let diff = x - mean
  let log_2pi = 1.8378770664093453
  let result = -0.5 * (log_2pi + @math.ln(variance)) -
    diff * diff /
    (2.0 * variance)
  result
}

///|
/// Run Viterbi algorithm to find most likely state sequence.
/// Returns array of state indices (one per bin).
pub fn hmmcopy_viterbi(
  bins : Array[HMMcopyBin],
  model : HMMcopyModel,
) -> Array[Int] {
  let n = bins.length()
  let k = model.num_states
  if n == 0 || k == 0 {
    return []
  }
  // Compute mean coverage of CN2 state (used to normalize)
  let mut mean_cov = 0.0
  for b in bins {
    mean_cov = mean_cov + b.corrected_reads
  }
  mean_cov = mean_cov / n.to_double()
  if mean_cov <= 0.0 {
    mean_cov = 1.0
  }
  // Viterbi matrix: viterbi[t][s] = best log-prob of being in state s at time t
  let viterbi : Array[Array[Double]] = Array::new()
  let backptr : Array[Array[Int]] = Array::new()
  let mut vi = 0
  while vi < n {
    viterbi.push(Array::make(k, -1.0e30))
    backptr.push(Array::make(k, 0))
    vi = vi + 1
  }
  // Initialize t=0
  let first_cov = bins[0].corrected_reads / mean_cov
  let mut s = 0
  while s < k {
    let expected = model.emissions[s]
    let log_emit = hmmcopy_log_gaussian(first_cov, expected, model.variances[s])
    viterbi[0][s] = @math.ln(model.initial[s]) + log_emit
    backptr[0][s] = 0
    s = s + 1
  }
  // Recursion
  let mut t = 1
  while t < n {
    let cov = bins[t].corrected_reads / mean_cov
    let mut s_to = 0
    while s_to < k {
      let log_emit = hmmcopy_log_gaussian(
        cov,
        model.emissions[s_to],
        model.variances[s_to],
      )
      let mut best_log = -1.0e30
      let mut best_from = 0
      let mut s_from = 0
      while s_from < k {
        let log_prob = viterbi[t - 1][s_from] +
          @math.ln(model.transition[s_from][s_to])
        if log_prob > best_log {
          best_log = log_prob
          best_from = s_from
        }
        s_from = s_from + 1
      }
      viterbi[t][s_to] = best_log + log_emit
      backptr[t][s_to] = best_from
      s_to = s_to + 1
    }
    t = t + 1
  }
  // Termination: find best final state
  let mut best_final = 0
  let mut best_log = -1.0e30
  s = 0
  while s < k {
    if viterbi[n - 1][s] > best_log {
      best_log = viterbi[n - 1][s]
      best_final = s
    }
    s = s + 1
  }
  // Backtrace
  let states = Array::make(n, 0)
  states[n - 1] = best_final
  t = n - 1
  while t > 0 {
    states[t - 1] = backptr[t][states[t]]
    t = t - 1
  }
  states
}

// ===========================================================================
// CNV calling and segmentation
// ===========================================================================

///|
/// Call copy number for each bin using Viterbi.
/// Modifies bins in-place: sets copy and state fields.
pub fn hmmcopy_call_cnv(
  bins : Array[HMMcopyBin],
  model : HMMcopyModel,
) -> Array[HMMcopyBin] {
  let states = hmmcopy_viterbi(bins, model)
  let result = Array::new()
  let mut i = 0
  while i < bins.length() {
    let b = HMMcopyBin::new(
      chr=bins[i].chr,
      start=bins[i].start,
      end=bins[i].end,
      reads=bins[i].reads,
      gc=bins[i].gc,
      mappability=bins[i].mappability,
    )
    b.corrected_reads = bins[i].corrected_reads
    b.state = states[i]
    b.copy = states[i]
    result.push(b)
    i = i + 1
  }
  result
}

///|
/// Segment bins by merging consecutive bins with the same copy number.
pub fn hmmcopy_segment(
  bins : Array[HMMcopyBin],
) -> Array[HMMcopySegment] {
  let segments = Array::new()
  let n = bins.length()
  if n == 0 {
    return segments
  }
  let mut i = 0
  while i < n {
    let cur_copy = bins[i].copy
    let cur_chr = bins[i].chr
    let start = bins[i].start
    let mut end = bins[i].end
    let mut sum_reads = bins[i].corrected_reads
    let mut count = 1
    let mut j = i + 1
    while j < n && bins[j].copy == cur_copy && bins[j].chr == cur_chr {
      end = bins[j].end
      sum_reads = sum_reads + bins[j].corrected_reads
      count = count + 1
      j = j + 1
    }
    segments.push(HMMcopySegment::new(
      chr=cur_chr,
      start=start,
      end=end,
      copy=cur_copy,
      state=bins[i].state,
      num_bins=count,
      mean_reads=sum_reads / count.to_double(),
    ))
    i = j
  }
  segments
}

///|
/// Compute log-likelihood of the observed data given the model.
pub fn hmmcopy_log_likelihood(
  bins : Array[HMMcopyBin],
  model : HMMcopyModel,
) -> Double {
  let n = bins.length()
  if n == 0 {
    return 0.0
  }
  let mut mean_cov = 0.0
  for b in bins {
    mean_cov = mean_cov + b.corrected_reads
  }
  mean_cov = mean_cov / n.to_double()
  if mean_cov <= 0.0 {
    mean_cov = 1.0
  }
  // Forward algorithm
  let k = model.num_states
  let alpha : Array[Array[Double]] = Array::new()
  let mut ai = 0
  while ai < n {
    alpha.push(Array::make(k, 0.0))
    ai = ai + 1
  }
  // Initialize
  let first_cov = bins[0].corrected_reads / mean_cov
  let mut s = 0
  while s < k {
    let log_emit = hmmcopy_log_gaussian(
      first_cov,
      model.emissions[s],
      model.variances[s],
    )
    alpha[0][s] = model.initial[s] * @math.exp(log_emit)
    s = s + 1
  }
  // Recursion
  let mut t = 1
  while t < n {
    let cov = bins[t].corrected_reads / mean_cov
    let mut s_to = 0
    while s_to < k {
      let log_emit = hmmcopy_log_gaussian(
        cov,
        model.emissions[s_to],
        model.variances[s_to],
      )
      let emit = @math.exp(log_emit)
      let mut sum = 0.0
      let mut s_from = 0
      while s_from < k {
        sum = sum + alpha[t - 1][s_from] * model.transition[s_from][s_to]
        s_from = s_from + 1
      }
      alpha[t][s_to] = sum * emit
      s_to = s_to + 1
    }
    // Scale to prevent underflow
    let mut scale_sum = 0.0
    let mut si = 0
    while si < k {
      scale_sum = scale_sum + alpha[t][si]
      si = si + 1
    }
    if scale_sum > 0.0 {
      si = 0
      while si < k {
        alpha[t][si] = alpha[t][si] / scale_sum
        si = si + 1
      }
    }
    t = t + 1
  }
  // Sum final
  let mut total = 0.0
  s = 0
  while s < k {
    total = total + alpha[n - 1][s]
    s = s + 1
  }
  if total > 0.0 {
    @math.ln(total)
  } else {
    -1.0e30
  }
}

///|
/// Generate sample bins for testing (simulated CNV data).
/// 100 bins on chr1, with a CN0 region (bins 10-19) and CN3 region (bins 60-69).
pub fn hmmcopy_sample_bins() -> Array[HMMcopyBin] {
  let bins = Array::new()
  let mut i = 0
  while i < 100 {
    let reads = if i >= 10 && i < 20 {
      // CN0 region: very low reads
      2 + (i % 3)
    } else if i >= 60 && i < 70 {
      // CN3 region: elevated reads
      90 + (i % 10)
    } else {
      // CN2 region: normal
      60 + (i % 10)
    }
    let gc = 0.4 + (i % 10).to_double() * 0.02
    let bin = HMMcopyBin::new(
      chr="chr1",
      start=i * 100000,
      end=(i + 1) * 100000,
      reads=reads,
      gc=gc,
      mappability=0.95,
    )
    bins.push(bin)
    i = i + 1
  }
  bins
}

///|
/// Generate a sample model for testing.
pub fn hmmcopy_sample_model() -> HMMcopyModel {
  let model = HMMcopyModel::new(num_states=7, ploidy=2)
  // Set emissions: CN0=0, CN1=0.5, CN2=1.0, CN3=1.5, CN4=2.0, CN5=2.5, CN6=3.0
  let means = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
  let mut i = 0
  while i < 7 {
    model.set_emission(i, means[i])
    model.set_variance(i, 0.05)
    i = i + 1
  }
  // Set initial: strongly prefer CN2
  i = 0
  while i < 7 {
    if i == 2 {
      model.set_initial(i, 0.9)
    } else {
      model.set_initial(i, 0.1 / 6.0)
    }
    i = i + 1
  }
  model
}
