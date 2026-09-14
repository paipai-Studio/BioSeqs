///|
/// dada2 - Divisive Amplicon Denoising Algorithm 2.
/// Infers exact amplicon sequence variants (ASVs) from high-throughput sequencing data.
/// Inspired by Bioconductor dada2 package (Callahan et al. 2016).
/// Core pipeline: dereplicate → learnErrors → dada → mergePairs → removeChimeras → makeTable → assignTaxonomy.

///|
/// Nucleotide encoding: A=0, C=1, G=2, T=3
pub fn nt_to_int(c : Char) -> Int {
  match c.to_int() {
    65 | 97 => 0 // A or a
    67 | 99 => 1 // C or c
    71 | 103 => 2 // G or g
    84 | 116 => 3 // T or t
    _ => -1
  }
}

///|
/// Int to nucleotide character.
pub fn int_to_nt(i : Int) -> String {
  match i {
    0 => "A"
    1 => "C"
    2 => "G"
    3 => "T"
    _ => "N"
  }
}

///|
/// Compute Hamming distance between two equal-length sequences.
pub fn dada_hamming_distance(a : String, b : String) -> Int {
  if a.length() != b.length() {
    return -1
  }
  let a_bytes = a.to_array()
  let b_bytes = b.to_array()
  let mut dist = 0
  let mut i = 0
  while i < a_bytes.length() {
    if a_bytes[i] != b_bytes[i] {
      dist = dist + 1
    }
    i = i + 1
  }
  dist
}

///|
/// Dereplicated sequence entry.
pub struct DerepEntry {
  sequence : String
  abundance : Int
  qual_scores : Array[Double] // average quality per position
}

///|
/// Dereplication result.
pub struct Dereplication {
  uniques : Array[DerepEntry]
  total_reads : Int
}

///|
/// Dereplicate sequences: collapse identical sequences, sum abundances, average qualities.
/// sequences: array of DNA strings
/// quals: array of quality arrays (Phred scores), can be empty
/// Returns sorted by abundance descending.
pub fn dereplicate(
  sequences : Array[String],
  quals : Array[Array[Double]],
) -> Dereplication {
  let seq_map : Map[String, (Int, Array[Double])] = Map([], capacity=sequences.length())
  let has_quals = quals.length() == sequences.length() && quals.length() > 0
  let mut i = 0
  while i < sequences.length() {
    let seq = sequences[i]
    if seq_map.contains(seq) {
      let v = seq_map[seq]
      let new_abund = v.0 + 1
      let new_qual : Array[Double] = []
      if has_quals {
        let mut j = 0
        let qlen = v.1.length()
        while j < qlen {
          new_qual.push((v.1[j] * v.0.to_double() + quals[i][j]) / new_abund.to_double())
          j = j + 1
        }
      }
      seq_map[seq] = (new_abund, new_qual)
    } else {
      let q : Array[Double] = []
      if has_quals {
        let mut j = 0
        while j < quals[i].length() {
          q.push(quals[i][j])
          j = j + 1
        }
      }
      seq_map[seq] = (1, q)
    }
    i = i + 1
  }
  let entries : Array[DerepEntry] = []
  let keys = seq_map.keys().collect()
  let mut ki = 0
  while ki < keys.length() {
    let k = keys[ki]
    let v = seq_map[k]
    entries.push(DerepEntry::{ sequence: k, abundance: v.0, qual_scores: v.1 })
    ki = ki + 1
  }
  // Sort by abundance descending
  entries.sort_by(fn(a, b) {
    if b.abundance > a.abundance {
      1
    } else if b.abundance < a.abundance {
      -1
    } else {
      0
    }
  })
  let mut total = 0
  let mut ei = 0
  while ei < entries.length() {
    total = total + entries[ei].abundance
    ei = ei + 1
  }
  Dereplication::{ uniques: entries, total_reads: total }
}

///|
/// Quality profile: per-position mean, min, max, and quartiles.
pub struct QualityProfile {
  positions : Int
  mean : Array[Double]
  min : Array[Double]
  max : Array[Double]
  q25 : Array[Double]
  q75 : Array[Double]
}

///|
/// Compute quality profile from a set of quality arrays.
pub fn quality_profile(quals : Array[Array[Double]]) -> QualityProfile {
  if quals.length() == 0 {
    return QualityProfile::{ positions: 0, mean: [], min: [], max: [], q25: [], q75: [] }
  }
  let n_pos = quals[0].length()
  let mean_arr : Array[Double] = Array::make(n_pos, 0.0)
  let min_arr : Array[Double] = Array::make(n_pos, 1.0e18)
  let max_arr : Array[Double] = Array::make(n_pos, -1.0e18)
  let q25_arr : Array[Double] = Array::make(n_pos, 0.0)
  let q75_arr : Array[Double] = Array::make(n_pos, 0.0)
  let mut i = 0
  while i < quals.length() {
    let mut j = 0
    while j < n_pos {
      let q = quals[i][j]
      mean_arr[j] = mean_arr[j] + q
      if q < min_arr[j] {
        min_arr[j] = q
      }
      if q > max_arr[j] {
        max_arr[j] = q
      }
      j = j + 1
    }
    i = i + 1
  }
  let n = quals.length().to_double()
  let mut j = 0
  while j < n_pos {
    mean_arr[j] = mean_arr[j] / n
    // Collect and sort per-position for quartiles
    let col : Array[Double] = []
    let mut k = 0
    while k < quals.length() {
      col.push(quals[k][j])
      k = k + 1
    }
    col.sort_by(fn(a, b) {
      if a < b {
        -1
      } else if a > b {
        1
      } else {
        0
      }
    })
    let idx25 = (col.length() - 1).to_double() * 0.25
    let idx75 = (col.length() - 1).to_double() * 0.75
    let lo25 = idx25.floor().to_int()
    let hi25 = idx25.ceil().to_int()
    let lo75 = idx75.floor().to_int()
    let hi75 = idx75.ceil().to_int()
    if lo25 == hi25 {
      q25_arr[j] = col[lo25]
    } else {
      q25_arr[j] = col[lo25] + (col[hi25] - col[lo25]) * (idx25 - lo25.to_double())
    }
    if lo75 == hi75 {
      q75_arr[j] = col[lo75]
    } else {
      q75_arr[j] = col[lo75] + (col[hi75] - col[lo75]) * (idx75 - lo75.to_double())
    }
    j = j + 1
  }
  QualityProfile::{ positions: n_pos, mean: mean_arr, min: min_arr, max: max_arr, q25: q25_arr, q75: q75_arr }
}

///|
/// Error model: 16 x (max_qual+1) transition rate matrix.
/// err_mat[transition][qual] = estimated error rate.
/// Transitions: AA(0), AC(1), AG(2), AT(3), CA(4), CC(5), CG(6), CT(7),
///              GA(8), GC(9), GG(10), GT(11), TA(12), TC(13), TG(14), TT(15)
pub struct ErrorModel {
  err_mat : Array[Array[Double]] // 16 x (max_qual+1)
  max_qual : Int
  n_reads : Int
}

///|
/// Get transition index from original nucleotide (0-3) to observed nucleotide (0-3).
pub fn transition_index(from_nt : Int, to_nt : Int) -> Int {
  from_nt * 4 + to_nt
}

///|
/// Initialize error model with default rates (from dada2 loessErr initial estimates).
/// Uses 41 quality levels (0-40).
pub fn init_error_model(max_qual : Int) -> ErrorModel {
  let err_mat : Array[Array[Double]] = []
  let mut t = 0
  while t < 16 {
    let row : Array[Double] = Array::make(max_qual + 1, 0.0)
    let from_nt = t / 4
    let to_nt = t % 4
    let mut q = 0
    while q <= max_qual {
      if from_nt == to_nt {
        row[q] = 1.0 - 1.0e-3
      } else {
        // Initial estimate: ~10^(-0.1*qual) for transitions
        let rate = @math.pow(10.0, -0.1 * q.to_double() - 2.0)
        if rate > 0.333 {
          row[q] = 0.333
        } else {
          row[q] = rate
        }
      }
      q = q + 1
    }
    err_mat.push(row)
    t = t + 1
  }
  ErrorModel::{ err_mat, max_qual, n_reads: 0 }
}

///|
/// loessErr: DADA2's parametric error rate fitting.
/// Given observed error counts, fit error rates as a function of quality.
/// Returns updated error model.
/// err_counts[t][q] = (n_errors, n_total) for transition t at quality q
pub fn loess_err(
  err_counts : Array[Array[(Int, Int)]],
  max_qual : Int,
) -> ErrorModel {
  let err_mat : Array[Array[Double]] = []
  let mut t = 0
  while t < 16 {
    let row : Array[Double] = Array::make(max_qual + 1, 0.0)
    let from_nt = t / 4
    let to_nt = t % 4
    let mut q = 0
    while q <= max_qual {
      let (n_err, n_total) = err_counts[t][q]
      if from_nt == to_nt {
        if n_total > 0 {
          row[q] = 1.0 - n_err.to_double() / n_total.to_double()
        } else {
          row[q] = 1.0 - 1.0e-3
        }
      } else {
        if n_total > 0 {
          let rate = n_err.to_double() / n_total.to_double()
          // Apply monotone constraint: error rate must decrease with quality
          if q > 0 && rate > row[q - 1] {
            row[q] = row[q - 1]
          } else {
            row[q] = rate
          }
        } else {
          // Interpolate from neighbors
          if q > 0 {
            row[q] = row[q - 1] * 0.5
          } else {
            row[q] = 1.0e-3
          }
        }
      }
      q = q + 1
    }
    err_mat.push(row)
    t = t + 1
  }
  ErrorModel::{ err_mat, max_qual, n_reads: 0 }
}

///|
/// Learn error model from data (simplified: one round of self-consistent learning).
/// Given dereplicated sequences, estimate error rates.
pub fn learn_errors(derep : Dereplication, max_qual : Int) -> ErrorModel {
  // Count observed errors: assume most abundant sequences are "true"
  // and less abundant ones within small Hamming distance are errors
  let err_counts : Array[Array[(Int, Int)]] = []
  let mut t = 0
  while t < 16 {
    err_counts.push(Array::make(max_qual + 1, (0, 0)))
    t = t + 1
  }
  let n_uniques = derep.uniques.length()
  let mut i = 0
  while i < n_uniques {
    let seq_i = derep.uniques[i]
    // Compare to more abundant sequences (potential "true" sequences)
    let mut j = 0
    while j < i {
      let seq_j = derep.uniques[j]
      let dist = dada_hamming_distance(seq_i.sequence, seq_j.sequence)
      if dist >= 1 && dist <= 3 {
        // Count this as potential error of seq_j
        let a_bytes = seq_j.sequence.to_array()
        let b_bytes = seq_i.sequence.to_array()
        let mut pos = 0
        while pos < a_bytes.length() {
          let from_nt = nt_to_int(a_bytes[pos])
          let to_nt = nt_to_int(b_bytes[pos])
          if from_nt >= 0 && to_nt >= 0 && from_nt != to_nt {
            let ti = transition_index(from_nt, to_nt)
            let qual = if seq_i.qual_scores.length() > pos {
              seq_i.qual_scores[pos].round().to_int()
            } else {
              20
            }
            let q_clamped = if qual > max_qual { max_qual } else { qual }
            let (ne, nt_c) = err_counts[ti][q_clamped]
            err_counts[ti][q_clamped] = (ne + seq_i.abundance, nt_c + seq_i.abundance)
          }
          pos = pos + 1
        }
        // Also count matches
        let mut pos2 = 0
        while pos2 < a_bytes.length() {
          let from_nt = nt_to_int(a_bytes[pos2])
          let to_nt = nt_to_int(b_bytes[pos2])
          if from_nt >= 0 && to_nt >= 0 && from_nt == to_nt {
            let ti = transition_index(from_nt, to_nt)
            let qual = if seq_i.qual_scores.length() > pos2 {
              seq_i.qual_scores[pos2].round().to_int()
            } else {
              20
            }
            let q_clamped = if qual > max_qual { max_qual } else { qual }
            let (ne, nt_c) = err_counts[ti][q_clamped]
            err_counts[ti][q_clamped] = (ne, nt_c + seq_i.abundance)
          }
          pos2 = pos2 + 1
        }
        break // Only count against the most abundant "true" sequence
      }
      j = j + 1
    }
    i = i + 1
  }
  let model = loess_err(err_counts, max_qual)
  ErrorModel::{ err_mat: model.err_mat, max_qual, n_reads: derep.total_reads }
}

///|
/// DADA partition: a cluster of sequences centered on a "center" sequence.
pub struct Partition {
  center : String
  abundance : Int
  members : Array[(String, Int)] // (sequence, abundance) of all members
}

///|
/// DADA inference result.
pub struct DadaResult {
  partitions : Array[Partition]
  error_model : ErrorModel
  n_clusters : Int
  n_input_reads : Int
}

///|
/// Compute expected abundance of an error sequence given the true sequence and error model.
/// Uses Poisson model: E[abund_error] = abund_true * prod(err_rate at each error position)
pub fn expected_abundance(
  true_seq : String,
  error_seq : String,
  true_abund : Int,
  model : ErrorModel,
) -> Double {
  let a_bytes = true_seq.to_array()
  let b_bytes = error_seq.to_array()
  let mut log_prob = 0.0
  let mut pos = 0
  while pos < a_bytes.length() {
    let from_nt = nt_to_int(a_bytes[pos])
    let to_nt = nt_to_int(b_bytes[pos])
    if from_nt >= 0 && to_nt >= 0 {
      let ti = transition_index(from_nt, to_nt)
      let q = 20 // default quality
      let q_clamped = if q > model.max_qual { model.max_qual } else { q }
      let rate = model.err_mat[ti][q_clamped]
      if rate > 0.0 {
        log_prob = log_prob + @math.ln(rate)
      } else {
        log_prob = log_prob + (-50.0)
      }
    }
    pos = pos + 1
  }
  true_abund.to_double() * @math.exp(log_prob)
}

///|
/// P-value test: is the observed abundance of a sequence consistent with being an error?
/// Uses the abundance p-value test from the DADA2 paper.
/// If p < threshold, the sequence is NOT an error (it's a real variant).
pub fn abundance_p_value(
  observed : Int,
  expected : Double,
) -> Double {
  if expected <= 0.0 {
    return 0.0
  }
  // P(X >= observed | Poisson(expected))
  // = 1 - CDF(observed - 1) = 1 - e^{-lambda} * sum_{k=0}^{observed-1} lambda^k / k!
  let lambda = expected
  let mut cdf = 0.0
  let mut log_term = -lambda // k=0 term: e^{-lambda}
  let mut k = 0
  while k < observed {
    cdf = cdf + @math.exp(log_term)
    k = k + 1
    log_term = log_term + @math.ln(lambda) - @math.ln(k.to_double())
  }
  let p = 1.0 - cdf
  if p < 0.0 {
    0.0
  } else if p > 1.0 {
    1.0
  } else {
    p
  }
}

///|
/// Run DADA inference on dereplicated sequences.
/// Core algorithm: greedily assign sequences to partitions, test for errors.
/// OMEGA_A: p-value threshold (default 1e-40 in dada2)
pub fn dada_inference(
  derep : Dereplication,
  model : ErrorModel,
  omega_a~ : Double = 1.0e-40,
) -> DadaResult {
  let uniques = derep.uniques
  if uniques.length() == 0 {
    return DadaResult::{ partitions: [], error_model: model, n_clusters: 0, n_input_reads: 0 }
  }
  // Initialize: most abundant sequence is first partition center
  let partitions : Array[Partition] = []
  let assigned : Array[Bool] = Array::make(uniques.length(), false)
  let mut i = 0
  while i < uniques.length() {
    if !assigned[i] {
      let center = uniques[i]
      assigned[i] = true
      let members : Array[(String, Int)] = [(center.sequence, center.abundance)]
      // Try to assign remaining sequences to this partition
      let mut j = i + 1
      while j < uniques.length() {
        if !assigned[j] {
          let candidate = uniques[j]
          let dist = dada_hamming_distance(center.sequence, candidate.sequence)
          if dist >= 1 && dist <= 4 {
            let exp_abund = expected_abundance(center.sequence, candidate.sequence, center.abundance, model)
            let pval = abundance_p_value(candidate.abundance, exp_abund)
            // If p-value is large, it's consistent with being an error
            if pval > omega_a {
              members.push((candidate.sequence, candidate.abundance))
              assigned[j] = true
            }
          }
        }
        j = j + 1
      }
      // Compute total abundance
      let mut total_abund = 0
      let mut mi = 0
      while mi < members.length() {
        total_abund = total_abund + members[mi].1
        mi = mi + 1
      }
      partitions.push(Partition::{ center: center.sequence, abundance: total_abund, members })
    }
    i = i + 1
  }
  DadaResult::{
    partitions,
    error_model: model,
    n_clusters: partitions.length(),
    n_input_reads: derep.total_reads,
  }
}

///|
/// Merged paired-end read.
pub struct MergedPair {
  sequence : String
  quality : Array[Double]
  forward : String
  reverse : String
  overlap : Int
  mismatch : Int
  accepted : Bool
}

///|
/// Reverse complement of a DNA sequence.
pub fn dada_reverse_complement(seq : String) -> String {
  let bytes = seq.to_array()
  let n = bytes.length()
  let sb = StringBuilder::new()
  let mut i = 0
  while i < n {
    let c = bytes[n - 1 - i]
    let comp = match c.to_int() {
      65 | 97 => 84 // A->T
      84 | 116 => 65 // T->A
      67 | 99 => 71 // C->G
      71 | 103 => 67 // G->C
      _ => 78 // N->N
    }
    sb.write_char(comp.unsafe_to_char())
    i = i + 1
  }
  sb.to_string()
}

///|
/// Merge paired-end reads by finding overlap.
/// forward: 5'→3' forward read
/// reverse: 5'→3' reverse read (will be reverse-complemented)
/// min_overlap: minimum overlap length (default 12)
/// max_mismatch_rate: max allowed mismatch rate in overlap (default 0.1)
pub fn merge_pair(
  forward : String,
  reverse : String,
  min_overlap~ : Int = 12,
  max_mismatch_rate~ : Double = 0.1,
) -> MergedPair {
  let rc_reverse = dada_reverse_complement(reverse)
  let fwd_len = forward.length()
  let rev_len = rc_reverse.length()
  let max_over = if fwd_len < rev_len { fwd_len } else { rev_len }
  let mut best_overlap = -1
  let mut best_mismatch = 999999
  let mut best_start = -1
  // Try different overlap lengths
  let mut overlap = max_over
  while overlap >= min_overlap {
    // forward[...forward_len - overlap:] should match rc_reverse[:overlap]
    let fwd_start = fwd_len - overlap
    let mut mismatches = 0
    let mut k = 0
    while k < overlap {
      let fc = forward[fwd_start + k]
      let rc = rc_reverse[k]
      if fc != rc {
        mismatches = mismatches + 1
      }
      k = k + 1
    }
    let max_mm = (overlap.to_double() * max_mismatch_rate).floor().to_int()
    if mismatches <= max_mm {
      best_overlap = overlap
      best_mismatch = mismatches
      best_start = fwd_start
      break
    }
    overlap = overlap - 1
  }
  if best_overlap < 0 {
    return MergedPair::{
      sequence: "",
      quality: [],
      forward,
      reverse: rc_reverse,
      overlap: 0,
      mismatch: 0,
      accepted: false,
    }
  }
  // Build merged sequence
  let sb = StringBuilder::new()
  // Part from forward only
  let mut i = 0
  while i < best_start {
    sb.write_char(forward[i].unsafe_to_char())
    i = i + 1
  }
  // Overlap region: use forward (or could use higher quality)
  let mut k = 0
  while k < best_overlap {
    sb.write_char(forward[best_start + k].unsafe_to_char())
    k = k + 1
  }
  // Part from rc_reverse only (after overlap)
  let mut j = best_overlap
  while j < rev_len {
    sb.write_char(rc_reverse[j].unsafe_to_char())
    j = j + 1
  }
  let merged = sb.to_string()
  let qual : Array[Double] = Array::make(merged.length(), 30.0)
  MergedPair::{
    sequence: merged,
    quality: qual,
    forward,
    reverse: rc_reverse,
    overlap: best_overlap,
    mismatch: best_mismatch,
    accepted: true,
  }
}

///|
/// Merge multiple paired-end reads.
pub fn merge_pairs(
  forwards : Array[String],
  reverses : Array[String],
  min_overlap~ : Int = 12,
  max_mismatch_rate~ : Double = 0.1,
) -> Array[MergedPair] {
  let results : Array[MergedPair] = []
  let n = if forwards.length() < reverses.length() { forwards.length() } else { reverses.length() }
  let mut i = 0
  while i < n {
    results.push(merge_pair(forwards[i], reverses[i], min_overlap=min_overlap, max_mismatch_rate=max_mismatch_rate))
    i = i + 1
  }
  results
}

///|
/// Check if a sequence is a bimera (chimeric) of two parent sequences.
/// A bimera has a breakpoint where it switches from matching parent A to parent B.
pub fn is_bimera(
  seq : String,
  parent_a : String,
  parent_b : String,
  min_one_parent_frac~ : Double = 0.8,
) -> Bool {
  if seq.length() != parent_a.length() || seq.length() != parent_b.length() {
    return false
  }
  let n = seq.length()
  let s_bytes = seq.to_array()
  let a_bytes = parent_a.to_array()
  let b_bytes = parent_b.to_array()
  // Try all possible breakpoints
  let mut bp = 1
  while bp < n {
    // Try: first part from A, second from B
    let mut match_a_first = 0
    let mut match_b_second = 0
    let mut k = 0
    while k < bp {
      if s_bytes[k] == a_bytes[k] {
        match_a_first = match_a_first + 1
      }
      k = k + 1
    }
    while k < n {
      if s_bytes[k] == b_bytes[k] {
        match_b_second = match_b_second + 1
      }
      k = k + 1
    }
    let frac_a = match_a_first.to_double() / bp.to_double()
    let frac_b = match_b_second.to_double() / (n - bp).to_double()
    if frac_a >= min_one_parent_frac && frac_b >= min_one_parent_frac {
      return true
    }
    // Try: first part from B, second from A
    let mut match_b_first = 0
    let mut match_a_second = 0
    let mut k2 = 0
    while k2 < bp {
      if s_bytes[k2] == b_bytes[k2] {
        match_b_first = match_b_first + 1
      }
      k2 = k2 + 1
    }
    while k2 < n {
      if s_bytes[k2] == a_bytes[k2] {
        match_a_second = match_a_second + 1
      }
      k2 = k2 + 1
    }
    let frac_b2 = match_b_first.to_double() / bp.to_double()
    let frac_a2 = match_a_second.to_double() / (n - bp).to_double()
    if frac_b2 >= min_one_parent_frac && frac_a2 >= min_one_parent_frac {
      return true
    }
    bp = bp + 1
  }
  false
}

///|
/// Remove chimeras (bimeras) from a set of ASVs.
/// Returns (non-chimeric ASVs, non-chimeric abundances, chimera flags).
pub fn remove_chimeras(
  asvs : Array[String],
  abundances : Array[Int],
) -> (Array[String], Array[Int], Array[Bool]) {
  let n = asvs.length()
  let is_chimera : Array[Bool] = Array::make(n, false)
  // Sort indices by abundance descending for parent candidates
  let indices : Array[Int] = []
  let mut i = 0
  while i < n {
    indices.push(i)
    i = i + 1
  }
  indices.sort_by(fn(a, b) {
    if abundances[b] > abundances[a] {
      1
    } else if abundances[b] < abundances[a] {
      -1
    } else {
      0
    }
  })
  // Check each ASV against more abundant parents
  let mut idx = 0
  while idx < n {
    let seq_idx = indices[idx]
    let seq = asvs[seq_idx]
    // Parents are more abundant sequences
    let mut p1 = 0
    while p1 < idx {
      let parent_a_idx = indices[p1]
      let parent_a = asvs[parent_a_idx]
      if parent_a.length() == seq.length() {
        let mut p2 = p1 + 1
        while p2 < idx {
          let parent_b_idx = indices[p2]
          let parent_b = asvs[parent_b_idx]
          if parent_b.length() == seq.length() {
            if is_bimera(seq, parent_a, parent_b) {
              is_chimera[seq_idx] = true
            }
          }
          p2 = p2 + 1
        }
      }
      p1 = p1 + 1
    }
    idx = idx + 1
  }
  let clean_asvs : Array[String] = []
  let clean_abund : Array[Int] = []
  let mut j = 0
  while j < n {
    if !is_chimera[j] {
      clean_asvs.push(asvs[j])
      clean_abund.push(abundances[j])
    }
    j = j + 1
  }
  (clean_asvs, clean_abund, is_chimera)
}

///|
/// Sequence table: ASVs x samples count matrix.
pub struct SequenceTable {
  asv_sequences : Array[String]
  sample_names : Array[String]
  counts : Array[Array[Int]] // [asv][sample]
  n_asvs : Int
  n_samples : Int
}

///|
/// Make sequence table from per-sample DADA results.
/// samples: array of (sample_name, DadaResult)
pub fn make_sequence_table(samples : Array[(String, DadaResult)]) -> SequenceTable {
  // Collect all unique ASV sequences across samples
  let asv_set : Map[String, Int] = Map([], capacity=256)
  let sample_names : Array[String] = []
  let mut si = 0
  while si < samples.length() {
    sample_names.push(samples[si].0)
    let result = samples[si].1
    let mut pi = 0
    while pi < result.partitions.length() {
      let center = result.partitions[pi].center
      if !asv_set.contains(center) {
        asv_set[center] = asv_set.length()
      }
      pi = pi + 1
    }
    si = si + 1
  }
  let asv_sequences : Array[String] = Array::make(asv_set.length(), "")
  let keys = asv_set.keys().collect()
  let mut ki = 0
  while ki < keys.length() {
    asv_sequences[asv_set[keys[ki]]] = keys[ki]
    ki = ki + 1
  }
  let n_asvs = asv_sequences.length()
  let n_samples = sample_names.length()
  let counts : Array[Array[Int]] = []
  let mut ai = 0
  while ai < n_asvs {
    counts.push(Array::make(n_samples, 0))
    ai = ai + 1
  }
  // Fill counts
  let mut sj = 0
  while sj < n_samples {
    let result = samples[sj].1
    let mut pi = 0
    while pi < result.partitions.length() {
      let center = result.partitions[pi].center
      let abund = result.partitions[pi].abundance
      if asv_set.contains(center) {
        counts[asv_set[center]][sj] = abund
      }
      pi = pi + 1
    }
    sj = sj + 1
  }
  SequenceTable::{ asv_sequences, sample_names, counts, n_asvs, n_samples }
}

///|
/// Get total counts per sample from sequence table.
pub fn seq_table_sample_totals(table : SequenceTable) -> Array[Int] {
  let totals : Array[Int] = Array::make(table.n_samples, 0)
  let mut ai = 0
  while ai < table.n_asvs {
    let mut si = 0
    while si < table.n_samples {
      totals[si] = totals[si] + table.counts[ai][si]
      si = si + 1
    }
    ai = ai + 1
  }
  totals
}

///|
/// Get total counts per ASV from sequence table.
pub fn seq_table_asv_totals(table : SequenceTable) -> Array[Int] {
  let totals : Array[Int] = Array::make(table.n_asvs, 0)
  let mut ai = 0
  while ai < table.n_asvs {
    let mut si = 0
    while si < table.n_samples {
      totals[ai] = totals[ai] + table.counts[ai][si]
      si = si + 1
    }
    ai = ai + 1
  }
  totals
}

///|
/// Taxonomy assignment result.
pub struct TaxonomyResult {
  asv : String
  taxonomy : Array[String] // [kingdom, phylum, class, order, family, genus, species]
  confidence : Double
}

///|
/// Reference sequence for taxonomy assignment.
pub struct RefEntry {
  sequence : String
  taxonomy : Array[String]
}

///|
/// Assign taxonomy using naive Bayesian classifier (simplified k-mer approach).
/// ref_seqs: reference database entries
/// k: k-mer size (default 8)
pub fn assign_taxonomy(
  asvs : Array[String],
  ref_seqs : Array[RefEntry],
  k~ : Int = 8,
) -> Array[TaxonomyResult] {
  let results : Array[TaxonomyResult] = []
  let n_levels = if ref_seqs.length() > 0 { ref_seqs[0].taxonomy.length() } else { 0 }
  // Build k-mer index for each reference sequence
  let ref_kmers : Array[Map[String, Int]] = []
  let mut ri = 0
  while ri < ref_seqs.length() {
    let kmers = extract_kmers(ref_seqs[ri].sequence, k)
    let km : Map[String, Int] = Map([], capacity=kmers.length())
    let mut ki = 0
    while ki < kmers.length() {
      let key = kmers[ki]
      if km.contains(key) {
        km[key] = km[key] + 1
      } else {
        km[key] = 1
      }
      ki = ki + 1
    }
    ref_kmers.push(km)
    ri = ri + 1
  }
  // For each ASV, find best matching reference
  let mut ai = 0
  while ai < asvs.length() {
    let asv = asvs[ai]
    let asv_kmers = extract_kmers(asv, k)
    let asv_kmer_set : Map[String, Int] = Map([], capacity=asv_kmers.length())
    let mut ki = 0
    while ki < asv_kmers.length() {
      let key = asv_kmers[ki]
      if asv_kmer_set.contains(key) {
        asv_kmer_set[key] = asv_kmer_set[key] + 1
      } else {
        asv_kmer_set[key] = 1
      }
      ki = ki + 1
    }
    // Compute k-mer overlap with each reference
    let mut best_score = -1.0
    let mut best_tax : Array[String] = []
    let mut rj = 0
    while rj < ref_seqs.length() {
      let mut score = 0.0
      let asv_keys = asv_kmer_set.keys().collect()
      let mut aki = 0
      while aki < asv_keys.length() {
        let key = asv_keys[aki]
        if ref_kmers[rj].contains(key) {
          score = score + 1.0
        }
        aki = aki + 1
      }
      // Normalize by unique ASV k-mer count
      let n_unique = asv_kmer_set.length()
      if n_unique > 0 {
        score = score / n_unique.to_double()
      }
      if score > best_score {
        best_score = score
        best_tax = []
        let mut ti = 0
        while ti < ref_seqs[rj].taxonomy.length() {
          best_tax.push(ref_seqs[rj].taxonomy[ti])
          ti = ti + 1
        }
      }
      rj = rj + 1
    }
    // Bootstrap confidence: fraction of k-mers that agree with best match
    let confidence = if best_score > 0.0 { best_score } else { 0.0 }
    let tax = if best_tax.length() == 0 {
      let dummy : Array[String] = []
      let mut t = 0
      while t < n_levels {
        dummy.push("Unknown")
        t = t + 1
      }
      dummy
    } else {
      best_tax
    }
    results.push(TaxonomyResult::{ asv, taxonomy: tax, confidence })
    ai = ai + 1
  }
  results
}

///|
/// Extract k-mers from a sequence.
pub fn extract_kmers(seq : String, k : Int) -> Array[String] {
  let kmers : Array[String] = []
  if seq.length() < k {
    return kmers
  }
  let n = seq.length() - k + 1
  let mut i = 0
  while i < n {
    let sb = StringBuilder::new()
    let mut j = 0
    while j < k {
      sb.write_char(seq[i + j].unsafe_to_char())
      j = j + 1
    }
    kmers.push(sb.to_string())
    i = i + 1
  }
  kmers
}

///|
/// Complete dada2 pipeline: from raw sequences to ASV table with taxonomy.
pub struct PipelineResult {
  seq_table : SequenceTable
  taxonomy : Array[TaxonomyResult]
  n_chimeras_removed : Int
  error_model : ErrorModel
}

///|
/// Run complete dada2 pipeline on paired-end reads.
/// samples: array of (sample_name, forward_reads, reverse_reads)
pub fn run_pipeline(
  samples : Array[(String, Array[String], Array[String])],
  ref_seqs : Array[RefEntry],
  max_qual~ : Int = 40,
  omega_a~ : Double = 1.0e-40,
  min_overlap~ : Int = 12,
) -> PipelineResult {
  // Step 1: Merge pairs and collect per-sample ASVs
  let dada_results : Array[(String, DadaResult)] = []
  let mut si = 0
  while si < samples.length() {
    let (name, fwd, rev) = samples[si]
    // Merge pairs
    let merged = merge_pairs(fwd, rev, min_overlap=min_overlap)
    let good_seqs : Array[String] = []
    let mut mi = 0
    while mi < merged.length() {
      if merged[mi].accepted {
        good_seqs.push(merged[mi].sequence)
      }
      mi = mi + 1
    }
    // Dereplicate
    let derep = dereplicate(good_seqs, [])
    // Learn errors
    let model = learn_errors(derep, max_qual)
    // DADA inference
    let dada = dada_inference(derep, model, omega_a=omega_a)
    dada_results.push((name, dada))
    si = si + 1
  }
  // Step 2: Make sequence table
  let table = make_sequence_table(dada_results)
  // Step 3: Remove chimeras
  let asv_totals = seq_table_asv_totals(table)
  let (clean_asvs, _, chimera_flags) = remove_chimeras(table.asv_sequences, asv_totals)
  let mut n_chim = 0
  let mut ci = 0
  while ci < chimera_flags.length() {
    if chimera_flags[ci] {
      n_chim = n_chim + 1
    }
    ci = ci + 1
  }
  // Rebuild table without chimeras
  let clean_table = filter_table(table, chimera_flags)
  // Step 4: Assign taxonomy
  let taxonomy = assign_taxonomy(clean_asvs, ref_seqs)
  PipelineResult::{
    seq_table: clean_table,
    taxonomy,
    n_chimeras_removed: n_chim,
    error_model: if dada_results.length() > 0 {
      dada_results[0].1.error_model
    } else {
      init_error_model(max_qual)
    },
  }
}

///|
/// Filter sequence table by chimera flags (keep where flag is false).
pub fn filter_table(table : SequenceTable, flags : Array[Bool]) -> SequenceTable {
  let asv_sequences : Array[String] = []
  let counts : Array[Array[Int]] = []
  let mut i = 0
  while i < table.n_asvs {
    if !flags[i] {
      asv_sequences.push(table.asv_sequences[i])
      counts.push(table.counts[i])
    }
    i = i + 1
  }
  SequenceTable::{
    asv_sequences,
    sample_names: table.sample_names,
    counts,
    n_asvs: asv_sequences.length(),
    n_samples: table.n_samples,
  }
}

///|
/// Track abundance: convert sequence table to relative abundance.
pub fn seq_table_to_relative(table : SequenceTable) -> Array[Array[Double]] {
  let totals = seq_table_sample_totals(table)
  let rel : Array[Array[Double]] = []
  let mut ai = 0
  while ai < table.n_asvs {
    let row : Array[Double] = Array::make(table.n_samples, 0.0)
    let mut si = 0
    while si < table.n_samples {
      if totals[si] > 0 {
        row[si] = table.counts[ai][si].to_double() / totals[si].to_double()
      }
      si = si + 1
    }
    rel.push(row)
    ai = ai + 1
  }
  rel
}

///|
/// Create example data for testing/demo.
pub fn create_example_derep() -> Dereplication {
  let seqs = [
    "ACGTACGTACGTACGT",
    "ACGTACGTACGTACGT",
    "ACGTACGTACGTACGT",
    "ACGTACGTACGTACGA",
    "ACGTACGTACGTACGA",
    "ACGTACGTACGTACGT",
    "TTGCGCGCATATATAT",
    "TTGCGCGCATATATAT",
    "TTGCGCGCATATATAT",
    "TTGCGCGCATATATAT",
    "TTGCGCGCATATATAC",
    "ACGTACGTACGTACGT",
  ]
  dereplicate(seqs, [])
}

///|
/// Create example reference database for taxonomy.
pub fn create_example_ref_db() -> Array[RefEntry] {
  [
    RefEntry::{ sequence: "ACGTACGTACGTACGTACGTACGTACGTACGT", taxonomy: ["Bacteria", "Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Escherichia", "coli"] },
    RefEntry::{ sequence: "TTGCGCGCATATATATGCGCGCGCATATATA", taxonomy: ["Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus", "plantarum"] },
    RefEntry::{ sequence: "GCGCATATATATGCGCGCGCATATATATGCG", taxonomy: ["Bacteria", "Actinobacteria", "Actinomycetia", "Bifidobacteriales", "Bifidobacteriaceae", "Bifidobacterium", "longum"] },
  ]
}
