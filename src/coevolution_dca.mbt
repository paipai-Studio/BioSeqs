// ===========================================================================
// coevolution_dca.mbt — MSA Coevolution & Direct Coupling Analysis (DCA)
//
// Algorithms:
//   1. Mutual Information (MI) between alignment columns
//   2. MI with Average Product Correction (MI-APC)
//   3. Mean-Field DCA (mfDCA): infer Potts model couplings via covariance
//      inversion, score residue pairs by Frobenius norm and Direct Information
//   4. Contact prediction from DCA scores
//
// References:
//   - Morcos et al., PNAS 2011 (mfDCA / DCA)
//   - Dunn et al., Bioinformatics 2008 (APC correction)
//   - Jones et al., FEBS Lett 2012 (PSICOV, similar mean-field approach)
// ===========================================================================

///|
/// Configuration for DCA computation.
pub(all) struct DCAConfig {
  /// Alphabet: "protein" (21 states: 20 AA + gap) or "dna" (5 states: 4 NT + gap)
  alphabet : String
  /// Pseudocount weight (typically 0.01–0.1 for mfDCA)
  pseudocount : Double
  /// Reweight threshold: sequences with > threshold * L identical positions
  /// are clustered. Typically 0.8.
  reweight_threshold : Double
}

///|
/// Result of DCA computation.
pub(all) struct DCAResult {
  /// L x L MI matrix
  mi_scores : Array[Array[Double]]
  /// L x L MI-APC matrix
  mi_apc_scores : Array[Array[Double]]
  /// L x L DCA Frobenius norm scores
  dca_scores : Array[Array[Double]]
  /// L x L Direct Information scores
  di_scores : Array[Array[Double]]
  /// Effective number of sequences after reweighting
  m_eff : Double
  /// Alignment length
  align_len : Int
  /// Number of sequences
  num_seqs : Int
}

///|
/// A predicted residue-residue contact.
pub(all) struct Contact {
  pos_i : Int
  pos_j : Int
  score : Double
}

///|
/// Default config for protein DCA.
pub fn DCAConfig::default_protein() -> DCAConfig {
  DCAConfig::{
    alphabet: "protein",
    pseudocount: 0.01,
    reweight_threshold: 0.8,
  }
}

///|
/// Default config for DNA DCA.
pub fn DCAConfig::default_dna() -> DCAConfig {
  DCAConfig::{
    alphabet: "dna",
    pseudocount: 0.01,
    reweight_threshold: 0.8,
  }
}

// ---------------------------------------------------------------------------
// Alphabet handling
// ---------------------------------------------------------------------------

///|
/// Get the alphabet string for a given type.

///|
/// Square root helper.
fn coev_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut guess = if x > 1.0 { x } else { 1.0 }
  let mut prev = 0.0
  for i = 0; i < 50; i = i + 1 {
    prev = guess
    guess = (guess + x / guess) * 0.5
    if (guess - prev).abs() < 1.0e-15 * guess {
      break
    }
  }
  guess
}

fn get_alphabet(alphabet : String) -> String {
  match alphabet {
    "protein" => "ACDEFGHIKLMNPQRSTVWY-"
    "dna" => "ACGT-"
    _ => alphabet
  }
}

///|
/// Get the number of states q for an alphabet.
fn get_q(alphabet : String) -> Int {
  get_alphabet(alphabet).length()
}

///|
/// Encode a character to an index in the alphabet.
/// Returns -1 if the character is not found.
fn encode_char(c : Char, alphabet : String) -> Int {
  let ab = get_alphabet(alphabet)
  let s = ab.to_array()
  for i in 0..<s.length() {
    if s[i] == c {
      return i
    }
  }
  // Unknown character → map to gap (last state)
  s.length() - 1
}

///|
/// Encode an alignment (array of sequences) into integer codes.
/// Returns a 2D array: encoded[seq_idx][pos_idx] = state index.
pub fn encode_alignment(
  alignment : Array[String],
  alphabet : String,
) -> Array[Array[Int]] {
  let result = Array(capacity=alignment.length())
  for seq in alignment {
    let encoded = Array(capacity=seq.length())
    for i in 0..<seq.length() {
      let c = seq.unsafe_get(i).unsafe_to_char()
      encoded.push(encode_char(c, alphabet))
    }
    result.push(encoded)
  }
  result
}

// ---------------------------------------------------------------------------
// Sequence reweighting
// ---------------------------------------------------------------------------

///|
/// Compute sequence weights for an alignment.
/// Sequences with more than threshold * L identical positions are clustered.
/// Each sequence gets weight 1 / cluster_size.
pub fn sequence_weights(
  alignment : Array[String],
  threshold : Double,
) -> Array[Double] {
  let n = alignment.length()
  if n == 0 {
    return []
  }
  let l = alignment[0].length()
  let threshold_count = (threshold * l.to_double()).to_int()
  let weights = Array::make(n, 1.0)
  // For each pair, count identical positions
  for i in 0..<n {
    let mut count = 1
    for j in (i + 1)..<n {
      let mut identical = 0
      for k in 0..<l {
        if alignment[i].unsafe_get(k) == alignment[j].unsafe_get(k) {
          identical = identical + 1
        }
      }
      if identical >= threshold_count {
        count = count + 1
      }
    }
    weights[i] = 1.0 / count.to_double()
  }
  weights
}

///|
/// Compute the effective number of sequences from weights.
pub fn effective_sequences(weights : Array[Double]) -> Double {
  let mut sum = 0.0
  for w in weights {
    sum = sum + w
  }
  sum
}

// ---------------------------------------------------------------------------
// Single-site and pair frequencies
// ---------------------------------------------------------------------------

///|
/// Compute single-site frequencies f_i(a) with pseudocount.
/// Returns a 2D array: freq[i][a] = frequency of state a at position i.
pub fn single_site_frequencies(
  encoded : Array[Array[Int]],
  q : Int,
  pseudocount : Double,
  weights : Array[Double],
  m_eff : Double,
) -> Array[Array[Double]] {
  let n = encoded.length()
  let l = if n > 0 { encoded[0].length() } else { 0 }
  let result = Array(capacity=l)
  // Regularization parameter
  let lambda = pseudocount
  let one_minus_lambda = 1.0 - lambda
  // Uniform background
  let uniform = 1.0 / q.to_double()
  for i in 0..<l {
    let counts = Array::make(q, 0.0)
    for s in 0..<n {
      let a = encoded[s][i]
      counts[a] = counts[a] + weights[s]
    }
    let freq = Array(capacity=q)
    for a in 0..<q {
      let empirical = counts[a] / m_eff
      freq.push(one_minus_lambda * empirical + lambda * uniform)
    }
    result.push(freq)
  }
  result
}

///|
/// Compute pair frequencies f_ij(a,b) with pseudocount.
/// Returns a 4D array: freq[i][j][a][b].
pub fn pair_frequencies(
  encoded : Array[Array[Int]],
  q : Int,
  pseudocount : Double,
  weights : Array[Double],
  m_eff : Double,
) -> Array[Array[Array[Array[Double]]]] {
  let n = encoded.length()
  let l = if n > 0 { encoded[0].length() } else { 0 }
  let lambda = pseudocount
  let one_minus_lambda = 1.0 - lambda
  let uniform = 1.0 / (q.to_double() * q.to_double())
  let result = Array(capacity=l)
  for i in 0..<l {
    let row_i = Array(capacity=l)
    for j in 0..<l {
      let counts = Array::make(q * q, 0.0)
      for s in 0..<n {
        let a = encoded[s][i]
        let b = encoded[s][j]
        counts[a * q + b] = counts[a * q + b] + weights[s]
      }
      let freq_ij = Array(capacity=q)
      for a in 0..<q {
        let freq_a = Array(capacity=q)
        for b in 0..<q {
          let empirical = counts[a * q + b] / m_eff
          freq_a.push(one_minus_lambda * empirical + lambda * uniform)
        }
        freq_ij.push(freq_a)
      }
      row_i.push(freq_ij)
    }
    result.push(row_i)
  }
  result
}

// ---------------------------------------------------------------------------
// Mutual Information
// ---------------------------------------------------------------------------

///|
/// Compute the Mutual Information matrix for an alignment.
/// MI[i][j] = sum_a sum_b f_ij(a,b) * log(f_ij(a,b) / (f_i(a) * f_j(b)))
pub fn compute_mi(
  alignment : Array[String],
  alphabet : String,
) -> Array[Array[Double]] {
  let encoded = encode_alignment(alignment, alphabet)
  let q = get_q(alphabet)
  let n = encoded.length()
  let l = if n > 0 { encoded[0].length() } else { 0 }
  let weights = Array::make(n, 1.0)
  let m_eff = n.to_double()
  let fi = single_site_frequencies(encoded, q, 0.0, weights, m_eff)
  let fij = pair_frequencies(encoded, q, 0.0, weights, m_eff)
  let mi = Array(capacity=l)
  for i in 0..<l {
    let row = Array::make(l, 0.0)
    for j in 0..<l {
      if i == j {
        row[j] = 0.0
        continue
      }
      let mut s = 0.0
      for a in 0..<q {
        for b in 0..<q {
          let p = fij[i][j][a][b]
          if p > 0.0 && fi[i][a] > 0.0 && fi[j][b] > 0.0 {
            s = s + p * (@math.ln(p / (fi[i][a] * fi[j][b])))
          }
        }
      }
      row[j] = s
    }
    mi.push(row)
  }
  mi
}

// ---------------------------------------------------------------------------
// Average Product Correction (APC)
// ---------------------------------------------------------------------------

///|
/// Apply Average Product Correction to a score matrix.
/// APC[i][j] = S[i][j] - (row_mean[i] * col_mean[j] / matrix_mean)
pub fn apc_correction(scores : Array[Array[Double]]) -> Array[Array[Double]] {
  let l = scores.length()
  if l == 0 {
    return []
  }
  // Compute row means, col means, and overall mean
  let row_sums = Array::make(l, 0.0)
  let col_sums = Array::make(l, 0.0)
  let mut total = 0.0
  let mut count = 0
  for i in 0..<l {
    for j in 0..<l {
      if i != j {
        row_sums[i] = row_sums[i] + scores[i][j]
        col_sums[j] = col_sums[j] + scores[i][j]
        total = total + scores[i][j]
        count = count + 1
      }
    }
  }
  let matrix_mean = if count > 0 { total / count.to_double() } else { 0.0 }
  let result = Array(capacity=l)
  for i in 0..<l {
    let row = Array::make(l, 0.0)
    let row_mean = row_sums[i] / (l - 1).to_double()
    for j in 0..<l {
      if i == j {
        row[j] = 0.0
      } else {
        let col_mean = col_sums[j] / (l - 1).to_double()
        let apc = if matrix_mean > 0.0 {
          row_mean * col_mean / matrix_mean
        } else {
          0.0
        }
        row[j] = scores[i][j] - apc
      }
    }
    result.push(row)
  }
  result
}

///|
/// Compute MI with APC correction.
pub fn compute_mi_apc(
  alignment : Array[String],
  alphabet : String,
) -> Array[Array[Double]] {
  let mi = compute_mi(alignment, alphabet)
  apc_correction(mi)
}

// ---------------------------------------------------------------------------
// Mean-Field DCA (mfDCA)
// ---------------------------------------------------------------------------

///|
/// Compute the covariance matrix for mfDCA.
/// The covariance is C_{(i,a),(j,b)} = f_ij(a,b) - f_i(a) * f_j(b)
/// for a, b in 0..q-2 (excluding the last state, which is dependent).
/// The matrix has size (q-1)*L x (q-1)*L.
fn compute_covariance(
  fi : Array[Array[Double]],
  fij : Array[Array[Array[Array[Double]]]],
  q : Int,
  l : Int,
) -> Matrix {
  let qm = q - 1
  let size = qm * l
  let data = Array::make(size * size, 0.0)
  for i in 0..<l {
    for a in 0..<qm {
      for j in 0..<l {
        for b in 0..<qm {
          let row = i * qm + a
          let col = j * qm + b
          data[row * size + col] = fij[i][j][a][b] - fi[i][a] * fi[j][b]
        }
      }
    }
  }
  Matrix::new(size, size, data)
}

///|
/// Extract coupling J_ij(a,b) from the inverse covariance matrix.
/// J_ij(a,b) = -inv_C_{(i,a),(j,b)} for a,b in 0..q-2.
/// For the last state (q-1), J is determined by the constraint
/// sum_a J_ij(a,b) = 0.
fn extract_couplings(
  inv_c : Matrix,
  q : Int,
  l : Int,
  i : Int,
  j : Int,
) -> Array[Array[Double]] {
  let qm = q - 1
  // J[a][b] for a,b in 0..q-1
  let j_mat = Array::make(q * q, 0.0)
  // Fill for a,b in 0..q-2
  for a in 0..<qm {
    for b in 0..<qm {
      j_mat[a * q + b] = -inv_c.get(i * qm + a, j * qm + b)
    }
  }
  // Fill last row/column using sum-to-zero constraint:
  // J[q-1][b] = -sum_{a=0}^{q-2} J[a][b]
  // J[a][q-1] = -sum_{b=0}^{q-2} J[a][b]
  // J[q-1][q-1] = sum_{a,b} J[a][b]
  for b in 0..<qm {
    let mut s = 0.0
    for a in 0..<qm {
      s = s + j_mat[a * q + b]
    }
    j_mat[(q - 1) * q + b] = -s
  }
  for a in 0..<qm {
    let mut s = 0.0
    for b in 0..<qm {
      s = s + j_mat[a * q + b]
    }
    j_mat[a * q + (q - 1)] = -s
  }
  let mut s = 0.0
  for a in 0..<qm {
    for b in 0..<qm {
      s = s + j_mat[a * q + b]
    }
  }
  j_mat[(q - 1) * q + (q - 1)] = s
  // Convert flat to 2D
  let result = Array(capacity=q)
  for a in 0..<q {
    let row = Array(capacity=q)
    for b in 0..<q {
      row.push(j_mat[a * q + b])
    }
    result.push(row)
  }
  result
}

///|
/// Compute Frobenius norm of a coupling matrix.
/// S_ij = sqrt(sum_{a,b} J_ij(a,b)^2)
pub fn frobenius_norm(j_mat : Array[Array[Double]]) -> Double {
  let mut s = 0.0
  for a in 0..<j_mat.length() {
    for b in 0..<j_mat[a].length() {
      let v = j_mat[a][b]
      s = s + v * v
    }
  }
  coev_sqrt(s)
}

///|
/// Compute Direct Information (DI) score for a residue pair.
/// DI = sum_{a,b} f_ij^dir(a,b) * log(f_ij^dir(a,b) / (f_i(a) * f_j(b)))
/// where f_ij^dir is the direct coupling distribution computed from
/// the Potts model couplings.
pub fn direct_information(
  fi_i : Array[Double],
  fi_j : Array[Double],
  j_mat : Array[Array[Double]],
  q : Int,
) -> Double {
  // Compute the direct-coupling probability distribution
  // p_dir(a,b) = exp(J(a,b) + h_i(a) + h_j(b)) / Z
  // where h_i, h_j are local fields determined by marginalization.
  // We use a simplified mean-field approximation:
  // p_dir(a,b) ∝ f_i(a) * f_j(b) * exp(J(a,b))
  let mut z = 0.0
  let log_vals = Array::make(q * q, 0.0)
  for a in 0..<q {
    for b in 0..<q {
      let val = @math.ln(fi_i[a].max(1.0e-10)) +
        @math.ln(fi_j[b].max(1.0e-10)) +
        j_mat[a][b]
      log_vals[a * q + b] = val
    }
  }
  // Log-sum-exp for numerical stability
  let mut max_val = log_vals[0]
  let qq = q * q
  for k in (1)..<qq {
    if log_vals[k] > max_val {
      max_val = log_vals[k]
    }
  }
  let p_dir = Array::make(q * q, 0.0)
  for k in 0..<qq {
    p_dir[k] = @math.exp(log_vals[k] - max_val)
    z = z + p_dir[k]
  }
  for k in 0..<qq {
    p_dir[k] = p_dir[k] / z
  }
  // DI = sum p_dir(a,b) * ln(p_dir(a,b) / (f_i(a) * f_j(b)))
  let mut di = 0.0
  for a in 0..<q {
    for b in 0..<q {
      let p = p_dir[a * q + b]
      if p > 0.0 {
        di = di + p * @math.ln(p / (fi_i[a] * fi_j[b]).max(1.0e-10))
      }
    }
  }
  di
}

///|
/// Run full Mean-Field DCA on an alignment.
/// Returns MI, MI-APC, DCA Frobenius norm scores, and DI scores.
pub fn mf_dca(
  alignment : Array[String],
  config : DCAConfig,
) -> DCAResult {
  let q = get_q(config.alphabet)
  let n = alignment.length()
  let l = if n > 0 { alignment[0].length() } else { 0 }
  // Encode alignment
  let encoded = encode_alignment(alignment, config.alphabet)
  // Compute sequence weights and effective M
  let weights = sequence_weights(alignment, config.reweight_threshold)
  let m_eff = effective_sequences(weights)
  // Compute frequencies
  let fi = single_site_frequencies(encoded, q, config.pseudocount, weights, m_eff)
  let fij = pair_frequencies(encoded, q, config.pseudocount, weights, m_eff)
  // Compute MI (using empirical frequencies, no pseudocount)
  let mi = compute_mi(alignment, config.alphabet)
  let mi_apc = apc_correction(mi)
  // Compute covariance matrix
  let cov = compute_covariance(fi, fij, q, l)
  // Add small regularization to diagonal for numerical stability
  let qm = q - 1
  let size = qm * l
  let reg = 1.0e-4
  for i in 0..<size {
    cov.set(i, i, cov.get(i, i) + reg)
  }
  // Invert covariance matrix
  let inv_c = mat_inverse(cov)
  // Compute DCA scores and DI scores
  let dca_scores = Array(capacity=l)
  let di_scores = Array(capacity=l)
  for i in 0..<l {
    let dca_row = Array::make(l, 0.0)
    let di_row = Array::make(l, 0.0)
    for j in 0..<l {
      if i == j {
        continue
      }
      let j_mat = extract_couplings(inv_c, q, l, i, j)
      dca_row[j] = frobenius_norm(j_mat)
      di_row[j] = direct_information(fi[i], fi[j], j_mat, q)
    }
    dca_scores.push(dca_row)
    di_scores.push(di_row)
  }
  DCAResult::{
    mi_scores: mi,
    mi_apc_scores: mi_apc,
    dca_scores,
    di_scores,
    m_eff,
    align_len: l,
    num_seqs: n,
  }
}

// ---------------------------------------------------------------------------
// Contact prediction and ranking
// ---------------------------------------------------------------------------

///|
/// Rank residue pairs by coupling score (descending).
/// Only includes pairs (i, j) with i < j and |i - j| >= min_seq_dist.
pub fn rank_couplings(
  scores : Array[Array[Double]],
  min_seq_dist : Int,
) -> Array[Contact] {
  let l = scores.length()
  let pairs = Array(capacity=0)
  for i in 0..<l {
    for j in (i + 1)..<l {
      if j - i >= min_seq_dist {
        pairs.push(Contact::{ pos_i: i, pos_j: j, score: scores[i][j] })
      }
    }
  }
  // Sort by score descending (simple insertion sort for stability)
  sort_contacts_desc(pairs)
}

///|
/// Sort contacts by score descending (in-place).
fn sort_contacts_desc(arr : Array[Contact]) -> Array[Contact] {
  let n = arr.length()
  // Simple selection sort
  for i in 0..<n {
    let mut max_idx = i
    for j in (i + 1)..<n {
      if arr[j].score > arr[max_idx].score {
        max_idx = j
      }
    }
    if max_idx != i {
      let tmp = arr[i]
      arr[i] = arr[max_idx]
      arr[max_idx] = tmp
    }
  }
  arr
}

///|
/// Predict top-k contacts from a score matrix.
/// Only includes pairs with |i - j| >= min_seq_dist.
pub fn predict_contacts(
  scores : Array[Array[Double]],
  top_k : Int,
  min_seq_dist : Int,
) -> Array[Contact] {
  let ranked = rank_couplings(scores, min_seq_dist)
  let k = if top_k < ranked.length() { top_k } else { ranked.length() }
  let result = Array(capacity=k)
  for i in 0..<k {
    result.push(ranked[i])
  }
  result
}

///|
/// Evaluate predicted contacts against known contacts.
/// Returns (precision, recall, f1).
pub fn evaluate_contacts(
  predicted : Array[Contact],
  true_contacts : Array[(Int, Int)],
) -> (Double, Double, Double) {
  let true_map : Map[Int, Unit] = Map([], capacity=0)
  for pair in true_contacts {
    let (i, j) = pair
    true_map.set(i * 100000 + j, ())
  }
  let mut tp = 0
  for c in predicted {
    if true_map.contains(c.pos_i * 100000 + c.pos_j) {
      tp = tp + 1
    }
  }
  let precision = if predicted.length() > 0 {
    tp.to_double() / predicted.length().to_double()
  } else {
    0.0
  }
  let recall = if true_contacts.length() > 0 {
    tp.to_double() / true_contacts.length().to_double()
  } else {
    0.0
  }
  let f1 = if precision + recall > 0.0 {
    2.0 * precision * recall / (precision + recall)
  } else {
    0.0
  }
  (precision, recall, f1)
}

// ---------------------------------------------------------------------------
// DCAResult accessors
// ---------------------------------------------------------------------------

///|
/// Get the L x L MI matrix.
pub fn DCAResult::get_mi(self : DCAResult) -> Array[Array[Double]] {
  self.mi_scores
}

///|
/// Get the L x L MI-APC matrix.
pub fn DCAResult::get_mi_apc(self : DCAResult) -> Array[Array[Double]] {
  self.mi_apc_scores
}

///|
/// Get the L x L DCA (Frobenius norm) matrix.
pub fn DCAResult::get_dca(self : DCAResult) -> Array[Array[Double]] {
  self.dca_scores
}

///|
/// Get the L x L Direct Information matrix.
pub fn DCAResult::get_di(self : DCAResult) -> Array[Array[Double]] {
  self.di_scores
}

///|
/// Get the effective number of sequences.
pub fn DCAResult::get_m_eff(self : DCAResult) -> Double {
  self.m_eff
}

///|
/// Predict top-k contacts from DCA scores.
pub fn DCAResult::predict_contacts(
  self : DCAResult,
  top_k : Int,
  min_seq_dist : Int,
) -> Array[Contact] {
  predict_contacts(self.dca_scores, top_k, min_seq_dist)
}

///|
/// Rank all residue pairs by DCA score.
pub fn DCAResult::ranked_pairs(
  self : DCAResult,
  min_seq_dist : Int,
) -> Array[Contact] {
  rank_couplings(self.dca_scores, min_seq_dist)
}

///|
/// Get the top-scoring pair.
pub fn DCAResult::top_pair(self : DCAResult, min_seq_dist : Int) -> Contact? {
  let ranked = rank_couplings(self.dca_scores, min_seq_dist)
  if ranked.length() > 0 {
    Some(ranked[0])
  } else {
    None
  }
}
