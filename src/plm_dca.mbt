// ===========================================================================
// plm_dca.mbt — Pseudo-Likelihood Maximization DCA (plmDCA)
//
// plmDCA infers Potts model parameters (fields h_i and couplings J_ij) by
// maximizing the pseudo-likelihood of the MSA. Unlike mean-field DCA (mfDCA),
// which approximates the inference via covariance inversion, plmDCA directly
// optimizes the conditional log-likelihood at each site — a convex problem
// solved via gradient descent. This yields higher-accuracy couplings and
// better contact prediction, at higher computational cost.
//
// Potts model:
//   P(x) = exp(Σ_i h_i(x_i) + Σ_{i<j} J_ij(x_i, x_j)) / Z
//
// Pseudo-likelihood (per site i):
//   PL_i = (1/M_eff) Σ_s w_s log P(x_i^s | x_{\i}^s)
//   P(x_i=a | x_{\i}) = exp(g_i(a)) / Σ_b exp(g_i(b))
//   g_i(a) = h_i(a) + Σ_{j≠i} J_ij(a, x_j)
//
// Gradients:
//   ∂PL_i/∂h_i(a)     = (1/M_eff) Σ_s w_s [δ(x_i^s,a) - P_i^s(a)] - λ_h h_i(a)
//   ∂PL_i/∂J_ij(a,b)  = (1/M_eff) Σ_s w_s [δ(x_i^s,a) - P_i^s(a)] δ(x_j^s,b) - λ_J J_ij(a,b)
//
// Gauge: zero-sum (Σ_a J_ij(a,b)=0, Σ_b J_ij(a,b)=0, Σ_a h_i(a)=0)
//
// Reference:
//   Ekeberg, Lövkvist, Lan, Weigt, Aurell (2013).
//   "Improved contact prediction in proteins: Using Potts models with
//    pseudo-likelihood maximization." Phys Rev E 87, 012707.
// ===========================================================================

///|
/// Configuration for plmDCA computation.
pub(all) struct PlmDCAConfig {
  /// Alphabet: "protein" (21 states) or "dna" (5 states)
  alphabet : String
  /// Reweight threshold for sequence clustering (typically 0.8)
  reweight_threshold : Double
  /// L2 regularization strength for fields h_i (typically 0.01)
  lambda_h : Double
  /// L2 regularization strength for couplings J_ij (typically 0.01–0.05)
  lambda_j : Double
  /// Maximum gradient descent iterations per site (typically 300–500)
  max_iter : Int
  /// Learning rate / step size (typically 0.01–0.05)
  learning_rate : Double
  /// Convergence tolerance on gradient norm (typically 1e-5)
  tol : Double
}

///|
/// Result of plmDCA computation.
pub(all) struct PlmDCAResult {
  /// Coupling matrices: couplings[i][j][a][b] = J_ij(a,b)
  /// Symmetrized: J_ij = (J_ij^from_site_i + J_ij^from_site_j) / 2
  couplings : Array[Array[Array[Array[Double]]]]
  /// Fields: fields[i][a] = h_i(a)
  fields : Array[Array[Double]]
  /// L × L Frobenius norm scores: fn_scores[i][j] = sqrt(Σ_{a,b} J_ij(a,b)²)
  fn_scores : Array[Array[Double]]
  /// L × L Direct Information scores
  di_scores : Array[Array[Double]]
  /// Effective number of sequences after reweighting
  m_eff : Double
  /// Alignment length
  align_len : Int
  /// Number of sequences
  num_seqs : Int
}

///|
/// Default config for protein plmDCA.
pub fn PlmDCAConfig::default_protein() -> PlmDCAConfig {
  PlmDCAConfig::{
    alphabet: "protein",
    reweight_threshold: 0.8,
    lambda_h: 0.01,
    lambda_j: 0.01,
    max_iter: 300,
    learning_rate: 0.02,
    tol: 1.0e-5,
  }
}

///|
/// Default config for DNA plmDCA.
pub fn PlmDCAConfig::default_dna() -> PlmDCAConfig {
  PlmDCAConfig::{
    alphabet: "dna",
    reweight_threshold: 0.8,
    lambda_h: 0.01,
    lambda_j: 0.01,
    max_iter: 300,
    learning_rate: 0.05,
    tol: 1.0e-5,
  }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

///|
fn plm_sqrt(x : Double) -> Double {
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

///|
fn plm_get_alphabet(alphabet : String) -> String {
  match alphabet {
    "protein" => "ACDEFGHIKLMNPQRSTVWY-"
    "dna" => "ACGT-"
    _ => alphabet
  }
}

///|
fn plm_get_q(alphabet : String) -> Int {
  plm_get_alphabet(alphabet).length()
}

///|
fn plm_encode_char(c : Char, alphabet : String) -> Int {
  let ab = plm_get_alphabet(alphabet)
  let s = ab.to_array()
  for i in 0..<s.length() {
    if s[i] == c {
      return i
    }
  }
  s.length() - 1
}

///|
fn plm_encode_alignment(
  alignment : Array[String],
  alphabet : String,
) -> Array[Array[Int]] {
  let result = Array(capacity=alignment.length())
  for seq in alignment {
    let encoded = Array(capacity=seq.length())
    for i in 0..<seq.length() {
      let c = seq.unsafe_get(i).unsafe_to_char()
      encoded.push(plm_encode_char(c, alphabet))
    }
    result.push(encoded)
  }
  result
}

///|
/// Compute sequence weights (same logic as coevolution_dca).
fn plm_sequence_weights(
  alignment : Array[String],
  threshold : Double,
) -> Array[Double] {
  let n = alignment.length()
  let l = if n > 0 { alignment[0].length() } else { 0 }
  let weights = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let mut cluster_size = 0
    for j = 0; j < n; j = j + 1 {
      let mut identical = 0
      for k = 0; k < l; k = k + 1 {
        if alignment[i].unsafe_get(k) == alignment[j].unsafe_get(k) {
          identical = identical + 1
        }
      }
      if identical.to_double() > threshold * l.to_double() {
        cluster_size = cluster_size + 1
      }
    }
    weights[i] = 1.0 / cluster_size.to_double()
  }
  weights
}

///|
/// Sum of weights = effective number of sequences.
fn plm_effective_sequences(weights : Array[Double]) -> Double {
  let mut s = 0.0
  for w in weights {
    s = s + w
  }
  s
}

// ---------------------------------------------------------------------------
// Softmax with log-sum-exp stability
// ---------------------------------------------------------------------------

///|
/// Compute softmax of a vector, returning a new array.
/// Uses log-sum-exp for numerical stability.
fn softmax(logits : Array[Double]) -> Array[Double] {
  let q = logits.length()
  let mut max_val = logits[0]
  for i in 1..<q {
    if logits[i] > max_val {
      max_val = logits[i]
    }
  }
  let result = Array(capacity=q)
  let mut z = 0.0
  for i = 0; i < q; i = i + 1 {
    let e = @math.exp(logits[i] - max_val)
    result.push(e)
    z = z + e
  }
  for i = 0; i < q; i = i + 1 {
    result[i] = result[i] / z
  }
  result
}

// ---------------------------------------------------------------------------
// Gauge fixing
// ---------------------------------------------------------------------------

///|
/// Center a vector so its elements sum to zero (subtract mean).
fn center_vector(v : Array[Double]) -> Unit {
  let n = v.length()
  let mut s = 0.0
  for i = 0; i < n; i = i + 1 {
    s = s + v[i]
  }
  let mean = s / n.to_double()
  for i = 0; i < n; i = i + 1 {
    v[i] = v[i] - mean
  }
}

///|
/// Center a matrix so each row and each column sums to zero.
/// Applied iteratively until convergence (row/col centering don't commute).
fn plm_center_matrix(m : Array[Array[Double]]) -> Unit {
  let q = m.length()
  for _ in 0..<10 {
    // Center rows
    for a = 0; a < q; a = a + 1 {
      let mut s = 0.0
      for b = 0; b < q; b = b + 1 {
        s = s + m[a][b]
      }
      let mean = s / q.to_double()
      for b = 0; b < q; b = b + 1 {
        m[a][b] = m[a][b] - mean
      }
    }
    // Center columns
    for b = 0; b < q; b = b + 1 {
      let mut s = 0.0
      for a = 0; a < q; a = a + 1 {
        s = s + m[a][b]
      }
      let mean = s / q.to_double()
      for a = 0; a < q; a = a + 1 {
        m[a][b] = m[a][b] - mean
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Per-site optimization
// ---------------------------------------------------------------------------

///|
/// Optimize Potts model parameters for a single site i via gradient descent.
///
/// Returns (h_i, J_ij) where:
///   h_i : Array[Double] of size q (fields, zero-sum gauge)
///   J_ij : Array[Array[Array[Double]]] indexed [j][a][b] (couplings to site j)
///
/// Only stores j ≠ i (j == site_i entries are zero).
fn optimize_site(
  encoded : Array[Array[Int]],
  weights : Array[Double],
  m_eff : Double,
  q : Int,
  l : Int,
  site_i : Int,
  config : PlmDCAConfig,
) -> (Array[Double], Array[Array[Array[Double]]]) {
  let n = encoded.length()
  let lr = config.learning_rate
  let lam_h = config.lambda_h
  let lam_j = config.lambda_j
  // Initialize fields: h_i(a) = log(f_i(a)) - mean(log(f_i))
  let counts_i = Array::make(q, 0.0)
  for s = 0; s < n; s = s + 1 {
    counts_i[encoded[s][site_i]] = counts_i[encoded[s][site_i]] + weights[s]
  }
  let h = Array(capacity=q)
  for a = 0; a < q; a = a + 1 {
    let freq = counts_i[a] / m_eff
    h.push(@math.ln(freq.max(1.0e-10)))
  }
  center_vector(h)
  // Initialize couplings: J_ij(a,b) = 0 for all j, a, b
  let j_arr = Array(capacity=l)
  for j = 0; j < l; j = j + 1 {
    let jmat = Array(capacity=q)
    for a = 0; a < q; a = a + 1 {
      jmat.push(Array::make(q, 0.0))
    }
    j_arr.push(jmat)
  }
  // Gradient descent
  for iter = 0; iter < config.max_iter; iter = iter + 1 {
    // Accumulate gradients
    let grad_h = Array::make(q, 0.0)
    let grad_j = Array(capacity=l)
    for j = 0; j < l; j = j + 1 {
      let gj = Array(capacity=q)
      for a = 0; a < q; a = a + 1 {
        gj.push(Array::make(q, 0.0))
      }
      grad_j.push(gj)
    }
    // For each sequence, compute conditional probabilities and errors
    for s = 0; s < n; s = s + 1 {
      // Compute local field g(a) = h(a) + Σ_{j≠i} J_ij(a, x_j^s)
      let g = Array::make(q, 0.0)
      for a = 0; a < q; a = a + 1 {
        g[a] = h[a]
      }
      for j = 0; j < l; j = j + 1 {
        if j == site_i {
          continue
        }
        let xj = encoded[s][j]
        for a = 0; a < q; a = a + 1 {
          g[a] = g[a] + j_arr[j][a][xj]
        }
      }
      // Conditional probabilities
      let p = softmax(g)
      // Error: e(a) = δ(x_i^s, a) - P(a)
      let xi = encoded[s][site_i]
      let error = Array::make(q, 0.0)
      for a = 0; a < q; a = a + 1 {
        error[a] = (if a == xi { 1.0 } else { 0.0 }) - p[a]
      }
      // Accumulate gradients (weighted)
      let ws = weights[s]
      for a = 0; a < q; a = a + 1 {
        grad_h[a] = grad_h[a] + ws * error[a]
      }
      for j = 0; j < l; j = j + 1 {
        if j == site_i {
          continue
        }
        let xj = encoded[s][j]
        for a = 0; a < q; a = a + 1 {
          grad_j[j][a][xj] = grad_j[j][a][xj] + ws * error[a]
        }
      }
    }
    // Scale by 1/M_eff and add regularization
    let inv_m = 1.0 / m_eff
    let mut grad_norm = 0.0
    for a = 0; a < q; a = a + 1 {
      grad_h[a] = grad_h[a] * inv_m - lam_h * h[a]
      grad_norm = grad_norm + grad_h[a] * grad_h[a]
    }
    for j = 0; j < l; j = j + 1 {
      if j == site_i {
        continue
      }
      for a = 0; a < q; a = a + 1 {
        for b = 0; b < q; b = b + 1 {
          grad_j[j][a][b] = grad_j[j][a][b] * inv_m - lam_j * j_arr[j][a][b]
          grad_norm = grad_norm + grad_j[j][a][b] * grad_j[j][a][b]
        }
      }
    }
    grad_norm = plm_sqrt(grad_norm)
    // Check convergence
    if grad_norm < config.tol {
      break
    }
    // Update parameters (gradient ascent)
    for a = 0; a < q; a = a + 1 {
      h[a] = h[a] + lr * grad_h[a]
    }
    for j = 0; j < l; j = j + 1 {
      if j == site_i {
        continue
      }
      for a = 0; a < q; a = a + 1 {
        for b = 0; b < q; b = b + 1 {
          j_arr[j][a][b] = j_arr[j][a][b] + lr * grad_j[j][a][b]
        }
      }
    }
    // Re-apply gauge fixing
    center_vector(h)
    for j = 0; j < l; j = j + 1 {
      if j == site_i {
        continue
      }
      plm_center_matrix(j_arr[j])
    }
  }
  (h, j_arr)
}

// ---------------------------------------------------------------------------
// Scoring: Frobenius norm and Direct Information
// ---------------------------------------------------------------------------

///|
/// Compute Frobenius norm of a q×q coupling matrix.
/// F_ij = sqrt(Σ_{a,b} J_ij(a,b)²)
fn plm_frobenius_norm(j_mat : Array[Array[Double]]) -> Double {
  let mut s = 0.0
  for a in 0..<j_mat.length() {
    for b in 0..<j_mat[a].length() {
      let v = j_mat[a][b]
      s = s + v * v
    }
  }
  plm_sqrt(s)
}

///|
/// Compute Direct Information (DI) score for a residue pair.
/// DI = Σ_{a,b} p_dir(a,b) * ln(p_dir(a,b) / (f_i(a) * f_j(b)))
/// where p_dir(a,b) ∝ f_i(a) * f_j(b) * exp(J(a,b))
fn plm_direct_information(
  fi_i : Array[Double],
  fi_j : Array[Double],
  j_mat : Array[Array[Double]],
  q : Int,
) -> Double {
  let log_vals = Array::make(q * q, 0.0)
  for a = 0; a < q; a = a + 1 {
    for b = 0; b < q; b = b + 1 {
      log_vals[a * q + b] = @math.ln(fi_i[a].max(1.0e-10)) +
        @math.ln(fi_j[b].max(1.0e-10)) +
        j_mat[a][b]
    }
  }
  // Log-sum-exp
  let mut max_val = log_vals[0]
  let qq = q * q
  for k in 1..<qq {
    if log_vals[k] > max_val {
      max_val = log_vals[k]
    }
  }
  let p_dir = Array::make(qq, 0.0)
  let mut z = 0.0
  for k = 0; k < qq; k = k + 1 {
    p_dir[k] = @math.exp(log_vals[k] - max_val)
    z = z + p_dir[k]
  }
  for k = 0; k < qq; k = k + 1 {
    p_dir[k] = p_dir[k] / z
  }
  let mut di = 0.0
  for a = 0; a < q; a = a + 1 {
    for b = 0; b < q; b = b + 1 {
      let p = p_dir[a * q + b]
      if p > 0.0 {
        di = di + p * @math.ln(p / (fi_i[a] * fi_j[b]).max(1.0e-10))
      }
    }
  }
  di
}

///|
/// Compute empirical single-site frequencies (no pseudocount).
fn plm_single_site_freqs(
  encoded : Array[Array[Int]],
  weights : Array[Double],
  m_eff : Double,
  q : Int,
  l : Int,
) -> Array[Array[Double]] {
  let result = Array(capacity=l)
  for i = 0; i < l; i = i + 1 {
    let counts = Array::make(q, 0.0)
    for s in 0..<encoded.length() {
      counts[encoded[s][i]] = counts[encoded[s][i]] + weights[s]
    }
    let freq = Array(capacity=q)
    for a = 0; a < q; a = a + 1 {
      freq.push(counts[a] / m_eff)
    }
    result.push(freq)
  }
  result
}

// ---------------------------------------------------------------------------
// Main plmDCA entry point
// ---------------------------------------------------------------------------

///|
/// Run pseudo-likelihood maximization DCA on an alignment.
///
/// For each site i, optimizes the conditional log-likelihood of observing
/// residue x_i given all other positions, using gradient descent with L2
/// regularization. Couplings from both directions are symmetrized:
///   J_ij = (J_ij^{from site i} + J_ij^{from site j}) / 2
///
/// Returns fields, couplings, Frobenius norm scores, and DI scores.
pub fn plm_dca(
  alignment : Array[String],
  config : PlmDCAConfig,
) -> PlmDCAResult {
  let q = plm_get_q(config.alphabet)
  let n = alignment.length()
  let l = if n > 0 { alignment[0].length() } else { 0 }
  // Encode alignment
  let encoded = plm_encode_alignment(alignment, config.alphabet)
  // Compute sequence weights and effective M
  let weights = plm_sequence_weights(alignment, config.reweight_threshold)
  let m_eff = plm_effective_sequences(weights)
  // Compute empirical single-site frequencies for DI
  let fi = plm_single_site_freqs(encoded, weights, m_eff, q, l)
  // Optimize each site independently
  // all_fields[i] = h_i, all_couplings[i][j][a][b] = J_ij from site i's perspective
  let all_fields = Array(capacity=l)
  let all_couplings = Array(capacity=l)
  for i = 0; i < l; i = i + 1 {
    let (h_i, j_i) = optimize_site(encoded, weights, m_eff, q, l, i, config)
    all_fields.push(h_i)
    all_couplings.push(j_i)
  }
  // Symmetrize couplings: J_ij = (J_ij^i + J_ji^j^T) / 2
  // all_couplings[i][j][a][b] = J_ij(a,b) from optimizing site i
  // all_couplings[j][i][b][a] = J_ji(b,a) from optimizing site j
  // Symmetrized: J_ij(a,b) = (all_couplings[i][j][a][b] + all_couplings[j][i][b][a]) / 2
  let sym_couplings = Array(capacity=l)
  for i = 0; i < l; i = i + 1 {
    let row_i = Array(capacity=l)
    for j = 0; j < l; j = j + 1 {
      let jmat = Array(capacity=q)
      if i == j {
        // Diagonal: zero coupling
        for a = 0; a < q; a = a + 1 {
          jmat.push(Array::make(q, 0.0))
        }
      } else {
        for a = 0; a < q; a = a + 1 {
          let row_a = Array(capacity=q)
          for b = 0; b < q; b = b + 1 {
            let j_from_i = all_couplings[i][j][a][b]
            let j_from_j = all_couplings[j][i][b][a]
            row_a.push((j_from_i + j_from_j) * 0.5)
          }
          jmat.push(row_a)
        }
      }
      row_i.push(jmat)
    }
    sym_couplings.push(row_i)
  }
  // Compute Frobenius norm and DI scores
  let fn_scores = Array(capacity=l)
  let di_scores = Array(capacity=l)
  for i = 0; i < l; i = i + 1 {
    let fn_row = Array::make(l, 0.0)
    let di_row = Array::make(l, 0.0)
    for j = 0; j < l; j = j + 1 {
      if i == j {
        continue
      }
      fn_row[j] = plm_frobenius_norm(sym_couplings[i][j])
      di_row[j] = plm_direct_information(fi[i], fi[j], sym_couplings[i][j], q)
    }
    fn_scores.push(fn_row)
    di_scores.push(di_row)
  }
  PlmDCAResult::{
    couplings: sym_couplings,
    fields: all_fields,
    fn_scores,
    di_scores,
    m_eff,
    align_len: l,
    num_seqs: n,
  }
}

// ---------------------------------------------------------------------------
// Contact prediction
// ---------------------------------------------------------------------------

///|
/// A predicted residue-residue contact from plmDCA.
pub(all) struct PlmContact {
  pos_i : Int
  pos_j : Int
  score : Double
}

///|
/// Rank residue pairs by Frobenius norm score (descending).
/// Only includes pairs (i, j) with i < j and |i - j| >= min_seq_dist.
pub fn plm_rank_contacts(
  result : PlmDCAResult,
  min_seq_dist : Int,
) -> Array[PlmContact] {
  let l = result.align_len
  let pairs = Array(capacity=0)
  for i = 0; i < l; i = i + 1 {
    for j in (i + 1)..<l {
      if j - i >= min_seq_dist {
        pairs.push(PlmContact::{
          pos_i: i,
          pos_j: j,
          score: result.fn_scores[i][j],
        })
      }
    }
  }
  // Sort by score descending (selection sort)
  let n = pairs.length()
  for i = 0; i < n; i = i + 1 {
    let mut max_idx = i
    for j in (i + 1)..<n {
      if pairs[j].score > pairs[max_idx].score {
        max_idx = j
      }
    }
    if max_idx != i {
      let tmp = pairs[i]
      pairs[i] = pairs[max_idx]
      pairs[max_idx] = tmp
    }
  }
  pairs
}
