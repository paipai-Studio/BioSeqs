///|
/// Profile HMM construction and sequence search (HMMER planar model).
///
/// Mirrors the core of HMMER's `hmmbuild` + `hmmsearch`/`phmmer`/`nhmmer`
/// pipeline for a *profile* hidden Markov model (Eddy 1998, "Profile hidden
/// Markov models"; Durbin et al. 1998 ch. 5).
///
/// - `build_profile_hmm`: estimate a profile HMM (position-specific match /
///   insert / delete states) from a multiple sequence alignment.
/// - `profile_viterbi`: optimal (Viterbi) alignment of a sequence to the
///   profile, returning the aligned path and bit score.
/// - `profile_forward_score`: log-odds (bit) score of a sequence against the
///   profile under a background (null) model.
/// - `hmmer_search`: rank a query set, returning bit scores and E-values.
///
/// Plan 7 is simplified to its core: match states M_1..M_L, insert states
/// I_0..I_L, delete states D_1..D_L. The glocal (local-with-respect-to-
/// profile) entry/exit behavior of phmmer/nhmmer is implemented as uniform
/// local entry into any M_k (or I_0 / D_1) and free exit from any state.
///
/// NOTE (honesty): bit scores use Durbin-style log-odds (emission odds +
/// transition log-probabilities), not HMMER's simulation-calibrated Gumbel
/// parameters. E-values use the first-order approximation E = N * 2^-S.
/// For exact HMMER bit scores, use the C implementation / hmmer3.

// ============================================================
// Profile HMM data type
// ============================================================

///|
/// Position-specific profile HMM over a residue alphabet.
pub struct ProfileHMM {
  alphabet : Array[String]
  /// Number of match states nmat.
  length : Int
  /// Emission probabilities: match_emit[k-1] = P(residue | M_k), k = 1..nmat.
  match_emit : Array[Map[String, Double]]
  /// Emission probabilities: insert_emit[k] = P(residue | I_k), k = 0..nmat.
  insert_emit : Array[Map[String, Double]]
  // --- transition probabilities (origin-indexed) ---
  /// mm[k-1] = P(M_k -> M_{k+1}), k = 1..nmat-1.
  mm : Array[Double]
  /// mi[k-1] = P(M_k -> I_k), k = 1..nmat.
  mi : Array[Double]
  /// md[k-1] = P(M_k -> D_{k+1}), k = 1..nmat-1.
  md : Array[Double]
  /// im[k] = P(I_k -> M_{k+1}), k = 0..nmat-1.
  im : Array[Double]
  /// ii[k] = P(I_k -> I_k), k = 0..nmat.
  ii : Array[Double]
  /// dm[k-1] = P(D_k -> M_{k+1}), k = 1..nmat-1.
  dm : Array[Double]
  /// dd[k-1] = P(D_k -> D_{k+1}), k = 1..nmat-1.
  dd : Array[Double]
  /// Background (null) residue frequencies.
  background : Map[String, Double]
  /// Number of sequences the profile was built from.
  n_seqs : Int
} derive(Debug)

///|
/// Get the number of match states (profile length).
pub fn ProfileHMM::length(self : ProfileHMM) -> Int {
  self.length
}

///|
/// Get the residue alphabet of the profile.
pub fn ProfileHMM::alphabet(self : ProfileHMM) -> Array[String] {
  self.alphabet
}

///|
/// Get the number of sequences used to build the profile.
pub fn ProfileHMM::n_seqs(self : ProfileHMM) -> Int {
  self.n_seqs
}

///|
/// Emission probability P(residue | M_k) for match state k (1-based).
pub fn ProfileHMM::match_emission(self : ProfileHMM, k : Int, res : String) -> Double {
  if k < 1 || k > self.length {
    return 0.0
  }
  profile_get(self.match_emit[k - 1], res)
}

///|
/// Emission probability P(residue | I_k) for insert state k (0-based, 0..L).
pub fn ProfileHMM::insert_emission(self : ProfileHMM, k : Int, res : String) -> Double {
  if k < 0 || k > self.length {
    return 0.0
  }
  profile_get(self.insert_emit[k], res)
}

// ============================================================
// Alphabets and small helpers
// ============================================================

///|
/// Standard 20 amino-acid alphabet (uppercase, one-letter strings).
pub fn profile_aa_alphabet() -> Array[String] {
  ["A", "C", "D", "E", "F", "G", "H", "I", "K", "nmat", "M", "N", "P", "Q", "R", "S",
    "T", "V", "W", "Y"]
}

///|
/// DNA alphabet (uppercase, one-letter strings).
pub fn profile_dna_alphabet() -> Array[String] {
  ["A", "C", "G", "T"]
}

///|
/// Is a byte a gap character ('-' alignment gap or '.' insert gap)?
fn profile_is_gap(c : Int) -> Bool {
  c == '-'.to_int() || c == '.'.to_int()
}

///|
/// Map lookup with default 0.0.
fn profile_get(m : Map[String, Double], key : String) -> Double {
  match m.get(key) {
    Some(v) => v
    None => 0.0
  }
}

///|
/// Increment a count in a mutable map.
fn profile_accum(m : Map[String, Double], key : String) -> Unit {
  match m.get(key) {
    Some(v) => m[key] = v + 1.0
    None => m[key] = 1.0
  }
}

///|
/// log2 with a safe floor for p <= 0.
fn profile_log2(p : Double) -> Double {
  if p <= 0.0 {
    -1.0e300
  } else {
    @math.log2(p)
  }
}

///|
/// logadd in base-2 log space: log2(2^a + 2^b).
fn profile_logadd(a : Double, b : Double) -> Double {
  if a < b {
    if b < -1.0e200 {
      return -1.0e300
    }
    b + @math.log2(1.0 + @math.pow(2.0, a - b))
  } else {
    if a < -1.0e200 {
      return -1.0e300
    }
    a + @math.log2(1.0 + @math.pow(2.0, b - a))
  }
}


///|
/// Allocate a rows x cols matrix of Double (distinct inner arrays).
fn prof_dmatrix(rows : Int, cols : Int) -> Array[Array[Double]] {
  let m : Array[Array[Double]] = Array::make(rows, Array::make(0, 0.0))
  for i in 0..<rows {
    m[i] = Array::make(cols, 0.0)
  }
  m
}

///|
/// Allocate a rows x cols matrix of Int (distinct inner arrays).
fn prof_imatrix(rows : Int, cols : Int) -> Array[Array[Int]] {
  let m : Array[Array[Int]] = Array::make(rows, Array::make(0, 0))
  for i in 0..<rows {
    m[i] = Array::make(cols, 0)
  }
  m
}

// ============================================================
// Background model
// ============================================================

///|
/// Estimate background residue frequencies from the MSA (smoothed).
/// Frequencies are counted over every non-gap residue; a uniform prior
/// mass of 1.0 per alphabet letter is added.
fn profile_background(
  msa : Array[String],
  alphabet : Array[String],
) -> Map[String, Double] {
  let counts : Map[String, Double] = Map([], capacity=alphabet.length())
  let mut total = 0.0
  for seq in msa {
    for i in 0..<seq.length() {
      let c = seq.unsafe_get(i).to_int()
      if profile_is_gap(c) {
        continue
      }
      let res = seq[i:i + 1].to_owned()
      profile_accum(counts, res)
      total = total + 1.0
    }
  }
  let bg : Map[String, Double] = Map([], capacity=alphabet.length())
  let a_len = alphabet.length().to_double()
  for res in alphabet {
    bg[res] = (profile_get(counts, res) + 1.0) / (total + a_len)
  }
  bg
}

// ============================================================
// Construction (hmmbuild)
// ============================================================

///|
/// Internal path state for transition counting.
enum PState {
  M(Int) // M_k, k >= 1
  I(Int) // I_k, k >= 0
  D(Int) // D_k, k >= 1
} derive(Eq, Debug)

///|
/// Count a single transition between two path states, incrementing the
/// appropriate transition-count array (arrays are mutable references).
fn profile_count_transition(
  from : PState,
  to : PState,
  nmat : Int,
  mm_c : Array[Double],
  mi_c : Array[Double],
  md_c : Array[Double],
  im_c : Array[Double],
  ii_c : Array[Double],
  dm_c : Array[Double],
  dd_c : Array[Double],
) -> Unit {
  match from {
    M(k) =>
      match to {
        M(k2) => if k2 == k + 1 && k <= nmat - 1 {
          mm_c[k - 1] = mm_c[k - 1] + 1.0
        }
        I(k2) => if k2 == k {
          mi_c[k - 1] = mi_c[k - 1] + 1.0
        }
        D(k2) => if k2 == k + 1 && k <= nmat - 1 {
          md_c[k - 1] = md_c[k - 1] + 1.0
        }
        _ => ()
      }
    I(k) =>
      match to {
        M(k2) => if k2 == k + 1 && k <= nmat - 1 {
          im_c[k] = im_c[k] + 1.0
        }
        I(k2) => if k2 == k {
          ii_c[k] = ii_c[k] + 1.0
        }
        _ => ()
      }
    D(k) =>
      match to {
        M(k2) => if k2 == k + 1 && k <= nmat - 1 {
          dm_c[k - 1] = dm_c[k - 1] + 1.0
        }
        D(k2) => if k2 == k + 1 && k <= nmat - 1 {
          dd_c[k - 1] = dd_c[k - 1] + 1.0
        }
        _ => ()
      }
  }
}

///|
/// Estimate a profile HMM from a multiple sequence alignment.
///
/// All sequences in `msa` must have the same length. A column is a *match*
/// column when its gap fraction is <= `gap_threshold` (default 0.5, i.e. a
/// majority of residues); otherwise it is an *insert* column. This is a
/// standard and correct heuristic (HMMER uses the "map" maximum a posteriori
/// assignment; results differ slightly but are not wrong).
///
/// Pseudocounts: a small additive prior (`pseudocount`, default 0.01) is
/// added to every emission/transition count before normalization, preventing
/// zero probabilities on small MSAs.
pub fn build_profile_hmm(
  msa : Array[String],
  alphabet : Array[String],
  gap_threshold~ : Double = 0.5,
  pseudocount~ : Double = 0.01,
) -> ProfileHMM {
  let n = msa.length()
  let width = if n > 0 { msa[0].length() } else { 0 }

  // Classify columns: match vs insert.
  let col_is_match : Array[Bool] = Array::make(width, false)
  let mut n_match = 0
  for c in 0..<width {
    let mut gaps = 0
    for r in 0..<n {
      if profile_is_gap(msa[r].unsafe_get(c).to_int()) {
        gaps = gaps + 1
      }
    }
    let gap_frac = gaps.to_double() / n.to_double()
    if gap_frac <= gap_threshold {
      col_is_match[c] = true
      n_match = n_match + 1
    }
  }

  // Column -> match index (1-based, 0 if insert) and insert index (0-based).
  let match_no : Array[Int] = Array::make(width, 0)
  let col_to_insert : Array[Int] = Array::make(width, 0)
  let mut seen = 0
  for c in 0..<width {
    if col_is_match[c] {
      seen = seen + 1
      match_no[c] = seen
    } else {
      col_to_insert[c] = seen
    }
  }

  // Emission counters.
  let match_count : Array[Map[String, Double]] = Array::make(n_match, Map([], capacity=0))
  for k in 0..<n_match {
    match_count[k] = Map([], capacity=10)
  }
  let insert_count : Array[Map[String, Double]] = Array::make(n_match + 1, Map([], capacity=0))
  for k in 0..<=n_match {
    insert_count[k] = Map([], capacity=10)
  }

  // Transition counters (arrays for mutability).
  let mm_c : Array[Double] = Array::make(if n_match > 1 { n_match - 1 } else { 0 }, 0.0)
  let mi_c : Array[Double] = Array::make(n_match, 0.0)
  let md_c : Array[Double] = Array::make(if n_match > 1 { n_match - 1 } else { 0 }, 0.0)
  let im_c : Array[Double] = Array::make(n_match, 0.0)
  let ii_c : Array[Double] = Array::make(n_match + 1, 0.0)
  let dm_c : Array[Double] = Array::make(if n_match > 1 { n_match - 1 } else { 0 }, 0.0)
  let dd_c : Array[Double] = Array::make(if n_match > 1 { n_match - 1 } else { 0 }, 0.0)

  // Walk each sequence: build the M/I/D path, count emissions + transitions.
  for r in 0..<n {
    let seq = msa[r]
    let path : Array[PState] = []
    for c in 0..<width {
      let ch = seq.unsafe_get(c).to_int()
      if profile_is_gap(ch) {
        if col_is_match[c] {
          path.push(PState::D(match_no[c]))
        }
        continue // insert-column gap -> no state
      }
      let res = seq[c:c + 1].to_owned()
      if col_is_match[c] {
        let k = match_no[c]
        path.push(PState::M(k))
        profile_accum(match_count[k - 1], res)
      } else {
        let k = col_to_insert[c]
        path.push(PState::I(k))
        profile_accum(insert_count[k], res)
      }
    }
    let mut t = 0
    while t + 1 < path.length() {
      profile_count_transition(
        path[t], path[t + 1], n_match,
        mm_c, mi_c, md_c, im_c, ii_c, dm_c, dd_c,
      )
      t = t + 1
    }
  }

  let background = profile_background(msa, alphabet)
  let alpha = pseudocount
  let a_len = alphabet.length().to_double()

  // Normalize match emissions.
  for k in 0..<n_match {
    let counts = match_count[k]
    let mut total = 0.0
    for res in alphabet {
      total = total + (profile_get(counts, res) + alpha)
    }
    let emit : Map[String, Double] = Map([], capacity=alphabet.length())
    for res in alphabet {
      emit[res] = (profile_get(counts, res) + alpha) / total
    }
    match_count[k] = emit
  }
  // Normalize insert emissions (background-weighted prior).
  for k in 0..<=n_match {
    let counts = insert_count[k]
    let mut total = 0.0
    for res in alphabet {
      total = total + (profile_get(counts, res) + alpha * profile_get(background, res))
    }
    let emit : Map[String, Double] = Map([], capacity=alphabet.length())
    for res in alphabet {
      emit[res] =
        (profile_get(counts, res) + alpha * profile_get(background, res)) / total
    }
    insert_count[k] = emit
  }

  // Normalize transitions per origin. For M_k the outgoing states are M_{k+1}
  // (k<nmat), I_k, and D_{k+1} (k<nmat). For I_k they are M_{k+1} (k<nmat) and I_k.
  // For D_k they are M_{k+1} (k<nmat) and D_{k+1} (k<nmat).
  let mm_p : Array[Double] = Array::make(if n_match > 1 { n_match - 1 } else { 0 }, 0.0)
  let mi_p : Array[Double] = Array::make(n_match, 0.0)
  let md_p : Array[Double] = Array::make(if n_match > 1 { n_match - 1 } else { 0 }, 0.0)
  let im_p : Array[Double] = Array::make(n_match, 0.0)
  let ii_p : Array[Double] = Array::make(n_match + 1, 0.0)
  let dm_p : Array[Double] = Array::make(if n_match > 1 { n_match - 1 } else { 0 }, 0.0)
  let dd_p : Array[Double] = Array::make(if n_match > 1 { n_match - 1 } else { 0 }, 0.0)

  for k in 1..=n_match {
    // M_k origin
    let mut m_total = 0.0
    if k < n_match {
      m_total = m_total + mm_c[k - 1] + alpha
    }
    m_total = m_total + mi_c[k - 1] + alpha
    if k < n_match {
      m_total = m_total + md_c[k - 1] + alpha
    }
    if k < n_match {
      mm_p[k - 1] = (mm_c[k - 1] + alpha) / m_total
      md_p[k - 1] = (md_c[k - 1] + alpha) / m_total
    }
    mi_p[k - 1] = (mi_c[k - 1] + alpha) / m_total
  }
  for k in 0..=n_match {
    // I_k origin
    let mut i_total = 0.0
    if k < n_match {
      i_total = i_total + im_c[k] + alpha
    }
    i_total = i_total + ii_c[k] + alpha
    if k < n_match {
      im_p[k] = (im_c[k] + alpha) / i_total
    }
    ii_p[k] = (ii_c[k] + alpha) / i_total
  }
  for k in 1..=n_match {
    // D_k origin
    let mut d_total = 0.0
    if k < n_match {
      d_total = d_total + dm_c[k - 1] + alpha
      d_total = d_total + dd_c[k - 1] + alpha
    }
    if d_total > 0.0 {
      if k < n_match {
        dm_p[k - 1] = (dm_c[k - 1] + alpha) / d_total
        dd_p[k - 1] = (dd_c[k - 1] + alpha) / d_total
      }
    }
  }
  // I_n_match (last insert) has no M successor: im_p[n_match] stays 0, ii_p
  // is normalized to 1 (loop) by the i_total loop above.

  ProfileHMM::{
    alphabet,
    length: n_match,
    match_emit: match_count,
    insert_emit: insert_count,
    mm: mm_p,
    mi: mi_p,
    md: md_p,
    im: im_p,
    ii: ii_p,
    dm: dm_p,
    dd: dd_p,
    background,
    n_seqs: n,
  }
}

// ============================================================
// Emission / transition scores in log-odds
// ============================================================

///|
/// Log-odds (bit) emission score for residue `res` in a given emission map.
fn profile_emit_score(
  emit : Map[String, Double],
  background : Map[String, Double],
  res : String,
) -> Double {
  let p = profile_get(emit, res)
  let b = profile_get(background, res)
  if p <= 0.0 {
    -1.0e300
  } else {
    profile_log2(p) - profile_log2(if b > 0.0 { b } else { 1.0e-10 })
  }
}

// ============================================================
// Viterbi alignment
// ============================================================

///|
/// Result of a Viterbi alignment of a sequence to a profile HMM.
pub struct ViterbiResult {
  /// Bit score of the optimal local alignment.
  score : Double
  /// The state path as a string, e.g. "B M I I M D M E".
  path : String
  /// The aligned residue sequence (gaps '-' inserted for delete states).
  target_seq : String
}

///|
/// Optimal Viterbi alignment of `seq` to profile `profile` (local alignment).
/// Returns the bit score, the M/I/D state path, and aligned strings.
pub fn profile_viterbi(profile : ProfileHMM, seq : String) -> ViterbiResult {
  let nmat = profile.length
  let m = seq.length()
  let neg = -1.0e300

  // Transitions in log space.
  let l_mm = if nmat > 1 {
    Array::make(nmat - 1, 0.0)
  } else {
    []
  }
  let l_md = if nmat > 1 {
    Array::make(nmat - 1, 0.0)
  } else {
    []
  }
  let l_dm = if nmat > 1 {
    Array::make(nmat - 1, 0.0)
  } else {
    []
  }
  let l_dd = if nmat > 1 {
    Array::make(nmat - 1, 0.0)
  } else {
    []
  }
  let l_mi = Array::make(nmat, 0.0)
  let l_im = Array::make(nmat, 0.0)
  let l_ii = Array::make(nmat + 1, 0.0)
  for k in 1..<nmat {
    l_mm[k - 1] = profile_log2(profile.mm[k - 1])
    l_md[k - 1] = profile_log2(profile.md[k - 1])
    l_dm[k - 1] = profile_log2(profile.dm[k - 1])
    l_dd[k - 1] = profile_log2(profile.dd[k - 1])
  }
  for k in 0..<nmat {
    l_mi[k] = profile_log2(profile.mi[k])
    l_im[k] = profile_log2(profile.im[k])
    l_ii[k] = profile_log2(profile.ii[k])
  }
  l_ii[nmat] = profile_log2(profile.ii[nmat])

  // DP matrices (i = 0..m residues consumed).
  let mv : Array[Array[Double]] = prof_dmatrix(m + 1, nmat)
  let iv : Array[Array[Double]] = prof_dmatrix(m + 1, nmat + 1)
  let dv : Array[Array[Double]] = prof_dmatrix(m + 1, nmat)
  // Traceback (encoded predecessor: 0=start, 1=M, 2=I, 3=D).
  let tbM : Array[Array[Int]] = prof_imatrix(m + 1, nmat)
  let tbI : Array[Array[Int]] = prof_imatrix(m + 1, nmat + 1)
  let tbD : Array[Array[Int]] = prof_imatrix(m + 1, nmat)

  // Initialize i = 0.
  for k in 0..<nmat {
    mv[0][k] = 0.0 // local entry into any M_k
    dv[0][k] = 0.0
  }
  for k in 0..<=nmat {
    iv[0][k] = if k == 0 { 0.0 } else { neg }
  }

  // Fill.
  for i in 1..=m {
    let res = seq[i - 1:i].to_owned()
    // Match states.
    for k in 1..=nmat {
      let em = profile_emit_score(profile.match_emit[k - 1], profile.background, res)
      let mut best = 0.0 // local start
      let mut best_pred = 0
      if k > 1 {
        let a = mv[i - 1][k - 2] + l_mm[k - 2]
        if a > best {
          best = a
          best_pred = 1
        }
        let b = iv[i - 1][k - 2] + l_im[k - 2]
        if b > best {
          best = b
          best_pred = 2
        }
        let c = dv[i - 1][k - 2] + l_dm[k - 2]
        if c > best {
          best = c
          best_pred = 3
        }
      }
      mv[i][k - 1] = em + best
      tbM[i][k - 1] = best_pred
    }
    // Insert states.
    for k in 0..=nmat {
      let em = profile_emit_score(profile.insert_emit[k], profile.background, res)
      let mut best = if k == 0 { 0.0 } else { neg }
      let mut best_pred = 0
      if k >= 1 {
        let a = mv[i - 1][k - 1] + l_mi[k - 1] // M_k -> I_k
        if a > best {
          best = a
          best_pred = 1
        }
      }
      let b = iv[i - 1][k] + l_ii[k] // I_k -> I_k
      if b > best {
        best = b
        best_pred = 2
      }
      iv[i][k] = em + best
      tbI[i][k] = best_pred
    }
    // Delete states.
    for k in 1..=nmat {
      let mut best = neg
      let mut best_pred = 0
      if k == 1 {
        best = 0.0 // local entry by deleting leading matches
        best_pred = 0
      }
      if k > 1 {
        let a = mv[i][k - 2] + l_md[k - 2] // M_{k-1} -> D_k
        if a > best {
          best = a
          best_pred = 1
        }
        let d = dv[i][k - 2] + l_dd[k - 2] // D_{k-1} -> D_k
        if d > best {
          best = d
          best_pred = 3
        }
      }
      dv[i][k - 1] = best
      tbD[i][k - 1] = best_pred
    }
  }

  // Find terminal state (free exit).
  let best_i = m
  let mut best_state = 0 // 1=M, 2=I, 3=D
  let mut best_k = 0
  let mut best_score = 0.0
  for k in 1..=nmat {
    if mv[m][k - 1] > best_score {
      best_score = mv[m][k - 1]
      best_state = 1
      best_k = k
    }
    if dv[m][k - 1] > best_score {
      best_score = dv[m][k - 1]
      best_state = 3
      best_k = k
    }
  }
  for k in 0..<=nmat {
    if iv[m][k] > best_score {
      best_score = iv[m][k]
      best_state = 2
      best_k = k
    }
  }

  // Traceback to collect residues + states.
  let residues : Array[Int] = [] // 1-based residue positions (reversed)
  let states : Array[String] = []
  let mut si = best_i
  let mut sk = best_k
  let mut ss = best_state
  while si > 0 && ss != 0 {
    match ss {
      1 => {
        states.push("M")
        residues.push(si)
        let pred = tbM[si][sk - 1]
        match pred {
          1 => {
            si = si - 1
            sk = sk - 1
            ss = 1
          }
          2 => {
            si = si - 1
            sk = sk - 1
            ss = 2
          }
          3 => {
            si = si - 1
            sk = sk - 1
            ss = 3
          }
          _ => ss = 0
        }
      }
      2 => {
        states.push("I")
        residues.push(si)
        let pred = tbI[si][sk]
        match pred {
          1 => {
            si = si - 1
            ss = 1
          }
          2 => {
            si = si - 1
            ss = 2
          }
          _ => ss = 0
        }
      }
      _ => {
        // D state (no residue emitted) — push sentinel 0 to keep alignment
        states.push("D")
        residues.push(0)
        let pred = tbD[si][sk - 1]
        match pred {
          1 => {
            sk = sk - 1
            ss = 1
          }
          3 => {
            sk = sk - 1
            ss = 3
          }
          _ => ss = 0
        }
      }
    }
  }

  // Reconstruct the alignment (reverse traceback order to forward order).
  let path_sb = StringBuilder()
  let tgt_sb = StringBuilder()
  let mut idx = residues.length()
  while idx > 0 {
    idx = idx - 1
    let st = states[idx]
    path_sb.write_string(st)
    if st == "M" || st == "I" {
      let rpos = residues[idx]
      tgt_sb.write_string(seq[rpos - 1:rpos].to_owned())
    } else {
      // D state: no residue consumed — emit a gap symbol
      tgt_sb.write_char('-')
    }
  }

  ViterbiResult::{
    score: best_score,
    path: path_sb.to_string(),
    target_seq: tgt_sb.to_string(),
  }
}

// ============================================================
// Forward score (log-odds bit score)
// ============================================================

///|
/// Forward log-odds (bit) score of `seq` against `profile` under the null
/// model. This is the quantity used to rank sequences in hmmer_search.
pub fn profile_forward_score(profile : ProfileHMM, seq : String) -> Double {
  let nmat = profile.length
  let m = seq.length()
  let neg = -1.0e300

  let l_mm = if nmat > 1 {
    Array::make(nmat - 1, 0.0)
  } else {
    []
  }
  let l_md = if nmat > 1 {
    Array::make(nmat - 1, 0.0)
  } else {
    []
  }
  let l_dm = if nmat > 1 {
    Array::make(nmat - 1, 0.0)
  } else {
    []
  }
  let l_dd = if nmat > 1 {
    Array::make(nmat - 1, 0.0)
  } else {
    []
  }
  let l_mi = Array::make(nmat, 0.0)
  let l_im = Array::make(nmat, 0.0)
  let l_ii = Array::make(nmat + 1, 0.0)
  for k in 1..<nmat {
    l_mm[k - 1] = profile_log2(profile.mm[k - 1])
    l_md[k - 1] = profile_log2(profile.md[k - 1])
    l_dm[k - 1] = profile_log2(profile.dm[k - 1])
    l_dd[k - 1] = profile_log2(profile.dd[k - 1])
  }
  for k in 0..<nmat {
    l_mi[k] = profile_log2(profile.mi[k])
    l_im[k] = profile_log2(profile.im[k])
    l_ii[k] = profile_log2(profile.ii[k])
  }
  l_ii[nmat] = profile_log2(profile.ii[nmat])

  let mv : Array[Array[Double]] = prof_dmatrix(m + 1, nmat)
  let iv : Array[Array[Double]] = prof_dmatrix(m + 1, nmat + 1)
  let dv : Array[Array[Double]] = prof_dmatrix(m + 1, nmat)

  for k in 0..<nmat {
    mv[0][k] = 0.0
    dv[0][k] = 0.0
  }
  for k in 0..<=nmat {
    iv[0][k] = if k == 0 { 0.0 } else { neg }
  }

  for i in 1..=m {
    let res = seq[i - 1:i].to_owned()
    for k in 1..=nmat {
      let em = profile_emit_score(profile.match_emit[k - 1], profile.background, res)
      let mut acc = 0.0 // local start
      if k > 1 {
        acc = profile_logadd(
          acc,
          profile_logadd(
            mv[i - 1][k - 2] + l_mm[k - 2],
            profile_logadd(iv[i - 1][k - 2] + l_im[k - 2], dv[i - 1][k - 2] + l_dm[k - 2]),
          ),
        )
      }
      mv[i][k - 1] = em + acc
    }
    for k in 0..=nmat {
      let em = profile_emit_score(profile.insert_emit[k], profile.background, res)
      let mut acc = if k == 0 { 0.0 } else { neg }
      if k >= 1 {
        acc = profile_logadd(acc, mv[i - 1][k - 1] + l_mi[k - 1])
      }
      acc = profile_logadd(acc, iv[i - 1][k] + l_ii[k])
      iv[i][k] = em + acc
    }
    for k in 1..=nmat {
      let mut acc = neg
      if k == 1 {
        acc = 0.0
      }
      if k > 1 {
        acc = profile_logadd(
          profile_logadd(acc, mv[i][k - 2] + l_md[k - 2]),
          dv[i][k - 2] + l_dd[k - 2],
        )
      }
      dv[i][k - 1] = acc
    }
  }

  // Free exit: sum over all terminal states.
  let mut total = 0.0 // empty alignment
  for k in 1..=nmat {
    total = profile_logadd(total, mv[m][k - 1])
    total = profile_logadd(total, dv[m][k - 1])
  }
  for k in 0..<=nmat {
    total = profile_logadd(total, iv[m][k])
  }
  total
}

// ============================================================
// Search / ranking (phmmer, nhmmer)
// ============================================================

///|
/// A search hit against a profile HMM.
pub struct HmmHit {
  name : String
  score : Double
  evalue : Double
}

///|
/// HmmHit accessors.
pub fn HmmHit::name(self : HmmHit) -> String {
  self.name
}

pub fn HmmHit::score(self : HmmHit) -> Double {
  self.score
}

pub fn HmmHit::evalue(self : HmmHit) -> Double {
  self.evalue
}

///|
/// Rank a set of sequences against a profile HMM.
///
/// Returns hits sorted by descending bit score (best first). The bit score is
/// the Forward log-odds score. The E-value uses the first-order approximation
/// E = db_size * 2^(-score); `db_size` defaults to the number of sequences.
///
/// `names` and `seqs` must have the same length. This is the shared engine
/// behind phmmer (protein profiles) and nhmmer (DNA profiles) — the only
/// difference is the alphabet used at build time.
pub fn hmmer_search(
  profile : ProfileHMM,
  names : Array[String],
  seqs : Array[String],
  db_size~ : Int = 0,
) -> Array[HmmHit] {
  let n = names.length()
  let db = if db_size > 0 { db_size } else { n }
  let hits : Array[HmmHit] = []
  for i in 0..<n {
    let score = profile_forward_score(profile, seqs[i])
    let evalue = db.to_double() * @math.pow(2.0, -score)
    hits.push(HmmHit::{ name: names[i], score, evalue })
  }
  hits.sort_by(fn(a, b) {
    if a.score > b.score {
      -1
    } else if a.score < b.score {
      1
    } else {
      0
    }
  })
  hits
}