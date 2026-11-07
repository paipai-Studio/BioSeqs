/// biopython_seq_distance.mbt — Sequence distance metrics
/// Inspired by Biopython: Hamming, Levenshtein, LCS, Jukes-Cantor, Kimura

///|
/// Hamming distance: count positions where two sequences differ.
/// Sequences must be equal length.
pub fn biopython_hamming_distance(seq1 : String, seq2 : String) -> Int32? {
  if seq1.length() != seq2.length() {
    return None
  }
  let mut dist = 0
  for i = 0; i < seq1.length(); i = i + 1 {
    if seq1[i] != seq2[i] {
      dist = dist + 1
    }
  }
  Some(dist)
}

///|
/// Levenshtein (edit) distance: minimum insertions/deletions/substitutions.
pub fn biopython_levenshtein_distance(seq1 : String, seq2 : String) -> Int32 {
  let m = seq1.length()
  let n = seq2.length()
  if m == 0 {
    return n
  }
  if n == 0 {
    return m
  }
  // DP table (two rows)
  let prev : Array[Int32] = Array::make(n + 1, 0)
  let curr : Array[Int32] = Array::make(n + 1, 0)
  for j = 0; j <= n; j = j + 1 {
    prev[j] = j
  }
  for i = 1; i <= m; i = i + 1 {
    curr[0] = i
    for j = 1; j <= n; j = j + 1 {
      let cost = if seq1[i - 1] == seq2[j - 1] { 0 } else { 1 }
      let del = prev[j] + 1
      let ins = curr[j - 1] + 1
      let sub = prev[j - 1] + cost
      curr[j] = if del <= ins && del <= sub {
        del
      } else if ins <= sub {
        ins
      } else {
        sub
      }
    }
    for j = 0; j <= n; j = j + 1 {
      prev[j] = curr[j]
    }
  }
  prev[n]
}

///|
/// Longest common subsequence length.
pub fn biopython_lcs_length(seq1 : String, seq2 : String) -> Int32 {
  let m = seq1.length()
  let n = seq2.length()
  if m == 0 || n == 0 {
    return 0
  }
  // DP table (two rows)
  let prev : Array[Int32] = Array::make(n + 1, 0)
  let curr : Array[Int32] = Array::make(n + 1, 0)
  for i = 1; i <= m; i = i + 1 {
    for j = 1; j <= n; j = j + 1 {
      if seq1[i - 1] == seq2[j - 1] {
        curr[j] = prev[j - 1] + 1
      } else {
        curr[j] = if prev[j] >= curr[j - 1] { prev[j] } else { curr[j - 1] }
      }
    }
    for j = 0; j <= n; j = j + 1 {
      prev[j] = curr[j]
    }
  }
  prev[n]
}

///|
/// Sequence similarity (fraction of matching positions).
pub fn biopython_sequence_similarity(seq1 : String, seq2 : String) -> Double {
  let m = seq1.length()
  let n = seq2.length()
  if m == 0 && n == 0 {
    return 1.0
  }
  let lcs = biopython_lcs_length(seq1, seq2)
  let max_len = if m >= n { m } else { n }
  lcs.to_double() / max_len.to_double()
}

///|
/// Jukes-Cantor distance for DNA sequences.
/// d = -3/4 * ln(1 - 4/3 * p)
/// where p is the proportion of different sites.
pub fn biopython_jukes_cantor_distance(seq1 : String, seq2 : String) -> Double? {
  if seq1.length() != seq2.length() || seq1.length() == 0 {
    return None
  }
  let mut diff = 0
  for i = 0; i < seq1.length(); i = i + 1 {
    if seq1[i] != seq2[i] {
      diff = diff + 1
    }
  }
  let p = diff.to_double() / seq1.length().to_double()
  let val = 1.0 - 4.0 / 3.0 * p
  if val <= 0.0 {
    return None
  }
  Some(-0.75 * @math.ln(val))
}

///|
/// Kimura 2-parameter distance for DNA sequences.
/// d = -1/2 * ln((1 - 2P - Q) * sqrt(1 - 2Q))
/// where P = transitions/n, Q = transversions/n
pub fn biopython_kimura2p_distance(seq1 : String, seq2 : String) -> Double? {
  let n = seq1.length()
  if n != seq2.length() || n == 0 {
    return None
  }
  let mut transitions = 0
  let mut transversions = 0
  for i = 0; i < n; i = i + 1 {
    if seq1[i] != seq2[i] {
      let a = seq1[i].to_int()
      let b = seq2[i].to_int()
      // Transitions: A<->G, C<->T
      let is_transition = (a == 65 && b == 71) ||
        (a == 71 && b == 65) ||
        (a == 97 && b == 103) ||
        (a == 103 && b == 97) ||
        (a == 67 && b == 84) ||
        (a == 84 && b == 67) ||
        (a == 99 && b == 116) ||
        (a == 116 && b == 99)
      if is_transition {
        transitions = transitions + 1
      } else {
        transversions = transversions + 1
      }
    }
  }
  let p = transitions.to_double() / n.to_double()
  let q = transversions.to_double() / n.to_double()
  let val1 = 1.0 - 2.0 * p - q
  let val2 = 1.0 - 2.0 * q
  if val1 <= 0.0 || val2 <= 0.0 {
    return None
  }
  Some(-0.5 * @math.ln(val1 * @math.pow(val2, 0.5)))
}
